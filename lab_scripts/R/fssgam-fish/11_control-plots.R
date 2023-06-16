###
# Project: Parks - Abrolhos
# Data:    BOSS fish, habitat
# Task:    Plotting species richness and greater than legal targeted species
# author:  Claude
# date:    Nov-Dec 2021
##

# Clear memory----
rm(list=ls())
gc()

# Libraries required
library(GlobalArchive)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

#load theme
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    #legend.title = element_blank(),
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=10),
    axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

## Set your working directory ----
working.dir <- getwd()
setwd(working.dir)
#OR set manually once

#standard error function
se <- function(x) sd(x)/sqrt(length(x))

# read in maxn and lengths
maxn <- readRDS("data/Tidy/dat.maxn.rds")%>%
  glimpse()

length <- readRDS("data/Tidy/dat.length.rds")%>%
  glimpse()

# read in raw maxn data 
boss <- read.csv("data/Tidy/2021-05_Abrolhos_BOSS.complete.maxn.csv")%>%
  dplyr::filter(maxn>0)%>%
  glimpse()

bruv <- read.csv("data/Tidy/2021-05_Abrolhos_stereo-BRUVs.complete.maxn.csv")%>%
  dplyr::filter(maxn>0)%>%
  glimpse()

full.maxn <- bind_rows(boss,bruv)

length(unique(full.maxn$id))

#read in SST
sst <- readRDS("data/spatial/oceanography/Abrolhos_SST_winter.rds")%>%
  ungroup()%>%
  dplyr::mutate(year=as.numeric(year))%>%
  glimpse()

locations <-  read.csv("data/spatial/oceanography/network_scale_boundaries.csv", 
                       header = TRUE) %>%
  glimpse()

sst.npz6 <- sst %>%
  dplyr::filter(Lat <= -27.916 & Lat >= -28.238,Lon <= 113.666 & Lon >= 113.112)%>%
  dplyr::group_by(year)%>%
  dplyr::summarise(sst.mean=mean(sst), sd = mean(sd))%>%
  glimpse()

str(sst.npz6)
sst.npz9 <- sst %>%
  dplyr::filter(Lat <= -26.991 & Lat >= -27.302,Lon <= 113.402 & Lon >= 112.903)%>%
  dplyr::group_by(year)%>%
  dplyr::summarise(sst.mean=mean(sst), sd = mean(sd))%>%
  glimpse()

# get rls thermal niche values ----
url <- "https://docs.google.com/spreadsheets/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?ts=5e6f36e2#gid=825736197"

master<-googlesheets4::read_sheet(url)%>%
  ga.clean.names()%>%
  filter(grepl('Australia', global.region))%>% # Change country here
  dplyr::select(family,genus,species,rls.thermal.niche)%>%
  distinct()%>%
  glimpse()

cti <- full.maxn %>%
  left_join(master)%>%
  dplyr::filter(!is.na(rls.thermal.niche))%>%
  dplyr::mutate(log.maxn=log1p(maxn),weightedSTI=log.maxn*rls.thermal.niche)%>%
  dplyr::group_by(id,sample,location,status)%>%
  dplyr::summarise(log.maxn=sum(log.maxn),w.STI = sum(weightedSTI),CTI=w.STI/log.maxn)%>%
  dplyr::ungroup()%>%
  glimpse()

#need to make a new dataframe - year, species richness (plus SE), greater than legal (plus SE)
year <- c("2017","2017","2018","2018","2019","2019","2020","2020","2021","2021","2022","2022")
status <- c("Fished","No-take")
dat <- data.frame(year,status)
dat$year <- as.numeric(dat$year)
str(dat)

#data for npz6
#species richness
spr.npz6.sr <- maxn %>%
  dplyr::filter(location%in%"NPZ6",scientific%in%"species.richness")%>%
  dplyr::group_by(status)%>%
  dplyr::summarise(species.richness = mean(maxn),species.richness.se=se(maxn))%>%
  dplyr::mutate(year=as.numeric("2021"))

#greater than legal
spr.npz6.l <- length %>%
  dplyr::filter(location%in%"NPZ6",scientific%in%"greater than legal size")%>%
  dplyr::group_by(status)%>%
  summarise(legal = mean(number),legal.se=se(number))%>%
  dplyr::mutate(year=as.numeric("2021"))

#thermal index
spr.npz6.cti <- cti %>%
  dplyr::filter(location%in%"NPZ6")%>%
  dplyr::group_by(status)%>%
  summarise(cti = mean(CTI),cti.se=se(CTI))%>%
  dplyr::mutate(year=as.numeric("2021"))

npz6 <- dat %>%
  left_join(spr.npz6.sr)%>%
  left_join(spr.npz6.l)%>%
  left_join(spr.npz6.cti)%>%
  left_join(sst.npz6)%>%
  dplyr::filter(!year=="2022")%>%
  glimpse()

#data for npz9
#no data for status
#species richness
spr.npz9.sr <- maxn %>%
  dplyr::filter(location%in%"NPZ9",scientific%in%"species.richness")%>%
  dplyr::group_by(status)%>%
  summarise(species.richness = mean(maxn),species.richness.se=se(maxn))%>%
  dplyr::mutate(year=as.numeric("2021"))

#greater than legal
spr.npz9.l <- length %>%
  dplyr::filter(location%in%"NPZ9",scientific%in%"greater than legal size")%>%
  dplyr::group_by(status)%>%
  summarise(legal = mean(number),legal.se=se(number))%>%
  dplyr::mutate(year=as.numeric("2021"))

#thermal index
spr.npz9.cti <- cti %>%
  dplyr::filter(location%in%"NPZ9")%>%
  dplyr::group_by(status)%>%
  summarise(cti = mean(CTI),cti.se=se(CTI))%>%
  dplyr::mutate(year=as.numeric("2021"))

npz9 <- dat %>%
  left_join(spr.npz9.sr)%>%
  left_join(spr.npz9.l)%>%
  left_join(spr.npz9.cti)%>%
  left_join(sst.npz9)%>%
  dplyr::filter(!year=="2022")

#NPZ6
# plot year by species richness - plus a line for MPA gazetting time ---
gg.npz6.sr <- ggplot(data = npz6, aes(x = year, y = species.richness, fill = status))+
  geom_errorbar(data = npz6,aes(ymin=species.richness-species.richness.se,
                                ymax= species.richness+species.richness.se), 
                width = 0.2,position=position_dodge(width=0.1))+
  geom_point(shape = 21,size = 2, position=position_dodge(width=0.1),stroke = 1, color = "black")+
  theme_classic()+
  scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(limits = c(2017,2021.5))+
  geom_vline(xintercept = 2018, linetype="dashed",color = "black", size=0.5,alpha = 0.5)+
  ylab("Species richness")+
  xlab("Year")+
  labs(title = "a)")+
  scale_fill_manual(labels = c("Special Purpose Zone*", "National Park Zone"),values=c("#6daff4", "#7bbc63"))+
  guides(fill=guide_legend(title = "Marine Park Zone"))+
  Theme1
gg.npz6.sr

#greater than legal - including traffic light bands
gg.npz6.l <- ggplot(data = npz6, aes(x = year, y = legal, fill = status))+
  scale_fill_manual(labels = c("Special Purpose Zone*", "National Park Zone"),values=c("#6daff4", "#7bbc63"))+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.25, ymax = 1.5),fill = "#ffeec7")+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2),fill = "#c7d6ff")+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf),fill = "#caffc7")+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.25),fill = "#ffc7c7")+
  geom_errorbar(data = npz6,aes(ymin=legal-legal.se,ymax= legal+legal.se), width = 0.2,position=position_dodge(width=0.1))+
  geom_point(shape = 21,size = 2, position=position_dodge(width=0.1),stroke = 1, color = "black")+
  theme_classic()+
  scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(limits = c(2017,2021.5))+
  geom_vline(xintercept = 2018, linetype="dashed",color = "black", size=0.5,alpha = 0.5)+
  ylab("Greater than legal size")+
  xlab("Year")+
  labs(title = "b)")+
  guides(fill=guide_legend(title = "Marine Park Zone"))+
  Theme1
gg.npz6.l

# plot year by community thermal index - plus a line for MPA gazetting time ---

gg.npz6.cti <- ggplot()+ 
  geom_line(data = npz6,aes(group = 1, x = year, y = sst.mean))+
  geom_ribbon(data = npz6,aes(group = 1, x = year, y = sst.mean, 
                              ymin = sst.mean - sd, ymax = sst.mean+sd), 
              alpha = 0.2)+
  geom_errorbar(data = npz6,aes(x = year, y = cti,ymin=cti-cti.se,
                                ymax= cti+cti.se, fill = status), 
                width = 0.2, position = position_dodge(width = 0.1))+
  geom_point(data = npz6, aes(x = year, y = cti, fill = status),shape = 21,size = 2,
             stroke = 1, color = "black", position = position_dodge(width = 0.1))+
  theme_classic()+
  scale_y_continuous(limits = c(20,23))+
  scale_x_continuous(limits = c(2017,2021.5))+
  geom_vline(xintercept = 2018, linetype="dashed",color = "black", 
             size=0.5,alpha = 0.5)+
  ylab(expression(paste("Temperature (",degree~C,")")))+
  xlab("Year")+
  scale_fill_manual(labels = c("Special Purpose Zone*", "National Park Zone"),
                    values=c("#6daff4", "#7bbc63"))+
  guides(fill=guide_legend(title = "Marine Park Zone"))+
  labs(title = "c)")+
  Theme1

gg.npz6.cti

#NPZ9
# plot year by species richness - plus a line for MPA gazetting time ---
gg.npz9.sr <- npz9 %>%
  dplyr::filter((status%in%"No-take"))%>%
  ggplot(aes(x = year, y = species.richness, fill = status))+
  geom_errorbar(aes(ymin=species.richness-species.richness.se,ymax= species.richness+species.richness.se), width = 0.1,position=position_dodge(width=0.3))+
  geom_point(shape = 21,size = 2,stroke = 1, color = "black")+
  theme_classic()+
  scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(limits = c(2017,2021.5))+
  geom_vline(xintercept = 2018, linetype="dashed",color = "black", size=0.5,alpha = 0.5)+
  ylab("Species richness")+
  xlab("Year")+
  scale_fill_manual(labels = c("National Park Zone"),values=c( "#7bbc63"))+
  guides(fill=guide_legend(title = "Marine Park Zone"))+
  labs(title = "a)")+
  Theme1
gg.npz9.sr

#greater than legal - including traffic light bands
gg.npz9.l <- npz9 %>%
  dplyr::filter((status%in%"No-take"))%>%
  ggplot(aes(x = year, y = legal))+ 
  scale_fill_manual(labels = c("National Park Zone"),values=c("#7bbc63"))+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1.5),fill = "#ffeec7")+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2),fill = "#c7d6ff")+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf),fill = "#caffc7")+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.25),fill = "#ffc7c7")+
  geom_errorbar(aes(ymin=legal-legal.se,ymax= legal+legal.se), width = 0.1,position=position_dodge(width=0.3))+
  geom_point(shape = 21,size = 2,stroke = 1, color = "black", aes(fill = status))+
  theme_classic()+
  scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(limits = c(2017,2021.5))+
  geom_vline(xintercept = 2018, linetype="dashed",color = "black", size=0.5,alpha = 0.5)+
  ylab("Greater than legal size")+
  xlab("Year")+
  guides(fill=guide_legend(title = "Marine Park Zone"))+
  labs(title = "b)")+
  Theme1
gg.npz9.l

gg.npz9.cti <- ggplot(data = npz9%>%filter((status%in%"No-take")))+ 
  geom_line(aes(group = 1, x = year, y = sst.mean))+
  geom_ribbon(aes(group = 1, x = year, y = sst.mean, 
                  ymin = sst.mean - sd, ymax = sst.mean+sd), 
              alpha = 0.2)+
  geom_errorbar(aes(x = year, y = cti,ymin=cti-cti.se,
                                ymax= cti+cti.se, fill = status), 
                width = 0.1)+
  geom_point(aes(x = year, y = cti, fill = status),shape = 21,size = 2,
             stroke = 1, color = "black")+
  theme_classic()+
  scale_y_continuous(limits = c(20,23))+
  scale_x_continuous(limits = c(2017,2021.5), expand = c(0,0))+
  geom_vline(xintercept = 2018, linetype="dashed",color = "black", 
             size=0.5,alpha = 0.5)+
  ylab(expression(paste("Temperature (",degree~C,")")))+
  xlab("Year")+
  scale_fill_manual(labels = c("National Park Zone"),values=c( "#7bbc63"))+
  guides(fill=guide_legend(title = "Marine Park Zone"))+
  labs(title = "c)")+
  Theme1

gg.npz9.cti

grid.npz6 <- gg.npz6.sr/gg.npz6.l/gg.npz6.cti+plot_layout(guides = 'collect')
grid.npz6

grid.npz9 <- gg.npz9.sr/gg.npz9.l/gg.npz9.cti+plot_layout(guides = 'collect')
grid.npz9

#save out plot
save_plot("plots/fish/control-plots.npz6.png",grid.npz6,base_height = 6,base_width = 8)
save_plot("plots/fish/control-plots.npz9.png",grid.npz9,base_height = 6,base_width = 8)


