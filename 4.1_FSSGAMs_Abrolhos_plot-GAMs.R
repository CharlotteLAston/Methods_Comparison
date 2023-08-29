###
# Project: Methods Comparison
# Data:    BOSS & BRUV fish
# Task:    Plotting fish GAM relationships for Abrolhos
# author:  Charlotte (Claude)
# date:    July 2023
##

rm(list=ls())

library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(GlobalArchive)
library(stringr)
library(ggplot2)
library(gamm4)
library(ggmap)
library(rgdal)
library(raster)
library(png)
library(cowplot)
library(forcats)
library(grid)
library(gridExtra)
library(ggnewscale)

# set theme
# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    legend.title = element_blank(),
    #legend.position = c(0.2, 0.8),
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=12),
    axis.title.y=element_text(vjust=0.6, angle=90, size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())
a4.width <- 160

rug.1.colour <- c("#FFFFFF", "#88CCEE")
avalues= c(0,1)
rug.1.colour = sapply(1:2,function(i)alpha(rug.1.colour[i],avalues[i]))

rug.2.colour <- c("#117733", "#FFFFFF")
avalues= c(1,0)
rug.2.colour = sapply(1:2,function(i)alpha(rug.2.colour[i],avalues[i]))

rug.3.colour <- c("white")


# Set the study name
study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

## Set working directory----
## Set your working directory ----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")


setwd(data_dir)
dat <- readRDS(paste0(name, sep="_", "dat_length.rds"))

dat.white = dat %>% 
  filter(method %in% "BRUV")

# Format data

dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity2) %>% 
  summarise(Abundance = length(Maturity2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity2", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(less_mat = ifelse(is.na(less_mat), 0, less_mat),
         greater_mat = ifelse(is.na(greater_mat), 0, greater_mat)) %>% 
  pivot_longer(cols= c(less_mat, greater_mat), names_to="Maturity2", values_to="Abundance") %>% 
  mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
  group_by(method, sample, Maturity2) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat.greater.less <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))

#### ALL SPECIES GREATER LESS THAN ####

#* MODEL < LM (macroalgae) ----

dat.small <- dat.greater.less %>% 
  filter(Maturity2 %in% c("< Length Maturity"))

mod <- gam(Abundance~ s(macroalgae, k=3, bs='cr') + method, family=tw, data=dat.small)
summary(mod)

gam.check(mod, pch=19,cex=0.8)
# predict - relief ----
testdata <- expand.grid(method=(mod$model$method),
                        macroalgae=mean(mod$model$macroalgae)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.all.less = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

## Plot for < LM 
# Method ----
ggmod.all.less<- ggplot(data=dat.small, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nfish < length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), data=predicts.all.less, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.all.less,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x=0.6, y=2.4, label="(a)", size = 4, fontface=1)
ggmod.all.less


#* MODEL > LM ----

use.dat <- dat.greater.less %>% 
  filter(Maturity2 %in% c("> Length Maturity"))

mod <- mod <- gam(Abundance~s(biog,k=3,bs='cr') + s(detrended,k=3,bs='cr') + method, family=tw, data=dat.small)
summary(mod)
gam.check(mod, pch=19,cex=0.8)

testdata <- expand.grid(method=(mod$model$method),
                        detrended=mean(mod$model$detrended),
                        biog=mean(mod$model$biog)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.all.greater = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

## Plot for > LM 

ggmod.all.greater<- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nfish > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.all.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.all.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.all.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x=0.6, y=2.4, label="(b)", size = 4, fontface=1)
ggmod.all.greater


## kde plot of length 

kde.all <- dat %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.06, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.03, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.03, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  ylim(-0.05, NA)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (mm)")+
  ggplot2::annotate("text", x=140, y=1.1, label="(c)", size = 4, fontface=1)
kde.all

setwd(fig_dir)
plot.layout <- matrix(c(1,2,
                        3,3), ncol=2, byrow=TRUE)

less.all <-grid.arrange(arrangeGrob(ggmod.all.less,
                                    ggmod.all.greater,
                                    kde.all,
                                    layout_matrix = plot.layout))

ggsave(less.all, filename="All_species_greater_less_LM.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



#### MODEL GREATER OR LESS THAN BY SPECIES ####
setwd(data_dir)
dat <- readRDS(paste0(name, sep="_", "dat_length.rds"))

dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity2) %>% 
  summarise(Abundance = length(Maturity2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity2", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(less_mat = ifelse(is.na(less_mat), 0, less_mat),
         greater_mat = ifelse(is.na(greater_mat), 0, greater_mat)) %>% 
  pivot_longer(cols= c(less_mat, greater_mat), names_to="Maturity2", values_to="Abundance") %>% 
  mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
  mutate(Mat.Species = paste0(Maturity2, sep="_", scientific)) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended", "sd.relief") %>% 
  distinct()

dat.species <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))

#* Baldchin Grouper > Maturity ####
 
use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Labridae Choerodon rubescens"))

mod=gam(Abundance~s(roughness,k=3,bs='cr') + s(sd.relief,k=3,bs='cr')  + method, family=tw, data=use.dat)
summary(mod)
gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        sd.relief=mean(mod$model$sd.relief),
                        roughness=mean(mod$model$roughness)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.baldchin.greater = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.baldchin.greater <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nC. rubescens > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8))+
  scale_x_discrete(limits = levels(predicts.baldchin.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), data=predicts.baldchin.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.baldchin.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x=0.5, y=0.78, label="(a)", size = 4, fontface=1)
ggmod.baldchin.greater

kde.baldchin <- dat %>% 
  filter(scientific %in% "Labridae Choerodon rubescens") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.06, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.03, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.03, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(100,750)+
  ylim(-0.005, NA)+
  geom_segment(x=279, y=0, xend=279, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (mm)")+
  ggplot2::annotate("text", x=100, y=0.12, label="(b)", size = 4, fontface=1) 
kde.baldchin

setwd(fig_dir)

greater.baldchin <-grid.arrange(arrangeGrob(ggmod.baldchin.greater,
                                    kde.baldchin,
                                    nrow=2))

ggsave(greater.baldchin, filename="Baldchin_greater_maturity.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#* Red throat < Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("< Length Maturity_Lethrinidae Lethrinus miniatus"))

mod=gam(Abundance~s(macroalgae,k=3,bs='cr') + s(sd.relief,k=3,bs='cr') + method, family=tw, data=use.dat)
summary(mod)
gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        macroalgae=mean(mod$model$macroalgae),
                        sd.relief = mean(mod$model$sd.relief)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.redthroat.less = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.redthroat.less <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nL. miniatus < length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))+
  scale_x_discrete(limits = levels(predicts.redthroat.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.redthroat.less, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.redthroat.less,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x=0.6, y=0.49, label="(a)", size = 4, fontface=1)
ggmod.redthroat.less

#* Red throat > Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Lethrinidae Lethrinus miniatus"))

mod=gam(Abundance~s(biog,k=3,bs='cr') + s(detrended,k=3,bs='cr') + method, family=tw, data=use.dat)
summary(mod)
gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        biog=mean(mod$model$biog),
                        detrended = mean(mod$model$detrended)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.redthroat.greater = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.redthroat.greater <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nL. miniatus > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.redthroat.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.redthroat.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.redthroat.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x=0.6, y=2.45, label="(b)", size = 4, fontface=1)
ggmod.redthroat.greater


kde.redthroat <- dat %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.06, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.03, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.03, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(200,750)+
  ylim(-0.05, NA)+
  geom_segment(x=372, y=0, xend=372, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (mm)")+
  ggplot2::annotate("text", x=210, y=0.6, label="(c)", size = 4, fontface=1) 
kde.redthroat

setwd(fig_dir)
plot.layout <- matrix(c(1,2,
                        3,3), ncol=2, byrow=TRUE)

greater.less.redthroat <-grid.arrange(arrangeGrob(ggmod.redthroat.less,
                                             ggmod.redthroat.greater,
                                            kde.redthroat,
                                            layout_matrix = plot.layout))

ggsave(greater.less.redthroat, filename="Redthroat_greater_less_maturity.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#* Pink Snapper < Maturity ####

# use.dat <- dat.species %>% 
#   filter(Mat.Species %in% ("< Length Maturity_Sparidae Chrysophrys auratus"))
# 
# mod=gam(Abundance~s(biog, k=3, bs='cr') + s(roughness, k=3, bs='cr')  + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         biog=mean(mod$model$biog),
#                         roughness=mean(mod$model$roughness)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.snapper.less = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.snapper.less <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of\nC.auratus < length at maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.snapper.less$method))+
#   geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.snapper.less, alpha=0.75)+
#   scale_fill_manual( values = c("#117733", "#88CCEE"))+
#   scale_colour_manual(values=c("#117733", "#88CCEE"))+
#   geom_errorbar(data=predicts.snapper.less, aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")+
#   ggplot2::annotate("text", x=0.5, y=1.45, label="(a)", size = 4, fontface=1)
# ggmod.snapper.less

#* Predict snapper > length maturity ####
use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Sparidae Chrysophrys auratus"))

mod=gam(Abundance~s(biog, k=3, bs='cr') + s(roughness, k=3, bs='cr') + method, family=tw, data=use.dat)
summary(mod)
gam.check(mod, pch=19, cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        biog=mean(mod$model$biog),
                        roughness=mean(mod$model$roughness)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.snapper.greater = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.snapper.greater <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nC.auratus > length at maturity")+
  xlab("Method")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.25))+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_x_discrete(limits = levels(predicts.snapper.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.snapper.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.snapper.greater, aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x=0.5, y=1.2, label="(a)", size = 4, fontface=1)
ggmod.snapper.greater

kde.snapper <- dat %>% 
  filter(scientific %in% "Sparidae Chrysophrys auratus") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.06, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.03, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.03, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(250,700)+
  ylim(-0.05, NA)+
  geom_segment(x=365, y=0, xend=365, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (mm)")+
  ggplot2::annotate("text", x=265, y=0.8, label="(b)", size = 4, fontface=1) 
kde.snapper 

setwd(fig_dir)

greater.snapper <-grid.arrange(arrangeGrob(ggmod.snapper.greater,
                                             kde.snapper,
                                             nrow=2))

ggsave(greater.snapper, filename="Snapper_greater_maturity.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

# ##### PLOTTING MODELS BY GROUPS #####
# setwd(data_dir)
# dat <- readRDS(paste0(name, sep="_", "dat_length.rds"))
# 
# # Format data
# dat.response <- dat %>% 
#   filter(status %in% "No-take") %>% 
#   group_by(method, sample, scientific, Maturity) %>% 
#   summarise(Abundance = length(Maturity)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = "Maturity", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
#   mutate(greater_mat_125 = ifelse(is.na(greater_mat_125), 0, greater_mat_125),
#          greater_mat_less_125 = ifelse(is.na(greater_mat_less_125), 0, greater_mat_less_125),
#          greater_50_less_mat = ifelse(is.na(greater_50_less_mat), 0, greater_50_less_mat)) %>% 
#   pivot_longer(cols= c(greater_mat_125, greater_mat_less_125, greater_50_less_mat), names_to="Maturity", values_to="Abundance") %>% 
#   mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
#   mutate(Maturity = fct_recode(Maturity, "< 50 Length Maturity" = "less_50", ">50 Maturity but < Maturity"="greater_50_less_mat",
#                                "> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
#                                "> 1.25x Maturity"="greater_mat_125")) %>% 
#   group_by(method, sample, Maturity) %>% 
#   summarise(Abundance = sum(Abundance)) %>% 
#   glimpse()
# 
# dat.preds <- dat %>% 
#   dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
#   distinct()
# 
# dat.groups.all <- dat.response %>% 
#   inner_join(., dat.preds, by=c("sample", "method")) %>% 
#   mutate(method = as.factor(method),
#          sample = as.factor(sample))
# 
# #* ALL SPECIES > 125% Maturity #####
# 
# use.dat <- dat.groups.all %>% 
#   filter(Maturity %in% ("> 1.25x Maturity"))
# 
# mod=gam(Abundance~s(tpi, k=3, bs='cr') + s(macroalgae, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         tpi=mean(mod$model$tpi),
#                         macroalgae=mean(mod$model$macroalgae)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.greater.125 = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.greater.125 <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of indicator species > 125% length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(-0.1, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.greater.125$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.greater.125, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.greater.125,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.greater.125
# 
# #* ALL SPECIES > Maturity but <125% Maturity #####
# 
# use.dat <- dat.groups.all %>% 
#   filter(Maturity %in% ("> Length Maturity but < 1.25x Maturity"))
# 
# mod=gam(Abundance ~  method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.greater.mat = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.greater.mat <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of indicator species > length maturity but < 125% length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.greater.mat$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.greater.mat, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.greater.mat,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.greater.mat
# 
# #* ALL SPECIES > 50% but < Maturity #####
# 
# use.dat <- dat.groups.all %>% 
#   filter(Maturity %in% (">50 Maturity but < Maturity"))
# 
# mod=gam(Abundance~s(macroalgae, k=3, bs='cr') + s(biog, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         macroalgae=mean(mod$model$macroalgae),
#                         biog=mean(mod$model$biog)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.greater.50 = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.greater.50 <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of indicator species > 50% length maturity but < length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.greater.50$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.greater.50, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.greater.50,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.greater.50
# 
# #### LENGTH GROUPS BY SPECIES ####
# 
# dat.response <- dat %>% 
#   filter(status %in% "No-take") %>% 
#   group_by(method, sample, scientific, Maturity) %>% 
#   summarise(Abundance = length(Maturity)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = "Maturity", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
#   mutate(greater_mat_125 = ifelse(is.na(greater_mat_125), 0, greater_mat_125),
#          greater_mat_less_125 = ifelse(is.na(greater_mat_less_125), 0, greater_mat_less_125),
#          greater_50_less_mat = ifelse(is.na(greater_50_less_mat), 0, greater_50_less_mat)) %>% 
#   pivot_longer(cols= c(greater_mat_125, greater_mat_less_125, greater_50_less_mat), names_to="Maturity", values_to="Abundance") %>% 
#   mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
#   mutate(Maturity = fct_recode(Maturity, "< 50 Length Maturity" = "less_50", ">50 Maturity but < Maturity"="greater_50_less_mat",
#                                "> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
#                                "> 1.25x Maturity"="greater_mat_125")) %>%
#   mutate(Mat.Species = paste0(Maturity, sep="_", scientific)) %>% 
#   glimpse()
# 
# dat.preds <- dat %>% 
#   dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
#   distinct()
# 
# dat.groups.species <- dat.response %>% 
#   inner_join(., dat.preds, by=c("sample", "method")) %>% 
#   mutate(method = as.factor(method),
#          sample = as.factor(sample))
# 
# 
# #* Baldchin > 50% Maturity but < LM #####
# 
# use.dat <- dat.groups.species %>% 
#   filter(Mat.Species %in% (">50 Maturity but < Maturity_Labridae Choerodon rubescens"))
# 
# mod=gam(Abundance~s(tpi, k=3, bs='cr') + s(depth, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         tpi=mean(mod$model$tpi),
#                         depth=mean(mod$model$depth)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.baldchin.greater.50 = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.baldchin.greater.50<- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of C. rubescens > 50% length maturity but < length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(-0.05, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.baldchin.greater.50$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.baldchin.greater.50, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.baldchin.greater.50,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.baldchin.greater.50
# 
# #* Redthroat >LM but < 125% Maturity #####
# 
# use.dat <- dat.groups.species %>% 
#   filter(Mat.Species %in% ("> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus"))
# 
# mod=gam(Abundance~s(detrended, k=3, bs='cr') + s(biog, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         detrended=mean(mod$model$detrended),
#                         biog=mean(mod$model$biog)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.redthroat.greater.mat = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.redthroat.greater.mat<- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of L. miniatus > length maturity but < 125% length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.redthroat.greater.mat$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.redthroat.greater.mat, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.redthroat.greater.mat,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.redthroat.greater.mat
# 
# #* Redthroat < 50% Maturity but < Maturity #####
# 
# use.dat <- dat.groups.species %>% 
#   filter(Mat.Species %in% (">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus"))
# 
# mod=gam(Abundance~s(detrended, k=3, bs='cr') + s(macroalgae, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         detrended=mean(mod$model$detrended),
#                         macroalgae=mean(mod$model$macroalgae)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.redthroat.greater.50 = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.redthroat.greater.50<- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of L. miniatus > 50% length maturity but < length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.redthroat.greater.50$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.redthroat.greater.50, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.redthroat.greater.50,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.redthroat.greater.50
# 
# #* Goldband > 125% Maturity #####
# 
# use.dat <- dat.groups.species %>% 
#   filter(Mat.Species %in% ("> Length Maturity but < 1.25x Maturity_Lutjanidae Pristipomoides multidens"))
# 
# mod=gam(Abundance~s(depth, k=3, bs='cr') + s(biog, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         depth=mean(mod$model$depth),
#                         biog=mean(mod$model$biog)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.goldband.greater.lm = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.goldband.greater.lm<- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of P. multidens > length maturity but < 125% length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(-0.002, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.goldband.greater.lm$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.goldband.greater.lm, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.goldband.greater.lm,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.goldband.greater.lm
# 
# #* Pink Snapper > 50% Maturity but < length maturity #####
# 
# use.dat <- dat.groups.species %>% 
#   filter(Mat.Species %in% (">50 Maturity but < Maturity_Sparidae Chrysophrys auratus"))
# 
# mod=gam(Abundance~s(biog, k=3, bs='cr') + s(macroalgae, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)
# 
# # predict method
# testdata <- expand.grid(method=(mod$model$method),
#                         macroalgae=mean(mod$model$macroalgae),
#                         biog=mean(mod$model$biog)) %>%
#   
#   distinct()%>%
#   glimpse()
# 
# fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
# 
# predicts.snapper.greater.50 = testdata%>%data.frame(fits)%>%
#   group_by(method)%>% #only change here
#   summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
#   ungroup()
# 
# ggmod.snapper.greater.50 <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
#   ylab("Predicated abundance of C. auratus > 50% length maturity but < length maturity")+
#   xlab("Method")+
#   scale_y_continuous(expand = c(0, 0), limits = c(-0.002, NA))+
#   #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
#   scale_x_discrete(limits = levels(predicts.snapper.greater.50$method))+
#   geom_point(aes(x=method, y=Abundance, color="#C77CFF", fill="#C77CFF", size=2), data=predicts.snapper.greater.50, alpha=0.75)+
#   #geom_bar(stat = "identity")+
#   geom_errorbar(data=predicts.snapper.greater.50,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
#   theme_classic()+
#   Theme1+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(legend.position = "none")
# ggmod.snapper.greater.50
# 
