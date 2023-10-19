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
library(gtable)

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
a4.width <- 210
a4.height <- 297

rug.1.colour <- c("#FFFFFF", "#88CCEE")
avalues= c(0,1)
rug.1.colour = sapply(1:2,function(i)alpha(rug.1.colour[i],avalues[i]))

rug.2.colour <- c("#117733", "#FFFFFF")
avalues= c(1,0)
rug.2.colour = sapply(1:2,function(i)alpha(rug.2.colour[i],avalues[i]))

rug.3.colour <- c("white")



## Set working directory----
## Set your working directory ----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

# Set the study name
study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

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
# summary(mod)

# gam.check(mod, pch=19,cex=0.8)
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
ggmod.all.less.ab<- ggplot(data=dat.small, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nfish < length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), data=predicts.all.less, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.all.less,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.all.less.ab


#* MODEL > LM ----

use.dat <- dat.greater.less %>% 
  filter(Maturity2 %in% c("> Length Maturity"))

mod <- mod <- gam(Abundance~s(biog,k=3,bs='cr') + s(detrended,k=3,bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19,cex=0.8)

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

ggmod.all.greater.ab <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nfish > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.all.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.all.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.all.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
  
ggmod.all.greater.ab

## kde plot of length 

kde.all.ab <- dat %>% 
  mutate(length=length/10) %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  ylim(-0.05, NA)+
  xlim(10,80)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (mm)")
kde.all.ab

# setwd(fig_dir)
# plot.layout <- matrix(c(1,2,
#                         3,3), ncol=2, byrow=TRUE)
# 
# less.all <-grid.arrange(arrangeGrob(ggmod.all.less,
#                                     ggmod.all.greater,
#                                     kde.all,
#                                     layout_matrix = plot.layout))
# 
# ggsave(less.all, filename="All_species_greater_less_LM.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



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

dat.species.ab <- dat.species
#* Baldchin Grouper > Maturity ####
 
use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Labridae Choerodon rubescens"))

mod=gam(Abundance~s(roughness,k=3,bs='cr') + s(sd.relief,k=3,bs='cr')  + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19,cex=0.8)

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
  theme(legend.position = "none")#+
  #ggplot2::annotate("text", x=0.5, y=0.78, label="(a)", size = 4, fontface=1)
ggmod.baldchin.greater

kde.baldchin <- dat %>% 
  mutate(length = length/10) %>% 
  filter(scientific %in% "Labridae Choerodon rubescens") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(5,80)+
  ylim(-0.005, NA)+
  geom_segment(x=27.9, y=0, xend=27.9, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)") #+
  # ggplot2::annotate("text", x=100, y=0.12, label="(b)", size = 4, fontface=1) 
kde.baldchin

# setwd(fig_dir)
# 
# greater.baldchin <-grid.arrange(arrangeGrob(ggmod.baldchin.greater,
#                                     kde.baldchin,
#                                     nrow=2))
# 
# ggsave(greater.baldchin, filename="Baldchin_greater_maturity.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
# 

#* Red throat < Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("< Length Maturity_Lethrinidae Lethrinus miniatus"))

mod=gam(Abundance~s(macroalgae,k=3,bs='cr') + s(sd.relief,k=3,bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        macroalgae=mean(mod$model$macroalgae),
                        sd.relief = mean(mod$model$sd.relief)) %>%
  
  distinct()%>%
  glimpse()

testdata <- expand.grid(method=(mod$model$method),
                        sd.relief=seq(min(mod$model$sd.relief), max(mod$model$sd.relief), length=20),
                        macroalgae = mean(mod$model$macroalgae)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.redthroat.less.ab = testdata%>%data.frame(fits)%>%
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
  theme(legend.position = "none") #+
  # ggplot2::annotate("text", x=0.6, y=0.49, label="(a)", size = 4, fontface=1)
ggmod.redthroat.less

#* Red throat > Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Lethrinidae Lethrinus miniatus"))

mod=gam(Abundance~s(biog,k=3,bs='cr') + s(detrended,k=3,bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        biog=mean(mod$model$biog),
                        detrended = mean(mod$model$detrended)) %>%
  
  distinct()%>%
  glimpse()


fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.redthroat.greater.ab = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.redthroat.greater <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nL. miniatus > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6))+
  scale_x_discrete(limits = levels(predicts.redthroat.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.redthroat.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.redthroat.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") #+
  # ggplot2::annotate("text", x=0.6, y=2.45, label="(b)", size = 4, fontface=1)
ggmod.redthroat.greater


kde.redthroat.ab <- dat %>% 
  mutate(length=length/10) %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(20,80)+
  ylim(-0.05, NA)+
  geom_segment(x=37.2, y=0, xend=37.2, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)") #+
  # ggplot2::annotate("text", x=210, y=0.6, label="(c)", size = 4, fontface=1) 
kde.redthroat.ab

# setwd(fig_dir)
# plot.layout <- matrix(c(1,2,
#                         3,3), ncol=2, byrow=TRUE)
# 
# greater.less.redthroat <-grid.arrange(arrangeGrob(ggmod.redthroat.less,
#                                              ggmod.redthroat.greater,
#                                             kde.redthroat,
#                                             layout_matrix = plot.layout))
# 
# ggsave(greater.less.redthroat, filename="Redthroat_greater_less_maturity.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

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
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)

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
  theme(legend.position = "none") #+
  # ggplot2::annotate("text", x=0.5, y=1.2, label="(a)", size = 4, fontface=1)
ggmod.snapper.greater

kde.snapper <- dat %>% 
  mutate(length=length/10) %>% 
  filter(scientific %in% "Sparidae Chrysophrys auratus") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(25,70)+
  ylim(-0.05, NA)+
  geom_segment(x=36.5, y=0, xend=36.5, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)") #+
  # ggplot2::annotate("text", x=265, y=0.8, label="(b)", size = 4, fontface=1) 
kde.snapper 

# setwd(fig_dir)
# 
# greater.snapper <-grid.arrange(arrangeGrob(ggmod.snapper.greater,
#                                              kde.snapper,
#                                              nrow=2))
# 
# ggsave(greater.snapper, filename="Snapper_greater_maturity.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### NINGALOO PLOTS ####

# Set the study name
study <- "Ningaloo_PtCloates_BOSS-BRUV" 
name <- study

setwd(data_dir)
dat <-  readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  mutate(method = as.factor(method)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  mutate(scientific = as.factor(scientific)) %>% 
  mutate(Maturity2 = as.factor(Maturity2)) %>% 
  mutate(status = as.factor(status)) %>% 
  glimpse()

dat.white <- dat %>% 
  filter(method %in% "BRUV")

# Format data

dat.response <- dat %>% 
  dplyr::filter(status %in% "No-Take") %>% 
  dplyr::group_by(method, sample, scientific, Maturity2) %>% 
  dplyr::summarise(Abundance = length(Maturity2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity2", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(less_mat = ifelse(is.na(less_mat), 0, less_mat),
         greater_mat = ifelse(is.na(greater_mat), 0, greater_mat)) %>% 
  pivot_longer(cols= c(less_mat, greater_mat), names_to="Maturity2", values_to="Abundance") %>% 
  mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
  dplyr::group_by(method, sample, Maturity2) %>% 
  dplyr::summarise(Abundance = sum(Abundance)) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "sand", "biog", "relief","tpi", "roughness", "detrended", "sdrel") %>% 
  distinct()

dat.greater.less <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))

#### ALL SPECIES GREATER LESS THAN ####

#* MODEL < LM  ----

dat.small <- dat.greater.less %>% 
  filter(Maturity2 %in% c("< Length Maturity"))

mod <- gam(Abundance~s(depth, k=3, bs='cr') + s(roughness, k=3, bs='cr') + method, family=tw, data=dat.small)
# summary(mod)

# gam.check(mod, pch=19,cex=0.8)
# predict - relief ----
testdata <- expand.grid(method=(mod$model$method),
                        depth=mean(mod$model$depth),
                        roughness=mean(mod$model$roughness)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.all.less.ni = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

## Plot for < LM 
# Method ----
ggmod.all.less.ni<- ggplot(data=dat.small, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nfish < length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), data=predicts.all.less, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  #geom_bar(stat = "identity")+
  geom_errorbar(data=predicts.all.less,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.all.less.ni


#* MODEL > LM ----

use.dat <- dat.greater.less %>% 
  filter(Maturity2 %in% c("> Length Maturity"))

mod <- mod <- gam(Abundance~s(depth, k=3,bs='cr') + s(sdrel, k=3,bs='cr') + method, family=tw, data=dat.small)
# summary(mod)
# gam.check(mod, pch=19,cex=0.8)

testdata <- expand.grid(method=(mod$model$method),
                        depth=mean(mod$model$depth),
                        sdrel=mean(mod$model$sdrel)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.all.greater.ni = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

## Plot for > LM 

ggmod.all.greater.ni<- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nfish > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5))+
  scale_x_discrete(limits = levels(predicts.all.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.all.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.all.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.all.greater.ni


## kde plot of length 

kde.all.ni <- dat %>% 
  mutate(length = length/10) %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  ylim(-0.05, NA)+
  xlim(5,80)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)")
kde.all.ni

#* Save whole indicator plots ####


## Combine the greater and less than plots into 1 plot
# Abrolhos
predicts.all.less <- predicts.all.less %>% 
  mutate(Model = "< LM")

predicts.all.greater <- predicts.all.greater %>% 
  mutate(Model = "> LM")

predicts.all <- rbind(predicts.all.less, predicts.all.greater) %>% 
  mutate(across(Model, factor, levels=c("< LM","> LM")))

ggmod.all.ab <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.all.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.all, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.all,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  facet_grid(.~Model)+
  theme_classic()+
  Theme1+
  theme(strip.text.x = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.all.ab

# Ningaloo
predicts.all.less.ni <- predicts.all.less.ni %>% 
  mutate(Model = "< LM")

predicts.all.greater.ni <- predicts.all.greater.ni %>% 
  mutate(Model = "> LM")

predicts.all.ni <- rbind(predicts.all.less.ni, predicts.all.greater.ni) %>% 
  mutate(across(Model, factor, levels=c("< LM","> LM")))

ggmod.all.ni <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.all.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.all.ni, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.all.ni,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  facet_grid(.~Model)+
  theme_classic()+
  Theme1+
  theme(strip.text.x = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.all.ni

 

setwd(fig_dir)
legend <- gtable_filter(ggplotGrob(kde.all.ni), "guide-box")
plot_layout <- rbind(c(1,1),
                     c(2,2),
                     c(3,4))


indicator.all <-grid.arrange(arrangeGrob(ggmod.all.ni + 
                                           theme(legend.position="none") +
                                           ggtitle("(a) Ningaloo") +
                                           theme(plot.title = element_text(hjust=0.05)),
                                           #ggplot2::annotate("text", x=0.5, y=4, label="(a)", size = 4, fontface=1),
                                         ggmod.all.ab + 
                                           theme(legend.position="none")+
                                           ggtitle("(b) Abrolhos") +
                                           theme(plot.title = element_text(hjust=0.05)),#+
                                           #ggplot2::annotate("text", x=13, y=11, label="(b)", size = 4, fontface=1),
                                         kde.all.ni + 
                                           theme(legend.position="none") +
                                           ggplot2::annotate("text", x=8, y=60, label="(c)", size = 4, fontface=1),
                                         kde.all.ab + 
                                           theme(legend.position="none") +
                                           ggplot2::annotate("text", x=13, y=11.75, label="(d)", size = 4, fontface=1),
                                         layout_matrix = plot_layout), right=legend)

ggsave(indicator.all, filename=paste0("All_indicator_greater_less_LM.png"), height = a4.width*1.2, width = a4.width, units  ="mm", dpi = 300 )



#### MODEL GREATER OR LESS THAN BY SPECIES ####
setwd(data_dir)

dat <- readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  mutate(method = as.factor(method)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  mutate(scientific = as.factor(scientific)) %>% 
  mutate(Maturity2 = as.factor(Maturity2)) %>% 
  mutate(status = as.factor(status)) %>% 
  glimpse()

dat.response <- dat %>% 
  filter(status %in% "No-Take") %>% 
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
  dplyr::select("method","sample","depth", "sand", "biog", "relief","tpi", "roughness", "detrended", "sdrel") %>% 
  distinct()

dat.species <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))

dat.species.ni <- dat.species

#* Redthroat > Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Lethrinidae Lethrinus miniatus"))

mod=gam(Abundance~s(detrended,k=3,bs='cr') + s(sdrel,k=3,bs='cr')  + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        sdrel=mean(mod$model$sdrel),
                        detrended=mean(mod$model$detrended)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.redthroat.greater.ni = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.redthroat.greater.ni <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nL. miniatus > length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6))+
  scale_x_discrete(limits = levels(predicts.redthroat.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), data=predicts.redthroat.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.redthroat.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") #+
  #ggplot2::annotate("text", x=0.5, y=5.9, label="(a)", size = 4, fontface=1)
ggmod.redthroat.greater.ni

kde.redthroat.ni <- dat %>% 
  mutate(length=length/10) %>% 
  filter(scientific %in% "Lethrinidae Lethrinus miniatus") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(20,75)+
  ylim(-0.05, NA)+
  geom_segment(x=37.2, y=0, xend=37.2, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)")#+
  #ggplot2::annotate("text", x=210, y=2.9, label="(b)", size = 4, fontface=1) 
kde.redthroat.ni

# setwd(fig_dir)
# 
# greater.redthroat <-grid.arrange(arrangeGrob(ggmod.redthroat.greater,
#                                              kde.redthroat,
#                                              nrow=2))
# 
# ggsave(greater.redthroat, filename=paste0(name, sep="_", "Redthroat_greater_maturity.png"), height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
# 

#* Spanglie > Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Lethrinidae Lethrinus nebulosus"))

mod=gam(Abundance~s(depth,k=3,bs='cr') + s(relief,k=3,bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
#gam.check(mod, pch=19,cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        depth=mean(mod$model$depth),
                        relief = mean(mod$model$relief)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.spanglie.greater = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.spanglie.greater <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nL. nebulosus < length at maturity")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.6))+
  scale_x_discrete(limits = levels(predicts.spanglie.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.spanglie.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.spanglie.greater,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") #+
  # ggplot2::annotate("text", x=0.6, y=1.55, label="(a)", size = 4, fontface=1)
ggmod.spanglie.greater


kde.spanglie <- dat %>%
  mutate(length = length/10) %>% 
  filter(scientific %in% "Lethrinidae Lethrinus nebulosus") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(30,65)+
  ylim(-0.05, NA)+
  geom_segment(x=35.0, y=0, xend=35.0, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)")#+
  #ggplot2::annotate("text", x=310, y=0.95, label="(b)", size = 4, fontface=1) 
kde.spanglie

# setwd(fig_dir)
# 
# greater.spanglie <-grid.arrange(arrangeGrob(ggmod.spanglie.greater,
#                                             kde.spanglie))
# 
# ggsave(greater.spanglie, filename=paste0(name, sep="_", "Spanglie_greater_less_maturity.png"),height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#* Goldband < Maturity ####

use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("< Length Maturity_Lutjanidae Pristipomoides multidens"))

mod=gam(Abundance~s(depth, k=3, bs='cr') + s(roughness, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        depth=mean(mod$model$depth),
                        roughness=mean(mod$model$roughness)) %>%
  
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.goldband.less = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.goldband.less <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nPristipomoides spp. < length at maturity")+
  xlab("Method")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2))+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_x_discrete(limits = levels(predicts.goldband.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.goldband.less, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.goldband.less, aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") #+
  # theme(axis.title.x=element_text(vjust=0.3, size=11)) #+
  # theme(axis.title.y=element_text(vjust=0.6, angle=90, size=11)) #+
  #ggplot2::annotate("text", x=0.6, y=1.9, label="(a)", size = 4, fontface=1)
ggmod.goldband.less

#* Predict snapper > length maturity ####
use.dat <- dat.species %>% 
  filter(Mat.Species %in% ("> Length Maturity_Lutjanidae Pristipomoides multidens"))

mod=gam(Abundance~s(depth, k=3, bs='cr') + s(roughness, k=3, bs='cr') + method, family=tw, data=use.dat)
# summary(mod)
# gam.check(mod, pch=19, cex=0.8)

# predict method
testdata <- expand.grid(method=(mod$model$method),
                        depth=mean(mod$model$depth),
                        roughness=mean(mod$model$roughness))

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.goldband.greater = testdata%>%data.frame(fits)%>%
  group_by(method)%>% #only change here
  summarise(Abundance=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

ggmod.goldband.greater <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance of\nPristipomoides spp. > length at maturity")+
  xlab("Method")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  scale_x_discrete(limits = levels(predicts.goldband.greater$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill="#C77CFF", size=2), predicts.goldband.greater, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.goldband.greater, aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") #+
  # theme(axis.title.x=element_text(vjust=0.3, size=11)) +
  # theme(axis.title.y=element_text(vjust=0.6, angle=90, size=11)) #+
  # ggplot2::annotate("text", x=0.6, y=0.475, label="(b)", size = 4, fontface=1)
ggmod.goldband.greater

kde.goldband <- dat %>% 
  mutate(length = length/10) %>% 
  filter(scientific %in% "Lutjanidae Pristipomoides multidens") %>% 
  ggplot() +
  geom_density(aes(length, y = stat(count), color = method, fill = method), alpha = 0.5)+
  scale_fill_manual( values = c("#117733", "#88CCEE"), name="Method")+
  scale_colour_manual(values=c("#117733", "#88CCEE"), name="Method")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method, group=method), sides = 'b', outside = F, length=unit(0.04, "npc")) +
  scale_colour_manual(values=rug.1.colour, guide="none")+
  new_scale_colour()+
  geom_rug(data=dat.white, aes(x = length, colour=method), sides = 'b', outside = F, length=unit(0.02, "npc"), size=1.1) +
  scale_colour_manual(values=rug.3.colour, guide="none")+
  new_scale_colour()+
  geom_rug(aes(x = length, color = method), sides = 'b', outside = F, length=unit(0.02, "npc")) +
  scale_colour_manual(values=rug.2.colour, guide="none")+
  xlim(15,85)+
  ylim(-0.05, 20)+
  geom_segment(x=52.6, y=0, xend=52.6, yend=Inf)+
  #geom_vline(xintercept = 600*1.25, linetype="dotted")+
  #geom_vline(xintercept = 600*0.5, linetype="dotted")+
  theme_classic()+
  Theme1+
  ylab("Count\n")+
  xlab("Length (cm)") #+
  # ggplot2::annotate("text", x=155, y=1.9, label="(c)", size = 4, fontface=1) 
kde.goldband 

# setwd(fig_dir)
# plot.layout <- matrix(c(1,2,
#                         3,3), ncol=2, byrow=TRUE)
# 
# greater.goldband <-grid.arrange(arrangeGrob(ggmod.goldband.less,
#                                             ggmod.goldband.greater,
#                                             kde.goldband,
#                                             layout_matrix = plot.layout))
# 
# ggsave(greater.goldband, filename=paste0(name, sep="_", "Goldband_greater_maturity.png"),height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### Save grouped plots ####
#* L miniatus ####
predicts.redthroat.less.ab <- predicts.redthroat.less.ab %>% 
  mutate(Model = "< LM")

predicts.redthroat.greater.ab <- predicts.redthroat.greater.ab %>% 
  mutate(Model = "> LM")

predicts.redthroat.ab <- rbind(predicts.redthroat.less.ab, predicts.redthroat.greater.ab) %>% 
  mutate(across(Model, factor, levels=c("< LM","> LM")))

use.dat.redthroat <- dat.species.ab %>% 
  filter(scientific %in% c("Lethrinidae Lethrinus miniatus"))

ggmod.redthroat.ab <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.all.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.redthroat.ab, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.redthroat.ab,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  facet_grid(.~Model)+
  theme_classic()+
  Theme1+
  theme(strip.text.x = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.redthroat.ab

# Ningaloo
# predicts.redthroat.less.ni <- predicts.redthroat.less.ni %>% 
#   mutate(Model = "< LM")

predicts.redthroat.greater.ni <- predicts.redthroat.greater.ni %>% 
  mutate(Model = "> LM")

predicts.redthroat.ni <- predicts.redthroat.greater.ni%>% 
  mutate(across(Model, factor, levels=c("< LM","> LM")))

use.dat.redthroat.ni <- dat.species.ni %>% 
  filter(scientific %in% c("Lethrinidae Lethrinus miniatus"))

ggmod.redthroat.ni <- ggplot(data=use.dat, aes(x=method, y=Abundance)) +
  ylab("Predicated abundance")+
  xlab("Method")+
  #geom_bar(aes(fill=method, colour=method), stat = "summary", fun = "mean", alpha=0.2)+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_x_discrete(limits = levels(predicts.redthroat.less$method))+
  geom_point(aes(x=method, y=Abundance, color=method, fill=method, size=2), data=predicts.redthroat.ni, alpha=0.75)+
  scale_fill_manual( values = c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  geom_errorbar(data=predicts.redthroat.ni,aes(ymin =Abundance-se.fit,ymax = Abundance+se.fit), colour="grey20",width = 0.5) +
  facet_grid(.~Model, drop=FALSE)+
  theme_classic()+
  Theme1+
  theme(strip.text.x = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none")
ggmod.redthroat.ni
 

setwd(fig_dir)
legend <- gtable_filter(ggplotGrob(kde.all.ni), "guide-box")
plot_layout <- rbind(c(1,1),
                     c(2,2),
                     c(3,4))

miniatus <-grid.arrange(arrangeGrob(ggmod.redthroat.ni + 
                                      theme(legend.position="none") +
                                      ggtitle("(a) Ningaloo") +
                                      theme(plot.title = element_text(hjust=0.05)),
                                    ggmod.redthroat.ab + 
                                      theme(legend.position="none") +
                                      ggtitle("(b) Abrolhos") +
                                      theme(plot.title = element_text(hjust=0.05)),
                                    kde.redthroat.ni+ 
                                      theme(legend.position="none") +
                                      ggplot2::annotate("text", x=23, y=29, label="(c)", size = 4, fontface=1),
                                    kde.redthroat.ab+ 
                                      theme(legend.position="none") +
                                      ggplot2::annotate("text", x=23, y=5, label="(d)", size = 4, fontface=1),
                                    layout_matrix = plot_layout,
                                    nrow=3, ncol=2), right=legend)

ggsave(miniatus, filename=paste0("Redthroat_greater_less_LM.png"), height = a4.width*1.2, width = a4.width, units  ="mm", dpi = 300 )

#* All other species put together ####
setwd(fig_dir)
legend <- gtable_filter(ggplotGrob(kde.all.ni), "guide-box")
plot_layout <- rbind(c(1,2,2),
                     c(3,4,4),
                     c(5,6,6),
                     c(7,8,9))

other.species <-grid.arrange(arrangeGrob(ggmod.baldchin.greater + 
                                      theme(legend.position="none") +
                                        theme(axis.title = element_text(size=8.5))+
                                      ggplot2::annotate("text", x=0.6, y=0.475, label="(a)", size = 4, fontface=1),
                                    kde.baldchin +
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) +
                                      ggplot2::annotate("text", x=0.6, y=0.475, label="(a)", size = 4, fontface=1),
                                    ggmod.snapper.greater + 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) +
                                      ggplot2::annotate("text", x=0.6, y=5.5, label="(c)", size = 4, fontface=1),
                                    kde.snapper+ 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) +
                                      ggplot2::annotate("text", x=23, y=5, label="(d)", size = 4, fontface=1),
                                    ggmod.spanglie.greater + 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) + 
                                      ggplot2::annotate("text", x=0.6, y=2.35, label="(b)", size = 4, fontface=1),
                                    kde.spanglie + 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) +
                                      ggplot2::annotate("text", x=23, y=27, label="(e)", size = 4, fontface=1),
                                    ggmod.goldband.less + 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) +
                                      ggplot2::annotate("text", x=0.6, y=2.35, label="(b)", size = 4, fontface=1),
                                    ggmod.goldband.greater + 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) + 
                                      ggplot2::annotate("text", x=0.6, y=2.35, label="(b)", size = 4, fontface=1),
                                    kde.goldband + 
                                      theme(legend.position="none") +
                                      theme(axis.title = element_text(size=8.5)) +
                                      ggplot2::annotate("text", x=23, y=27, label="(e)", size = 4, fontface=1),
                                    layout_matrix = plot_layout,
                                    nrow=4, ncol=3), right=legend)

ggsave(other.species, filename=paste0("All_other_greater_less_LM.png"), height = a4.width, width = a4.height, units  ="mm", dpi = 300 )

