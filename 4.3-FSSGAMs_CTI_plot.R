###
# Project: Methods Comparison
# Data:    BOSS & BRUV fish
# Task:    Plotting the CTI models for both Abrolhos and Ningaloo
# author:  Charlotte (Claude)
# date:    September 2023
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
library(gtable)
library(ggnewscale)
library(rcartocolor)
library(cowplot)
library(patchwork)

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


## Set your working directory ----
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")
sp_dir <- paste(working.dir, "spatial_data", sep="/")


# Set the study name
study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

setwd(data_dir)
dat.abrolhos <- readRDS(paste0(name, sep="_", "dat_cti.rds"))

#### ABROLHOS ####

## Get SST
setwd(sp_dir)
zone <- "Abrolhos"
sst <- readRDS(paste0(zone, "_SST_winter.rds")) %>%
  ungroup()%>%
  dplyr::filter(!is.na(sst)) %>%
  dplyr::mutate(year = as.numeric(year))%>%
  dplyr::filter(year %in% c(2017,2018,2019,2022,2021, 2022)) %>% 
  group_by(year) %>% 
  dplyr::summarise(sst.mean = mean(sst), se = (sd(sst)/sqrt(5)))%>%
  mutate(method = "SST") %>% 
  glimpse()

#* MODEL CTI ----

mod.abrolhos <- gam(CTI~ s(sd.relief, k=3, bs='cr') + s(detrended, k=3, bs='cr') + method, family=gaussian, data=dat.abrolhos)
summary(mod.abrolhos)

gam.check(mod.abrolhos, pch=19,cex=0.8)
# predict - relief ----
testdata.abrolhos <- expand.grid(method=(mod.abrolhos$model$method),
                        biog=mean(mod.abrolhos$model$biog),
                        detrended=mean(mod.abrolhos$model$detrended)) %>%
  
  distinct()%>%
  glimpse()

fits.abrolhos <- predict.gam(mod.abrolhos, newdata=testdata.abrolhos, type='response', se.fit=T)

predicts.CTI.abrolhos = testdata.abrolhos%>%data.frame(fits.abrolhos)%>%
  group_by(method)%>% #only change here
  summarise(CTI=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>% 
  mutate(year = 2021)

## Plot for Abrolhos 

# Method ----
ggmod.CTI.Abrolhos <- ggplot() +
  ylab(NULL)+
  xlab(NULL)+
  geom_errorbar(data=predicts.CTI.abrolhos, aes(ymin =CTI-se.fit,ymax = CTI+se.fit, x=year), colour="grey20",width = 0.3) +
  geom_point(aes(x=year, y=CTI, color=method, fill=method),size=5,data=predicts.CTI.abrolhos, alpha=0.75)+
  geom_line(aes(x=year, y=sst.mean), data=sst)+
  geom_ribbon(data=sst, aes(x=year, ymin=sst.mean-se, ymax=sst.mean+se), fill="grey20", alpha=0.2)+
  ylim(19,26)+
  xlim(2017,2021.5)+ 
  scale_fill_manual(values=c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  geom_vline(xintercept=2018, linetype="dashed")+
  ylab(expression(paste("Temperature (",degree~C,")")))+
  xlab("Year")+
  #theme(legend.position = "none")+
ggmod.CTI.Abrolhos

# predict - sd.relief ----
testdata.abrolhos <- expand.grid(method=(mod.abrolhos$model$method),
                                 sd.relief=seq(min(mod.abrolhos$model$sd.relief), max(mod.abrolhos$model$sd.relief), length=20),
                                 detrended=mean(mod.abrolhos$model$detrended)) %>%
  
  distinct()%>%
  glimpse()

fits.abrolhos <- predict.gam(mod.abrolhos, newdata=testdata.abrolhos, type='response', se.fit=T)

predicts.CTI.abrolhos.sd.relief = testdata.abrolhos%>%data.frame(fits.abrolhos)%>%
  group_by(sd.relief)%>% #only change here
  summarise(CTI=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>% 
  mutate(year = 2021) 


ggmod.CTI.abrolhos.sd.relief <- ggplot() +
  ylab("CTI")+
  xlab("Standard Deviation of Relief")+
  geom_jitter(width = 0.25,height = 0)+
  #geom_point(data=use.dat, aes(x=sd.relief.relief, y=Abundance), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.CTI.abrolhos.sd.relief, aes(x=sd.relief, y=CTI), colour='grey20', alpha=0.75)+
  geom_line(data=predicts.CTI.abrolhos.sd.relief, aes(x=sd.relief, y=CTI - se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.CTI.abrolhos.sd.relief, aes(x=sd.relief, y=CTI + se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  #xlim(0, 0.9)+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") +
  ggplot2::annotate("text", x=0.01, y=25, label="(a)", size = 4, fontface=1)
ggmod.CTI.abrolhos.sd.relief

# predict - detrended ----
testdata.abrolhos <- expand.grid(method=(mod.abrolhos$model$method),
                                 detrended=seq(min(mod.abrolhos$model$detrended), max(mod.abrolhos$model$detrended), length=20),
                                 sd.relief=mean(mod.abrolhos$model$sd.relief)) %>%
  
  distinct()%>%
  glimpse()

fits.abrolhos <- predict.gam(mod.abrolhos, newdata=testdata.abrolhos, type='response', se.fit=T)

predicts.CTI.abrolhos.detrended = testdata.abrolhos%>%data.frame(fits.abrolhos)%>%
  group_by(detrended)%>% #only change here
  summarise(CTI=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>% 
  mutate(year = 2021) 


ggmod.CTI.abrolhos.detrended <- ggplot() +
  ylab("CTI")+
  xlab("Detrended Bathymetry (m)")+
  geom_jitter(width = 0.25,height = 0)+
  #geom_point(data=use.dat, aes(x=biog.relief, y=Abundance), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.CTI.abrolhos.detrended, aes(x=detrended, y=CTI), colour='grey20', alpha=0.75)+
  geom_line(data=predicts.CTI.abrolhos.detrended, aes(x=detrended, y=CTI - se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.CTI.abrolhos.detrended, aes(x=detrended, y=CTI + se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  # xlim(0, 30)+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") +
  ggplot2::annotate("text", x=-30, y=25, label="(b)", size = 4, fontface=1)
ggmod.CTI.abrolhos.detrended

setwd(fig_dir)

y.label <- textGrob(expression(paste("Temperature (",degree~C,")")), gp=gpar(fontsize=13), rot=90)

cti.other <-grid.arrange(arrangeGrob(ggmod.CTI.abrolhos.sd.relief + ylab(NULL),
                                     ggmod.CTI.abrolhos.detrended + ylab(NULL)),
                                            left=y.label,
                                            # ncol=2, 
                                            nrow=1)
ggsave(cti.other, filename="Abrolhos_CTI_other_predictors.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



#### NINGALOO #####
# Set the study name
study <- "Ningaloo_PtCloates_BOSS-BRUV"
name <- study

setwd(data_dir)
dat.ningaloo <- readRDS(paste0(name, sep="_", "dat_cti.rds")) %>% 
  filter(!campaignid %in% "2021-05_PtCloates_stereo-BRUVS" )

## Get SST
setwd(sp_dir)
zone <- "Ningaloo"
sst.ningaloo <- readRDS(paste0(zone, "_SST_winter.rds")) %>%
  ungroup()%>%
  dplyr::filter(!is.na(sst)) %>%
  dplyr::mutate(year = as.numeric(year))%>%
  dplyr::filter(year %in% c(2017,2018,2019,2022,2021,2022)) %>% 
  group_by(year) %>% 
  dplyr::summarise(sst.mean = mean(sst), se = (sd(sst)/sqrt(5)))%>%
  mutate(method = "SST") %>% 
  glimpse()

#* MODEL CTI ----

mod.ningaloo <- gam(CTI ~ s(depth, k=3, bs='cr') + s(detrended, k=3, bs='cr') + s(sdrel, k=3, bs='cr') + method, family=gaussian, data=dat.ningaloo)
summary(mod.ningaloo)

gam.check(mod.ningaloo, pch=19,cex=0.8)
# predict - relief ----
testdata.ningaloo <- expand.grid(method=(mod.ningaloo$model$method),
                                 depth=mean(mod.ningaloo$model$depth),
                                 detrended=mean(mod.ningaloo$model$detrended),
                                 sdrel=mean(mod.ningaloo$model$sdrel)) %>%
  
  distinct()%>%
  glimpse()

fits.ningaloo <- predict.gam(mod.ningaloo, newdata=testdata.ningaloo, type='response', se.fit=T)

predicts.CTI.ningaloo = testdata.ningaloo%>%data.frame(fits.ningaloo)%>%
  group_by(method)%>% #only change here
  summarise(CTI=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>% 
  mutate(across(method, factor, levels=c("BOSS","BRUV"))) %>% 
  mutate(year=2021)

## Plot for ningaloo 

# Method ----
ggmod.CTI.ningaloo <- ggplot() +
  geom_errorbar(data=predicts.CTI.ningaloo, aes(ymin =CTI-se.fit,ymax = CTI+se.fit, x=year), colour="grey20",width = 0.3) +
  geom_point(aes(x=year, y=CTI, color=method, fill=method),size=5,data=predicts.CTI.ningaloo, alpha=0.75)+
  geom_line(aes(x=year, y=sst.mean), data=sst.ningaloo)+
  geom_ribbon(data=sst.ningaloo, aes(x=year, ymin=sst.mean-se, ymax=sst.mean+se), fill="grey20", alpha=0.2)+
  ylim(19,26)+
  xlim(2017,2021.5)+ 
  scale_fill_manual(values=c("#117733", "#88CCEE"))+
  scale_colour_manual(values=c("#117733", "#88CCEE"))+
  theme_classic()+
  Theme1+
  ylab(expression(paste("Temperature (",degree~C,")")))+
  xlab("Year")+
  theme(plot.title = element_text(hjust = 0))+
  geom_vline(xintercept=2018, linetype="dashed")
  
ggmod.CTI.ningaloo

# predict - detrended ----
testdata.ningaloo <- expand.grid(method=(mod.ningaloo$model$method),
                                 depth=mean(mod.ningaloo$model$depth),
                                 detrended=seq(min(mod.ningaloo$model$detrended), max(mod.ningaloo$model$detrended), length=20),
                                 sdrel=mean(mod.ningaloo$model$sdrel)) %>%
  distinct()%>%
  glimpse()


fits.ningaloo <- predict.gam(mod.ningaloo, newdata=testdata.ningaloo, type='response', se.fit=T)

predicts.CTI.ningaloo.detrended = testdata.ningaloo%>%data.frame(fits.ningaloo)%>%
  group_by(detrended)%>% #only change here
  summarise(CTI=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>% 
  mutate(year = 2021) 


ggmod.CTI.ningaloo.detrended <- ggplot() +
  ylab("CTI")+
  xlab("Detrended Bathymetry (m)")+
  geom_jitter(width = 0.25,height = 0)+
  #geom_point(data=use.dat, aes(x=detrended.relief, y=Abundance), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.CTI.ningaloo.detrended, aes(x=detrended, y=CTI), colour='grey20', alpha=0.75)+
  geom_line(data=predicts.CTI.ningaloo.detrended, aes(x=detrended, y=CTI - se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.CTI.ningaloo.detrended, aes(x=detrended, y=CTI + se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  xlim(10, 45)+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") +
  ggplot2::annotate("text", x=10, y=25.5, label="(b)", size = 4, fontface=1)
ggmod.CTI.ningaloo.detrended

# predict - sd rel ----
testdata.ningaloo <- expand.grid(method=(mod.ningaloo$model$method),
                                 depth=mean(mod.ningaloo$model$depth),
                                 sdrel=seq(min(mod.ningaloo$model$sdrel), max(mod.ningaloo$model$sdrel), length=20),
                                 detrended=mean(mod.ningaloo$model$detrended)) %>%
  distinct()%>%
  glimpse()


fits.ningaloo <- predict.gam(mod.ningaloo, newdata=testdata.ningaloo, type='response', se.fit=T)

predicts.CTI.ningaloo.sdrel = testdata.ningaloo%>%data.frame(fits.ningaloo)%>%
  group_by(sdrel)%>% #only change here
  summarise(CTI=mean(fit),se.fit=mean(se.fit))%>%
  ungroup() %>% 
  mutate(year = 2021) 


ggmod.CTI.ningaloo.sdrel <- ggplot() +
  ylab("CTI")+
  xlab("Standard Deviation of Relief")+
  geom_jitter(width = 0.25,height = 0)+
  #geom_point(data=use.dat, aes(x=sdrel.relief, y=Abundance), colour="lightblue", alpha=0.75, size=2,show.legend=F)+
  geom_line(data=predicts.CTI.ningaloo.sdrel, aes(x=sdrel, y=CTI), colour='grey20', alpha=0.75)+
  geom_line(data=predicts.CTI.ningaloo.sdrel, aes(x=sdrel, y=CTI - se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  geom_line(data=predicts.CTI.ningaloo.sdrel, aes(x=sdrel, y=CTI + se.fit), colour='grey20', linetype="dashed",alpha=0.75)+
  #xlim(10, 45)+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(hjust = 0))+
  theme(legend.position = "none") +
  ggplot2::annotate("text", x=0.01, y=25.5, label="(a)", size = 4, fontface=1)
ggmod.CTI.ningaloo.sdrel

setwd(fig_dir)

y.label <- textGrob(expression(paste("Temperature (",degree~C,")")), gp=gpar(fontsize=13), rot=90)

cti.other <-grid.arrange(arrangeGrob(ggmod.CTI.abrolhos.sdrel + ylab(NULL),
                                     ggmod.CTI.abrolhos.detrended + ylab(NULL)),
                         left=y.label,
                         # ncol=2, 
                         nrow=1)
ggsave(cti.other, filename="Abrolhos_CTI_other_predictors.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#### PUT PLOTS TOGETHER AND SAVE #####
setwd(fig_dir)
y.label <- textGrob(expression(paste("Temperature (",degree~C,")")), gp=gpar(fontsize=13), rot=90)
x.label <- textGrob("Method", gp=gpar(fontsize=13))
legend <- gtable_filter(ggplotGrob(ggmod.CTI.Abrolhos), "guide-box")

CTI.Plots <-grid.arrange(arrangeGrob(ggmod.CTI.ningaloo + 
                                             theme(legend.position="none") +
                                             ggplot2::annotate("text", x=0.675, y=26, label="(a) Ningaloo", size = 4, fontface=1),
                                           ggmod.CTI.Abrolhos +
                                             theme(legend.position="none") + 
                                             ggplot2::annotate("text", x=0.675, y=26, label="(b) Abrolhos", size = 4, fontface=1)),
                               right=legend,
                               left=y.label,
                               bottom=x.label)


ggsave(CTI.Plots, filename=paste0("CTI_plots_both_areas.png"),height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

CTI.Plots.Other <-grid.arrange(arrangeGrob(ggmod.CTI.ningaloo.sdrel + ylab(NULL),
                                     ggmod.CTI.ningaloo.detrended + ylab(NULL)),
                         right=legend,
                         left=y.label,
                         bottom=x.label)

ggsave(CTI.Plots.Other, filename=paste0("CTI_plotsother_predictors.png"),height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



#### IMPORTANCE SCORE PLOTS ####

##### READ IN FORMATTED DATA ####
setwd(out_dir)
#read in data 
dat1 <- read.csv("2021-05_Abrolhos_BOSS-BRUV_CTI_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  mutate(resp.var = "Abrolhos") %>% 
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

dat2 <- read.csv("Ningaloo_PtCloates_BOSS-BRUV_CTI_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  mutate(resp.var = "Ningaloo") %>% 
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

CTI_imp <- bind_rows(dat1,dat2)%>%
  mutate(importance = ifelse(predictor=="depth", NA, importance)) %>% 
  glimpse()

dat.CTI.Imp <- CTI_imp %>%
  mutate(label=NA)%>%
  mutate(resp.var=factor(resp.var, levels = c("Abrolhos", "Ningaloo")))%>%
  # greater or less than length maturity
  #mutate(label=ifelse(predictor=="method"&resp.var=="Abrolhos","X",label))%>%
  mutate(label=ifelse(predictor=="sd.relief"&resp.var=="Abrolhos","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"&resp.var=="Abrolhos","X",label))%>%
  
  #mutate(label=ifelse(predictor=="depth"&resp.var=="Ningaloo","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"&resp.var=="Ningaloo","X",label))%>%
  mutate(label=ifelse(predictor=="sdrel"&resp.var=="Ningaloo","X",label))%>%

  mutate(predictor = fct_relevel(predictor,c("biog", "depth", "detrended","macroalgae", "mean.relief","relief","roughness", "sd.relief", "sdrel","sand","tpi","method"))) %>% 
  glimpse()

# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=8),
    legend.title = element_text(size=8, face="bold"),
    legend.position = "top",
    legend.direction="horizontal",
    text=element_text(size=10),
    strip.text.y = element_text(size = 6,angle = 0),
    axis.title.x=element_text(vjust=0.3, hjust = 0.5, size=6),
    axis.title.y=element_text(vjust=0.6, hjust=0.5 ,angle=90, size=8),
    axis.text.x=element_text(size=8,angle = 45, hjust=1,vjust=1),
    axis.text.y=element_text(size=7,face="italic"),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# colour ramps-
my_colors <- carto_pal(7,"Teal")
re <- colorRampPalette(my_colors)(200)

# Labels-
legend_title<-"Importance"

# Plot gg.importance.scores ----
# gg.importance.npz6 <- ggplot(dat.taxa.npz6, 
#                                aes(x=predictor,y=resp.var,fill=importance)) +
#    geom_tile(show.legend=T) +
#    scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
#                          limits = c(-1, 1))+
#    scale_y_discrete(labels=c("Smaller than legal size*","Greater than legal size*","Species richness","Total abundance"))+
#    scale_x_discrete(labels = c("Biogenic", "Depth", "Detrended", "Macroalgae", "Relief","Roughness", 'TPI'))+
#    labs(x = NULL, y = NULL) +
#    theme_classic()+
#    Theme1+
#    geom_text(aes(label=label))
# gg.importance.npz6

imp.cti.abrolhos<- ggplot(dat.CTI.Imp%>%dplyr::filter(resp.var%in%c("Abrolhos")), 
                               aes(x=predictor,y=resp.var,fill=importance)) +
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
                       limits = c(0, 1))+
  scale_y_discrete(labels=c("CTI"))+
  labs(x = NULL, y = NULL, title = "Abrolhos") +
  theme_classic()+
  Theme1+
  geom_text(aes(label=label)) +
  geom_vline(xintercept=8.5, colour="white", linewidth=3) +
  scale_x_discrete(labels = c("% Biogenic\nReef", "Depth", "Detrended\nbathymetry","% Macroalgae", "Mean Relief","Roughness","SD Relief", "TPI","Method"))+
  theme(plot.title = element_text(hjust = -0.05, vjust=-30, size=8)) # Looks crap here but title comes back in exported version
imp.cti.abrolhos

imp.cti.ningaloo <- ggplot(dat.CTI.Imp%>%dplyr::filter(resp.var%in%c("Ningaloo")), 
                                   aes(x=predictor,y=resp.var,fill=importance)) +
  geom_tile(show.legend=F) +
  scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
                       limits = c(0, 1))+
  scale_y_discrete(labels=c("CTI"))+
  scale_x_discrete(labels = c("% Biogenic\nReef","Depth","Detrended\nBathymetry", "Mean\nRelief", "Roughness", "SD Relief", "% Sand", "TPI","Method"))+
  labs(x = NULL, y = NULL, title = "Ningaloo") +
  theme_classic()+
  Theme1+
  geom_text(aes(label=label)) +
  geom_vline(xintercept=8.5, colour="white", linewidth=3) +
  theme(plot.title = element_text(hjust = -0.05, size=8)) # Looks crap here but title comes back in exported version
imp.cti.ningaloo

gg.importance.CTI <- imp.cti.abrolhos / imp.cti.ningaloo


#save output - changed dimensions for larger text in report
setwd(fig_dir)
save_plot("CTI_importance.png", gg.importance.CTI, base_height = 5, base_width = 6.275)





