###
# Project: Parks - Abrolhos
# Data:    BOSS & BRUV fish, habitat
# Task:    Plotting fish GAM relationships at npz6
# author:  Claude
# date:    Nov-Dec 2021
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
    legend.position = c(0.2, 0.8),
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=10),
    axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# Set the study name
name <- "2021-05_Abrolhos_npz6" # for the study

## Set working directory----
working.dir <- getwd()
setwd(working.dir)
#OR Set manually once

# Bring in and format the  data----
dat.maxn <- readRDS("data/Tidy/dat.maxn.rds")%>%
  dplyr::filter(location%in%"NPZ6")%>%
  dplyr::rename(number = maxn)%>%
  glimpse()

dat.length <- readRDS("data/Tidy/dat.length.rds")%>%
  dplyr::filter(location%in%"NPZ6")%>%
  glimpse()

dat <- bind_rows(dat.maxn,dat.length)

# Manually make the most parsimonious GAM models for each taxa ----
#### Abrolhos MaxN ####
unique(dat$scientific)

# MODEL Total abundance (relief) ----
dat.total <- dat %>% filter(scientific=="total.abundance")

mod=gam(number~s(mean.relief,k=3,bs='cr'), family=tw,data=dat.total)

# predict - relief ----
testdata <- expand.grid(mean.relief=seq(min(dat$mean.relief),max(dat$mean.relief),length.out = 20)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.total.relief = testdata%>%data.frame(fits)%>%
  group_by(mean.relief)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for Total abundance ----
# detrended bathy ----
ggmod.total.relief<- ggplot() +
  ylab("")+
  xlab("Relief")+
  geom_point(data=dat.total,aes(x=mean.relief,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.total.relief,aes(x=mean.relief,y=number),alpha=0.5)+
  geom_line(data=predicts.total.relief,aes(x=mean.relief,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.total.relief,aes(x=mean.relief,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Total abundance") +
  theme(plot.title = element_text(hjust = 0))
ggmod.total.relief

# MODEL Species richness (biog + depth + tpi) ----
dat.species <- dat %>% filter(scientific=="species.richness")

mod=gam(number~s(biog,k=3,bs='cr') + s(depth,k=3,bs='cr') + s(tpi,k=3,bs='cr'), family=tw,data=dat.species)

# predict - bioge ----
testdata <- expand.grid(biog=seq(min(dat$biog),max(dat$biog),length.out = 20),
                        depth=mean(mod$model$depth),
                        tpi=mean(mod$model$tpi)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.species.biog = testdata%>%data.frame(fits)%>%
  group_by(biog)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# predict - depth ----
testdata <- expand.grid(depth=seq(min(dat$depth),max(dat$depth),length.out = 20),
                        biog=mean(mod$model$biog),
                        tpi=mean(mod$model$tpi)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.species.depth = testdata%>%data.frame(fits)%>%
  group_by(depth)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# predict - tpi ----
testdata <- expand.grid(tpi=seq(min(dat$tpi),max(dat$tpi),length.out = 20),
                        biog=mean(mod$model$biog),
                        depth=mean(mod$model$depth)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.species.tpi = testdata%>%data.frame(fits)%>%
  group_by(tpi)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for Species richness ----
# Sessile inverts ----
ggmod.species.biog<- ggplot() +
  ylab("")+
  xlab("Sessile invertebrates")+
  geom_point(data=dat.species,aes(x=biog,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.species.biog,aes(x=biog,y=number),alpha=0.5)+
  geom_line(data=predicts.species.biog,aes(x=biog,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.species.biog,aes(x=biog,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Species richness") +
  theme(plot.title = element_text(hjust = 0))
ggmod.species.biog

# depth ----
ggmod.species.depth<- ggplot() +
  ylab("")+
  xlab("Depth")+
  geom_point(data=dat.species,aes(x=depth,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.species.depth,aes(x=depth,y=number),alpha=0.5)+
  geom_line(data=predicts.species.depth,aes(x=depth,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.species.depth,aes(x=depth,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.species.depth

# tpi ----
ggmod.species.tpi<- ggplot() +
  ylab("")+
  xlab("TPI")+
  geom_point(data=dat.species,aes(x=tpi,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.species.tpi,aes(x=tpi,y=number),alpha=0.5)+
  geom_line(data=predicts.species.tpi,aes(x=tpi,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.species.tpi,aes(x=tpi,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.species.tpi

# MODEL Greater than legal size (detrended) ----
dat.legal <- dat %>% filter(scientific=="greater than legal size")

mod=gam(number~s(detrended,k=3,bs='cr'), family=tw,data=dat.legal)

# predict - detrended ----
testdata <- expand.grid(detrended=seq(min(dat$detrended),max(dat$detrended),length.out = 20)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.detrended = testdata%>%data.frame(fits)%>%
  group_by(detrended)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for Greater than legal size ----
# detrended bathymetry ----
ggmod.legal.detrended<- ggplot() +
  ylab("")+
  xlab("Detrended bathymetry")+
  geom_point(data=dat.legal,aes(x=detrended,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.legal.detrended,aes(x=detrended,y=number),alpha=0.5)+
  geom_line(data=predicts.legal.detrended,aes(x=detrended,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.detrended,aes(x=detrended,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Greater than legal size") +
  theme(plot.title = element_text(hjust = 0))
ggmod.legal.detrended

# MODEL Smaller than legal size (tpi) ----
dat.sublegal <- dat %>% filter(scientific=="smaller than legal size")

mod=gam(number~s(tpi,k=3,bs='cr'), family=tw,data=dat.sublegal)

# predict - depth ----
testdata <- expand.grid(tpi=seq(min(dat$tpi),max(dat$tpi),length.out = 20)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.sublegal.tpi = testdata%>%data.frame(fits)%>%
  group_by(tpi)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for Smaller than legal size ----
# tpi ----
ggmod.sublegal.tpi<- ggplot() +
  ylab("")+
  xlab("TPI")+
  geom_point(data=dat.sublegal,aes(x=tpi,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.sublegal.tpi,aes(x=tpi,y=number),alpha=0.5)+
  geom_line(data=predicts.sublegal.tpi,aes(x=tpi,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.sublegal.tpi,aes(x=tpi,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Smaller than legal size") +
  theme(plot.title = element_text(hjust = 0))
ggmod.sublegal.tpi

# Combine with cowplot
library(cowplot)

# view plots
plot.grid.npz6 <- plot_grid(ggmod.total.relief, NULL, NULL,
                       ggmod.species.biog, ggmod.species.depth, ggmod.species.tpi,
                       ggmod.legal.detrended, NULL, NULL,
                       ggmod.sublegal.tpi, NULL, NULL,
                       ncol = 3, labels = c('a','','','b','c','d','e','', '', 'f', '', ''),align = "vh")
plot.grid.npz6

#Save plots
save_plot("plots/fish/abrolhos.npz6.gam.png", plot.grid.npz6,base_height = 9,base_width = 8.5)
