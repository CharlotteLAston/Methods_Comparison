###
# Project: Parks - Abrolhos
# Data:    BOSS & BRUV fish, habitat
# Task:    Plotting fish GAM relationships at npz9
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
name <- "2021-05_Abrolhos_npz9" # for the study

## Set working directory----
working.dir <- getwd()
setwd(working.dir)
#OR Set manually once

# Bring in and format the  data----
dat.maxn <- readRDS("data/Tidy/dat.maxn.rds")%>%
  dplyr::filter(location%in%"NPZ9")%>%
  dplyr::rename(number = maxn)%>%
  glimpse()

dat.length <- readRDS("data/Tidy/dat.length.rds")%>%
  dplyr::filter(location%in%"NPZ9")%>%
  glimpse()

dat <- bind_rows(dat.maxn,dat.length)

# Manually make the most parsimonious GAM models for each taxa ----
#### Abrolhos MaxN ####
unique(dat$scientific)

# MODEL Total abundance (relief + roughness) ----
dat.total <- dat %>% filter(scientific=="total.abundance")

mod=gam(number~s(mean.relief,k=3,bs='cr')+s(roughness,k=3,bs='cr'), family=tw,data=dat.total)

# predict - relief ----
testdata <- expand.grid(mean.relief=seq(min(dat$mean.relief),max(dat$mean.relief),length.out = 20),
                        roughness=mean(mod$model$roughness)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.total.relief = testdata%>%data.frame(fits)%>%
  group_by(mean.relief)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# predict - roughness ----
testdata <- expand.grid(roughness=seq(min(dat$roughness),max(dat$roughness),length.out = 20),
                        mean.relief=mean(mod$model$mean.relief)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.total.roughness = testdata%>%data.frame(fits)%>%
  group_by(roughness)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for Total abundance ----
# relief ----
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

# roughness ----
ggmod.total.roughness<- ggplot() +
  ylab("")+
  xlab("Roughness")+
  geom_point(data=dat.total,aes(x=roughness,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.total.roughness,aes(x=roughness,y=number),alpha=0.5)+
  geom_line(data=predicts.total.roughness,aes(x=roughness,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.total.roughness,aes(x=roughness,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.total.roughness

# MODEL species richness (mean.relief + roughness) ----
dat.species <- dat %>% filter(scientific=="species.richness")

mod=gam(number~s(mean.relief,k=3,bs='cr') + s(roughness,k=3,bs='cr'), family=tw,data=dat.species)

# predict - mean.relief ----
testdata <- expand.grid(mean.relief=seq(min(dat$mean.relief),max(dat$mean.relief),length.out = 20),
                        roughness=mean(mod$model$roughness)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.species.relief = testdata%>%data.frame(fits)%>%
  group_by(mean.relief)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# predict - roughness ----
testdata <- expand.grid(roughness=seq(min(dat$roughness),max(dat$roughness),length.out = 20),
                        mean.relief=mean(mod$model$mean.relief)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.species.roughness = testdata%>%data.frame(fits)%>%
  group_by(roughness)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for species richness ----
# relief ----
ggmod.species.relief<- ggplot() +
  ylab("")+
  xlab("Relief")+
  geom_point(data=dat.species,aes(x=mean.relief,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.species.relief,aes(x=mean.relief,y=number),alpha=0.5)+
  geom_line(data=predicts.species.relief,aes(x=mean.relief,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.species.relief,aes(x=mean.relief,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Species richness") +
  theme(plot.title = element_text(hjust = 0))
ggmod.species.relief

# roughness ----
ggmod.species.roughness<- ggplot() +
  ylab("")+
  xlab("Roughness")+
  geom_point(data=dat.species,aes(x=roughness,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.species.roughness,aes(x=roughness,y=number),alpha=0.5)+
  geom_line(data=predicts.species.roughness,aes(x=roughness,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.species.roughness,aes(x=roughness,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.species.roughness

# MODEL Greater than legal size (roughness) ----
dat.legal <- dat %>% filter(scientific=="greater than legal size")

mod=gam(number~s(roughness,k=3,bs='cr') + s(tpi, k = 3, bs = "cr"), family=tw,data=dat.legal)

# predict - roughness ----
testdata <- expand.grid(roughness=seq(min(dat$roughness),max(dat$roughness),length.out = 20),
                        tpi = mean(mod$model$tpi)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.roughness = testdata%>%data.frame(fits)%>%
  group_by(roughness)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# predict - tpi ----
testdata <- expand.grid(tpi=seq(min(dat$tpi),max(dat$tpi),length.out = 20),
                        roughness = mean(mod$model$roughness)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.legal.tpi = testdata%>%data.frame(fits)%>%
  group_by(tpi)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for Greater than legal size ----
# roughness ----
ggmod.legal.roughness<- ggplot() +
  ylab("")+
  xlab("Roughness")+
  geom_point(data=dat.legal,aes(x=roughness,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.legal.roughness,aes(x=roughness,y=number),alpha=0.5)+
  geom_line(data=predicts.legal.roughness,aes(x=roughness,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.roughness,aes(x=roughness,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Greater than legal size") +
  theme(plot.title = element_text(hjust = 0))
ggmod.legal.roughness

# tpi ----
ggmod.legal.tpi<- ggplot() +
  ylab("")+
  xlab("TPI")+
  geom_point(data=dat.legal,aes(x=tpi,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.legal.tpi,aes(x=tpi,y=number),alpha=0.5)+
  geom_line(data=predicts.legal.tpi,aes(x=tpi,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.legal.tpi,aes(x=tpi,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.legal.tpi

# MODEL Smaller than legal size (depth + roughness) ----
dat.sublegal <- dat %>% filter(scientific=="smaller than legal size")

mod=gam(number~s(depth,k=3,bs='cr')+s(roughness,k=3,bs='cr'), family=tw,data=dat.sublegal)

# predict - depth ----
testdata <- expand.grid(depth=seq(min(dat$depth),max(dat$depth),length.out = 20),
                        roughness=mean(mod$model$roughness)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.sublegal.depth = testdata%>%data.frame(fits)%>%
  group_by(depth)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# predict - roughness ----
testdata <- expand.grid(roughness=seq(min(dat$roughness),max(dat$roughness),length.out = 20),
                        depth=mean(mod$model$depth)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.sublegal.roughness = testdata%>%data.frame(fits)%>%
  group_by(roughness)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

# PLOTS for smaller than legal size ----
# depth ----
ggmod.sublegal.depth<- ggplot() +
  ylab("")+
  xlab("Depth")+
  geom_point(data=dat.sublegal,aes(x=depth,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.sublegal.depth,aes(x=depth,y=number),alpha=0.5)+
  geom_line(data=predicts.sublegal.depth,aes(x=depth,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.sublegal.depth,aes(x=depth,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  ggtitle("Smaller than legal size") +
  theme(plot.title = element_text(hjust = 0))
ggmod.sublegal.depth

# roughness ----
ggmod.sublegal.roughness<- ggplot() +
  ylab("")+
  xlab("Roughness")+
  geom_point(data=dat.sublegal,aes(x=roughness,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.sublegal.roughness,aes(x=roughness,y=number),alpha=0.5)+
  geom_line(data=predicts.sublegal.roughness,aes(x=roughness,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.sublegal.roughness,aes(x=roughness,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.sublegal.roughness

# Combine with cowplot
library(cowplot)

# view plots
plot.grid.npz9 <- plot_grid(ggmod.total.relief, ggmod.total.roughness,
                            ggmod.species.relief,ggmod.species.roughness,
                            ggmod.legal.roughness, ggmod.legal.tpi,
                            ggmod.sublegal.depth, ggmod.sublegal.roughness,
                            ncol = 2, labels = c('a','b','c','d','e','f','g','h'),align = "vh")
plot.grid.npz9

#Save plots
save_plot("plots/fish/abrolhos.npz9.gam.png", plot.grid.npz9,base_height = 9,base_width = 8.5)
