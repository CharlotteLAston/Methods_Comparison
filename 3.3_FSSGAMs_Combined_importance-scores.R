###
# Project: Methods Comparison
# Data:    BOSS & BRUV fish Abrolhos
# Task:    Plotting fish importance scores
# Author:  Charlotte (Claude)
# Date:    Jan 2024
##

rm(list=ls())

# libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(patchwork)
library(forcats)
library(rcartocolor)

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

##### READ IN FORMATTED DATA ####
setwd(out_dir)
#read in data 
dat1 <- read.csv("2021-05_Abrolhos_BOSS-BRUV_greater_less_by_species_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  mutate(location = "Abrolhos") %>% 
  glimpse()

dat2 <- read.csv("Ningaloo_PtCloates_BOSS-BRUV_greater_less_by_species_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  mutate(location = "Ningaloo") %>% 
  glimpse()

greater_less <- bind_rows(dat1,dat2)%>%
  glimpse()

dat.greater.less <- greater_less %>%
  mutate(label=NA)%>%
  mutate(resp.var=factor(resp.var, levels = c("> Length Maturity_Lethrinidae Lethrinus miniatus", 
                                              "< Length Maturity_Lethrinidae Lethrinus miniatus"
                                              )))%>%

  # greater or less than length maturity by species

  #Abrolhos
  mutate(label=ifelse(predictor=="macroalgae"& location=="Abrolhos" & resp.var=="< Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="sd.relief"& location=="Abrolhos" & resp.var=="< Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="method"& location=="Abrolhos" & resp.var=="< Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  
  mutate(label=ifelse(predictor=="biog"& location=="Abrolhos" & resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"& location=="Abrolhos" & resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="method"& location=="Abrolhos" & resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  
  # Ningaloo
  mutate(label=ifelse(predictor=="sdrel"& location=="Ningaloo" & resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"& location=="Ningaloo" & resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="method"& location=="Ningaloo" & resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  filter(!is.na(resp.var)) %>% 
 
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

imp.species.greater.less.ningaloo <- ggplot(dat.greater.less%>%dplyr::filter(location%in%c("Ningaloo")) %>% 
                                            mutate(predictor = fct_relevel(predictor,c("biog", "depth", "detrended", "relief", "roughness", "sand","sdrel", "tpi","method"))), 
                                   aes(x=predictor,y=resp.var,fill=importance)) +
  geom_tile(show.legend=F) +
  scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
                       limits = c(0, 1))+
  scale_y_discrete(labels=c("> Length Maturity\nLethrinus miniatus"))+
  scale_x_discrete(labels = c("% Biogenic\nReef","Depth","Detrended\nBathymetry", "Mean\nRelief", "Roughness", "SD Relief", "% Sand", "TPI","Method"))+
  labs(x = NULL, y = NULL, title = "Ningaloo") +
  theme_classic()+
  Theme1+
  geom_text(aes(label=label)) +
  geom_vline(xintercept=8.5, colour="white", linewidth=3) +
  theme(plot.title = element_text(hjust = -0.2, size=8)) # Looks crap here but title comes back in exported version
imp.species.greater.less.ningaloo

imp.species.greater.less.abrolhos <- ggplot(dat.greater.less%>%dplyr::filter(location %in% c("Abrolhos")) %>% 
                                            mutate(predictor = fct_relevel(predictor,c("biog", "depth", "detrended","macroalgae", "mean.relief","roughness", "sd.relief", "tpi","method"))), 
                                   aes(x=predictor,y=resp.var,fill=importance)) +
  geom_tile(show.legend=F) +
  scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
                       limits = c(0, 1))+
  scale_y_discrete(labels=c("> Length Maturity\nLethrinus miniatus","< Length Maturity\nLethrinus miniatus"))+
  scale_x_discrete(labels = c("% Biogenic\nReef", "Depth", "Detrended\nbathymetry","% Macroalgae", "Mean Relief","Roughness","SD Relief", "TPI","Method"))+
  labs(x = NULL, y = NULL, title = "Abrolhos") +
  theme_classic()+
  Theme1+
  geom_text(aes(label=label)) +
  geom_vline(xintercept=8.5, colour="white", linewidth=3) +
  theme(plot.title = element_text(hjust = -0.2, size=8)) # Looks crap here but title comes back in exported version
imp.species.greater.less.abrolhos

gg.importance.greater.less <- imp.species.greater.less.ningaloo / imp.species.greater.less.abrolhos

#save output - changed dimensions for larger text in report
setwd(fig_dir)
save_plot("combined.importance.redthroat.png", gg.importance.greater.less, base_height = 5, base_width = 6.275)

