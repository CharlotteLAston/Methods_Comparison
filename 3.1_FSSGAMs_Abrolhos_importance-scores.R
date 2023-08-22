###
# Project: Methods Comparison
# Data:    BOSS & BRUV fish Abrolhos
# Task:    Plotting fish importance scores
# Author:  Charlotte (Claude)
# Date:    July 2023
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
dat1 <- read.csv("2021-05_Abrolhos_BOSS-BRUV_greater_less_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

dat2 <- read.csv("2021-05_Abrolhos_BOSS-BRUV_greater_less_by_species_all.var.imp.csv")%>% #from local copy
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  glimpse()

greater_less <- bind_rows(dat1,dat2)%>%
  glimpse()

dat.greater.less <- greater_less %>%
  mutate(label=NA)%>%
  mutate(resp.var=factor(resp.var, levels = c("< Length Maturity","> Length Maturity",	
                                              "> Length Maturity_Labridae Choerodon rubescens",
                                              "> Length Maturity_Lethrinidae Lethrinus miniatus", 
                                              "< Length Maturity_Lethrinidae Lethrinus miniatus",
                                              "> Length Maturity_Sparidae Chrysophrys auratus",
                                              "< Length Maturity_Sparidae Chrysophrys auratus")))%>%
  # greater or less than length maturity
  mutate(label=ifelse(predictor=="method"&resp.var=="> Length Maturity","X",label))%>%
  mutate(label=ifelse(predictor=="biog"&resp.var=="> Length Maturity","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"&resp.var=="> Length Maturity","X",label))%>%
  mutate(label=ifelse(predictor=="method"&resp.var=="< Length Maturity","X",label))%>%
  mutate(label=ifelse(predictor=="macroalgae"&resp.var=="< Length Maturity","X",label))%>%
  # greater or less than length maturity by species
  mutate(label=ifelse(predictor=="method"&resp.var=="> Length Maturity_Labridae Choerodon rubescens","X",label))%>%
  
  mutate(label=ifelse(predictor=="macroalgae"&resp.var=="< Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"&resp.var=="< Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="method"&resp.var=="< Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  
  mutate(label=ifelse(predictor=="biog"&resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="detrended"&resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  mutate(label=ifelse(predictor=="method"&resp.var=="> Length Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
  
  mutate(label=ifelse(predictor=="macroalgae"&resp.var=="< Length Maturity_Sparidae Chrysophrys auratus","X",label))%>%
  mutate(label=ifelse(predictor=="method"&resp.var=="< Length Maturity_Sparidae Chrysophrys auratus","X",label))%>%
  
  mutate(label=ifelse(predictor=="method"&resp.var=="> Length Maturity_Sparidae Chrysophrys auratus","X",label))%>%
  mutate(label=ifelse(predictor=="biog"&resp.var=="> Length Maturity_Sparidae Chrysophrys auratus","X",label))%>%
  mutate(label=ifelse(predictor=="roughness"&resp.var=="> Length Maturity_Sparidae Chrysophrys auratus","X",label))%>%
  mutate(predictor = fct_relevel(predictor,c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended", "method"))) %>% 
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

imp.full.greater.less<- ggplot(dat.greater.less%>%dplyr::filter(resp.var%in%c("< Length Maturity", "> Length Maturity")), 
                   aes(x=predictor,y=resp.var,fill=importance)) +
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
                       limits = c(0, 1))+
  scale_y_discrete(labels=c("< Length at maturity","> Length at maturity"))+
  labs(x = NULL, y = NULL, title = "All indicator species") +
  theme_classic()+
  Theme1+
  geom_text(aes(label=label)) +
  geom_vline(xintercept=7.5, colour="white", linewidth=3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = -0.2, vjust = -15, size=8)) # Looks crap here but title comes back in exported version
imp.full.greater.less

imp.species.greater.less <- ggplot(dat.greater.less%>%dplyr::filter(resp.var%in%c("> Length Maturity_Labridae Choerodon rubescens",
                                                                               "> Length Maturity_Lethrinidae Lethrinus miniatus",
                                                                               "< Length Maturity_Lethrinidae Lethrinus miniatus",
                                                                               "> Length Maturity_Sparidae Chrysophrys auratus",
                                                                               "< Length Maturity_Sparidae Chrysophrys auratus")), 
                   aes(x=predictor,y=resp.var,fill=importance)) +
  geom_tile(show.legend=F) +
  scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
                       limits = c(0, 1))+
  scale_y_discrete(labels=c("> Length Maturity\nChoerodon rubescens","> Length Maturity\nLethrinus miniatus","< Length Maturity\nLethrinus miniatus",
                            "< Length Maturity\nChrysophrys auratus", "> Length Maturity\nChrysophrys auratus"))+
  scale_x_discrete(labels = c("Depth", "Macroalgae", "Biogenic\nReef", "Mean\nRelief","TPI", "Roughness", "Detrended\nBathymetry", "Method"))+
  labs(x = NULL, y = NULL, title = "By indicator species") +
  theme_classic()+
  Theme1+
  geom_text(aes(label=label)) +
  geom_vline(xintercept=7.5, colour="white", linewidth=3) +
  theme(plot.title = element_text(hjust = -0.2, size=8)) # Looks crap here but title comes back in exported version
imp.species.greater.less

gg.importance.greater.less <- imp.full.greater.less / imp.species.greater.less


#save output - changed dimensions for larger text in report
setwd(fig_dir)
save_plot("abrolhos.importance.greater.less.png", gg.importance.greater.less, base_height = 5, base_width = 6.275)

#read in data - negative values manually added
# setwd(out_dir)
# dat3 <- read.csv("2021-05_Abrolhos_BOSS-BRUV_mat_groups_all.var.imp.csv")%>% #from local copy
#   rename(resp.var=X)%>%
#   gather(key=predictor,value=importance,2:ncol(.))%>%
#   glimpse()
# 
# dat4 <- read.csv("2021-05_Abrolhos_BOSS-BRUV_mat_groups_by_species_all.var.imp.csv")%>% #from local copy
#   rename(resp.var=X)%>%
#   gather(key=predictor,value=importance,2:ncol(.))%>%
#   glimpse()
# 
# dat.groups <- bind_rows(dat3,dat4)%>%
#   glimpse()
# 
# dat.groups.all <- dat.groups %>%
#   mutate(label=NA)%>%
#   mutate(resp.var=factor(resp.var, levels = c("> 1.25x Maturity", "> Length Maturity but < 1.25x Maturity",
#                                               ">50 Maturity but < Maturity",
#                                               ">50 Maturity but < Maturity_Labridae Choerodon rubescens",
#                                               "> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus",
#                                               ">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus",
#                                               "> Length Maturity but < 1.25x Maturity_Lutjanidae Pristipomoides multidens",
#                                               ">50 Maturity but < Maturity_Sparidae Chrysophrys auratus")))%>%
#   mutate(label=ifelse(predictor=="macroalgae"&resp.var==">50 Maturity but < Maturity","X",label))%>%
#   mutate(label=ifelse(predictor=="biog"&resp.var==">50 Maturity but < Maturity","X",label))%>%
#   mutate(label=ifelse(predictor=="method"&resp.var==">50 Maturity but < Maturity","X",label))%>%
#   mutate(label=ifelse(predictor=="method"&resp.var=="> Length Maturity but < 1.25x Maturity","X",label))%>%
#   mutate(label=ifelse(predictor=="tpi"&resp.var=="> 1.25x Maturity","X",label))%>%
#   mutate(label=ifelse(predictor=="macroalgae"&resp.var=="> 1.25x Maturity","X",label))%>%
#   mutate(label=ifelse(predictor=="depth"&resp.var==">50 Maturity but < Maturity_Labridae Choerodon rubescens","X",label))%>%
#   mutate(label=ifelse(predictor=="tpi"&resp.var==">50 Maturity but < Maturity_Labridae Choerodon rubescens","X",label))%>%
#   mutate(label=ifelse(predictor=="method"&resp.var==">50 Maturity but < Maturity_Labridae Choerodon rubescens","X",label))%>%
#   mutate(label=ifelse(predictor=="method"&resp.var=="> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
#   mutate(label=ifelse(predictor=="detrended"&resp.var=="> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
#   mutate(label=ifelse(predictor=="biog"&resp.var=="> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
#   mutate(label=ifelse(predictor=="macroalgae"&resp.var==">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
#   mutate(label=ifelse(predictor=="detrended"&resp.var==">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
#   mutate(label=ifelse(predictor=="method"&resp.var==">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus","X",label))%>%
#   mutate(label=ifelse(predictor=="biog"&resp.var=="> Length Maturity but < 1.25x Maturity_Lutjanidae Pristipomoides multidens","X",label))%>%
#   mutate(label=ifelse(predictor=="depth"&resp.var=="> Length Maturity but < 1.25x Maturity_Lutjanidae Pristipomoides multidens","X",label))%>%
#   mutate(label=ifelse(predictor=="roughness"&resp.var=="> Length Maturity but < 1.25x Maturity_Lutjanidae Pristipomoides multidens","X",label))%>%
#   mutate(label=ifelse(predictor=="macroalgae"&resp.var==">50 Maturity but < Maturity_Sparidae Chrysophrys auratus","X",label))%>%
#   mutate(label=ifelse(predictor=="method"&resp.var==">50 Maturity but < Maturity_Sparidae Chrysophrys auratus","X",label))%>%
#   mutate(label=ifelse(predictor=="biog"&resp.var==">50 Maturity but < Maturity_Sparidae Chrysophrys auratus","X",label))%>%
#   mutate(predictor = fct_relevel(predictor,c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended", "method"))) %>% 
#   glimpse()
# 
# # Plot gg.importance.scores ----
# 
# imp.all.groups <- ggplot(dat.groups.all%>%dplyr::filter(resp.var%in%c("> 1.25x Maturity", "> Length Maturity but < 1.25x Maturity",
#                                                                      ">50 Maturity but < Maturity")), 
#                         aes(x=predictor,y=resp.var,fill=importance)) +
#   geom_tile(show.legend=T) +
#   scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
#                        limits = c(0, 1))+
#   scale_y_discrete(labels=c("> 1.25x Maturity", "> Length Maturity but < 1.25x Maturity",
#                             ">50 Maturity but < Maturity"))+
#   labs(x = NULL, y = NULL, title = "All indicator species") +
#   theme_classic()+
#   Theme1+
#   geom_text(aes(label=label)) +
#   geom_vline(xintercept=7.5, colour="white", linewidth=3) +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         axis.line.x = element_blank(),
#         plot.title = element_text(hjust = -0.4, vjust = -15, size=8)) # Looks crap here but title comes back in exported version
# imp.all.groups
# 
# imp.species.groups <- ggplot(dat.groups.all %>%
#                                dplyr::filter(resp.var%in%c(">50 Maturity but < Maturity_Labridae Choerodon rubescens",
#                                                                           "> Length Maturity but < 1.25x Maturity_Lethrinidae Lethrinus miniatus",
#                                                                           ">50 Maturity but < Maturity_Lethrinidae Lethrinus miniatus",
#                                                                           "> Length Maturity but < 1.25x Maturity_Lutjanidae Pristipomoides multidens",
#                                                                           ">50 Maturity but < Maturity_Sparidae Chrysophrys auratus")), 
#                         aes(x=predictor,y=resp.var,fill=importance)) +
#   geom_tile(show.legend=F) +
#   scale_fill_gradientn(legend_title, colours=c(re), na.value = "grey98",
#                        limits = c(0, 1))+
#   scale_y_discrete(labels=c(">50 Maturity but < Maturity\nChoerodon rubescens",
#                             "> Length Maturity but < 1.25x Maturity\nLethrinus miniatus",
#                             ">50 Maturity but < Maturity\nLethrinus miniatus",
#                             "> Length Maturity but < 1.25x Maturity\nPristipomoides multidens",
#                             ">50 Maturity but < Maturity\nChrysophrys auratus"))+
#   scale_x_discrete(labels = c("Depth", "Macroalgae", "Biogenic\nReef", "Mean\nRelief","TPI", "Roughness", "Detrended\nBathymetry", "Method"))+
#   labs(x = NULL, y = NULL, title = "By indicator species") +
#   theme_classic()+
#   geom_vline(xintercept=7.5, colour="white", linewidth=3) +
#   Theme1+
#   geom_text(aes(label=label)) +
#   theme(plot.title = element_text(hjust = -0.4, size=8)) # Looks crap here but title comes back in exported version
# imp.species.groups
# 
# gg.importance.groups <- imp.all.groups / imp.species.groups
# gg.importance.groups
# 
# #save output - changed dimensions for larger text in report
# setwd(fig_dir)
# save_plot("abrolhos.importance.groups.png", gg.importance.groups,base_height = 5,base_width = 6.275)





