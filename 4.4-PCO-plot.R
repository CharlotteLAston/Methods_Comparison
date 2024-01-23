###############################
# Project: Methods Comparison
# Data:    BOSS & BRUV fish
# Task:    Plotting PCO of the 
#          communities observed 
#          by the BOSS/BRUV
# Author:  Charlotte (Sam)
# Date:    September 2023
###############################

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

#### Import data for PCO of method (using presence/absence data) ####
setwd(data_dir)

pco.abrolhos <- read.csv("Abrolhos_pco.csv")                 #import CAP.xlsx data sheet (manually exported from primer)
vectors.abrolhos <- read.csv("Abrolhos_vectors.csv") %>%  #import CAP.xlsx vectors sheet (manually exported from primer)
  filter(r>=0.5|taxa %in% "Lethrinidae Lethrinus miniatus") %>% 
  mutate(taxa.short = str_extract(taxa, "(?<=\\s)[:alpha:]+\\s[:alpha:]+")) %>%
  mutate(genus = str_extract(taxa.short, "[:alpha:]+(?=\\s)"),
         species = str_extract(taxa.short, "(?<=\\s)[:alpha:]+")) %>% 
  mutate(taxa.short = str_replace(taxa.short, "([:alpha:]{1})(.*)\\s", "\\1\\. ")) %>% 
  mutate(taxa.short = ifelse(species %in% "spp", str_replace(taxa.short,"([:alpha:]{1})(.)", str_extract(genus, "[:alpha:]{1,}")), taxa.short))

pco.ningaloo <- read.csv("ningaloo_PCO_values.csv")              #import CAP.xlsx data sheet (manually exported from primer)
vectors.ningaloo <- read.csv("ningaloo_vector_values.csv") %>%  #import CAP.xlsx vectors sheet (manually exported from primer)
  filter(r>=0.5) %>% 
  mutate(taxa.short = str_extract(taxa, "(?<=\\s)[:alpha:]+\\s[:alnum:]+")) %>%
  mutate(genus = str_extract(taxa.short, "[:alpha:]+(?=\\s)"),
         species = str_extract(taxa.short, "(?<=\\s)[:alnum:]+")) %>% 
  mutate(taxa.short = str_replace(taxa.short, "([:alpha:]{1})(.*)\\s", "\\1\\. "))%>% 
  mutate(taxa.short = ifelse(species %in% "sp1", str_replace(taxa.short,"([:alpha:]{1})(.)", str_extract(genus, "[:alpha:]{1,}")), taxa.short))

#04.3.2 : Set theme and palette for PCO plot ----

Palette<- c("#117733", "#88CCEE")

#04.3.3 : Import fish pics for PCO plot ----
# pic.p.n <- readJPEG("G:/My Drive/meg_graphics/Parapercis nebulosa 5cmL.jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.c.c <- readJPEG("G:/My Drive/meg_graphics/Carangoides chrysophrys.jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.s.f <- readJPEG("G:/My Drive/meg_graphics/Sufflamen fraenatum 5cmL.jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.x.l <- readJPEG("G:/My Drive/meg_graphics/Xanthichthys lineopunctatus (online - fishbase).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.e.b <- readJPEG("G:/My Drive/meg_graphics/Elagatis bipinnulata (online - fishbase).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.p.f <- readJPEG("G:/My Drive/meg_graphics/Pristipomoides filamentosus (online fishbase).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.e <- readJPEG("G:/My Drive/meg_graphics/Echeneidae (fishes of australia).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.l.r <- readJPEG("G:/My Drive/meg_graphics/Lethrinus rubrioperculatus (fishbase online).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.g.e <- readJPEG("G:/My Drive/meg_graphics/Gymnocranius_euanus(online).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.p.d <- readJPEG("G:/My Drive/meg_graphics/Pterocaesio digramma (online).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.n.m <- readJPEG("G:/My Drive/meg_graphics/naso tuberosus (online).jpg")%>%
#   rasterGrob(interpolate=TRUE)
# pic.s.l <- readJPEG("G:/My Drive/meg_graphics/Seriola lalandi (online).jpg")%>%
#   rasterGrob(interpolate=TRUE)

#### Check labels/vectors exist in multiple sites ####
# 
# vectors.abrolhos <- vectors.abrolhos %>%
#   mutate(vectors = case_when(
#       vectors == "yes" & (BRUV + eDNA + `Drop-Camera`) < 2 ~ "NA",
#       vectors == NA & (BRUV + eDNA + `Drop-Camera`) > 2 & r >= 0.4 ~ "yes",
#       TRUE ~ vectors  # Keep the original value if none of the conditions match
#     )
#   )

#### Set custom label positioning for PCO plot ####
# labels <- vectors %>% filter(vectors =="yes")
vectors.abrolhos <- vectors.abrolhos %>% 
  mutate(PCO1.lab = ifelse(taxa.short %in% "C. auratus", PCO1-0.23, 
                          ifelse(taxa.short %in% "C. rubescens", PCO1, 
                                 ifelse(taxa.short %in% "P. spilurus", PCO1-0.35,
                                        ifelse(taxa.short %in% "S. cyanolaemus", PCO1+0.63,
                                               ifelse(taxa.short %in% "Chromis spp", PCO1+0.32, 
                                                      ifelse(taxa.short %in% "L. miniatus", PCO1+0, NA))))))) %>% 
  mutate(PCO2.lab = ifelse(taxa.short %in% "C. auratus", PCO2+0.037, 
                           ifelse(taxa.short %in% "C. rubescens", PCO2+0.11, 
                                  ifelse(taxa.short %in% "P. spilurus", PCO2-0.03,
                                         ifelse(taxa.short %in% "S. cyanolaemus", PCO2+0.04,
                                                ifelse(taxa.short %in% "Chromis spp", PCO2+0.01, 
                                                       ifelse(taxa.short %in% "L. miniatus", PCO2+0.03, NA)))))))
#### Plot PCO ####
gg.cap.abrolhos <- ggplot()+
  # annotation_custom(pic.p.n, xmin=-0.65, xmax=0.65, ymin=0.5, ymax=0.56)+
  # annotation_custom(pic.c.c, xmin=0.05, xmax=0.175, ymin=0.55, ymax=0.675)+
  # annotation_custom(pic.s.f, xmin=0.345, xmax=0.545, ymin=-0.175, ymax=-0.3)+
  # annotation_custom(pic.p.f, xmin=-0.325, xmax=-0.125, ymin=-0.125, ymax=-0.325)+
  # annotation_custom(pic.e, xmin=-0.2, xmax=0.06, ymin=0.345, ymax=0.445)+
  # annotation_custom(pic.l.r, xmin=0.35, xmax=0.55, ymin=0.625, ymax=0.75)+
  # annotation_custom(pic.g.e, xmin=0.3, xmax=0.5, ymin=0.325, ymax=0.425)+
  # annotation_custom(pic.g.e, xmin=0.3, xmax=0.5, ymin=0.325, ymax=0.425)+
  # annotation_custom(pic.p.d, xmin=0.385, xmax=0.535, ymin=-0.395, ymax=-0.475)+
  # annotation_custom(pic.n.m, xmin=0.05, xmax=0.205, ymin=-0.425, ymax=-0.555)+
  # annotation_custom(pic.s.l, xmin=-0.15, xmax=0.05, ymin=-0.325, ymax=-0.455)+
  #annotation_custom(pic.x.l, xmin=-0.275, xmax=-0.455, ymin=0.115, ymax=0.3)+
  #annotation_custom(pic.e.b, xmin=-0.325, xmax=-0.15, ymin=-0.175, ymax=-0.325)+
  geom_point(data=pco.abrolhos, aes(x=PCO1, y=PCO2, fill=method),alpha=1, size=5,shape=21)+#fill
  scale_fill_manual(values=c("#117733", "#88CCEE"))+
  theme_classic() +
  Theme1+
  # guides(fill=guide_legend(override.aes = list(size = 6)))+
  # guides(shape=guide_legend(override.aes = list(size = 6)))+
  #scale_x_continuous(limits = c(-1.5, 1.5)) +
  #scale_y_continuous(limits = c(-2, 2)) +
  theme(legend.position = c(0.925, 0.896)) +
  xlab ("PCO1 (9.1%)")+                
  ylab ("PCO2 (36.5%)") +
  geom_segment(data=vectors.abrolhos,aes(x=0,y=0,xend=(PCO1),yend=(PCO2)),size = 0.5,colour="darkgrey",arrow = arrow(angle=25,length=unit(0.25,"cm"))) +
  geom_text(data= vectors.abrolhos ,aes(x=(PCO1.lab),y=((PCO2.lab)),label = taxa.short),size = 4,colour="black",fontface='italic',angle=0)+
  ggplot2::annotate("text", x=-1.65, y=1, label="(b) Abrolhos", size = 4, fontface=1)
gg.cap.abrolhos
#ggsave(paste(wd, "04-out", paste0("04-PCO-method_", primer, ".png"), sep="/"), gg.cap,width = 20, height = 20,units = "cm")

#### CAP Plot for Ningaloo ####
# labels <- vectors %>% filter(vectors =="yes")
vectors.ningaloo <- vectors.ningaloo %>% 
  mutate(PCO1.lab = ifelse(taxa.short %in% "Gymnocranius sp1", PCO1, 
                           ifelse(taxa.short %in% "L. miniatus", PCO1-0.175, 
                                  ifelse(taxa.short %in% "L. rubrioperculatus", PCO1-0.22,
                                         ifelse(taxa.short %in% "P. nebulosa", PCO1, NA))))) %>% 
  mutate(PCO2.lab = ifelse(taxa.short %in% "Gymnocranius sp1", PCO2+0.03, 
                           ifelse(taxa.short %in% "L. miniatus", PCO2-0.031, 
                                  ifelse(taxa.short %in% "L. rubrioperculatus", PCO2,
                                         ifelse(taxa.short %in% "P. nebulosa", PCO2-0.025, NA))))) %>% 
  mutate(taxa.short = ifelse(taxa.short %in% "Gymnocranius sp1", "Gymnocranius spp.", taxa.short))

gg.cap.ningaloo <- ggplot()+
  # annotation_custom(pic.p.n, xmin=-0.65, xmax=0.65, ymin=0.5, ymax=0.56)+
  # annotation_custom(pic.c.c, xmin=0.05, xmax=0.175, ymin=0.55, ymax=0.675)+
  # annotation_custom(pic.s.f, xmin=0.345, xmax=0.545, ymin=-0.175, ymax=-0.3)+
  # annotation_custom(pic.p.f, xmin=-0.325, xmax=-0.125, ymin=-0.125, ymax=-0.325)+
  # annotation_custom(pic.e, xmin=-0.2, xmax=0.06, ymin=0.345, ymax=0.445)+
  # annotation_custom(pic.l.r, xmin=0.35, xmax=0.55, ymin=0.625, ymax=0.75)+
  # annotation_custom(pic.g.e, xmin=0.3, xmax=0.5, ymin=0.325, ymax=0.425)+
  # annotation_custom(pic.g.e, xmin=0.3, xmax=0.5, ymin=0.325, ymax=0.425)+
  # annotation_custom(pic.p.d, xmin=0.385, xmax=0.535, ymin=-0.395, ymax=-0.475)+
  # annotation_custom(pic.n.m, xmin=0.05, xmax=0.205, ymin=-0.425, ymax=-0.555)+
  # annotation_custom(pic.s.l, xmin=-0.15, xmax=0.05, ymin=-0.325, ymax=-0.455)+
#annotation_custom(pic.x.l, xmin=-0.275, xmax=-0.455, ymin=0.115, ymax=0.3)+
#annotation_custom(pic.e.b, xmin=-0.325, xmax=-0.15, ymin=-0.175, ymax=-0.325)+
geom_point(data=pco.ningaloo, aes(x=PCO1, y=PCO2, fill=method),alpha=1, size=5,shape=21)+#fill
  scale_fill_manual(values=c("#117733", "#88CCEE"))+
  theme_classic() +
  Theme1+
  # guides(fill=guide_legend(override.aes = list(size = 6)))+
  # guides(shape=guide_legend(override.aes = list(size = 6)))+
  scale_x_continuous(limits = c(-1, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.7)) +
  theme(legend.position = "none") +
  xlab ("PCO1 (24.5%)")+                
  ylab ("PCO2 (6.1%)") +
  geom_segment(data=vectors.ningaloo,aes(x=0,y=0,xend=(PCO1),yend=(PCO2)),size = 0.5,colour="darkgrey",arrow = arrow(angle=25,length=unit(0.25,"cm")))+
  geom_text(data= vectors.ningaloo, aes(x=(PCO1.lab),y=((PCO2.lab)),label = taxa.short),size = 4,colour="black",fontface='italic',angle=0)+
  ggplot2::annotate("text", x=-0.95, y=0.7, label="(a) Ningaloo", size = 4, fontface=1)
gg.cap.ningaloo

## Put the plots together 
setwd(fig_dir)
#y.label <- textGrob(expression(paste("Temperature (",degree~C,")")), gp=gpar(fontsize=13), rot=90)
#x.label <- textGrob("Method", gp=gpar(fontsize=13))
legend <- gtable_filter(ggplotGrob(gg.cap.ningaloo), "guide-box")

PCO.Plots <-grid.arrange(arrangeGrob(gg.cap.ningaloo + theme(legend.position="none"),
                                     gg.cap.abrolhos))

ggsave(PCO.Plots , filename=paste0("PCO_plots_both_areas.png"),height = a4.width*1.35, width = a4.width*1.2, units  ="mm", dpi = 300 )

