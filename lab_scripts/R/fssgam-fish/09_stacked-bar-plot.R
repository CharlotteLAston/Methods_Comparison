###
# Project: Parks - Abrolhos
# Data:    BOSS fish, habitat
# Task:    Plotting 10 most abundant species w/ cute pics
# author:  Claude
# date:    Nov-Dec 2021
##

# Set directories----
rm(list=ls())

# Libraries required
library(GlobalArchive)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(rgdal)
library(raster)
library(png)
library(cowplot)

## Set your working directory ----
working.dir <- getwd()
setwd(working.dir)
#OR set manually once
theme_collapse<-theme(      
  panel.grid.major=element_line(colour = "white"), 
  panel.grid.minor=element_line(colour = "white", size = 0.25), 
  plot.margin= grid::unit(c(0, 0, 0, 0), "in"))

theme.larger.text<-theme(
  strip.text.x = element_text(size = 5,angle = 0),
  strip.text.y = element_text(size = 5),
  axis.title.x=element_text(vjust=-0.0, size=10),
  axis.title.y=element_text(vjust=0.0,size=10),
  axis.text.x=element_text(size=8),
  axis.text.y=element_text(size=8),
  legend.title = element_text(family="TN",size=8),
  legend.text = element_text(family="TN",size=8))

# Load fish pictures for plotting ----
#anampses geographicus
a.g <- readPNG("data/images/Anampses_geographicus_nb_TAYLOR.png")
a.g <- as.raster(a.g)

#chaetodon assarius
c.ass <- readPNG("data/images/Chaetodon assarius-3cmL.png")
c.ass <- as.raster(c.ass)

#choerodon rubescens
c.r <- readPNG("data/images/Choerodon rubescens 3cm.png")
c.r <- as.raster(c.r)

#chromis westaustralis
c.w <- readPNG("data/images/Chromis westaustralis-3cmL.png")
c.w <- as.raster(c.w)

#coris auricularis
c.a <- readPNG("data/images/Coris auricularis-3cmL.png")
c.a <- as.raster(c.a)

#suzeicthys cyanolaemus
s.c <- readPNG("data/images/Labridae-Dark.png")
s.c <- as.raster(s.c)

#lethrinus miniatus
l.m <- readPNG("data/images/Lethrinus miniatus 3cm.png")
l.m <- as.raster(l.m)

#neatypus obliquus
n.o <- readPNG("data/images/Neatypus obliquus-3cmL.png")
n.o <- as.raster(n.o)

#parupeneus spilurus
p.s <- readPNG("data/images/Parupeneus_spilurus_nb_TAYLOR.png")
p.s <- as.raster(p.s)

#chrysophrys auratus
c.au <- readPNG("data/images/Chrysophrys auratus 3cm.png")
c.au <- as.raster(c.au)

#pentapous nagasakiensis
p.n <- readPNG("data/images/Pentapodus porosus-3cmL.png")
p.n <- as.raster(p.n)

#lethrinus nebulosus
l.n <- readPNG("data/images/lethrinus nebulosus 3cm.png")
l.n <- as.raster(l.n)

#seriola dumerili
s.d <- readPNG("data/images/seriola_dumerili_nb.png")
s.d <- as.raster(s.d)

#pristipomoides multidens
p.m <- readPNG("data/images/Pristipomoides multidens 3cm.png")
p.m <- as.raster(p.m)

#pseudocaranx spp
p.spp <- readPNG("data/images/Pseudocaranx dentex-3cm.png")
p.spp <- as.raster(p.spp)

#gymnothorax woodwardi
#aint no pic for this one

# #pseudanthias
# p.spp <- readPNG("data/images/Pseudanthias rubrizonatus.png")
# p.spp <- as.raster(p.spp)

# read in maxn
maxn.boss <- read.csv("data/Tidy/2021-05_Abrolhos_BOSS.complete.maxn.csv")%>%
  glimpse()

maxn.bruv <- read.csv("data/Tidy/2021-05_Abrolhos_stereo-BRUVs.complete.maxn.csv")%>%
  glimpse()

maxn <- bind_rows(maxn.boss,maxn.bruv)

#NPZ6

# workout total maxn for each species ---
maxn.npz6.10<-maxn%>%
  dplyr::filter(location%in%"NPZ6")%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  ungroup()%>%
  top_n(11)%>%
  dplyr::filter(!scientific %in% c('Unknown spp', 'SUS sus'))%>%
  glimpse()

test1<-maxn%>%
  dplyr::filter(location%in%"NPZ6"&status%in%"Fished")%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  ungroup()%>%
  top_n(11)%>%
  dplyr::filter(!scientific %in% c('Unknown spp', 'SUS sus'))%>%
  glimpse()

test2<-maxn%>%
  dplyr::filter(location%in%"NPZ6"&status%in%"No-take")%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  ungroup()%>%
  top_n(11)%>%
  dplyr::filter(!scientific %in% c('Unknown spp', 'SUS sus'))%>%
  glimpse()

## Total frequency of occurance 
# I think we could remove this section - but maybe good to see sometimes
bar.npz6<-ggplot(maxn.npz6.10, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",position=position_dodge())+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  #scale_x_discrete(limits = rev(levels(scientific)))+
  #annotation_custom(lcpic, xmin=0.5, xmax=2.5, ymin=.75, ymax=1.5)+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_collapse
bar.npz6

## Top ten plot ----
bar.npz6.top.10<-ggplot(maxn.npz6.10%>%mutate(scientific=str_replace_all(.$scientific,
  c("miniatus"="miniatus*","auratus"="auratus*","rubescens"="rubescens*","nebulosus"="nebulosus*"))), aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",colour="black",fill="lightgrey",position=position_dodge())+
  ylim (0, 1180)+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"))+
  theme_collapse+
  theme.larger.text+
  annotation_raster(c.w, xmin=9.85,xmax=10.15,ymin=1100, ymax=1180)+
  annotation_raster(c.a, xmin=8.7,xmax=9.3,ymin=805, ymax=1005)+
  annotation_raster(l.m, xmin=7.5, xmax=8.5, ymin=120, ymax=400)+
  annotation_raster(c.au, xmin=6.5,xmax=7.5,ymin=110, ymax=405)+
  annotation_raster(p.s, xmin=5.7,xmax=6.3,ymin=80, ymax=265)+
  annotation_raster(c.r, xmin=4.6,xmax=5.4,ymin=75, ymax=330)+
  annotation_raster(s.c, xmin=3.8,xmax=4.2,ymin=65, ymax=185)+
  annotation_raster(n.o, xmin=2.75,xmax=3.25,ymin=65, ymax=190)+
  annotation_raster(p.n, xmin=1.72,xmax=2.25,ymin=55, ymax=225)+
  annotation_raster(l.n, xmin=0.5,xmax=1.5,ymin=55, ymax=330)
# bar.npz6.top.10

#save out plot
ggsave("plots/fish/stacked.bar.plot.npz6.png",bar.npz6.top.10,dpi=600,width=6.0, height = 6)

#NPZ9

# workout total maxn for each species ---
maxn.npz9.10<-maxn%>%
  dplyr::filter(location%in%"NPZ9")%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  ungroup()%>%
  top_n(11)%>%
  dplyr::filter(!scientific %in% c('Unknown spp', 'SUS sus'))%>%
  glimpse()

## Total frequency of occurance 
# I think we could remove this section - but maybe good to see sometimes
bar.npz9<-ggplot(maxn.npz9.10, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",position=position_dodge())+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  #scale_x_discrete(limits = rev(levels(scientific)))+
  #annotation_custom(lcpic, xmin=0.5, xmax=2.5, ymin=.75, ymax=1.5)+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_collapse
bar.npz9

## Top ten plot ----
bar.npz9.top.10<-ggplot(maxn.npz9.10%>%mutate(scientific=str_replace_all(.$scientific,
  c("miniatus"="miniatus*","auratus"="auratus*","dumerili"="dumerili*","multidens"="multidens*","Pseudocaranx spp"="Pseudocaranx spp*"))), 
  aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",colour="black",fill="lightgrey",position=position_dodge())+
  ylim (0, 175)+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"))+
  theme_collapse+
  theme.larger.text+
  annotation_raster(c.ass, xmin=9.75,xmax=10.25,ymin=168, ymax=183)+
  annotation_raster(c.w, xmin=8.8,xmax=9.2,ymin=138, ymax=153)+
  annotation_raster(l.m, xmin=7.5, xmax=8.5, ymin=75, ymax=115)+
  annotation_raster(c.au, xmin=6.5,xmax=7.5,ymin=73, ymax=117)+
  annotation_raster(s.d, xmin=5.5,xmax=6.5,ymin=30, ymax=80)+
  annotation_raster(p.m, xmin=4.5,xmax=5.5,ymin=29, ymax=75)+          
  annotation_raster(c.a, xmin=3.65,xmax=4.25,ymin=20, ymax=50)+
  annotation_raster(p.spp, xmin=2.65,xmax=3.35,ymin=14, ymax=46)+
  annotation_raster(p.s, xmin=1.65,xmax=2.35,ymin=13, ymax=42)
# bar.npz9.top.10

#save out plot
ggsave("plots/fish/stacked.bar.plot.npz9.png",bar.npz9.top.10,dpi=600,width=6.0, height = 6)

#Recreationally targeted species
#targeted species top 10 abundance
# Read in life history
url <- "https://docs.google.com/spreadsheets/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?ts=5e6f36e2#gid=825736197"

master <- googlesheets4::read_sheet(url)%>%
  ga.clean.names()%>%
  filter(grepl('Australia', global.region))%>%
  # filter(grepl('NW', marine.region))%>%
  dplyr::select(family,genus,species,iucn.ranking,fishing.mortality,fishing.type,australian.common.name,minlegal.wa)%>% 
  distinct()%>%
  glimpse()

fished.species <- maxn %>%
  dplyr::left_join(master) %>%
  dplyr::mutate(fishing.type = ifelse(scientific %in%c("Serranidae Plectropomus spp","Scombridae Scomberomorus spp","Lethrinidae Gymnocranius spp",
                                                       "Lethrinidae Lethrinus spp","Lethrinidae Unknown spp","Platycephalidae Platycephalus spp")
                                      ,"R",fishing.type))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Serranidae Plectropomus spp"), "450", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Scombridae Scomberomorus spp"), "900", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Gymnocranius spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Lethrinus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Unknown spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Platycephalidae Platycephalus spp"), "280", minlegal.wa))%>%
  dplyr::filter(fishing.type %in% c("B/R","B/C/R","R","C/R","C"))%>%
  dplyr::filter(!family%in%c("Monacanthidae", "Scorpididae", "Mullidae"))%>%    # Brooke removed leatherjackets, sea sweeps and goat fish
  dplyr::filter(!species%in%c("albimarginatus","longimanus"))%>%
  glimpse()

# workout total maxn for each species ---
#npz6
maxn.fished.npz6<-fished.species %>%
  dplyr::filter(location%in%"NPZ6")%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  # dplyr::filter(!scientific%in%"Lethrinus spp")%>%
  ungroup()%>%
  top_n(10)%>%
  glimpse()

#npz9
maxn.fished.npz9<-fished.species %>%
  dplyr::filter(location%in%"NPZ9")%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  group_by(scientific)%>%
  dplyr::summarise(maxn=sum(maxn))%>%
  ungroup()%>%
  top_n(10)%>%
  glimpse()

#have a look
#npz6
bar.npz6<-ggplot(maxn.fished.npz6, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",position=position_dodge())+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_collapse
bar.npz6

#npz9
bar.npz9<-ggplot(maxn.fished.npz9, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",position=position_dodge())+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_collapse
bar.npz9

#load fish pics
#NPZ6
#1 Lethrinus miniatus
#already loaded

#2 Chrysophrys auratus
#already loaded

#3 Choerodon rubescens
#already loaded

#4 Lethrinus nebulosus
#already loaded

#5 Seriola hippos
s.h <- as.raster(readPNG("data/images/Seriola_hippos_nb_HQ_TAYLOR.png"))

#6 Scomberomorus spp
s.spp <- as.raster(readPNG("data/images/Scombridae-Dark.png"))

#7 Lethrinus spp
l.spp <- as.raster(readPNG("data/images/Lethrinidae-Dark.png"))

#8 Lethrinus ravus
#as above

#9 Epinephelus rivulatus
e.r <- as.raster(readPNG("data/images/Serranidae-Dark.png"))

#10 Carcharhinus plumbeus
c.p <- as.raster(readPNG("data/images/Carcharinus plumbeus 5cmL adapt.png"))



#plot final bar plot
#npz6
bar.fished.npz6<-ggplot(maxn.fished.npz6, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",colour="black",fill="lightgrey",position=position_dodge())+
  ylim (0, 110)+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"))+
  theme_collapse+
  theme.larger.text+
  annotation_raster(l.m, xmin=9.6,xmax=10.4,ymin=90, ymax=115)+
  annotation_raster(c.au, xmin=8.6,xmax=9.4,ymin=80, ymax=105)+
  annotation_raster(c.r, xmin=7.65, xmax=8.35, ymin=50, ymax=75)+
  annotation_raster(l.n, xmin=6.55,xmax=7.45,ymin=28, ymax=55)+
  annotation_raster(s.h, xmin=5.55,xmax=6.45,ymin=20, ymax=50)+
  annotation_raster(s.spp, xmin=4.55,xmax=5.45,ymin=10, ymax=45)+
  annotation_raster(l.spp, xmin=3.7,xmax=4.3,ymin=7, ymax=25)+
  annotation_raster(l.spp, xmin=2.7,xmax=3.3,ymin=7, ymax=25)+
  annotation_raster(e.r, xmin=1.75,xmax=2.25,ymin=5, ymax=20)+
  annotation_raster(c.p, xmin=0.4,xmax=1.6,ymin=5, ymax=50)
# ggtitle("10 most abundant species") +
# theme(plot.title = element_text(hjust = 0))
bar.fished.npz6
  
#save out plot
ggsave("plots/fish/abundant.targets.npz6.png",bar.fished.npz6,dpi=600,width=6.0, height = 6.0)

#NPZ9
#1 Lethrinus miniatus
#already loaded
  
#2 Chrysophrus auratus
#already loaded
  
#3 Pristipomoides multidens
#already loaded
  
#4 Gymnocranius grandoculis
g.g <- as.raster(readPNG("data/images/Gymnocranius_grandoculis_nb_TAYLOR.png"))

#5 Seriola hippos
#already loaded
  
#6 Pristipomoides typus
p.t <- as.raster(readPNG("data/images/Pristipomoides typus 3cm.png"))
  
#7 Lethrinus nebulosus
 #already loaded
  
#8 Gymnocranius spp
#same as g.g
  
#9 Carcharhinus plumbeus
#already loaded
  
#10 Glaucosoma buergeri
g.b <- as.raster(readPNG("data/images/glaucosoma_buergeri_nb.png"))
  

#npz9
bar.fished.npz9<-ggplot(maxn.fished.npz9, aes(x=reorder(scientific,maxn), y=maxn)) +   
  geom_bar(stat="identity",colour="black",fill="lightgrey",position=position_dodge())+
  ylim (0, 92)+
  coord_flip()+
  xlab("Species")+
  ylab(expression(Overall~abundance~(Sigma~MaxN)))+
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"))+
  theme_collapse+
  theme.larger.text+
  annotation_raster(l.m, xmin=9.55,xmax=10.45,ymin=73, ymax=93)+          #1
  annotation_raster(c.au, xmin=8.55,xmax=9.45,ymin=73, ymax=95)+               #2
  annotation_raster(p.m, xmin=7.55, xmax=8.45, ymin=30, ymax=53)+         #3
  annotation_raster(g.g, xmin=6.55,xmax=7.45,ymin=8, ymax=25)+               #4
  annotation_raster(s.h, xmin=5.6,xmax=6.3,ymin=7, ymax=33)+                #5
  annotation_raster(p.t, xmin=4.6,xmax=5.4,ymin=7, ymax=27)+                 #6
  annotation_raster(l.n, xmin=3.55,xmax=4.45,ymin=5, ymax=25)+                 #7
  annotation_raster(g.g, xmin=2.55,xmax=3.45,ymin=5, ymax=20)+              #8
  annotation_raster(c.p, xmin=1.5,xmax=2.5,ymin=5, ymax=35)+                #9
  annotation_raster(g.b, xmin=0.55,xmax=1.45,ymin=4, ymax=20)                 #10
# ggtitle("10 most abundant species") +
# theme(plot.title = element_text(hjust = 0))
# bar.fished.npz9

# Save out plot
ggsave("plots/fish/abundant.targets.npz9.png",bar.fished.npz9,dpi=600,width=6.0, height = 6.0)






