###
# Project: Methods Comparison Paper
# Data:    BOSS & BRUV fish from Abrolhos and Ningaloo 
# Task:    Join BOSS and BRUV, format data for fssGAM
# Author:  Charlotte (Claude, Brooke, Kingsley)
# Date:    July 2023
##

rm(list=ls())

# libraries----
#detach("package:plyr", unload=TRUE)#will error - don't worry
library(tidyr)
library(dplyr)
options(dplyr.width = Inf) #enables head() to display all coloums
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
# library(doParallel) #this can removed?
library(doSNOW)
library(gamm4)
library(RCurl) #needed to download data from GitHub
library(FSSgam)
library(GlobalArchive)
library(ggplot2)
library(stringr)
library(sf)

#### SET DIRECTORIES AND READ IN DATA ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")
sg_dir <- paste(working.dir, "staging", sep="/")
sp_dir <- paste(working.dir, "spatial_data", sep="/")

#### ABROLHOS ####

setwd(working.dir)
name <- "2021-05_Abrolhos_BOSS-BRUV"  # set study name

# load and join datasets
#MaxN
setwd(data_dir)
boss.maxn   <- read.csv("2021-05_Abrolhos_BOSS.complete.maxn.csv")%>%
  dplyr::mutate(method = "BOSS")%>%
  glimpse()
bruv.maxn <- read.csv("2021-05_Abrolhos_stereo-BRUVs.complete.maxn.csv")%>%
  dplyr::mutate(method = "BRUV")%>%
  glimpse()
#join
maxn <- bind_rows(boss.maxn,bruv.maxn)%>%
  glimpse()


#length
boss.length <- read.csv("2021-05_Abrolhos_BOSS.complete.length.csv")%>%
  dplyr::mutate(method = "BOSS")%>%
  glimpse()
bruv.length <- read.csv("2021-05_Abrolhos_stereo-BRUVs.complete.length.csv")%>%
  dplyr::mutate(method = "BRUV")%>%
  glimpse()
#join
length <- bind_rows(boss.length,bruv.length)%>%
  dplyr::mutate(scientific = paste(family,genus,species, sep = " "))%>%
  glimpse()

#habitat
allhab <- readRDS("merged_habitat.rds")%>%
  dplyr::select(-status) %>%
  ga.clean.names()%>%
  glimpse()

allhab <- allhab %>%
  transform(kelps = kelps / broad.total.points.annotated) %>%
  transform(macroalgae = macroalgae / broad.total.points.annotated) %>%
  transform(sand = sand / broad.total.points.annotated) %>%
  transform(rock = rock / broad.total.points.annotated) %>%
  transform(biog = biog / broad.total.points.annotated) %>%
  glimpse()

names(maxn)

metadata <- maxn %>%
  distinct(sample, method,latitude, longitude, date, time, location, status, site, 
           depth, observer, successful.count, successful.length)

# look at top species ----
maxn.sum <- maxn %>%
  mutate(scientific = paste(genus, species, sep = " ")) %>%
  group_by(scientific) %>%
  dplyr::summarise(maxn = sum(maxn)) %>%
  top_n(10)%>%
  ungroup()

## Total frequency of occurrence
ggplot(maxn.sum, aes(x = reorder(scientific, maxn), y = maxn)) +   
  geom_bar(stat="identity",position = position_dodge()) +
  coord_flip() +
  xlab("Species") +
  ylab(expression(Overall ~ abundance ~ (Sigma ~ MaxN))) +
  #Theme1+
  theme(axis.text.y = element_text(face = "italic"))+
  #theme_collapse+
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))#+

# Create total abundance and species richness ----
ta.sr <- maxn %>%
  dplyr::ungroup() %>%
  dplyr::group_by(scientific,sample,method) %>%
  dplyr::summarise(maxn = sum(maxn)) %>%
  tidyr::spread(scientific, maxn, fill = 0) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(total.abundance = rowSums(.[, 3:(ncol(.))], na.rm = TRUE )) %>% #Add in Totals
  dplyr::mutate(species.richness = rowSums(.[, 3:(ncol(.))] > 0)) %>% # double check these
  dplyr::select(sample, total.abundance, species.richness,method) %>%
  tidyr::gather(., "scientific", "maxn", 2:3) %>%
  dplyr::rename(metric = "scientific") %>% 
  dplyr::glimpse()

dat.maxn <- bind_rows(ta.sr) %>%
  left_join(allhab) %>%
  left_join(metadata) %>%
  distinct()

# Set predictor variables---
names(dat.maxn)
names(allhab)

pred.vars = c("depth", 
              "macroalgae", 
              "sand", 
              "biog", 
              "mean.relief",
              "tpi",
              "roughness",
              "detrended") 

# predictor variables Removed at first pass---
# broad.Sponges and broad.Octocoral.Black and broad.Consolidated 

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat.maxn[,pred.vars]), 2)
# nothing is highly correlated 

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow = c(3, 2))
for (i in pred.vars) {
  x <- dat.maxn[ , i]
  x = as.numeric(unlist(x)) 
  hist((x)) #Looks best
  plot((x), main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x + 1))
  plot(log(x + 1))
}

#all looks fine
#write data to load in to next script
setwd(data_dir)
saveRDS(dat.maxn, "dat_maxn.rds")

#### GET LEGAL LENGTHS FOR FISHING ####
url <- "https://docs.google.com/spreadsheets/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?ts=5e6f36e2#gid=825736197"

master<-googlesheets4::read_sheet(url)%>%
  ga.clean.names()%>%
  filter(grepl('Australia', global.region))%>% # Change country here
  dplyr::select(family,genus,species,fishing.type,australian.common.name,minlegal.wa)%>%
  distinct()%>%
  glimpse()

unique(master$fishing.type)

all.species <- length %>%
  dplyr::left_join(., master) %>%
  dplyr::mutate(fishing.type = ifelse(scientific %in%c("Serranidae Plectropomus spp","Scombridae Scomberomorus spp","Lethrinidae Gymnocranius spp",
                                                       "Lethrinidae Lethrinus spp","Lethrinidae Unknown spp","Platycephalidae Platycephalus spp")
                                      ,"R",fishing.type))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Serranidae Plectropomus spp"), "450", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Scombridae Scomberomorus spp"), "900", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Gymnocranius spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Lethrinus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Unknown spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Platycephalidae Platycephalus spp"), "280", minlegal.wa))%>%
  #dplyr::filter(fishing.type %in% c("B/R","B/C/R","R","C/R","C"))%>%
  #dplyr::filter(!family%in%c("Monacanthidae", "Scorpididae", "Mullidae"))%>%    # Brooke removed leatherjackets, sea sweeps and goat fish
  #dplyr::filter(!species%in%c("albimarginatus","longimanus")) %>%
  dplyr::mutate(minlegal.wa = as.double(minlegal.wa)) %>%
  # dplyr::select(-method) %>% 
  glimpse()

complete.length <- all.species %>%
  dplyr::right_join(metadata, by = c("sample", "method")) %>% # add in all samples
  dplyr::select(sample,scientific,length,number,method) %>%
  tidyr::complete(nesting(sample,method), scientific) %>%
  replace_na(list(number = 0)) %>% #we add in zeros - in case we want to calculate abundance of species based on a length rule (e.g. greater than legal size)
  dplyr::ungroup()%>%
  dplyr::filter(!is.na(scientific)) %>% # this should not do anything
  dplyr::left_join(.,metadata, by = c("sample", "method")) %>%
  dplyr::left_join(.,allhab, by=c("sample", "method", "latitude", "longitude", "date", "location", "site", "depth")) %>%
  dplyr::filter(successful.length%in%c("Y", "Yes")) %>%
  dplyr::mutate(scientific=as.character(scientific)) %>%
  dplyr::glimpse()

testboss <- complete.length %>%
  dplyr::filter(method%in%"BOSS")

testbruv <- complete.length %>%
  dplyr::filter(method%in%"BRUV")

length(unique(testboss$sample))
75*2
length(unique(testbruv$sample))
47*2

#### GET LENGTH AT MATURITY ####
# Read in simple life history and synonyms ----
maturity <- readRDS("maturity.RDS")

mat.dat <- maturity %>% 
  mutate(scientific = paste0(Family, sep=" ", Genus, sep=" ", Species)) %>% 
  filter(scientific %in% c("Sparidae Chrysophrys auratus", "Labridae Choerodon rubescens", "Lethrinidae Lethrinus miniatus",
                           "Glaucosomatidae Glaucosoma hebraicum", "NA Centroberyx gerrardi", "Carangidae Seriola hippos",
                           "Lutjanidae Pristipomoides multidens")) %>% 
  # Tim has decided he wants the values from the fish.indicator.species sheet so this is me manually changing them #
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Labridae Choerodon rubescens"), 37.9, FB.length.at.maturity.cm)) %>%   # Got this value from somewhere else on the internet, probs a fisheries paper
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Sparidae Chrysophrys auratus"), 36.5, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lethrinidae Lethrinus miniatus"), 37.2, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Glaucosomatidae Glaucosoma hebraicum"), 31.1, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("NA Centroberyx gerrardi"), 26.9, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lutjanidae Pristipomoides multidens"), 52.6, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lethrinidae Lethrinus nebulosus"), 35, FB.length.at.maturity.cm)) %>%
  filter(!is.na(FB.length.at.maturity.cm)) %>% 
  mutate(scientific = str_replace(scientific, "NA", "Berycidae"))

length.dat <- complete.length %>% 
  filter(scientific %in% c("Sparidae Chrysophrys auratus", "Labridae Choerodon rubescens", "Lethrinidae Lethrinus miniatus",
                           "Glaucosomatidae Glaucosoma hebraicum", "Berycidae Centroberyx gerrardi",
                           "Lutjanidae Pristipomoides multidens", "Lethrinidae Lethrinus nebulosus")) %>% 
  left_join(., mat.dat, by="scientific") %>% 
  rename(length.mat = "FB.length.at.maturity.cm") %>% 
  mutate(length.mat = length.mat*10) %>% 
  mutate(Maturity = ifelse(length < (length.mat*0.5), "less_50", 
                          ifelse(length>(length.mat*0.5) & length < length.mat, "greater_50_less_mat",
                                 ifelse(length > length.mat & length < (length.mat*1.25), "greater_mat_less_125","greater_mat_125")))) %>% 
  mutate(Maturity2 = ifelse(length<length.mat, "less_mat", "greater_mat")) 

# Set predictor variables---
names(complete.length)
names(allhab)

pred.vars = c("depth", 
              "macroalgae", 
              "sand", 
              "biog", 
              "mean.relief",
              "tpi",
              "roughness",
              "detrended") 

# predictor variables Removed at first pass---
# broad.Sponges and broad.Octocoral.Black and broad.Consolidated , "InPreds","BioTurb" are too rare

dat.length <- complete.length

# Check for correalation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat.length[,pred.vars]),2)
# nothing is highly correlated 

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in pred.vars) {
  x<-dat.length[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

# all looks fine
# write data to load in to next script
setwd(data_dir)
filename <- paste0(name, sep="_", "dat_length.rds")
saveRDS(length.dat, filename)

#### NINGALOO ####

setwd(working.dir)
name <- "Ningaloo_PtCloates_BOSS-BRUV"  # set study name

# load and join datasets
#MaxN
setwd(data_dir)
boss.maxn   <- read.csv("Ningaloo_PtCloates_BOSS.complete.maxn.csv")%>%
  dplyr::mutate(method = "BOSS")%>%
  glimpse() %>% 
  dplyr::select(!date)
bruv.maxn <- read.csv("Ningaloo_PtCloates_stereo-BRUVs.complete.maxn.csv")%>%
  dplyr::mutate(method = "BRUV")%>%
  glimpse() %>% 
  dplyr::select(!date)
#join
maxn <- bind_rows(boss.maxn,bruv.maxn)%>%
  glimpse()


#length
boss.length <- read.csv("Ningaloo_PtCloates_BOSS.complete.length.csv")%>%
  dplyr::mutate(method = "BOSS")%>%
  glimpse() %>% 
  dplyr::select(!date)
bruv.length <- read.csv("Ningaloo_PtCloates_stereo-BRUVs.complete.length.csv")%>%
  dplyr::mutate(method = "BRUV")%>%
  glimpse()%>% 
  dplyr::select(!date)
#join
length <- bind_rows(boss.length,bruv.length)%>%
  dplyr::mutate(scientific = paste(family,genus,species, sep = " "))%>%
  glimpse()

#habitat
allhab <- readRDS("merged_habitat.rds")%>%
  dplyr::select(-status) %>%
  ga.clean.names()%>%
  glimpse()

allhab <- allhab %>%
  transform(kelps = kelps / broad.total.points.annotated) %>%
  transform(macroalgae = macroalgae / broad.total.points.annotated) %>%
  transform(sand = sand / broad.total.points.annotated) %>%
  transform(rock = rock / broad.total.points.annotated) %>%
  transform(biog = biog / broad.total.points.annotated) %>%
  glimpse()

names(maxn)

metadata <- maxn %>%
  distinct(id, method,latitude, longitude, time, location, status, site, 
           depth, observer, successful.count, successful.length)

# look at top species ----
maxn.sum <- maxn %>%
  mutate(scientific = paste(genus, species, sep = " ")) %>%
  group_by(scientific) %>%
  dplyr::summarise(maxn = sum(maxn)) %>%
  top_n(10)%>%
  ungroup()

## Total frequency of occurrence
ggplot(maxn.sum, aes(x = reorder(scientific, maxn), y = maxn)) +   
  geom_bar(stat="identity",position = position_dodge()) +
  coord_flip() +
  xlab("Species") +
  ylab(expression(Overall ~ abundance ~ (Sigma ~ MaxN))) +
  #Theme1+
  theme(axis.text.y = element_text(face = "italic"))+
  #theme_collapse+
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))#+

# Create total abundance and species richness ----
ta.sr <- maxn %>%
  dplyr::ungroup() %>%
  dplyr::group_by(scientific,sample,method) %>%
  dplyr::summarise(maxn = sum(maxn)) %>%
  tidyr::spread(scientific, maxn, fill = 0) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(total.abundance = rowSums(.[, 3:125], na.rm = TRUE )) %>% #Add in Totals
  dplyr::mutate(species.richness = rowSums(.[, 3:(ncol(.))] > 0)) %>% # double check these
  dplyr::select(sample, total.abundance, species.richness,method) %>%
  tidyr::gather(., "scientific", "maxn", 2:3) %>%
  dplyr::rename(metric = "scientific") %>% 
  dplyr::glimpse()

dat.maxn <- bind_rows(ta.sr) %>%
  left_join(allhab) %>%
  left_join(metadata) %>%
  distinct()

# Set predictor variables---
names(dat.maxn)
names(allhab)

pred.vars = c("depth", 
              "macroalgae", 
              "sand", 
              "biog", 
              "mean.relief",
              "tpi",
              "roughness",
              "detrended") 

# predictor variables Removed at first pass---
# broad.Sponges and broad.Octocoral.Black and broad.Consolidated 

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat.maxn[,pred.vars]), 2)
# nothing is highly correlated 

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow = c(3, 2))
for (i in pred.vars) {
  x <- dat.maxn[ , i]
  x = as.numeric(unlist(x)) 
  hist((x)) #Looks best
  plot((x), main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x + 1))
  plot(log(x + 1))
}

#all looks fine
#write data to load in to next script
setwd(data_dir)
saveRDS(dat.maxn, "dat_maxn.rds")

#### GET LEGAL LENGTHS FOR FISHING ####

all.species <- length %>%
  dplyr::left_join(., master) %>%
  dplyr::mutate(fishing.type = ifelse(scientific %in%c("Serranidae Plectropomus spp","Scombridae Scomberomorus spp","Lethrinidae Gymnocranius spp",
                                                       "Lethrinidae Lethrinus spp","Lethrinidae Unknown spp","Platycephalidae Platycephalus spp")
                                      ,"R",fishing.type))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Serranidae Plectropomus spp"), "450", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Scombridae Scomberomorus spp"), "900", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Gymnocranius spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Lethrinus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Unknown spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Platycephalidae Platycephalus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = as.double(minlegal.wa)) %>%
  dplyr::select(-method) %>% 
  glimpse()

complete.length <- all.species %>%
  dplyr::right_join(metadata, by = c("id")) %>% # add in all samples
  mutate(scientific = ifelse(scientific %in% c("Lutjanidae Pristipomoides sp1"), "Lutjanidae Pristipomoides multidens", scientific)) %>% # Changes it to multidens to make it easy but it's really sp1
  dplyr::select(sample,id,campaignid,scientific,length,number,method) %>%
  tidyr::complete(nesting(sample,method), scientific) %>%
  replace_na(list(number = 0)) %>% #we add in zeros - in case we want to calculate abundance of species based on a length rule (e.g. greater than legal size)
  dplyr::ungroup()%>%
  dplyr::filter(!is.na(scientific)) %>% # this should not do anything
  dplyr::left_join(.,metadata, by="id") %>%
  dplyr::left_join(.,allhab) %>%
  dplyr::filter(successful.length%in%c("Yes", "Y")) %>%
  dplyr::mutate(scientific=as.character(scientific)) %>%
  dplyr::select(-method.y, -method) %>% 
  dplyr::rename(method = "method.x") %>% 
  dplyr::glimpse()

testboss <- complete.length %>%
  dplyr::filter(method%in%"BOSS")

testbruv <- complete.length %>%
  dplyr::filter(method%in%"BRUV")

length(unique(testboss$id))

length(unique(testbruv$id))


#### GET LENGTH AT MATURITY ####
# Read in simple life history and synonyms ----
maturity <- readRDS("maturity.RDS")

mat.dat <- maturity %>% 
  mutate(scientific = paste0(Family, sep=" ", Genus, sep=" ", Species)) %>% 
  filter(scientific %in% c("Sparidae Chrysophrys auratus", "Lethrinidae Lethrinus nebulosus",
                           "Lutjanidae Pristipomoides multidens", "Lethrinidae Lethrinus miniatus",
                           "Labridae Choerodon rubescens", "Lethrinidae Lethrinus punctulatus", 
                           "Serranidae Epinephelus multinotatus", "Serranidae Epinephelides armatus")) %>% 
  # Tim has decided he wants to use the values from the fish.indicator.species sheet so this is me going through and changing all the values #
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Sparidae Chrysophrys auratus"), 36.5, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lethrinidae Lethrinus nebulosus"), 35, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lutjanidae Pristipomoides multidens"), 52.6, FB.length.at.maturity.cm)) %>%
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lethrinidae Lethrinus miniatus"), 37.2, FB.length.at.maturity.cm)) %>% 
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Lethrinidae Lethrinus punctulatus"), 20.6, FB.length.at.maturity.cm)) %>%   # Got this value from somewhere else on the internet, probs a fisheries paper
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Labridae Choerodon rubescens"), 37.9, FB.length.at.maturity.cm)) %>%   # Got this value from somewhere else on the internet, probs a fisheries paper
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Serranidae Epinephelus multinotatus"), 39.1, FB.length.at.maturity.cm)) %>% 
  mutate(FB.length.at.maturity.cm = ifelse(scientific %in% c("Serranidae Epinephelides armatus"), 28.1, FB.length.at.maturity.cm)) %>% 
  filter(!is.na(FB.length.at.maturity.cm)) 

length.dat <- complete.length %>% 
  filter(scientific %in% c("Sparidae Chrysophrys auratus", "Lethrinidae Lethrinus nebulosus",
                           "Lutjanidae Pristipomoides multidens", "Lethrinidae Lethrinus miniatus",
                           "Labridae Choerodon rubescens", "Lethrinidae Lethrinus punctulatus", 
                           "Serranidae Epinephelus multinotatus", "Serranidae Epinephelides armatus")) %>% 
  left_join(., mat.dat, by="scientific") %>% 
  rename(length.mat = "FB.length.at.maturity.cm") %>% 
  mutate(length.mat = length.mat*10) %>% 
  mutate(Maturity = ifelse(length < (length.mat*0.5), "less_50", 
                           ifelse(length>(length.mat*0.5) & length < length.mat, "greater_50_less_mat",
                                  ifelse(length > length.mat & length < (length.mat*1.25), "greater_mat_less_125","greater_mat_125")))) %>% 
  mutate(Maturity2 = ifelse(length<length.mat, "less_mat", "greater_mat")) 

# Set predictor variables---
names(complete.length)
names(allhab)

pred.vars = c("depth", 
              "macroalgae", 
              "sand", 
              "biog", 
              "mean.relief",
              "tpi",
              "roughness",
              "detrended") 

# predictor variables Removed at first pass---
# broad.Sponges and broad.Octocoral.Black and broad.Consolidated , "InPreds","BioTurb" are too rare

dat.length <- length.dat

# Check for correalation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat.length[,pred.vars]),2)
# nothing is highly correlated 

# Plot of likely transformations - thanks to Anna Cresswell for this loop!
par(mfrow=c(3,2))
for (i in pred.vars) {
  x<-dat.length[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

#### GET STATUS FOR SAMPLES THAT ARE CURRENTLY NA ####
setwd(sp_dir)
NTZ <- st_read("NTZ and Fished areas for status.shp")%>% # Read in the Commonwealth NPZ so that we can work out which points are inside
  filter(ZoneName %in% c("National Park Zone")) %>% 
  st_transform(4283)
plot(NTZ$geometry)

points <- metadata %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs=4283) # Convert points to spatial data so we can see where they overlap
plot(points$geometry, add=T)

inside.points <- st_intersects(points,NTZ) %>% 
  as.data.frame() %>% 
  dplyr::select(row.id)

metadata.ntz <- metadata[as.numeric(inside.points$row.id), ]

length.dat <- dat.length %>% 
  mutate(status = ifelse(is.na(status) & id %in% c(metadata.ntz$id), "No-Take", "Fished"))

setwd(data_dir)
filename <- paste0(name, sep="_", "dat_length.rds") 
saveRDS(length.dat, filename)










