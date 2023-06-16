###############################
# Script to bring in data from
# Abrolhos and Ningaloo and
# format it for plotting
###############################

## For data wrangling/cleaning
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
## For plotting
library(RColorBrewer)
## For life history info
library(cli)
library(devtools)
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)
library(httpuv)
library(googlesheets4)


rm(list = ls())


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")

#### LOAD FILES ####
setwd(data_dir)

## Abrolhos

# BOSS
ab_boss_lengths <- read.csv("2021-05_Abrolhos_BOSS.complete.length.csv") %>% 
  mutate(method = "BOSS")
ab_boss_meta <- read.csv("2021-05_Abrolhos_BOSS.checked.metadata.csv") %>% 
  mutate(method = "BOSS")


#BRUV
ab_bruv_lengths <- read.csv("2021-05_Abrolhos_stereo-BRUVs.complete.length.csv") %>% 
  mutate(method = "BRUV")
ab_bruv_meta <- read.csv("2021-05_Abrolhos_stereo-BRUVs.checked.metadata.csv") %>% 
  mutate(method = "BRUV")

## Ningaloo
# BRUV
ni_bruv_lengths <- read.csv("Ningaloo_PtCloates_stereo-BRUVs.complete.length.csv") %>% 
  mutate(method = "BRUV") %>% 
  ga.clean.names()

ni_bruv_meta<- read.csv("Ningaloo_PtCloates_stereo-BRUVs.checked.metadata.csv") %>% 
  mutate(method = "BRUV") %>% 
  ga.clean.names()

#BOSS
ni_boss_lengths <- read.csv("Ningaloo_PtCloates_BOSS.complete.length.csv") %>% 
  mutate(method = "BRUV") %>% 
  ga.clean.names()

ni_boss_meta <- read.csv("Ningaloo_PtCloates_BOSS.checked.metadata.csv") %>% 
  mutate(method = "BRUV") %>% 
  ga.clean.names()

## Length info from fishbase 
url <- "https://docs.google.com/spreadsheets/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?ts=5e6f36e2#gid=825736197"

master<-googlesheets4::read_sheet(url)%>%ga.clean.names()%>%
  filter(grepl('Australia', global.region))%>% # Change country here
  filter(grepl('NW', marine.region))%>% # Select marine region (currently this is only for Australia)
  dplyr::mutate(all=as.numeric(all))%>%
  dplyr::mutate(bll=as.numeric(bll))%>%
  dplyr::mutate(a=as.numeric(a))%>%
  dplyr::mutate(b=as.numeric(b))%>%
  dplyr::select(family,genus,species,marine.region,length.measure,a,b,all,bll,fb.length_max,fb.ltypemaxm, minlegal.wa, maxlegal.wa, fishing.type)%>%
  distinct()%>%
  glimpse()

##### FORMAT DATA ####

## Abrolhos
ab_metadata <- rbind(ab_boss_meta, ab_bruv_meta)

ab_all_species <- rbind(ab_boss_lengths, ab_bruv_lengths) %>% 
  filter(successful.length %in% c("Y")) %>% 
  mutate(scientific = paste(genus, species, sep = " ")) %>% 
  full_join(., master, by=c("family","genus","species")) %>% 
  dplyr::mutate(fishing.type = ifelse(scientific %in%c("Serranidae Plectropomus spp","Scombridae Scomberomorus spp","Lethrinidae Gymnocranius spp",
                                                       "Lethrinidae Lethrinus spp","Lethrinidae Unknown spp","Platycephalidae Platycephalus spp"),"R",fishing.type))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Serranidae Plectropomus spp"), "450", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Scombridae Scomberomorus spp"), "900", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Gymnocranius spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Lethrinus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Unknown spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Platycephalidae Platycephalus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(target = ifelse(fishing.type %in% c("B/R","B/C/R","R","C/R","C"), "Y", "N")) %>% 
  dplyr::mutate(target = ifelse(family %in% c("Monacanthidae", "Scorpididae", "Mullidae"), "N", target)) %>% 
  dplyr::mutate(target = ifelse(species %in% c("albimarginatus","longimanus"), "N", target)) %>% 
  dplyr::mutate(minlegal.wa = as.double(minlegal.wa)) %>%
  glimpse()

ab_all_species <- ab_all_species %>%
  tidyr::replace_na(list(minlegal.wa=0)) %>%
  dplyr::mutate(size = ifelse(length > minlegal.wa, "greater than legal size", "smaller than legal size")) %>%
  filter(!is.na(length)) %>% 
  dplyr::select(campaignid, sample, family, genus, species, scientific, length, number, location, status, method, target, size, fishing.type, minlegal.wa) %>% 
  dplyr::glimpse()

ab_species_lengths <- ab_all_species 


## Ningaloo
ni_metadata <- rbind(ni_boss_meta, ni_bruv_meta)

ni_all_species <- rbind(ni_boss_lengths, ni_bruv_lengths) %>% 
  filter(successful.length %in% c("Y")) %>% 
  mutate(scientific = paste(genus, species, sep = " ")) %>% 
  dplyr::left_join(master) %>%
  dplyr::mutate(fishing.type = ifelse(scientific %in% c("Serranidae Plectropomus spp","Scombridae Scomberomorus spp",
                                                        "Lethrinidae Gymnocranius spp","Lethrinidae Lethrinus spp",
                                                        "Lethrinidae Unknown spp","Platycephalidae Platycephalus spp", 
                                                        "Lutjanidae Pristipomoides spp", "Lutjanidae Pristipomoides sp1",
                                                        "Lethrinidae Gymnocranius sp1")
                                      ,"R",fishing.type))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Serranidae Plectropomus spp"), "450", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Scombridae Scomberomorus spp"), "900", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Gymnocranius spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Gymnocranius sp1"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Lethrinus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Lethrinidae Unknown spp"), "280", minlegal.wa))%>%
  dplyr::mutate(minlegal.wa = ifelse(scientific %in% c("Platycephalidae Platycephalus spp"), "280", minlegal.wa))%>%
  dplyr::mutate(target = ifelse(fishing.type %in% c("B/R","B/C/R","R","C/R","C"), "Y", "N"))%>%
  dplyr::mutate(target = ifelse(family %in% c("Monacanthidae", "Scorpididae", "Mullidae", 
                             "Carcharhinidae", "Sphyrnidae", "Pomacanthidae"), "N", target))%>%    # Remove non-targeted families   
  dplyr::mutate(minlegal.wa = as.double(minlegal.wa)) %>%
  glimpse()

ni_all_species <- ni_all_species %>%
  tidyr::replace_na(list(minlegal.wa=0)) %>%
  dplyr::mutate(size = ifelse(length > minlegal.wa, "greater than legal size", "smaller than legal size")) %>%
  filter(!is.na(length)) %>% 
  dplyr::select(campaignid, sample, family, genus, species, scientific, length, number, location, status, method, target, size, fishing.type, minlegal.wa) %>% 
  dplyr::glimpse()


