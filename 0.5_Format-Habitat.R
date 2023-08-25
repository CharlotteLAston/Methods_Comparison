###
# Project: Parks - Abrolhos Post-Survey
# Data:    BRUVS, BOSS Habitat data
# Task:    Merging habitat data
# author:  Charlotte (Kingsley Griffin, Brooke Gibbons & Claude)
# date:    July-Oct 2021
##

rm(list = ls())

library(reshape2)
library(dplyr)
library(raster)
library(sp)
library(ggplot2)

#### SET DIRECTORIES AND READ IN DATA ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "EM Export", sep="/")
sg_dir <- paste(working.dir, "staging", sep="/")
sp_dir <- paste(working.dir, "spatial_data", sep="/")

name <- "Ningaloo_PtCloates_BOSS-BRUV" 

# read in data
setwd(data_dir)

bosmet<- read.csv("Ningaloo_PtCloates_BOSS.checked.metadata.csv")                              # boss metadata
bruvmet <- read.csv("Ningaloo_PtCloates_stereo-BRUVs.checked.metadata.csv")

setwd(out_dir)

# BOSS Habitat
boshab.sq.2021 <- read.table('2021-08_PtCloates_Squid-BOSS_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2021-08_PtCloates_BOSS")
boshab.fl.2021 <- read.table('2021-08_PtCloates_Flasher-BOSS_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2021-08_PtCloates_BOSS")
boshab.sq.rel.2021 <- read.table('2021-08_PtCloates_Squid-BOSS_Relief_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2021-08_PtCloates_BOSS")
boshab.fl.rel.2021 <- read.table('2021-08_PtCloates_Flasher-BOSS_Relief_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2021-08_PtCloates_BOSS")

boshab.sq.2022 <- read.table('2022-05_PtCloates_Squid-BOSS_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2022-05_PtCloates_BOSS")
boshab.na.2022 <- read.table('2022-05_PtCloates_Naked-BOSS_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2022-05_PtCloates_BOSS")
boshab.sq.rel.2022 <- read.table('2022-05_PtCloates_Squid-BOSS_Relief_Dot Point Measurements.txt',skip = 5, sep = "\t") %>% 
  mutate(campaignid = "2022-05_PtCloates_BOSS")
boshab.na.rel.2022 <- read.table('2022-05_PtCloates_Naked-BOSS_Relief_Dot Point Measurements.txt',skip = 5, sep = "\t")%>% 
  mutate(campaignid = "2022-05_PtCloates_BOSS")

fbrhab.2021 <- read.table('2021-08-PC-BRUVs-Forwards_Dot Point Measurements.txt',skip = 5, sep = "\t")                                      # habitat notes bruv forwards
bbrhab.2021 <- read.table('2021-08-PC-BRUVs-Backwards_Dot Point Measurements.txt', skip = 5, sep = "\t")                                      # habitat notes bruv backwards
fbrrel.2021 <- read.table('2021-08-PC-BRUVs-Forwards_Relief_Dot Point Measurements.txt', skip = 5, sep = "\t")                                      # habitat notes bruv forwards
bbrrel.2021 <- read.table('2021-08-PC-BRUVs-Backwards_Relief_Dot Point Measurements.txt',skip = 5, sep = "\t")                                      # habitat notes bruv backwards

fbrhab.2022 <- read.table('2022-05_PtCloates_stereo-BRUVS_Forwards_Dot Point Measurements.txt',skip = 5, sep = "\t")                                      # habitat notes bruv forwards
bbrhab.2022 <- read.table('2022-05_PtCloates_stereo-BRUVS_Backwards_Dot Point Measurements.txt', skip = 5, sep = "\t")                                      # habitat notes bruv backwards
fbrrel.2022 <- read.table('2022-05_PtCloates_stereo-BRUVS_Forwards_Relief_Dot Point Measurements.txt', skip = 5, sep = "\t")                                      # habitat notes bruv forwards
bbrrel.2022 <- read.table('2022-05_PtCloates_stereo-BRUVS_Backwards_Relief_Dot Point Measurements.txt',skip = 5, sep = "\t")  

#### BOSS DATA ####
# clean all and merge into a single habitat dataframe
# wrangling boss data first

#* Relief ####
bosrel <- rbind.fill(boshab.sq.rel.2021, boshab.fl.rel.2021, boshab.sq.rel.2022, boshab.fl.rel.2022)

bosrel <- bosrel %>% 
  dplyr::select(campaignid,V1, V21, V27)

colnames(bosrel) <- c("campaignid", "filename","relief.old", "relief.new")                                      
bosrel$relief.new <- substr(bosrel$relief.new, start = 1, stop = 3)                  # extract relief value
bosrel$relief.new <- as.numeric(gsub("\\.", "", bosrel$relief.new))  

bosrel <- bosrel %>% 
  mutate(relief.old = ifelse(is.na(relief.old), relief.new, relief.old)) %>% 
  dplyr::rename(relief = "relief.old") %>% 
  dplyr::select(-relief.new) 

bosrel <- bosrel %>% 
  mutate(filename = as.factor(filename))  
 
bosrel$filename <- gsub(".jpg", "", bosrel$filename)                             # extract sample id

bosrel <- bosrel %>% 
  dplyr::group_by(campaignid, filename) %>% 
  dplyr::summarise(sdrel  = sd(relief, na.rm = TRUE),
            relief = mean(relief, na.rm = TRUE)) %>% 
  mutate(filename = as.character(filename)) %>% 
  mutate(filename = str_replace(filename, "6-FLASH_PCFB15", "6-FLASH_PCFB15.2")) %>% 
  mutate(filename = str_replace(filename, "3-FLASH_PCFB15", "3-FLASH_PCFB15.1")) %>% 
  mutate(filename = str_replace(filename, "3-FLASH_PCFB1_", "3-FLASH_PCFB1.1_")) %>% 
  mutate(filename = str_replace(filename, "4-FLASH_PCFB1_", "4-FLASH_PCFB1.2_")) %>% 
  mutate(sample = str_replace(filename, "SQUID-1", "SQUID")) %>% 
  mutate(sample = str_extract(sample, "(?<=(SQUID|FLASH|NAKED)\\_)[^\\_]+"))

#* Habitat data ####

summary(boshab)
boshab.2021 <- rbind.fill(boshab.sq.2021, boshab.fl.2021)

boshab.2021 <- boshab.2021[ , c(1, 4, 5, 18:21, 24, 26, 27)]                                  # omit bare columns
colnames(boshab.2021) <- c("filename", "image row", "image col", "broad",
                      "morphology", "type", "FOV", "CODE", "radius", "campaignid")            # fix colnames

boshab.2021 <- boshab.2021 %>%
  mutate(morphology = ifelse(broad %in% "Sessile invertebrates", "Mixed sessile invertebrates", morphology)) %>% 
  mutate(broad  = ifelse(broad %in% "Sessile invertebrates", "Matrix", broad)) %>% 
  mutate(morphology = ifelse(broad %in% "Unscorable", "Unscorable", morphology)) %>% 
  mutate(filename = str_replace(filename, "6-FLASH_PCFB15", "6-FLASH_PCFB15.2")) %>% 
  mutate(filename = str_replace(filename, "3-FLASH_PCFB15", "3-FLASH_PCFB15.1")) %>% 
  mutate(filename = str_replace(filename, "3-FLASH_PCFB1_", "3-FLASH_PCFB1.1_")) %>% 
  mutate(filename = str_replace(filename, "4-FLASH_PCFB1_", "4-FLASH_PCFB1.2_")) %>% 
  mutate(sample = str_replace(filename, "SQUID-1", "SQUID")) %>% 
  mutate(sample = str_extract(filename, "(?<=(SQUID|FLASH|NAKED)\\_)[^\\_]+"))

boshab.2022 <- rbind.fill(boshab.sq.2022, boshab.na.2022)

boshab.2022 <- boshab.2022[ , c(1, 4, 5, 24:26, 28, 32, 33)]                                  # omit bare columns
colnames(boshab.2022) <- c("filename", "image row", "image col",
                           "morphology", "type", "FOV", "CODE","radius", "campaignid")  

boshab.2022 <- boshab.2022 %>% 
  mutate(broad = str_extract(morphology, "^[:alpha:]+(?=[:blank:])")) %>% 
  relocate(broad, .before=morphology) %>% 
  mutate(morphology = str_remove(morphology, "^[:alpha:]+[:blank:]\\W[:blank:]")) %>% 
  mutate(broad = ifelse(is.na(broad) & morphology %in% "Unconsolidated"|is.na(broad) & morphology %in% "Consolidated", "Substrate", broad)) %>% 
  mutate(broad = ifelse(is.na(broad) & morphology %in% "Hydroids", "Cnidaria", broad)) %>% 
  mutate(FOV = ifelse(morphology %in% "Octocoral/Black", type, FOV)) %>% 
  mutate(type = ifelse(morphology %in% "Octocoral/Black", "Black & Octocorals", type)) %>% 
  mutate(morphology = ifelse(morphology %in% "Octocoral/Black", "Corals", morphology)) %>% 
  mutate(broad = ifelse(morphology %in% "Corals", "Cnidaria", broad)) %>% 
  mutate(broad = ifelse(is.na(broad) & morphology %in% "Sponges", "Sponges", broad)) %>% 
  mutate(morphology = ifelse(broad %in% "Sponges" & morphology %in% "Sponges", type, morphology)) %>% 
  mutate(morphology = ifelse(morphology %in% ""|is.na(morphology), "Unscorable", morphology)) %>% 
  mutate(sample = str_replace(filename, "SQUID-1", "SQUID")) %>% 
  mutate(sample = str_extract(sample, "(?<=(SQUID|FLASH|NAKED)\\_)[^\\_]+"))

boshab <- rbind(boshab.2021, boshab.2022)

#* Put together and merge with metadata ####

bosmet <- bosmet[, colnames(bosmet) %in% c("sample", "sate", "time", "latitude",
                                           "longitude", "site", "sample",
                                           "location", "status", "depth",
                                           "type")] 

allbos <- merge(bosmet, boshab, by = "sample")                                  # metadata and habitat measures
allbos <- merge(allbos, bosrel, by = "sample")                                  # add relief
allbos$pa <- c(1)                                                           # presence column for summing
allbos$method <- c("BOSS")                                                      # method id
head(allbos)

#### BRUV ####
#* forwards relief - see notes above for explanation #####

head(fbrrel.2021)
summary(fbrrel.2021)
fbrrel.2021 <- fbrrel.2021[ , c(1, 21)]
colnames(fbrrel.2021) <- c("filename", "relief")
fbrrel.2021$relief <- as.numeric(gsub("\\.", "", fbrrel.2021$relief))
fbrrel.2021 <- fbrrel.2021 %>% 
  mutate(sample = str_extract(filename, "^[:alnum:]+(?=\\_)")) %>% 
  mutate(sample = ifelse(is.na(sample), filename, sample))

head(fbrrel.2022)
summary(fbrrel.2022)
fbrrel.2022 <- fbrrel.2022[ , c(1, 27)]
colnames(fbrrel.2022) <- c("filename", "relief")
fbrrel.2022$relief    <- substr(fbrrel.2022$relief, start = 1, stop = 3)  
fbrrel.2022$relief <- as.numeric(gsub("\\.", "", fbrrel.2022$relief))
fbrrel.2022 <- fbrrel.2022 %>% 
  mutate(sample = str_extract(filename, "^[:alnum:]+(?=\\_)")) %>% 
  mutate(sample = ifelse(is.na(sample), filename, sample))

fbrrel <- rbind(fbrrel.2021, fbrrel.2022)

#* backwards relief ####
head(bbrrel)
bbrrel.2021 <- bbrrel.2021[ , c(1, 21)]
colnames(bbrrel.2021) <- c("filename", "relief")
bbrrel.2021$sample <- gsub(".jpg", "", bbrrel.2021$filename, ignore.case = TRUE)

bbrrel.2022 <- bbrrel.2022[ , c(1, 27)]
colnames(bbrrel.2022) <- c("filename", "relief")
bbrrel.2022$relief    <- substr(bbrrel.2022$relief, start = 1, stop = 3)  
bbrrel.2022$relief <- as.numeric(gsub("\\.", "", bbrrel.2022$relief))
bbrrel.2022$sample <- gsub(".jpg", "", bbrrel.2022$filename, ignore.case = TRUE)

bbrrel <- rbind(bbrrel.2021, bbrrel.2022)

#* merge and calc sample mean relief ####
buvrel <- rbind(fbrrel, bbrrel) %>% 
  mutate(filename = as.factor(filename))
buvrel <- buvrel %>% 
  dplyr::group_by(filename) %>% 
  dplyr::summarise(sdrel  = sd(relief, na.rm = TRUE),
            relief = mean(relief, na.rm = TRUE)) %>% 
  mutate(sample = str_extract(filename, "^[:alnum:]+(?=\\_|\\.)"))
  
head(buvrel)

#* bruv habitat ####
## forwards

fbrhab.old <- fbrhab.2021[ , c(1, 4, 5, 18:21,24, 26)]
colnames(fbrhab.old) <- c("filename", "image row", "image col", "broad",
                      "morphology", "type", "FOV", "CODE", "radius")
fbrhab.old$filename <- gsub(".jpg", "", fbrhab.old$filename, ignore.case = TRUE)

fbrhab.old <- fbrhab.old %>% 
  mutate(broad = ifelse(morphology %in% "Hydroids", "Cnidaria", broad)) %>% 
  mutate(morphology = ifelse(broad %in% "Sessile invertebrates", "Mixed sessile invertebrates", morphology)) %>% 
  mutate(broad  = ifelse(broad %in% "Sessile invertebrates", "Matrix", broad)) %>% 
  mutate(broad = ifelse(broad %in% "Sessile invertebrates", "Matrix", broad))

fbrhab.new <- fbrhab.2022[ , c(1, 4, 5, 23,24:26,28, 32)]
colnames(fbrhab.new) <- c("filename", "image row", "image col","broad",
                          "morphology", "type", "FOV", "CODE", "radius")
fbrhab.new$filename <- gsub(".jpg", "", fbrhab.new$filename, ignore.case = TRUE)

fbrhab.new <- fbrhab.new %>% 
  mutate(broad = ifelse(morphology %in% "Hydroids", "Cnidaria", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Invertebrate Complex", "Matrix", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Unconsolidated", "Substrate", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Sponges", "Sponges", broad)) %>% 
  mutate(morphology = ifelse(morphology %in% "Sponges", type, morphology)) %>% 
  mutate(type = ifelse(morphology %in% "Octocoral/Black", morphology, type)) %>% 
  mutate(morphology = ifelse(type %in% "Octocoral/Black", "Coral", morphology)) %>% 
  mutate(broad = ifelse(morphology %in% "Coral", "Cnidaria", broad))

fbrhab <- rbind.fill(fbrhab.old, fbrhab.new) %>% 
  mutate(sample = str_extract(filename, "^[:alnum:]+(?=\\_)")) %>% 
  mutate(sample = ifelse(is.na(sample), filename, sample))

## backwards
summary(bbrhab.2021)
bbrhab.old <- bbrhab.2021[ , c(1, 4, 5, 18:21, 24, 26)]
colnames(bbrhab.old) <- c("filename", "image row", "image col", "broad",
                      "morphology", "type", "FOV", "CODE", "radius")
bbrhab.old$sample   <- gsub(".jpg", "", bbrhab.old$filename, ignore.case = TRUE)

bbrhab.old <- bbrhab.old %>% 
  mutate(morphology = ifelse(broad %in% "Sessile invertebrates", "Mixed sessile invertebrates", morphology)) %>% 
  mutate(broad  = ifelse(broad %in% "Sessile invertebrates", "Matrix", broad))
  
summary(bbrhab.2022)
bbrhab.new <- bbrhab.2022[ , c(1, 4, 5, 23,24:26, 28, 32)]
colnames(bbrhab.new) <- c("filename", "image row", "image col", "broad",
                          "morphology", "type", "FOV", "CODE", "radius")
bbrhab.new$sample    <- gsub(".jpg", "", bbrhab.new$filename, ignore.case = TRUE)

bbrhab.new <- bbrhab.new %>% 
  mutate(broad = ifelse(morphology %in% "Hydroids", "Cnidaria", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Invertebrate Complex", "Matrix", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Unconsolidated", "Substrate", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Sponges", "Sponges", broad)) %>% 
  mutate(morphology = ifelse(morphology %in% "Sponges", type, morphology)) %>% 
  mutate(type = ifelse(morphology %in% "Octocoral/Black", morphology, type)) %>% 
  mutate(morphology = ifelse(type %in% "Octocoral/Black", "Coral", morphology)) %>% 
  mutate(broad = ifelse(morphology %in% "Coral", "Cnidaria", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Bryozoa", "Matrix", broad))


bbrhab <- rbind.fill(bbrhab.old, bbrhab.new)

#* join both bruv views and merge with metadata ####
buvhab <- rbind(bbrhab, fbrhab)
bruvmet <- bruvmet[, colnames(bruvmet) %in% c("campaignid","sample", "date", "time", "latitude",
                                           "longitude", "site", "sample",
                                           "location", "status", "depth",
                                           "type")]
allbuv <- merge(bruvmet, buvhab, by = "sample")
allbuv <- merge(allbuv, buvrel, by = "sample")
allbuv$pa     <- c(1)
allbuv$method <- c("BRUV")
head(allbuv)

#### Join both methods ####
allbos <- allbos %>% 
  dplyr::select(sample, latitude, longitude, location, status, depth, broad, morphology, type, FOV,
                campaignid.x, sdrel, relief, pa, method) %>% 
  dplyr::rename(campaignid = "campaignid.x")
allbuv <- allbuv %>% 
  dplyr::select(sample, latitude, longitude, location, status, depth, broad, morphology, type, FOV,
                campaignid, sdrel, relief, pa, method)

allhab <- rbind(allbos, allbuv) %>% 
  mutate(morphology = ifelse(morphology %in% "", "Unknown", morphology)) %>% 
  mutate(broad = ifelse(morphology %in% "Unscorable", "Unscorable", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Unknown", "Unknown", broad)) %>% 
  mutate(broad = ifelse(morphology %in% "Open Water", "Open Water", broad))


#### Long to wide and summarise ####
allhabw <- reshape2::dcast(allhab, campaignid + sample + method + location + latitude + longitude + depth + relief + sdrel ~ broad + morphology + type,
                           value.var = "pa", fun.aggregate = sum, drop = TRUE)
allhabw$totalpts <- rowSums(allhabw[ ,10:48]) - c(allhabw$Unknown_Unknown_ + allhabw$`Open Water_Open Water_` + allhabw$Unknown_Unknown_Unknown)
head(allhabw)
allhabw[(allhabw$totalpts - allhabw$Unconsolidated_Sand) < 0, ]

allhabw

# data checks (Brooke)
# check for habitat data that is missing metadata
allmet <- rbind.fill(bosmet, bruvmet)

t1 <- dplyr::anti_join(allhabw, allmet, by="sample") # none

# check for samples in metadata missing habitat
t2 <- dplyr::anti_join(allmet, allhabw, by="sample") # none


## extract bathy derivatives for modelling
# spatial setup
wgscrs <- CRS("+proj=longlat +datum=WGS84 +south")
utmcrs <- CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs")

setwd(sp_dir)
preds  <- readRDS("spatial_covariates.rds")
crs(preds)

# align crs and check samples over bathy and extract terrain data
allhab_sp <- SpatialPointsDataFrame(coords = allhabw[6:5], data = allhabw,
                                    proj4string = wgscrs) 
allhab_t  <- spTransform(allhab_sp, CRS = utmcrs)
plot(preds[[1]])
plot(allhab_t, add=T)
habt_df   <- as.data.frame(allhab_t, xy = T)
habi_df   <- cbind(habt_df, raster::extract(preds, allhab_t))

# Get combined habitat columns for modelling, same as section below
allhab <- habi_df %>%
  mutate(biog = base::rowSums(dplyr::select(., dplyr::starts_with(c("Bryozoa", "Cnidaria", "Sponges", "Matrix", "Invertebrate"))))) %>% 
  mutate(rock = rowSums(dplyr::select(., dplyr::contains(c("Consolidated"))))) %>% 
  mutate(sand = rowSums(dplyr::select(., dplyr::contains(c("Unconsolidated"))))) %>% 
  dplyr::select(campaignid, sample, method, latitude, longitude, depth, relief, sdrel, slope,
                aspect, TPI, TRI, roughness, detrended, lineartrend, biog, sand, rock, totalpts)

# # check all this malarkey has lined up (should get zero rows here)
# habi_df[(habi_df$totalpts - habi_df$Consolidated_Rock) < 0,]
# 
# 
# # collate generalised habitat tags
# habi_df <- habi_df %>%
#   mutate(kelps = Macroalgae_Large.canopy.forming_Ecklonia.radiata) %>%
#   mutate(macroalgae = rowSums(habi_df[ , grep("Macroalgae", colnames(habi_df))])) %>%
#   mutate(biog = rowSums(habi_df[ , c(grep("Sponge", colnames(habi_df)),
#                                       grep("Invertebrate", colnames(habi_df)),
#                                       grep("coral", colnames(habi_df)),
#                                       grep("Ascidian", colnames(habi_df)),
#                                       grep("Bryozoa", colnames(habi_df)),
#                                       grep("Hydroid", colnames(habi_df)))]
#                           )) %>%
#   mutate(sand = rowSums(habi_df[ , grep("Unconsolidated", colnames(habi_df))])) %>%
#   # mutate(turf = rowSums(habi_df[ , grep("Turf", colnames(habi_df))])) %>%
#   mutate(rock = rowSums(habi_df[ , grep("Consolidated", colnames(habi_df))])) # Should be everything including rock w/ turf mat
# # # biogenic reef - superceded above
# # brfc <- colnames(habi_df[ , - c(1:10, 30, 38:ncol(habi_df),
# #                                grep("Unconsolidated", colnames(habi_df)),
# #                                grep("Consolidated", colnames(habi_df))
# #                                )])
# # habi_df <- habi_df %>%
# #   mutate(biog = rowSums(habi_df[ , colnames(habi_df) %in% brfc])) %>%
# #   mutate(sample = Sample)
# 
# habi_df[(habi_df$totalpts - habi_df$biog) < 0,]
setwd(data_dir)
saveRDS(allhab, paste0(name, sep="_", "merged_habitat.rds"))
