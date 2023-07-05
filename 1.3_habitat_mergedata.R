###
# Project: SW Habitat & Fish 
# Data:    BRUVS, BOSS Habitat data
# Task:    Merging habitat data
# Author:  Claude Spencer
# Date:    October 2022
##

# This script takes formatted habitat data from TransectMeasure and joins it with bathymetry derivatives for modelling

# Clear your environment
rm(list = ls())

# Load libraries
library(reshape2)
library(dplyr)
library(terra)
library(sf)
library(ggplot2)
library(magrittr)
library(RNetCDF)
library(weathermetrics)
library(lubridate)
library(stars)
library(stringr)

# Set your study name
name <- "sw-network"                                                  

# Load in tidy data from the formatting scripts
hab  <- read.csv(paste0("data/tidy/", name, "_broad-habitat.csv")) %>%
  glimpse()

test <- hab %>%
  dplyr::mutate(reef.prop = reef / total.pts) %>%
  dplyr::filter(reef.prop > 1) # 0 obs - v nice

# Set up CRS and load spatial covariates from 02_spatial_layers.R 
wgscrs <- "+proj=longlat +datum=WGS84 +south"                                   # Latlong projection

preds <- readRDS(paste(paste0('data/spatial/rasters/raw bathymetry/', name),    # This is ignored - too big!
                       'spatial_covariates.rds', sep = "_"))
# Extract bathy derivatives for modelling
# Align crs and check samples over bathy and extract terrain data
allhab_sp <- vect(hab, geom = c("longitude", "latitude"), crs = wgscrs)         # Convert the habitat data to a terra vector
plot(preds[[13]])                                                                # Plot the first bathymetry derivative
plot(allhab_sp, add = T)                                                        # Add the sampling points to check if they align
habt_df   <- as.data.frame(allhab_sp) %>%
  left_join(hab) %>%
  dplyr::select(longitude, latitude, everything()) %>%
  dplyr::rename(x = longitude, y = latitude)
habi_df   <- cbind(habt_df, terra::extract(preds, allhab_sp)) %>% 
  dplyr::filter(!is.na(Z),
                !is.na(PROD),
                !is.na(SST),
                !is.na(roughness))

# Save the output
saveRDS(habi_df, paste(paste0('data/tidy/', name), 
                       'habitat-bathy-derivatives.rds', sep = "_"))