###
# Project: Parks - Abrolhos Post-Survey
# Data:    BRUVS, BOSS Habitat data
# Task:    Creating spatial covariates from bathy 
# author:  Charlotte (Kingsley Griffin, Brooke Gibbons & Claude)
# date:    July-Oct 2021
##

rm(list = ls())

# Load libraries
library(tidyverse)
library(sp)
library(raster)
library(terra)
library(sf)
library(stars)
library(starsExtra)
library(magrittr)
library(RNetCDF)
library(weathermetrics)
library(lubridate)
library(gstat)
library(rgdal)

#### SET DIRECTORIES AND READ IN DATA ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "EM Export", sep="/")
sg_dir <- paste(working.dir, "staging", sep="/")
sp_dir <- paste(working.dir, "spatial_data", sep="/")

name <- "Ningaloo_PtCloates_BOSS-BRUV" 

#### READ IN BATHY ####
setwd(sp_dir)

e <- ext(112.5, 114.7, -24, -20.5)

# Read in the bathymetry
setwd(sp_dir)
bathy <- rast("bath_250_good.tif") %>%
  crop(e) %>%
  clamp(upper = 0, lower = -300, values = F) %>%
  trim()
plot(bathy)

# Save the bathymetry data out as a dataframe
fbath <- bathy

fbath_df <- as.data.frame(bathy, xy = TRUE, na.rm = T) %>%                      # Convert to a dataframe
  dplyr::rename(Z = bath_250_good)

# transform bathy to projected coords for modelling
wgscrs  <- "+proj=longlat +datum=WGS84"
sppcrs  <- "+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs"     # crs for sp objects
proj4string(fbath) <- wgscrs
fbath_t <- fbath %>%  project("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs")
plot(fbath_t)

# calculate terrain on fine bathy
preds <- terrain(fbath_t, neighbors = 8,
                 v = c("slope", "aspect", "TPI", "TRI", "roughness"),           # Remove here as necessary
                 unit = "degrees")
preds <- rast(list(fbath_t, preds))  

# Calculate detrended bathymetry
zstar <- st_as_stars(fbath_t)                                                     # Convert to a stars object
detre <- detrend(zstar, parallel = 8)                                           # Detrend bathymetry - This usually runs quite slow!
detre <- as(object = detre, Class = "SpatRaster")                               # Convert it to a terra raster
names(detre) <- c("detrended", "lineartrend")
preds <- rast(list(preds, detre))                                               # Make a rasterstack
plot(preds)
names(preds)[1] <- "Z"
preds <- as(object = preds, Class = "Raster")

setwd(sp_dir)
saveRDS(preds, "spatial_covariates.rds")
