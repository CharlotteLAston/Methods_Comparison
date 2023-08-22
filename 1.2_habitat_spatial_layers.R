###
# Project: Methods Comparison
# Data:    Bathymetry Data
# Task:    Prepare spatial layers for modelling - Just Ningaloo as this has already been done for Abrolhos
# Author:  Charlotte Aston (Claude Spencer)
# Date:    July 2023
## 

# This script formats bathymetry data and extracts bathymetry derivatives for modelling habitat and fish
# As this raw bathymetry data is often too large for GitHub, the raw files are hidden in the .gitnore
# You have to download data and create this folder directory yourself for the script to run!

# Clear your environment
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
out_dir <- paste(working.dir, "fssgam_output", sep="/")
sg_dir <- paste(working.dir, "staging", sep="/")
sp_dir <- paste(working.dir, "spatial_data", sep="/")

# Set your study name
name <- "Ningaloo_PtCloates_BOSS-BRUV"                        

e <- ext(112.5, 114.7, -24, -20.5)

# Read in the bathymetry
setwd(sp_dir)
bathy <- rast("bath_250_good.tif") %>%
  crop(e) %>%
  clamp(upper = 0, lower = -300, values = F) %>%
  trim()
plot(bathy)
summary(bathy)

# Save the bathymetry data out as a dataframe
fbath_df <- as.data.frame(bathy, xy = TRUE, na.rm = T) %>%                      # Convert to a dataframe
  dplyr::rename(Z = bath_250_good)

saveRDS(fbath_df, paste0(name, 'bathymetry.rds', sep = "_"))

#### CALCULATE BATHY DERIVATIVES ####
aus <- st_read("cstauscd_r.mif", crs = 4283) %>% #Shapefile of Australia
  dplyr::filter(!FEAT_CODE %in% "sea") %>%
  st_crop(e) %>%
  st_union(.) %>%
  st_cast("MULTILINESTRING") %>%
  st_transform(crs = 4326)
plot(aus)

openness <- st_as_sf(fbath_df, coords = c("x", "y"), crs = 4326) %>%
  dplyr::mutate(openness = st_distance(., aus)) %>%
  bind_cols(st_coordinates(.)) %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  dplyr::mutate(openness = as.numeric(openness)) %>%
  glimpse()

openness$openness <- as.numeric(openness$openness)
openness.rast <- rast(openness %>% dplyr::select(X, Y, openness), type = "xyz",
                      crs = "epsg:4326")
plot(openness.rast)
# Calculate TERRA terrain derivatives
preds <- terrain(bathy, neighbors = 8,
                 v = c("slope", "aspect", "TPI", "TRI", "roughness"),           # Remove here as necessary
                 unit = "degrees")              
preds <- rast(list(bathy, preds, openness.rast))                                # Stack the derivatives with the bathymetry

# Calculate detrended bathymetry
zstar <- st_as_stars(bathy)                                                     # Convert to a stars object
detre <- detrend(zstar, parallel = 8)                                           # Detrend bathymetry - This usually runs quite slow!
detre <- as(object = detre, Class = "SpatRaster")                               # Convert it to a terra raster
names(detre) <- c("detrended", "lineartrend")
preds <- rast(list(preds, detre))                                               # Make a rasterstack
plot(preds)
names(preds)[1] <- "Z"
preds <- as(object = preds, Class = "Raster")

# Oceanography data
# Sea level anomaly and currents
nc_sla <- open.nc("data/spatial/oceanography/IMOS_aggregation_20230418T233649Z.nc", write = TRUE)
print.nc(nc_sla) # shows you all the file details

time_nc <- var.get.nc(nc_sla, 'TIME')  # NC_CHAR time:units = "days since 1981-01-01 00:00:00" ;
time_nc_sla <- utcal.nc("days since 1985-01-01 00:00:00 UTC", time_nc,type = "c")
dates_sla <- as.Date(time_nc_sla)

lat <- var.get.nc(nc_sla, 'LATITUDE') # Some latitude, some lat -> watch for spelling
lon <- var.get.nc(nc_sla, 'LONGITUDE')
summary(lat)
summary(lon)
# Get lats of sst file which correspond to the lats of the zone
ext(preds)

lat_i <- which(lat <= ext(preds)[4] & lat >= ext(preds)[3])
check_lat <- lat[lat_i]
lon_i <- which(lon <= ext(preds)[2] & lon >= ext(preds)[1])
check_lon <- lon[lon_i]

#load only the subset of data
#get all sst
sla_all  <- var.get.nc(nc_sla,'GSLA', start = c(lon_i[1], lat_i[1],1),
                       count = c(length(lon_i), length(lat_i), length(dates_sla)));
ucur_all <- var.get.nc(nc_sla,'UCUR', start = c(lon_i[1], lat_i[1],1),
                       count = c(length(lon_i), length(lat_i), length(dates_sla)));
vcur_all <- var.get.nc(nc_sla,'VCUR', start = c(lon_i[1], lat_i[1],1),
                       count = c(length(lon_i), length(lat_i), length(dates_sla)));

time_data <- list()
time_data$dates <- as.character(dates_sla)
time_data$month <- lubridate::month(as.POSIXlt(time_data$dates, format="%Y-%m-%d"))
time_data <- as.data.frame(time_data)

lat_sla <- check_lat
lon_sla <- check_lon

# SEA LEVEL ANOMALY
# As an array
arr.sla = array(sla_all, dim = c(length(lon_i), length(lat_i), length(time_data$dates)),
                dimnames = list(lon_sla, lat_sla, time_data$dates))

# To a dataframe - 10 year average
sla.df <- arr.sla %>%
  reshape2::melt(varnames = c("lon","lat","date")) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::summarise(sla = mean(value, na.rm = T)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(sla)) %>%
  dplyr::select(lon, lat, sla) %>%
  glimpse()

# Rasterise
wgscrs <- "+proj=longlat +datum=WGS84 +no_defs"

sla.rast <- rasterFromXYZ(sla.df, crs = wgscrs) 
# To a sf object - for some reason not working as a dataframe
sla.sf <- st_as_sf(sla.df, coords = c("lon", "lat"), crs = 4326)

# Gstat Nearest Neighbour interpolation
gs.sla <- gstat(id = "sla", formula = sla~1, data = sla.sf, 
                nmax = 5, set = list(idp = 0))
sla.nn <- raster::interpolate(object = sla.rast, model = gs.sla, na.rm = T, debug.level = 0) %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::crop(preds[[1]]) %>%
  raster::mask(preds[[1]])
plot(sla.nn)

# Cover the NAs of the original dataframe with the interpolated values
sla.rast <- sla.rast %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::cover(sla.nn) %>%
  raster::crop(sla.nn) %>%
  raster::mask(sla.nn)
names(sla.rast) <- "SLA"
plot(sla.rast)

# VCUR
# As an array
arr.vcur = array(vcur_all, dim = c(length(lon_i), length(lat_i), length(time_data$dates)),
                 dimnames = list(lon_sla, lat_sla, time_data$dates))

# To a dataframe - 10 year average
vcur.df <- arr.vcur %>%
  reshape2::melt(varnames = c("lon","lat","date")) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::summarise(vcur = mean(value, na.rm = T)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(vcur)) %>%
  dplyr::select(lon, lat, vcur) %>%
  glimpse()

# Rasterise
vcur.rast <- rasterFromXYZ(vcur.df, crs = wgscrs) 
# To a sf object - for some reason not working as a dataframe
vcur.sf <- st_as_sf(vcur.df, coords = c("lon", "lat"), crs = 4326)

# Gstat Nearest Neighbour interpolation
gs.vcur <- gstat(id = "vcur", formula = vcur~1, data = vcur.sf, 
                 nmax = 5, set = list(idp = 0))
vcur.nn <- raster::interpolate(object = vcur.rast, model = gs.vcur, na.rm = T, debug.level = 0) %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::crop(preds[[1]]) %>%
  raster::mask(preds[[1]])
plot(vcur.nn)

# Cover the NAs of the original dataframe with the interpolated values
vcur.rast <- vcur.rast %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::cover(vcur.nn) %>%
  raster::crop(vcur.nn) %>%
  raster::mask(vcur.nn)
names(vcur.rast) <- "VCUR"
plot(vcur.rast)

# UCUR
# As an array
arr.ucur = array(ucur_all, dim = c(length(lon_i), length(lat_i), length(time_data$dates)),
                 dimnames = list(lon_sla, lat_sla,time_data$dates))

# To a dataframe - 10 year average
ucur.df <- arr.ucur %>%
  reshape2::melt(varnames = c("lon","lat","date")) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::summarise(ucur = mean(value, na.rm = T)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(ucur)) %>%
  dplyr::select(lon, lat, ucur) %>%
  glimpse()

# Rasterise
ucur.rast <- rasterFromXYZ(ucur.df, crs = wgscrs) 
# To a sf object - for some reason not working as a dataframe
ucur.sf <- st_as_sf(ucur.df, coords = c("lon", "lat"), crs = 4326)

# Gstat Nearest Neighbour interpolation
gs.ucur <- gstat(id = "ucur", formula = ucur~1, data = ucur.sf, 
                 nmax = 5, set = list(idp = 0))
ucur.nn <- raster::interpolate(object = ucur.rast, model = gs.ucur, na.rm = T, debug.level = 0) %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::crop(preds[[1]]) %>%
  raster::mask(preds[[1]])
plot(ucur.nn)

# Cover the NAs of the original dataframe with the interpolated values
ucur.rast <- ucur.rast %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::cover(ucur.nn) %>%
  raster::crop(ucur.nn) %>%
  raster::mask(ucur.nn)
names(ucur.rast) <- "UCUR"
plot(ucur.rast)

# Sea surface temperature
nc_sst <- open.nc("data/spatial/oceanography/IMOS_aggregation_20230418T233102Z.nc", write = TRUE)

print.nc(nc_sst) # shows you all the file details

time_nc <- var.get.nc(nc_sst, 'time')  # NC_CHAR time:units = "days since 1981-01-01 00:00:00" ;
time_nc_sst <- utcal.nc("seconds since 1981-01-01 00:00:00", time_nc,type = "c")
dates_sst <- as.Date(time_nc_sst)

lat <- var.get.nc(nc_sst, 'lat') # Some latitude, some lat -> watch for spelling
lon <- var.get.nc(nc_sst, 'lon')
summary(lat)
summary(lon)
# Get lats of sst file which correspond to the lats of the zone
ext(preds)

lat_i <- which(lat <= ext(preds)[4] & lat >= ext(preds)[3])
check_lat <- lat[lat_i]
lon_i <- which(lon <= ext(preds)[2] & lon >= ext(preds)[1])
check_lon <- lon[lon_i]

#load only the subset of data
#get all sst
sst_all  <- var.get.nc(nc_sst,'sea_surface_temperature', start = c(lon_i[1], lat_i[1],1),
                       count = c(length(lon_i), length(lat_i), length(dates_sst)));

time_data <- list()
time_data$dates <- as.character(dates_sst)
time_data$month <- lubridate::month(as.POSIXlt(time_data$dates, format="%Y-%m-%d"))
time_data <- as.data.frame(time_data)

lat_sst <- check_lat
lon_sst <- check_lon

# SEA SURFACE TEMPERATURE
# As an array
arr.sst = array(sst_all, dim = c(length(lon_i), length(lat_i), length(time_data$dates)),
                dimnames = list(lon_sst, lat_sst, time_data$dates))

# To a dataframe - 10 year average - winter average
sst.df <- arr.sst %>%
  reshape2::melt(varnames = c("lon","lat","date")) %>%
  dplyr::mutate(month = as.numeric(month(date))) %>%
  dplyr::filter(between(month, 6, 8)) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::summarise(sst = mean(value, na.rm = T)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(sst)) %>%
  dplyr::select(lon, lat, sst) %>%
  dplyr::mutate(sst = kelvin.to.celsius(sst)) %>%
  glimpse()

# Determine the extent of the raster
xmin <- min(sst.df$lon)
xmax <- max(sst.df$lon)
ymin <- min(sst.df$lat)
ymax <- max(sst.df$lat)

# Convert the meter resolution to degree resolution
# Function to convert degree distance to meter distance
lat <- mean(sst.df$lat)

# Chat GPT represent
degree_to_meter <- function(latitude, degree_distance) {
  # Radius of the Earth in meters
  earth_radius <- 6371000
  
  # Calculate the conversion factor based on the latitude
  conversion_factor <- 2 * pi * earth_radius * cos(latitude * pi / 180) / 360
  
  # Convert the degree distance to meters
  meter_distance <- degree_distance * conversion_factor
  
  # Return the result
  return(meter_distance)
}

# Example usage
degree_distance <- 0.2 # Distance in degrees

meter_distance <- degree_to_meter(lat, degree_distance)
print(meter_distance)

resolution_degrees <- meter_distance / 111000

# Calculate the number of rows and columns based on the resolution
ncols <- ceiling((xmax - xmin) / resolution_degrees)
nrows <- ceiling((ymax - ymin) / resolution_degrees)

# Create an empty raster with the specified resolution and extent
r <- raster(nrows = nrows, ncols = ncols, xmn = xmin, xmx = xmax, 
            ymn = ymin, ymx = ymax, crs = "+proj=longlat +datum=WGS84")

# Rasterize
sst.rast <- rasterize(cbind(sst.df$lon, sst.df$lat), r, field = sst.df$sst, fun = mean)
plot(sst.rast)

# To a sf object - for some reason not working as a dataframe
sst.sf <- st_as_sf(sst.df, coords = c("lon", "lat"), crs = 4326)

# Gstat Nearest Neighbour interpolation
gs.sst <- gstat(id = "sst", formula = sst~1, data = sst.sf, 
                nmax = 5, set = list(idp = 0))
sst.nn <- raster::interpolate(object = sst.rast, model = gs.sst, na.rm = T, debug.level = 0) %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::crop(preds[[1]]) %>%
  raster::mask(preds[[1]])
plot(sst.nn)

# Cover the NAs of the original dataframe with the interpolated values
sst.rast <- sst.rast %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::cover(sst.nn) %>%
  raster::crop(sst.nn) %>%
  raster::mask(sst.nn)
names(sst.rast) <- "SST"
plot(sst.rast)

# Primary productivity - only 1 variable in the .nc so just average across all layers RIP computer
prod.rast <- rast("data/spatial/oceanography/IMOS_aggregation_20230501T023239Z.nc") %>%
  mean(na.rm = T)
prod.rast <- as(object = prod.rast, Class = "Raster")
plot(prod.rast)
extent(preds)
prod.rast <- crop(prod.rast, extent(preds))
names(prod.rast) <- "PROD"
# To a dataframe
prod.df <- as.data.frame(prod.rast, xy = T, na.rm = T)

# To a sf object - for some reason not working as a dataframe
prod.sf <- st_as_sf(prod.df, coords = c("x", "y"), crs = 4326)

# Gstat Nearest Neighbour interpolation
gs.prod <- gstat(id = "PROD", formula = PROD~1, data = prod.sf, 
                 nmax = 5, set = list(idp = 0))
prod.nn <- prod.rast %>% # raster::interpolate(object = prod.rast, model = gs.prod, na.rm = T, debug.level = 0)
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::crop(preds[[1]]) %>%
  raster::mask(preds[[1]])

# Cover the NAs of the original dataframe with the interpolated values
# prod.rast <- prod.rast %>%
#   raster::resample(preds[[1]], method = "bilinear") %>%
#   raster::cover(prod.nn) %>%
#   raster::crop(prod.nn) %>%
#   raster::mask(prod.nn)

# names(prod.rast) <- "PROD"
plot(prod.rast)

# Recreational fishing 
recfish <- st_read("data/spatial/shapefiles/NESP_D6_recboating_modelledPolygon.shp") %>%
  dplyr::select(geometry, meanTripsY) %>%
  st_rasterize() %>%                                                            # Stars                 
  rast()
recfish <- as(object = recfish, Class = "Raster")
names(recfish) <- "recfish"

# To a dataframe
recfish.df <- as.data.frame(recfish, xy = T, na.rm = T)

# To a sf object - for some reason not working as a dataframe
recfish.sf <- st_as_sf(recfish.df, coords = c("x", "y"), crs = 4326)

# Gstat Nearest Neighbour interpolation
gs.recfish <- gstat(id = "recfish", formula = recfish~1, data = recfish.sf, 
                    nmax = 5, set = list(idp = 0))
recfish.nn <- raster::interpolate(object = recfish, model = gs.recfish, na.rm = T, debug.level = 0) %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::crop(preds[[1]]) %>%
  raster::mask(preds[[1]])

# Cover the NAs of the original dataframe with the interpolated values
recfish <- recfish %>%
  raster::resample(preds[[1]], method = "bilinear") %>%
  raster::cover(recfish.nn) %>%
  raster::crop(recfish.nn) %>%
  raster::mask(recfish.nn)

names(recfish) <- "recfish"
plot(recfish)

# Join all the layers
preds <- stack(preds, sla.rast, ucur.rast, vcur.rast, sst.rast, recfish, prod.nn)
plot(preds)

preds <- rast(preds)
plot(preds)
# Save the output
preds <- terra::wrap(preds)
saveRDS(preds, paste(paste0('data/spatial/rasters/raw bathymetry/', name),      # This is ignored - too big!
                     'spatial_covariates.rds', sep = "_")) 