###
# Project: Parks Abrolhos
# Data:    Geoscience Australia 250m res bathy
# Task:    Bathymetry cross section figure
# author:  Claude Spencer
# date:    September 2022
##

# Clear your environment
rm(list = ls())

# Load libraries
library(dplyr)
library(sf)
library(rgeos)
library(rnaturalearth)
library(ggplot2)
library(metR)
library(stringr)
library(patchwork)
library(terra)
library(ggnewscale)
library(GlobalArchive)
library(tidyverse)
library(viridis)

# Set your study name
name <- "Abrolhos"                                                              # Change here

# Set CRS for transformations
wgscrs <- "+proj=longlat +datum=WGS84"
gdacrs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

# Australian outline and state and commonwealth marine parks
aus    <- st_read("data/spatial/shp/cstauscd_r.mif") %>%                        # Geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/
  dplyr::filter(FEAT_CODE %in% c("mainland", "island"))
st_crs(aus) <- gdacrs

# read in and merge GA coarse bathy tiles from https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/67703
cbaths <- list.files("data/spatial/raster", "*tile", full.names = TRUE)
cbathy <- lapply(cbaths, function(x){read.table(file = x, header = TRUE, sep = ",")})
cbathy <- do.call("rbind", lapply(cbathy, as.data.frame))

# 1. Bathymetry cross section
# Shallow bank = 28.07
sf_use_s2(T)
auss <- st_transform(aus, wgscrs)
auss <- auss[auss$FEAT_CODE %in% "mainland", ]
auss <- st_union(auss)
ausout <- st_cast(auss, "MULTILINESTRING")

bath_cross <- cbathy %>%
  dplyr::filter(abs(Y - -28.07) == min(abs(Y - -28.07)),
                 Z > -220) %>% 
  st_as_sf(coords = c("X", "Y"), crs = wgscrs)

bath_sf <- bath_cross %>%
  dplyr::mutate("distance.from.coast" = st_distance(bath_cross, ausout),
                x = unlist(map(bath_cross$geometry, 1)),
                y = unlist(map(bath_cross$geometry, 2)),
                land = lengths(st_intersects(bath_cross, auss)) > 0) %>%
  glimpse()

bath_slice1 <- as.data.frame(bath_sf) %>%
  dplyr::select(-geometry) %>%
  dplyr::rename(depth = "Z") %>%
  dplyr::mutate(distance.from.coast = as.numeric(distance.from.coast/1000)) %>%
  dplyr::mutate(distance.from.coast = ifelse(land %in% "FALSE", distance.from.coast*-1, distance.from.coast)) %>%
  glimpse()

paleo <- data.frame(depth = c(-118, -94, -63, -41),
                    label = c("18-20 Ka", "15-17 Ka", "12-13 Ka", "9-10 Ka"))

for (i in 1:nrow(paleo)) {
  temp <- bath_slice1 %>%
    dplyr::filter(abs(bath_slice1$depth - paleo$depth[i]) == min(abs(bath_slice1$depth - paleo$depth[i]))) %>%
    dplyr::select(depth, distance.from.coast) %>%
    slice(1)
  
  if (i == 1) {
    dat <- temp
  } 
  else {
    dat <- bind_rows(dat, temp)
  }
}

paleo$distance.from.coast <- dat$distance.from.coast
rm("temp", "dat")

p1 <- ggplot() +
  geom_rect(aes(xmin = min(bath_slice1$distance.from.coast), xmax = 9, ymin =-Inf, ymax = 0), fill = "#12a5db", alpha = 0.5) +
  annotate("segment", x = -5.556, xend = -5.556, y = -29, yend = 0, colour = "red") +
  annotate("segment", x = -57.6, xend = -57.6, y = -49, yend = 0, colour = "#7bbc63") + # Inside
  annotate("segment", x = -88, xend = -88, y = -210, yend = 0, colour = "#7bbc63") + # Outside
  geom_line(data = bath_slice1, aes(y = depth, x = distance.from.coast)) +
  geom_ribbon(data = bath_slice1, aes(ymin = -Inf, ymax = depth, x = distance.from.coast), fill = "tan") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(min(bath_slice1$distance.from.coast), 15)) +
  labs(x = "Distance from coast (km)", y = "Elevation (m)", title = "b) Shallow Bank") +
  geom_segment(data = paleo, aes(x = distance.from.coast, xend = distance.from.coast + 10, 
                                 y = depth, yend = depth), linetype = 2, alpha = 0.5) +
  geom_text(data = paleo, aes(x = distance.from.coast + 13, y = depth, label = label), size = 2)
p1

# Big Bank = 27.13
bath_cross <- cbathy %>%
  dplyr::filter(abs(Y - -27.13) == min(abs(Y - -27.13)),
                Z > -220) %>% 
  st_as_sf(coords = c("X", "Y"), crs = wgscrs)

bath_sf <- bath_cross %>%
  dplyr::mutate("distance.from.coast" = st_distance(bath_cross, ausout),
                x = unlist(map(bath_cross$geometry, 1)),
                y = unlist(map(bath_cross$geometry, 2)),
                land = lengths(st_intersects(bath_cross, auss)) > 0) %>%
  glimpse()

bath_slice2 <- as.data.frame(bath_sf) %>%
  dplyr::select(-geometry) %>%
  dplyr::rename(depth = "Z") %>%
  dplyr::mutate(distance.from.coast = as.numeric(distance.from.coast/1000)) %>%
  dplyr::mutate(distance.from.coast = ifelse(land %in% "FALSE", distance.from.coast*-1, distance.from.coast)) %>%
  glimpse()

paleo <- data.frame(depth = c(-118, -94, -63, -41),
                    label = c("18-20 Ka", "15-17 Ka", "12-13 Ka", "9-10 Ka"))

for (i in 1:nrow(paleo)) {
  temp <- bath_slice2 %>%
    dplyr::filter(abs(bath_slice2$depth - paleo$depth[i]) == min(abs(bath_slice2$depth - paleo$depth[i]))) %>%
    dplyr::select(depth, distance.from.coast) %>%
    slice(1)
  
  if (i == 1) {
    dat <- temp
  } 
  else {
    dat <- bind_rows(dat, temp)
  }
}

paleo$distance.from.coast <- dat$distance.from.coast
rm("temp", "dat")

p2 <- ggplot() +
  geom_rect(aes(xmin = min(bath_slice2$distance.from.coast), 
                xmax = 9, ymin =-Inf, ymax = 0), 
            fill = "#12a5db", alpha = 0.5) +
  annotate("segment", x = -5.556, xend = -5.556, y = -66, yend = 0, colour = "red") +
  annotate("segment", x = -52.2, xend = -52.2, y = -104, yend = 0, colour = "#7bbc63") + # Inside
  annotate("segment", x = -73.5, xend = -73.5, y = -143, yend = 0, colour = "#7bbc63") + # Outside
  geom_line(data = bath_slice2, aes(y = depth, x = distance.from.coast)) +
  geom_ribbon(data = bath_slice2, 
              aes(ymin = -Inf, ymax = depth, x = distance.from.coast), 
              fill = "tan") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(min(bath_slice2$distance.from.coast), 15)) +
  labs(x = "Distance from coast (km)", y = "Elevation (m)", title = "a) Big Bank") +
  geom_segment(data = paleo, aes(x = distance.from.coast, xend = distance.from.coast + 10, 
                                 y = depth, yend = depth), linetype = 2, alpha = 0.5) +
  geom_text(data = paleo, aes(x = distance.from.coast + 13, y = depth, label = label), size = 2)
p2

# Southern group = 28.91 - doesnt really go through gero!
bath_cross <- cbathy %>%
  dplyr::filter(abs(Y - -28.91) == min(abs(Y - -28.91)),
                Z > -250) %>% 
  st_as_sf(coords = c("X", "Y"), crs = wgscrs)

bath_sf <- bath_cross %>%
  dplyr::mutate("distance.from.coast" = st_distance(bath_cross, ausout),
                x = unlist(map(bath_cross$geometry, 1)),
                y = unlist(map(bath_cross$geometry, 2)),
                land = lengths(st_intersects(bath_cross, auss)) > 0) %>%
  glimpse()

bath_slice3 <- as.data.frame(bath_sf) %>%
  dplyr::select(-geometry) %>%
  dplyr::rename(depth = "Z") %>%
  dplyr::mutate(distance.from.coast = as.numeric(distance.from.coast/1000)) %>%
  dplyr::mutate(distance.from.coast = ifelse(land %in% "FALSE", distance.from.coast*-1, distance.from.coast)) %>%
  glimpse()

paleo <- data.frame(depth = c(-118, -94, -63, -45),
                    label = c("18-20 Ka", "15-17 Ka", "12-13 Ka", "9-10 Ka"))

for (i in 1:nrow(paleo)) {
  temp <- bath_slice3 %>%
    dplyr::filter(abs(bath_slice3$depth - paleo$depth[i]) == min(abs(bath_slice3$depth - paleo$depth[i]))) %>%
    dplyr::select(depth, distance.from.coast) %>%
    arrange(distance.from.coast) %>%
    slice(1)
  
  if (i == 1) {
    dat <- temp
  } 
  else {
    dat <- bind_rows(dat, temp)
  }
}

paleo$distance.from.coast <- dat$distance.from.coast
rm("temp", "dat")

p3 <- ggplot() +
  geom_rect(aes(xmin = min(bath_slice3$distance.from.coast), xmax = 9, ymin =-Inf, ymax = 0), fill = "#12a5db", alpha = 0.5) +
  geom_line(data = bath_slice3, aes(y = depth, x = distance.from.coast)) +
  geom_ribbon(data = bath_slice3, aes(ymin = -Inf, ymax = depth, x = distance.from.coast), fill = "tan") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(min(bath_slice3$distance.from.coast), 15)) +
  labs(x = "Distance from coast (km)", y = "Elevation (m)", title = "c) Southern Group") +
  geom_segment(data = paleo, aes(x = distance.from.coast, xend = distance.from.coast + 10, 
                                 y = depth, yend = depth), linetype = 2, alpha = 0.5) +
  annotate("text", x = -70, y = 30, label = "Southern Group", size = 2.5) +
  annotate("text", x = 0, y = 80, label = "Geraldton", size = 2.5) +
  annotate("segment", x = -5.556, xend = -5.556, y = -21, yend = 0, colour = "red") +
  annotate("segment", x = -79.13277, xend = -79.13277, y = -69, yend = 0, colour = "red") +
  annotate("segment", x = -52.07578, xend = -52.07578, y = -39, yend = 0, colour = "red") +
  geom_text(data = paleo, aes(x = distance.from.coast + 13, y = depth, label = label), size = 2)
p3

plots <- p2 / p1 / p3
plots

png(filename = paste(paste0('plots/spatial/', name) , 'bathymetry-cross-section.png', 
                     sep = "-"), units = "in", res = 200, width = 8, height = 10)
plots
dev.off()
