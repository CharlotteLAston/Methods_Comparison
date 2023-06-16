###
# Project: Parks Abrolhos
# Data:    BRUVS, BOSS
# Task:    Scatterpies
# author:  Claude
# date:    July 2022
##

rm(list=ls())

library(dplyr)
library(ggplot2)
library(scatterpie)
library(viridis)
library(sf)
library(raster)
library(ggnewscale)
library(metR)
library(cowplot)

working.dir <- getwd()
setwd(working.dir)

#define crs
wgscrs <- CRS("+proj=longlat +datum=WGS84")
sppcrs <- CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs")     # crs for sp objects

# bring in marine parks
aumpa  <- st_read("data/spatial/shp/AustraliaNetworkMarineParks.shp")%>%           # all aus mpas
  dplyr::filter(ResName%in%"Abrolhos",ZoneIUCN%in%"II")
#bring in state MP
wampa  <- st_read("data/spatial/shp/WA_MPA_2018.shp")                           # all wa mpas
ab_mpa <- wampa[wampa$NAME %in% c("Abrolhos Islands", #"Jurien Bay", "Ningaloo",
                                  "Hamelin Pool", "Shark Bay"), ]               # just wa parks nearby
# simplify state parks names
ab_mpa$waname <- gsub("( \\().+(\\))", "", ab_mpa$ZONE_TYPE)
ab_mpa$waname <- gsub(" [1-4]", "", ab_mpa$waname)
ab_mpa$waname[ab_mpa$NAME == "Hamelin Pool"]     <- "Marine Nature Reserve"
ab_mpa$waname[ab_mpa$NAME == "Abrolhos Islands"] <- "Fish Habitat Protection Area"

ab_mpa$waname <- dplyr::recode(ab_mpa$waname, 
                               "General Use" = "General Use Zone",
                               "Special Purpose Zone (Shore Based Activities)" = 
                                 "Special Purpose Zone\n(Shore Based Activities)",
                               "Special Purpose Zone (Seagrass Protection) (IUCN IV)" = "Special Purpose Zone",
)

# get aus outline data
aus    <- st_read("data/spatial/shp/cstauscd_r.mif")                            # geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/
aus    <- aus[aus$FEAT_CODE == "mainland", ]
st_crs(aus)         <- st_crs(aumpa)

cwatr  <- readRDS('output/coastal_waters_limit_trimmed.rds')                    # coastal waters line trimmed in 'R/GA_coast_trim.R'

dat <- readRDS("data/Tidy/merged_habitat.rds") %>%
  dplyr::rename("Sessile invertebrates" = biog,
                "Rock" = rock,
                "Macroalgae" = macroalgae,
                "Sand" = sand,
                "Kelp" = kelps) %>%
  dplyr::mutate(grouping = factor(1:125)) %>%
  glimpse()

#bring in bathy for contour lines
bathdf <- readRDS("output/ga_bathy_fine.rds")                                   # bathymetry trimmed in 'R/GA_coast_trim.R'
colnames(bathdf)[3] <- "Depth"
# assign commonwealth zone colours
nmpa_cols <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
                                          "Habitat Protection Zone" = "#fff8a3",
                                          "Multiple Use Zone" = "#b9e6fb",
                                          "Special Purpose Zone\n(Mining Exclusion)" = "#368ac1")) 

# state colours
wampa_cols <- scale_fill_manual(values = c("Fish Habitat Protection Area" = "#fac86b",
                                           "Reef Observation Area" = "#ddccff",
                                           "Sanctuary Zone" = "#bfd054",
                                           "General Use Zone" = "#bddde1",
                                           "Recreation Zone" = "#f4e952",
                                           "Special Purpose Zone" = "#c5bcc9",
                                           "Marine Nature Reserve" = "#bfd054"
))

#class colours 
hab_cols <- scale_fill_manual(values = c("Rock" = "grey40",
                                         "Sessile invertebrates" = "plum",
                                         "Macroalgae" = "darkgoldenrod4",
                                         "Kelp" = "goldenrod1",
                                         "Sand" = "wheat"))

#depth colours 
depth_cols <- scale_fill_manual(values = c("#a7cfe0","#9acbec","#98c4f7", "#a3bbff"),guide = "none")
#shallow to deep

#make the plot
test <- dat%>%dplyr::filter(location%in%"NPZ6")
min(test$longitude)
max(test$longitude)
min(test$latitude)
max(test$latitude)

gg.scatterpie.npz6 <- ggplot() + 
  geom_contour_filled(data = bathdf, aes(x, y, z = Depth, fill = after_stat(level)), color = "black",
                      breaks = c(-30, -70, -200,-700,-10000), size = 0.1) +
  annotate("text", x = c(114.40,114.467,114.72,114.945), y = -33.85, label = c("700m","200m","70m","30m"), size = 2)+
  depth_cols+
  new_scale_fill()+
  geom_sf(data = aus, fill = "seashell2", colour = "black", size = 0.1) +
  geom_sf(data = wampa,fill = "#bfd054", alpha = 2/5, color = NA)+
  wampa_cols+
  labs(fill = "State Marine Parks")+
  new_scale_fill()+
  geom_sf(data = aumpa, fill = "#7bbc63",alpha = 2/5, color = NA) +
  labs(fill = "Australian Marine Parks")+
  nmpa_cols+
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.3) +
  new_scale_fill()+
  geom_scatterpie(aes(x=longitude, y=latitude, group=grouping), data=dat,
                  cols = c("Rock","Sessile invertebrates","Macroalgae",
                           "Kelp", "Sand"),
                  pie_scale = 0.45, color = NA) +
  labs(fill = "Habitat",x = 'Longitude', y = 'Latitude', title = "Shallow Bank")+
  hab_cols + 
  annotate("text", x = c(113.47, 113.405, 113.278), y = c(-28.13, -28.13, -28.13), label = c("30m", "70m", "200m"),
           size = 1.5, colour = "black")+
  coord_sf(xlim = c(113.169637818, 113.592952023), ylim = c(-28.147530871, -27.951387524))+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "#b9d1d6"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

gg.scatterpie.npz9 <- ggplot() + 
  geom_contour_filled(data = bathdf, aes(x, y, z = Depth, fill = after_stat(level)), color = "black",
                      breaks = c(-30, -70, -200,-700,-10000), size = 0.1) +
  annotate("text", x = c(114.40,114.467,114.72,114.945), y = -33.85, label = c("700m","200m","70m","30m"), size = 2)+
  depth_cols+
  new_scale_fill()+
  geom_sf(data = aus, fill = "seashell2", colour = "black", size = 0.1) +
  geom_sf(data = wampa,fill = "#bfd054", alpha = 2/5, color = NA)+
  wampa_cols+
  labs(fill = "State Marine Parks")+
  new_scale_fill()+
  geom_sf(data = aumpa, fill = "#7bbc63",alpha = 2/5, color = NA) +
  labs(fill = "Australian Marine Parks")+
  nmpa_cols+
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.3) +
  new_scale_fill()+
  geom_scatterpie(aes(x=longitude, y=latitude, group=grouping), data=dat,
                  cols = c("Rock","Sessile invertebrates","Macroalgae",
                           "Kelp", "Sand"),
                  pie_scale = 0.45, color = NA) +
  labs(fill = "Habitat",x = 'Longitude', y = 'Latitude', title = "Big Bank")+
  hab_cols+
  coord_sf(xlim = c(113.016535099, 113.316210677), ylim = c(-27.249154771, -27.077529621))+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "#b9d1d6"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

gg.scatterpie <- gg.scatterpie.npz9 / gg.scatterpie.npz6 + plot_layout(guide = "collect")

save_plot("plots/habitat/scatterpies.png", gg.scatterpie,base_height = 6.5,base_width = 7)

