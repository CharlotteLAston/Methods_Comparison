##########################################
# Project: Methods Comparison 
# Data:    Ningaloo and Abrolhos BOSS/BRUV
# Task:    Deployment Maps
# Author:  Charlotte
# date:    September-October 2023
##########################################

# Libraries required
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(terra)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(forcats)
library(ggridges)
library(grid)
library(gridExtra)
library(gtable)
library(purrr)
library(matrixStats)
library(sfnetworks)
library(rcartocolor)
library(ggnewscale)
library(rgeos)
library(rnaturalearth)
library(ggpattern)
library(patchwork)
library(metR)

# Clear memory----
rm(list=ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")
sp_dir <- paste(working.dir, "spatial_data", sep="/")

a4.width <- 160

#### ABROLHOS DATA ####
#* Read in sampling points ####
study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

# Set CRS for transformations
wgscrs <- "+proj=longlat +datum=WGS84"
gdacrs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
sppcrs <- CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs")       # crs for sp objects

setwd(data_dir)

abrolhos.meta <-  readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  dplyr::select(sample, method, latitude, longitude) %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs=gdacrs)

# get and sort spatial boundaries
setwd(sp_dir)
aus    <- st_read("cstauscd_r.mif")                            # geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/
dirkh  <- aus[aus$ISLAND_NAME%in%c("DIRK HARTOG ISLAND", "DORRE ISLAND", "BERNIER ISLAND"), ]                        # just dirk hartog island
abro    <- aus[aus$GROUP_NAME %in% c("WALLABI GROUP/HOUTMAN ABROLHOS", "EASTERN IS/EASTER GP/HOUTMAN ABROLHOS",
                                     "EASTER GROUP/HOUTMAN ABROLHOS", "PELSAERT GROUP/HOUTMAN ABROLHOS",
                                     "MANGROVE GP/PELSAERT GP/HOUTMAN ABROLHOS"), ]
aus    <- aus[aus$FEAT_CODE == "mainland", ]
aumpa  <- st_read("AustraliaNetworkMarineParks.shp")           # all aus mpas
wampa  <- st_read("WA_MPA_2018.shp")                           # all wa mpas
ab_mpa <- wampa[wampa$NAME %in% c("Abrolhos Islands", #"Jurien Bay", "Ningaloo",
                                  "Hamelin Pool", "Shark Bay"), ]               # just wa parks nearby
sw_mpa <- aumpa[aumpa$NetName %in% c("South-west"), ]             # just W nat parks
nw_mpa <- aumpa[aumpa$NetName %in% c("North-west"), ] 
ab_nmp <- sw_mpa[sw_mpa$ResName %in% c("Abrolhos", "Jurien", "Shark Bay"), ]    # just nat parks nearby
cwatr  <- st_read('amb_coastal_waters_limit.shp')                    # coastal waters line trimmed in 'R/GA_coast_trim.R'


fisher.ceo <- st_read("Fisheries_Guide_CEO_Notices_Determinations_DPIRD_060.shp") %>% 
  filter(legal==6) %>% 
  mutate(waname = "Reef Observation Area") 


cbathy <- raster("bath_250_good.tif") 
bath_r <- rast(cbathy)
crs(bath_r) <- wgscrs
bath_r <- crop(bath_r, ext(109.4, 115.0607, -29.25, -24.4))
bath_df <- as.data.frame(bath_r, xy = T, na.rm = T)                             # Dataframe - cropped and above 0 use for bath cross section
bath_r <- clamp(bath_r, upper = 0, value = F)                               # Only data below 0
bathy <- as.data.frame(bath_r, xy = T, na.rm = T)

st_crs(aus)         <- st_crs(aumpa)
st_crs(dirkh)       <- st_crs(aumpa) 
st_crs(abro)       <- st_crs(aumpa) 



# simplify state parks names
ab_mpa$waname <- gsub("( \\().+(\\))", "", ab_mpa$ZONE_TYPE)
ab_mpa$waname <- gsub(" [1-4]", "", ab_mpa$waname)
# ab_mpa$waname[ab_mpa$ZONE_TYPE == unique(ab_mpa$ZONE_TYPE)[14]] <- 
#   c("Special Purpose Zone\n(Habitat Protection)")

ab_mpa$waname[ab_mpa$NAME == "Hamelin Pool"]     <- "Sanctuary Zone"
ab_mpa$waname[ab_mpa$NAME == "Abrolhos Islands"] <- "Fish Habitat Protection Area"

ab_mpa$waname <- dplyr::recode(ab_mpa$waname, 
                               "General Use" = "General Use Zone",
                               "Special Purpose Zone (Shore Based Activities)" = 
                                 "Special Purpose Zone\n(Shore Based Activities)",
                               "Special Purpose Zone (Seagrass Protection) (IUCN IV)" = "Special Purpose Zone",
)


ab_fpa <- ab_mpa %>% 
  filter(waname %in% "Fish Habitat Protection Area")

ab_mpa <- ab_mpa %>% 
  filter(waname != "Fish Habitat Protection Area")

# # reduce terrestrial parks
# terrnp <- terrnp[terrnp$leg_catego %in% c("Nature Reserve", "National Park"), ] # exclude state forests etc
# terrnp <- st_crop(terrnp, xmin = 113, ymin = -30, xmax = 116, ymax = -24)       # just abrolhos area
# plot(terrnp["leg_catego"])

#Key Ecological Features
# kef <- st_read("data/spatial/shp/AU_DOEE_KEF_2015.shp")
# sf_use_s2(F)                                                              
# kef <- st_crop(kef, c(xmin = 109, xmax = 116, ymin = -30, ymax = -24))  
# kef$NAME <- dplyr::recode(kef$NAME,"Perth Canyon and adjacent shelf break, and other west coast canyons" = "West coast canyons",                 
#                           "Commonwealth marine environment within and adjacent to the west coast inshore lagoons" = "West coast lagoons",                
#                           "Ancient coastline at 90-120m depth" = "Ancient coastline",                                                   
#                           "Western demersal slope and associated fish communities" = "Western demersal fish",                               
#                           "Western rock lobster" = "Western rock lobster",
#                           "Commonwealth marine environment surrounding the Houtman Abrolhos Islands" = "Abrolhos Islands")

# assign mpa colours - full levels are saved at end of script for future ref
nmpa_cols <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
                                          "Multiple Use Zone" = "#b9e6fb",
                                          "Special Purpose Zone" = "#6daff4",
                                          "Habitat Protection Zone" = "#fff8a3",
                                          "Recreational Use Zone" = "#b9e6fb"
))

wampa_cols <- scale_fill_manual(values = c("Fish Habitat Protection Area" = "#fac86b",
                                           "Reef Observation Area" = "#ddccff",
                                           "Sanctuary Zone" = "#bfd054",
                                           "General Use Zone" = "#bddde1",
                                           "Recreation Zone" = "#f4e952",
                                           "Special Purpose Zone" = "#c5bcc9"
                                           #"Marine Nature Reserve" = "#bfd054"
))

# WA terrestrial parks colours
# waterr_cols <- scale_fill_manual(values = c("National Park" = "#c4cea6",
#                                             "Nature Reserve" = "#e4d0bb"))

# # build basic plot elements
p1 <- ggplot() +
  geom_contour_filled(data = bathy, aes(x = x, y = y, z = bath_250_good,
                                        fill = after_stat(level)),
                      breaks = c(-30, -70, -200, - 700, -2000 , -4000,-6000)
                      ) +
  geom_contour(data = bathy, aes(x = x, y = y, z = bath_250_good),
               breaks = c(-30, -70, -200, - 700, -2000 , -4000,-6000),
               colour = "white", alpha = 3/5, size = 0.1) +
  scale_fill_grey(start = 1, end = 0.5, guide = "none") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = dirkh, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  geom_sf(data = ab_mpa, aes(fill = waname), alpha = 2/5, colour = NA) +
  # geom_sf(data = wampa, aes(fill = waname), alpha = 2/5, colour = NA) +
  wampa_cols +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  geom_sf(data = ab_nmp, aes(fill = ZoneName), alpha = 2/5, colour = NA) +
  nmpa_cols +
  labs(fill = "Australian Marine Parks") +
  new_scale_fill() +
  geom_sf(data = ab_fpa, aes(fill = waname), alpha = 2/5, colour = NA) +
  geom_sf(data = fisher.ceo, aes(fill = waname), alpha = 2/5, colour=NA)+
  wampa_cols +
  labs(fill = "Fisheries Protection Area") +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  labs(x = NULL, y = NULL, fill = "Fisheries Protection Area") +
  guides(fill = guide_legend(order = 1)) +
  annotate("rect", xmin = 113.02, xmax = 113.29, ymin = -27.19, ymax = -27.08,
           colour = "goldenrod1", fill = "white", alpha = 1/5, size = 0.4) +
  annotate("rect", xmin = 113.24, xmax = 113.58, ymin = -28.13, ymax = -28.02,
           colour = "goldenrod1", fill = "white", alpha = 1/5, size = 0.4) +
  annotate("point", y = c(-28.7761, -27.7115, -24.8838), x = c(114.6113, 114.1714, 113.6571), size = 0.75) +
  annotate("text", y = c(-28.7761, - 27.7115, -24.8838), x = c(114.95, 114.48, 114.075),
           label = c("Geraldton", "Kalbarri", "Carnarvon"), size = 3) +
  coord_sf(xlim = c(112, 115.061), ylim = c(-29.25, -24.4)) +
  theme_minimal()
p1

setwd(fig_dir)
ggsave(p1, filename=paste0("Abrolhos_site_map.png"),height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )


# site zoom plots

# reduce mpa levels for these plots
snmpa_cols <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
                                           "Multiple Use Zone" = "#b9e6fb",
                                           "Special Purpose Zone" = "#6daff4"
))

swampa_cols <- scale_fill_manual(values = c(
  "Fish Habitat Protection Area" = "#fbff85"
))

# closer plot
# sitebathy <- readRDS('output/ga_bathy_fine.rds')                                # finer bathy
# colnames(sitebathy)[3] <- "Depth"
# # sitebathy <- sitebathy[sitebathy$Depth > -1000, ]                               # trim to reduce legend
# sitebathy <- sitebathy[sitebathy$x > 112.7 & sitebathy$x < 115, ] 
# sitebathy <- sitebathy[sitebathy$y > -30 & sitebathy$y < -26.6, ]

p3 <- ggplot() +
  # geom_raster(data = sitebathy, aes(x, y, fill = Depth), alpha = 4/5) +
  geom_contour_filled(data = bathy, aes(x = x, y = y, z = bath_250_good,
                                            fill = after_stat(level)),
                      breaks = c(0, -30, -70, -200, -700, -2000, -4000, -10000), alpha = 4/5) +
  scale_fill_grey(start = 1, end = 0.5 , guide = "none") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  # geom_sf(data = terrnp, aes(fill = leg_catego), alpha = 4/5, colour = NA) +
  # labs(fill = "State Managed Areas") +
  # waterr_cols +
  # new_scale_fill() +
  geom_sf(data = ab_nmp, aes(fill = ZoneName), alpha = 3/5, colour = NA) +
  snmpa_cols +
  labs(x = NULL, y = NULL, fill = "Australian Marine Park") +
  geom_contour(data = bathy, aes(x = x, y = y, z = bath_250_good),
               breaks = c(0, -30, -70, -200, - 700, - 9000), colour = "white", alpha = 1, size = 0.2) +
  new_scale_colour()+
  geom_sf(data = abrolhos.meta, aes(colour = method),
             alpha = 3/5, shape = 10, size=1) +
  scale_colour_manual(values = c("#117733", "#88CCEE")) +
  new_scale_fill()+
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  annotate("rect", xmin = 113.02, xmax = 113.29, ymin = -27.19, ymax = -27.08,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("text", x = 113.15, y = -27.05, size = 2, 
           colour = "grey20", label = "Big Bank") +
  annotate("rect", xmin = 113.24, xmax = 113.58, ymin = -28.13, ymax = -28.02,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("text", x = 113.42, y = -27.99, size = 2,
           colour = "grey20", label = "Shallow Bank") +
  coord_sf(xlim = c(112.8, 115), ylim = c(-28.9, -26.7)) +
  annotate("text", y = c(-27.875,-27.875,-27.875,-27.875, -26.87), 
           x = c(114.07,113.32, 113.15, 112.79, 113.167), 
           label = c("30m", "70m", "200m", "700m", "70m"), size = 2) +
  annotate("point", y = c(-27.7115, -28.7761), x = c(114.1714, 114.6113), size = 0.5) +
  annotate("text", y = c(-27.7115, -28.7761), x = c(114.275 + 0.05, 114.6113 + 0.2),
           label = c("Kalbarri", "Geraldton"), size = 3) +
  # annotate(geom = "segment", linetype = "dashed", 
  #          x = c(113.2362, 112.8363, 113.7188), 
  #          xend = c(114.3337, 114.0763, 114.8787), 
  #          y = c(-28.07, -27.13, -28.91),                                       # Shallow Bank, Big Bank, Southern Group
  #          yend = c(-28.07, -27.13, -28.91), alpha = 0.5, size = 0.2) +
  # annotate(geom = "text", x = c(114.3337 + 0.05, 114.0763 + 0.05, 114.8787 + 0.05), 
  #          y = c(-28.07, -27.13, -28.91),
  #          label = c("b)", "a)", "c)"), size = 2) +
  labs(colour = "Sample", x = NULL, y = NULL) +
  theme_minimal()
p3

setwd(fig_dir)
ggsave(p3, filename=paste0("Abrolhos_larger_deployment.png"),height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


# ggsave("plots/spatial/siteplot.png", dpi = 200, width = 8, height = 6)

## single site zoom plots
snmpa_cols <- scale_colour_manual(values = c("National Park Zone" = "#7bbc63"))

nsitebathy <- bathy[bathy$bath_250_good > -160, ]                               # trim to reduce legend
nsitebathy <- nsitebathy[nsitebathy$x > 113 & nsitebathy$x < 113.3, ]
nsitebathy <- nsitebathy[nsitebathy$y > -27.2 & nsitebathy$y < -27.05, ]


p4 <- ggplot() +
  geom_raster(data = nsitebathy, aes(x, y, fill = bath_250_good), alpha = 4/5) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_contour(data = nsitebathy, aes(x = x, y = y, z = bath_250_good), 
               binwidth = 10, colour = "white", alpha = 1, size = 0.1) +
  geom_text_contour(data = nsitebathy, aes(x = x, y = y, z = bath_250_good), 
                    binwidth = 10, size = 2.5,
                    label.placer = label_placer_n(1)) +
  geom_sf(data = ab_nmp, aes(colour = ZoneName), alpha = 4/5, fill = NA, show.legend = F) +
  snmpa_cols + 
  labs(x = NULL, y = NULL, colour = NULL) +
  new_scale_colour() +
  geom_sf(data = abrolhos.meta, aes(colour = method),
          alpha = 3/5, shape = 10, size=3) +
  scale_colour_manual(values = c("#117733", "#88CCEE")) +
  # geom_point(data = bossd, aes(Longitude, Latitude, colour = "Drop Camera"), 
  #            alpha = 3/5, shape = 10) +
  # scale_colour_manual(values = c("BRUV" = "indianred4",
  #                                "Drop Camera" = "navyblue")) +
  coord_sf(xlim = c(113.02, 113.28), ylim = c(-27.18, -27.08)) +
  labs(colour = NULL, x = NULL, y = NULL) +
  theme_minimal()
p4

setwd(fig_dir)
ggsave(p4, filename="Abrolhos_big_bank_deployment.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


snmpa_cols <- scale_colour_manual(values = c("National Park Zone" = "#7bbc63"
                                             #  "Multiple Use Zone" = "#b9e6fb",
                                             # "Special Purpose Zone" = "#6daff4"
))
sab_nmp <- ab_nmp[ab_nmp$ZoneName == "National Park Zone", ]

ssitebathy <- bathy[bathy$bath_250_good < -10 & bathy$bath_250_good > -380, ]                               # trim to reduce legend
ssitebathy <- ssitebathy[ssitebathy$x > 113.23 & ssitebathy$x < 113.6, ]
ssitebathy <- ssitebathy[ssitebathy$y > -28.15 & ssitebathy$y < -28, ]

p5 <- ggplot() +
  geom_raster(data = ssitebathy, aes(x, y, fill = bath_250_good), alpha = 4/5) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_contour(data = ssitebathy, aes(x = x, y = y, z = bath_250_good), 
               binwidth = 20, colour = "white", alpha = 4/5, size = 0.1) +
  geom_text_contour(data = ssitebathy, aes(x = x, y = y, z = bath_250_good), 
                    binwidth = 20, size = 2.5, label.placer = label_placer_n(1)) +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = sab_nmp, aes(colour = ZoneName), alpha = 1, fill = NA, show.legend = F) +
  snmpa_cols + 
  labs(x = NULL, y = NULL, colour = NULL) +
  new_scale_colour() +
  geom_sf(data = abrolhos.meta, aes(colour = method),
          alpha = 3/5, shape = 10, size=3) +
  scale_colour_manual(values = c("#117733", "#88CCEE")) +
  coord_sf(xlim = c(113.24, 113.58), ylim = c(-28.125, -28.03)) +
  labs(colour = NULL, x = NULL, y = NULL) +
  theme_minimal()
p5

setwd(fig_dir)
ggsave(p5, filename="Abrolhos_shallow_bank_deployment.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### NINGALOO PLOTs ####
#* Read in sampling points ####
study <- "Ningaloo_PtCloates_BOSS-BRUV"
name <- study

setwd(data_dir)
ningaloo.meta <-  readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  dplyr::select(campaignid,sample, method, latitude, longitude) %>% 
  as.data.frame()# %>% 
  #st_as_sf(coords = c("longitude", "latitude"), crs=gdacrs)


setwd(sp_dir)
aus <- st_read("cstauscd_r.mif")   
aus <- aus[aus$FEAT_CODE == "mainland", ]
st_crs(aus) <- st_crs(aumpa)

cbathy <- raster("bath_250_good.tif") 
bath_r <- rast(cbathy)
crs(bath_r) <- wgscrs
bath_r <- crop(bath_r, ext(113.3, 114.5, -24, -21.6))
bath_df <- as.data.frame(bath_r, xy = T, na.rm = T)                             # Dataframe - cropped and above 0 use for bath cross section
bath_r <- clamp(bath_r, upper = 0, value = F)                               # Only data below 0
bathy <- as.data.frame(bath_r, xy = T, na.rm = T)

nw_mpa <- nw_mpa %>% 
  filter(ResName %in% "Ningaloo")

wampa <- st_read("WA_MPA_2020.shp")
st_crs(wampa) <- gdacrs
# Simplify names for plot legend
wampa$waname <- gsub("( \\().+(\\))", "", wampa$ZONE_TYPE)
wampa$waname <- gsub(" [1-4]", "", wampa$waname)
wampa$waname[wampa$NAME == "Hamelin Pool"]     <- "Marine Nature Reserve"
wampa$waname[wampa$NAME == "Abrolhos Islands"] <- "Fish Habitat Protection Area"
wampa$waname <- dplyr::recode(wampa$waname, 
                              "General Use" = "General Use Zone",
                              "Special Purpose Zone (Shore Based Activities)" = 
                                "Special Purpose Zone\n(Shore Based Activities)",
                              "Special Purpose Zone (Seagrass Protection) (IUCN IV)" = 
                                "Special Purpose Zone",
                              "MMA" = 'Marine Management Area' ) 

wampa_fills <- scale_fill_manual(values = c("Marine Management Area" = "#b7cfe1",
                                            "Conservation Area" = "#b3a63d",
                                            "Sanctuary Zone" = "#bfd054",
                                            "General Use Zone" = "#bddde1",
                                            "Recreation Area" = "#f4e952",
                                            "Special Purpose Zone" = "#c5bcc9",
                                            "Marine Nature Reserve" = "#bfd054"
),
name = "State Marine Parks")

ningaloo.site.plot <- ggplot() +
  # geom_raster(data = bathdf, aes(x, y, fill = Depth), alpha = 0.9) +
  # scale_fill_gradient(low = "black", high = "grey70") +
  # geom_contour_filled(data = bathy, aes(x = x, y = y, z = bath_250_good,
  #                                        fill = after_stat(level)),
  #                     breaks = c(0, -40, -70, -120, -7000)) +
  # geom_contour(data = bathy, aes(x = x, y = y, z = bath_250_good),
  # binwidth = 250, colour = "white", alpha = 3/5, size = 0.1) +
  # scale_fill_grey(start = 1, end = 0.5, guide = "none") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  #geom_sf(data = nw_mpa, aes(fill = waname), alpha = 2/5, colour = NA) +
  geom_sf(data = wampa, aes(fill = waname), alpha = 2/5, colour = NA) +
  wampa_fills +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  # geom_sf(data = terrnp, aes(fill = leg_catego), alpha = 4/5, colour = NA) +
  # labs(fill = "State Managed Areas") +
  # waterr_cols +
  # new_scale_fill() +
  geom_sf(data = nw_mpa, aes(fill = ZoneName), alpha = 4/5, colour = NA) +
  nmpa_cols +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  labs(x = NULL, y = NULL, fill = "Australian Marine Parks") +
  guides(fill = guide_legend(order = 1)) +
  annotate("rect", xmin = min(ningaloo.meta$longitude), xmax = max(ningaloo.meta$longitude),
           ymin = min(ningaloo.meta$latitude), ymax = max(ningaloo.meta$latitude),
           colour = "goldenrod1", fill = "white", alpha = 0.2, size = 0.6) +
  coord_sf(xlim = c(113.3, 114.5), ylim = c(-23.5, -21.6)) +
  annotate(geom = "text", x = c((114.1279 + 0.12), (113.6775 + 0.125)), 
           y = c(-21.9323, (-22.7212+0.03)), label = c("Exmouth", "Pt Cloates"),
           size = 3) +
  annotate(geom = "point", x = c(114.1279, 113.6775), 
           y = c(-21.9323, -22.7212)) +
  theme_minimal()
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank())
ningaloo.site.plot

inset.map <- ggplot() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey90", size = 0.075) +
  #geom_sf(data = aumpa, alpha = 5/6, colour = "grey85", size = 0.02) +
  coord_sf(xlim = c(108, 125), ylim = c(-37, -13)) +
  annotate(geom = "text", x=c(110.2), y=c(-29.4), label = c("Indian\nOcean"), size=3)+
  annotate(geom = "text", x=c(120.8), y=c(-25.94), label = c("Western\nAustralia"), size=3.5)+
  annotate("rect", xmin = 113, xmax = 114.35, ymin = -22.8, ymax = -21.5,   # Change here 
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("rect", xmin = 113, xmax = 115.061, ymin = -29.25, ymax = -28, 
            colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey70"))+
  ylab(NULL)+
  xlab(NULL)
inset.map

ningaloo.site.plot <- ningaloo.site.plot + inset_element(inset.map, left = -1.29, right = 2.874, top = 0.4, bottom = 0)  

setwd(fig_dir)
ggsave(ningaloo.site.plot, filename=paste0("Ningaloo_site_map.png"),height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )

#* Ningaloo sample plots ####
ningaloo.meta.point <-  ningaloo.meta %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs=gdacrs) 

ningaloo.sample.big <- ggplot() +
  # geom_raster(data = sitebathy, aes(x, y, fill = Depth), alpha = 4/5) +
  geom_contour_filled(data = bathy, aes(x = x, y = y, z = bath_250_good,
                                        fill = after_stat(level)),
                      breaks = c(0, -30, -70, -200, -700, -2000, -4000, -10000), alpha = 4/5) +
  scale_fill_grey(start = 1, end = 0.5 , guide = "none") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  geom_sf(data = wampa, aes(fill = waname), alpha = 2/5, colour = NA) +
  wampa_fills +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  geom_sf(data = nw_mpa, aes(fill = ZoneName), alpha = 4/5, colour = NA) +
  nmpa_cols +
  labs(x = NULL, y = NULL, fill = "Australian Marine Park") +
  geom_contour(data = bathy, aes(x = x, y = y, z = bath_250_good),
               breaks = c(0, -30, -70, -200, - 700, - 9000), colour = "white", alpha = 1, size = 0.2) +
  new_scale_colour()+
  geom_sf(data = ningaloo.meta.point, aes(colour = method),
          alpha = 3/5, shape = 10, size=1) +
  scale_colour_manual(values = c("#117733", "#88CCEE")) +
  new_scale_fill()+
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  coord_sf(xlim = c(113.3, 114), ylim = c(-23, -22.5)) +
  annotate(geom = "text", x = c((114.1279 + 0.12), (113.6775 + 0.035)), 
           y = c(-21.9323, (-22.7212+0.015)), label = c("Exmouth", "Pt Cloates"),
           size = 3) +
  annotate(geom = "point", x = c(114.1279, 113.6775), 
           y = c(-21.9323, -22.7212)) +
  labs(colour = "Sample", x = NULL, y = NULL) +
  theme_minimal()
ningaloo.sample.big

setwd(fig_dir)
ggsave(ningaloo.sample.big, filename=paste0("Ningaloo_site_samples.png"),height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )

#* 2021 Samples ####

sitebathy <- bathy[bathy$bath_250_good > -160, ]                               # trim to reduce legend
sitebathy <- sitebathy[sitebathy$x > 113.3 & sitebathy$x < 114, ]
sitebathy <- sitebathy[sitebathy$y > -23 & sitebathy$y < -22.5, ]

ningaloo.meta.2021 <- ningaloo.meta.point %>% 
  filter(campaignid %in% c("2021-08_Pt-Cloates_stereo-BRUVs", "2021-08_PtCloates_BOSS" )) %>% 
  distinct(sample, .keep_all = T)

sample.2021 <- ggplot() +
  geom_raster(data = sitebathy, aes(x, y, fill = bath_250_good), alpha = 4/5) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_contour(data = sitebathy, aes(x = x, y = y, z = bath_250_good), 
               binwidth = 10, colour = "white", alpha = 1, size = 0.1) +
  geom_text_contour(data = nsitebathy, aes(x = x, y = y, z = bath_250_good), 
                    binwidth = 10, size = 2.5,
                    label.placer = label_placer_n(1)) +
  new_scale_fill()+
  geom_sf(data = wampa, aes(fill = waname), alpha = 2/5, colour = NA) +
  wampa_fills +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  geom_sf(data = nw_mpa, aes(fill = ZoneName), alpha = 4/5, colour = NA) +
  nmpa_cols +
  labs(x = NULL, y = NULL, colour = NULL) +
  new_scale_colour() +
  geom_sf(data = ningaloo.meta.2021, aes(colour = method),
          alpha = 3/5, shape = 10, size=3) +
  scale_colour_manual(values = c("#117733", "#88CCEE")) +
  # geom_point(data = bossd, aes(Longitude, Latitude, colour = "Drop Camera"), 
  #            alpha = 3/5, shape = 10) +
  # scale_colour_manual(values = c("BRUV" = "indianred4",
  #                                "Drop Camera" = "navyblue")) +
  coord_sf(xlim = c(113.5, 113.65), ylim = c(-22.775, -22.66)) +
  labs(colour = NULL, x = NULL, y = NULL) +
  theme_minimal()
sample.2021

setwd(fig_dir)
ggsave(sample.2021, filename=paste0("Ningaloo_2021_samples.png"),height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )


#* 2022 Samples ####
ningaloo.meta.2022 <- ningaloo.meta.point %>% 
  filter(campaignid %in% c("2022-05_PtCloates_stereo-BRUVS", "2022-05_PtCloates_BOSS" )) %>% 
  distinct(sample, .keep_all = T)

sample.2022 <- ggplot() +
  geom_raster(data = sitebathy, aes(x, y, fill = bath_250_good), alpha = 4/5) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_contour(data = sitebathy, aes(x = x, y = y, z = bath_250_good), 
               binwidth = 10, colour = "white", alpha = 1, size = 0.1) +
  geom_text_contour(data = nsitebathy, aes(x = x, y = y, z = bath_250_good), 
                    binwidth = 10, size = 2.5,
                    label.placer = label_placer_n(1)) +
  new_scale_fill()+
  geom_sf(data = wampa, aes(fill = waname), alpha = 2/5, colour = NA) +
  wampa_fills +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  geom_sf(data = nw_mpa, aes(fill = ZoneName), alpha = 2/5, colour = NA) +
  nmpa_cols +
  labs(x = NULL, y = NULL, colour = NULL) +
  new_scale_colour() +
  geom_sf(data = ningaloo.meta.2022, aes(colour = method),
          alpha = 4/5, shape = 10, size=3) +
  scale_colour_manual(values = c("#117733", "#88CCEE")) +
  coord_sf(xlim = c(113.5, 113.65), ylim = c(-22.775, -22.66)) +
  labs(colour = NULL, x = NULL, y = NULL) +
  theme_minimal()
sample.2022

setwd(fig_dir)
ggsave(sample.2022, filename=paste0("Ningaloo_2022_samples.png"),height = a4.width*1, width = a4.width*1.5, units  ="mm", dpi = 300 )
