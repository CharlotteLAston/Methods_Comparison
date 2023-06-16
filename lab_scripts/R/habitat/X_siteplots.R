
###
# Project: Parks - Abrolhos Post-Survey
# Data:    BRUVS, BOSS
# Task:    Overview maps
# author:  Kingsley Griffin & Claude Spencer
# date:    Jun 2021
##

rm(list = ls())

library(sf)
library(rgeos)
library(rnaturalearth)
library(ggplot2)
library(metR)
library(googlesheets4)
library(patchwork)
library(raster)
library(ggnewscale)
library(maptools)

# get and sort spatial boundaries
aus    <- st_read("data/spatial/shp/cstauscd_r.mif")                            # geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/
dirkh  <- aus[aus$ISLAND_NAME%in%c("DIRK HARTOG ISLAND", "DORRE ISLAND", "BERNIER ISLAND"), ]                        # just dirk hartog island
abro    <- aus[aus$GROUP_NAME %in% c("WALLABI GROUP/HOUTMAN ABROLHOS", "EASTERN IS/EASTER GP/HOUTMAN ABROLHOS",
                                   "EASTER GROUP/HOUTMAN ABROLHOS", "PELSAERT GROUP/HOUTMAN ABROLHOS",
                                   "MANGROVE GP/PELSAERT GP/HOUTMAN ABROLHOS"), ]
aus    <- aus[aus$FEAT_CODE == "mainland", ]
aumpa  <- st_read("data/spatial/shp/AustraliaNetworkMarineParks.shp")           # all aus mpas
wampa  <- st_read("data/spatial/shp/WA_MPA_2018.shp")                           # all wa mpas
ab_mpa <- wampa[wampa$NAME %in% c("Abrolhos Islands", #"Jurien Bay", "Ningaloo",
                                  "Hamelin Pool", "Shark Bay"), ]               # just wa parks nearby
sw_mpa <- aumpa[aumpa$NetName %in% c("South-west", "North-west"), ]             # just W nat parks
ab_nmp <- sw_mpa[sw_mpa$ResName %in% c("Abrolhos", "Jurien", "Shark Bay"), ]    # just nat parks nearby
cwatr  <- readRDS('output/coastal_waters_limit_trimmed.rds')                    # coastal waters line trimmed in 'R/GA_coast_trim.R'
bathdf <- readRDS("output/ga_bathy_trim.rds")                                   # bathymetry trimmed in 'R/GA_coast_trim.R'
colnames(bathdf)[3] <- "Depth"
terrnp <- st_read("data/spatial/shp/Legislated_Lands_and_Waters_DBCA_011.shp") %>%  # terrestrial reserves
  dplyr::filter(leg_catego %in% c("Nature Reserve", "National Park"))
st_crs(aus)         <- st_crs(aumpa)
st_crs(dirkh)       <- st_crs(aumpa) 
st_crs(abro)       <- st_crs(aumpa) 

roas <- st_read("data/spatial/shp/Abrolhos_ROAs.shp")                           # from matt's state reserve shapefile 

# simplify state parks names
ab_mpa$waname <- gsub("( \\().+(\\))", "", ab_mpa$ZONE_TYPE)
ab_mpa$waname <- gsub(" [1-4]", "", ab_mpa$waname)
# ab_mpa$waname[ab_mpa$ZONE_TYPE == unique(ab_mpa$ZONE_TYPE)[14]] <- 
#   c("Special Purpose Zone\n(Habitat Protection)")

ab_mpa$waname[ab_mpa$NAME == "Hamelin Pool"]     <- "Marine Nature Reserve"
ab_mpa$waname[ab_mpa$NAME == "Abrolhos Islands"] <- "Fish Habitat Protection Area"

ab_mpa$waname <- dplyr::recode(ab_mpa$waname, 
                               "General Use" = "General Use Zone",
                               "Special Purpose Zone (Shore Based Activities)" = 
                                 "Special Purpose Zone\n(Shore Based Activities)",
                               "Special Purpose Zone (Seagrass Protection) (IUCN IV)" = "Special Purpose Zone",
                               )

# reduce terrestrial parks
terrnp <- terrnp[terrnp$leg_catego %in% c("Nature Reserve", "National Park"), ] # exclude state forests etc
terrnp <- st_crop(terrnp, xmin = 113, ymin = -30, xmax = 116, ymax = -24)       # just abrolhos area
plot(terrnp["leg_catego"])

#Key Ecological Features
kef <- st_read("data/spatial/shp/AU_DOEE_KEF_2015.shp")
sf_use_s2(F)                                                              
kef <- st_crop(kef, c(xmin = 109, xmax = 116, ymin = -30, ymax = -24))  
kef$NAME <- dplyr::recode(kef$NAME,"Perth Canyon and adjacent shelf break, and other west coast canyons" = "West coast canyons",                 
                          "Commonwealth marine environment within and adjacent to the west coast inshore lagoons" = "West coast lagoons",                
                          "Ancient coastline at 90-120m depth" = "Ancient coastline",                                                   
                          "Western demersal slope and associated fish communities" = "Western demersal fish",                               
                          "Western rock lobster" = "Western rock lobster",
                          "Commonwealth marine environment surrounding the Houtman Abrolhos Islands" = "Abrolhos Islands")

# assign mpa colours - full levels are saved at end of script for future ref
nmpa_cols <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
                                          "Multiple Use Zone" = "#b9e6fb",
                                          "Special Purpose Zone" = "#6daff4",
                                          "Habitat Protection Zone" = "#fff8a3"
))

wampa_cols <- scale_fill_manual(values = c("Fish Habitat Protection Area" = "#fac86b",
                                           "Reef Observation Area" = "#ddccff",
                                           "Sanctuary Zone" = "#bfd054",
                                           "General Use Zone" = "#bddde1",
                                           "Recreation Zone" = "#f4e952",
                                           "Special Purpose Zone" = "#c5bcc9",
                                           "Marine Nature Reserve" = "#bfd054"
))

# WA terrestrial parks colours
waterr_cols <- scale_fill_manual(values = c("National Park" = "#c4cea6",
                                            "Nature Reserve" = "#e4d0bb"))

# build basic plot elements
p1 <- ggplot() +
  geom_contour_filled(data = bathdf, aes(x = x, y = y, z = Depth,
                                         fill = after_stat(level)),
                      breaks = c(-30, -70, -200, - 700, -2000 , -4000,-6000)) +
  geom_contour(data = bathdf, aes(x = x, y = y, z = Depth),
  breaks = c(-30, -70, -200, - 700, -2000 , -4000,-6000), colour = "white", alpha = 3/5, size = 0.1) +
  scale_fill_grey(start = 1, end = 0.5, guide = "none") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = dirkh, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  geom_sf(data = ab_mpa, aes(fill = waname), alpha = 2/5, colour = NA) +
  geom_sf(data = roas, aes(fill = Type), alpha = 2/5, colour = NA) +
  wampa_cols +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  geom_sf(data = terrnp, aes(fill = leg_catego), alpha = 4/5, colour = NA) +
  labs(fill = "State Managed Areas") +
  waterr_cols +
  new_scale_fill() +
  geom_sf(data = ab_nmp, aes(fill = ZoneName), alpha = 4/5, colour = NA) +
  nmpa_cols +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  labs(x = NULL, y = NULL, fill = "Australian Marine Parks") +
  guides(fill = guide_legend(order = 1)) +
  annotate("rect", xmin = 113, xmax = 113.7, ymin = -28.15, ymax = -27.05,
           colour = "goldenrod1", fill = "white", alpha = 0.2, size = 0.6) +
  annotate("point", y = c(-28.7761, -27.7115, -24.8838), x = c(114.6113, 114.1714, 113.6571), size = 0.75) +
  annotate("text", y = c(-28.7761, - 27.7115, -24.8838), x = c(115, 114.48, 114.075),
           label = c("Geraldton", "Kalbarri", "Carnarvon"), size = 3) +
  coord_sf(xlim = c(109.4, 115.0607), ylim = c(-29.25, -24.4)) +
  theme_minimal()
p1

# inset map
p2 <- ggplot(data = aus) +
  geom_sf(fill = "seashell1", colour = "grey90", size = 0.05, alpha = 4/5) +
  geom_sf(data = sw_mpa, alpha = 5/6, colour = "grey85", size = 0.02) +
  # geom_sf(data = ab_mpa, alpha = 4/5, colour = "grey85") +
  coord_sf(xlim = c(108, 125), ylim = c(-37, -13)) +
  annotate("rect", xmin = 108.9, xmax = 115.0607, ymin = -29.4, ymax = -24.2, 
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey70"))
p2

# plot both 
p2 + p1 + plot_layout(widths = c(0.8, 2.2))

ggsave("plots/spatial/locplot.png", dpi = 200, width = 10, height = 6)


# site zoom plots

# get sampling data
bossd <- readRDS('data/2105_abrolhos_boss.rds')
bruvd <- readRDS('data/2105_abrolhos_bruv.rds')

# reduce mpa levels for these plots
snmpa_cols <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
                                          "Multiple Use Zone" = "#b9e6fb",
                                          "Special Purpose Zone" = "#6daff4"
))

swampa_cols <- scale_fill_manual(values = c(
  "Fish Habitat Protection Area" = "#fbff85"
))

# closer plot
sitebathy <- readRDS('output/ga_bathy_fine.rds')                                # finer bathy
colnames(sitebathy)[3] <- "Depth"
# sitebathy <- sitebathy[sitebathy$Depth > -1000, ]                               # trim to reduce legend
sitebathy <- sitebathy[sitebathy$x > 112.7 & sitebathy$x < 115, ] 
sitebathy <- sitebathy[sitebathy$y > -30 & sitebathy$y < -26.6, ]

p3 <- ggplot() +
  # geom_raster(data = sitebathy, aes(x, y, fill = Depth), alpha = 4/5) +
  geom_contour_filled(data = sitebathy, aes(x = x, y = y, z = Depth,
                                             fill = after_stat(level)),
                      breaks = c(0, -30, -70, -200, -700, -2000, -4000, -10000), alpha = 4/5) +
  scale_fill_grey(start = 1, end = 0.5 , guide = "none") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +  
  geom_sf(data = terrnp, aes(fill = leg_catego), alpha = 4/5, colour = NA) +
  labs(fill = "State Managed Areas") +
  waterr_cols +
  new_scale_fill() +
  geom_sf(data = ab_nmp, aes(fill = ZoneName), alpha = 3/5, colour = NA) +
  snmpa_cols +
  labs(x = NULL, y = NULL, fill = "Australian Marine Park") +
  geom_contour(data = sitebathy, aes(x = x, y = y, z = Depth), 
               breaks = c(0, -30, -70, -200, - 700, - 9000), colour = "white", alpha = 1, size = 0.2) +
  geom_point(data = bruvd, aes(Longitude, Latitude, colour = "indianred4"), 
             alpha = 3/5, shape = 10) +
  geom_point(data = bossd, aes(Longitude, Latitude, colour = "Drop Camera"),
             alpha = 3/5, shape = 10) +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  scale_colour_manual(values = c("BRUV" = "indianred4","Drop Camera" = "seagreen4")) +
  annotate("rect", xmin = 113.02, xmax = 113.29, ymin = -27.19, ymax = -27.08,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("text", x = 113.15, y = -27.05, size = 3, 
           colour = "grey20", label = "Big Bank") +
  annotate("rect", xmin = 113.24, xmax = 113.58, ymin = -28.13, ymax = -28.02,
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  annotate("text", x = 113.42, y = -27.99, size = 3,
           colour = "grey20", label = "Shallow Bank") +
  coord_sf(xlim = c(112.8, 115), ylim = c(-28.9, -26.7)) +
  annotate("text", y = c(-27.875,-27.875,-27.875,-27.875, -26.87), 
           x = c(114.07,113.32, 113.15, 112.79, 113.167), 
           label = c("30m", "70m", "200m", "700m", "70m"), size = 2) +
  annotate("point", y = c(-27.7115, -28.7761), x = c(114.1714, 114.6113), size = 0.75) +
  annotate("text", y = c(-27.7115, -28.7761), x = c(114.275 + 0.05, 114.6113 + 0.16),
           label = c("Kalbarri", "Geraldton"), size = 3) +
  annotate(geom = "segment", linetype = "dashed", 
           x = c(113.2362, 112.8363, 113.7188), 
           xend = c(114.3337, 114.0763, 114.8787), 
           y = c(-28.07, -27.13, -28.91),                                       # Shallow Bank, Big Bank, Southern Group
           yend = c(-28.07, -27.13, -28.91), alpha = 0.5, size = 0.2) +
  annotate(geom = "text", x = c(114.3337 + 0.05, 114.0763 + 0.05, 114.8787 + 0.05), 
           y = c(-28.07, -27.13, -28.91),
           label = c("b)", "a)", "c)"), size = 2) +
  labs(colour = "Sample", x = NULL, y = NULL) +
  theme_minimal()

png(filename = "plots/spatial/siteplot.png", units = "in", res = 200, width = 8, height = 6)
p3
dev.off()

# ggsave("plots/spatial/siteplot.png", dpi = 200, width = 8, height = 6)

## single site zoom plots
snmpa_cols <- scale_colour_manual(values = c("National Park Zone" = "#7bbc63"))

nsitebathy <- sitebathy[sitebathy$Depth > -160, ]                               # trim to reduce legend
nsitebathy <- nsitebathy[nsitebathy$x > 113 & nsitebathy$x < 113.3, ]
nsitebathy <- nsitebathy[nsitebathy$y > -27.2 & nsitebathy$y < -27.05, ]


p4 <- ggplot() +
  geom_raster(data = nsitebathy, aes(x, y, fill = Depth), alpha = 4/5) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_contour(data = nsitebathy, aes(x = x, y = y, z = Depth), 
               binwidth = 10, colour = "white", alpha = 1, size = 0.1) +
  geom_text_contour(data = nsitebathy, aes(x = x, y = y, z = Depth), 
                    binwidth = 10, size = 2.5,
                    label.placer = label_placer_n(1)) +
  geom_sf(data = ab_nmp, aes(colour = ZoneName), alpha = 4/5, fill = NA, show.legend = F) +
  snmpa_cols + 
  labs(x = NULL, y = NULL, colour = NULL) +
  new_scale_colour() +
  geom_point(data = bruvd, aes(Longitude, Latitude, colour = "indianred4"), 
             alpha = 3/5, shape = 10, show.legend = F) +
  # geom_point(data = bossd, aes(Longitude, Latitude, colour = "Drop Camera"), 
  #            alpha = 3/5, shape = 10) +
  # scale_colour_manual(values = c("BRUV" = "indianred4",
  #                                "Drop Camera" = "navyblue")) +
  coord_sf(xlim = c(113.02, 113.28), ylim = c(-27.18, -27.08)) +
  labs(colour = NULL, x = NULL, y = NULL) +
  theme_minimal()
p4
ggsave("plots/nthsite.png", dpi = 200, width = 7, height = 4)


snmpa_cols <- scale_colour_manual(values = c("National Park Zone" = "#7bbc63"
                                          #  "Multiple Use Zone" = "#b9e6fb",
                                          # "Special Purpose Zone" = "#6daff4"
))
sab_nmp <- ab_nmp[ab_nmp$ZoneName == "National Park Zone", ]

ssitebathy <- sitebathy[sitebathy$Depth < -10 & sitebathy$Depth > -380, ]                               # trim to reduce legend
ssitebathy <- ssitebathy[ssitebathy$x > 113.23 & ssitebathy$x < 113.6, ]
ssitebathy <- ssitebathy[ssitebathy$y > -28.15 & ssitebathy$y < -28, ]

p5 <- ggplot() +
  geom_raster(data = ssitebathy, aes(x, y, fill = Depth), alpha = 4/5) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_contour(data = ssitebathy, aes(x = x, y = y, z = Depth), 
               binwidth = 20, colour = "white", alpha = 4/5, size = 0.1) +
  geom_text_contour(data = ssitebathy, aes(x = x, y = y, z = Depth), 
                    binwidth = 20, size = 2.5, label.placer = label_placer_n(1)) +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = sab_nmp, aes(colour = ZoneName), alpha = 1, fill = NA, show.legend = F) +
  snmpa_cols + 
  labs(x = NULL, y = NULL, colour = NULL) +
  new_scale_colour() +
  geom_point(data = bruvd, aes(Longitude, Latitude, colour = "indianred4"), 
             alpha = 3/5, shape = 10, show.legend = F) +
  # geom_point(data = bossd, aes(Longitude, Latitude, colour = "Drop Camera"), 
  #            alpha = 3/5, shape = 10) +
  # scale_colour_manual(values = c("BRUV" = "indianred4",
  #                                "Drop Camera" = "navyblue")) +
  coord_sf(xlim = c(113.24, 113.58), ylim = c(-28.125, -28.03)) +
  labs(colour = NULL, x = NULL, y = NULL) +
  theme_minimal()
p5

ggsave("plots/sthsite.png", dpi = 200, width = 7, height = 4)

formolly <- p4/p5 + plot_annotation(tag_levels = "a")
formolly

ggsave("plots/mollyforthesis.png", dpi = 200, width = 7, height = 8)

# Key Ecological Features plot
# KEF colours
kef_cols <- scale_fill_manual(values = c("Ancient coastline" = "#ff6db6",                             
                                         "Western rock lobster" = "#6db6ff",
                                         "West coast canyons" = "#dbb865", # Too similar to west coast lagoons, changed
                                         "West coast lagoons" = "#188e8e",
                                         "Abrolhos Islands" = "#2bf446",
                                         "Western demersal fish" = "#016dda",
                                         "Wallaby Saddle" = "#940000"))

nmpa_colors <- scale_color_manual(values = c("National Park Zone" = "#7bbc63",
                                          "Multiple Use Zone" = "#b9e6fb",
                                          "Special Purpose Zone" = "#6daff4",
                                          "Habitat Protection Zone" = "#fff8a3"
))

# Reorder KEFs so they plot in a sensible order
kef$NAME <- factor(kef$NAME, levels = c("Western rock lobster", "Western demersal fish", "Wallaby Saddle", 
                                        "Abrolhos Islands", "Ancient coastline", 
                                        "West coast canyons", "West coast lagoons"))

# Reorder aus national MPs so they plot in a sensible order
ab_nmp$ZoneName <- factor(ab_nmp$ZoneName, levels = c("Multiple Use Zone", "Special Purpose Zone",
                                                      "Habitat Protection Zone", "National Park Zone"))

# Merge the zones into one polygon for each park
# ab_nmpspat <- as_Spatial(ab_nmp)
# ab_nmpmerge <- unionSpatialPolygons(ab_nmpspat, ab_nmpspat$ResName)
# ab_nmpmerge <- as(ab_nmpmerge, "sf")

# build basic plot elements
p7 <- ggplot() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = dirkh, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  geom_sf(data = terrnp, aes(fill = leg_catego), alpha = 4/5, colour = NA, show.legend = F) +
  labs(fill = "Terrestrial Managed Areas") +
  waterr_cols +
  new_scale_fill() +
  geom_sf(data = kef, aes(fill = NAME), alpha = 0.7, color = NA) +
  kef_cols+
  geom_sf(data = ab_nmp, fill = NA, alpha = 1, aes(color = ZoneName), show.legend = F, size = 0.4) +
  nmpa_colors +
  # geom_sf(data = ab_nmpmerge, fill = NA, alpha = 2/5, color = "black", show.legend = F, size = 0.2) +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.2) +
  labs(x = NULL, y = NULL,  fill = "Key Ecological Features") +
  guides(fill = guide_legend(order = 1)) +
  coord_sf(xlim = c(109.4, 115.0607), ylim = c(-29.25, -24.4)) +
  annotate("point", y = c(-28.7761, -27.7115, -24.8838), x = c(114.6113, 114.1714, 113.6571), size = 0.75) +
  annotate("text", y = c(-28.7761, - 27.7115, -24.8838), x = c(115, 114.48, 114.075),
           label = c("Geraldton", "Kalbarri", "Carnarvon"), size = 2.3) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png(filename = "plots/spatial/key-ecological-features.png", height = 4, width = 6,
    units = "in", res = 200)
p7
dev.off()


# Old sea level map
depth_cols <- scale_fill_manual(values = c("#b8d9a9","#8dbc80", "#5d9d52"),
                                labels = c("9-10 Ka", "15-17 Ka", "20-30 Ka"),
                                name = "Shoreline age")

# read in and merge GA coarse bathy tiles from https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/67703
cbaths <- list.files("data/spatial/raster", "*tile", full.names = TRUE)
cbathy <- lapply(cbaths, function(x){read.table(file = x, header = TRUE, sep = ",")})
bathdf <- do.call("rbind", lapply(cbathy, as.data.frame))

# Convert back to a raster and smooth it out
bath_r <- rasterFromXYZ(bathdf)
e <- extent(111, 115.2,-29.75, -24)
bath_r <- crop(bath_r, e) # Crop to the map area
bath_r <- disaggregate(bath_r, 5, method = 'bilinear')
# Back to dataframe to plot
bathdf <- as.data.frame(bath_r, xy = T)
colnames(bathdf)[3] <- "Depth"
# Cut out shallow bits of the depth tile so it doesn't clash with the old coastlines
bathdf_f <- bathdf %>%
  dplyr::filter(Depth < -100)

# build basic plot elements
p8 <- ggplot() +
  geom_tile(data = bathdf_f, aes(x = x, y = y, fill = Depth)) +
  scale_fill_gradient2(low = "royalblue4", mid = "lightskyblue1", high = "white", name = "Depth (m)") +
  new_scale_fill() +
  geom_contour_filled(data = bathdf, aes(x = x, y = y, z = Depth,
                                         fill = after_stat(level)),
                      breaks = c(0, -40, -70, -125)) +
  depth_cols +
  new_scale_fill() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey62", size = 0.2) +
  geom_sf(data = dirkh, fill = "seashell2", colour = "grey62", size = 0.2) +
  geom_sf(data = abro, fill = "seashell2", colour = "grey62", size = 0.2) +
  new_scale_fill() +
  geom_sf(data = ab_nmp%>%dplyr::filter(!ZoneName %in% "National Park Zone"), colour = "grey61", size = 0.2, fill = NA) +
  geom_sf(data = ab_nmp%>%dplyr::filter(ZoneName %in% "National Park Zone"), colour = "#7bbc63", size = 0.55, fill = NA) +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 0.7, size = 0.3) +
  annotate("point", y = c(-28.7761, -27.7115, -24.8838), x = c(114.6113, 114.1714, 113.6571), size = 0.75) +
  annotate("text", y = c(-28.7761, - 27.7115, -24.8838), x = c(114.9, 114.38, 113.975),
           label = c("Geraldton", "Kalbarri", "Carnarvon"), size = 3) +
  coord_sf(xlim = c(111.5, 115.0607), ylim = c(-29.25, -24.4)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "#b8d9a9"))
png("plots/spatial/old-sealevel.png", res = 200, width = 7, height = 9, units = "in")
p8
dev.off()


## saving all colours for later

# 
# # assign mpa colours
# nmpa_cols <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
#                                           "Habitat Protection Zone" = "#fff8a3",# Commonwealth MPA colours
#                                           # "Habitat Protection Zone (Reefs)" = "#fbff85",
#                                           "Multiple Use Zone" = "#b9e6fb",
#                                           # "Recreational Use Zone" = "#ffb36b",
#                                           # "Sanctuary Zone" = "#f7c0d8",
#                                           # "Special Purpose Zone (Mining Exclusion)" = "#368ac1",
#                                           # "Special Purpose Zone (Trawl)" = "#3e8ec4",
#                                           "Special Purpose Zone" = "#6daff4"
# ))
# 
# wampa_cols <- scale_fill_manual(values = c("Fish Habitat Protection Area" = "#fbff85",
#                                            # "Sanctuary Zone" = "#bfd054",
#                                            # "Marine Nature Reserve" = "#bfd054",
#                                            # "Conservation Area" = "#b3a63d",
#                                            # "Habitat Protection Zone" = "#fffbcc",# State MPA colours
#                                            # "National Park Zone" = "#a4d194",
#                                            # "General Use Zone" = "#bddde1",
#                                            # "Recreation Zone" = "#f4e952",
#                                            # "Special Purpose Zone" = "#c5bcc9",
#                                            # "Special Purpose Zone\n(Shore Based Activities)" = "#ba3030"
#                                            # "Special Purpose Zone\n(Habitat Protection)" = "#f0ac41",
#                                            # "Marine Management Area" = "#b7cfe1",
#                                            "Reef Observation Area" = "#ddccff"
# ))
