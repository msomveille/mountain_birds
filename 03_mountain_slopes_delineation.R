## This script contains the code for delineating mountain slopes ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

# Load the the shapefile of mountain ranges prepared using the script 02_mountain_ranges_ebird.R
mountain_ranges <- st_read("mountain_ranges_2.shp")

# Two mountain ranges were removed from the dataset: Dinaric Alps and Siwalik Hills, because ebird data was too spread out and it wasn't possible to clearly delineate mountain slopes with sufficient ebird data
mountain_ranges <- mountain_ranges[!(mountain_ranges$NAME %in% c("Dinaric Alps", "Siwalik Hills")),]

mountain_slopes_polygons <- ebrd2_slopes <- list()

## ALPS 
## Swiss Alps
j = 1
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split the mountain range longitudinally
qq <- seq(ext(mountain_range)[1], ext(mountain_range)[2], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(qq[i], qq[i+1]), ymin(mountain_range), ymax(mountain_range)))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,2], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("latitude", "elevation")
  ss <- smooth.spline(elev$latitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(mean(c(qq[i], qq[i+1])), ss$x[which.max(ss$y)])
}
coords_highestPoint <- append(coords_highestPoint, list(c(xmax(mountain_range), coords_highestPoint[[i]][2])), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(xmin(mountain_range), coords_highestPoint[[1]][2])), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = 6.75, ymin = ymin(mountain_range), xmax = 8.75, ymax = ymax(mountain_range))
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"Y"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "north", "south"))
sf_use_s2(FALSE)
mountain_range_slope_north <- bbox_splitted %>% 
  filter(position == "north") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Swiss_Alps <- mountain_range_slope_north
# Extract ebird data in slope polygons
ebrd2_summer_sf_alps <- st_filter(ebrd2_summer_sf, mountain_range_slope_north)
ebrd2_winter_sf_alps <- st_filter(ebrd2_winter_sf, mountain_range_slope_north)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_alps) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Swiss_Alps")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_alps) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Swiss_Alps")
ebrd2_slopes$Swiss_Alps <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_alps, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_alps, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 7, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## ATLAS MOUNTAINS 
## Morocco
j = 2
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split the mountain range longitudinally
qq <- seq(ext(mountain_range)[1], ext(mountain_range)[2], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(qq[i], qq[i+1]), ymin(mountain_range), ymax(mountain_range)))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,2], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("latitude", "elevation")
  ss <- smooth.spline(elev$latitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(mean(c(qq[i], qq[i+1])), ss$x[which.max(ss$y)])
}
coords_highestPoint <- append(coords_highestPoint, list(c(xmax(mountain_range), coords_highestPoint[[i]][2])), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(xmin(mountain_range), coords_highestPoint[[1]][2])), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 30.5, xmax = xmax(mountain_range), ymax = 32.5)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"Y"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "north", "south"))
sf_use_s2(FALSE)
mountain_range_slope_north <- bbox_splitted[1,] %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Morocco <- mountain_range_slope_north
# Extract ebird data in slope polygons
ebrd2_summer_sf_atlas <- st_filter(ebrd2_summer_sf, mountain_range_slope_north)
ebrd2_winter_sf_atlas <- st_filter(ebrd2_winter_sf, mountain_range_slope_north)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_atlas) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Morocco")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_atlas) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Morocco")
ebrd2_slopes$Morocco <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() +
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_atlas, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_atlas, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 7, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## Australian Alps
j = 3
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
mountain_slopes_polygons$Australian_Alps <- mountain_ranges[j,]
# Extract ebird data in slope polygons
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Australian_Alps")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Australian_Alps")
ebrd2_slopes$Australian_Alps <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() +
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf, col="red", size=0.05) +
  geom_sf(data=mountain_range, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf, col="red", size=0.05) +
  geom_sf(data=mountain_range, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## Blue Ridge (North Carolina)
j = 4
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split the mountain range latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 35, xmax = xmax(mountain_range), ymax = 37)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub) %>% st_geometry()
mountain_slopes_polygons$Blue_Ridge <- mountain_range_slope_east
# Extract ebird data in slope polygons
ebrd2_summer_sf_blueRidge <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_blueRidge <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_blueRidge) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Blue_Ridge")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_blueRidge) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Blue_Ridge")
ebrd2_slopes$Blue_Ridge <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_blueRidge, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_blueRidge, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 9, height = 5, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## CASCADE RANGE
j = 5
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 43.8, xmax = xmax(mountain_range), ymax = 45.8)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub) %>% st_geometry()
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub) %>% st_geometry()
mountain_slopes_polygons$Cascade_East <- mountain_range_slope_east
mountain_slopes_polygons$Cascade_West <- mountain_range_slope_west
# Extract ebird data in slope polygons
# East
ebrd2_summer_sf_east <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_east <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cascade_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cascade_East")
ebrd2_slopes$Cascade_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# West
ebrd2_summer_sf_west <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_west <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cascade_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cascade_West")
ebrd2_slopes$Cascade_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 10, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## COAST RANGES
j=6
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 39.4, xmax = xmax(mountain_range), ymax = 41.4)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Coastal_California <- mountain_range_slope_west
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Coastal_California")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Coastal_California")
ebrd2_slopes$Coastal_California <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 11, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## COLUMBIA MTS. (eastern British Columbia)
j = 7
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 48.8, xmax = xmax(mountain_range), ymax = 50.8)
mountain_slopes_polygons$Eastern_BC <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Eastern_BC")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Eastern_BC")
ebrd2_slopes$Eastern_BC <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## Cord. Cantábrica
j = 8
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges longitudinally
qq <- seq(ext(mountain_range)[1], ext(mountain_range)[2], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(qq[i], qq[i+1]), ymin(mountain_range), ymax(mountain_range)))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,2], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("latitude", "elevation")
  ss <- smooth.spline(elev$latitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(mean(c(qq[i], qq[i+1])), ss$x[which.max(ss$y)])
}
coords_highestPoint <- append(coords_highestPoint, list(c(xmax(mountain_range), coords_highestPoint[[i]][2])), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(xmin(mountain_range), coords_highestPoint[[1]][2])), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = -6, ymin = ymin(mountain_range), xmax = -4, ymax = ymax(mountain_range))
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"Y"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "north", "south"))
sf_use_s2(FALSE)
mountain_range_slope_north <- bbox_splitted %>% 
  filter(position == "north") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Cord_Cantabrica <- mountain_range_slope_north
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_north)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_north)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cord_Cantabrica")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cord_Cantabrica")
ebrd2_slopes$Cord_Cantabrica <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 8, height = 5, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## CORD. CENTRAL (Central Peru)
j = 9
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -11.4, xmax = xmax(mountain_range), ymax = -9.4)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Peru_Central <- mountain_range_slope_east
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Peru_Central")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Peru_Central")
ebrd2_slopes$Peru_Central <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## Cord. Merida
j = 10
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 8, xmax = xmax(mountain_range), ymax = 10)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Cord_Merida <- mountain_range_slope_west
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cord_Merida")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Cord_Merida")
ebrd2_slopes$Cord_Merida <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## CORDILLERA OCCIDENTAL
j = 11
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -17.4, xmax = xmax(mountain_range), ymax = -15.4)
mountain_slopes_polygons$Peru_SouthWest <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Peru_SouthWest")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Peru_SouthWest")
ebrd2_slopes$Peru_SouthWest <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 10, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()



## CORDILLERA ORIENTAL
j = 12
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree bands
mountain_range_sub_1 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -12.5, xmax = xmax(mountain_range), ymax = -14.5) # Cuzco 
mountain_range_sub_2 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -22.5, xmax = xmax(mountain_range), ymax = -24.5) # Northern Argentina
mountain_slopes_polygons$Peru_East <- mountain_range_sub_1
mountain_slopes_polygons$Argentina_North <- mountain_range_sub_2
# Extract ebird data in slope polygons
# Peru East
ebrd2_summer_sf_peru <- st_filter(ebrd2_summer_sf, mountain_range_sub_1)
ebrd2_winter_sf_peru <- st_filter(ebrd2_winter_sf, mountain_range_sub_1)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_peru) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Peru_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_peru) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Peru_East")
ebrd2_slopes$Peru_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# Argentina North
ebrd2_summer_sf_argen <- st_filter(ebrd2_summer_sf, mountain_range_sub_2)
ebrd2_winter_sf_argen <- st_filter(ebrd2_winter_sf, mountain_range_sub_2)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_argen) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Argentina_North")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_argen) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Argentina_North")
ebrd2_slopes$Argentina_North <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_peru, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_argen, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub_1, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_sub_2, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_peru, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_argen, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub_1, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_sub_2, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7.5, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()



## CORD. OCCIDENTAL
j = 13
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree bands
mountain_range_sub_1 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -1, xmax = xmax(mountain_range), ymax = 1) 
mountain_range_sub_2 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 4, xmax = xmax(mountain_range), ymax = 6) 
mountain_slopes_polygons$Ecuador_West <- mountain_range_sub_1
mountain_slopes_polygons$Colombia_West <- mountain_range_sub_2
# Extract ebird data in slope polygons
# Ecuador_West
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub_1)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub_1)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Ecuador_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Ecuador_West")
ebrd2_slopes$Ecuador_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# Colombia_West
ebrd2_summer_sf_3 <- st_filter(ebrd2_summer_sf, mountain_range_sub_2)
ebrd2_winter_sf_3 <- st_filter(ebrd2_winter_sf, mountain_range_sub_2)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_3) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Colombia_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_3) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Colombia_West")
ebrd2_slopes$Colombia_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_3, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub_1, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_sub_2, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_3, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub_1, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_sub_2, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 7, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()



## DRAKENSBERG
j = 14
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -30.2, xmax = xmax(mountain_range), ymax = -28.2) # Cuzco 
mountain_slopes_polygons$Drakensberg <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Drakensberg")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Drakensberg")
ebrd2_slopes$Drakensberg <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7.5, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## Khasi Hills
j = 15
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = 90.8, ymin = ymin(mountain_range), xmax = 92.8, ymax = ymax(mountain_range)) 
mountain_slopes_polygons$Khasi_Hills <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Khasi_Hills")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Khasi_Hills")
ebrd2_slopes$Khasi_Hills <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 8, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()



## PYRENEES
j = 16
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges longitudinally
qq <- seq(ext(mountain_range)[1], ext(mountain_range)[2], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(qq[i], qq[i+1]), ymin(mountain_range), ymax(mountain_range)))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,2], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("latitude", "elevation")
  ss <- smooth.spline(elev$latitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(mean(c(qq[i], qq[i+1])), ss$x[which.max(ss$y)])
}
coords_highestPoint <- append(coords_highestPoint, list(c(xmax(mountain_range), coords_highestPoint[[i]][2])), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(xmin(mountain_range), coords_highestPoint[[1]][2])), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree bands
mountain_range_sub_1 <- st_crop(mountain_ranges[j,], xmin = -1.9, ymin = ymin(mountain_range), xmax = 0.1, ymax = ymax(mountain_range))
mountain_range_sub_2 <- st_crop(mountain_ranges[j,], xmin = 1.2, ymin = ymin(mountain_range), xmax = 3.2, ymax = ymax(mountain_range))
# Create polygons for mountain slopes 
bbox_splitted_1 <- mountain_range_sub_1 %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"Y"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "north", "south"))
bbox_splitted_2 <- mountain_range_sub_2 %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"Y"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "north", "south"))
sf_use_s2(FALSE)
mountain_range_slope_north <- bbox_splitted_1 %>% 
  filter(position == "north") %>% st_intersection(mountain_range_sub_1) %>% st_geometry()
mountain_range_slope_south <- bbox_splitted_2 %>% 
  filter(position == "south") %>% st_intersection(mountain_range_sub_2) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Pyrenees_Atlantic <- mountain_range_slope_north
mountain_slopes_polygons$Pyrenees_Catalonia <- mountain_range_slope_south
# Extract ebird data in slope polygons
# Pyrenees Atlantic
ebrd2_summer_sf_Atlantic <- st_filter(ebrd2_summer_sf, mountain_range_slope_north)
ebrd2_winter_sf_Atlantic <- st_filter(ebrd2_winter_sf, mountain_range_slope_north)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_Atlantic) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Pyrenees_Atlantic")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_Atlantic) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Pyrenees_Atlantic")
ebrd2_slopes$Pyrenees_Atlantic <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# Pyrenees Catalonia
ebrd2_summer_sf_Catalonia <- st_filter(ebrd2_summer_sf, mountain_range_slope_south)
ebrd2_winter_sf_Catalonia <- st_filter(ebrd2_winter_sf, mountain_range_slope_south)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_Catalonia) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Pyrenees_Catalonia")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_Catalonia) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Pyrenees_Catalonia")
ebrd2_slopes$Pyrenees_Catalonia <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_Atlantic, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_Catalonia, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_south, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_Atlantic, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_Catalonia, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_south, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## Serra da Mantiqueira
j = 17
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -24, xmax = xmax(mountain_range), ymax = -22) 
mountain_slopes_polygons$Serra_Mantiqueira <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Serra_Mantiqueira")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Serra_Mantiqueira")
ebrd2_slopes$Serra_Mantiqueira <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 7, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## Sierra Chiapas
j = 18
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Custom split
coords_highestPoint <- data.frame(longitude = c(-91.6, -91.6, -94.2),
                                  latitude = c(ymin(mountain_range), 16, ymax(mountain_range)))
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = -93, ymin = ymin(mountain_range), xmax = -91, ymax = ymax(mountain_range))
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Chiapas_East <- mountain_range_slope_east
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Chiapas_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Chiapas_East")
ebrd2_slopes$Chiapas_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 10, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## Sierra Madre del Sur
j = 19
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Custom split
coords_highestPoint <- data.frame(longitude = c(xmax(mountain_range), -96, -98.5),
                                  latitude = c(16.9, 16.9, ymax(mountain_range)))
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree bands
mountain_range_sub_1 <- st_crop(mountain_ranges[j,], xmin = -97.5, ymin = ymin(mountain_range), xmax = -95.5, ymax = ymax(mountain_range)) # Veracruz 
mountain_range_sub_2 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = ymin(mountain_range), xmax = xmin(mountain_range)+2, ymax = ymax(mountain_range)) # Jalisco 
# Create polygons for mountain slopes 
bbox_splitted_1 <- mountain_range_sub_1 %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_east <- bbox_splitted_1 %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub_1) %>% st_geometry()
mountain_range_slope_west <- mountain_range_sub_2
sf_use_s2(TRUE)
mountain_slopes_polygons$Veracruz <- mountain_range_slope_east
mountain_slopes_polygons$Jalisco <- mountain_range_slope_west
# Extract ebird data in slope polygons
# Veracruz
ebrd2_summer_sf_east <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_east <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Veracruz")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Veracruz")
ebrd2_slopes$Veracruz <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# Jalisco
ebrd2_summer_sf_west <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_west <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Jalisco")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Jalisco")
ebrd2_slopes$Jalisco <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 7, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## SIERRA MADRE OCCIDENTAL
j = 20
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 23, xmax = xmax(mountain_range), ymax = 25)
mountain_slopes_polygons$Sierra_Occidental <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Sierra_Occidental")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Sierra_Occidental")
ebrd2_slopes$Sierra_Occidental <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## SIERRA MADRE ORIENTAL
j = 21
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Custom split
coords_highestPoint <- data.frame(longitude = c(-101.2, -101.2),
                                  latitude = c(ymin(mountain_range), ymax(mountain_range)))
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 24, xmax = xmax(mountain_range), ymax = 26) 
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub) %>% st_geometry()
mountain_slopes_polygons$Nuevo_Leon <- mountain_range_slope_east
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Nuevo_Leon")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Nuevo_Leon")
ebrd2_slopes$Nuevo_Leon <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## SIERRA NEVADA
j = 22
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 36.2, xmax = xmax(mountain_range), ymax = 38.2)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$California_Sierra_Nevada <- mountain_range_slope_west
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="California_Sierra_Nevada")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="California_Sierra_Nevada")
ebrd2_slopes$California_Sierra_Nevada <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## S. Nevada
j = 23
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
mountain_slopes_polygons$S_Nevada <- mountain_ranges[j,]
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="S_Nevada")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="S_Nevada")
ebrd2_slopes$S_Nevada <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 8, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## SOUTHERN ALPS
j = 24
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges longitudinally
qq <- seq(ext(mountain_range)[1], ext(mountain_range)[2], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(qq[i], qq[i+1]), ymin(mountain_range), ymax(mountain_range)))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,2], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("latitude", "elevation")
  ss <- smooth.spline(elev$latitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(mean(c(qq[i], qq[i+1])), ss$x[which.max(ss$y)])
}
coords_highestPoint <- append(coords_highestPoint, list(c(xmax(mountain_range), coords_highestPoint[[i]][2])), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(xmin(mountain_range), coords_highestPoint[[1]][2])), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = 169.25, ymin = ymin(mountain_range), xmax = 171.25, ymax = ymax(mountain_range))
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"Y"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "north", "south"))
sf_use_s2(FALSE)
mountain_range_slope_north <- bbox_splitted %>% 
  filter(position == "north") %>% st_intersection(mountain_range_sub) %>% st_geometry()
mountain_range_slope_south <- bbox_splitted %>% 
  filter(position == "south") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Southern_Alps_West <- mountain_range_slope_north
mountain_slopes_polygons$Southern_Alps_East <- mountain_range_slope_south
# Extract ebird data in slope polygons
# East
ebrd2_summer_sf_east <- st_filter(ebrd2_summer_sf, mountain_range_slope_south)
ebrd2_winter_sf_east <- st_filter(ebrd2_winter_sf, mountain_range_slope_south)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Alps_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Alps_East")
ebrd2_slopes$Southern_Alps_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# West
ebrd2_summer_sf_west <- st_filter(ebrd2_summer_sf, mountain_range_slope_north)
ebrd2_winter_sf_west <- st_filter(ebrd2_winter_sf, mountain_range_slope_north)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Alps_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Alps_West")
ebrd2_slopes$Southern_Alps_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_south, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_north, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_south, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 8, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## Southern Ghats
j = 25
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_ranges[j,]) %>% st_geometry()
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_ranges[j,]) %>% st_geometry()
mountain_slopes_polygons$Southern_Ghats_East <- mountain_range_slope_east
mountain_slopes_polygons$Southern_Ghats_West <- mountain_range_slope_west
# Extract ebird data in slope polygons
# East
ebrd2_summer_sf_east <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_east <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Ghats_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Ghats_East")
ebrd2_slopes$Southern_Ghats_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# West
ebrd2_summer_sf_west <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_west <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Ghats_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Southern_Ghats_West")
ebrd2_slopes$Southern_Ghats_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 7, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## SOUTHERN ANDES
j = 26
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = -34.5, xmax = xmax(mountain_range), ymax = -32.5)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Chile_Central <- mountain_range_slope_west
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Chile_Central")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Chile_Central")
ebrd2_slopes$Chile_Central <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## EASTERN HIMALAYAS (Bhutan)
j = 27
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = 88.2, ymin = ymin(mountain_range), xmax = 90.2, ymax = ymax(mountain_range))
mountain_slopes_polygons$Eastern_Himalayas <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Eastern_Himalayas")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Eastern_Himalayas")
ebrd2_slopes$Eastern_Himalayas <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## WESTERN HIMALAYAS
j = 28
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 31.5, xmax = xmax(mountain_range), ymax = 33.5)
mountain_slopes_polygons$Western_Himalayas <- mountain_range_sub
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_range_sub)
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_range_sub)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Western_Himalayas")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Western_Himalayas")
ebrd2_slopes$Western_Himalayas <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range_sub, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 8, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


## Hawaii
j = 29
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
mountain_slopes_polygons$Hawaii <- mountain_ranges[j,]
# Extract ebird data in slope polygons
ebrd2_summer_sf_2 <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf_2 <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Hawaii")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_2) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Hawaii")
ebrd2_slopes$Hawaii <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_2, col="red", size=0.05) +
  geom_sf(data=mountain_range, fill=NA, linewidth=2, col="orange2") +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 9, height = 5, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## Taiwan
j = 30
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree bands
mountain_range_sub_1 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 22.5, xmax = xmax(mountain_range), ymax = 24.5)
mountain_range_sub_2 <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 23, xmax = xmax(mountain_range), ymax = 25)
# Create polygons for mountain slopes 
bbox_splitted_1 <- mountain_range_sub_1 %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
bbox_splitted_2 <- mountain_range_sub_2 %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_west <- bbox_splitted_1 %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub_1) %>% st_geometry()
mountain_range_slope_east <- bbox_splitted_2 %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub_2) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Taiwan_East <- mountain_range_slope_east
mountain_slopes_polygons$Taiwan_West <- mountain_range_slope_west
# Extract ebird data in slope polygons
# East 
ebrd2_summer_sf_east <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_east <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Taiwan_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Taiwan_East")
ebrd2_slopes$Taiwan_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# West
ebrd2_summer_sf_west <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_west <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Taiwan_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Taiwan_West")
ebrd2_slopes$Taiwan_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 6, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=1, ncol=2, legend="bottom", common.legend = T)
dev.off()


## Cord. Centroamericana Sur
j = 31
mountain_range <- terra::vect(mountain_ranges[j,])
# Get ebird data
ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[j],".csv")) %>%
  filter(duplicated(sampling_event_identifier) == F)
ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[j,])
ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[j,])
elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[j,])
elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[j,])
# Split mountain ranges latitudinally
qq <- seq(ext(mountain_range)[3], ext(mountain_range)[4], 0.25)
coords_highestPoint <- list()
for(i in 1:(length(qq)-1)){
  # crop elevation raster to the elevational bin
  elevation_raster_2 <- terra::crop(elevation_raster_1, terra::ext(c(xmin(mountain_range), xmax(mountain_range), qq[i], qq[i+1])))
  # smooth spline of elevation against longitude within the elevational bin
  elev <- cbind(crds(terra::as.points(elevation_raster_2))[,1], values(terra::as.points(elevation_raster_2)))
  colnames(elev) <- c("longitude", "elevation")
  ss <- smooth.spline(elev$longitude, elev$elevation, spar=0.5)
  # Get coordinates of the mode of the smooth spline
  coords_highestPoint[[i]] <- c(ss$x[which.max(ss$y)], mean(c(qq[i], qq[i+1])))
}
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[i]][1], ymax(mountain_range))), after=i)
coords_highestPoint <- append(coords_highestPoint, list(c(coords_highestPoint[[1]][1], ymin(mountain_range))), after=0)
coords_highestPoint <- as.data.frame(do.call(rbind, coords_highestPoint))
colnames(coords_highestPoint) <- c("longitude", "latitude")
coords_highestPoint_sf <- coords_highestPoint %>% sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(crs(mountain_range))
coords_highestPoint_sf_line <- coords_highestPoint_sf %>% st_coordinates() %>% st_linestring() %>% st_sfc() %>% sf::st_set_crs(crs(mountain_range))
# Crop mountain range to 2-degree band
mountain_range_sub <- st_crop(mountain_ranges[j,], xmin = xmin(mountain_range), ymin = 8.7, xmax = -81.5, ymax = 10.7)
# Create polygons for mountain slopes 
bbox_splitted <- mountain_range_sub %>% st_bbox() %>% st_as_sfc() %>% 
  lwgeom::st_split(coords_highestPoint_sf_line) %>% 
  st_collection_extract("POLYGON") %>% st_as_sf() %>%
  mutate(xpos = st_coordinates(st_centroid(.))[,"X"]) %>% 
  mutate(position = ifelse(xpos == max(xpos), "east", "west"))
sf_use_s2(FALSE)
mountain_range_slope_west <- bbox_splitted %>% 
  filter(position == "west") %>% st_intersection(mountain_range_sub) %>% st_geometry()
mountain_range_slope_east <- bbox_splitted %>% 
  filter(position == "east") %>% st_intersection(mountain_range_sub) %>% st_geometry()
sf_use_s2(TRUE)
mountain_slopes_polygons$Costa_Rica_East <- mountain_range_slope_east
mountain_slopes_polygons$Costa_Rica_West <- mountain_range_slope_west
# Extract ebird data in slope polygons
# East 
ebrd2_summer_sf_east <- st_filter(ebrd2_summer_sf, mountain_range_slope_east)
ebrd2_winter_sf_east <- st_filter(ebrd2_winter_sf, mountain_range_slope_east)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Costa_Rica_East")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_east) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Costa_Rica_East")
ebrd2_slopes$Costa_Rica_East <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
# West
ebrd2_summer_sf_west <- st_filter(ebrd2_summer_sf, mountain_range_slope_west)
ebrd2_winter_sf_west <- st_filter(ebrd2_winter_sf, mountain_range_slope_west)
ebrd2_slope_summer <- as.data.frame(ebrd2_summer_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Costa_Rica_West")
ebrd2_slope_winter <- as.data.frame(ebrd2_winter_sf_west) %>% dplyr::select(scientificname, latitude, longitude, sampling_event_identifier, year, day, mnt_id, season) %>% mutate(slope="Costa_Rica_West")
ebrd2_slopes$Costa_Rica_West <- rbind(ebrd2_slope_summer, ebrd2_slope_winter)
##  Plot
g_summer <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_summer_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_summer_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(A) Northern summer") + theme_void() + theme(plot.title = element_text(hjust=0.5))
g_winter <- ggplot() + 
  geom_spatraster(data = elevation_raster_1, aes(fill = SRTM_global)) + 
  scale_fill_viridis_c(na.value = "transparent", name = "Elevation") +
  geom_sf(data=ebrd2_winter_sf_east, col="red", size=0.05) +
  geom_sf(data=ebrd2_winter_sf_west, col="red", size=0.05) +
  geom_sf(data=mountain_range_slope_west, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=mountain_range_slope_east, fill=NA, linewidth=2, col="orange2") +
  geom_sf(data=coords_highestPoint_sf_line, linewidth=1.1) +
  ggtitle("(B) Northern winter") + theme_void() + theme(plot.title = element_text(hjust=0.5))
png(filename = paste0("results/figures/mountain_splits/", mountain_range$NAME, ".png"), width = 6, height = 8, units = "in", res=300, bg = "white")
ggarrange(g_summer, g_winter, nrow=2, ncol=1, legend="bottom", common.legend = T)
dev.off()


save(mountain_slopes_polygons, ebrd2_slopes, file="mountain_slopes.RData")

