## This script contains the code for initial processing of the mountain range data  ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)

# Upload mountain range polygons
mountain.key <- read.csv("resources/mountain-key.csv")
mountain_ranges <- st_read("resources/ne_10m_geography_regions_polys/ne_10m_geography_regions_polys.shp") %>% dplyr::select(NAME)
mountain_ranges <- mountain_ranges[match(mountain.key$NAME, mountain_ranges$NAME),] %>% left_join(mountain.key)

# Crop the Andes to only keep the southern part
southern_Andes <- st_crop(mountain_ranges[which(mountain_ranges$NAME == "ANDES"),], c(xmin = -80, ymin = -51.9, xmax = -60, ymax = -32))
southern_Andes$NAME <- "SOUTHERN ANDES"
mountain_ranges <- rbind(mountain_ranges, southern_Andes)
mountain_ranges <- mountain_ranges[-which(mountain_ranges$NAME == "ANDES"),]

# Split Himalayas into west and east
eastern_himalayas <- st_crop(mountain_ranges[which(mountain_ranges$NAME == "HIMALAYAS"),], c(xmin = 85, ymin = 20, xmax = 100, ymax = 40))
eastern_himalayas$NAME <- "EASTERN HIMALAYAS"
western_himalayas <- st_crop(mountain_ranges[which(mountain_ranges$NAME == "HIMALAYAS"),], c(xmin = 70, ymin = 20, xmax = 85, ymax = 40)) 
western_himalayas$NAME <- "WESTERN HIMALAYAS"
mountain_ranges <- rbind(mountain_ranges, eastern_himalayas, western_himalayas)
mountain_ranges <- mountain_ranges[-which(mountain_ranges$NAME == "HIMALAYAS"),]

# Projection
st_crs(mountain_ranges) <- "+proj=lonlat +datum=WGS84 +no_defs"
mountain_ranges <- sf::st_transform(mountain_ranges, "+proj=eqearth +datum=WGS84 +no_defs")

# Add Hawaii
mountain_ranges_hawaii <- st_read("resources/ne_10m_geography_regions_polys/OBJECTID-55.shp") %>% rename(NAME = Name) %>% dplyr::select(NAME)
mountain_ranges_hawaii <- st_transform(mountain_ranges_hawaii, st_crs(mountain_ranges)) 
mountain_ranges_hawaii <- mountain_ranges_hawaii %>% add_column(mnt.id = "NA", .after=1)
mountain_ranges <- rbind(mountain_ranges, mountain_ranges_hawaii)

# Add Taiwan
mountain_ranges_taiwan <- st_read("resources/ne_10m_geography_regions_polys/OBJECTID-398.shp") %>% rename(NAME = Name) %>% dplyr::select(NAME)
mountain_ranges_taiwan <- st_transform(mountain_ranges_taiwan, st_crs(mountain_ranges))
mountain_ranges_taiwan <- mountain_ranges_taiwan %>% add_column(mnt.id = "NA", .after=1)
mountain_ranges <- rbind(mountain_ranges, mountain_ranges_taiwan)

# Add Cordillera Centroamericana Sur
mountain_ranges_cordCentroamericana <- st_read("resources/ne_10m_geography_regions_polys/Cordillera_Centroamericana.shp") %>% dplyr::rename(NAME = MapName) %>% dplyr::select(NAME)
mountain_ranges_cordCentroamericana <- st_transform(mountain_ranges_cordCentroamericana, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
mountain_ranges_cordCentroamericana_sur <- st_crop(mountain_ranges_cordCentroamericana, c(xmin = -96, ymin = 0, xmax = -80, ymax = 11))
mountain_ranges_cordCentroamericana_sur <- st_transform(mountain_ranges_cordCentroamericana_sur, st_crs(mountain_ranges))
mountain_ranges_cordCentroamericana_sur$NAME <- "Cord. Centroamericana Sur"
mountain_ranges_cordCentroamericana_sur <- mountain_ranges_cordCentroamericana_sur %>% add_column(mnt.id = "NA", .after=1)
mountain_ranges <- rbind(mountain_ranges, mountain_ranges_cordCentroamericana_sur)

# Load elevation data (downloaded from NASA Shuttle Radar Topography Mission V003)
elevation_raster <- terra::rast("resources/SRTM_2km.tif")
elevation_raster <- terra::project(elevation_raster, terra::vect(mountain_ranges[1,]))

# Elevation bins
bins <- seq(200, 8000, 200)
mountain_ranges_elevation_bins <- list()
for(k in 1:nrow(mountain_ranges)){
  mountain_range <- terra::vect(mountain_ranges[k,])
  if(terra::ext(mountain_range)[4] < terra::ext(elevation_raster)[4] & terra::ext(mountain_range)[3] > terra::ext(elevation_raster)[3]){
    elevation_raster_1 <- terra::crop(elevation_raster, mountain_range)
    elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_range)
    elevation_bins_land <- vector()
    for(i in 1:length(bins)){
      elevation_bins_land[i] <- length(which(terra::values(elevation_raster_1) >= bins[i]-200 & terra::values(elevation_raster_1) < bins[i]))
    }
    mountain_ranges_elevation_bins[[k]] <- elevation_bins_land
  }
}
toremove <- which(lapply(mountain_ranges_elevation_bins, length) == 0)
mountain_ranges <- mountain_ranges[-toremove,]
mountain_ranges_elevation_bins <- do.call(cbind, mountain_ranges_elevation_bins)
colnames(mountain_ranges_elevation_bins) <- mountain_ranges$NAME
rownames(mountain_ranges_elevation_bins) <- bins

# Remove mountain ranges with no elevation data below 600m and/or no elevation data above 2000m
toremove <- which(apply(mountain_ranges_elevation_bins, 2, function(x) ifelse(sum(x[1:3]) == 0 | sum(x[9:40]) == 0, 1, 0)) == 1)
mountain_ranges_elevation_bins <- mountain_ranges_elevation_bins[,-toremove]

# Remove mountain ranges with less than 10% of area above 1000m
mountain_ranges_elevation_bins_2 <- apply(mountain_ranges_elevation_bins, 2, function(x) x / sum(x))
toremove <- which(apply(mountain_ranges_elevation_bins_2, 2, function(x) ifelse(sum(x[6:40]) < 0.1, 1, 0)) == 1)
mountain_ranges_elevation_bins <- mountain_ranges_elevation_bins[,-toremove]
mountain_ranges <- mountain_ranges[match(colnames(mountain_ranges_elevation_bins), mountain_ranges$NAME),]

# Redundancy among mountain ranges 
#mountain_ranges_intersects <- st_intersects(mountain_ranges)
#k=0
#k=k+1
#mountain_ranges$NAME[k]
#mountain_ranges$NAME[mountain_ranges_intersects[[k]]]
#ggplot(mountain_ranges[mountain_ranges_intersects[[k]],]) + geom_sf()

# Remove mountain ranges that overlap with others
mountains_toremove <- c("Anti Atlas", "ATLAS SAHARIEN", "Atlas Tellien", "Er Rif", "HAUT ATLAS", "Moyen Atlas",
                        "ROCKY MOUNTAINS", "Transylvanian Alps", "Central Highlands", "Klamath Mts.",
                        "Maoke Mts.", "Owen Stanley Ra.", "Ruwenzori Range")

mountain_ranges_elevation_bins <- mountain_ranges_elevation_bins[,-which(colnames(mountain_ranges_elevation_bins) %in% mountains_toremove)]
mountain_ranges <- mountain_ranges[match(colnames(mountain_ranges_elevation_bins), mountain_ranges$NAME),]

# Save shapefile containing the mountain range polygons
st_write(mountain_ranges, "resources/mountain_ranges.shp")


