## This script contains the code for processing eBird data within each mountain range and filtering out mountain ranges with insufficient ebird data ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

# Load the the shapefile of mountain ranges prepared using the script 01_mountain_ranges.R
mountain_ranges <- st_read("resources/mountain_ranges.shp")
mountain_ranges$mnt_id[85:87] <- c("220", "221", "222")

##  First set of filters to remove mountains with insufficient ebird data ##

# List of terrestrial bird species
# Remove the following habitat categories associated with water: Coastal, Marine, Riverine, Wetland (+ Human Modified)
# Remove the following primary lifestyle category: Aquatic
# Remove the following trophic niches: Aquatic predator, Herbivore aquatic
avonet_data <- read_csv("resources/avonet_data.csv") %>%
  rename(SCI_NAME = Species2) %>%
  left_join(read_csv("resources/eBird_Taxonomy_v2019.csv")) %>%
  filter(Habitat %in% c("Desert", "Grassland", "Rock", "Shrubland", "Forest", "Woodland") & Primary.Lifestyle %in% c("Aerial", "Generalist", "Insessorial", "Terrestrial") & Trophic.Niche %in% c("Invertivore", "Omnivore", "Frugivore", "Nectarivore", "Granivore", "Herbivore terrestrial", "Scavenger", "Vertivore"))


# Load eBird data for each mountain range
load("resources/ebird-mountain-ranges.RData")

# Filter eBird data and split by mountain range
for(i in 1:nrow(mountain_ranges)){
  ebrd2 <- dat %>% filter(mnt_id == mountain_ranges$mnt_id[i])
  ebrd2 <- ebrd2 %>% 
    filter(scientificname %in% avonet_data$SCI_NAME) %>%
    filter(samplingprotocol %in% c("Traveling", "Stationary", "Area")) %>%
    filter(all_species_reported == "1") %>%
    filter(duration_minutes < 300) %>%
    filter(number_observers <= 10) %>%
    filter(is.na(effort_distance_km) == T | effort_distance_km <= 2) %>%
    filter(is.na(effort_area_ha) == T | effort_area_ha <= 100) 
  # Seasonal ebird data
  ebrd2_summer <- ebrd2 %>% subset(day > 151 & day < 227) %>% mutate(season = "summer") # 1 June — 15 August
  ebrd2_winter <- ebrd2 %>% subset(day < 45 | day > 334) %>% mutate(season = "winter") # 1 Dec — 15 Feb
  ebrd2 <- rbind(ebrd2_summer, ebrd2_winter)
  write_csv(ebrd2, paste0("resources/ebird_data/", mountain_ranges$NAME[i],".csv"))
  rm(ebrd2_summer, ebrd2_winter, ebrd2)
  print(i)
}

# Calculate the number of checklists available per mountain range per season
ebrd2_summer_checklists <- ebrd2_winter_checklists <- mnt_names <- vector()
ebrd2_summer_checklists_elevation <- list()
for(i in 1:nrow(mountain_ranges)){
  ebrd2 <- read_csv(paste0("resources/ebird_data/", mountain_ranges$NAME[i],".csv")) %>%
    filter(duplicated(sampling_event_identifier) == F)
  mnt_names[i] <- mountain_ranges$NAME[i]
  ebrd2_summer_checklists[i] <- length(which(ebrd2$season == "summer"))
  ebrd2_winter_checklists[i] <- length(which(ebrd2$season == "winter"))
  rm(ebrd2)
}
mnt_checklists <- data.frame(mountain = mnt_names, 
                             checklists_summer = ebrd2_summer_checklists, 
                             checklists_winter = ebrd2_winter_checklists)

# Remove mountain ranges with less than 500 checklists at both seasons
mountain_ranges <- mountain_ranges[which(mnt_checklists$checklists_winter >= 1000 & mnt_checklists$checklists_summer >= 1000),]


##  Extract elevation for each checklist ##

# Load elevation data (downloaded from NASA Shuttle Radar Topography Mission V003)
elevation_raster <- terra::rast("resources/SRTM_2km.tif")
elevation_raster <- terra::project(elevation_raster, terra::vect(mountain_ranges[1,]))

ebrd2_summer_checklists_elevation <- ebrd2_winter_checklists_elevation <- list()
for(i in 1:nrow(mountain_ranges)){
  ebrd2 <- read_csv(paste0("resources/ebird_data/", mountain_ranges$NAME[i],".csv")) %>%
    filter(duplicated(sampling_event_identifier) == F)
  ebrd2_summer_sf <- st_as_sf(x = ebrd2 %>% filter(season == "summer"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
  ebrd2_winter_sf <- st_as_sf(x = ebrd2 %>% filter(season == "winter"), coords = c("longitude", "latitude"), crs = st_crs(mountain_ranges), remove=F)
  ebrd2_summer_sf <- st_filter(ebrd2_summer_sf, mountain_ranges[i,])
  ebrd2_winter_sf <- st_filter(ebrd2_winter_sf, mountain_ranges[i,])
  elevation_raster_1 <- terra::crop(elevation_raster, mountain_ranges[i,])
  elevation_raster_1 <- terra::mask(elevation_raster_1, mountain_ranges[i,])
  ebrd2_summer_checklists_elevation[[i]] <- ebrd2_summer_sf %>% mutate(elevation = terra::extract(elevation_raster_1, terra::vect(ebrd2_summer_sf))$SRTM_global) %>% select(latitude, longitude, elevation)
  ebrd2_winter_checklists_elevation[[i]] <- ebrd2_winter_sf %>% mutate(elevation = terra::extract(elevation_raster_1, terra::vect(ebrd2_winter_sf))$SRTM_global) %>% select(latitude, longitude, elevation)
  rm(ebrd2, elevation_raster_1, ebrd2_summer_sf, ebrd2_winter_sf)
}

# Remove mountain ranges with elevational range of checklists <1500m for at least one season
tokeep <- which(unlist(lapply(ebrd2_summer_checklists_elevation, function(x) diff(range(x$elevation, na.rm=T)))) > 1500 & unlist(lapply(ebrd2_winter_checklists_elevation, function(x) diff(range(x$elevation, na.rm=T)))) > 1500)
mountain_ranges <- mountain_ranges[tokeep,]

# Save mountain range polygons as shapefile
st_write(mountain_ranges, "resources/mountain_ranges_2.shp", append=F)

