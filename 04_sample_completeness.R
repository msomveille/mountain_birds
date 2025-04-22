## This script contains the code for filtrering mountain ramges based on sampling completeness ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(iNEXT)

# Load the shapefile of mountain ranges prepared using the script 02_mountain_ranges_ebird.R
mountain_ranges <- st_read("resources/mountain_ranges_2.shp")

# Load the mountain slopes polygons and ebird data prepared using the script 03_mountain_slopes_delineation.R
load("results/mountain_slopes_4.RData")

# Load elevation data (downloaded from NASA Shuttle Radar Topography Mission V003)
elevation_raster <- terra::rast("resources/SRTM_2km.tif")
elevation_raster <- terra::project(elevation_raster, terra::vect(mountain_ranges[1,]))

## Extract elevation for all eBird observation
ebrd2_slopes <- do.call(rbind, ebrd2_slopes)
ebrd2_slopes_sf <- st_as_sf(x = ebrd2_slopes, coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_slopes_elevation <- terra::extract(elevation_raster, ebrd2_slopes_sf %>% as.data.frame %>% dplyr::select(longitude, latitude))
ebrd2_slopes$elevation <- ebrd2_slopes_elevation$SRTM_global

# Which mountain range does each slope belong to?
slopes_mountain_name <- c(1,2,3,4,5,5,6,7,8,9,10,11,12,12,13,13,15,16,17,17,18,19,20,20,21,22,23,25,26,26,27,27,28,29,30,31,32,32,33,33)


## Sample completeness

elevation_bins <- seq(0, 7000, 200)
slopes <- unique(ebrd2_slopes$slope)
completeness_summer <- completeness_winter <- list()
for(k in 1:length(slopes_mountain_name)){
  
  # Get polygon for the mountain range
  mountain_range <- terra::vect(mountain_ranges[slopes_mountain_name[k],])
  
  # Load ebird data for the mountain range
  ebrd2 <- read_csv(paste0("resources/ebird_data/", mountain_ranges$NAME[slopes_mountain_name[k]],".csv")) 
  
  # filter ebird data for checklist that are within the 2 degrees gradient
  ebrd2_slopes_sub <- ebrd2_slopes %>% filter(slope == slopes[k])
  ebrd2 <- ebrd2 %>% filter(sampling_event_identifier %in% ebrd2_slopes_sub$sampling_event_identifier)
  
  ebrd2_sf <- st_as_sf(x = ebrd2, coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
  ebrd2_elevation <- terra::extract(elevation_raster, ebrd2_sf %>% as.data.frame %>% dplyr::select(longitude, latitude))
  ebrd2$elevation <- ebrd2_elevation$SRTM_global
  
  # Seasonal ebird data
  ebrd2_summer <- ebrd2 %>% filter(season == "summer")
  ebrd2_winter <- ebrd2 %>% filter(season == "winter")
  
  # Checklists per season per elevation bin
  ebrd2_summer_elevation_bins <- ebrd2_winter_elevation_bins <- list()
  for(i in 1:(length(elevation_bins)-1)){
    ebrd2_summer_elevation_bins[[i]] <- ebrd2_summer %>%
      subset(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
    ebrd2_winter_elevation_bins[[i]] <- ebrd2_winter %>%
      subset(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
  }
  
  # Incidence matrix for the whole mountain range
  incidence_mat <- table(ebrd2$scientificname, ebrd2$sampling_event_identifier)
  incidence_mat[which(incidence_mat > 1)] <- 1
  
  # Create incidence matrices for each elevation bin and each season
  # summer
  incidence_matrices_summer <- list()
  for(i in 1:(length(elevation_bins)-1)){
    if(nrow(ebrd2_summer_elevation_bins[[i]]) > 0){
      if(length(unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)) == 1){
        if(unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier) %in% colnames(incidence_mat)){
          incidence_matrices_summer[[i]] <- as.matrix(incidence_mat[,colnames(incidence_mat) %in% unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)])
          colnames(incidence_matrices_summer[[i]]) <- unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)
        }else{
          incidence_matrices_summer[[i]] <- NA
        }
      }else{
        incidence_matrices_summer[[i]] <- incidence_mat[,colnames(incidence_mat) %in% unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)]
      }
    }else{
      incidence_matrices_summer[[i]] <- NA
    }
    incidence_matrices_summer[[i]][which(incidence_matrices_summer[[i]] > 1)] <- 1
    if(sum(is.na(incidence_matrices_summer[[i]])) == 0){
      if(sum(apply(incidence_matrices_summer[[i]], 1, mean)) == 1){
        incidence_matrices_summer[[i]] <- NA
      }
    } 
  }
  # winter
  incidence_matrices_winter <- list()
  for(i in 1:(length(elevation_bins)-1)){
    if(nrow(ebrd2_winter_elevation_bins[[i]]) > 0){
      if(length(unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)) == 1){
        if(unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier) %in% colnames(incidence_mat)){
          incidence_matrices_winter[[i]] <- as.matrix(incidence_mat[,colnames(incidence_mat) %in% unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)])
          colnames(incidence_matrices_winter[[i]]) <- unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)
        }else{
          incidence_matrices_winter[[i]] <- NA
        }
      }else{
        incidence_matrices_winter[[i]] <- incidence_mat[,colnames(incidence_mat) %in% unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)]
      }
    }else{
      incidence_matrices_winter[[i]] <- NA
    }
    incidence_matrices_winter[[i]][which(incidence_matrices_winter[[i]] > 1)] <- 1
    if(sum(is.na(incidence_matrices_winter[[i]])) == 0){
      if(sum(apply(incidence_matrices_winter[[i]], 1, mean)) == 1){
        incidence_matrices_winter[[i]] <- NA
      }
    } 
  }
  
  # Calculate sample completeness
  names(incidence_matrices_summer) <- as.character(elevation_bins + 100)[-length(elevation_bins)]
  incidence_matrices_summer <- incidence_matrices_summer[names(incidence_matrices_summer) %in% names(unlist(lapply(incidence_matrices_summer, function(x) dim(x)[1])))]
  to_keep <- which(lapply(incidence_matrices_summer, ncol) >= 10)
  incidence_matrices_summer <- lapply(incidence_matrices_summer[to_keep], function(x) t(apply(x, 1, as.numeric)))
  names(incidence_matrices_winter) <- as.character(elevation_bins + 100)[-length(elevation_bins)]
  incidence_matrices_winter <- incidence_matrices_winter[names(incidence_matrices_winter) %in% names(unlist(lapply(incidence_matrices_winter, function(x) dim(x)[1])))]
  to_keep <- which(lapply(incidence_matrices_winter, ncol) >= 10)
  incidence_matrices_winter <- lapply(incidence_matrices_winter[to_keep], function(x) t(apply(x, 1, as.numeric)))
  
  if(length(incidence_matrices_summer) >= 7 & length(incidence_matrices_winter) >= 7){
    # summer
    out.raw.summer <- iNEXT(incidence_matrices_summer, datatype="incidence_raw")
    completeness_summer[[k]] <- out.raw.summer$AsyEst %>% subset(Diversity == "Species richness") %>% mutate(Season = "summer") %>% 
      dplyr::select(Season, Assemblage, Observed, Estimator, LCL, UCL) %>% 
      mutate(completeness = Observed/Estimator) %>%
      mutate(Assemblage = as.numeric(as.character(Assemblage))) %>%
      arrange(Assemblage)
    # winter
    out.raw.winter <- iNEXT(incidence_matrices_winter, datatype="incidence_raw")
    completeness_winter[[k]] <- out.raw.winter$AsyEst %>% subset(Diversity == "Species richness") %>% mutate(Season = "winter") %>% 
      dplyr::select(Season, Assemblage, Observed, Estimator, LCL, UCL) %>% 
      mutate(completeness = Observed/Estimator) %>%
      mutate(Assemblage = as.numeric(as.character(Assemblage))) %>%
      arrange(Assemblage)
  }
  print(k)
}


# Remove mountain slopes with less than 7 elevational bins with at least 10 checklists for at least one season
toremove <- c(2,18,29,30) # "Morocco","Khasi_Hills","Southern_Alps_West","Southern_Alps_East"
completeness_summer <- completeness_summer[-toremove]
completeness_winter <- completeness_winter[-toremove]
mountain_slopes_polygons <- mountain_slopes_polygons[-toremove]
slopes_mountain_name <- slopes_mountain_name[-toremove]
sc <- data.frame(slope = names(mountain_slopes_polygons),
                 completeness_summer = unlist(lapply(completeness_summer, function(x) median(x$completeness))),
                 completeness_winter = unlist(lapply(completeness_winter, function(x) median(x$completeness))))

# Remove mountains with at least one season with sample completeness < 0.7
toremove <- which(sc$completeness_summer < 0.7 | sc$completeness_winter < 0.7) # "Peru_Central", "Eastern_Himalayas"
mountain_slopes_polygons <- mountain_slopes_polygons[-toremove]
slopes_mountain_name <- slopes_mountain_name[-toremove]

ggplot() +
  geom_point(data=sc, aes(x=Assemblage, y=completeness, col=Season)) + 
  geom_line(data=sc, aes(x=Assemblage, y=completeness, col=Season)) + 
  ylim(0,1) + theme_classic() + xlab("Elevation") + ylab("Sampling completeness (obs richness / Chao2)")


# Put all the mountain slopes into one sf object 
mountain_slopes_polygons$Australian_Alps <- st_as_sfc(mountain_slopes_polygons$Australian_Alps)
mountain_slopes_polygons$Eastern_BC <- st_as_sfc(mountain_slopes_polygons$Eastern_BC)
mountain_slopes_polygons$Peru_SouthWest <- st_as_sfc(mountain_slopes_polygons$Peru_SouthWest)
mountain_slopes_polygons$Peru_East <- st_as_sfc(mountain_slopes_polygons$Peru_East)
mountain_slopes_polygons$Argentina_North <- st_as_sfc(mountain_slopes_polygons$Argentina_North)
mountain_slopes_polygons$Ecuador_West <- st_as_sfc(mountain_slopes_polygons$Ecuador_West)
mountain_slopes_polygons$Colombia_West <- st_as_sfc(mountain_slopes_polygons$Colombia_West)
mountain_slopes_polygons$Drakensberg <- st_as_sfc(mountain_slopes_polygons$Drakensberg)
mountain_slopes_polygons$Serra_Mantiqueira <- st_as_sfc(mountain_slopes_polygons$Serra_Mantiqueira)
mountain_slopes_polygons$Jalisco <- st_as_sfc(mountain_slopes_polygons$Jalisco)
mountain_slopes_polygons$Sierra_Occidental <- st_as_sfc(mountain_slopes_polygons$Sierra_Occidental)
mountain_slopes_polygons$S_Nevada <- st_as_sfc(mountain_slopes_polygons$S_Nevada)
mountain_slopes_polygons$Western_Himalayas <- st_as_sfc(mountain_slopes_polygons$Western_Himalayas)
mountain_slopes_polygons$Hawaii <- st_as_sfc(mountain_slopes_polygons$Hawaii)
mountain_slopes_polys = mountain_slopes_polygons[[1]]
for (i in 2:length(mountain_slopes_polygons)) {
  mountain_slopes_polys <- c(mountain_slopes_polys, mountain_slopes_polygons[[i]])
}
mountain_slopes_polys <- st_as_sf(mountain_slopes_polys) %>% mutate(name = names(mountain_slopes_polygons))

save(mountain_slopes_polys, slopes_mountain_name, ebrd2_slopes, file="results/mnt_slopes_data_new.RData")
