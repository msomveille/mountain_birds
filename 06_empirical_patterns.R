## This script contains the code for computing and plotting the empirical patterns ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ncdf4)
library(raster)
library(Hmisc)
library(rworldmap)

# Get environmental data for the present prepared in the script 05_environmental_data.R
env_data <- read.csv("env_data.csv")

# Load the data prepared in the script 04_sample_completeness.R
load("mnt_slopes_data_new.RData")
mountain_ranges <- vect("mountain_ranges_2.shp")


elevation_bins <- seq(0, 7000, 200)
slopes <- mountain_slopes_polys$name

## Species elevational ranges and climatic niches ##

species_elevational_range_winter <- species_elevational_range_summer <- list()
for(k in 1:nrow(mountain_slopes_polys)){
  
  mountain_range <- mountain_ranges[slopes_mountain_name[k],]
  
  # Load ebird data for the mountain range
  ebrd2 <- read_csv(paste0("ebird_data/", mountain_ranges$NAME[slopes_mountain_name[k]],".csv")) 
  
  # filter ebird data for checklist that are within the 2 degrees gradient
  ebrd2_slopes_sub <- ebrd2_slopes %>% filter(slope == slopes[k])
  ebrd2 <- ebrd2 %>% filter(sampling_event_identifier %in% ebrd2_slopes_sub$sampling_event_identifier)
  
  ebrd2_sf <- st_as_sf(x = ebrd2, coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
  ebrd2_elevation <- terra::extract(elevation_raster, ebrd2_sf %>% as.data.frame %>% dplyr::select(longitude, latitude))
  ebrd2$elevation <- ebrd2_elevation$SRTM_global
  
  # Seasonal ebird data
  ebrd2_summer <- ebrd2 %>% filter(season == "summer")
  ebrd2_winter <- ebrd2 %>% filter(season == "winter")
  
  # filter environmental data
  env_data_sub <- env_data %>% filter(slope == slopes[k])
  
  # Seasonal incidence matrix for the whole mountain range
  incidence_mat_w <- table(ebrd2_winter$scientificname, ebrd2_winter$sampling_event_identifier) # winter incidence matrix
  incidence_mat_w[which(incidence_mat_w > 1)] <- 1
  incidence_mat_s <- table(ebrd2_summer$scientificname, ebrd2_summer$sampling_event_identifier) # summer incidence matrix
  incidence_mat_s[which(incidence_mat_s > 1)] <- 1
  
  # Remove species observed in less than 20 checklist over the entire mountain range in a given season
  species_incidence_w <- apply(incidence_mat_w, 1, function(x) which(x == 1))
  species_selected_w <- which(unlist(lapply(species_incidence_w, length)) >= 20)
  incidence_mat_w <- incidence_mat_w[species_selected_w,]
  species_incidence_s <- apply(incidence_mat_s, 1, function(x) which(x == 1))
  species_selected_s <- which(unlist(lapply(species_incidence_s, length)) >= 20)
  incidence_mat_s <- incidence_mat_s[species_selected_s,]
  
  # Checklists per season per elevation bin
  ebrd2_summer_elevation_bins <- ebrd2_winter_elevation_bins <- list()
  for(i in 1:(length(elevation_bins)-1)){
    ebrd2_summer_elevation_bins[[i]] <- ebrd2_summer %>%
      subset(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
    ebrd2_winter_elevation_bins[[i]] <- ebrd2_winter %>%
      subset(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
  }
  
  # Summer
  # Create incidence matrices for each elevation bin and each season
  incidence_matrices_summer <- list()
  for(i in 1:(length(elevation_bins)-1)){
    if(nrow(ebrd2_summer_elevation_bins[[i]]) > 0){
      if(length(unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)) == 1){
        if(unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier) %in% colnames(incidence_mat_s)){
          incidence_matrices_summer[[i]] <- as.matrix(incidence_mat_s[,colnames(incidence_mat_s) %in% unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)])
          colnames(incidence_matrices_summer[[i]]) <- unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)
        }else{
          incidence_matrices_summer[[i]] <- NA
        }
      }else{
        incidence_matrices_summer[[i]] <- incidence_mat_s[,colnames(incidence_mat_s) %in% unique(ebrd2_summer_elevation_bins[[i]]$sampling_event_identifier)]
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
  # Estimate species elevational ranges
  names(incidence_matrices_summer) <- elevation_bins[1:(length(elevation_bins)-1)] + 100
  incidence_matrices_summer <- incidence_matrices_summer[which(is.na(incidence_matrices_summer)==F)]
  species_fraction_checklist_summer <- vector()
  for(i in 1:length(incidence_matrices_summer)){
    species_fraction_checklist_summer <- cbind(species_fraction_checklist_summer, apply(incidence_matrices_summer[[i]], 1, sum) / ncol(incidence_matrices_summer[[i]]))
  }
  species_fraction_checklist_summer <- t(apply(species_fraction_checklist_summer, 1, function(x) x/max(x)))
  ranges.min <- apply(species_fraction_checklist_summer, 1, function(x) as.numeric(names(incidence_matrices_summer))[min(which(x > 0.05))] - 100)
  ranges.max <- apply(species_fraction_checklist_summer, 1, function(x) as.numeric(names(incidence_matrices_summer))[max(which(x > 0.05))] + 100)
  ranges.mean <- apply(species_fraction_checklist_summer, 1, function(x) wtd.mean(as.numeric(names(incidence_matrices_summer)), x))
  ranges.sd <- apply(species_fraction_checklist_summer, 1, function(x) sqrt(wtd.var(as.numeric(names(incidence_matrices_summer)), x)))
  ranges.sd[which(is.na(ranges.sd) == T)] <- 100
  # Estimate summer realised thermal niche
  env_data_sub_summer <- env_data_sub %>% filter(season == "summer")
  temperature_elevation_bin <- vector()
  for(i in 1:length(as.numeric(names(incidence_matrices_summer)))){
    temperature_elevation_bin[i] <- mean(env_data_sub_summer %>% filter(env_data_sub_summer$elevation >= as.numeric(names(incidence_matrices_summer))[i]-100 & env_data_sub_summer$elevation < as.numeric(names(incidence_matrices_summer))[i]+100) %>% dplyr::select(temperature) %>% as.vector() %>% unlist(), na.rm=T)
  }
  temp.mean <- apply(species_fraction_checklist_summer, 1, function(x) wtd.mean(temperature_elevation_bin, x))
  temp.sd <- apply(species_fraction_checklist_summer, 1, function(x) sqrt(wtd.var(temperature_elevation_bin, x)))
  species_elevational_range_summer[[k]] <- data.frame(slope = slopes[k],
                                                      season = "summer",
                                                      species = names(ranges.min),
                                                      range_min = ranges.min,
                                                      range_max = ranges.max,
                                                      range_mean = ranges.mean,
                                                      range_sd = ranges.sd,
                                                      temp_mean = temp.mean,
                                                      temp_sd = temp.sd)
  rm(incidence_matrices_summer, incidence_mat_s)
  
  # Winter
  # Create incidence matrices for each elevation bin and each season
  incidence_matrices_winter <- list()
  for(i in 1:(length(elevation_bins)-1)){
    if(nrow(ebrd2_winter_elevation_bins[[i]]) > 0){
      if(length(unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)) == 1){
        if(unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier) %in% colnames(incidence_mat_w)){
          incidence_matrices_winter[[i]] <- as.matrix(incidence_mat_w[,colnames(incidence_mat_w) %in% unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)])
          colnames(incidence_matrices_winter[[i]]) <- unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)
        }else{
          incidence_matrices_winter[[i]] <- NA
        }
      }else{
        incidence_matrices_winter[[i]] <- incidence_mat_w[,colnames(incidence_mat_w) %in% unique(ebrd2_winter_elevation_bins[[i]]$sampling_event_identifier)]
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
  # Estimate species elevational ranges
  names(incidence_matrices_winter) <- elevation_bins[1:(length(elevation_bins)-1)] + 100
  incidence_matrices_winter <- incidence_matrices_winter[which(is.na(incidence_matrices_winter)==F)]
  species_fraction_checklist_winter <- vector()
  for(i in 1:length(incidence_matrices_winter)){
    species_fraction_checklist_winter <- cbind(species_fraction_checklist_winter, apply(incidence_matrices_winter[[i]], 1, sum) / ncol(incidence_matrices_winter[[i]]))
  }
  species_fraction_checklist_winter <- t(apply(species_fraction_checklist_winter, 1, function(x) x/max(x)))
  ranges.min <- apply(species_fraction_checklist_winter, 1, function(x) as.numeric(names(incidence_matrices_winter))[min(which(x > 0.05))] - 100)
  ranges.max <- apply(species_fraction_checklist_winter, 1, function(x) as.numeric(names(incidence_matrices_winter))[max(which(x > 0.05))] + 100)
  ranges.mean <- apply(species_fraction_checklist_winter, 1, function(x) wtd.mean(as.numeric(names(incidence_matrices_winter)), x))
  ranges.sd <- apply(species_fraction_checklist_winter, 1, function(x) sqrt(wtd.var(as.numeric(names(incidence_matrices_winter)), x)))
  ranges.sd[which(is.na(ranges.sd) == T)] <- 100
  # Estimate wintering realised thermal niche
  env_data_sub_winter <- env_data_sub %>% filter(season == "winter")
  temperature_elevation_bin <- vector()
  for(i in 1:length(as.numeric(names(incidence_matrices_winter)))){
    temperature_elevation_bin[i] <- mean(env_data_sub_winter %>% filter(env_data_sub_winter$elevation >= as.numeric(names(incidence_matrices_winter))[i]-100 & env_data_sub_winter$elevation < as.numeric(names(incidence_matrices_winter))[i]+100) %>% dplyr::select(temperature) %>% as.vector() %>% unlist(), na.rm=T)
  }
  temp.mean <- apply(species_fraction_checklist_winter, 1, function(x) wtd.mean(temperature_elevation_bin, x))
  temp.sd <- apply(species_fraction_checklist_winter, 1, function(x) sqrt(wtd.var(temperature_elevation_bin, x)))
  species_elevational_range_winter[[k]] <- data.frame(slope = slopes[k],
                                                      season = "winter",
                                                      species = names(ranges.min),
                                                      range_min = ranges.min,
                                                      range_max = ranges.max,
                                                      range_mean = ranges.mean,
                                                      range_sd = ranges.sd,
                                                      temp_mean = temp.mean,
                                                      temp_sd = temp.sd)
  
  rm(incidence_matrices_winter, incidence_mat_w)
  
  print(k) 
}

species_elevational_range_summer <- do.call(rbind, species_elevational_range_summer)
species_elevational_range_winter <- do.call(rbind, species_elevational_range_winter)
species_elevational_ranges <- rbind(species_elevational_range_summer, species_elevational_range_winter)
write.csv(species_elevational_ranges, "species_elevational_ranges.csv")




############################################
## Plot Taiwan Yuhina for illustration
############################################

# Run the loop above for k = 32 and without removing incidence_matrices_summer and incidence_matrices_winter

aa <- data.frame(elev = as.numeric(c(names(incidence_matrices_summer), names(incidence_matrices_winter))),
                 freq = c(species_fraction_checklist_summer[which(rownames(species_fraction_checklist_summer) == "Yuhina brunneiceps"),], species_fraction_checklist_winter[which(rownames(species_fraction_checklist_winter) == "Yuhina brunneiceps"),]),
                 season = c(rep("summer", ncol(species_fraction_checklist_summer)), rep("winter", ncol(species_fraction_checklist_winter))))
weighted.mean(aa$elev[aa$season == "summer"], aa$freq[aa$season == "summer"])
weighted.mean(aa$elev[aa$season == "winter"], aa$freq[aa$season == "winter"])
ggplot() +
  geom_point(data=aa, aes(x=freq, y=elev, col=season)) +
  theme_bw() + geom_vline(xintercept=0.05) + xlab("Rescaled frequency of observation") + ylab("Elevation") + theme(legend.position = "bottom") +
  geom_vline(xintercept=1.1, col="red") 

################################################



# Plot species' seasonal elevational ranges for each mountain range
k=0
k=k+1
species_elevational_ranges_sub <- species_elevational_ranges %>% filter(slope == slopes[k])
pdf(file=paste0("results/figures/elevational_ranges/", slopes[k], ".pdf"), width = 8, height = 8)
ggplot(species_elevational_ranges_sub, aes(x = reorder(species, -range_mean), y = range_mean)) + ylab("Elevation") + xlab("Species") +
  geom_pointrange(aes(ymin = range_min, ymax = range_max, colour = season), fatten = 1, size=1.3, alpha=0.7, position = position_dodge2(width=0.7)) + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=20), axis.text.x=element_blank())
dev.off()


##  FIGURE 1  ##

newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, '+proj=longlat +datum=WGS84')
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- vect(newmap)
newmap <- terra::buffer(newmap, 0)
newmap <- terra::aggregate(newmap)
newmap <- terra::crop(newmap, ext(-180,180,-60,90))
g_slopes <- ggplot() +
  geom_sf(data=newmap, col="grey80", fill="grey80") +
  geom_sf(data = mountain_slopes_polys, col="black", fill="black")
pdf(file = paste0("results/figures/mountain_slopes_map.pdf"), width = 12, height = 7)
g_slopes
dev.off()








