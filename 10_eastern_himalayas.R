## This script contains the code for predicting Eastern Himalayas ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

# Present climate

# Load environmental data prepared in 05_environmental_data.R
env_data_E_Him <- read_csv("env_data_E_Him.csv")[,-1]

# Load the mountain slopes polygons and ebird data prepared using the script 03_mountain_slopes_delineation.R
load("mountain_slopes.RData")

## Extract elevation for all eBird observation
ebrd2_slopes_E_Him <- ebrd2_slopes %>% filter(slope == "Eastern_Himalayas")
ebrd2_slopes_E_Him_sf <- st_as_sf(x = ebrd2_slopes_E_Him, coords = c("longitude", "latitude"), crs = crs(mountain_ranges), remove=F)
ebrd2_slopes_E_Him_elevation <- terra::extract(elevation_raster, ebrd2_slopes_E_Him_sf %>% as.data.frame %>% dplyr::select(longitude, latitude))
ebrd2_slopes_E_Him$elevation <- ebrd2_slopes_E_Him_elevation$SRTM_global


# Energy supply per elevational band
env_data_sub_summer <- env_data_E_Him %>% dplyr::filter(season == "summer")
env_data_sub_winter <- env_data_E_Him %>% dplyr::filter(season == "winter")
env_data_elevation_bins_summer <- env_data_elevation_bins_winter <- list()
for(i in 1:(length(elevation_bins)-1)){
  env_data_elevation_bins_summer[[i]] <- env_data_sub_summer %>%
    filter(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
  env_data_elevation_bins_winter[[i]] <- env_data_sub_winter %>%
    filter(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
}
env_data_elevation_bins_summer <- data.frame(slope = "Bhutan",
                                             season = "summer",
                                             elevation_bin = elevation_bins[-length(elevation_bins)] + 100,
                                             ndvi = unlist(lapply(env_data_elevation_bins_summer, function(x) mean(x$ndvi, na.rm=T))),
                                             temperature = unlist(lapply(env_data_elevation_bins_summer, function(x) mean(x$temperature, na.rm=T))),
                                             area = unlist(lapply(env_data_elevation_bins_summer, function(x) nrow(x))))
env_data_elevation_bins_winter <- data.frame(slope = "Bhutan",
                                             season = "winter",
                                             elevation_bin = elevation_bins[-length(elevation_bins)] + 100,
                                             ndvi = unlist(lapply(env_data_elevation_bins_winter, function(x) mean(x$ndvi, na.rm=T))),
                                             temperature = unlist(lapply(env_data_elevation_bins_winter, function(x) mean(x$temperature, na.rm=T))),
                                             area = unlist(lapply(env_data_elevation_bins_winter, function(x) nrow(x))))
env_data_elevation_bins_E_Him <- rbind(env_data_elevation_bins_summer, env_data_elevation_bins_winter)
env_data_elevation_bins_E_Him$ndvi[is.na(env_data_elevation_bins_E_Him$ndvi) == T] <- 0

ggplot() +
  geom_smooth(data=env_data_elevation_bins_E_Him, aes(x=elevation_bin, y=ndvi, col=season))


# Regional proportion of migrants 
Eastern_Himalayas_polygon <- mountain_slopes_polygons$Eastern_Himalayas
data_global_rich_E_Him <- data_global_rich[st_intersects(Eastern_Himalayas_polygon, data_global_rich_sf)[[1]],]
prop_long_dist_migr_E_Him = mean(data_global_rich_E_Him$DIFFERENCE / (data_global_rich_E_Him$RESIDENT + data_global_rich_E_Him$MIGRA_JAN))

## Mechanistic model
# Reproduction cost
lat_mnt <- Eastern_Himalayas_polygon %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y)
repro_cost = 0.41 #0.19 #0.41 # 0.1 # 1
reproduction_cost_summer = repro_cost 
reproduction_cost_winter = 0
# Long distance migration cost
migra_dist <- (lat_mnt$latitude - winter_peak_east_asia) * 111
migra_cost <- migra_dist * 6.45e-5
# Thermoregulation cost for long distance migrants
beta = 23.6 # 33.6 # 13.6
temp_peak_summer_east_asia <- 11.84833
temp_peak_winter_east_asia <- 21.29227
thermo_cost <- ifelse(temp_peak_winter_east_asia < (40 - beta), (40 - temp_peak_winter_east_asia - beta) / beta, 0)

# Prep environmental data
df <- env_data_elevation_bins_E_Him %>% mutate(energy_supply = area * ndvi) #%>% filter(is.na(energy_supply) == F)
energy_supply_summer <- df$energy_supply[df$season == "summer"]
energy_supply_winter <- df$energy_supply[df$season == "winter"]
total_energy_supply_summer <- sum(energy_supply_summer)
total_energy_supply_winter <- sum(energy_supply_winter)
# Thermoregulation cost  
thermoregulation_cost <- ifelse(df$temperature < (40 - beta), (40 - df$temperature - beta) / beta, 0)
thermoregulation_cost_summer <- thermoregulation_cost[df$season == "summer"]
thermoregulation_cost_winter <- thermoregulation_cost[df$season == "winter"]
# number of elevational bins
elevbins <- nrow(df) / 2
# Add long-distance migration
energy_supply_diff <- (total_energy_supply_summer / (1 + reproduction_cost_summer + mean(thermoregulation_cost_summer))) - (total_energy_supply_winter / (1 + reproduction_cost_winter + mean(thermoregulation_cost_winter)))
#energy_supply_winter_migration <- (total_energy_supply_summer / (1+reproduction_cost_summer)) - (total_energy_supply_winter / (1+reproduction_cost_winter))
energy_supply_winter <- c(energy_supply_winter, prop_long_dist_migr_E_Him*total_energy_supply_winter)
energy_supply_summer <- c(energy_supply_summer, 0)
migration_cost_winter = c(rep(0, elevbins), migra_cost) # 0.05 # 0.25
migration_cost_summer = c(rep(0, elevbins), 0) # 0.05 # 0.25
thermoregulation_cost_winter <- c(thermoregulation_cost_winter, thermo_cost)
thermoregulation_cost_summer <- c(thermoregulation_cost_summer, 0)
reproduction_cost_summer <- c(rep(reproduction_cost_summer, elevbins), repro_cost)
reproduction_cost_winter <- c(rep(reproduction_cost_winter, elevbins), 0)

# Simulate first population
energy_supply_yearound <- energy_efficiency_score <- list()
for(j in 1:length(energy_supply_winter)){
  energy_supply_yearound[[j]] <- energy_supply_summer + energy_supply_winter[j]
  energy_efficiency_score[[j]] <- (energy_supply_summer/(1 + reproduction_cost_summer + thermoregulation_cost_summer + migration_cost_summer)) + (energy_supply_winter[j] / (1 + reproduction_cost_winter + thermoregulation_cost_winter[j] + migration_cost_winter[j]))
}
energy_supply_yearound <- do.call(cbind, energy_supply_yearound)
energy_efficiency_score <- do.call(cbind, energy_efficiency_score)
energy_efficiency_score[nrow(energy_efficiency_score), ncol(energy_efficiency_score)] <- 0

selected_distribution <- vector()
# Run the model until 95% of the energy supply is consumed
while(sum(energy_supply_summer, na.rm=T) > total_energy_supply_summer*0.05 & sum(energy_supply_winter, na.rm=T) > total_energy_supply_winter*0.05){
  # Select most energy efficient population
  sel <- which(energy_efficiency_score == max(energy_efficiency_score, na.rm=T), arr.ind=T)
  if(nrow(sel) > 1){
    sel <- sel[sample(1:nrow(sel), 1),]
  }
  sel <- as.vector(sel)
  selected_distribution <- rbind(selected_distribution, sel)
  
  # Adjust energy supply
  energy_supply_summer[sel[1]] <- energy_supply_summer[sel[1]] - (1 + reproduction_cost_summer[sel[1]] + thermoregulation_cost_summer[sel[1]] + migration_cost_summer[sel[1]])
  energy_supply_winter[sel[2]] <- energy_supply_winter[sel[2]] - (1 + reproduction_cost_winter[sel[2]] + thermoregulation_cost_winter[sel[2]] + migration_cost_winter[sel[2]])
  
  # Calculate energy efficiency scores
  energy_supply_yearound <- energy_efficiency_score <- list()
  for(j in 1:length(energy_supply_winter)){
    energy_supply_yearound[[j]] <- energy_supply_summer + energy_supply_winter[j]
    energy_efficiency_score[[j]] <- (energy_supply_summer /(1 + reproduction_cost_summer + thermoregulation_cost_summer + migration_cost_summer)) + (energy_supply_winter[j] / (1 + reproduction_cost_winter[j] + thermoregulation_cost_winter[j] + migration_cost_winter[j]))
  }
  energy_supply_yearound <- do.call(cbind, energy_supply_yearound)
  energy_efficiency_score <- do.call(cbind, energy_efficiency_score)
  energy_efficiency_score[nrow(energy_efficiency_score), ncol(energy_efficiency_score)] <- 0
}
selected_distribution <- as.data.frame(selected_distribution)
colnames(selected_distribution) <- c("summer", "winter")
summer_richness_sim <- table(selected_distribution$summer)
winter_richness_sim <- table(selected_distribution$winter)
if(energy_supply_diff > 0){
  summer_richness_simu <- rep(0, sum(df$season == "summer"))
  winter_richness_simu <- c(rep(0, sum(df$season == "winter")), 0)
}else{
  summer_richness_simu <- c(rep(0, sum(df$season == "summer")), 0)
  winter_richness_simu <- rep(0, sum(df$season == "winter"))
}
summer_richness_simu[as.numeric(names(summer_richness_sim))] <- summer_richness_sim
winter_richness_simu[as.numeric(names(winter_richness_sim))] <- winter_richness_sim

# Estimate species richness from results based on the Species-Abundance Distribution — logseries 
species_from_individuals_logseries <- function(alpha, N){
  S = alpha * log(1 + (N/alpha))
  return(S)
}
summer_richness_simu_S <- species_from_individuals_logseries(100, summer_richness_simu)
winter_richness_simu_S <- species_from_individuals_logseries(100, winter_richness_simu)

# Prepare data to plot
SEDS_output_E_Him <- df %>% mutate(richness_simu = c(summer_richness_simu_S[1:(nrow(df)/2)], winter_richness_simu_S[1:(nrow(df)/2)]))
richness_E_Him_df <- SEDS_output_E_Him %>% filter(season == "summer") %>% dplyr::select(richness_simu) %>% rename(richness_simu_summer = richness_simu) %>%
  mutate(richness_simu_winter = SEDS_output_E_Him$richness_simu[SEDS_output_E_Him$season == "winter"]) %>%
  mutate(elevation = SEDS_output_E_Him$elevation_bin[SEDS_output_E_Him$season == "winter"]) %>%
  mutate(latitude = lat_mnt$latitude) %>%
  mutate(richness_difference_simu = richness_simu_summer - richness_simu_winter)



###  Future climate  ###

env_data_future_E_Him <- read_csv("results/env_data_future_E_Him.csv")[,-1]
# Energy supply per elevational band
env_data_sub_summer <- env_data_future_E_Him %>% dplyr::filter(season == "summer")
env_data_sub_winter <- env_data_future_E_Him %>% dplyr::filter(season == "winter")
env_data_elevation_bins_summer <- env_data_elevation_bins_winter <- list()
for(i in 1:(length(elevation_bins)-1)){
  env_data_elevation_bins_summer[[i]] <- env_data_sub_summer %>%
    filter(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
  env_data_elevation_bins_winter[[i]] <- env_data_sub_winter %>%
    filter(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
}
env_data_elevation_bins_summer <- data.frame(slope = "Bhutan",
                                             season = "summer",
                                             elevation_bin = elevation_bins[-length(elevation_bins)] + 100,
                                             ndvi = unlist(lapply(env_data_elevation_bins_summer, function(x) mean(x$ndvi, na.rm=T))),
                                             temperature = unlist(lapply(env_data_elevation_bins_summer, function(x) mean(x$temperature, na.rm=T))),
                                             area = unlist(lapply(env_data_elevation_bins_summer, function(x) nrow(x))))
env_data_elevation_bins_winter <- data.frame(slope = "Bhutan",
                                             season = "winter",
                                             elevation_bin = elevation_bins[-length(elevation_bins)] + 100,
                                             ndvi = unlist(lapply(env_data_elevation_bins_winter, function(x) mean(x$ndvi, na.rm=T))),
                                             temperature = unlist(lapply(env_data_elevation_bins_winter, function(x) mean(x$temperature, na.rm=T))),
                                             area = unlist(lapply(env_data_elevation_bins_winter, function(x) nrow(x))))
env_data_elevation_bins_E_Him <- rbind(env_data_elevation_bins_summer, env_data_elevation_bins_winter)
env_data_elevation_bins_E_Him$ndvi[is.na(env_data_elevation_bins_E_Him$ndvi) == T] <- 0

ggplot() +
  geom_smooth(data=env_data_elevation_bins_E_Him, aes(x=elevation_bin, y=ndvi, col=season))


# Regional proportion of migrants 
data_global_rich_E_Him <- data_global_rich[st_intersects(Eastern_Himalayas_polygon, data_global_rich_sf)[[1]],]
prop_long_dist_migr_E_Him = mean(data_global_rich_E_Him$DIFFERENCE / (data_global_rich_E_Him$RESIDENT + data_global_rich_E_Him$MIGRA_JAN))

## Mechanistic model
# Reproduction cost
lat_mnt <- Eastern_Himalayas_polygon %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y)
repro_cost = 0.41 #0.19 #0.41 # 0.1 # 1
reproduction_cost_summer = repro_cost 
reproduction_cost_winter = 0
# Long distance migration cost
migra_dist <- (lat_mnt$latitude - winter_peak_east_asia) * 111
migra_cost <- migra_dist * 6.45e-5
# Thermoregulation cost for long distance migrants
beta = 23.6 # 33.6 # 13.6
temp_peak_summer_east_asia <- 11.84833
temp_peak_winter_east_asia <- 21.29227
thermo_cost <- ifelse(temp_peak_winter_east_asia < (40 - beta), (40 - temp_peak_winter_east_asia - beta) / beta, 0)

# Prep environmental data
df <- env_data_elevation_bins_E_Him %>% mutate(energy_supply = area * ndvi) #%>% filter(is.na(energy_supply) == F)
energy_supply_summer <- df$energy_supply[df$season == "summer"]
energy_supply_winter <- df$energy_supply[df$season == "winter"]
total_energy_supply_summer <- sum(energy_supply_summer)
total_energy_supply_winter <- sum(energy_supply_winter)
# Thermoregulation cost  
thermoregulation_cost <- ifelse(df$temperature < (40 - beta), (40 - df$temperature - beta) / beta, 0)
thermoregulation_cost_summer <- thermoregulation_cost[df$season == "summer"]
thermoregulation_cost_winter <- thermoregulation_cost[df$season == "winter"]
# number of elevational bins
elevbins <- nrow(df) / 2
# Add long-distance migration
#energy_supply_diff <- (total_energy_supply_summer / (1 + reproduction_cost_summer + mean(thermoregulation_cost_summer))) - (total_energy_supply_winter / (1 + reproduction_cost_winter + mean(thermoregulation_cost_winter)))
#energy_supply_winter_migration <- (total_energy_supply_summer / (1+reproduction_cost_summer)) - (total_energy_supply_winter / (1+reproduction_cost_winter))
energy_supply_winter <- c(energy_supply_winter, prop_long_dist_migr_E_Him*total_energy_supply_winter)
energy_supply_summer <- c(energy_supply_summer, 0)
migration_cost_winter = c(rep(0, elevbins), migra_cost) # 0.05 # 0.25
migration_cost_summer = c(rep(0, elevbins), 0) # 0.05 # 0.25
thermoregulation_cost_winter <- c(thermoregulation_cost_winter, thermo_cost)
thermoregulation_cost_summer <- c(thermoregulation_cost_summer, 0)
reproduction_cost_summer <- c(rep(reproduction_cost_summer, elevbins), repro_cost)
reproduction_cost_winter <- c(rep(reproduction_cost_winter, elevbins), 0)

# Simulate first population
energy_supply_yearound <- energy_efficiency_score <- list()
for(j in 1:length(energy_supply_winter)){
  energy_supply_yearound[[j]] <- energy_supply_summer + energy_supply_winter[j]
  energy_efficiency_score[[j]] <- (energy_supply_summer/(1 + reproduction_cost_summer + thermoregulation_cost_summer + migration_cost_summer)) + (energy_supply_winter[j] / (1 + reproduction_cost_winter + thermoregulation_cost_winter[j] + migration_cost_winter[j]))
}
energy_supply_yearound <- do.call(cbind, energy_supply_yearound)
energy_efficiency_score <- do.call(cbind, energy_efficiency_score)
energy_efficiency_score[nrow(energy_efficiency_score), ncol(energy_efficiency_score)] <- 0

selected_distribution <- vector()
# Run the model until 95% of the energy supply is consumed
while(sum(energy_supply_summer, na.rm=T) > total_energy_supply_summer*0.05 & sum(energy_supply_winter, na.rm=T) > total_energy_supply_winter*0.05){
  # Select most energy efficient population
  sel <- which(energy_efficiency_score == max(energy_efficiency_score, na.rm=T), arr.ind=T)
  if(nrow(sel) > 1){
    sel <- sel[sample(1:nrow(sel), 1),]
  }
  sel <- as.vector(sel)
  selected_distribution <- rbind(selected_distribution, sel)
  
  # Adjust energy supply
  energy_supply_summer[sel[1]] <- energy_supply_summer[sel[1]] - (1 + reproduction_cost_summer[sel[1]] + thermoregulation_cost_summer[sel[1]] + migration_cost_summer[sel[1]])
  energy_supply_winter[sel[2]] <- energy_supply_winter[sel[2]] - (1 + reproduction_cost_winter[sel[2]] + thermoregulation_cost_winter[sel[2]] + migration_cost_winter[sel[2]])
  
  # Calculate energy efficiency scores
  energy_supply_yearound <- energy_efficiency_score <- list()
  for(j in 1:length(energy_supply_winter)){
    energy_supply_yearound[[j]] <- energy_supply_summer + energy_supply_winter[j]
    energy_efficiency_score[[j]] <- (energy_supply_summer /(1 + reproduction_cost_summer + thermoregulation_cost_summer + migration_cost_summer)) + (energy_supply_winter[j] / (1 + reproduction_cost_winter[j] + thermoregulation_cost_winter[j] + migration_cost_winter[j]))
  }
  energy_supply_yearound <- do.call(cbind, energy_supply_yearound)
  energy_efficiency_score <- do.call(cbind, energy_efficiency_score)
  energy_efficiency_score[nrow(energy_efficiency_score), ncol(energy_efficiency_score)] <- 0
}
selected_distribution <- as.data.frame(selected_distribution)
colnames(selected_distribution) <- c("summer", "winter")
summer_richness_sim <- table(selected_distribution$summer)
winter_richness_sim <- table(selected_distribution$winter)
if(energy_supply_diff > 0){
  summer_richness_simu <- rep(0, sum(df$season == "summer"))
  winter_richness_simu <- c(rep(0, sum(df$season == "winter")), 0)
}else{
  summer_richness_simu <- c(rep(0, sum(df$season == "summer")), 0)
  winter_richness_simu <- rep(0, sum(df$season == "winter"))
}
summer_richness_simu[as.numeric(names(summer_richness_sim))] <- summer_richness_sim
winter_richness_simu[as.numeric(names(winter_richness_sim))] <- winter_richness_sim

# Estimate species richness from results based on the Species-Abundance Distribution — logseries 
species_from_individuals_logseries <- function(alpha, N){
  S = alpha * log(1 + (N/alpha))
  return(S)
}
summer_richness_simu_S <- species_from_individuals_logseries(100, summer_richness_simu)
winter_richness_simu_S <- species_from_individuals_logseries(100, winter_richness_simu)

# Prepare data to plot
SEDS_output_E_Him_future <- df %>% mutate(richness_simu = c(summer_richness_simu_S[1:(nrow(df)/2)], winter_richness_simu_S[1:(nrow(df)/2)]))
richness_E_Him_df_future <- SEDS_output_E_Him_future %>% filter(season == "summer") %>% dplyr::select(richness_simu) %>% rename(richness_simu_summer = richness_simu) %>%
  mutate(richness_simu_winter = SEDS_output_E_Him_future$richness_simu[SEDS_output_E_Him_future$season == "winter"]) %>%
  mutate(elevation = SEDS_output_E_Him_future$elevation_bin[SEDS_output_E_Him_future$season == "winter"]) %>%
  mutate(latitude = lat_mnt$latitude) %>%
  mutate(richness_difference_simu = richness_simu_summer - richness_simu_winter)


g_E_Him <- ggplot() +
  geom_hline(yintercept = 0, col="grey50") +
  geom_point(data=richness_E_Him_df, aes(x=elevation, y=richness_difference_simu), col="black", size=0.8) +
  geom_smooth(data=richness_E_Him_df, aes(x=elevation, y=richness_difference_simu), col="black", linewidth = 1, se=T, method="loess", alpha=0.1, span=0.5) +
  geom_point(data=richness_E_Him_df_future, aes(x=elevation, y=richness_difference_simu), col="purple", size=0.8) +
  geom_smooth(data=richness_E_Him_df_future, aes(x=elevation, y=richness_difference_simu), col="purple", linewidth = 1, se=T, method="loess", alpha=0.1, span=0.5) +
  xlab("Elevation") + ylab("Seasonal difference in richness") + theme_classic() + labs(tag="B")
