## This script contains the code for predicting future patterns under climate warming ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)


# Load environmental data prepared in 05_environmental_data.R
env_data_future <- read_csv("env_data_future.csv")[,-1]

# Load empirical patterns in the seasonal difference in richness
richness_all_df <- read_csv("richness_all_df.csv")

slopes <- unique(richness_all_df$slope)
elevation_bins <- seq(0, 7000, 200)

# Energy supply per elevational band
env_data_elevation_bins_future <- list()
for(k in 1:length(slopes)){
  env_data_sub_summer <- env_data_future %>% dplyr::filter(slope == slopes[k] & season == "summer")
  env_data_sub_winter <- env_data_future %>% dplyr::filter(slope == slopes[k] & season == "winter")
  env_data_elevation_bins_summer <- env_data_elevation_bins_winter <- list()
  for(i in 1:(length(elevation_bins)-1)){
    env_data_elevation_bins_summer[[i]] <- env_data_sub_summer %>%
      filter(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
    env_data_elevation_bins_winter[[i]] <- env_data_sub_winter %>%
      filter(elevation > elevation_bins[i] & elevation <= elevation_bins[i+1]) 
  }
  env_data_elevation_bins_summer <- data.frame(slope = slopes[k],
                                               season = "summer",
                                               elevation_bin = elevation_bins[-length(elevation_bins)] + 100,
                                               ndvi = unlist(lapply(env_data_elevation_bins_summer, function(x) mean(x$ndvi, na.rm=T))),
                                               temperature = unlist(lapply(env_data_elevation_bins_summer, function(x) mean(x$temperature, na.rm=T))),
                                               area = unlist(lapply(env_data_elevation_bins_summer, function(x) nrow(x))))
  env_data_elevation_bins_winter <- data.frame(slope = slopes[k],
                                               season = "winter",
                                               elevation_bin = elevation_bins[-length(elevation_bins)] + 100,
                                               ndvi = unlist(lapply(env_data_elevation_bins_winter, function(x) mean(x$ndvi, na.rm=T))),
                                               temperature = unlist(lapply(env_data_elevation_bins_winter, function(x) mean(x$temperature, na.rm=T))),
                                               area = unlist(lapply(env_data_elevation_bins_winter, function(x) nrow(x))))
  env_data_elevation_bins_future[[k]] <- rbind(env_data_elevation_bins_summer, env_data_elevation_bins_winter)
}
env_data_elevation_bins_future <- do.call(rbind, env_data_elevation_bins_future)
env_data_elevation_bins_future$ndvi[is.na(env_data_elevation_bins_future$ndvi) == T] <- 0
env_data_elevation_bins_future <- env_data_elevation_bins_future[is.na(env_data_elevation_bins_future$temperature) == F,]

ggplot() +
  geom_smooth(data=env_data_elevation_bins %>% filter(slope == "Swiss_Alps"), aes(x=elevation_bin, y=temperature, col=season)) +
  geom_smooth(data=env_data_elevation_bins_future %>% filter(slope == "Swiss_Alps"), aes(x=elevation_bin, y=temperature, col=season))


## Mechanistic model

# Reproduction cost
lat_mnt <- mountain_slopes_polys %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y)
reproduction_proba_summer <- rep(0, length(slopes))
reproduction_proba_summer[which(lat_mnt > 15)] <- 1
reproduction_proba_summer[which(lat_mnt > 0 & lat_mnt < 15)] <- 0.75
reproduction_proba_summer[which(lat_mnt > -15 & lat_mnt < 0)] <- 0.25
repro_cost = 0.41 #0.19 #0.41 # 0.1 # 1
repro_cost_summer = repro_cost * reproduction_proba_summer 
repro_cost_winter = repro_cost * (1-reproduction_proba_summer)

# Long distance migration cost
migra_dist <- rep(0, length(slopes))
migra_dist[prop_long_dist_migr > 0 & lat_mnt$latitude > 0] <- (lat_mnt$latitude[prop_long_dist_migr > 0 & lat_mnt$latitude > 0] - winter_peak_global[prop_long_dist_migr > 0 & lat_mnt$latitude > 0]) * 111
migra_dist[prop_long_dist_migr < 0 & lat_mnt$latitude > 0] <- abs(lat_mnt$latitude[prop_long_dist_migr < 0 & lat_mnt$latitude > 0] - summer_peak_global[prop_long_dist_migr < 0 & lat_mnt$latitude > 0]) * 111
migra_dist[lat_mnt$latitude < 0] <- 3000
migra_cost <- migra_dist * 6.45e-5

# Thermoregulation cost for long distance migrants
beta = 23.6 # 33.6 # 13.6
temp_peak_summer_americas <- 17.71086 #12.72042
temp_peak_winter_americas <- 20.29994 #18.32428
temp_peak_summer_euro_africa <- 15.53568 #11.65147
temp_peak_winter_euro_africa <- 28.86751 #26.66647
temp_peak_summer_east_asia <- 16.13605 #11.84833
temp_peak_winter_east_asia <- 24.29982 #21.29227
thermo_cost <- rep(0, length(slopes))
thermo_cost[prop_long_dist_migr > 0 & lat_mnt$latitude > 0 & flyway=="americas"] <- ifelse(temp_peak_winter_americas < (40 - beta), (40 - temp_peak_winter_americas - beta) / beta, 0)
thermo_cost[prop_long_dist_migr > 0 & lat_mnt$latitude > 0 & flyway=="euro_africa"] <- ifelse(temp_peak_winter_euro_africa < (40 - beta), (40 - temp_peak_winter_euro_africa - beta) / beta, 0)
thermo_cost[prop_long_dist_migr > 0 & lat_mnt$latitude > 0 & flyway=="east_asia"] <- ifelse(temp_peak_winter_east_asia < (40 - beta), (40 - temp_peak_winter_east_asia - beta) / beta, 0)
thermo_cost[prop_long_dist_migr < 0 & lat_mnt$latitude > 0 & flyway=="americas"] <- ifelse(temp_peak_summer_americas < (40 - beta), (40 - temp_peak_summer_americas - beta) / beta, 0)
thermo_cost[prop_long_dist_migr < 0 & lat_mnt$latitude > 0 & flyway=="euro_africa"] <- ifelse(temp_peak_summer_euro_africa < (40 - beta), (40 - temp_peak_summer_euro_africa - beta) / beta, 0)
thermo_cost[prop_long_dist_migr < 0 & lat_mnt$latitude > 0 & flyway=="east_asia"] <- ifelse(temp_peak_summer_east_asia < (40 - beta), (40 - temp_peak_summer_east_asia - beta) / beta, 0)
thermo_cost[lat_mnt$latitude < 0] <- ifelse(20 < (40 - beta), (40 - 20 - beta) / beta, 0)

SEDS_output_future <- list()
for(k in 1:length(slopes)){
  # Prep environmental data
  df <- env_data_elevation_bins_future %>% filter(slope == slopes[k]) %>%
    mutate(energy_supply = area * ndvi) #%>% filter(is.na(energy_supply) == F)
  energy_supply_summer <- df$energy_supply[df$season == "summer"]
  energy_supply_winter <- df$energy_supply[df$season == "winter"]
  total_energy_supply_summer <- sum(energy_supply_summer, na.rm=T)
  total_energy_supply_winter <- sum(energy_supply_winter, na.rm=T)
  
  print(table(df$season)) # check that there are the same number of elevational bins at both seasons
  
  # Reproduction cost
  reproduction_cost_summer <- repro_cost_summer[k]
  reproduction_cost_winter <- repro_cost_winter[k]
  
  # Thermoregulation cost  
  thermoregulation_cost <- ifelse(df$temperature < (40 - beta), (40 - df$temperature - beta) / beta, 0)
  thermoregulation_cost_summer <- thermoregulation_cost[df$season == "summer"]
  thermoregulation_cost_winter <- thermoregulation_cost[df$season == "winter"]
  
  # number of elevational bins
  elevbins <- nrow(df) / 2
  
  # Add long-distance migration
  #energy_supply_diff <- (total_energy_supply_summer / (1 + reproduction_cost_summer + mean(thermoregulation_cost_summer))) - (total_energy_supply_winter / (1 + reproduction_cost_winter + mean(thermoregulation_cost_winter)))
  #energy_supply_winter_migration <- (total_energy_supply_summer / (1+reproduction_cost_summer)) - (total_energy_supply_winter / (1+reproduction_cost_winter))
  if(prop_long_dist_migr[k] > 0){
    energy_supply_winter <- c(energy_supply_winter, prop_long_dist_migr[k]*total_energy_supply_winter)
    energy_supply_summer <- c(energy_supply_summer, 0)
    migration_cost_winter = c(rep(0, elevbins), migra_cost[k]) # 0.05 # 0.25
    migration_cost_summer = c(rep(0, elevbins), 0) # 0.05 # 0.25
    thermoregulation_cost_winter <- c(thermoregulation_cost_winter, thermo_cost[k])
    thermoregulation_cost_summer <- c(thermoregulation_cost_summer, 0)
  }else{
    energy_supply_winter <- c(energy_supply_winter, 0)
    energy_supply_summer <- c(energy_supply_summer, abs(prop_long_dist_migr[k])*total_energy_supply_summer)
    migration_cost_winter = c(rep(0, elevbins), 0) # 0.05 # 0.25
    migration_cost_summer = c(rep(0, elevbins), migra_cost[k]) # 0.05 # 0.25
    thermoregulation_cost_winter <- c(thermoregulation_cost_winter, 0)
    thermoregulation_cost_summer <- c(thermoregulation_cost_summer, thermo_cost[k])
  }
  
  if(lat_mnt$latitude[k] > 0){
    reproduction_cost_summer <- c(rep(reproduction_cost_summer, elevbins), repro_cost)
    reproduction_cost_winter <- c(rep(reproduction_cost_winter, elevbins), 0)
  }else{
    reproduction_cost_summer <- c(rep(reproduction_cost_summer, elevbins), reproduction_cost_summer)
    reproduction_cost_winter <- c(rep(reproduction_cost_winter, elevbins), reproduction_cost_winter)
  }
  
  # Simulate first population
  energy_supply_yearound <- energy_efficiency_score <- list()
  for(j in 1:length(energy_supply_winter)){
    energy_supply_yearound[[j]] <- energy_supply_summer + energy_supply_winter[j]
    energy_efficiency_score[[j]] <- (energy_supply_summer/(1 + reproduction_cost_summer + thermoregulation_cost_summer + migration_cost_summer)) + (energy_supply_winter[j] / (1 + reproduction_cost_winter[j] + thermoregulation_cost_winter[j] + migration_cost_winter[j]))
  }
  energy_supply_yearound <- do.call(cbind, energy_supply_yearound)
  energy_efficiency_score <- do.call(cbind, energy_efficiency_score)
  energy_efficiency_score[nrow(energy_efficiency_score), ncol(energy_efficiency_score)] <- 0
  
  selected_distribution <- vector()
  # Run the model until 99% of the energy supply is consumed
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
  
  SEDS_output_future[[k]] <- df %>%
    mutate(richness_simu = c(summer_richness_simu_S[1:(nrow(df)/2)], winter_richness_simu_S[1:(nrow(df)/2)]))
  
  print(k)
}
SEDS_output_future <- do.call(rbind, SEDS_output_future)

# Putting all the data into one data frame for plotting
richness_all_df_future <- list()
for(k in 1:length(slopes)){
  richness_all_df_future[[k]] <- SEDS_output_future %>% filter(slope == slopes[k] & season == "summer") %>% dplyr::select(richness_simu) %>% rename(richness_simu_summer = richness_simu) %>%
    mutate(richness_simu_winter = SEDS_output_future$richness_simu[SEDS_output_future$slope == slopes[k] & SEDS_output_future$season == "winter"]) %>%
    mutate(elevation = SEDS_output_future$elevation_bin[SEDS_output_future$slope == slopes[k] & SEDS_output_future$season == "winter"]) %>%
    mutate(slope = slopes[k]) %>%
    mutate(latitude = as.numeric(mountain_slopes_polys[k,] %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y))) %>%
    mutate(richness_difference_simu = richness_simu_summer - richness_simu_winter) %>% 
    left_join(richness_patterns_observed %>% filter(slope == slopes[k]))
}
richness_all_df_future <- do.call(rbind, richness_all_df_future)

richness_all_df_future$slope <- gsub("_"," ",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Coastal California","Coastal Range",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Eastern BC","Columbia Mountains",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Peru SouthWest","Peru Pacific",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Peru East","Peru Amazon",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Ecuador West","Ecuador Pacific",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Colombia West","Colombia Pacific",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Chiapas East","Chiapas",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("S Nevada","Sierra Nevada Spain",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Costa Rica West","Costa Rica Pacific",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Costa Rica East","Costa Rica Caribbean",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Cord Merida","Cord. de Mérida",richness_all_df_future$slope)
richness_all_df_future$slope <- gsub("Nuevo Leon","Nuevo León",richness_all_df_future$slope)

pdf(file=paste0("results/figures/Fig_S3_future.pdf"), width = 10, height = 10)
ggplot() +
  geom_hline(yintercept = 0, col="grey50") +
  geom_point(data=richness_all_df, aes(x=elevation, y=richness_difference_simu), col="grey35", size=0.8) +
  geom_smooth(data=richness_all_df, aes(x=elevation, y=richness_difference_simu), col="grey20", linewidth = 1, se=T, method="loess", alpha=0.1) +
  geom_point(data=richness_all_df_future, aes(x=elevation, y=richness_difference_simu), col="red", size=0.8) +
  geom_smooth(data=richness_all_df_future, aes(x=elevation, y=richness_difference_simu), col="red", linewidth = 1, se=T, method="loess", alpha=0.1) +
  xlab("Elevation") + ylab("Seasonal difference in richness") + 
  facet_wrap(~reorder(slope, -latitude), scales="free_y") + theme_classic()
dev.off()


# Estimate features of the future predicted patterns

slopes2 <- unique(richness_all_df$slope)
ss_diff_simu_future <- list()
for(k in 1:length(slopes2)){
  richness_all_df_future_sub <- richness_all_df_future %>% filter(slope == slopes2[k])
  richness_all_df_future_sub <- richness_all_df_future_sub[which(is.na(richness_all_df_future_sub$richness_difference)==F),]
  ss_diff_simu_future[[k]] <- smooth.spline(x=richness_all_df_future_sub$elevation, y=richness_all_df_future_sub$richness_difference_simu, spar=0.5)
}
names(ss_diff_simu_future) <- slopes

# Simulated patterns for the future
transition_elev_simu_future <- unlist(list(
  Swiss_Alps = NA,
  Australian_Alps = NA,
  Blue_Ridge = NA,
  Cascade_East = NA,
  Cascade_West = NA,
  Coastal_California = seq(min(ss_diff_simu_future$Coastal_California$x), 2000, 1)[which.min(abs(predict(ss_diff_simu_future$Coastal_California, seq(min(ss_diff_simu_future$Coastal_California$x), 2000, 1))$y))],
  Eastern_BC = NA,
  Cord_Cantabrica = NA,
  Cord_Merida = NA,
  Peru_SouthWest = NA,
  Peru_East = seq(min(ss_diff_simu_future$Peru_East$x), 4500, 1)[which.min(abs(predict(ss_diff_simu_future$Peru_East, seq(min(ss_diff_simu_future$Peru_East$x), 4500, 1))$y))],
  Argentina_North = NA, 
  Ecuador_West = NA,
  Colombia_West = NA,
  Drakensberg = NA, 
  Pyrenees_Atlantic = NA,
  Pyrenees_Catalonia = seq(min(ss_diff_simu_future$Pyrenees_Catalonia$x), max(ss_diff_simu_future$Pyrenees_Catalonia$x), 1)[which.min(abs(predict(ss_diff_simu_future$Pyrenees_Catalonia, seq(min(ss_diff_simu_future$Pyrenees_Catalonia$x), max(ss_diff_simu_future$Pyrenees_Catalonia$x), 1))$y))],
  Serra_Mantiqueira = NA, 
  Chiapas_East = seq(min(ss_diff_simu_future$Chiapas_East$x), 3000, 1)[which.min(abs(predict(ss_diff_simu_future$Chiapas_East, seq(min(ss_diff_simu_future$Chiapas_East$x), 3000, 1))$y))],
  Veracruz = seq(min(ss_diff_simu_future$Veracruz$x), 3000, 1)[which.min(abs(predict(ss_diff_simu_future$Veracruz, seq(min(ss_diff_simu_future$Veracruz$x), 3000, 1))$y))], 
  Jalisco = NA,
  Sierra_Occidental = NA, 
  Nuevo_Leon = seq(min(ss_diff_simu_future$Nuevo_Leon$x), 2000, 1)[which.min(abs(predict(ss_diff_simu_future$Nuevo_Leon, seq(min(ss_diff_simu_future$Nuevo_Leon$x), 2000, 1))$y))], 
  California_Sierra_Nevada = NA,
  S_Nevada = NA, 
  Southern_Ghats_East = NA, 
  Southern_Ghats_West = NA, 
  Chile_Central = seq(min(ss_diff_simu_future$Chile_Central$x), 2000, 1)[which.min(abs(predict(ss_diff_simu_future$Chile_Central, seq(min(ss_diff_simu_future$Chile_Central$x), 2000, 1))$y))], 
  Western_Himalayas = seq(min(ss_diff_simu_future$Western_Himalayas$x), 2000, 1)[which.min(abs(predict(ss_diff_simu_future$Western_Himalayas, seq(min(ss_diff_simu_future$Western_Himalayas$x), 2000, 1))$y))], 
  Hawaii = seq(min(ss_diff_simu_future$Hawaii$x), 2500, 1)[which.min(abs(predict(ss_diff_simu_future$Hawaii, seq(min(ss_diff_simu_future$Hawaii$x), 2500, 1))$y))],
  Taiwan_East = seq(min(ss_diff_simu_future$Taiwan_East$x), max(ss_diff_simu_future$Taiwan_East$x), 1)[which.min(abs(predict(ss_diff_simu_future$Taiwan_East, seq(min(ss_diff_simu_future$Taiwan_East$x), max(ss_diff_simu_future$Taiwan_East$x), 1))$y))], 
  Taiwan_West = NA, 
  Costa_Rica_East = seq(min(ss_diff_simu_future$Costa_Rica_East$x), 2500, 1)[which.min(abs(predict(ss_diff_simu_future$Costa_Rica_East, seq(min(ss_diff_simu_future$Costa_Rica_East$x), 2500, 1))$y))],
  Costa_Rica_West = NA)) 

summer_peak_simu_future <- unlist(list(
  Swiss_Alps = seq(min(ss_diff_simu_future$Swiss_Alps$x), max(ss_diff_simu_future$Swiss_Alps$x), 1)[which.max(predict(ss_diff_simu_future$Swiss_Alps, seq(min(ss_diff_simu_future$Swiss_Alps$x), max(ss_diff_simu_future$Swiss_Alps$x), 1))$y)], 
  Australian_Alps = NA,
  Blue_Ridge = seq(300, max(ss_diff_simu_future$Blue_Ridge$x), 1)[which.max(predict(ss_diff_simu_future$Blue_Ridge, seq(300, max(ss_diff_simu_future$Blue_Ridge$x), 1))$y)], 
  Cascade_East = seq(min(ss_diff_simu_future$Cascade_East$x), max(ss_diff_simu_future$Cascade_East$x), 1)[which.max(predict(ss_diff_simu_future$Cascade_East, seq(min(ss_diff_simu_future$Cascade_East$x), max(ss_diff_simu_future$Cascade_East$x), 1))$y)], 
  Cascade_West = seq(min(ss_diff_simu_future$Cascade_West$x), max(ss_diff_simu_future$Cascade_West$x), 1)[which.max(predict(ss_diff_simu_future$Cascade_West, seq(min(ss_diff_simu_future$Cascade_West$x), max(ss_diff_simu_future$Cascade_West$x), 1))$y)], 
  Coastal_California = seq(min(ss_diff_simu_future$Coastal_California$x), max(ss_diff_simu_future$Coastal_California$x), 1)[which.max(predict(ss_diff_simu_future$Coastal_California, seq(min(ss_diff_simu_future$Coastal_California$x), max(ss_diff_simu_future$Coastal_California$x), 1))$y)],
  Eastern_BC = seq(min(ss_diff_simu_future$Eastern_BC$x), max(ss_diff_simu_future$Eastern_BC$x), 1)[which.max(predict(ss_diff_simu_future$Eastern_BC, seq(min(ss_diff_simu_future$Eastern_BC$x), max(ss_diff_simu_future$Eastern_BC$x), 1))$y)], 
  Cord_Cantabrica = seq(min(ss_diff_simu_future$Cord_Cantabrica$x), max(ss_diff_simu_future$Cord_Cantabrica$x), 1)[which.max(predict(ss_diff_simu_future$Cord_Cantabrica, seq(min(ss_diff_simu_future$Cord_Cantabrica$x), max(ss_diff_simu_future$Cord_Cantabrica$x), 1))$y)],
  Cord_Merida = NA,
  Peru_SouthWest = NA, 
  Peru_East = seq(min(ss_diff_simu_future$Peru_East$x), max(ss_diff_simu_future$Peru_East$x), 1)[which.max(predict(ss_diff_simu_future$Peru_East, seq(min(ss_diff_simu_future$Peru_East$x), max(ss_diff_simu_future$Peru_East$x), 1))$y)],
  Argentina_North = NA,
  Ecuador_West = NA,
  Colombia_West = NA,
  Drakensberg = NA,
  Pyrenees_Atlantic = seq(min(ss_diff_simu_future$Pyrenees_Atlantic$x), max(ss_diff_simu_future$Pyrenees_Atlantic$x), 1)[which.max(predict(ss_diff_simu_future$Pyrenees_Atlantic, seq(min(ss_diff_simu_future$Pyrenees_Atlantic$x), max(ss_diff_simu_future$Pyrenees_Atlantic$x), 1))$y)],
  Pyrenees_Catalonia = seq(min(ss_diff_simu_future$Pyrenees_Catalonia$x), max(ss_diff_simu_future$Pyrenees_Catalonia$x), 1)[which.max(predict(ss_diff_simu_future$Pyrenees_Catalonia, seq(min(ss_diff_simu_future$Pyrenees_Catalonia$x), max(ss_diff_simu_future$Pyrenees_Catalonia$x), 1))$y)], 
  Serra_Mantiqueira = NA,
  Chiapas_East = seq(min(ss_diff_simu_future$Chiapas_East$x), max(ss_diff_simu_future$Chiapas_East$x), 1)[which.max(predict(ss_diff_simu_future$Chiapas_East, seq(min(ss_diff_simu_future$Chiapas_East$x), max(ss_diff_simu_future$Chiapas_East$x), 1))$y)],
  Veracruz = NA, 
  Jalisco = NA, 
  Sierra_Occidental = NA,
  Nuevo_Leon = seq(min(ss_diff_simu_future$Nuevo_Leon$x), max(ss_diff_simu_future$Nuevo_Leon$x), 1)[which.max(predict(ss_diff_simu_future$Nuevo_Leon, seq(min(ss_diff_simu_future$Nuevo_Leon$x), max(ss_diff_simu_future$Nuevo_Leon$x), 1))$y)], 
  California_Sierra_Nevada = seq(min(ss_diff_simu_future$California_Sierra_Nevada$x), max(ss_diff_simu_future$California_Sierra_Nevada$x), 1)[which.max(predict(ss_diff_simu_future$California_Sierra_Nevada, seq(min(ss_diff_simu_future$California_Sierra_Nevada$x), max(ss_diff_simu_future$California_Sierra_Nevada$x), 1))$y)], 
  S_Nevada = seq(min(ss_diff_simu_future$S_Nevada$x), max(ss_diff_simu_future$S_Nevada$x), 1)[which.max(predict(ss_diff_simu_future$S_Nevada, seq(min(ss_diff_simu_future$S_Nevada$x), max(ss_diff_simu_future$S_Nevada$x), 1))$y)],
  Southern_Ghats_East = NA, 
  Southern_Ghats_West = NA,
  Chile_Central = NA,
  Western_Himalayas = seq(min(ss_diff_simu_future$Western_Himalayas$x), max(ss_diff_simu_future$Western_Himalayas$x), 1)[which.max(predict(ss_diff_simu_future$Western_Himalayas, seq(min(ss_diff_simu_future$Western_Himalayas$x), max(ss_diff_simu_future$Western_Himalayas$x), 1))$y)], 
  Hawaii = seq(min(ss_diff_simu_future$Hawaii$x), max(ss_diff_simu_future$Hawaii$x), 1)[which.max(predict(ss_diff_simu_future$Hawaii, seq(min(ss_diff_simu_future$Hawaii$x), max(ss_diff_simu_future$Hawaii$x), 1))$y)],
  Taiwan_East = seq(min(ss_diff_simu_future$Taiwan_East$x), max(ss_diff_simu_future$Taiwan_East$x), 1)[which.max(predict(ss_diff_simu_future$Taiwan_East, seq(min(ss_diff_simu_future$Taiwan_East$x), max(ss_diff_simu_future$Taiwan_East$x), 1))$y)], 
  Taiwan_West = NA, 
  Costa_Rica_East = seq(min(ss_diff_simu_future$Costa_Rica_East$x), max(ss_diff_simu_future$Costa_Rica_East$x), 1)[which.max(predict(ss_diff_simu_future$Costa_Rica_East, seq(min(ss_diff_simu_future$Costa_Rica_East$x), max(ss_diff_simu_future$Costa_Rica_East$x), 1))$y)],
  Costa_Rica_West = NA)) 

winter_peak_simu_future <- unlist(list(
  Swiss_Alps = NA,
  Australian_Alps = seq(min(ss_diff_simu_future$Australian_Alps$x), max(ss_diff_simu_future$Australian_Alps$x), 1)[which.min(predict(ss_diff_simu_future$Australian_Alps, seq(min(ss_diff_simu_future$Australian_Alps$x), max(ss_diff_simu_future$Australian_Alps$x), 1))$y)],
  Blue_Ridge = NA,
  Cascade_East = NA,
  Cascade_West = NA,
  Coastal_California = NA,
  Eastern_BC = NA,
  Cord_Cantabrica = NA,
  Cord_Merida = seq(min(ss_diff_simu_future$Cord_Merida$x), max(ss_diff_simu_future$Cord_Merida$x), 1)[which.min(predict(ss_diff_simu_future$Cord_Merida, seq(min(ss_diff_simu_future$Cord_Merida$x), max(ss_diff_simu_future$Cord_Merida$x), 1))$y)],
  Peru_SouthWest = seq(min(ss_diff_simu_future$Peru_SouthWest$x), max(ss_diff_simu_future$Peru_SouthWest$x), 1)[which.min(predict(ss_diff_simu_future$Peru_SouthWest, seq(min(ss_diff_simu_future$Peru_SouthWest$x), max(ss_diff_simu_future$Peru_SouthWest$x), 1))$y)], 
  Peru_East = NA,
  Argentina_North = seq(min(ss_diff_simu_future$Argentina_North$x), max(ss_diff_simu_future$Argentina_North$x), 1)[which.min(predict(ss_diff_simu_future$Argentina_North, seq(min(ss_diff_simu_future$Argentina_North$x), max(ss_diff_simu_future$Argentina_North$x), 1))$y)],
  Ecuador_West = NA,
  Colombia_West = seq(min(ss_diff_simu_future$Colombia_West$x), max(ss_diff_simu_future$Colombia_West$x), 1)[which.min(predict(ss_diff_simu_future$Colombia_West, seq(min(ss_diff_simu_future$Colombia_West$x), max(ss_diff_simu_future$Colombia_West$x), 1))$y)],
  Drakensberg = seq(min(ss_diff_simu_future$Drakensberg$x), max(ss_diff_simu_future$Drakensberg$x), 1)[which.min(predict(ss_diff_simu_future$Drakensberg, seq(min(ss_diff_simu_future$Drakensberg$x), max(ss_diff_simu_future$Drakensberg$x), 1))$y)], 
  Pyrenees_Atlantic = NA,
  Pyrenees_Catalonia = NA,
  Serra_Mantiqueira = NA, 
  Chiapas_East = seq(min(ss_diff_simu_future$Chiapas_East$x), max(ss_diff_simu_future$Chiapas_East$x), 1)[which.min(predict(ss_diff_simu_future$Chiapas_East, seq(min(ss_diff_simu_future$Chiapas_East$x), max(ss_diff_simu_future$Chiapas_East$x), 1))$y)], 
  Veracruz = NA, 
  Jalisco = seq(min(ss_diff_simu_future$Jalisco$x), max(ss_diff_simu_future$Jalisco$x), 1)[which.min(predict(ss_diff_simu_future$Jalisco, seq(min(ss_diff_simu_future$Jalisco$x), max(ss_diff_simu_future$Jalisco$x), 1))$y)], 
  Sierra_Occidental = seq(min(ss_diff_simu_future$Sierra_Occidental$x), max(ss_diff_simu_future$Sierra_Occidental$x), 1)[which.min(predict(ss_diff_simu_future$Sierra_Occidental, seq(min(ss_diff_simu_future$Sierra_Occidental$x), max(ss_diff_simu_future$Sierra_Occidental$x), 1))$y)],
  Nuevo_Leon = NA, 
  California_Sierra_Nevada = NA,
  S_Nevada = NA,
  Southern_Ghats_East = NA, 
  Southern_Ghats_West = seq(min(ss_diff_simu_future$Southern_Ghats_West$x), max(ss_diff_simu_future$Southern_Ghats_West$x), 1)[which.min(predict(ss_diff_simu_future$Southern_Ghats_West, seq(min(ss_diff_simu_future$Southern_Ghats_West$x), max(ss_diff_simu_future$Southern_Ghats_West$x), 1))$y)], 
  Chile_Central = seq(min(ss_diff_simu_future$Chile_Central$x), max(ss_diff_simu_future$Chile_Central$x), 1)[which.min(predict(ss_diff_simu_future$Chile_Central, seq(min(ss_diff_simu_future$Chile_Central$x), max(ss_diff_simu_future$Chile_Central$x), 1))$y)], 
  Western_Himalayas = NA, 
  Hawaii = NA, 
  Taiwan_East = seq(min(ss_diff_simu_future$Taiwan_East$x), max(ss_diff_simu_future$Taiwan_East$x), 1)[which.min(predict(ss_diff_simu_future$Taiwan_East, seq(min(ss_diff_simu_future$Taiwan_East$x), max(ss_diff_simu_future$Taiwan_East$x), 1))$y)], 
  Taiwan_West = seq(min(ss_diff_simu_future$Taiwan_West$x), max(ss_diff_simu_future$Taiwan_West$x), 1)[which.min(predict(ss_diff_simu_future$Taiwan_West, seq(min(ss_diff_simu_future$Taiwan_West$x), max(ss_diff_simu_future$Taiwan_West$x), 1))$y)],
  Costa_Rica_East = NA, 
  Costa_Rica_West = NA)) 


k = 34
ggplot() + 
  geom_smooth(data=richness_all_df %>% filter(slope == slopes2[k]), aes(x=elevation, y=richness_difference_simu), col="black", se=T, method="loess") +
  geom_point(data=richness_all_df %>% filter(slope == slopes2[k]), aes(x=elevation, y=richness_difference_simu)) + theme_bw() +
  geom_smooth(data=richness_all_df_future %>% filter(slope == slopes2[k]), aes(x=elevation, y=richness_difference_simu), col="grey45", se=T, method="loess") +
  geom_point(data=richness_all_df_future %>% filter(slope == slopes2[k]), aes(x=elevation, y=richness_difference_simu), col="grey45") + theme_bw() +
  geom_hline(yintercept = 0, col="grey50") +
  xlab("") + ylab("") + xlim(c(0,5650))


# Put all the pattern features in the same data frame
data_diff_mnts <- data.frame(
  mountain_slope = slopes,
  latitude = mountain_slopes_polys %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y),
  lower_limit_obs = unlist(lapply(ss_diff_obs, function(x) min(x$x))),
  winter_peak_obs = winter_peak_obs,
  transition_obs = transition_elev_obs,
  summer_peak_obs = summer_peak_obs,
  upper_limit_obs = unlist(lapply(ss_diff_obs, function(x) max(x$x))),
  lower_limit_simu = unlist(lapply(ss_diff_simu, function(x) min(x$x))),
  winter_peak_simu = winter_peak_simu,
  transition_simu = transition_elev_simu,
  summer_peak_simu = summer_peak_simu,
  upper_limit_simu = unlist(lapply(ss_diff_simu, function(x) max(x$x))),
  lower_limit_simu_future = unlist(lapply(ss_diff_simu_future, function(x) min(x$x))),
  winter_peak_simu_future = winter_peak_simu_future,
  transition_simu_future = transition_elev_simu_future,
  summer_peak_simu_future = summer_peak_simu_future,
  upper_limit_simu_future = unlist(lapply(ss_diff_simu_future, function(x) max(x$x))))
rownames(data_diff_mnts) <- NULL

# Estimating the degree of upslope shift predicted by the model in the future compared to the present
winter_peak_shift <- data_diff_mnts$winter_peak_simu_future - data_diff_mnts$winter_peak_simu
transition_shift <- data_diff_mnts$transition_simu_future - data_diff_mnts$transition_simu
summer_peak_shift <- data_diff_mnts$summer_peak_simu_future - data_diff_mnts$summer_peak_simu
mean(winter_peak_shift, na.rm=T)
mean(transition_shift, na.rm=T)
mean(summer_peak_shift, na.rm=T)

toplot <- data.frame(feature = factor(c(rep("Boreal winter peak", length(winter_peak_shift)), rep("Transition", length(winter_peak_shift)), rep("Boreal summer peak", length(winter_peak_shift))), levels=c("Boreal winter peak", "Transition", "Boreal summer peak"), ordered = T),
                     val = c(winter_peak_shift, transition_shift, summer_peak_shift))

g_future <- ggplot(data=toplot, aes(y=val, fill=feature)) +
  geom_boxplot(alpha=0.7) +
  xlab("") + ylab("Elevation change") + theme_classic() +
  scale_fill_manual(values=c("blue","black","red"), name="") + theme(legend.position = c(0.8, 0.12), legend.box.background = element_rect(colour = "black"), legend.title = element_blank(), axis.text.x=element_blank()) + 
  labs(tag="A") + ylim(c(-300, 600))


t.test(winter_peak_shift)
t.test(transition_shift)
t.test(summer_peak_shift)
t.test(c(winter_peak_shift, transition_shift, summer_peak_shift))


## PLOT FIGURE 5
# g_E_Him is plotted in the script 10_eastern_himalayas.R
pdf(file = "results/figures/Fig_5.pdf", bg = "white", width = 11, height = 5)
ggarrange(g_future, g_E_Him, nrow=1, ncol=2, widths = c(0.75, 1), align="h")
dev.off()


