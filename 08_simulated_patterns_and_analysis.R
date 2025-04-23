## This script contains the code for simulating the patterns based on energy efficiency ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(ggpubr)
library(lme4)
library(lmerTest)
library(MuMIn)

# Load the data prepared in the previous scripts
load("mnt_slopes_data_new.RData")
mountain_ranges <- vect("mountain_ranges_2.shp")
species_elevational_ranges <- read.csv("species_elevational_ranges.csv")


# Global latitudinal patterns of bird migration

data_global_rich <- read.csv("PresAbs_global/isea3h7_niche_apet.csv")
data_global_rich_sf <- data_global_rich %>% st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs='+proj=longlat +datum=WGS84')

# Americas
ss <- smooth.spline(data_global_rich$LATITUDE[data_global_rich$LONGITUDE <= -30], data_global_rich$DIFFERENCE[data_global_rich$LONGITUDE <= -30])
pp <- predict(ss, seq(min(data_global_rich$LATITUDE[data_global_rich$LONGITUDE <= -30]), max(data_global_rich$LATITUDE[data_global_rich$LONGITUDE <= -30]), 0.1))
summer_peak_americas <- pp$x[which.max(pp$y)]
summer_peak_americas_magnitude <- pp$y[which.max(pp$y)]
winter_peak_americas <- pp$x[which.min(pp$y)]
winter_peak_americas_magnitude <- pp$y[which.min(pp$y)]
transition_americas <- mean(pp$x[c(min(which(pp$y > 0)), min(which(pp$y > 0)) - 1)])
# Central-East Asia
ss <- smooth.spline(data_global_rich$LATITUDE[data_global_rich$LONGITUDE > 50], data_global_rich$DIFFERENCE[data_global_rich$LONGITUDE > 50])
pp <- predict(ss, seq(-40, max(data_global_rich$LATITUDE[data_global_rich$LONGITUDE > 50]), 0.1))
summer_peak_east_asia <- pp$x[which.max(pp$y)]
summer_peak_east_asia_magnitude <- pp$y[which.max(pp$y)]
winter_peak_east_asia <- pp$x[which.min(pp$y)]
winter_peak_east_asia_magnitude <- pp$y[which.min(pp$y)]
transition_east_asia <- mean(pp$x[c(min(which(pp$y > 0)), min(which(pp$y > 0)) - 1)])
# Europe Africa
ss <- smooth.spline(data_global_rich$LATITUDE[data_global_rich$LONGITUDE > -30 & data_global_rich$LONGITUDE <= 50], data_global_rich$DIFFERENCE[data_global_rich$LONGITUDE > -30 & data_global_rich$LONGITUDE <= 50])
pp <- predict(ss, seq(min(data_global_rich$LATITUDE[data_global_rich$LONGITUDE > -30 & data_global_rich$LONGITUDE <= 50]), max(data_global_rich$LATITUDE[data_global_rich$LONGITUDE > -30 & data_global_rich$LONGITUDE <= 50]), 0.1))
summer_peak_euro_africa <- pp$x[which.max(pp$y)]
summer_peak_euro_africa_magnitude <- pp$y[which.max(pp$y)]
winter_peak_euro_africa <- pp$x[which.min(pp$y)]
winter_peak_euro_africa_magnitude <- pp$y[which.min(pp$y)]
transition_euro_africa <- mean(pp$x[c(min(which(pp$y > 0 & pp$x > 0)), min(which(pp$y > 0 & pp$x > 0)) - 1)])

# Which global flyway does each mountain slope belongs to?
flyway <- c("euro_africa","east_asia","americas","americas","americas","americas","americas","euro_africa","americas","americas","americas","americas","americas","americas","euro_africa","euro_africa","euro_africa","americas","americas","americas","americas","americas","americas", "americas","euro_africa","east_asia","east_asia","americas","east_asia","americas","east_asia","east_asia","americas","americas")

# Latitudinal summer peaks corresponding to each mountain slope
summer_peak_global <- c(summer_peak_euro_africa,summer_peak_east_asia,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_euro_africa,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_euro_africa,summer_peak_euro_africa,summer_peak_euro_africa,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas,summer_peak_americas, summer_peak_americas,summer_peak_euro_africa,summer_peak_east_asia,summer_peak_east_asia,summer_peak_americas,summer_peak_east_asia,summer_peak_americas,summer_peak_east_asia,summer_peak_east_asia,summer_peak_americas,summer_peak_americas)

# Latitudinal winter peaks corresponding to each mountain slope
winter_peak_global <- c(winter_peak_euro_africa,winter_peak_east_asia,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_euro_africa,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_euro_africa,winter_peak_euro_africa,winter_peak_euro_africa,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas,winter_peak_americas, winter_peak_americas,winter_peak_euro_africa,winter_peak_east_asia,winter_peak_east_asia,winter_peak_americas,winter_peak_east_asia,winter_peak_americas,winter_peak_east_asia,winter_peak_east_asia,winter_peak_americas,winter_peak_americas)

# Latitudinal transition corresponding to each mountain slope
transition_global <- c(transition_euro_africa,transition_east_asia,transition_americas,transition_americas,transition_americas,transition_americas,transition_americas,transition_euro_africa,transition_americas,transition_americas,transition_americas,transition_americas,transition_americas,transition_americas,transition_euro_africa,transition_euro_africa,transition_euro_africa,transition_americas,transition_americas,transition_americas,transition_americas,transition_americas,transition_americas, transition_americas,transition_euro_africa,transition_east_asia,transition_east_asia,transition_americas,transition_east_asia,transition_americas,transition_east_asia,transition_east_asia,transition_americas,transition_americas)

# Plot global map of seasonal difference in richness with the location of the mountain slopes
newmap2 <- ne_coastline()
theme_set(theme_void())
g_slopes2 <- ggplot() +
  geom_sf(data = data_global_rich_sf, aes(col=DIFFERENCE), size=1.3) +
  scale_colour_gradient2("Richness difference", trans = "reverse") +
  geom_sf(data = mountain_slopes_polys, col="black", fill="black") +
  ylim(c(-55,80)) + theme(legend.position='bottom') +
  geom_sf(data=newmap2, col="grey30", fill="grey30")
pdf(file = paste0("results/figures/mountain_slopes_map_diff.pdf"), width = 12, height = 7)
g_slopes2
dev.off()



##  Simulating the elevational seasonal distribution of birds  ##

# Energy supply per elevational band
env_data <- read_csv("env_data.csv")[,-1]
elevation_bins <- seq(0, 7000, 200)
slopes <- mountain_slopes_polys$name
env_data_elevation_bins <- list()
for(k in 1:length(slopes)){
  env_data_sub_summer <- env_data %>% dplyr::filter(slope == slopes[k] & season == "summer")
  env_data_sub_winter <- env_data %>% dplyr::filter(slope == slopes[k] & season == "winter")
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
  env_data_elevation_bins[[k]] <- rbind(env_data_elevation_bins_summer, env_data_elevation_bins_winter)
}
env_data_elevation_bins <- do.call(rbind, env_data_elevation_bins)
env_data_elevation_bins$ndvi[is.na(env_data_elevation_bins$ndvi) == T] <- 0
env_data_elevation_bins <- env_data_elevation_bins[is.na(env_data_elevation_bins$temperature) == F,]

ggplot() +
  geom_smooth(data=env_data_elevation_bins %>% filter(slope == "Swiss_Alps"), aes(x=elevation_bin, y=ndvi, col=season))
ggplot() +
  geom_smooth(data=env_data_elevation_bins %>% filter(slope == "California_Sierra_Nevada"), aes(x=elevation_bin, y=ndvi, col=season))


# Regional proportion of migrants for each mountain slope

prop_long_dist_migr <- prop_long_dist_migr_summer <- prop_long_dist_migr_winter <- vector()
for(k in 1:length(slopes)){
  species_elevational_ranges_sub <- species_elevational_ranges %>% filter(slope == slopes[k])
  sum(table(species_elevational_ranges_sub$species) - 1) # total number of species along the slope
  aa <- table(species_elevational_ranges_sub$species, species_elevational_ranges_sub$season)
  aa_residents = sum(apply(aa, 1, sum) - 1)
  aa_summer_migrants = length(which(aa[,1] == 1 & aa[,2] == 0))
  aa_winter_migrants = length(which(aa[,1] == 0 & aa[,2] == 1))
  aa_diff <- aa_summer_migrants - aa_winter_migrants
  if(aa_diff > 0){
    prop_long_dist_migr[k] = aa_diff / (aa_residents + aa_winter_migrants)
  }else{
    prop_long_dist_migr[k] = aa_diff / (aa_residents + aa_summer_migrants)
  }
  prop_long_dist_migr_summer[k] = aa_summer_migrants / (aa_residents + aa_summer_migrants)
  prop_long_dist_migr_winter[k] = aa_winter_migrants / (aa_residents + aa_winter_migrants)
}




## SEDS model

# Reproduction cost
lat_mnt <- mountain_slopes_polys %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y)
reproduction_proba_summer <- rep(0, length(slopes))
reproduction_proba_summer[which(lat_mnt > 15)] <- 1
reproduction_proba_summer[which(lat_mnt > 0 & lat_mnt < 15)] <- 0.75
reproduction_proba_summer[which(lat_mnt > -15 & lat_mnt < 0)] <- 0.25
repro_cost = 0.41 #0.19 #0.41 # 0.1 # 1
repro_cost_summer = repro_cost * reproduction_proba_summer 
repro_cost_winter = repro_cost * (1-reproduction_proba_summer)

# Long distance (latitudinal) migration cost
migra_dist <- rep(0, length(slopes))
migra_dist[prop_long_dist_migr > 0 & lat_mnt$latitude > 0] <- (lat_mnt$latitude[prop_long_dist_migr > 0 & lat_mnt$latitude > 0] - winter_peak_global[prop_long_dist_migr > 0 & lat_mnt$latitude > 0]) * 111
migra_dist[prop_long_dist_migr < 0 & lat_mnt$latitude > 0] <- abs(lat_mnt$latitude[prop_long_dist_migr < 0 & lat_mnt$latitude > 0] - summer_peak_global[prop_long_dist_migr < 0 & lat_mnt$latitude > 0]) * 111
migra_dist[lat_mnt$latitude < 0] <- 3000
migra_cost <- migra_dist * 6.45e-5

# Thermoregulation cost for long distance migrants
beta = 23.6 # 33.6 # 13.6
temp_peak_summer_americas <- 12.72042
temp_peak_winter_americas <- 18.32428
temp_peak_summer_euro_africa <- 11.65147
temp_peak_winter_euro_africa <- 26.66647
temp_peak_summer_east_asia <- 11.84833
temp_peak_winter_east_asia <- 21.29227
thermo_cost <- rep(0, length(slopes))
thermo_cost[prop_long_dist_migr > 0 & lat_mnt$latitude > 0 & flyway=="americas"] <- ifelse(temp_peak_winter_americas < (40 - beta), (40 - temp_peak_winter_americas - beta) / beta, 0)
thermo_cost[prop_long_dist_migr > 0 & lat_mnt$latitude > 0 & flyway=="euro_africa"] <- ifelse(temp_peak_winter_euro_africa < (40 - beta), (40 - temp_peak_winter_euro_africa - beta) / beta, 0)
thermo_cost[prop_long_dist_migr > 0 & lat_mnt$latitude > 0 & flyway=="east_asia"] <- ifelse(temp_peak_winter_east_asia < (40 - beta), (40 - temp_peak_winter_east_asia - beta) / beta, 0)
thermo_cost[prop_long_dist_migr < 0 & lat_mnt$latitude > 0 & flyway=="americas"] <- ifelse(temp_peak_summer_americas < (40 - beta), (40 - temp_peak_summer_americas - beta) / beta, 0)
thermo_cost[prop_long_dist_migr < 0 & lat_mnt$latitude > 0 & flyway=="euro_africa"] <- ifelse(temp_peak_summer_euro_africa < (40 - beta), (40 - temp_peak_summer_euro_africa - beta) / beta, 0)
thermo_cost[prop_long_dist_migr < 0 & lat_mnt$latitude > 0 & flyway=="east_asia"] <- ifelse(temp_peak_summer_east_asia < (40 - beta), (40 - temp_peak_summer_east_asia - beta) / beta, 0)
thermo_cost[lat_mnt$latitude < 0] <- ifelse(20 < (40 - beta), (40 - 20 - beta) / beta, 0)

SEDS_output <- simulated_distributions <- list()
for(k in 1:length(slopes)){
  # Prep environmental data
  df <- env_data_elevation_bins %>% filter(slope == slopes[k]) %>%
    mutate(energy_supply = area * ndvi) #%>% filter(is.na(energy_supply) == F)
  energy_supply_summer <- df$energy_supply[df$season == "summer"]
  energy_supply_winter <- df$energy_supply[df$season == "winter"]
  total_energy_supply_summer <- sum(energy_supply_summer, na.rm=T)
  total_energy_supply_winter <- sum(energy_supply_winter, na.rm=T)
  
  #print(table(df$season)) # check that there are the same number of elevational bins at both seasons
  
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
  simulated_distributions[[k]] <- selected_distribution %>% mutate(slope = slopes[k])
  summer_richness_sim <- table(selected_distribution$summer)
  winter_richness_sim <- table(selected_distribution$winter)
  if(prop_long_dist_migr[k] > 0){
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
  
  SEDS_output[[k]] <- df %>%
    mutate(richness_simu = c(summer_richness_simu_S[1:(nrow(df)/2)], winter_richness_simu_S[1:(nrow(df)/2)]))
  
  print(k)
}
SEDS_output <- do.call(rbind, SEDS_output)
simulated_distributions <- do.call(rbind, simulated_distributions)
rownames(simulated_distributions) <- NULL


## Seasonal richness patterns ##

# Empirical patterns
richness_patterns_observed <- list()
for(k in 1:nrow(mountain_slopes_polys)){
  species_elevational_ranges_sub_summer <- species_elevational_ranges %>% filter(slope == slopes[k] & season == "summer")
  species_elevational_ranges_sub_winter <- species_elevational_ranges %>% filter(slope == slopes[k] & season == "winter") 
  richness_summer <- richness_winter <- vector()
  for(i in 1:(length(elevation_bins)-1)){
    richness_summer[i] <- length(which(species_elevational_ranges_sub_summer$range_min <= elevation_bins[i] & species_elevational_ranges_sub_summer$range_max >= elevation_bins[i+1]))
    richness_winter[i] <- length(which(species_elevational_ranges_sub_winter$range_min <= elevation_bins[i] & species_elevational_ranges_sub_winter$range_max >= elevation_bins[i+1]))
  }
  richness_patterns_observed[[k]] <- data.frame(slope = slopes[k],
                                                elevation = elevation_bins[-length(elevation_bins)] + 100,
                                                richness_summer = richness_summer,
                                                richness_winter = richness_winter, 
                                                richness_difference = richness_summer - richness_winter)
  bins_without_species <- which(elevation_bins[-length(elevation_bins)] < min(species_elevational_ranges_sub_summer$range_min, species_elevational_ranges_sub_winter$range_min) | elevation_bins[-length(elevation_bins)] >= max(species_elevational_ranges_sub_summer$range_max, species_elevational_ranges_sub_winter$range_max))
  richness_patterns_observed[[k]][bins_without_species,3:5] <- NA
}
richness_patterns_observed <- do.call(rbind, richness_patterns_observed) 


# Putting empirical and simulated patterns into one data frame
richness_all_df <- list()
for(k in 1:length(slopes)){
  richness_all_df[[k]] <- SEDS_output %>% filter(slope == slopes[k] & season == "summer") %>% dplyr::select(richness_simu) %>% rename(richness_simu_summer = richness_simu) %>%
    mutate(richness_simu_winter = SEDS_output$richness_simu[SEDS_output$slope == slopes[k] & SEDS_output$season == "winter"]) %>%
    mutate(elevation = SEDS_output$elevation_bin[SEDS_output$slope == slopes[k] & SEDS_output$season == "winter"]) %>%
    mutate(slope = slopes[k]) %>%
    mutate(latitude = as.numeric(mountain_slopes_polys[k,] %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y))) %>%
    mutate(richness_difference_simu = richness_simu_summer - richness_simu_winter) %>% 
    left_join(richness_patterns_observed %>% filter(slope == slopes[k]))
}
richness_all_df <- do.call(rbind, richness_all_df)
#write_csv(richness_all_df, file="results/richness_all_df.csv")


## Correlation between observed and simulated

correlations_simu_obs_pval <- correlations_simu_obs_cor <- vector()
for(k in 1:length(slopes)){
  correlations_simu_obs <- cor.test(richness_all_df$richness_difference[richness_all_df$slope == slopes[k]], richness_all_df$richness_difference_simu[richness_all_df$slope == slopes[k]], use="pairwise.complete.obs")
  correlations_simu_obs_pval[k] <- correlations_simu_obs$p.value
  correlations_simu_obs_cor[k] <- correlations_simu_obs$estimate
}
correlations_simu_obs <- data.frame(slope = slopes, cor = correlations_simu_obs_cor, pval = correlations_simu_obs_pval, pval_signif = ifelse(correlations_simu_obs_pval < 0.05, "yes", "no"))

g_dens_cor <- ggplot(correlations_simu_obs) +
  geom_density(aes(x=cor)) + theme_classic() + xlim(c(-1,1)) + ylab("Density") + xlab("Correlation coefficient between simulated and observed") +
  geom_point(aes(x=cor, y=0, col=pval_signif)) + theme(legend.position="none") + scale_colour_manual(values=c("yes"="black", "no"="grey50")) + labs(tag="B")

g_cor <- ggplot() +
  geom_point(data=richness_all_df, aes(x=richness_difference, y=richness_difference_simu, by=slope, col=latitude), size=0.8) + theme_bw() +
  geom_smooth(data=richness_all_df, aes(x=richness_difference, y=richness_difference_simu, by=slope, col=latitude), se=F, method="lm", alpha=0.8) +
  geom_smooth(data=richness_all_df, aes(x=richness_difference, y=richness_difference_simu), col="grey20", linewidth = 1, se=F, method="lm") +
  xlab("Observed seasonal difference in richness") + ylab("Simulated seasonal difference in richness") + 
  scale_color_terrain_c(limits=c(-60,60), name="Latitude") + theme(legend.position = "right") + labs(tag="A") + ylim(c(-150, 150))


## Model predictions for thermal tracking
seas_tracking <- prop_migr <- vector()
for(k in 1:length(slopes)){
  simulated_distributions_slope <- simulated_distributions %>% filter(slope == slopes[k])
  SEDS_output_slope_winter <- SEDS_output %>% filter(slope == slopes[k] & season == "winter")
  SEDS_output_slope_summer <- SEDS_output %>% filter(slope == slopes[k] & season == "summer")
  # birds not migrating latitutinally
  simulated_distributions_slope <- simulated_distributions_slope %>% filter(winter <= nrow(SEDS_output_slope_winter) & summer <= nrow(SEDS_output_slope_summer))
  # altitudinal migrant populations
  simulated_distributions_slope_migrants <- simulated_distributions_slope %>% filter(abs(winter - summer) > 1)
  simulated_distributions_slope_sedentary <- simulated_distributions_slope %>% filter(abs(winter - summer) <= 1)
  # seasonal thermal tracking: seasonal thermal difference for altitudinal migrants compared to sedentary
  seas_tracking_migr <- abs(SEDS_output_slope_summer$temperature[simulated_distributions_slope_migrants$summer] - SEDS_output_slope_winter$temperature[simulated_distributions_slope_migrants$winter])
  seas_tracking_sed <- mean(abs(SEDS_output_slope_summer$temperature[simulated_distributions_slope_sedentary$summer] - SEDS_output_slope_winter$temperature[simulated_distributions_slope_sedentary$winter]))
  seas_tracking <- c(seas_tracking, seas_tracking_sed - seas_tracking_migr)
  prop_migr <- c(prop_migr, nrow(simulated_distributions_slope_migrants) / nrow(simulated_distributions_slope))
}    

g_hist_tracking <- ggplot() + 
  geom_histogram(aes(x=seas_tracking)) +
  theme_classic() + xlab("Seasonal thermal tracking") + ylab("Number of simulated avian populations") + labs(tag="C") +
  geom_vline(xintercept = 0,linetype="dotted")


pdf(file=paste0("results/figures/Fig_3_new3.pdf"), width = 15, height = 5)
ggarrange(g_cor, g_dens_cor, g_hist_tracking, widths = c(1, 0.75, 0.75), nrow=1)
dev.off()

# Median correlation between empirical and simulated patterns across all mountain slopes 
median(correlations_simu_obs$cor)
# How many mountain slopes have a positive correlation
length(which(correlations_simu_obs$cor > 0))
# How many mountain slopes have a correlation smaller than 0.15
length(which(correlations_simu_obs$cor < 0.15))

# Rename some mountain slopes
richness_all_df$slope <- gsub("_"," ",richness_all_df$slope)
richness_all_df$slope <- gsub("Coastal California","Coastal Range",richness_all_df$slope)
richness_all_df$slope <- gsub("Eastern BC","Rocky Mountains N",richness_all_df$slope)
richness_all_df$slope <- gsub("Peru SouthWest","Peru Pacific",richness_all_df$slope)
richness_all_df$slope <- gsub("Peru East","Peru Amazon",richness_all_df$slope)
richness_all_df$slope <- gsub("Ecuador West","Ecuador Pacific",richness_all_df$slope)
richness_all_df$slope <- gsub("Colombia West","Colombia Pacific",richness_all_df$slope)
richness_all_df$slope <- gsub("Chiapas East","Chiapas",richness_all_df$slope)
richness_all_df$slope <- gsub("S Nevada","Sierra Nevada Spain",richness_all_df$slope)
richness_all_df$slope <- gsub("Costa Rica West","Costa Rica Pacific",richness_all_df$slope)
richness_all_df$slope <- gsub("Costa Rica East","Costa Rica Caribbean",richness_all_df$slope)
richness_all_df$slope <- gsub("Cord Merida","Cord. de Mérida",richness_all_df$slope)
richness_all_df$slope <- gsub("Nuevo Leon","Nuevo León",richness_all_df$slope)

pdf(file=paste0("results/figures/Fig_S3_new.pdf"), width = 10, height = 10)
ggplot() +
  geom_hline(yintercept = 0, col="grey50") +
  geom_point(data=richness_all_df, aes(x=elevation, y=richness_difference), col="grey35", size=0.8) +
  geom_smooth(data=richness_all_df, aes(x=elevation, y=richness_difference), col="grey20", linewidth = 1, se=T, method="loess", alpha=0.1) +
  geom_point(data=richness_all_df, aes(x=elevation, y=richness_difference_simu), col="red", size=0.8) +
  geom_smooth(data=richness_all_df, aes(x=elevation, y=richness_difference_simu), col="red", linewidth = 1, se=T, method="loess", alpha=0.1) +
  xlab("Elevation") + ylab("Seasonal difference in richness") + 
  facet_wrap(~reorder(slope, -latitude), scales="free_y") + theme_classic()
dev.off()

write_csv(richness_all_df, file="richness_all_df.csv")

# Random slope model to test the correlation between simulated and empirical richness patterns globally
mod <- lmer(richness_difference ~ 1 + richness_difference_simu + (1 + richness_difference_simu | slope), data=richness_all_df)
summary(mod)
coef(mod)$slope
anova(mod)
r.squaredGLMM(mod)

cor(richness_all_df$richness_difference, richness_all_df$richness_difference_simu, use="pairwise.complete.obs")




## NULL MODEL: random seasonal redistributions

SEDS_output_null1 <- list()
for(k in 1:length(slopes)){
  SEDS_output_sub <- SEDS_output %>% filter(slope==slopes[k])
  species_total <- round(max(sum(SEDS_output_sub$richness_simu[SEDS_output_sub$season == "winter"]), sum(SEDS_output_sub$richness_simu[SEDS_output_sub$season == "summer"])))
  df <- env_data_elevation_bins %>% filter(slope == slopes[k]) %>%
    mutate(energy_supply = area * ndvi) %>%
    filter(is.na(energy_supply) == F)
  selected_distribution <- vector()
  for(i in 1:species_total){
    sel <- sample(size=2, 1:(nrow(df)/2))
    selected_distribution <- rbind(selected_distribution, sel)
  }
  selected_distribution <- as.data.frame(selected_distribution)
  colnames(selected_distribution) <- c("summer", "winter")
  summer_richness_sim <- table(selected_distribution$summer)
  winter_richness_sim <- table(selected_distribution$winter)
  SEDS_output_null1[[k]] <- df %>%
    mutate(richness_simu = c(summer_richness_sim, winter_richness_sim))
}
SEDS_output_null1 <- do.call(rbind, SEDS_output_null1)

richness_all_df_null1 <- list()
for(k in 1:length(slopes)){
  richness_all_df_null1[[k]] <- SEDS_output_null1 %>% filter(slope == slopes[k] & season == "summer") %>% dplyr::select(richness_simu) %>% rename(richness_simu_summer = richness_simu) %>%
    mutate(richness_simu_winter = SEDS_output_null1$richness_simu[SEDS_output_null1$slope == slopes[k] & SEDS_output_null1$season == "winter"]) %>%
    mutate(elevation = SEDS_output_null1$elevation_bin[SEDS_output_null1$slope == slopes[k] & SEDS_output_null1$season == "winter"]) %>%
    mutate(slope = slopes[k]) %>%
    mutate(latitude = as.numeric(mountain_slopes_polys[k,] %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y))) %>%
    mutate(richness_difference_simu = richness_simu_summer - richness_simu_winter) %>% 
    left_join(richness_patterns_observed %>% filter(slope == slopes[k]))
}
richness_all_df_null1 <- do.call(rbind, richness_all_df_null1)

correlations_simu_obs_pval <- correlations_simu_obs_cor <- vector()
for(k in 1:length(slopes)){
  correlations_simu_obs <- cor.test(richness_all_df_null1$richness_difference[richness_all_df_null1$slope == slopes[k]], richness_all_df_null1$richness_difference_simu[richness_all_df_null1$slope == slopes[k]], use="pairwise.complete.obs")
  correlations_simu_obs_pval[k] <- correlations_simu_obs$p.value
  correlations_simu_obs_cor[k] <- correlations_simu_obs$estimate
}
correlations_simu_obs_null1 <- data.frame(slope = slopes, cor = correlations_simu_obs_cor, pval = correlations_simu_obs_pval, pval_signif = ifelse(correlations_simu_obs_pval < 0.05, "yes", "no"))

g_dens_cor <- ggplot(correlations_simu_obs_null1) +
  geom_density(aes(x=cor)) + theme_classic() + xlim(c(-1,1)) + ylab("Density") + xlab("Correlation coefficient between simulated and observed") +
  geom_point(aes(x=cor, y=0, col=pval_signif)) + theme(legend.position="none") + scale_colour_manual(values=c("yes"="black", "no"="grey50")) + labs(tag="B")

g_cor <- ggplot() +
  geom_point(data=richness_all_df_null1, aes(x=richness_difference, y=richness_difference_simu, by=slope, col=latitude), size=0.8) + theme_bw() +
  geom_smooth(data=richness_all_df_null1, aes(x=richness_difference, y=richness_difference_simu, by=slope, col=latitude), se=F, method="lm", alpha=0.8) +
  geom_smooth(data=richness_all_df_null1, aes(x=richness_difference, y=richness_difference_simu), col="grey20", linewidth = 1, se=F, method="lm") +
  xlab("Observed seasonal difference in richness") + ylab("Simulated seasonal difference in richness") + 
  scale_color_terrain_c(limits=c(-60,60), name="Latitude") + theme(legend.position = "right") + labs(tag="A")

pdf(file=paste0("results/figures/Fig_S4.pdf"), width = 11, height = 5)
ggarrange(g_cor, g_dens_cor, widths = c(1, 0.7))
dev.off()

mod <- lmer(richness_difference ~ 1 + richness_difference_simu + (1 + richness_difference_simu | slope), data=richness_all_df_null1)
summary(mod)
coef(mod)$slope
anova(mod)
r.squaredGLMM(mod)

ks.test(correlations_simu_obs$cor, correlations_simu_obs_null1$cor, alternative="less")



######################################
##  Measuring peaks and transition  ##
######################################

slopes3 <- unique(richness_all_df$slope)
ss_diff_obs <- ss_diff_simu <- list()
for(k in 1:length(slopes3)){
  richness_all_df_sub <- richness_all_df %>% filter(slope == slopes3[k])
  richness_all_df_sub <- richness_all_df_sub[which(is.na(richness_all_df_sub$richness_difference)==F),]
  ss_diff_obs[[k]] <- smooth.spline(x=richness_all_df_sub$elevation, y=richness_all_df_sub$richness_difference, spar=0.5)
  ss_diff_simu[[k]] <- smooth.spline(x=richness_all_df_sub$elevation, y=richness_all_df_sub$richness_difference_simu, spar=0.5)
}
names(ss_diff_obs) <- names(ss_diff_simu) <- slopes


###  Features in observed patterns  ###

## Transition between net richness surpluses during the boreal winter and summer 
transition_elev_obs <- unlist(list(
  Swiss_Alps = NA,
  Australian_Alps = NA,
  Blue_Ridge = NA,
  Cascade_East = NA,
  Cascade_West = NA,
  Coastal_California = seq(min(ss_diff_obs$Coastal_California$x), max(ss_diff_obs$Coastal_California$x), 1)[which.min(abs(predict(ss_diff_obs$Coastal_California, seq(min(ss_diff_obs$Coastal_California$x), max(ss_diff_obs$Coastal_California$x), 1))$y))],
  Eastern_BC = NA,
  Cord_Cantabrica = NA,
  Cord_Merida = NA,
  Peru_SouthWest = seq(min(ss_diff_obs$Peru_SouthWest$x), max(ss_diff_obs$Peru_SouthWest$x), 1)[which.min(abs(predict(ss_diff_obs$Peru_SouthWest, seq(min(ss_diff_obs$Peru_SouthWest$x), max(ss_diff_obs$Peru_SouthWest$x), 1))$y))],
  Peru_East = NA,
  Argentina_North = NA,
  Ecuador_West = NA,
  Colombia_West = seq(min(ss_diff_obs$Colombia_West$x), max(ss_diff_obs$Colombia_West$x), 1)[which.min(abs(predict(ss_diff_obs$Colombia_West, seq(min(ss_diff_obs$Colombia_West$x), max(ss_diff_obs$Colombia_West$x), 1))$y))],
  Drakensberg = NA,
  Pyrenees_Atlantic = NA,
  Pyrenees_Catalonia = NA,
  Serra_Mantiqueira = NA,
  Chiapas_East = NA,
  Veracruz = NA, 
  Jalisco = seq(min(ss_diff_obs$Jalisco$x), 4000, 1)[which.min(abs(predict(ss_diff_obs$Jalisco, seq(min(ss_diff_obs$Jalisco$x), 4000, 1))$y))],
  Sierra_Occidental = NA,
  Nuevo_Leon = seq(min(ss_diff_obs$Nuevo_Leon$x), 2000, 1)[which.min(abs(predict(ss_diff_obs$Nuevo_Leon, seq(min(ss_diff_obs$Nuevo_Leon$x), 2000, 1))$y))],
  California_Sierra_Nevada = seq(min(ss_diff_obs$California_Sierra_Nevada$x), max(ss_diff_obs$California_Sierra_Nevada$x), 1)[which.min(abs(predict(ss_diff_obs$California_Sierra_Nevada, seq(min(ss_diff_obs$California_Sierra_Nevada$x), max(ss_diff_obs$California_Sierra_Nevada$x), 1))$y))], 
  S_Nevada = NA,
  Southern_Ghats_East = NA,
  Southern_Ghats_West = NA,
  Chile_Central = seq(min(ss_diff_obs$Chile_Central$x), max(ss_diff_obs$Chile_Central$x), 1)[which.min(abs(predict(ss_diff_obs$Chile_Central, seq(min(ss_diff_obs$Chile_Central$x), max(ss_diff_obs$Chile_Central$x), 1))$y))], 
  Western_Himalayas = seq(min(ss_diff_obs$Western_Himalayas$x), max(ss_diff_obs$Western_Himalayas$x), 1)[which.min(abs(predict(ss_diff_obs$Western_Himalayas, seq(min(ss_diff_obs$Western_Himalayas$x), max(ss_diff_obs$Western_Himalayas$x), 1))$y))],
  Hawaii = seq(min(ss_diff_obs$Hawaii$x), max(ss_diff_obs$Hawaii$x), 1)[which.min(abs(predict(ss_diff_obs$Hawaii, seq(min(ss_diff_obs$Hawaii$x), max(ss_diff_obs$Hawaii$x), 1))$y))],
  Taiwan_East = seq(min(ss_diff_obs$Taiwan_East$x), max(ss_diff_obs$Taiwan_East$x), 1)[which.min(abs(predict(ss_diff_obs$Taiwan_East, seq(min(ss_diff_obs$Taiwan_East$x), max(ss_diff_obs$Taiwan_East$x), 1))$y))], 
  Taiwan_West = seq(min(ss_diff_obs$Taiwan_West$x), max(ss_diff_obs$Taiwan_West$x), 1)[which.min(abs(predict(ss_diff_obs$Taiwan_West, seq(min(ss_diff_obs$Taiwan_West$x), max(ss_diff_obs$Taiwan_West$x), 1))$y))], 
  Costa_Rica_East = NA,
  Costa_Rica_West = NA))

## Boreal summer peak in richness
summer_peak_obs <- unlist(list(
  Swiss_Alps = seq(min(ss_diff_obs$Swiss_Alps$x), max(ss_diff_obs$Swiss_Alps$x), 1)[which.max(predict(ss_diff_obs$Swiss_Alps, seq(min(ss_diff_obs$Swiss_Alps$x), max(ss_diff_obs$Swiss_Alps$x), 1))$y)],
  Australian_Alps = NA,
  Blue_Ridge = seq(min(ss_diff_obs$Blue_Ridge$x), max(ss_diff_obs$Blue_Ridge$x), 1)[which.max(predict(ss_diff_obs$Blue_Ridge, seq(min(ss_diff_obs$Blue_Ridge$x), max(ss_diff_obs$Blue_Ridge$x), 1))$y)],
  Cascade_East = seq(min(ss_diff_obs$Cascade_East$x), max(ss_diff_obs$Cascade_East$x), 1)[which.max(predict(ss_diff_obs$Cascade_East, seq(min(ss_diff_obs$Cascade_East$x), max(ss_diff_obs$Cascade_East$x), 1))$y)],
  Cascade_West = seq(min(ss_diff_obs$Cascade_West$x), max(ss_diff_obs$Cascade_West$x), 1)[which.max(predict(ss_diff_obs$Cascade_West, seq(min(ss_diff_obs$Cascade_West$x), max(ss_diff_obs$Cascade_West$x), 1))$y)], 
  Coastal_California = seq(min(ss_diff_obs$Coastal_California$x), max(ss_diff_obs$Coastal_California$x), 1)[which.max(predict(ss_diff_obs$Coastal_California, seq(min(ss_diff_obs$Coastal_California$x), max(ss_diff_obs$Coastal_California$x), 1))$y)],
  Eastern_BC = seq(min(ss_diff_obs$Eastern_BC$x), max(ss_diff_obs$Eastern_BC$x), 1)[which.max(predict(ss_diff_obs$Eastern_BC, seq(min(ss_diff_obs$Eastern_BC$x), max(ss_diff_obs$Eastern_BC$x), 1))$y)], 
  Cord_Cantabrica = seq(min(ss_diff_obs$Cord_Cantabrica$x), max(ss_diff_obs$Cord_Cantabrica$x), 1)[which.max(predict(ss_diff_obs$Cord_Cantabrica, seq(min(ss_diff_obs$Cord_Cantabrica$x), max(ss_diff_obs$Cord_Cantabrica$x), 1))$y)],
  Cord_Merida = NA,
  Peru_SouthWest = NA, 
  Peru_East = NA,
  Argentina_North = NA,
  Ecuador_West = NA,
  Colombia_West = NA,
  Drakensberg = NA,
  Pyrenees_Atlantic = seq(min(ss_diff_obs$Pyrenees_Atlantic$x), max(ss_diff_obs$Pyrenees_Atlantic$x), 1)[which.max(predict(ss_diff_obs$Pyrenees_Atlantic, seq(min(ss_diff_obs$Pyrenees_Atlantic$x), max(ss_diff_obs$Pyrenees_Atlantic$x), 1))$y)],
  Pyrenees_Catalonia = seq(min(ss_diff_obs$Pyrenees_Catalonia$x), max(ss_diff_obs$Pyrenees_Catalonia$x), 1)[which.max(predict(ss_diff_obs$Pyrenees_Catalonia, seq(min(ss_diff_obs$Pyrenees_Catalonia$x), max(ss_diff_obs$Pyrenees_Catalonia$x), 1))$y)], 
  Serra_Mantiqueira = NA,
  Chiapas_East = NA,
  Veracruz = NA,
  Jalisco = NA,
  Sierra_Occidental = NA,
  Nuevo_Leon = seq(min(ss_diff_obs$Nuevo_Leon$x), max(ss_diff_obs$Nuevo_Leon$x), 1)[which.max(predict(ss_diff_obs$Nuevo_Leon, seq(min(ss_diff_obs$Nuevo_Leon$x), max(ss_diff_obs$Nuevo_Leon$x), 1))$y)], 
  California_Sierra_Nevada = seq(min(ss_diff_obs$California_Sierra_Nevada$x), max(ss_diff_obs$California_Sierra_Nevada$x), 1)[which.max(predict(ss_diff_obs$California_Sierra_Nevada, seq(min(ss_diff_obs$California_Sierra_Nevada$x), max(ss_diff_obs$California_Sierra_Nevada$x), 1))$y)], 
  S_Nevada = seq(min(ss_diff_obs$S_Nevada$x), max(ss_diff_obs$S_Nevada$x), 1)[which.max(predict(ss_diff_obs$S_Nevada, seq(min(ss_diff_obs$S_Nevada$x), max(ss_diff_obs$S_Nevada$x), 1))$y)],
  Southern_Ghats_East = NA,
  Southern_Ghats_West = NA,
  Chile_Central = NA, 
  Western_Himalayas = seq(min(ss_diff_obs$Western_Himalayas$x), max(ss_diff_obs$Western_Himalayas$x), 1)[which.max(predict(ss_diff_obs$Western_Himalayas, seq(min(ss_diff_obs$Western_Himalayas$x), max(ss_diff_obs$Western_Himalayas$x), 1))$y)], 
  Hawaii = NA,
  Taiwan_East = seq(min(ss_diff_obs$Taiwan_East$x), max(ss_diff_obs$Taiwan_East$x), 1)[which.max(predict(ss_diff_obs$Taiwan_East, seq(min(ss_diff_obs$Taiwan_East$x), max(ss_diff_obs$Taiwan_East$x), 1))$y)], 
  Taiwan_West = NA,
  Costa_Rica_East = NA,
  Costa_Rica_West = NA))  

## Boreal winter peak in richness
winter_peak_obs <- unlist(list(
  Swiss_Alps = NA,
  Australian_Alps = seq(min(ss_diff_obs$Australian_Alps$x), max(ss_diff_obs$Australian_Alps$x), 1)[which.min(predict(ss_diff_obs$Australian_Alps, seq(min(ss_diff_obs$Australian_Alps$x), max(ss_diff_obs$Australian_Alps$x), 1))$y)],
  Blue_Ridge = NA,
  Cascade_East = NA,
  Cascade_West = NA,
  Coastal_California = NA,
  Eastern_BC = NA,
  Cord_Cantabrica = NA,
  Cord_Merida = seq(min(ss_diff_obs$Cord_Merida$x), max(ss_diff_obs$Cord_Merida$x), 1)[which.min(predict(ss_diff_obs$Cord_Merida, seq(min(ss_diff_obs$Cord_Merida$x), max(ss_diff_obs$Cord_Merida$x), 1))$y)],
  Peru_SouthWest = seq(min(ss_diff_obs$Peru_SouthWest$x), max(ss_diff_obs$Peru_SouthWest$x), 1)[which.min(predict(ss_diff_obs$Peru_SouthWest, seq(min(ss_diff_obs$Peru_SouthWest$x), max(ss_diff_obs$Peru_SouthWest$x), 1))$y)], 
  Peru_East = NA,
  Argentina_North = seq(min(ss_diff_obs$Argentina_North$x), max(ss_diff_obs$Argentina_North$x), 1)[which.min(predict(ss_diff_obs$Argentina_North, seq(min(ss_diff_obs$Argentina_North$x), max(ss_diff_obs$Argentina_North$x), 1))$y)],
  Ecuador_West = NA,
  Colombia_West = seq(min(ss_diff_obs$Colombia_West$x), max(ss_diff_obs$Colombia_West$x), 1)[which.min(predict(ss_diff_obs$Colombia_West, seq(min(ss_diff_obs$Colombia_West$x), max(ss_diff_obs$Colombia_West$x), 1))$y)],
  Drakensberg = seq(min(ss_diff_obs$Drakensberg$x), max(ss_diff_obs$Drakensberg$x), 1)[which.min(predict(ss_diff_obs$Drakensberg, seq(min(ss_diff_obs$Drakensberg$x), max(ss_diff_obs$Drakensberg$x), 1))$y)], 
  Pyrenees_Atlantic = NA,
  Pyrenees_Catalonia = NA,
  Serra_Mantiqueira = NA,
  Chiapas_East = seq(min(ss_diff_obs$Chiapas_East$x), max(ss_diff_obs$Chiapas_East$x), 1)[which.min(predict(ss_diff_obs$Chiapas_East, seq(min(ss_diff_obs$Chiapas_East$x), max(ss_diff_obs$Chiapas_East$x), 1))$y)], 
  Veracruz = NA,
  Jalisco = seq(min(ss_diff_obs$Jalisco$x), max(ss_diff_obs$Jalisco$x), 1)[which.min(predict(ss_diff_obs$Jalisco, seq(min(ss_diff_obs$Jalisco$x), max(ss_diff_obs$Jalisco$x), 1))$y)], 
  Sierra_Occidental = seq(min(ss_diff_obs$Sierra_Occidental$x), max(ss_diff_obs$Sierra_Occidental$x), 1)[which.min(predict(ss_diff_obs$Sierra_Occidental, seq(min(ss_diff_obs$Sierra_Occidental$x), max(ss_diff_obs$Sierra_Occidental$x), 1))$y)],
  Nuevo_Leon = seq(min(ss_diff_obs$Nuevo_Leon$x), max(ss_diff_obs$Nuevo_Leon$x), 1)[which.min(predict(ss_diff_obs$Nuevo_Leon, seq(min(ss_diff_obs$Nuevo_Leon$x), max(ss_diff_obs$Nuevo_Leon$x), 1))$y)],
  California_Sierra_Nevada = NA,
  S_Nevada = NA,
  Southern_Ghats_East = seq(min(ss_diff_obs$Southern_Ghats_East$x), max(ss_diff_obs$Southern_Ghats_East$x), 1)[which.min(predict(ss_diff_obs$Southern_Ghats_East, seq(min(ss_diff_obs$Southern_Ghats_East$x), max(ss_diff_obs$Southern_Ghats_East$x), 1))$y)], 
  Southern_Ghats_West = seq(min(ss_diff_obs$Southern_Ghats_West$x), max(ss_diff_obs$Southern_Ghats_West$x), 1)[which.min(predict(ss_diff_obs$Southern_Ghats_West, seq(min(ss_diff_obs$Southern_Ghats_West$x), max(ss_diff_obs$Southern_Ghats_West$x), 1))$y)],
  Chile_Central = seq(min(ss_diff_obs$Chile_Central$x), max(ss_diff_obs$Chile_Central$x), 1)[which.min(predict(ss_diff_obs$Chile_Central, seq(min(ss_diff_obs$Chile_Central$x), max(ss_diff_obs$Chile_Central$x), 1))$y)],
  Western_Himalayas = NA, 
  Hawaii = NA,
  Taiwan_East = seq(min(ss_diff_obs$Taiwan_East$x), max(ss_diff_obs$Taiwan_East$x), 1)[which.min(predict(ss_diff_obs$Taiwan_East, seq(min(ss_diff_obs$Taiwan_East$x), max(ss_diff_obs$Taiwan_East$x), 1))$y)],
  Taiwan_West = NA,
  Costa_Rica_East = seq(min(ss_diff_obs$Costa_Rica_East$x), max(ss_diff_obs$Costa_Rica_East$x), 1)[which.min(predict(ss_diff_obs$Costa_Rica_East, seq(min(ss_diff_obs$Costa_Rica_East$x), max(ss_diff_obs$Costa_Rica_East$x), 1))$y)], 
  Costa_Rica_West = seq(min(ss_diff_obs$Costa_Rica_West$x), max(ss_diff_obs$Costa_Rica_West$x), 1)[which.min(predict(ss_diff_obs$Costa_Rica_West, seq(min(ss_diff_obs$Costa_Rica_West$x), max(ss_diff_obs$Costa_Rica_West$x), 1))$y)]))


###  Features in simulated patterns  ###

## Transition between net richness surpluses during the boreal winter and summer 
transition_elev_simu <- unlist(list(
  Swiss_Alps = NA,
  Australian_Alps = NA,
  Blue_Ridge = NA,
  Cascade_East = NA,
  Cascade_West = NA,
  Coastal_California = seq(min(ss_diff_simu$Coastal_California$x), 2000, 1)[which.min(abs(predict(ss_diff_simu$Coastal_California, seq(min(ss_diff_simu$Coastal_California$x), 2000, 1))$y))],
  Eastern_BC = NA,
  Cord_Cantabrica = NA,
  Cord_Merida = NA,
  Peru_SouthWest = NA,
  Peru_East = seq(min(ss_diff_simu$Peru_East$x), 4500, 1)[which.min(abs(predict(ss_diff_simu$Peru_East, seq(min(ss_diff_simu$Peru_East$x), 4500, 1))$y))],
  Argentina_North = NA, 
  Ecuador_West = NA,
  Colombia_West = NA,
  Drakensberg = NA, 
  Pyrenees_Atlantic = NA,
  Pyrenees_Catalonia = seq(min(ss_diff_simu$Pyrenees_Catalonia$x), max(ss_diff_simu$Pyrenees_Catalonia$x), 1)[which.min(abs(predict(ss_diff_simu$Pyrenees_Catalonia, seq(min(ss_diff_simu$Pyrenees_Catalonia$x), max(ss_diff_simu$Pyrenees_Catalonia$x), 1))$y))],
  Serra_Mantiqueira = NA, 
  Chiapas_East = seq(min(ss_diff_simu$Chiapas_East$x), 3000, 1)[which.min(abs(predict(ss_diff_simu$Chiapas_East, seq(min(ss_diff_simu$Chiapas_East$x), 3000, 1))$y))], 
  Veracruz = seq(min(ss_diff_simu$Veracruz$x), 3000, 1)[which.min(abs(predict(ss_diff_simu$Veracruz, seq(min(ss_diff_simu$Veracruz$x), 3000, 1))$y))], 
  Jalisco = NA, 
  Sierra_Occidental = NA, 
  Nuevo_Leon = seq(min(ss_diff_simu$Nuevo_Leon$x), 2000, 1)[which.min(abs(predict(ss_diff_simu$Nuevo_Leon, seq(min(ss_diff_simu$Nuevo_Leon$x), 2000, 1))$y))], 
  California_Sierra_Nevada = NA, 
  S_Nevada = NA, 
  Southern_Ghats_East = NA, 
  Southern_Ghats_West = NA, 
  Chile_Central = NA, 
  Western_Himalayas = seq(min(ss_diff_simu$Western_Himalayas$x), 2000, 1)[which.min(abs(predict(ss_diff_simu$Western_Himalayas, seq(min(ss_diff_simu$Western_Himalayas$x), 2000, 1))$y))], 
  Hawaii = seq(min(ss_diff_simu$Hawaii$x), 2500, 1)[which.min(abs(predict(ss_diff_simu$Hawaii, seq(min(ss_diff_simu$Hawaii$x), 2500, 1))$y))],
  Taiwan_East = seq(min(ss_diff_simu$Taiwan_East$x), max(ss_diff_simu$Taiwan_East$x), 1)[which.min(abs(predict(ss_diff_simu$Taiwan_East, seq(min(ss_diff_simu$Taiwan_East$x), max(ss_diff_simu$Taiwan_East$x), 1))$y))], 
  Taiwan_West = NA, 
  Costa_Rica_East = seq(min(ss_diff_simu$Costa_Rica_East$x), 2500, 1)[which.min(abs(predict(ss_diff_simu$Costa_Rica_East, seq(min(ss_diff_simu$Costa_Rica_East$x), 2500, 1))$y))],
  Costa_Rica_West = NA))

## Boreal summer peak in richness
summer_peak_simu <- unlist(list(
  Swiss_Alps = seq(min(ss_diff_simu$Swiss_Alps$x), max(ss_diff_simu$Swiss_Alps$x), 1)[which.max(predict(ss_diff_simu$Swiss_Alps, seq(min(ss_diff_simu$Swiss_Alps$x), max(ss_diff_simu$Swiss_Alps$x), 1))$y)],
  Australian_Alps = NA,
  Blue_Ridge = seq(300, max(ss_diff_simu$Blue_Ridge$x), 1)[which.max(predict(ss_diff_simu$Blue_Ridge, seq(300, max(ss_diff_simu$Blue_Ridge$x), 1))$y)], 
  Cascade_East = seq(min(ss_diff_simu$Cascade_East$x), max(ss_diff_simu$Cascade_East$x), 1)[which.max(predict(ss_diff_simu$Cascade_East, seq(min(ss_diff_simu$Cascade_East$x), max(ss_diff_simu$Cascade_East$x), 1))$y)], 
  Cascade_West = seq(min(ss_diff_simu$Cascade_West$x), max(ss_diff_simu$Cascade_West$x), 1)[which.max(predict(ss_diff_simu$Cascade_West, seq(min(ss_diff_simu$Cascade_West$x), max(ss_diff_simu$Cascade_West$x), 1))$y)], 
  Coastal_California = seq(min(ss_diff_simu$Coastal_California$x), max(ss_diff_simu$Coastal_California$x), 1)[which.max(predict(ss_diff_simu$Coastal_California, seq(min(ss_diff_simu$Coastal_California$x), max(ss_diff_simu$Coastal_California$x), 1))$y)],
  Eastern_BC = seq(min(ss_diff_simu$Eastern_BC$x), max(ss_diff_simu$Eastern_BC$x), 1)[which.max(predict(ss_diff_simu$Eastern_BC, seq(min(ss_diff_simu$Eastern_BC$x), max(ss_diff_simu$Eastern_BC$x), 1))$y)], 
  Cord_Cantabrica = seq(min(ss_diff_simu$Cord_Cantabrica$x), max(ss_diff_simu$Cord_Cantabrica$x), 1)[which.max(predict(ss_diff_simu$Cord_Cantabrica, seq(min(ss_diff_simu$Cord_Cantabrica$x), max(ss_diff_simu$Cord_Cantabrica$x), 1))$y)],
  Cord_Merida = NA,
  Peru_SouthWest = NA, 
  Peru_East = seq(min(ss_diff_simu$Peru_East$x), max(ss_diff_simu$Peru_East$x), 1)[which.max(predict(ss_diff_simu$Peru_East, seq(min(ss_diff_simu$Peru_East$x), max(ss_diff_simu$Peru_East$x), 1))$y)],
  Argentina_North = NA,
  Ecuador_West = NA,
  Colombia_West = NA,
  Drakensberg = NA,
  Pyrenees_Atlantic = seq(min(ss_diff_simu$Pyrenees_Atlantic$x), max(ss_diff_simu$Pyrenees_Atlantic$x), 1)[which.max(predict(ss_diff_simu$Pyrenees_Atlantic, seq(min(ss_diff_simu$Pyrenees_Atlantic$x), max(ss_diff_simu$Pyrenees_Atlantic$x), 1))$y)],
  Pyrenees_Catalonia = seq(min(ss_diff_simu$Pyrenees_Catalonia$x), max(ss_diff_simu$Pyrenees_Catalonia$x), 1)[which.max(predict(ss_diff_simu$Pyrenees_Catalonia, seq(min(ss_diff_simu$Pyrenees_Catalonia$x), max(ss_diff_simu$Pyrenees_Catalonia$x), 1))$y)], 
  Serra_Mantiqueira = NA,
  Chiapas_East = seq(min(ss_diff_simu$Chiapas_East$x), max(ss_diff_simu$Chiapas_East$x), 1)[which.max(predict(ss_diff_simu$Chiapas_East, seq(min(ss_diff_simu$Chiapas_East$x), max(ss_diff_simu$Chiapas_East$x), 1))$y)],
  Veracruz = NA, 
  Jalisco = NA, 
  Sierra_Occidental = NA,
  Nuevo_Leon = seq(min(ss_diff_simu$Nuevo_Leon$x), max(ss_diff_simu$Nuevo_Leon$x), 1)[which.max(predict(ss_diff_simu$Nuevo_Leon, seq(min(ss_diff_simu$Nuevo_Leon$x), max(ss_diff_simu$Nuevo_Leon$x), 1))$y)],
  California_Sierra_Nevada = seq(min(ss_diff_simu$California_Sierra_Nevada$x), max(ss_diff_simu$California_Sierra_Nevada$x), 1)[which.max(predict(ss_diff_simu$California_Sierra_Nevada, seq(min(ss_diff_simu$California_Sierra_Nevada$x), max(ss_diff_simu$California_Sierra_Nevada$x), 1))$y)], 
  S_Nevada = seq(min(ss_diff_simu$S_Nevada$x), max(ss_diff_simu$S_Nevada$x), 1)[which.max(predict(ss_diff_simu$S_Nevada, seq(min(ss_diff_simu$S_Nevada$x), max(ss_diff_simu$S_Nevada$x), 1))$y)],
  Southern_Ghats_East = NA, 
  Southern_Ghats_West = NA,
  Chile_Central = NA,
  Western_Himalayas = seq(min(ss_diff_simu$Western_Himalayas$x), max(ss_diff_simu$Western_Himalayas$x), 1)[which.max(predict(ss_diff_simu$Western_Himalayas, seq(min(ss_diff_simu$Western_Himalayas$x), max(ss_diff_simu$Western_Himalayas$x), 1))$y)], 
  Hawaii = seq(min(ss_diff_simu$Hawaii$x), max(ss_diff_simu$Hawaii$x), 1)[which.max(predict(ss_diff_simu$Hawaii, seq(min(ss_diff_simu$Hawaii$x), max(ss_diff_simu$Hawaii$x), 1))$y)],
  Taiwan_East = seq(min(ss_diff_simu$Taiwan_East$x), max(ss_diff_simu$Taiwan_East$x), 1)[which.max(predict(ss_diff_simu$Taiwan_East, seq(min(ss_diff_simu$Taiwan_East$x), max(ss_diff_simu$Taiwan_East$x), 1))$y)],
  Taiwan_West = NA, 
  Costa_Rica_East = seq(min(ss_diff_simu$Costa_Rica_East$x), max(ss_diff_simu$Costa_Rica_East$x), 1)[which.max(predict(ss_diff_simu$Costa_Rica_East, seq(min(ss_diff_simu$Costa_Rica_East$x), max(ss_diff_simu$Costa_Rica_East$x), 1))$y)],
  Costa_Rica_West = NA)) 

## Boreal winter peak in richness
winter_peak_simu <- unlist(list(
  Swiss_Alps = NA,
  Australian_Alps = seq(min(ss_diff_simu$Australian_Alps$x), max(ss_diff_simu$Australian_Alps$x), 1)[which.min(predict(ss_diff_simu$Australian_Alps, seq(min(ss_diff_simu$Australian_Alps$x), max(ss_diff_simu$Australian_Alps$x), 1))$y)],
  Blue_Ridge = NA,
  Cascade_East = NA,
  Cascade_West = NA,
  Coastal_California = NA,
  Eastern_BC = NA,
  Cord_Cantabrica = NA,
  Cord_Merida = seq(min(ss_diff_simu$Cord_Merida$x), max(ss_diff_simu$Cord_Merida$x), 1)[which.min(predict(ss_diff_simu$Cord_Merida, seq(min(ss_diff_simu$Cord_Merida$x), max(ss_diff_simu$Cord_Merida$x), 1))$y)],
  Peru_SouthWest = seq(min(ss_diff_simu$Peru_SouthWest$x), max(ss_diff_simu$Peru_SouthWest$x), 1)[which.min(predict(ss_diff_simu$Peru_SouthWest, seq(min(ss_diff_simu$Peru_SouthWest$x), max(ss_diff_simu$Peru_SouthWest$x), 1))$y)], 
  Peru_East = NA,
  Argentina_North = seq(min(ss_diff_simu$Argentina_North$x), max(ss_diff_simu$Argentina_North$x), 1)[which.min(predict(ss_diff_simu$Argentina_North, seq(min(ss_diff_simu$Argentina_North$x), max(ss_diff_simu$Argentina_North$x), 1))$y)],
  Ecuador_West = NA,
  Colombia_West = seq(min(ss_diff_simu$Colombia_West$x), max(ss_diff_simu$Colombia_West$x), 1)[which.min(predict(ss_diff_simu$Colombia_West, seq(min(ss_diff_simu$Colombia_West$x), max(ss_diff_simu$Colombia_West$x), 1))$y)],
  Drakensberg = seq(min(ss_diff_simu$Drakensberg$x), max(ss_diff_simu$Drakensberg$x), 1)[which.min(predict(ss_diff_simu$Drakensberg, seq(min(ss_diff_simu$Drakensberg$x), max(ss_diff_simu$Drakensberg$x), 1))$y)], 
  Pyrenees_Atlantic = NA,
  Pyrenees_Catalonia = NA,
  Serra_Mantiqueira = NA, 
  Chiapas_East = seq(min(ss_diff_simu$Chiapas_East$x), max(ss_diff_simu$Chiapas_East$x), 1)[which.min(predict(ss_diff_simu$Chiapas_East, seq(min(ss_diff_simu$Chiapas_East$x), max(ss_diff_simu$Chiapas_East$x), 1))$y)], 
  Veracruz = NA, 
  Jalisco = seq(min(ss_diff_simu$Jalisco$x), max(ss_diff_simu$Jalisco$x), 1)[which.min(predict(ss_diff_simu$Jalisco, seq(min(ss_diff_simu$Jalisco$x), max(ss_diff_simu$Jalisco$x), 1))$y)], 
  Sierra_Occidental = seq(min(ss_diff_simu$Sierra_Occidental$x), max(ss_diff_simu$Sierra_Occidental$x), 1)[which.min(predict(ss_diff_simu$Sierra_Occidental, seq(min(ss_diff_simu$Sierra_Occidental$x), max(ss_diff_simu$Sierra_Occidental$x), 1))$y)],
  Nuevo_Leon = NA, 
  California_Sierra_Nevada = NA,
  S_Nevada = NA,
  Southern_Ghats_East = NA,
  Southern_Ghats_West = seq(min(ss_diff_simu$Southern_Ghats_West$x), max(ss_diff_simu$Southern_Ghats_West$x), 1)[which.min(predict(ss_diff_simu$Southern_Ghats_West, seq(min(ss_diff_simu$Southern_Ghats_West$x), max(ss_diff_simu$Southern_Ghats_West$x), 1))$y)], 
  Chile_Central = seq(min(ss_diff_simu$Chile_Central$x), max(ss_diff_simu$Chile_Central$x), 1)[which.min(predict(ss_diff_simu$Chile_Central, seq(min(ss_diff_simu$Chile_Central$x), max(ss_diff_simu$Chile_Central$x), 1))$y)], 
  Western_Himalayas = NA, 
  Hawaii = NA, 
  Taiwan_East = seq(min(ss_diff_simu$Taiwan_East$x), max(ss_diff_simu$Taiwan_East$x), 1)[which.min(predict(ss_diff_simu$Taiwan_East, seq(min(ss_diff_simu$Taiwan_East$x), max(ss_diff_simu$Taiwan_East$x), 1))$y)], 
  Taiwan_West = seq(min(ss_diff_simu$Taiwan_West$x), max(ss_diff_simu$Taiwan_West$x), 1)[which.min(predict(ss_diff_simu$Taiwan_West, seq(min(ss_diff_simu$Taiwan_West$x), max(ss_diff_simu$Taiwan_West$x), 1))$y)],
  Costa_Rica_East = NA, 
  Costa_Rica_West = NA)) 


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
  upper_limit_simu = unlist(lapply(ss_diff_simu, function(x) max(x$x))))
rownames(data_diff_mnts) <- NULL

# Mountains located in the Northern Hemisphere
mountains_diff_data_north <- data_diff_mnts %>% filter(latitude > 1)
# Temperate mountains in the Southern Hemisphere
mountains_diff_data_south <- data_diff_mnts %>% filter(latitude <= -23)



## FIGURE S2

# Example of seasonal elevational range using Taiwan Yuhina

aa <- data.frame(elev = as.numeric(c(names(incidence_matrices_summer), names(incidence_matrices_winter))),
                 freq = c(species_fraction_checklist_summer[which(rownames(species_fraction_checklist_summer) == "Yuhina brunneiceps"),], species_fraction_checklist_winter[which(rownames(species_fraction_checklist_winter) == "Yuhina brunneiceps"),]),
                 season = c(rep("summer", ncol(species_fraction_checklist_summer)), rep("winter", ncol(species_fraction_checklist_winter))))
weighted.mean(aa$elev[aa$season == "summer"], aa$freq[aa$season == "summer"])
weighted.mean(aa$elev[aa$season == "winter"], aa$freq[aa$season == "winter"])
g_yuhina <- ggplot() +
  geom_point(data=aa, aes(x=freq, y=elev, col=season)) + scale_color_manual(values=c("coral2","deepskyblue2")) +
  theme_bw() + geom_vline(xintercept=0.05, linetype="dotted") + xlab("Rescaled frequency of observation") + ylab("Elevation") +
  geom_segment(aes(x=1.12, y=aa$elev[min(which(aa$season == "summer" & aa$freq > 0.05))] - 100, xend=1.12, yend=aa$elev[max(which(aa$season == "summer" & aa$freq > 0.05))] + 100), col="coral2", linewidth = 1) +
  geom_segment(aes(x=1.15, y=aa$elev[min(which(aa$season == "winter" & aa$freq > 0.05))] - 100, xend=1.15, yend=aa$elev[max(which(aa$season == "winter" & aa$freq > 0.05))] + 100), col="deepskyblue2", linewidth = 1) +
  geom_point(aes(x=1.12, y=wtd.mean(as.numeric(names(incidence_matrices_summer)), aa$freq[aa$season == "summer"])), col="coral2", size=3) +
  geom_point(aes(x=1.15, y=wtd.mean(as.numeric(names(incidence_matrices_winter)), aa$freq[aa$season == "winter"])), col="deepskyblue2", size=3) +
  theme_classic() + xlim(c(0,1.15)) + theme(legend.position = "none", text = element_text(size=15)) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) + ylim(c(0,3400))

# Plot species' seasonal elevational ranges for Taiwan East
species_elevational_ranges_sub <- species_elevational_ranges %>% filter(slope == slopes[k])
g_ranges_taiwan <- ggplot(species_elevational_ranges_sub, aes(x = reorder(species, -range_mean), y = range_mean)) + ylab("") + xlab("Species") +
  geom_pointrange(aes(ymin = range_min, ymax = range_max, colour = season),fatten = 1, size=1.5, alpha=0.7, position = position_dodge2(width=0.7)) + theme_bw() + scale_color_manual(values=c("coral2","deepskyblue2")) +
  theme(legend.position = c(0.85, 0.9), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=13), axis.text.x=element_blank(), axis.text.y=element_blank()) +
  ylim(c(0,3400))

# Plot seasonal difference in richness
g_richness_diff_taiwan <- ggplot() + 
  geom_hline(yintercept = 0, col="grey50") +
  geom_smooth(data=richness_all_df %>% filter(slope == slopes[k]), aes(x=elevation, y=richness_difference), col="black", se=F, method="loess") +
  geom_point(data=richness_all_df %>% filter(slope == slopes[k]), aes(x=elevation, y=richness_difference)) + theme_bw() +
  ylab("Seasonal difference in richness") + xlab("") + xlim(c(0,3400)) + coord_flip() + theme_classic() + theme(legend.position = "none", text = element_text(size=15), axis.text.y=element_blank()) +
  geom_segment(aes(y=30, x=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$winter_peak_obs, yend=30, xend=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$transition_obs), col="blue", linewidth = 1) +
  geom_segment(aes(y=30, x=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$transition_obs, yend=30, xend=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$summer_peak_obs), col="red", linewidth = 1) +
  geom_segment(aes(y=30, x=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$summer_peak_obs, yend=30, xend=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$upper_limit_obs), col="orange", linewidth = 1) +
  geom_point(aes(y=30, x=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$winter_peak_obs), col="blue", size=3) +
  geom_point(aes(y=30, x=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$transition_obs), col="black", size=3) +
  geom_point(aes(y=30, x=data_diff_mnts[data_diff_mnts$mountain_slope == slopes[k],]$summer_peak_obs), col="red", size=3)

pdf(file="results/figures/Fig_S2.pdf", width = 15, height = 5.8)
ggarrange(g_yuhina, g_ranges_taiwan, g_richness_diff_taiwan, labels = c("A", "B", "C"), align = "h", ncol=3, nrow=1, widths=c(0.7,1, 0.7))
dev.off()



############################################################################
##  Correspondence between elevational gradients and latitudinal gradient  ##
############################################################################

# Slope of transition line fit
# Observed
summary(lm(mountains_diff_data_north$transition_obs ~ mountains_diff_data_north$latitude)) # estimate = -100.17, R2 = 0.747  # estimate = -101.7, R2 = 0.79 # estimate = -111.04, R2 = 0.495  # estimate = -69.483, R2 = 0.7813. 

# Simulated
summary(lm(mountains_diff_data_north$transition_simu ~ mountains_diff_data_north$latitude)) # estimate = -58.55, R2 = 0.6268 # estimate = -84.21, R2 = 0.83 # estimate = -73.66, R2 = 0.79  # estimate = -69.483, R2 = 0.7813. 

cor.test(mountains_diff_data_north$transition_obs, mountains_diff_data_north$transition_simu)
cor.test(mountains_diff_data_north$winter_peak_obs, mountains_diff_data_north$winter_peak_simu)
cor.test(mountains_diff_data_north$summer_peak_obs, mountains_diff_data_north$summer_peak_simu)

# Distance between transition and peaks in mountains
# Observed
mean_elevation_distance_winterPeak_obs <- mean(mountains_diff_data_north$transition_obs - mountains_diff_data_north$winter_peak_obs, na.rm=T)
mean_elevation_distance_summerPeak_obs <- mean(mountains_diff_data_north$summer_peak_obs - mountains_diff_data_north$transition_obs, na.rm=T)
# Simulated
mean_elevation_distance_winterPeak_simu <- mean(mountains_diff_data_north$transition_simu - mountains_diff_data_north$winter_peak_simu, na.rm=T)
mean_elevation_distance_summerPeak_simu <- mean(mountains_diff_data_north$summer_peak_simu - mountains_diff_data_north$transition_simu, na.rm=T)

# Distance between transition and peaks in latitudinal gradient
lat_distance_winterPeak <- transition_global - winter_peak_global
lat_distance_summerPeak <- summer_peak_global - transition_global

# observed
mean_elevation_distance_winterPeak_obs / mean(lat_distance_winterPeak, na.rm=T) # for 1 degree lat -> 127.8892 metres elevation 
mean_elevation_distance_summerPeak_obs / mean(lat_distance_summerPeak, na.rm=T) # for 1 degree lat -> 86.61066 metres elevation

# simulated
mean_elevation_distance_winterPeak_simu / mean(lat_distance_winterPeak, na.rm=T) # for 1 degree lat -> 141.9008 metres elevation
mean_elevation_distance_summerPeak_simu / mean(lat_distance_summerPeak, na.rm=T) # for 1 degree lat -> 71.78396 metres elevation


# 1 degree along the latitudinal gradient = ~100m along the elevational gradient 



## FIGURE 4

g_ranges_obs <- ggplot(data_diff_mnts, aes(x = latitude, y = lower_limit_obs)) + ylab("Elevation") + xlab("Latitude") + ylim(c(0,6350)) + xlim(c(-40,50)) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_obs) == F), aes(ymin = lower_limit_obs, ymax = winter_peak_obs), col="blue", size=0.4) + 
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_obs) == F & is.na(transition_obs) == F), aes(ymin = winter_peak_obs, ymax = transition_obs), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_obs) == F & is.na(transition_obs) == T), aes(ymin = winter_peak_obs, ymax = upper_limit_obs), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_obs) == T & is.na(transition_obs) == F), aes(ymin = lower_limit_obs, ymax = transition_obs), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(transition_obs) == F & is.na(summer_peak_obs) == F), aes(ymin = transition_obs, ymax = summer_peak_obs), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(transition_obs) == F & is.na(summer_peak_obs) == T), aes(ymin = transition_obs, ymax = upper_limit_obs), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(transition_obs) == T & is.na(summer_peak_obs) == F), aes(ymin = lower_limit_obs, ymax = summer_peak_obs), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(summer_peak_obs) == F), aes(ymin = summer_peak_obs, ymax = upper_limit_obs), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(winter_peak_obs) == F), aes(ymin = winter_peak_obs, ymax = upper_limit_obs), col="blue", size=0.4) + 
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(winter_peak_obs) == F & is.na(transition_obs) == F), aes(ymin = transition_obs, ymax = winter_peak_obs), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(winter_peak_obs) == F & is.na(transition_obs) == T), aes(ymin = lower_limit_obs, ymax = winter_peak_obs), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(transition_obs) == F & is.na(summer_peak_obs) == F), aes(ymin = summer_peak_obs, ymax = transition_obs), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(summer_peak_obs) == F), aes(ymin = lower_limit_obs, ymax = summer_peak_obs), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(transition_obs) == F & is.na(summer_peak_obs) == T), aes(ymin = lower_limit_obs, ymax = transition_obs), col="red", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope %in% c("Serra_Mantiqueira", "Ecuador_West")), aes(ymin = lower_limit_obs, ymax = upper_limit_obs), col="grey40", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_East"), aes(ymin = lower_limit_obs, ymax = upper_limit_obs), col="red", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_SouthWest"), aes(ymin = lower_limit_obs, ymax = transition_obs), col="blue", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_SouthWest"), aes(ymin = transition_obs, ymax = upper_limit_obs), col="red", size=0.4) +
  geom_point(data=data_diff_mnts, aes(x=latitude, y=transition_obs), size=1.7) + 
  geom_point(data=data_diff_mnts, aes(x=latitude, y=winter_peak_obs), col="blue", size=1.2) + 
  geom_point(data=data_diff_mnts, aes(x=latitude, y=summer_peak_obs), col="red", size=1.2) + 
  geom_smooth(data=mountains_diff_data_north, aes(x=latitude, y=transition_obs), se=F, method = "lm", col="black", span=1) + theme_classic() +
  coord_cartesian(ylim = c(0,5900)) +
  labs(tag="A")

g_ranges_simu <- ggplot(data_diff_mnts, aes(x = latitude, y = lower_limit_simu)) + ylab("Elevation") + xlab("Latitude") + ylim(c(0,6350)) + xlim(c(-40,50)) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_simu) == F), aes(ymin = lower_limit_simu, ymax = winter_peak_simu), col="blue", size=0.4) + 
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_simu) == F & is.na(transition_simu) == F), aes(ymin = winter_peak_simu, ymax = transition_simu), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_simu) == F & is.na(transition_simu) == T), aes(ymin = winter_peak_simu, ymax = upper_limit_simu), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(winter_peak_simu) == T & is.na(transition_simu) == F), aes(ymin = lower_limit_simu, ymax = transition_simu), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(transition_simu) == F & is.na(summer_peak_simu) == F), aes(ymin = transition_simu, ymax = summer_peak_simu), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(transition_simu) == F & is.na(summer_peak_simu) == T), aes(ymin = transition_simu, ymax = upper_limit_simu), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(transition_simu) == T & is.na(summer_peak_simu) == F), aes(ymin = lower_limit_simu, ymax = summer_peak_simu), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_north %>% filter(is.na(summer_peak_simu) == F), aes(ymin = summer_peak_simu, ymax = upper_limit_simu), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(winter_peak_simu) == F), aes(ymin = winter_peak_simu, ymax = upper_limit_simu), col="blue", size=0.4) + 
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(winter_peak_simu) == F & is.na(transition_simu) == F), aes(ymin = transition_simu, ymax = winter_peak_simu), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(winter_peak_simu) == F & is.na(transition_simu) == T), aes(ymin = lower_limit_simu, ymax = winter_peak_simu), col="blue", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(transition_simu) == F & is.na(summer_peak_simu) == F), aes(ymin = summer_peak_simu, ymax = transition_simu), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(summer_peak_simu) == F), aes(ymin = lower_limit_simu, ymax = summer_peak_simu), col="red", size=0.4) +
  geom_linerange(data = mountains_diff_data_south %>% filter(is.na(transition_simu) == F & is.na(summer_peak_simu) == T), aes(ymin = lower_limit_simu, ymax = transition_simu), col="red", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope %in% c("Serra_Mantiqueira")), aes(ymin = lower_limit_simu, ymax = upper_limit_simu), col="grey40", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope %in% c("Ecuador_West")), aes(ymin = lower_limit_simu, ymax = upper_limit_simu), col="blue", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_East"), aes(ymin = lower_limit_simu, ymax = transition_simu), col="red", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_East"), aes(ymin = transition_simu, ymax = upper_limit_simu), col="blue", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_SouthWest"), aes(ymin = lower_limit_simu, ymax = winter_peak_simu), col="blue", size=0.4) +
  geom_linerange(data = data_diff_mnts %>% filter(mountain_slope == "Peru_SouthWest"), aes(ymin = winter_peak_simu, ymax = upper_limit_simu), col="blue", size=0.4) +
  geom_point(data=data_diff_mnts, aes(x=latitude, y=transition_simu), size=1.7) + 
  geom_point(data=data_diff_mnts, aes(x=latitude, y=winter_peak_simu), col="blue", size=1.2) + 
  geom_point(data=data_diff_mnts, aes(x=latitude, y=summer_peak_simu), col="red", size=1.2) + 
  geom_smooth(data=mountains_diff_data_north, aes(x=latitude, y=transition_simu), se=F, method = "lm", col="black", span=1) + theme_classic() +
  coord_cartesian(ylim = c(0,5900)) +
  labs(tag="B")

## Plot adjusted seasonal difference in richness 
adjustment <- (4627/100) - mean(c(transition_euro_africa, transition_americas, transition_east_asia))
richness_all_df <- richness_all_df %>% 
  mutate(lat_translated = latitude + (elevation / 100) - adjustment)
g_translated <- ggplot() + 
  geom_hline(yintercept = 0, col="grey50") + theme_classic() +
  geom_smooth(data=richness_all_df %>% filter(latitude > 3 & slope != "Hawaii"), aes(x=lat_translated, y=richness_difference, by=slope), col="grey60", se=F, alpha=0.5, linewidth = 0.8) +
  geom_smooth(data=richness_all_df %>% filter(latitude > 3 & slope != "Hawaii"), aes(x=lat_translated, y=richness_difference), col="black", se=F, linewidth = 1.8, alpha=0.5) +
  geom_smooth(data=data_global_rich, aes(x=LATITUDE, y=DIFFERENCE), col="green4", se=F, linewidth = 1.8) +
  coord_cartesian(ylim = c(-100,110)) +
  #geom_point(data=richness_all_df %>% filter(latitude > 3 & slope != "Hawaii"), aes(x=lat_translated, y=richness_difference, col=slope), alpha=0.1) +
  #scale_color_viridis_d("Mountain slopes") + 
  theme(legend.position = "none") + 
  xlab("Latitude") + ylab("Seadonal difference in richness") + xlim(c(-12,80)) + ylim(c(-136, 180)) + labs(tag="C")

pdf(file = "results/figures/Fig_4.pdf", bg = "white", width = 12, height = 4)
ggarrange(g_ranges_obs, g_ranges_simu, g_translated, nrow=1, ncol=3)
dev.off()


