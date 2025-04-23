## This script contains the code for quantifying empirical seasonal thermal tracking ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggpubr)

# Load the data prepared in the previous scripts
load("mnt_slopes_data_new.RData")
mountain_ranges <- vect("mountain_ranges_2.shp")
species_elevational_ranges <- read.csv("species_elevational_ranges.csv")


species_per_slope <- table(species_elevational_ranges$species, species_elevational_ranges$slope)
slopes <- unique(species_elevational_ranges$slope)
slopes2 <- colnames(species_per_slope)

# Calculate altitudinal migration distance and seasonal thermal difference
migration_thermal_distances <- list()
for(k in 1:length(slopes2)){
  # Select resident species
  resident_species <- names(which(species_per_slope[,k] == 2))
  species_elevational_ranges_sel <- species_elevational_ranges %>% filter(slope == slopes2[k]) %>% filter(species %in% resident_species)
  # Migration distance
  migra_dist <- abs(species_elevational_ranges_sel$range_mean[species_elevational_ranges_sel$season == "summer"] - species_elevational_ranges_sel$range_mean[species_elevational_ranges_sel$season == "winter"])
  # Thermal distance
  thermal_dist <- abs(species_elevational_ranges_sel$temp_mean[species_elevational_ranges_sel$season == "summer"] - species_elevational_ranges_sel$temp_mean[species_elevational_ranges_sel$season == "winter"])
  
  migration_thermal_distances[[k]] <- data.frame(slope = slopes2[k], 
                                                 species = resident_species, 
                                                 migration_distance = migra_dist,
                                                 thermal_distance = thermal_dist)
}
migration_thermal_distances_all <- do.call(rbind, migration_thermal_distances)

# Latitude of the mountain slopes
lat_mnt <- mountain_slopes_polys %>% st_centroid() %>% st_coordinates() %>% as.data.frame() %>% dplyr::select("Y") %>% rename(latitude = Y)
lat_mnt <- cbind(slopes, lat_mnt) %>% rename(slope = slopes)
migration_thermal_distances_all <- migration_thermal_distances_all %>% left_join(lat_mnt, by="slope")

g_hist_migr <- ggplot(migration_thermal_distances_all) + 
  geom_histogram(aes(x=sqrt(migration_distance))) +
  theme_classic() + xlim(c(0,60)) + xlab("Altitudinal migration distance (square-root)") + ylab("Number of avian populations") + labs(tag="A") +
  geom_vline(xintercept = sqrt(200),linetype="dotted")


# How many altitudinal migrants are doing worst than if they stayed resident in terms of seasonal thermal tracking?

alt_migr_results <- list()
for(k in 1:length(slopes2)){
  mean_thermal_dist_residents <- mean(migration_thermal_distances[[k]]$thermal_distance[migration_thermal_distances[[k]]$migration_distance < 200])
  alt_migr_results[[k]] <- data.frame(slope = slopes2[k],
                                      richness_total = nrow(migration_thermal_distances[[k]]),
                                      richness_alt_migr = length(which(migration_thermal_distances[[k]]$migration_distance >= 200)),
                                      richness_alt_migr_upslope = length(which(migration_thermal_distances[[k]]$thermal_distance[migration_thermal_distances[[k]]$migration_distance >= 200] >= mean_thermal_dist_residents)))
}
alt_migr_results <- do.call(rbind, alt_migr_results) %>%
  mutate(prop_alt_migr = richness_alt_migr / richness_total) %>%
  mutate(prop_uplope = richness_alt_migr_upslope / richness_alt_migr) %>% 
  left_join(lat_mnt)

# Total proportion of migrants
sum(alt_migr_results$richness_alt_migr) / sum(alt_migr_results$richness_total)
# Proportion of migrants versus latitude
ggplot(data = alt_migr_results) +
  geom_point(aes(x=abs(latitude), y=prop_alt_migr))
cor.test(alt_migr_results$latitude, alt_migr_results$prop_alt_migr)
# Proportion of altitudinal migrants in the equatorial tropics
sum(alt_migr_results$richness_alt_migr[abs(alt_migr_results$latitude) < 10]) / sum(alt_migr_results$richness_total[abs(alt_migr_results$latitude) < 10])
# Proportion of altitudinal migrants that are doing worse than residents in the equatorial tropics
sum(alt_migr_results$richness_alt_migr_upslope[abs(alt_migr_results$latitude) < 10]) / sum(alt_migr_results$richness_alt_migr[abs(alt_migr_results$latitude) < 10])
# Proportion of altitudinal migrants at mid-latitude
sum(alt_migr_results$richness_alt_migr[abs(alt_migr_results$latitude) > 15 & abs(alt_migr_results$latitude) < 35]) / sum(alt_migr_results$richness_total[abs(alt_migr_results$latitude) > 15 & abs(alt_migr_results$latitude) < 35])
# Proportion of altitudinal migrants that are doing worse than residents at mid-latitude
sum(alt_migr_results$richness_alt_migr_upslope[abs(alt_migr_results$latitude) > 15 & abs(alt_migr_results$latitude) < 35]) / sum(alt_migr_results$richness_alt_migr[abs(alt_migr_results$latitude) > 15 & abs(alt_migr_results$latitude) < 35])
# Proportion of altitudinal migrants in temperate regions
sum(alt_migr_results$richness_alt_migr[abs(alt_migr_results$latitude) > 35]) / sum(alt_migr_results$richness_total[abs(alt_migr_results$latitude) > 35])

# Proportion of altitudinal migrants that are doing worse than residents globally
sum(alt_migr_results$richness_alt_migr_upslope) / sum(alt_migr_results$richness_alt_migr)



# Frequency distribution of the seasonal thermal tracking compared to sedentary for altitudinal migrants
seasonal_climate_tracking_spp <- vector()
for(k in 1:length(slopes2)){
  mean_thermal_dist_residents <- mean(migration_thermal_distances[[k]]$thermal_distance[migration_thermal_distances[[k]]$migration_distance < 200])
  seasonal_climate_tracking_spp <-c(seasonal_climate_tracking_spp, mean_thermal_dist_residents - migration_thermal_distances[[k]]$thermal_distance[migration_thermal_distances[[k]]$migration_distance > 200])
}
g_hist_tracking <- ggplot() + 
  geom_histogram(aes(x=seasonal_climate_tracking_spp)) +
  theme_classic() + xlab("Seasonal thermal tracking") + ylab("Number of avian populations") + labs(tag="B") +
  geom_vline(xintercept = 0,linetype="dotted")


# Plot Figure 2

pdf("results/figures/Fig_2.pdf", width=9, height=4)
ggarrange(g_hist_migr, g_hist_tracking, nrow=1, ncol=2)
dev.off()

