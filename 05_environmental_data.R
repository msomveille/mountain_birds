## This script contains the code for extracting and processing environmental data ##

# Load relevant packages

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ncdf4)

# Load the data prepared in the script 04_sample_completeness.R
load("results/mnt_slopes_data_new.RData")
mountain_ranges <- vect("resources/mountain_ranges_2.shp")

# Load elevation data (downloaded from NASA Shuttle Radar Topography Mission V003)
elevation_raster <- terra::rast("resources/SRTM_2km.tif")
elevation_raster <- terra::project(elevation_raster, terra::vect(mountain_ranges[1,]))

elevation_bins <- seq(0, 7000, 200)
slopes <- mountain_slopes_polys$name

# Extract NDVI per elevational band
mois <- c("01","02","03","04","05","06","07","08","09","10","11","12")
jours <- c("01","11","21")
env_data <- vector()
for(k in 1:length(slopes)){
  # prep elevation
  elevation_raster_sub <- terra::as.polygons(crop(elevation_raster, ext(mountain_slopes_polys[k,])), values=T, aggregate=F)
  # load and prep ndvi
  # Summer
  ndvi.raster.summer <- rast()
  for(i in 4:9){
    for(j in 1:length(jours)){
      ndvi_file <- nc_open(paste0("/Users/msomveille/Desktop/Projects/global-migration-model/resources/NDVI/copernicus/c_gls_NDVI-LTS_1999-2017-", mois[i], jours[j], "_GLOBE_VGT-PROBAV_V2.2.1.nc"))
      ndvi.lat <- ncvar_get(nc = ndvi_file, varid = "lat")
      ndvi.lon <- ncvar_get(nc = ndvi_file, varid = "lon")
      ndvi.array <- ncvar_get(nc = ndvi_file, varid = "mean")
      nc_close(ndvi_file)
      ndvi.array[ndvi.array < 0] <- 0
      ndvi.array[ndvi.array > 0.92] <- NA
      ndvi.rast <- terra::rast(raster(t(ndvi.array), xmn=min(ndvi.lon), xmx=max(ndvi.lon), ymn=min(ndvi.lat), ymx=max(ndvi.lat), crs="+proj=longlat +datum=WGS84 +no_defs"))
      ndvi.rast <- crop(ndvi.rast, ext(mountain_slopes_polys[k,]))
      crs(ndvi.rast) <- "+proj=longlat +datum=WGS84 +no_defs"
      ndvi.raster.summer <- c(ndvi.raster.summer, ndvi.rast)
      rm(ndvi.array, ndvi.rast)
    }
  }
  ndvi.raster.summer <- mean(ndvi.raster.summer)
  # Winter
  ndvi.raster.winter <- rast()
  for(i in c(1:3, 10:12)){
    for(j in 1:length(jours)){
      ndvi_file <- nc_open(paste0("/Users/msomveille/Desktop/Projects/global-migration-model/resources/NDVI/copernicus/c_gls_NDVI-LTS_1999-2017-", mois[i], jours[j], "_GLOBE_VGT-PROBAV_V2.2.1.nc"))
      ndvi.lat <- ncvar_get(nc = ndvi_file, varid = "lat")
      ndvi.lon <- ncvar_get(nc = ndvi_file, varid = "lon")
      ndvi.array <- ncvar_get(nc = ndvi_file, varid = "mean")
      nc_close(ndvi_file)
      ndvi.array[ndvi.array < 0] <- 0
      ndvi.array[ndvi.array > 0.92] <- NA
      ndvi.rast <- terra::rast(raster(t(ndvi.array), xmn=min(ndvi.lon), xmx=max(ndvi.lon), ymn=min(ndvi.lat), ymx=max(ndvi.lat), crs="+proj=longlat +datum=WGS84 +no_defs"))
      ndvi.rast <- crop(ndvi.rast, ext(mountain_slopes_polys[k,]))
      crs(ndvi.rast) <- "+proj=longlat +datum=WGS84 +no_defs"
      ndvi.raster.winter <- c(ndvi.raster.winter, ndvi.rast)
      rm(ndvi.array, ndvi.rast)
    }
  }
  ndvi.raster.winter <- mean(ndvi.raster.winter)
  
  # load and prep temperature
  ## Extract climate for each ecoregion. Climate was dowloaded from Chelsa (freely available)
  years <- 2000:2018
  months_days <- c("001", "032", "060", "091", "121", "152", "182", "213", "244", "274", "305", "335")
  Temp_files <- list.files("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas")
  Temp_files <- Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[4])) %in% years]
  
  mois_summer <- mois[4:9]
  mois_winter <- mois[c(1:3, 10:12)]
  temp_summer <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% mois_summer])
  temp.raster.summer <- terra::rast(apply(temp_summer, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas/", x)))) # load temperature rasters
  temp.raster.summer <- (mean(temp.raster.summer) / 10) - 273.15
  temp_winter <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% mois_winter])
  temp.raster.winter <- terra::rast(apply(temp_winter, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas/", x)))) # load temperature rasters
  temp.raster.winter <- (mean(temp.raster.winter) / 10) - 273.15
  crs(temp.raster.summer) <- crs(temp.raster.winter) <- "+proj=longlat +datum=WGS84"
  
  # put together the environmental data
  env_data <- rbind(env_data, 
                    data.frame(slope = slopes[k],
                               season = "summer",
                               elevation = elevation_raster_sub$SRTM_global,
                               ndvi = terra::extract(ndvi.raster.summer, elevation_raster_sub, fun=mean, na.rm=T)$layer,
                               temperature = terra::extract(temp.raster.summer, elevation_raster_sub, fun=mean, na.rm-T)$mean)
  )
  env_data <- rbind(env_data, 
                    data.frame(slope = slopes[k],
                               season = "winter",
                               elevation = elevation_raster_sub$SRTM_global,
                               ndvi = terra::extract(ndvi.raster.winter, elevation_raster_sub, fun=mean, na.rm=T)$layer,
                               temperature = terra::extract(temp.raster.winter, elevation_raster_sub, fun=mean, na.rm-T)$mean)
  )
  rm(ndvi.array, ndvi.raster.summer, ndvi.raster.winter, temp.raster.summer, temp.raster.winter)
  print(k)
}
env_data$ndvi[is.na(env_data$ndvi) == T] <- 0

write.csv(env_data, "results/env_data.csv")



#####  Eastern Himalayas  #####

mountain_range_E_Him <- mountain_ranges[mountain_ranges$NAME == "EASTERN HIMALAYAS",]

# Load the mountain slopes polygons and ebird data prepared using the script 03_mountain_slopes_delineation.R
load("results/mountain_slopes.RData")

# Extract environmental conditions per elevational band

elevation_bins <- seq(0, 7000, 200)
mois <- c("01","02","03","04","05","06","07","08","09","10","11","12")
jours <- c("01","11","21")
env_data_E_Him <- vector()
# prep elevation
elevation_raster_sub <- terra::as.polygons(crop(elevation_raster, ext(mountain_slopes_polygons$Eastern_Himalayas)), values=T, aggregate=F)
# load and prep ndvi
# Summer
ndvi.raster.summer <- rast()
for(i in 4:9){
  for(j in 1:length(jours)){
    ndvi_file <- nc_open(paste0("/Users/msomveille/Desktop/Projects/global-migration-model/resources/NDVI/copernicus/c_gls_NDVI-LTS_1999-2017-", mois[i], jours[j], "_GLOBE_VGT-PROBAV_V2.2.1.nc"))
    ndvi.lat <- ncvar_get(nc = ndvi_file, varid = "lat")
    ndvi.lon <- ncvar_get(nc = ndvi_file, varid = "lon")
    ndvi.array <- ncvar_get(nc = ndvi_file, varid = "mean")
    nc_close(ndvi_file)
    ndvi.array[ndvi.array < 0] <- 0
    ndvi.array[ndvi.array > 0.92] <- NA
    ndvi.rast <- terra::rast(raster(t(ndvi.array), xmn=min(ndvi.lon), xmx=max(ndvi.lon), ymn=min(ndvi.lat), ymx=max(ndvi.lat), crs="+proj=longlat +datum=WGS84 +no_defs"))
    ndvi.rast <- crop(ndvi.rast, ext(mountain_slopes_polygons$Eastern_Himalayas))
    crs(ndvi.rast) <- "+proj=longlat +datum=WGS84 +no_defs"
    ndvi.raster.summer <- c(ndvi.raster.summer, ndvi.rast)
    rm(ndvi.array, ndvi.rast)
  }
}
ndvi.raster.summer <- mean(ndvi.raster.summer)
# Winter
ndvi.raster.winter <- rast()
for(i in c(1:3, 10:12)){
  for(j in 1:length(jours)){
    ndvi_file <- nc_open(paste0("/Users/msomveille/Desktop/Projects/global-migration-model/resources/NDVI/copernicus/c_gls_NDVI-LTS_1999-2017-", mois[i], jours[j], "_GLOBE_VGT-PROBAV_V2.2.1.nc"))
    ndvi.lat <- ncvar_get(nc = ndvi_file, varid = "lat")
    ndvi.lon <- ncvar_get(nc = ndvi_file, varid = "lon")
    ndvi.array <- ncvar_get(nc = ndvi_file, varid = "mean")
    nc_close(ndvi_file)
    ndvi.array[ndvi.array < 0] <- 0
    ndvi.array[ndvi.array > 0.92] <- NA
    ndvi.rast <- terra::rast(raster(t(ndvi.array), xmn=min(ndvi.lon), xmx=max(ndvi.lon), ymn=min(ndvi.lat), ymx=max(ndvi.lat), crs="+proj=longlat +datum=WGS84 +no_defs"))
    ndvi.rast <- crop(ndvi.rast, ext(mountain_slopes_polygons$Eastern_Himalayas))
    crs(ndvi.rast) <- "+proj=longlat +datum=WGS84 +no_defs"
    ndvi.raster.winter <- c(ndvi.raster.winter, ndvi.rast)
    rm(ndvi.array, ndvi.rast)
  }
}
ndvi.raster.winter <- mean(ndvi.raster.winter)

# load and prep temperature
## Extract climate for each ecoregion. Climate was dowloaded from Chelsa (freely available)
years <- 2000:2018
months_days <- c("001", "032", "060", "091", "121", "152", "182", "213", "244", "274", "305", "335")
Temp_files <- list.files("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas")
Temp_files <- Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[4])) %in% years]
mois_summer <- mois[4:9]
mois_winter <- mois[c(1:3, 10:12)]
temp_summer <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% mois_summer])
temp.raster.summer <- terra::rast(apply(temp_summer, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas/", x)))) # load temperature rasters
temp.raster.summer <- (mean(temp.raster.summer) / 10) - 273.15
temp_winter <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% mois_winter])
temp.raster.winter <- terra::rast(apply(temp_winter, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas/", x)))) # load temperature rasters
temp.raster.winter <- (mean(temp.raster.winter) / 10) - 273.15
crs(temp.raster.summer) <- crs(temp.raster.winter) <- "+proj=longlat +datum=WGS84"

# put together the environmental data
env_data_E_Him <- rbind(env_data_E_Him, 
                        data.frame(slope = slopes[k],
                                   season = "summer",
                                   elevation = elevation_raster_sub$SRTM_global,
                                   ndvi = terra::extract(ndvi.raster.summer, elevation_raster_sub, fun=mean, na.rm=T)$layer,
                                   temperature = terra::extract(temp.raster.summer, elevation_raster_sub, fun=mean, na.rm-T)$mean)
)
env_data_E_Him <- rbind(env_data_E_Him, 
                        data.frame(slope = slopes[k],
                                   season = "winter",
                                   elevation = elevation_raster_sub$SRTM_global,
                                   ndvi = terra::extract(ndvi.raster.winter, elevation_raster_sub, fun=mean, na.rm=T)$layer,
                                   temperature = terra::extract(temp.raster.winter, elevation_raster_sub, fun=mean, na.rm-T)$mean)
)
rm(ndvi.array, ndvi.raster.summer, ndvi.raster.winter, temp.raster.summer, temp.raster.winter)

env_data_E_Him$ndvi[is.na(env_data$ndvi) == T] <- 0

write.csv(env_data_E_Him, "results/env_data_E_Him.csv")


# For future climate, the script above can be repeated with simply changing the temperature rasters that are loaded, and instead loading the ones downlowded for a future scenario


