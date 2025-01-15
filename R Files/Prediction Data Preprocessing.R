################################################################################
#
# Produce a dataframe of raster data for predicting across the study area
#
################################################################################

# Load necessary libraries
library(sf)
library(tidyverse)
library(terra)
library(raster)

# Read in raster files
covars.name <- read.csv("Data/BoF_mw_sh_covariate_list_rasterfolder.csv")
covars.name$FilenameAll <- paste0(covars.name$Name, covars.name$file.extension)
ras_list <- paste0(getwd(), "/Data/Rasters/", covars.name$FilenameAll)
r <- lapply(ras_list, terra::rast)

# Produce an average bottom temp. raster for each year (by averaging across
# individual rasters for June, July, August)
mean_bt_rasters <- c()
years <- 2011:2023
for(i in 1:length(years))
{
  year_bt_stack <- stack(c(r[[6+3*(i-1)]], r[[7+3*(i-1)]], r[[8+3*(i-1)]]))
  mean_bt_raster <- calc(year_bt_stack, fun = mean)
  mean_bt_rasters <- c(mean_bt_rasters, mean_bt_raster)
}

# Crop bottom temp. average rasters and static rasters to the area covered by the 
# benthoscape & 2021+ bottom temp. data, and then aggregate the resulting rasters
# to approx. 1 km resolution
benthoscape <- raster(r[[5]])
btm_temp_june_2021 <- raster(r[[37]])
unprocessed_rasters <- c(raster(r[[1]]), raster(r[[2]]), raster(r[[3]]),
                         raster(r[[4]]), raster(r[[5]]), mean_bt_rasters)
processed_rasters <- list()
for (i in 1:length(unprocessed_rasters))
{
  processed_rasters[[i]] <- unprocessed_rasters[[i]] %>%
    crop(benthoscape) %>% mask(benthoscape) %>%
    resample(btm_temp_june_2021) %>%
    crop(btm_temp_june_2021) %>% mask(btm_temp_june_2021) %>%
    aggregate(fact = 1000/res(r[[1]])[1], na.rm = FALSE, fun = modal)
}

# Produce a dataframe from the processed raster files
raster_df_init <- data.frame()
for (i in 1:length(years))
{
  pred_raster_data <- stack(c(processed_rasters[1:5], processed_rasters[5+i]))
  names(pred_raster_data) <- c("Bathymetry", "SedMobFreq", "ShearVelocity",
                               "Backscatter", "Benthoscape", "BottomTempAtSurvey")
  raster_df_new <- as.data.frame(na.omit(rasterToPoints(pred_raster_data)))
  raster_df_new$year <- years[i]
  raster_df_new$f_year <- as.factor(raster_df_new$year)
  raster_df_init <- rbind(raster_df_init, raster_df_new)
}

# Filter out 2020 rows (as no survey was conducted that year), filter out bedrocks 
# and boulders part of benthoscape (i.e., benthoscape = 1, as area is tiny relative 
# to scale of survey), and then create log versions of some covariates and remove any 
# NAs
raster_df <- raster_df_init %>% 
  filter(year != 2020) %>%
  mutate(LogBathymetry = log(-Bathymetry), 
         LogBottomTempAtSurvey = log(BottomTempAtSurvey),
         LogShearVelocity = log(ShearVelocity),
         Bathymetry = -Bathymetry, # Make depths positive
         Benthoscape = round(Benthoscape)) %>%
  filter(!Benthoscape %in% c(1,4)) %>% na.omit()

# Benthoscape: merge silt with silty gravel with anemones, then make benthoscape 
# a factor
raster_df$Benthoscape <- ifelse(raster_df$Benthoscape == 7, 6, raster_df$Benthoscape)
benthoscape_labels <- c("Gravelly sand", "Mixed sediment", 
                        "Sand", "Silt", "Tidal scoured mixed sediments")
raster_df$Benthoscape <- factor(raster_df$Benthoscape, levels = c(2:3, 5, 6, 8),
                                benthoscape_labels)

# Convert to sf object
raster_df[c("x.copy", "y.copy")] <- raster_df[c("x","y")]
raster_df <- st_as_sf(raster_df, coords= c("x.copy","y.copy"),  crs = 32620)

# Read in the SPAs shapefiles
temp <- tempfile()
download.file(paste0("https://raw.githubusercontent.com/Mar-scal/GIS_layers/mast",
                     "er/inshore_boundaries/inshore_boundaries.zip"), temp)
temp2 <- tempfile()
unzip(zipfile=temp, exdir=temp2)
SPA1A <- st_read(paste0(temp2, "/SPA1A_polygon_NAD83.shp")) %>% 
  mutate(SPA = "SPA 1A") %>% dplyr::select(SPA)
SPA1B <- st_read(paste0(temp2, "/SPA1B_polygon_NAD83.shp")) %>% 
  mutate(SPA = "SPA 1B") %>% dplyr::select(SPA)
SPA4 <- st_read(paste0(temp2, "/SPA4_polygon_NAD83.shp"))  %>% 
  mutate(SPA = "SPA 4") %>% dplyr::select(SPA)

# Combine the SPA shapefiles into a single UTM Zone 20N object and then intersect
# with the raster data
SPA1ABand4 <- rbind(SPA1A, SPA1B, SPA4) %>% st_transform(crs = 32620)
raster_df_SPA <- st_intersection(SPA1ABand4, raster_df)

# Write processed dataframe to rds file so can use it during prediction
saveRDS(raster_df_SPA, paste(getwd(),"/Data/Processed_Data_For_Prediction.rds", sep = ""))