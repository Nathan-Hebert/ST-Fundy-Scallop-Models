################################################################################
#
# Prepare survey data and covariates for modelling
#
################################################################################

# Load necessary libraries
library(sf)
library(tidyverse)
library(terra)
library(mapview)
library(raster)
library(lubridate)

##################################Survey Data###################################

# Notes: 
# *** All cruise.tow (IDs) in meat weight - shell height (MWSH) sampling dataset 
# (mwsh.dat) should have a corresponding matching tow in the shell height 
# frequency (SHF) dataset (shf.dat). 
# *** Only a subsample of tows are sampled for meat weight-shell height, so 
# therefore there will be more tows with shell height data than meat weight - shell 
# height data. 
# *** Only a subset of scallops in a tow is sampled for the meat weight-shell 
# height data within a sampled tow.
# ***If one wanted to get commercial size (fishery-targeted size), one would use 
# shf.dat bins from BIN_ID_80 and up; i.e. sum fields BIN_ID_80... up to BIN_ID_195.

# Set seed
set.seed(123)

# Read in the SPA shapefiles for SPAs 1A, 1B, and 4
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

# Combine the SPA shapefiles into a single UTM Zone 20N object
SPA1ABand4 <- rbind(SPA1A, SPA1B, SPA4) %>% st_transform(crs = 4326)

####SHF Data####

# Load in dataset with abundance by 5mm shell height frequency bin 
shf.dat <- read.csv('Data/bof.shf.unlined.gear.2011to2023.JunetoAug.csv')
str(shf.dat)
dim(shf.dat)

# Convert to an sf object, and then crop to SPAs 1A, 1B, and 4... data comes in 
# as WGS 84 
shf.dat.sf <-  st_as_sf(shf.dat, coords= c("mid.lon","mid.lat"),  crs = 4326) %>%
  st_intersection(st_make_valid(SPA1ABand4))
mapview(shf.dat.sf)

# Convert cropped sf object to WGS 84 UTM zone 20N
shf.dat.sf.utm <- st_transform(shf.dat.sf, crs = 32620)
shf.dat.sf.utm[c("x","y")] <- st_coordinates(shf.dat.sf.utm$geometry)
st_crs(shf.dat.sf.utm)

# Calculate commercial numbers in thousands per km^2 by summing the commercial bins,
# applying the q correction, and converting units
commercial_bins <- names(shf.dat.sf.utm)[grep("BIN", 
                                              names(shf.dat.sf.utm))[-c(1:16)]]
shf.dat.sf.utm$nums_tow <- (rowSums(shf.dat.sf.utm[,commercial_bins, drop = T])/
  ((800*5.334)/(1000)))/(0.36) # 1 standardized tow = 800m*5.334m, q = 0.36

# Make year factor variable and day of year variable
shf.dat.sf.utm$f_year <- as.factor(shf.dat.sf.utm$year)
shf.dat.sf.utm$doy <- lubridate::yday(make_date(year = shf.dat.sf.utm$year, 
                                                month = shf.dat.sf.utm$month, 
                                                day = shf.dat.sf.utm$day))

####MWSH Data####

# Load in meat weight-shell height (MWSH) sampling data and subset to commercial 
# size scallops... the data is normalized by individual scallop sample within a 
# tow within a cruise 
mwsh.dat <- read.csv("Data/BFdetail.dat.2011to2023.JunetoMidAug.clean.csv") %>%
  filter(HEIGHT >= 80)
str(mwsh.dat)
dim(mwsh.dat)

# Convert to an sf object, and then crop to SPAs 1A, 1B, and 4... data comes in 
# as WGS 84 
mwsh.dat.sf <- st_as_sf(mwsh.dat, coords= c("mid.lon","mid.lat"),  crs = 4326) %>%
  st_intersection(st_make_valid(SPA1ABand4))
mapview(mwsh.dat.sf)

# Convert cropped sf object to WGS 84 UTM zone 20N
mwsh.dat.sf.utm <- st_transform(mwsh.dat.sf, crs = 32620)
mwsh.dat.sf.utm[c("x","y")] <- st_coordinates(mwsh.dat.sf.utm$geometry)
st_crs(mwsh.dat.sf.utm)

# Make year factor variable and day of year variable
mwsh.dat.sf.utm$f_year <- as.factor(mwsh.dat.sf.utm$year)
mwsh.dat.sf.utm$doy <- lubridate::yday(make_date(year = mwsh.dat.sf.utm$year, 
                                                 month = mwsh.dat.sf.utm$month, 
                                                 day = mwsh.dat.sf.utm$day))

#########################Add Environmental Covariates########################### 

# Read in covariate raster files
covars.name <- read.csv("Data/BoF_mw_sh_covariate_list_rasterfolder.csv")
covars.name$FilenameAll <- paste0(covars.name$Name, covars.name$file.extension)
ras_list <- paste0(getwd(), "/Data/Rasters/", covars.name$FilenameAll)
r <- lapply(ras_list, terra::rast)
names(r[[3]]) <- "ShearVelocity"
r

# Intersect static environmental data with MWSH and SHF data 
mwsh_df <- mwsh.dat.sf.utm 
shf_df <- shf.dat.sf.utm
for (i in 1:5)
{
  mwsh_df[, names(r[[i]]) ] <- raster::extract(r[[i]], mwsh_df)[,2]
  shf_df[, names(r[[i]]) ] <- raster::extract(r[[i]], shf_df)[,2]
}

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

# Intersect each MWSH and SHF data with the average bottom temperature over 
# June-August
mwsh_df$BottomTempAtSurvey <- NA
shf_df$BottomTempAtSurvey <- NA
for (j in 1:length(years)){
  temp.year <- 2010 + j
  if (temp.year %in% shf_df$year)
  {
    val_mwsh <- raster::extract(mean_bt_rasters[[j]], (mwsh_df %>% 
                                           filter(year == temp.year)))
    val_sh <- raster::extract(mean_bt_rasters[[j]], (shf_df %>% filter(year == temp.year)))
    mwsh_df$BottomTempAtSurvey[(mwsh_df$year == temp.year)] <- val_mwsh
    shf_df$BottomTempAtSurvey[(shf_df$year == temp.year)] <- val_sh
  }
}

# Make depth positive
mwsh_df <- mwsh_df %>% mutate(Bathymetry = -Bathymetry)
shf_df <- shf_df %>% mutate(Bathymetry = -Bathymetry)

# Remove tows that occurred where benthoscape or backscatter are NA, or where 
# temperature data is not available for the entire time series... also filter out 
# bedrocks and boulders part of benthoscape (i.e., benthoscape = 1, as area is tiny 
# relative to scale of survey)
sh_na_idx <- which(is.na(raster::extract(r[[37]], shf_df)[,2]))
mwsh_na_idx <- which(is.na(raster::extract(r[[37]], mwsh_df)[,2]))
shf_df <- shf_df %>% 
  filter(!Benthoscape %in% c(1,4) & is.na(Benthoscape)==FALSE & 
           is.na(Backscatter)==FALSE & !(row_number() %in% sh_na_idx))
mwsh_df <- mwsh_df %>%
  filter(!Benthoscape %in% c(1,4) & is.na(Benthoscape)==FALSE & 
           is.na(Backscatter)==FALSE & !(row_number() %in% mwsh_na_idx))

# Log-transform all covariates except for benthoscape, backscatter, and SedMobFreq
covar_cols <- c("Bathymetry", "Backscatter", "BottomTempAtSurvey",
                "ShearVelocity", "SedMobFreq", "Benthoscape",
                "LogBathymetry", "LogBottomTempAtSurvey",
                "LogShearVelocity")
log_cols <- covar_cols[grep("Log", covar_cols)]
shf_df[log_cols] <- log(shf_df[, substr(log_cols, 4, 22), drop = TRUE])
mwsh_df[log_cols] <- log(mwsh_df[, substr(log_cols, 4, 22), drop = TRUE])

# Scale all covariates except for the benthoscape... scale MWSH covariates using
# SHF mean and sd
non_benthoscape_cols <- covar_cols[-which(covar_cols == "Benthoscape")]
scaled_cols <- paste0(non_benthoscape_cols, "Scaled")
shf_df[scaled_cols] <- scale(shf_df[ ,non_benthoscape_cols, drop = T])
mwsh_df[scaled_cols] <- mapply(
  x = mwsh_df[, non_benthoscape_cols, drop = TRUE], 
  y = shf_df[, non_benthoscape_cols, drop = TRUE],
  FUN = function(x, y) scale(x, center = mean(y), scale = sd(y))
)

# Benthoscape: merge silt with silty gravel with anemones, then make benthoscape 
# a factor
shf_df$Benthoscape <- ifelse(shf_df$Benthoscape == 7, 6, shf_df$Benthoscape)
mwsh_df$Benthoscape <- ifelse(mwsh_df$Benthoscape == 7, 6, mwsh_df$Benthoscape)
benthoscape_labels <- c("Gravelly sand", "Mixed sediment", 
                        "Sand", "Silt", "Tidal scoured mixed sediments")
shf_df$Benthoscape <- factor(shf_df$Benthoscape, levels = c(2:3, 5, 6, 8),
                            benthoscape_labels) %>% relevel(ref = "Mixed sediment")
mwsh_df$Benthoscape <- factor(mwsh_df$Benthoscape, levels = c(2:3, 5, 6, 8),
                              benthoscape_labels) %>% relevel(ref = "Mixed sediment")

############Simulated Individual Shell Height Measurements Dataset##############

# Reset seed
set.seed(123)

# Use commercial size SHF bins to simulate shell height values for individuals... 
# store these individual values in a new dataframe
nbins_commercial <- length(commercial_bins)
binocount <- matrix(0, nrow(shf_df), nbins_commercial) 
sh_df <- data.frame()
# Loop through each tow
for(i in 1:nrow(shf_df))
{
  # Loop through each commercial size frequency bin
  for(k in 1:nbins_commercial)
  {
    # Simulate individual shell height values
    val <- shf_df[i, commercial_bins[k], drop = T]/6
    binocount[i,k] <- floor(val) + rbinom(1, 1, val %% 1)
    new_values <- runif(binocount[i,k], min = 75+5*k, max=79.99+5*k)
    
    # Store each value as a row in a dataframe
    if(length(new_values)>0)
    {
      new_rows <- do.call(rbind, replicate(length(new_values),
                                           shf_df[i, -grep("BIN", names(shf_df))], 
                                           simplify = F)) %>% 
        mutate(HEIGHT = new_values, LogHEIGHT = log(new_values))
      sh_df <- rbind(sh_df, new_rows)
    }
  }
}

# Scale height, and then scale the MWSH data's height values using the same mean 
# and sd
height_col_names <- c("LogHEIGHT","HEIGHT")
sh_df[paste0(height_col_names, "Scaled")] <- scale(sh_df[ ,height_col_names, 
                                                           drop = T])
mwsh_df$LogHEIGHT <- log(mwsh_df$HEIGHT)
mwsh_df[paste0(height_col_names,"Scaled")] <- mapply(
  x = mwsh_df[, height_col_names, drop = TRUE], 
  y = sh_df[, height_col_names, drop = TRUE],
  FUN = function(x, y) scale(x, center = mean(y), scale = sd(y))
)

#################################Save Datasets##################################

# Examine dataframes
head(shf_df)
head(mwsh_df)
head(sh_df)

# Write processed MWSH, SHF, and SH data to rds file
data_list <- list(mwsh_df, shf_df, sh_df)
saveRDS(data_list, paste(getwd(),"/Data/Processed_MWSH_SH_Data.rds", sep = ""))