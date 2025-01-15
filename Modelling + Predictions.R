################################################################################
#
# Fit the models and examine predictions
#
################################################################################

# Load necessary libraries
library(sf)
library(tidyverse)
library(fmesher)
library(inlabru)
library(sdmTMB)
library(sdmTMBextra)
library(rnaturalearth)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
source("Helper Functions.R")

######################Model-Fitting + Covariate Effects#########################

# Load in processed survey data, and convert UTM units from m to km for stability 
# during modelling
survey_data <- readRDS("Data/Processed_MWSH_SH_Data.rds")
mwsh_df <- survey_data[[1]] %>% convert_to_km()
shf_df <- survey_data[[2]] %>% convert_to_km()
sh_df <- survey_data[[3]] %>% convert_to_km()

# Make non-sf versions of mwsh_df, shf_df, and sh_df as some functions like cor 
# can be finicky
mwsh_df_nonsf <- as.data.frame(mwsh_df)
shf_df_nonsf <- as.data.frame(shf_df)
sh_df_nonsf <- as.data.frame(sh_df)

# Examine the correlations between the covariates to ensure all pairwise |r| < 0.7
cor(mwsh_df_nonsf[c("BathymetryScaled", "LogBottomTempAtSurveyScaled", 
                "BackscatterScaled", "SedMobFreqScaled", 
                "LogShearVelocityScaled", "LogHEIGHTScaled")])
cor(shf_df_nonsf[c("BathymetryScaled", "LogBottomTempAtSurveyScaled", 
                  "BackscatterScaled", "SedMobFreqScaled", 
                  "LogShearVelocityScaled")])

# Grab land data for barrier and figures
land <- ne_countries(scale = "large", returnclass = "sf", 
                     continent = "North America") %>%
  st_crop(c(xmin = -75, ymin = 40, xmax = -50, ymax = 50)) %>% 
  st_transform(crs = 32620)
# Convert units from m to km
land_km <- st_geometry(land)/1000
st_crs(land_km) <- st_crs(land)
st_geometry(land) <- land_km

# Since MWSH tows are a subset of SHF tows, use SHF tows to create one mesh to 
# be used throughout... have to add barrier later
non_convex_bdry <- fm_nonconvex_hull(shf_df$geometry, 0.03, resolution = c(100, 50))
mesh <- fm_mesh_2d(boundary = non_convex_bdry, max.edge=c(5, 13.7), cutoff = 2.5, 
                   loc = shf_df$geometry)

# Setup mesh with sdmTMB and add barrier
mesh_mwsh <- sdmTMB::make_mesh(data = mwsh_df_nonsf, 
                               xy_cols = c("x","y"), mesh = mesh) %>%
  add_barrier_mesh(land, range_fraction = 0.1)
mesh_shf <- sdmTMB::make_mesh(data = shf_df_nonsf, 
                               xy_cols = c("x","y"), mesh = mesh) %>%
  add_barrier_mesh(land, range_fraction = 0.1)
mesh_sh <- sdmTMB::make_mesh(data = sh_df_nonsf, 
                             xy_cols = c("x","y"), mesh = mesh) %>%
  add_barrier_mesh(land, range_fraction = 0.1)

# Plot barrier mesh with data points... green triangles for barrier triangles
mesh_df_land <- mesh_mwsh$mesh_sf[mesh_mwsh$barrier_triangles, ]
ggplot() + geom_sf(data = land, fill = "grey44") + 
  gg(mesh_mwsh$mesh, edge.color = "grey4", linewidth = 2) + 
  geom_sf(data = shf_df, size = 0.65, col = "orange") + 
  geom_sf(data = mwsh_df, size = 0.65, col = "blue") + 
  geom_sf(data = mesh_df_land, col = "green", shape = 2) + 
  coord_sf(xlim = c(200, 420), ylim = c(4890, 5070)) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), panel.grid = element_blank())
ggsave(paste0(getwd(),"/Figs/mesh.jpeg"), plot=last_plot(), 
       width=8, height=6, units="in")

###########Abundance Model###########

# Fit the abundance model... linear term for temperature to avoid smooth too 
# sensitive to high temp. values
abundance_fit <- sdmTMB(
  nums_tow ~ 0 + f_year + 
    Benthoscape + 
    s(BathymetryScaled, k = 3) +
    s(BackscatterScaled, k = 3) +
    s(SedMobFreqScaled, k = 3) + 
    s(LogShearVelocityScaled, k = 3) + 
    LogBottomTempAtSurveyScaled,
  family = tweedie(link = "log"), 
  data = shf_df_nonsf, 
  mesh = mesh_shf, spatial = "on", 
  spatiotemporal = "rw", time = "year", 
  extra_time = c(2020),
  share_range = T,
  control = sdmTMBcontrol(newton_loops = 2)
) %>% run_extra_optimization(newton_loops = 2)
sanity(abundance_fit)
table_fit(abundance_fit) # Summary and CIs (minus smoothers)

# Plot randomized quantile residuals from the abundance model
qq_residual_plot(abundance_fit)
ggsave(paste0(getwd(),"/Figs/abundance_residuals.jpeg"), plot=last_plot(), 
       width=8, height=6, units="in")

# Plot conditional environmental effects from the abundance model
covar_names <- c("Benthoscape", "BathymetryScaled", 
                 "BackscatterScaled", "LogShearVelocityScaled", 
                 "SedMobFreqScaled", "LogBottomTempAtSurveyScaled")
covar_plot_labels <- c("Benthoscape", "Depth (m)", "Backscatter (intensity)", 
                       "Shear velocity (m/s)", "SMF (log scale time percent)", 
                       "Bottom temperature (°C)")
plot_effects(abundance_fit, covar_names,
             xlab = covar_plot_labels,
             ylim = lapply(vector("list", length = 6), function(x) c(-11,250)), 
             ylab = expression("Predicted 80 mm+ density (thousands/km"^2*")"), 
             nrow = 2, ncol = 3, no_ylab = c(2,3,5,6), yr = 2023)
ggsave(paste0(getwd(),"/Figs/abundance_enviro_effects.jpeg"), plot=last_plot(), 
       width=6.5, height=4.5, units="in")

############Shell Height Model###########

# Fit the shell height model... linear terms used to replace smooths with an sd 
# estimate = 0 in an initial fit
SH_fit <- sdmTMB(
  HEIGHT ~ 0 + f_year + 
    Benthoscape + 
    s(BathymetryScaled, k = 3) +
    BackscatterScaled +
    s(SedMobFreqScaled, k = 3) + 
    s(LogShearVelocityScaled, k = 3) + 
    s(LogBottomTempAtSurveyScaled, k = 3),
  family = Gamma(link = "log"), 
  data = sh_df_nonsf, 
  mesh = mesh_sh, spatial = "on", 
  spatiotemporal = "rw", time = "year", 
  extra_time = c(2020),
  share_range = T,
  control = sdmTMBcontrol(newton_loops = 2)
) %>% run_extra_optimization(newton_loops = 2)
sanity(SH_fit)
table_fit(SH_fit) # Summary and CIs (minus smoothers)

# Plot randomized quantile residuals from the shell height model
qq_residual_plot(SH_fit)
ggsave(paste0(getwd(),"/Figs/SH_residuals.jpeg"), plot=last_plot(), 
       width=8, height=6, units="in")

# Plot conditional environmental effects from the shell height model
plot_effects(SH_fit, covar_names,
             xlab = covar_plot_labels,
             ylim = lapply(vector("list", length = 6), function(x) c(90,142)), 
             ylab = "Predicted shell height (if 80 mm+)", nrow = 2, ncol = 3,
             no_ylab = c(2,3,5,6), yr = 2023)
ggsave(paste0(getwd(),"/Figs/SH_enviro_effects.jpeg"), plot=last_plot(), 
       width=6.5, height=4.5, units="in")

############Meat Weight Model###########

# Fit the meat weight model... linear terms used to replace smooths with an sd 
# estimate = 0 in an initial fit
MWSH_fit <- sdmTMB(
  WET_MEAT_WGT ~ 0 + f_year + 
    Benthoscape + 
    s(BathymetryScaled, k = 3) +
    s(BackscatterScaled, k = 3) +
    SedMobFreqScaled + 
    s(LogShearVelocityScaled, k = 3) + 
    LogBottomTempAtSurveyScaled + 
    s(LogHEIGHTScaled, by = f_year, k = 3),
  family = Gamma(link = "log"), 
  data = mwsh_df_nonsf, 
  mesh = mesh_mwsh, spatial = "on", 
  time = "year", spatiotemporal = "iid",
  extra_time = c(2020),
  share_range = T,
  control = sdmTMBcontrol(newton_loops = 2)
) %>% run_extra_optimization(newton_loops = 2)
sanity(MWSH_fit)
table_fit(MWSH_fit) # Summary and CIs (minus smoothers)

# Plot randomized quantile residuals from the meat weight model
qq_residual_plot(MWSH_fit)
ggsave(paste0(getwd(),"/Figs/MWSH_residuals.jpeg"), plot=last_plot(), 
       width=8, height=6, units="in")

# Plot conditional environmental effects from the meat weight model... assume a
# 100 mm scallop
plot_effects(MWSH_fit, covar_names,
             xlab = covar_plot_labels,
             ylim = lapply(vector("list", length = 6), function(x) c(10,21)), 
             ylab = "Predicted meat weight of 100 mm scallop (g)", nrow = 2, ncol = 3,
             no_ylab = c(2,3,5,6), yr = 2023)
ggsave(paste0(getwd(),"/Figs/MWSH_enviro_effects.jpeg"), plot=last_plot(), 
       width=6.5, height=4.5, units="in")

# Estimate yearly conditional effect of shell height using the meat weight model
pred_dat <- expand.grid(year = c(2011:2019, 2021:2023), 
                        Benthoscape = "Mixed sediment", 
                        BathymetryScaled = 0, BackscatterScaled = 0, 
                        SedMobFreqScaled = 0,
                        LogShearVelocityScaled = 0, LogBottomTempAtSurveyScaled = 0, 
                        LogHEIGHTScaled = seq(-2.4,3, by = 0.1))
pred_dat$f_year <- as.factor(pred_dat$year)
prediction <- predict(MWSH_fit, pred_dat, re_form = NA, se_fit = TRUE)
# Generate plot of yearly conditional effect
prediction$HEIGHT <- exp(prediction$LogHEIGHTScaled*sd(sh_df$LogHEIGHT)+
                           mean(sh_df$LogHEIGHT))
cols <- rep(scales::hue_pal()(6), each = 2); shp_size <- 4
ggplot(prediction, aes(HEIGHT, exp(est),
                       ymin = exp(est - 1.96 * est_se),
                       ymax = exp(est + 1.96 * est_se),
                       col = f_year, fill = f_year, shape = f_year)) +
  geom_ribbon(aes(alpha = ifelse(f_year == "2023", 0.18, 0.09)), 
              col = NA, show.legend = F) + 
  geom_point(data = prediction[c(169:180, 385:396, 553:564),], size = shp_size) + 
  geom_line(aes(alpha = ifelse(f_year == "2023", 1, 0.5)), 
            lwd = 1.2, show.legend = F) +
  scale_alpha_identity() +
  coord_cartesian(expand = FALSE, ylim = c(0,NA), xlim = c(80,170)) + 
  theme_classic() + scale_x_continuous(breaks = c(80, 110, 140, 170)) +
  theme(text = element_text(size = 16)) +
  labs(x = "Shell height (mm)", y = "Predicted meat weight (g)",
       col = "Year", shape = "Year", fill = "Year") +  
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = rep(c(16,17), 6)) + 
  guides(shape = guide_legend(override.aes = list(size = shp_size)))
ggsave(paste0(getwd(),"/Figs/MWSH_SH_effect.jpeg"), plot=last_plot(), 
       width=6.5, height=6, units="in")

# Sensitivity analysis: fit the meat weight model without spatial or spatio-temporal 
# random effects...
MWSH_fit_nospatial <- sdmTMB(
  WET_MEAT_WGT ~ 0 + f_year + 
    Benthoscape + 
    s(BathymetryScaled, k = 3) +
    BackscatterScaled +
    s(SedMobFreqScaled, k = 3) + 
    s(LogShearVelocityScaled, k = 3) + 
    LogBottomTempAtSurveyScaled + 
    s(LogHEIGHTScaled, by = f_year, k = 3),
  family = Gamma(link = "log"), 
  data = mwsh_df_nonsf, 
  mesh = mesh_mwsh, spatial = "off", 
  control = sdmTMBcontrol(newton_loops = 2)
) %>% run_extra_optimization(newton_loops = 2)
sanity(MWSH_fit_nospatial)
table_fit(MWSH_fit_nospatial) # Summary and CIs (minus smoothers)

#... and plot conditional environmental effects
plot_effects(MWSH_fit_nospatial, covar_names,
             xlab = covar_plot_labels,
             ylim = lapply(vector("list", length = 6), function(x) c(10,23)), 
             ylab = "Predicted meat weight of 100 mm scallop (g)", nrow = 2, ncol = 3,
             no_ylab = c(2,3,5,6), yr = 2023)
ggsave(paste0(getwd(),"/Figs/MWSH_enviro_effects_sensitivity.jpeg"), plot=last_plot(), 
       width=6.5, height=4.5, units="in")

###########Biomass Model###########

# Determine mid-points of the SHF bins
commercial_bins <- names(shf_df_nonsf[grep("BIN", names(shf_df_nonsf))[-c(1:16)]])
bin_midpoints <- as.numeric(substr(commercial_bins, 8,11)) + (4.99/2)

# Use meat weight model to predict meat weight in kg for each mid-point+tow combination
expanded_df <- shf_df_nonsf %>%
  slice(rep(1:n(), each = length(bin_midpoints))) %>%
  mutate(HEIGHT = rep(bin_midpoints, times = nrow(shf_df_nonsf)),
         LogHEIGHT = rep(log(bin_midpoints), times = nrow(shf_df_nonsf)))
expanded_df$LogHEIGHTScaled <-  pred_cov_scaling(train = sh_df$LogHEIGHT, 
                                                 pred = expanded_df$LogHEIGHT) # Scaled log
tow_mw_pred <- predict(MWSH_fit, newdata = expanded_df, type = "response")$est/1000

# Convert SHF bins to total kg/km^2 using predicted meat weights... apply q correction
for(i in 1:nrow(shf_df_nonsf))
{
  tow_weights <- tow_mw_pred[(length(commercial_bins)*(i-1)+1):(length(commercial_bins)*i)]
  shf_df_nonsf[i, "wt_tow"] <- (sum(tow_weights*shf_df_nonsf[i,commercial_bins])/
    ((800*5.334)/(1000*1000)))/0.36 # 1 standardized tow = 800m*5.334m, q = 0.36
}

# Fit the biomass model... linear terms used to replace smooths with an sd 
# estimate = 0 in an initial fit
biomass_fit <- sdmTMB(
  wt_tow ~ 0 + f_year + 
    Benthoscape + 
    s(BathymetryScaled, k = 3) +
    s(BackscatterScaled, k = 3) +
    s(SedMobFreqScaled, k = 3) + 
    s(LogShearVelocityScaled, k = 3) + 
    LogBottomTempAtSurveyScaled,
  family = tweedie(link = "log"), 
  data = shf_df_nonsf, 
  mesh = mesh_shf, spatial = "on", 
  spatiotemporal = "rw", time = "year", 
  extra_time = c(2020),
  share_range = T,
  control = sdmTMBcontrol(newton_loops = 2)
) %>% run_extra_optimization(newton_loops = 2)
sanity(biomass_fit)
table_fit(biomass_fit) # Summary and CIs (minus smoothers)

# Plot randomized quantile residuals from the biomass model
qq_residual_plot(biomass_fit)
ggsave(paste0(getwd(),"/Figs/biomass_residuals.jpeg"), plot=last_plot(), 
       width=8, height=6, units="in")

# Plot conditional environmental effects from the biomass model
plot_effects(biomass_fit, covar_names,
             xlab = covar_plot_labels,
             ylim = lapply(vector("list", length = 6), function(x) c(-152,7222)), 
             ylab = expression("Predicted 80 mm+ biomass density (kg/km"^2*")"), 
             nrow = 2, ncol = 3, no_ylab = c(2,3,5,6), yr = 2023)
ggsave(paste0(getwd(),"/Figs/biomass_enviro_effects.jpeg"), plot=last_plot(), 
       width=6.5, height=4.5, units="in")

########################Spatial Predictions + Indices###########################

# Configure a few settings for the plots
xlimits <- c(220, 420); ylimits <- c(4910, 5050) # Limits for map axes
text_size <- 16 # For all figures below
map_width <- 10; map_height <- 8 # Dimensions of saved map figures
index_width <- 7; index_height <- 6 # Dimensions of saved index figures 
                                    # (including scatterplots)
nsims <- 500 # How many times to sample from joint precision matrices for SE maps

# Load in the processed raster data and then convert the UTM coordinates from m 
# to km to match the units used in modelling
raster_df <- readRDS("Data/Processed_Data_For_Prediction.rds") %>%
  convert_to_km()

# Clip the raster data using a convex hull generated from the SHF data locations
hulls <- shf_df %>%
  summarise( geometry = st_combine(geometry)) %>% st_convex_hull()
raster_df <- st_intersection(x = raster_df, y = hulls)

# Scale raster data where necessary, using the SHF data mean and sd for consistency
covariates <- c("Bathymetry", "LogBottomTempAtSurvey", "Backscatter",
                "SedMobFreq", "LogShearVelocity")
raster_df[paste0(covariates,"Scaled")] <- 
  pred_cov_scaling(train = shf_df_nonsf[covariates], 
                   pred = raster_df[,covariates, drop = T])

# Calculate the area of each cell for biomass and abundance indices
index_area <- (raster_df$y[2]-raster_df$y[1])^2

# To aid interpretation of the predictions, plot the data locations over the study 
# area (SHF data is in orange and MWSH data is in blue... MWSH locations are a 
# subset of the SHF locations)
raster_df$dummy <- 1 # Fake variable just to plot the study area
plot_predictions(raster_df, land, fill_var = "dummy", 
                 facet_vars = "year", xlimits = xlimits,
                 ylimits = ylimits, plot_title = "Tow locations", 
                 legend.pos = "none", font_size = text_size,
                 facet_ncol = 3) + 
  geom_point(data = shf_df, aes(x = x, y = y), col = "orange", size = 0.7) + 
  geom_point(data = mwsh_df, aes(x = x, y = y), col = "blue", size = 0.7)
ggsave(paste0(getwd(),"/Figs/study_area_data_ts.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Spatially plot summer bottom temps. by year to aid interpretation of predictions
plot_predictions(raster_df, land, fill_var = "BottomTempAtSurvey", 
                 facet_vars = "year", xlimits = xlimits,
                 ylimits = ylimits, plot_title = "Average summer bottom temperature (°C)", 
                 font_size = text_size)
ggsave(paste0(getwd(),"/Figs/temp_timeseries.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units = "in")

# Generate a plot and correlation to check that overall commercial size abundance is 
# a good proxy for 100 mm abundance (as assumed by the meat weight index)
correlation <- round(cor(shf_df$nums_tow, shf_df$BIN_ID_100, method = "spearman"),2)
ggplot(data = shf_df, aes(x = sqrt(nums_tow*0.36*1000/(5.334*800)), y = sqrt(BIN_ID_100))) + 
  geom_point(size = 0.75, alpha = 0.75) +
  xlab("Square root of 80 mm+ abundance per tow") + 
  ylab("Square root of 100-104 mm abundance per tow") + 
  theme(text = element_text(size = text_size)) 
ggsave(paste0(getwd(),"/Figs/corr_abundance_100_104.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")

###########Abundance Model###########

# Make predictions with standard errors using the abundance model
abundance_pred_df <- predict(abundance_fit, newdata = raster_df[drop = T], type = "response")
abundance_pred_df$se <- predict(abundance_fit, newdata = raster_df[drop = T], nsim = nsims, 
                             type = "response") %>% apply(1, sd)

# Calculate abundance index (in millions)
index_abundance <- predict(abundance_fit, newdata = raster_df[drop = T], 
                           return_tmb_object = T) %>% 
  get_index(bias_correct = T, area = index_area/1000)

# Plot the abundance model's predictions
plot_predictions(abundance_pred_df, land, fill_var = "est", facet_vars = "year", 
                 xlimits = xlimits, ylimits = ylimits, font_size = text_size,
                 plot_title = expression("Predicted 80 mm+ density (thousands/km"^2*")")) + 
  scale_fill_gradientn("", colors = viridis(100), 
                       trans = trans_new("log10", transform = sqrt, 
                                         inverse = function(x) x^2), 
                       breaks = c(0, 75, 350, 750), limits = c(0, 750))
ggsave(paste0(getwd(),"/Figs/abundance_pred.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the standard errors associated with the above predictions
plot_predictions(abundance_pred_df, land, fill_var = "se", xlimits = xlimits,
                 ylimits = ylimits, facet_vars = "year", font_size = text_size,
                 plot_title = expression("Predicted 80 mm+ density (thousands/km"^2*") - std. errors")) +
  scale_fill_gradientn("", colors = viridis(100), 
                       trans = trans_new("log10", transform = sqrt, 
                                         inverse = function(x) x^2), 
                       breaks = c(0, 50, 200, 700), limits = c(0, 700))
ggsave(paste0(getwd(),"/Figs/abundance_pred_std_errors.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the abundance model's random fields (spatial + spatio-temporal)
plot_predictions(abundance_pred_df, land, fill_var = "est_rf", 
                 facet_vars = "year", xlimits = xlimits, font_size = text_size,
                 ylimits = ylimits, 
                 plot_title = "Abundance model random effects")
ggsave(paste0(getwd(),"/Figs/abundance_RFs.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the abundance index
plot_index(df = index_abundance, year_col = "year", est_col = "est", lwr_ci_col = "lwr", 
           upr_ci_col = "upr", ylabel = "Total predicted 80 mm+\nabundance (millions)")
ggsave(paste0(getwd(),"/Figs/abundance_index.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")

############Shell Height Model###########

# Make predictions with standard errors using the shell height model
SH_pred_df <- predict(SH_fit, newdata = raster_df[drop = T], type = "response")
SH_pred_df$se <- predict(SH_fit, newdata = raster_df[drop = T], nsim = nsims, 
                         type = "response") %>% apply(1, sd)

# Calculate the shell height index, weighted by abundance
density_sum <- abundance_pred_df %>% group_by(year) %>% summarise(total = sum(est))
index_SH <- predict(SH_fit, newdata = raster_df[drop = T], return_tmb_object = T) %>% 
  get_index(bias_correct = T, area = abundance_pred_df$est) %>%
  left_join(density_sum) %>%
  mutate(est = est/total, lwr = lwr/total, upr = upr/total)

# Plot the shell height model's predictions
plot_predictions(SH_pred_df, land, fill_var = "est", facet_vars = "year", 
                 xlimits = xlimits, ylimits = ylimits, font_size = text_size,
                 plot_title = "Predicted shell height (if 80 mm+)")
ggsave(paste0(getwd(),"/Figs/SH_pred.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the standard errors associated with the above predictions
plot_predictions(SH_pred_df, land, fill_var = "se", xlimits = xlimits,
                 ylimits = ylimits, facet_vars = "year", font_size = text_size,
                 plot_title = "Predicted shell height (if 80 mm+) - std. errors",
                 scale_limits = c(0, NA))
ggsave(paste0(getwd(),"/Figs/SH_pred_std_errors.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the shell height model's random fields (spatial + spatio-temporal)
plot_predictions(SH_pred_df, land, fill_var = "est_rf", 
                 facet_vars = "year", xlimits = xlimits, font_size = text_size,
                 ylimits = ylimits, 
                 plot_title = "Shell height model random effects")
ggsave(paste0(getwd(),"/Figs/SH_RFs.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the shell height index
plot_index(df = index_SH, year_col = "year", est_col = "est", lwr_ci_col = "lwr", 
           upr_ci_col = "upr", ylabel = "Mean predicted shell\nheight (if 80 mm+)", 
           y_lowerlim = 105)
ggsave(paste0(getwd(),"/Figs/SH_index.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")

############Meat Weight Model###########

# Make predictions with standard errors using the meat weight model... assume a 
# 100 mm scallop
raster_df$LogHEIGHTScaled <- (log(100)-mean(sh_df$LogHEIGHT))/sd(sh_df$LogHEIGHT)
MWSH_pred_df <- predict(MWSH_fit, newdata = raster_df[drop = T], 
                        type = "response")
MWSH_pred_df$se <- predict(MWSH_fit, newdata = raster_df[drop = T], nsim = nsims, 
                           type = "response") %>% apply(1, sd)

# Calculate meat weight index (based on 100 mm scallop), weighted by abundance
index_MWSH <- predict(MWSH_fit, newdata = raster_df[drop = T], return_tmb_object = T) %>% 
  get_index(bias_correct = T, area = abundance_pred_df$est) %>%
  left_join(density_sum) %>%
  mutate(est = est/total, lwr = lwr/total, upr = upr/total)

# Plot the meat weight model's predictions for a 100 mm scallop
plot_predictions(MWSH_pred_df, land, fill_var = "est", 
                 facet_vars = "year", xlimits = xlimits, font_size = text_size,
                 ylimits = ylimits, 
                 plot_title = "Predicted meat weight of 100 mm scallop (g)")
ggsave(paste0(getwd(),"/Figs/MWSH_pred.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the standard errors associated with the above predictions
plot_predictions(MWSH_pred_df, land, fill_var = "se", xlimits = xlimits,
                 facet_vars = "year", ylimits = ylimits, font_size = text_size,
                 plot_title = "Predicted meat weight of 100 mm scallop (g) - std. errors",
                 scale_limits = c(0, NA))
ggsave(paste0(getwd(),"/Figs/MWSH_pred_std_errors.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the meat weight model's random fields (spatial + spatio-temporal)
plot_predictions(MWSH_pred_df, land, fill_var = "est_rf", 
                 facet_vars = "year", xlimits = xlimits, font_size = text_size,
                 ylimits = ylimits, 
                 plot_title = "Meat weight model random effects")
ggsave(paste0(getwd(),"/Figs/MWSH_RFs.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the meat weight index (based on 100 mm scallop)
plot_index(df = index_MWSH, year_col = "year", est_col = "est", lwr_ci_col = "lwr", 
           upr_ci_col = "upr", ylabel = "Mean predicted meat weight\nof 100 mm scallop (g)", 
           y_lowerlim = 9)
ggsave(paste0(getwd(),"/Figs/MWSH_index.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")

###########Biomass Model###########

# Make predictions with standard errors using the biomass model
biomass_pred_df <- predict(biomass_fit, newdata = raster_df[drop = T], type = "response")
biomass_pred_df$se <- predict(biomass_fit, newdata = raster_df[drop = T], nsim = nsims, 
                           type = "response") %>% apply(1, sd)

# Calculate biomass index (in metric tonnes)
index_biomass <- predict(biomass_fit, newdata = raster_df[drop = T], 
                         return_tmb_object = T) %>% 
  get_index(bias_correct = T, area = index_area/1000)

# Plot the biomass model's predictions
plot_predictions(biomass_pred_df, land, fill_var = "est", facet_vars = "year", 
                 xlimits = xlimits, ylimits = ylimits, font_size = text_size,
                 plot_title = expression("Predicted 80 mm+ biomass density (kg/km"^2*")")) + 
  scale_fill_gradientn("", colors = viridis(100), 
                       trans = trans_new("log10", transform = sqrt, 
                                         inverse = function(x) x^2), 
                       breaks = c(0, 2000, 8000, 16000), limits = c(0, 16000))
ggsave(paste0(getwd(),"/Figs/biomass_pred.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the standard errors associated with the above predictions
plot_predictions(biomass_pred_df, land, fill_var = "se", xlimits = xlimits,
                 ylimits = ylimits, facet_vars = "year", font_size = text_size,
                 plot_title = expression("Predicted 80 mm+ biomass density (kg/km"^2*") - std. errors")) +
  scale_fill_gradientn("", colors = viridis(100), 
                       trans = trans_new("log10", transform = sqrt, 
                                         inverse = function(x) x^2), 
                       breaks = c(0, 1000, 5000, 12000), limits = c(0, 12000))
ggsave(paste0(getwd(),"/Figs/biomass_pred_std_errors.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the biomass model's random fields (spatial + spatio-temporal)
plot_predictions(abundance_pred_df, land, fill_var = "est_rf", 
                 facet_vars = "year", xlimits = xlimits, font_size = text_size,
                 ylimits = ylimits, 
                 plot_title = "Biomass model random effects")
ggsave(paste0(getwd(),"/Figs/biomass_RFs.jpeg"), plot=last_plot(), 
       width=map_width, height=map_height, units="in")

# Plot the biomass index
plot_index(df = index_biomass, year_col = "year", est_col = "est", lwr_ci_col = "lwr", 
           upr_ci_col = "upr", ylabel = "Total predicted 80 mm+\nbiomass (metric tonnes)")
ggsave(paste0(getwd(),"/Figs/biomass_index.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")

# Plot the biomass and survey indices as a scatterplot to highlight 2023
biomass_abundance_indices <- merge(index_biomass, index_abundance, by = "year", 
                                   suffixes = c(".biomass",".abundance"))
biomass_abundance_indices_no2023 <- biomass_abundance_indices %>% filter(year < 2023)
ggplot(data = biomass_abundance_indices, aes(x = est.abundance, y = est.biomass, 
                                             label = paste0("'", substr(year,3,4)))) +
  geom_smooth(data = biomass_abundance_indices_no2023, method = "lm", se = F, 
              size = 2, col = "cornflowerblue") + 
  geom_point(size = 2) + 
  geom_linerange(aes(ymin = lwr.biomass, ymax = upr.biomass), width = 2, alpha = 0.45) +
  geom_linerange(aes(xmin = lwr.abundance, xmax = upr.abundance), height = 2, alpha = 0.45) + 
  geom_text(size = 4.5, nudge_x = 16, nudge_y = 500) +
  ylab("Total predicted 80 mm+ biomass (metric tonnes)") + 
  xlab("Total predicted 80 mm+ abundance (millions)") + theme_classic() +
  theme(text = element_text(size = text_size))
ggsave(paste0(getwd(),"/Figs/abundance_biomass_scatterplot.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")

######Comparing Biomass Model To Model Currently Used For Stock Assessment######

# Load output from the stock assessment model currently in use (separate output 
# for each management area)
Spa1A <- read.csv("Data/Spa1AModelOutput.csv")
Spa1B <- read.csv("Data/Spa1BModelOutput.csv")
Spa4 <- read.csv("Data/Spa4ModelOutput.csv")

# Combine the management areas into one dataframe
DF <- bind_rows(Spa1A %>% mutate(Dataset = 'Spa1A'),
                Spa1B %>% mutate(Dataset = 'Spa1B'),
                Spa4 %>% mutate(Dataset = 'Spa4'))

# Combine the biomass model index with an index calculated from the stock assessment
# output
plot_data <- data.frame(biomass = DF$X50.[grep("B", DF$X)], 
                        year = c(1997:2023, 1997:2023, 1983:2023)) %>%
  group_by(year) %>% filter(year>2010) %>% summarize(est = sum(biomass)) %>%
  bind_rows(index_biomass[c("year","est")]) %>% 
  mutate(model = c(rep("Current stock\nassessment\nmodel", 13), rep("Biomass\nmodel", 12))) %>%
  bind_rows(data.frame(year = 2020, est = NA, model = "Biomass\nmodel"))

# Plot the indices together
ggplot(plot_data, aes(x = year, y = est, col = model)) +
  geom_point(size = 2) + geom_line(lwd = 1) +
  scale_y_continuous(limits = c(0, 16600), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(2011, 2023, by = 2)) +
  theme_classic() + theme(text = element_text(size = text_size),
                          legend.text = element_text(margin = margin(b = 10))) +
  labs(y = "Total predicted 80 mm+ biomass (metric tonnes)", x = "Year", col = "")
ggsave(paste0(getwd(),"/Figs/biomass_SAM_indices_compare.jpeg"), plot=last_plot(), 
       width=index_width, height=index_height, units="in")