################################################################################
#
# Supplementary analysis involving a smaller spatial domain (and only July tows)
#
################################################################################

# Load necessary libraries
library(mgcv)
library(ggplot2)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(viridis)
source("Helper Functions.R")

# Load in the dataset restricted to a smaller spatial domain (and only July tows)
MWSH_dat <- read.csv("Data/July_MWSH.csv")
MWSH_dat$f_year <- factor(MWSH_dat$year)

# Also load in the larger survey datasets (to help with data processing etc.)
survey_data <- readRDS("Data/Processed_MWSH_SH_Data.rds")
mwsh_df <- survey_data[[1]] %>% convert_to_km()
shf_df <- survey_data[[2]] %>% convert_to_km()
sh_df <- survey_data[[3]] %>% convert_to_km()

# Plot the smaller spatial domain's raw July LWR data for all years, with a loess 
# smooth for each year
ggplot(data = MWSH_dat, aes(x = HEIGHT, y = WET_MEAT_WGT)) + 
  geom_point(alpha = 0.3) + facet_wrap(~year) + 
  geom_smooth(method = "loess", size = 1.2, se = F) + 
  xlab("Shell height (mm)") + ylab("Meat weight (g)") + 
  theme_classic() +
  theme(text = element_text(size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line.y = element_blank()) +
  coord_cartesian(xlim = c(80, 190), ylim = c(0, 100), expand = F) + 
  scale_x_continuous(breaks = c(80, 120, 160))
ggsave(paste0(getwd(),"/Figs/MWSH_yearly_data_supplement.jpeg"), plot=last_plot(), 
       width=6.5, height=4.5, units="in")

# Setup the shape file denoting the study area
raster_df <- readRDS("Data/Processed_Data_For_Prediction.rds") %>%
  convert_to_km()
hulls <- shf_df %>%
  summarise( geometry = st_combine(geometry)) %>% st_convex_hull()
raster_df <- st_intersection(x = raster_df, y = hulls)

# Grab land data for the spatial figure
land <- ne_countries(scale = "large", returnclass = "sf", 
                     continent = "North America") %>%
  st_crop(c(xmin = -75, ymin = 40, xmax = -50, ymax = 50)) %>% 
  st_transform(crs = 32620)
# Convert units from m to km
land_km <- st_geometry(land)/1000
st_crs(land_km) <- st_crs(land)
st_geometry(land) <- land_km

# Plot the tow locations from the smaller spatial domain (July only)
raster_df$dummy <- 1 # Fake variable just to plot the study area
xlimits <- c(220, 420); ylimits <- c(4910, 5050) # Limits for map axes
plot_predictions(raster_df, land, fill_var = "dummy", 
                 facet_vars = "year", xlimits = xlimits,
                 ylimits = ylimits, plot_title = "Tows for second sensitivity analysis", 
                 legend.pos = "none", font_size = 16,
                 facet_ncol = 4) + 
  geom_point(data = MWSH_dat, aes(x = x, y = y), col = "black", size = 0.7) +
  scale_fill_gradient(low = "darkgrey", high = "darkgrey") + guides(fill = "none")
ggsave(paste0(getwd(),"/Figs/study_area_data_ts_supplement.jpeg"), plot=last_plot(), 
       width=10, height=8, units="in")

# Create the scaled log shell height covariate for the smaller spatial domain dataset
MWSH_dat$LogHEIGHT <- log(MWSH_dat$HEIGHT)
MWSH_dat$LogHEIGHTScaled <- scale(
  MWSH_dat$LogHEIGHT,
  center = mean(sh_df$LogHEIGHT),
  scale = sd(sh_df$LogHEIGHT)
)

# Fit a GAM model with a time-varying shell height effect and random effects for tow
# to the smaller spatial domain dataset
MWSH_mgcv <- gam(WET_MEAT_WGT ~ 0 + f_year +  
    s(LogHEIGHTScaled, by = f_year, k = 3) +
    s(TOW_NO, bs = "re"),
  family = Gamma(link = "log"),
  data   = MWSH_dat,
  method = "REML")

# Estimate yearly conditional effect of shell height using this model
pred_dat <- expand.grid(year = c(2011:2019, 2021:2023), 
                        LogHEIGHTScaled = seq(-2.4,3, by = 0.1),
                        TOW_NO = 0)
pred_dat$f_year <- as.factor(pred_dat$year)
predictions <- predict(MWSH_mgcv, pred_dat, exclude = "s(TOW_NO)", se.fit = TRUE)
pred_dat$est <- predictions$fit 
pred_dat$est_se <- predictions$se.fit
# Generate plot of yearly conditional effect
pred_dat$HEIGHT <- exp(pred_dat$LogHEIGHTScaled*sd(sh_df$LogHEIGHT)+
                         mean(sh_df$LogHEIGHT))
cols <- rep(c("#DDBB00","#005BBB","black","firebrick2"), each = 3); shp_size <- 3.5
ggplot(pred_dat, aes(HEIGHT, exp(est),
                     ymin = exp(est - 1.96 * est_se),
                     ymax = exp(est + 1.96 * est_se),
                     col = f_year, fill = f_year, shape = f_year)) +
  geom_ribbon(col = NA, show.legend = F, alpha = 0.1) + 
  geom_line(aes(linetype = f_year), lwd = 1.3, alpha = 0.85) + 
  geom_point(data = pred_dat[c(169:180, 385:396, 553:564),], size = shp_size, 
             stroke = 1.5, alpha = 0.85)  +
  scale_alpha_identity() +
  coord_cartesian(expand = FALSE, ylim = c(0,75), xlim = c(80,170)) + 
  theme_classic() + scale_x_continuous(breaks = c(80, 110, 140, 170)) +
  theme(text = element_text(size = 16), legend.key.size = unit(2,"line")) +
  labs(x = "Shell height (mm)", y = "Predicted meat weight (g)",
       col = "Year", shape = "Year", fill = "Year", linetype = "Year") +  
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = rep(c(15,16,17,5), 3)) + 
  guides(shape = guide_legend(override.aes = list(size = shp_size))) + 
  scale_linetype_manual(values = rep(c("solid","dotted","dashed"), each = 4))
ggsave(paste0(getwd(),"/Figs/MWSH_SH_effect_supplement.jpeg"), plot=last_plot(), 
       width=6.5, height=6, units="in")