################################################################################
#
# Functions to help with fitting models and examining model output
#
################################################################################

# Takes a UTM sf object in m and converts the units to km
convert_to_km <- function(dat) {
  # Do conversion
  dat_km <- st_geometry(dat)/1000
  st_crs(dat_km) <- st_crs(dat)
  st_geometry(dat) <- dat_km
  
  # Save coordinates in a column
  dat[c("x", "y")] <- st_coordinates(dat)
  
  return(dat)
}

# Can take an sdmTMB model object and output a Q-Q plot of the randomized quantile 
# residuals
qq_residual_plot <- function(model, font_size = 16, ref_line_col = "red", 
                             ref_line_width = 1.2) {
  # Compute residuals (with arbitrary seed for reproducibility)
  set.seed(1)
  qq_data <- data.frame(Residuals = residuals(model))
  
  # Create Q-Q plot using ggplot
  ggplot(qq_data, aes(sample = Residuals)) +
    geom_qq() +
    geom_abline(intercept = 0, slope = 1, col = ref_line_col, linetype = "dashed", 
                linewidth = ref_line_width) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") +
    ggtitle("Q-Q plot of randomized quantile residuals") + 
    theme(text = element_text(size = font_size))
}

# For a given sdmTMB model fit, creates a table of parameter estimates and std. 
# errors, plus CIs... ignores smoothers
table_fit <- function(fit, rounding = 3, conf_level = 0.95) {
  fixed <- tidy(fit, "fixed", conf.int = T, conf.level = conf_level)[-1]
  ran_pars <- tidy(fit, "ran_pars", conf.int = T, conf.level = conf_level)[-1]
  parameters <- c(tidy(fit, "fixed")$term, tidy(fit, "ran_pars")$term)
  result <- data.frame(parameter = parameters, 
                       rbind(round(fixed, rounding), round(ran_pars, rounding)))
  return(result)
}

# Generates environmental covariate conditional effects plots with 95% confidence 
# bounds for a model fit... fit is the model object, xvars is a vector containing 
# the names of covariates to plot effects for, xlabs is a vector of labels 
# corresponding to those covariates, ylimits is a list containing a ylim vector 
# for each covariate, and ylab is a label for all y-axes to share... nrow and ncol 
# define how the plot facets are arranged... no_ylab is a vector of + integers 
# (e.g., c(1,2,3)) corresponding to the covariates to drop y-axes labels for... 
# yr is the year to set f_year/year to... exp = T will convert log-transformed 
# covariates back to their natural scale... function will always reverse mean-centering 
# and scaling
plot_effects <- function(fit, xvars, yr = 2023, xlabs, ylab, no_ylab = NULL, ylimits,
                         nrow = NULL, ncol = NULL,
                         font_size = 10, linewidth = 1, ptsize = 1.6, 
                         color_est = "deepskyblue", color_se = "grey",
                         alpha_est = 1, alpha_se = 0.5, exp = T)
{
  # Create a list to store the plots
  plot_list <- list()
  
  # Generate plots and store them in plot_list
  for(i in 1:length(xvars))
  {
    # Set up initial dataframe to expand... 100 mm scallop assumed for meat weight 
    # model
    pred_dat <- data.frame(year = yr, f_year = factor(yr),
                           Benthoscape = factor("Mixed sediment"), 
                           BathymetryScaled = 0, 
                           BackscatterScaled = 0, 
                           LogShearVelocityScaled = 0, 
                           LogBottomTempAtSurveyScaled = 0,
                           SedMobFreqScaled = 0)
    pred_dat$LogHEIGHTScaled <- (log(100)-mean(sh_df$LogHEIGHT))/sd(sh_df$LogHEIGHT)
    
    # Expand dataframe using covariate of interest
    if(xvars[i] != "Benthoscape")
    {
      xvar_data <- seq(min(fit$data[[xvars[i]]], na.rm = TRUE), 
                       max(fit$data[[xvars[i]]], na.rm = TRUE), by = 0.1)
      pred_dat <- pred_dat[rep(1, length(xvar_data)), ]
      pred_dat[[xvars[i]]] <- xvar_data
      pred_dat$unscaled_covar <- pred_dat[[xvars[i]]]*sd(shf_df[[gsub("Scaled", "", xvars[i])]]) + 
        mean(shf_df[[gsub("Scaled", "", xvars[i])]]) # Reverse mean-centering and scaling
      if(grepl("Log", xvars[i]) && exp == T) # Reverse log if required
      {
        pred_dat$unscaled_covar <- exp(pred_dat$unscaled_covar)
      }
    }else{
      xvar_data <- levels(fit$data$Benthoscape)
      pred_dat <- pred_dat[rep(1, length(xvar_data)), ]
      pred_dat[[xvars[i]]] <- factor(xvar_data)
    }
    
    # Make prediction
    prediction <- predict(fit, pred_dat, re_form = NA, se_fit = TRUE)
    
    # Generate plot
    if(xvars[i] != "Benthoscape")
    {
      rug_col_name <- if (grepl("Log", xvars[i]) && exp == TRUE) { # Reverse log if required
        gsub("Log|Scaled", "", xvars[i])
      } else {
        gsub("Scaled", "", xvars[i])
      }
      plot_list[[i]] <- ggplot(prediction, aes(x = unscaled_covar, exp(est),
                                               ymin = exp(est - 1.96 * est_se),
                                               ymax = exp(est + 1.96 * est_se))) +
        geom_rug(data = fit$data, inherit.aes = F, 
                 aes(x = .data[[rug_col_name]]), sides = "b", alpha = 0.5) +
        geom_ribbon(fill = color_se, col = NA, alpha = alpha_se) + 
        geom_line(lwd = linewidth, col = color_est, alpha = alpha_est) +
        coord_cartesian(expand = FALSE, ylim = ylimits[[i]]) + 
        theme_classic() +
        theme(text = element_text(size = font_size), 
              axis.text = element_text(size = font_size)) +
        labs(x = xlabs[i], y = ylab)
    }else{
      levels(prediction$Benthoscape) <- c("Gr","Mi","Sa","Si","Ti")
      plot_list[[i]] <- ggplot(prediction, aes(x = .data[[xvars[i]]], exp(est),
                                               ymin = exp(est - 1.96 * est_se),
                                               ymax = exp(est + 1.96 * est_se))) +
        geom_errorbar(col = color_se, lwd = linewidth*1.2, alpha = alpha_se) + 
        geom_point(size = ptsize, col = color_est, alpha = alpha_est) +
        coord_cartesian(expand = FALSE, ylim = ylimits[[i]]) + 
        theme_classic() +
        theme(text = element_text(size = font_size), 
              axis.text = element_text(size = font_size)) +
        labs(x = xlabs[i], y = ylab)
    }
    # If plot doesn't have a ylabel, remove
    if(i %in% c(no_ylab))
    {
      plot_list[[i]] <- plot_list[[i]] + theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())
    }
  }
  
  # Return the final combined plot
  wrap_plots(plot_list, ncol = ncol, nrow = nrow) + plot_layout(axis_titles = "collect")
}

# Takes one or more columns from a dataframe with data used in prediction, and 
# scales them with the associated mean and sd of the training data
pred_cov_scaling <- function(train, pred) {
  if (is.data.frame(train) && is.data.frame(pred)) {
    Map(x = pred, y = train, function(x, y) scale(x, center = mean(y), 
                                                  scale = sd(y)))
  } else {
    scale(pred, center = mean(train), scale = sd(train))
  }
}

# Plots a map of model output... pred_df is the df containing the model predictions, 
# landmass is an sf object containing the surrounding landmass, and x_var, y_var, and 
# fill_var are the names of the pred_df columns that contain the x coords, y coords, 
# and model output respectively. facet_vars is a vector of column names from pred_df 
# to facet over (such as a year column)
plot_predictions <- function(pred_df, landmass, x_var = "x", y_var = "y", 
                             fill_var = "est", font_size = 16, xlimits = NULL, 
                             ylimits = NULL, facet_vars = NULL, facet_ncol = 4,
                             facet_nrow = 4, plot_title = NULL, 
                             legend.pos = "right", scale_limits = NULL) {
  
  # Determine if facet_wrap or facet_grid should be used
  facet_type <- if (length(facet_vars) == 1) "wrap" else "grid"
  
  # Create the base ggplot object
  p <- ggplot() + 
    geom_raster(data = pred_df, aes(x = .data[[x_var]], y = .data[[y_var]], 
                                    fill = .data[[fill_var]])) +
    scale_fill_gradientn(name = "", colours = viridis(100),
                         limits = scale_limits) + 
    geom_sf(data = landmass, fill = "grey22") +
    theme(text = element_text(size = font_size), 
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = legend.pos) + 
    coord_sf(xlim = xlimits, ylim = ylimits) +
    ggtitle(plot_title)
  
  # Add facet based on determined facet_type
  if (facet_type == "wrap") {
    p <- p + facet_wrap(facet_vars, ncol = facet_ncol, nrow = facet_nrow)
  } else {
    p <- p + facet_grid(facet_vars)
  }
  return(p)
}

# Takes a dataframe df with a column for year (year_col) and columns for the index 
# estimate and lower and upper CI bounds (ie., est_col, lwr_ci_col, and upr_ci_col)
# and plots the time series... years_NA is a vector of years for which the index was not 
# estimated, which the function will use to create discontinuities in the plot
plot_index <- function(df, year_col, est_col, lwr_ci_col, upr_ci_col, 
                       years_NA = 2020, text_size = 16, xlabel = "Year", ylabel, 
                       y_lowerlim = 0, y_upperlim = NA, 
                       est_color = "black", ci_color = "grey", 
                       line_size = 1, point_size = 2)
{
  # Insert missing year NAs so proper discontinuities will appear in plot
  if (!is.null(years_NA)) {
    na_rows <- data.frame(
      year_col = years_NA,
      est_col = NA,
      lwr_ci_col = NA,
      upr_ci_col = NA
    )
    names(na_rows) <- c(year_col, est_col, lwr_ci_col, upr_ci_col)
    df <- rbind(df[c(year_col, est_col, lwr_ci_col, upr_ci_col)], na_rows)
  }
  
  # Generate the plot
  ggplot(df, aes(x = .data[[year_col]], y = .data[[est_col]])) +
    geom_ribbon(aes(ymin = .data[[lwr_ci_col]], ymax = .data[[upr_ci_col]]), 
                fill = ci_color) +
    geom_line(lwd = line_size, colour = est_color) +
    geom_point(size = point_size, colour = est_color) +
    labs(x = xlabel, y = ylabel) + theme_classic() +
    scale_y_continuous(limits = c(y_lowerlim, y_upperlim), expand = c(0,0)) +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size)) +
    scale_x_continuous(breaks = seq(min(df[year_col]), max(df[year_col]),
                                    by = floor((max(df[year_col])-1-min(df[year_col]))/4)))
}