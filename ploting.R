#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library(latex2exp)

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

load("improper1_noInt_fitted.RData")
load("improper1_typeI_fitted.RData")
load("improper1_typeII_fitted.RData")
load("improper1_typeIII_fitted.RData")
load("improper1_typeIV_fitted.RData")

################################################################################
#Load in simulation data (done stupidly, change!)
load("./Simulated_data/sc1_data_2.RData")
lambda_sc1.df <- lambda.df
lambda_sc1.df$space.time = 1:nrow(lambda_sc1.df)

load("./Simulated_data/sc2_data_2.RData")
lambda_sc2.df <- lambda.df
lambda_sc2.df$space.time = 1:nrow(lambda_sc2.df)

load("./Simulated_data/sc3_data_2.RData")
lambda_sc3.df <- lambda.df
lambda_sc3.df$space.time = 1:nrow(lambda_sc3.df)

load("./Simulated_data/sc4_data_2.RData")
lambda_sc4.df <- lambda.df
lambda_sc4.df$space.time = 1:nrow(lambda_sc4.df)

load("./Simulated_data/sc5_data_2.RData")
lambda_sc5.df <- lambda.df
lambda_sc5.df$space.time = 1:nrow(lambda_sc5.df)

load("./Simulated_data/sc6_data_2.RData")
lambda_sc6.df <- lambda.df
lambda_sc6.df$space.time = 1:nrow(lambda_sc6.df)

load("./Simulated_data/sc7_data_2.RData")
lambda_sc7.df <- lambda.df
lambda_sc7.df$space.time = 1:nrow(lambda_sc7.df)

load("./Simulated_data/sc8_data_2.RData")
lambda_sc8.df <- lambda.df
lambda_sc8.df$space.time = 1:nrow(lambda_sc8.df)

load("./Simulated_data/sc9_data_2.RData")
lambda_sc9.df <- lambda.df
lambda_sc9.df$space.time = 1:nrow(lambda_sc9.df)

load("./Simulated_data/sc10_data_2.RData")
lambda_sc10.df <- lambda.df
lambda_sc10.df$space.time = 1:nrow(lambda_sc10.df)

load("./Simulated_data/sc11_data_2.RData")
lambda_sc11.df <- lambda.df
lambda_sc11.df$space.time = 1:nrow(lambda_sc11.df)

load("./Simulated_data/sc12_data_2.RData")
lambda_sc12.df <- lambda.df
lambda_sc12.df$space.time = 1:nrow(lambda_sc12.df)

################################################################################
# Plot the maps themselves
plt_first_level <- ggplot(data = first_level_admin_map) + 
  geom_sf(aes(), 
          alpha = 1,
          color="black") + ggtitle("First-level administrative areas") + 
  theme(plot.title = element_text(size = 15,  hjust = 0.5),
        axis.title.x = element_blank(), #Remove axis and background grid
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin =  unit(c(0, 0, 0, 0), "inches"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.spacing = unit(1, 'lines')) +
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right"))

plt_second_level <- ggplot(data = second_level_admin_map) + 
  geom_sf(aes(), 
          alpha = 1,
          color="black") + ggtitle("Second-level administrative areas") + 
  theme(plot.title = element_text(size = 15,  hjust = 0.5),
        axis.title.x = element_blank(), #Remove axis and background grid
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin =  unit(c(0, 0, 0, 0), "inches"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.spacing = unit(1, 'lines')) +
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right"))

## Save as 10 by 6, name: germany_maps
ggarrange(plt_first_level, NULL, plt_second_level,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1))

################################################################################
# Plot two continuous risk surfaces at similar times w. different amount of knots
load("./Data/Simulated_risk_surfaces/sc3_risk_surface_2.RData")
risk_surface.list_sc3 = risk_surface.list

load("./Data/Simulated_risk_surfaces/sc9_risk_surface_2.RData")
risk_surface.list_sc9 = risk_surface.list

plt_sc3 <- heatmap_points(risk_surface.list_sc3,
                          polygon_grid2,
                          t_axis[1],
                          title = "20 X 20 Knots")


plt_sc9 <- heatmap_points(risk_surface.list_sc9,
                          polygon_grid2,
                          t_axis[1],
                          title = "10 X 10 Knots")

# Save to pdf as 10 by 6, name: continuous_risk_20_vs_10
ggarrange(plt_sc3, NULL, plt_sc9,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)

################################################################################
# Plot the temporal trends

# Save to pdf as 12 by 9: const_temporal_trend_fitted
plot_fitted_temporal_trends(model_on_first_level_20_knots = improper_typeI_sc1,
                            model_on_first_level_10_knots = improper_typeI_sc7,
                            model_on_second_level_20_knots = improper_typeI_sc2,
                            model_on_second_level_10_knots = improper_typeI_sc8,
                            tT)

# Save to pdf as 12 by 9: lin_increasing_temporal_trend_fitted
plot_fitted_temporal_trends(model_on_first_level_20_knots = improper_typeI_sc3,
                            model_on_first_level_10_knots = improper_typeI_sc9,
                            model_on_second_level_20_knots = improper_typeI_sc4,
                            model_on_second_level_10_knots = improper_typeI_sc10,
                            tT)

# Save to pdf as 12 by 9: change_point_temporal_trend_fitted
plot_fitted_temporal_trends(model_on_first_level_20_knots = improper_typeI_sc5,
                            model_on_first_level_10_knots = improper_typeI_sc11,
                            model_on_second_level_20_knots = improper_typeI_sc6,
                            model_on_second_level_10_knots = improper_typeI_sc12,
                            tT)



################################################################################
# Plot the continuous risk surface in one time interval (first and last)

plt_cont_risk_one_time_interval <- function(dir, risk_surface_filename,
                                            time_interval, t_axis,
                                            polygon_grid2){
  
  ## Load in a specific scenario
  filename = paste(dir, risk_surface_filename, sep = "")
  load(filename)
  tmp_ = risk_surface.list
  
  ## get the times
  t_axis_indices = which(time_interval - 1 < t_axis & 
                           t_axis <= time_interval)
  
  
  if(length(t_axis_indices) == 3){
    ### Plot the continuous true risk surface for times within time_interval
    plt_1 <- heatmap_points(tmp_,
                            polygon_grid2,
                            t_axis[t_axis_indices[1]],
                            title = paste("time: ", 
                                          round(t_axis[t_axis_indices[1]], 1)))
    
    plt_2 <- heatmap_points(tmp_,
                            polygon_grid2,
                            t_axis[t_axis_indices[2]],
                            title = paste("time: ", 
                                          round(t_axis[t_axis_indices[2]], 1)))
    
    plt_3 <- heatmap_points(tmp_,
                            polygon_grid2,
                            t_axis[t_axis_indices[3]],
                            title = paste("time: ", 
                                          round(t_axis[t_axis_indices[3]], 1)))
    
    
    
    ggarrange(plt_1, plt_2, plt_3,
              ncol = 3, nrow = 1,
              common.legend = T, legend = "right")
    
  } else {
    print("Unsuited")
  }
  
  
}

dir = "./Data/Simulated_risk_surfaces/"

##########
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc1_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc2_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc3_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc4_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc5_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
# Save as 15 by 6 name: continuous_scenario_6_time_1
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc6_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc7_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc8_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)

#Save to pdf as 15 by 6. name: continuous_scenario_9_time_1
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc9_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)

plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc10_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc11_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc12_risk_surface_2.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)


##########
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc1_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc2_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc3_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc4_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc5_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)

# Save as 15 by 6, name: continuous_scenario_6_time_12
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc6_risk_surface_2.RData",
                                time_interval = 12, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc7_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc8_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)

#Save to pdf as 15 by 6. name: continuous_scenario_9_time_12
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc9_risk_surface_2.RData",
                                time_interval = 12, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)

plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc10_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc11_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc12_risk_surface_2.RData",
                                time_interval = 13, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)










# Plot the continuous risk surface at three specific times (t = 0.3, 7, 12.5)




################################################################################

# the discrete rate resulting from this AND sampled_count/E_it for the regions
# Add the estimated rate for some of the models AND standard devation (or similar)


## Plot the discrete rate and the sampled-count/E_it for first-level
tmp_ = lambda_sc9.df[lambda_sc9.df$time_id == 1, ]
tmp_map_ = first_level_admin_map
tmp_map_$lambda = tmp_$lambda_it; tmp_map_$E_it = tmp_$E_it; tmp_map_$counts = tmp_$sampled_counts

plt_lambda <- heatmap_areas(tmp_map_, tmp_map_$lambda, title = "Rate")
plt_count_div_e_it <- heatmap_areas(tmp_map_, tmp_map_$counts / tmp_map_$E_it, title = TeX(r'($Y_{it}/E_{it}$)'))

# Save to pdf as 10 by 6, name: sc9_rate_and_counts_div_pop
ggarrange(plt_lambda, NULL, plt_count_div_e_it,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)

## Plot the discrete rate and the sampled-count/E_it for first-level
tmp_ = lambda_sc9.df[lambda_sc9.df$time_id == 12, ]
tmp_map_ = first_level_admin_map
tmp_map_$lambda = tmp_$lambda_it; tmp_map_$E_it = tmp_$E_it; tmp_map_$counts = tmp_$sampled_counts

plt_lambda <- heatmap_areas(tmp_map_, tmp_map_$lambda, title = "Rate")
plt_count_div_e_it <- heatmap_areas(tmp_map_, tmp_map_$counts / tmp_map_$E_it, title = TeX(r'($Y_{it}/E_{it}$)'))

# Save to pdf as 10 by 6, name: sc9_rate_and_counts_div_pop_2
ggarrange(plt_lambda, NULL, plt_count_div_e_it,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)


## Plot the estimated rate and std. deviation of rate for model on first-level
tmp_ = improper_typeIV_sc9$summary.fitted.values$mean[1:nrow(first_level_admin_map)]
tmp_map_ = first_level_admin_map
plt_fitted = heatmap_areas(tmp_map_, tmp_, title = "Posterior mean rate")

tmp_ = improper_typeIV_sc9$summary.fitted.values$sd[1:nrow(first_level_admin_map)]
plt_fitted_sd = heatmap_areas(tmp_map_, tmp_, title = "Posterior std deviation of rate")

# Save to pdf as 10 by 6, name: fitted_rate_n_sd_typeIV_sc9
ggarrange(plt_fitted, NULL, plt_fitted_sd,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)


## Plot the estimated rate and std. deviation of rate for model on first-level
tmp_ = improper_typeIV_sc9$summary.fitted.values$mean[(11 * nrow(first_level_admin_map) + 1):(12 * nrow(first_level_admin_map))]
tmp_map_ = first_level_admin_map
plt_fitted = heatmap_areas(tmp_map_, tmp_, title = "Posterior mean rate")

tmp_ = improper_typeIV_sc9$summary.fitted.values$sd[(11 * nrow(first_level_admin_map) + 1):(12 * nrow(first_level_admin_map))]
plt_fitted_sd = heatmap_areas(tmp_map_, tmp_, title = "Posterior std deviation of rate")

# Save to pdf as 10 by 6, name: fitted_rate_n_sd_typeIV_sc9_2
ggarrange(plt_fitted, NULL, plt_fitted_sd,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)




## Plot the discrete rate and the sampled-count/E_it for second_level
tmp_ = lambda_sc6.df[lambda_sc6.df$time_id == 1, ]
tmp_map_ = second_level_admin_map
tmp_map_$lambda = tmp_$lambda_it; tmp_map_$E_it = tmp_$E_it; tmp_map_$counts = tmp_$sampled_counts

plt_lambda <- heatmap_areas(tmp_map_, tmp_map_$lambda, title = "Rate")
plt_count_div_e_it <- heatmap_areas(tmp_map_, tmp_map_$counts / tmp_map_$E_it, title = TeX(r'($Y_{it}/E_{it}$)'))

# Save to pdf as 10 by 6, name: Rate_and_counts_div_pop_sc6
ggarrange(plt_lambda, NULL, plt_count_div_e_it,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)


## Plot fits
tmp_ = improper_typeIV_sc6$summary.fitted.values$mean[1:nrow(second_level_admin_map)]
tmp_map_ = second_level_admin_map
plt_fitted = heatmap_areas(tmp_map_, tmp_, title = "Posterior mean rate")

tmp_ = improper_typeIV_sc6$summary.fitted.values$sd[1:nrow(second_level_admin_map)]
plt_fitted_sd = heatmap_areas(tmp_map_, tmp_, title = "Posterior std deviation of rate")

# Save to pdf as 10 by 6, name: fitted_rate_n_sd_typeIV_sc6
ggarrange(plt_fitted, NULL, plt_fitted_sd,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)



## Plot the discrete rate and the sampled-count/E_it for second_level
tmp_ = lambda_sc6.df[lambda_sc6.df$time_id == 12, ]
tmp_map_ = second_level_admin_map
tmp_map_$lambda = tmp_$lambda_it; tmp_map_$E_it = tmp_$E_it; tmp_map_$counts = tmp_$sampled_counts

plt_lambda <- heatmap_areas(tmp_map_, tmp_map_$lambda, title = "Rate")
plt_count_div_e_it <- heatmap_areas(tmp_map_, tmp_map_$counts / tmp_map_$E_it, title = TeX(r'($Y_{it}/E_{it}$)'))

# Save to pdf as 10 by 6, name: Rate_and_counts_div_pop_sc6_2
ggarrange(plt_lambda, NULL, plt_count_div_e_it,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)


## Plot fits
tmp_ = improper_typeIV_sc6$summary.fitted.values$mean[(11 * nrow(second_level_admin_map) + 1):(12 * nrow(second_level_admin_map))]
tmp_map_ = second_level_admin_map
plt_fitted = heatmap_areas(tmp_map_, tmp_, title = "Posterior mean rate")

tmp_ = improper_typeIV_sc6$summary.fitted.values$sd[(11 * nrow(second_level_admin_map) + 1):(12 * nrow(second_level_admin_map))]
plt_fitted_sd = heatmap_areas(tmp_map_, tmp_, title = "Posterior std deviation of rate")

# Save to pdf as 10 by 6, name: fitted_rate_n_sd_typeIV_sc6_2
ggarrange(plt_fitted, NULL, plt_fitted_sd,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)


################################################################################
# Plot time-series plot for some regions w. fitted rate against
# The true rate AND sampled_count/E_{it}
regions = c(1, 2, 3, 4,
            5, 6, 7, 8)
# Save to pdf as 10 by 10   name: first_level_typeIV_time_series_sc5
select_regions_lin_pred_vs_true(true_risk = lambda_sc5.df,
                                model = improper_typeIV_sc5,
                                regions = regions,
                                n = nrow(first_level_admin_map),
                                tT = tT)

regions = c(9, 10, 11, 12,
            13, 14, 15, 16)
# Save to pdf as 10 by 10   name: first_level_typeIV_time_series_sc5_2
select_regions_lin_pred_vs_true(true_risk = lambda_sc5.df,
                                model = improper_typeIV_sc5,
                                regions = regions,
                                n = nrow(first_level_admin_map),
                                tT = tT)

## Plot the fitted linear pred. for 6 areas over time w. 95 %CIs against true values
regions = c(1, 50, 100,
            150, 200, 250,
            300, 350)

#region_time_series(true_risk = lambda.df, model = improper_noInt,
#                   region = regions[6], n = nrow(germany_map_2), tT = tT)

## Plot the fitted linear predictor vs the true values for the regions = regions
# Save as pdf as 10 by 10   name: second_level_typeIV_time_series_sc6
select_regions_lin_pred_vs_true(true_risk = lambda_sc6.df,
                                model = improper_typeIV_sc6,
                                regions = regions,
                                n = nrow(second_level_admin_map),
                                tT = tT)



################################################################################
# 












