#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("geofacet")
library(ggh4x)
library(ggridges)
library(latex2exp)

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)




#Save as B_spline_basis_for_time 7.5 by 3.5
#plot(Bt[, 1] ~ xt, type = "l", ylim = c(0, 1), ylab = TeX(r'($B_{t}(t')$)'),
#     xlab = "t'")
#for(i in 2:dim(Bt)[2]){lines(Bt[, i] ~ xt)}


################################################################################
#Load in simulation data (done stupidly, change!)
load("./Simulated_data/sc1/sc1_data.RData")
lambda_sc1.df <- lambda.df
lambda_sc1.df$space.time = 1:nrow(lambda_sc1.df)

load("./Simulated_data/sc2/sc2_data.RData")
lambda_sc2.df <- lambda.df
lambda_sc2.df$space.time = 1:nrow(lambda_sc2.df)

load("./Simulated_data/sc3/sc3_data.RData")
lambda_sc3.df <- lambda.df
lambda_sc3.df$space.time = 1:nrow(lambda_sc3.df)

load("./Simulated_data/sc4/sc4_data.RData")
lambda_sc4.df <- lambda.df
lambda_sc4.df$space.time = 1:nrow(lambda_sc4.df)

load("./Simulated_data/sc5/sc5_data.RData")
lambda_sc5.df <- lambda.df
lambda_sc5.df$space.time = 1:nrow(lambda_sc5.df)

load("./Simulated_data/sc6/sc6_data.RData")
lambda_sc6.df <- lambda.df
lambda_sc6.df$space.time = 1:nrow(lambda_sc6.df)

load("./Simulated_data/sc7/sc7_data.RData")
lambda_sc7.df <- lambda.df
lambda_sc7.df$space.time = 1:nrow(lambda_sc7.df)

load("./Simulated_data/sc8/sc8_data.RData")
lambda_sc8.df <- lambda.df
lambda_sc8.df$space.time = 1:nrow(lambda_sc8.df)

load("./Simulated_data/sc9/sc9_data.RData")
lambda_sc9.df <- lambda.df
lambda_sc9.df$space.time = 1:nrow(lambda_sc9.df)

load("./Simulated_data/sc10/sc10_data.RData")
lambda_sc10.df <- lambda.df
lambda_sc10.df$space.time = 1:nrow(lambda_sc10.df)

load("./Simulated_data/sc11/sc11_data.RData")
lambda_sc11.df <- lambda.df
lambda_sc11.df$space.time = 1:nrow(lambda_sc11.df)

load("./Simulated_data/sc12/sc12_data.RData")
lambda_sc12.df <- lambda.df
lambda_sc12.df$space.time = 1:nrow(lambda_sc12.df)

################################################################################
# Plot the maps themselves

exceptonial_areas_adm1_names = c("Bremen",
                                 "Hamburg",
                                 "Berlin",
                                 "Saarland",
                                 "Niedersachsen",
                                 "Brandenburg")

plt_first_level <- ggplot() + 
  geom_sf(data = first_level_admin_map,
          aes(), 
          alpha = 1,
          color="black") + 
  #geom_sf_label(data = first_level_admin_map[!(first_level_admin_map$NAME_1 %in%
  #                                               exceptonial_areas_adm1_names), ],
  #              aes(label = NAME_1), colour = "black", size = 2) + 
  #ggrepel::geom_label_repel(
  #  data = first_level_admin_map[first_level_admin_map$NAME_1 %in%
  #                                   exceptonial_areas_adm1_names, ],
  #  aes(label = NAME_1, geometry = geometry),
  #  stat = "sf_coordinates",
  #  min.segment.length = 0,
  #  colour = "black",
  #  segment.colour = "black",
  #  size = 2) + 
  ggrepel::geom_label_repel(
    data = first_level_admin_map,
    aes(label = NAME_1, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0.3,
    colour = "black",
    segment.colour = "black",
    size = 3.5) + 
  
  theme(plot.title = element_text(size = 15,  hjust = 0.5),
        axis.title.x = element_blank(), #Remove axis and background grid
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin =  unit(c(0, 0, 0, 0), "inches"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.spacing = unit(1, 'lines')) +
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right"))

plt_first_level


plt_second_level <- ggplot(data = second_level_admin_map) + 
  geom_sf(aes(), 
          alpha = 1,
          color="black") +  
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




plt_polygon_grid <- ggplot(data = polygon_grid2) + 
  geom_sf(aes(), 
          alpha = 1,
          color="black") +  
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

#Save as polygon_grid 7.5 by 3.5
plt_polygon_grid

################################################################################
# Plot two continuous risk surfaces at similar times w. different amount of knots
load("./Data/Simulated_risk_surfaces/sc1_risk_surfaces.RData")
risk_surface.list_sc1 = risk_surface.list[, c("t", "polygon_id", "geometry", "x", "y", 
                                              "unique_id", "time_id", 
                                              "first_level_area_id_mapping", 
                                              "second_level_area_id_mapping")]

risk_surface.list_sc1$values = risk_surface.list$values[, 1]

rm(risk_surface.list)
gc()

load("./Data/Simulated_risk_surfaces/sc7_risk_surfaces.RData")
risk_surface.list_sc7 = risk_surface.list[, c("t", "polygon_id", "geometry", "x", "y", 
                                              "unique_id", "time_id", 
                                              "first_level_area_id_mapping", 
                                              "second_level_area_id_mapping")]

risk_surface.list_sc7$values = risk_surface.list$values[, 1]

rm(risk_surface.list)
gc()

plt_sc1 <- heatmap_points(risk_surface.list_sc1,
                          polygon_grid2,
                          admin_map = germany_border,
                          t_axis[1],
                          title = "20 X 20 Knots")


plt_sc7 <- heatmap_points(risk_surface.list_sc7,
                          polygon_grid2,
                          admin_map = germany_border,
                          t_axis[1],
                          title = "10 X 10 Knots")

# Save to pdf as 10 by 6, name: continuous_risk_20_vs_10
ggarrange(plt_sc1, NULL, plt_sc7,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)

################################################################################
# Plot the temporal trends


xlab = "Year"
ylab = TeX(r'($\bar{Y}_{t}$)')

plt1 <- plot_temporal_trend_data_one_data_set("sc1", 3, "ADM1: const trend, short range", xlab = NULL, ylab = ylab)
plt2 <- plot_temporal_trend_data_one_data_set("sc2", 3, "ADM4: const trend, short range", xlab = NULL, ylab = ylab)
plt3 <- plot_temporal_trend_data_one_data_set("sc3", 3, "ADM1: linear trend, short range", xlab = NULL, ylab = NULL)
plt4 <- plot_temporal_trend_data_one_data_set("sc4", 3, "ADM4: linear trend, short range", xlab = NULL, ylab = NULL)
plt5 <- plot_temporal_trend_data_one_data_set("sc5", 3, "ADM1: change point, short range", xlab = NULL, ylab = NULL)
plt6 <- plot_temporal_trend_data_one_data_set("sc6", 3, "ADM4: change point, short range", xlab = NULL, ylab = NULL)
plt7 <- plot_temporal_trend_data_one_data_set("sc7", 3, "ADM1: const trend, long range", xlab = NULL, ylab = ylab)
plt8 <- plot_temporal_trend_data_one_data_set("sc8", 3, "ADM4: const trend, long range", xlab = xlab, ylab = ylab)
plt9 <- plot_temporal_trend_data_one_data_set("sc9", 3, "ADM1: linear trend, long range", xlab = NULL, ylab = NULL)
plt10 <- plot_temporal_trend_data_one_data_set("sc10", 3, "ADM4: linear trend, long range", xlab = xlab, ylab = NULL)
plt11 <- plot_temporal_trend_data_one_data_set("sc11", 3, "ADM1: change point, long range", xlab = NULL, ylab = NULL)
plt12 <- plot_temporal_trend_data_one_data_set("sc12", 3, "ADM4: change point, long range", xlab = xlab, ylab = NULL)


# Save to pdf as aggregated_temporal_trends: 10.5 by 10.5 
ggarrange(plt1, plt3, plt5, 
          plt7, plt9, plt11,
          plt2, plt4, plt6,
          plt8, plt10, plt12,
          ncol = 3, nrow = 4,
          common.legend = T, 
          legend = "top")

# Save to pdf as 12 by 9: const_temporal_trend_fitted
#plot_fitted_temporal_trends(model_on_first_level_20_knots = improper_typeI_sc1,
#                            model_on_first_level_10_knots = improper_typeI_sc7,
#                            model_on_second_level_20_knots = improper_typeI_sc2,
#                            model_on_second_level_10_knots = improper_typeI_sc8,
#                            tT)

# Save to pdf as 12 by 9: lin_increasing_temporal_trend_fitted
#plot_fitted_temporal_trends(model_on_first_level_20_knots = improper_typeI_sc3,
#                            model_on_first_level_10_knots = improper_typeI_sc9,
#                            model_on_second_level_20_knots = improper_typeI_sc4,
#                            model_on_second_level_10_knots = improper_typeI_sc10,
#                            tT)

# Save to pdf as 12 by 9: change_point_temporal_trend_fitted
#plot_fitted_temporal_trends(model_on_first_level_20_knots = improper_typeI_sc5,
#                            model_on_first_level_10_knots = improper_typeI_sc11,
#                            model_on_second_level_20_knots = improper_typeI_sc6,
#                            model_on_second_level_10_knots = improper_typeI_sc12,
#                            tT)



################################################################################
# Plot the continuous risk surface in one time interval (first and last)

plt_cont_risk_one_time_interval <- function(dir, risk_surface_filename,
                                            time_interval, t_axis,
                                            admin_map,
                                            polygon_grid2){
  
  ## Load in a specific scenario
  filename = paste(dir, risk_surface_filename, sep = "")
  load(filename)
  col_names <- colnames(risk_surface.list)[colnames(risk_surface.list) != "values"]
  tmp_ = risk_surface.list[, col_names]
  tmp_$values = risk_surface.list$values[, 1]
  rm(risk_surface.list)
  
  ## get the times
  t_axis_indices = which(time_interval - 1 < t_axis & 
                           t_axis <= time_interval)
  
  
  if(length(t_axis_indices) == 3){
    ### Plot the continuous true risk surface for times within time_interval
    plt_1 <- heatmap_points(tmp_,
                            polygon_grid2,
                            admin_map,
                            t_axis[t_axis_indices[1]],
                            title = paste("time: ", 
                                          round(t_axis[t_axis_indices[1]], 1)))
    
    plt_2 <- heatmap_points(tmp_,
                            polygon_grid2,
                            admin_map,
                            t_axis[t_axis_indices[2]],
                            title = paste("time: ", 
                                          round(t_axis[t_axis_indices[2]], 1)))
    
    plt_3 <- heatmap_points(tmp_,
                            polygon_grid2,
                            admin_map,
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
                                risk_surface_filename = "sc1_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                first_level_admin_map,
                                polygon_grid2 = polygon_grid2)

plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc2_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc3_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc4_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                polygon_grid2 = polygon_grid2)
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc5_risk_surfaces.RData",
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

tmp_map_$lambda = tmp_$lambda_it[, 1]; tmp_map_$E_it = tmp_$E_it; tmp_map_$counts = tmp_$sampled_counts[, 1]

plt_lambda <- heatmap_areas(tmp_map_, tmp_map_$lambda, title = "Rate")
plt_count_div_e_it <- heatmap_areas(tmp_map_, tmp_map_$counts / tmp_map_$E_it, title = TeX(r'($Y_{it}/E_{it}$)'))

# Save to pdf as 10 by 6, name: sc9_rate_and_counts_div_pop
ggarrange(plt_lambda, NULL, plt_count_div_e_it,
          ncol = 3, nrow = 1, widths = c(1, 0.05, 1),
          common.legend = F)

## Plot the discrete rate and the sampled-count/E_it for first-level
tmp_ = lambda_sc9.df[lambda_sc9.df$time_id == 12, ]
tmp_map_ = first_level_admin_map
tmp_map_$lambda = tmp_$lambda_it[, 1]; tmp_map_$E_it = tmp_$E_it; tmp_map_$counts = tmp_$sampled_counts[, 1]

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
# Plot the models fitted values over time for regions (for ADM1 plot all regions w. geofacet)

##############################
#ADM 1

# Dataset choosen
dataset_id = 3 


# Create ADM1 grid for geofacet
ADM1_grid <- data.frame(
  row = c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5),
  col = c(3, 4, 4, 6, 3, 2, 6, 3, 4, 3, 2, 6, 4, 1, 4, 3),
  code = c("15", "8", "4", "6", "9", "10", "5", "13", "14", "7", "11", "3", "16", "12", "2", "1"),
  name = c("Schleswig-H.", "Mecklenburg-V.", "Brandenburg", "Hamburg", "Niedersachsen", "Nordrhein-W.", "Bremen", "Sachsen-A.", "Sachsen", "Hessen", "Rheinland-P.", "Berlin", "Thuringen", "Saarland", "Bayern", "Baden-W."),
  stringsAsFactors = FALSE #olstein orpommern estfalen nhalt falz urttemberg
)
geofacet::grid_preview(ADM1_grid)


#####
#rate

plt_linpredictor_vs_true_rate <- function(geofacet_grid,
                                            pred_to_plot,
                                            title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  if(nrow(geofacet_grid) == 16){
    theme <- theme(axis.title=element_text(size=12),
                   plot.title = element_text(hjust = 0.5, size=12),
                   strip.text.x = element_text(size = 10),
                   legend.position = c(0.15, 1),
                   legend.justification = c("right", "top"),
                   legend.box.just = "right",
                   legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  } else{
    theme <- theme(axis.title=element_text(size=12),
                   plot.title = element_text(hjust = 0.5, size=12),
                   strip.text.x = element_text(size = 10),
                   legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  }
  
  plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = "Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_point(aes(x = time_id, y = lambda_it, col = "True rate")) + 
    geom_point((aes(x = time_id, y = count_div_pop, col = "sampled count/Eit"))) + # TeX(r'($Y_{it}/E_{it}$)')
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) + 
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue"),
                                  labels = unname(TeX(c("Posterior 95% CI",
                                                        "Posterior median risk",
                                                        "Simulated rate: $\\lambda_{it}$",
                                                        "Simulated count $\\frac{Y_{it}}{E_{it}}$")))) # 
  return(plt)
}

plt_linpredictor_vs_true_rate_improper <- function(ADM_grid, lambda_, model, title){
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = model$summary.fitted.values$'0.5quant', 
                             quantile_0.025 = model$summary.fitted.values$'0.025quant', 
                             quantile_0.975 = model$summary.fitted.values$'0.975quant', 
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it,
                             count_div_pop = lambda_$sampled_counts/lambda_$E_it)
  
  plt_linpredictor_vs_true_rate(ADM_grid, pred_to_plot, title)
  
}

plt_linpredictor_vs_true_rate_proper <- function(ADM_grid, lambda_, model, title){
  ### NB: Must sort the proper ones
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = sort_proper_fitted(model$summary.fitted.values$'0.5quant', length(unique(lambda_$area_id)), tT), # model$summary.fitted.values$'0.5quant',
                             quantile_0.025 = sort_proper_fitted(model$summary.fitted.values$'0.025quant', length(unique(lambda_$area_id)), tT),# model$summary.fitted.values$'0.025quant',
                             quantile_0.975 = sort_proper_fitted(model$summary.fitted.values$'0.975quant', length(unique(lambda_$area_id)), tT), #model$summary.fitted.values$'0.975quant',
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it,
                             count_div_pop = lambda_$sampled_counts/lambda_$E_it)
  
  plt_linpredictor_vs_true_rate(ADM_grid, pred_to_plot, title)
}


wrapper_plt_linpredictor_vs_true_rate <- function(ADM_grid, 
                                                scenario_name,
                                                dataset_id,
                                                model,
                                                improper = T,
                                                title){
  
  ### Load in simulated data for that scenario
  load(paste("./Simulated_data/", scenario_name, "/", 
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
  lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]
  
  #### Remove unessecary data
  rm(lambda.df)
  
  ## Actually plot it
  if(improper){
    plt_linpredictor_vs_true_rate_improper(ADM_grid, lambda_, model, title)
  } else{
    plt_linpredictor_vs_true_rate_proper(ADM_grid, lambda_, model, title)
  }
}


scenario_name = "sc1"
title = "Scenario: const trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_lin_pred_vs_true_sc1 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                    Improper1_typeI_ADM1, improper = T,
                                    title = paste("Model: Improper1_typeI", title, sep = "   "))


scenario_name = "sc3"
title = "Scenario: linear trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model select_regions_lin_pred_vs_true_sc3 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                    proper2_onlyInt_ADM1, improper = F,
                                    title = paste("Model: Proper2_onlyInt", title, sep = "   "))

scenario_name = "sc5"
title = "Scenario: change point, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                    Improper2_typeIV_ADM1, improper = T,
                                    title = paste("Model: Improper2_typeIV", title, sep = "   "))

scenario_name = "sc7"
title = "Scenario: const trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                    Improper1_typeIII_ADM1, improper = T,
                                    title = paste("Model: Improper1_typeIII", title, sep = "   "))


scenario_name = "sc9"
"Scenario: linear trend, short range (ADM1)"
title = "Scenario: linear trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                    proper1_full_ADM1, improper = F,
                                    title = paste("Model: proper1_full", title, sep = "   "))

scenario_name = "sc11"
title = "Scenario: change point, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                    Improper1_noInt_ADM1, improper = T,
                                    title = paste("Model: Improper1_noInt", title, sep = "   "))





#####
# Count

plt_predcount_vs_true_count <- function(geofacet_grid,
                                        pred_to_plot,
                                        title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true counts
  
  plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
    geom_point((aes(x = time_id, y = lambda_it, col = "True rate per 100"))) +
    geom_line(aes(x = time_id, y = post_mean, col = "Posterior mean risk")) + 
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=12),
          plot.title = element_text(hjust = 0.5, size=12),
          strip.text.x = element_text(size = 10),
          legend.position = c(0.15, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue"),
                                  labels = unname(TeX(c("Posterior 95% CI",
                                                        "Posterior mean $\\hat{Y}_{it}$",
                                                        "True $Y_{it}$",
                                                        "Expected true $Y_{it}$")))) # 
  return(plt)
}

wrapper_plt_predcount_vs_true_count <- function(ADM1_grid, 
                                                scenario_name,
                                                dataset_id,
                                                model,
                                                improper = T,
                                                title){
  
  ### Load in simulated data for that scenario
  load(paste("./Simulated_data/", scenario_name, "/", 
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
  lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]
  
  #### Remove unessecary data
  rm(lambda.df)
  
  ## If proper sort the marginals
  if(!improper){
    model$marginals.fitted.values <- sort_proper_fitted(model$marginals.fitted.values,
                                                        length(unique(lambda_$area_id)),
                                                        tT)
  } 
    
  ## Get the upper, lower, and median quantile for the pred. counts
  ul_each <- lapply(model$marginals.fitted.values, 
                    FUN = function(x){
                      return(find_ul_quants_counts_single_pred(x, 100)) #100 aka E_it
                    })
  
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = rep(NA, nrow(lambda_)),
                             post_mean = rep(NA, nrow(lambda_)),
                             quantile_0.025 = rep(NA, nrow(lambda_)),
                             quantile_0.975 = rep(NA, nrow(lambda_)),
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it * lambda_$E_it,
                             count_div_pop = lambda_$sampled_counts/lambda_$E_it)
  
  for(i in 1:nrow(pred_to_plot)){
    pred_to_plot[i, ]$median = ul_each[[i]]$median
    pred_to_plot[i, ]$post_mean = ul_each[[i]]$post_mean
    pred_to_plot[i, ]$quantile_0.025 = ul_each[[i]]$l
    pred_to_plot[i, ]$quantile_0.975 = ul_each[[i]]$u
  }
  
  plt_predcount_vs_true_count(ADM1_grid, pred_to_plot, title)
}


scenario_name = "sc1"
title = "Scenario: const trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_lin_pred_vs_true_sc1 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_typeI_ADM1, improper = T,
                                      title = paste("Model: Improper1_typeI", title, sep = "   "))


scenario_name = "sc3"
title = "Scenario: linear trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model select_regions_lin_pred_vs_true_sc3 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper2_onlyInt_ADM1, improper = F,
                                      title = paste("Model: Proper2_onlyInt", title, sep = "   "))

scenario_name = "sc5"
title = "Scenario: change point, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_typeIV_ADM1, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

scenario_name = "sc7"
title = "Scenario: const trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_typeIII_ADM1, improper = T,
                                      title = paste("Model: Improper1_typeIII", title, sep = "   "))


scenario_name = "sc9"
"Scenario: linear trend, short range (ADM1)"
title = "Scenario: linear trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_full_ADM1, improper = F,
                                      title = paste("Model: proper1_full", title, sep = "   "))

scenario_name = "sc11"
title = "Scenario: change point, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_noInt_ADM1, improper = T,
                                      title = paste("Model: Improper1_noInt", title, sep = "   "))



###################
# ADM4

dataset_id_2 = 4

ADM4_grid <- data.frame(
  code = c("1", "50", "100", "150", "200", "250", "300", "350", "400"),
  name = c("Alb-Donau-K.",
           "Ansbach",
           "Munchen",
           "Oberhavel",
           "Celle",
           "Dusseldorf",
           "Bad Kreuznach",
           "Stendal",
           "Wartburgkreis"),
  row = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
  col = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
  stringsAsFactors = FALSE)


geofacet::grid_preview(ADM4_grid)



#####
#rate

scenario_name = "sc2"
title = "Scenario: const trend, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_lin_pred_vs_true_sc1 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper1_typeIV_ADM4, improper = T,
                                      title = paste("Model: Improper1_typeIV", title, sep = "   "))


scenario_name = "sc4"
title = "Scenario: linear trend, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model select_regions_lin_pred_vs_true_sc3 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      proper2_onlyInt_ADM4, improper = F,
                                      title = paste("Model: Proper2_onlyInt", title, sep = "   "))

scenario_name = "sc6"
title = "Scenario: change point, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
wrapper_plt_linpredictor_vs_true_rate(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper2_typeIII_ADM4, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

scenario_name = "sc8"
title = "Scenario: const trend, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save
wrapper_plt_linpredictor_vs_true_rate(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper1_noInt_ADM4, improper = T,
                                      title = paste("Model: Improper1_noInt", title, sep = "   "))


scenario_name = "sc10"
title = "Scenario: linear trend, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model
wrapper_plt_linpredictor_vs_true_rate(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      proper1_onlyInt_ADM4, improper = F,
                                      title = paste("Model: proper1_onlyInt", title, sep = "   "))

scenario_name = "sc12"
title = "Scenario: change point, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
wrapper_plt_linpredictor_vs_true_rate(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper1_typeIV_ADM4, improper = T,
                                      title = paste("Model: Improper1_typeIV", title, sep = "   "))





#####
#count

scenario_name = "sc2"
title = "Scenario: const trend, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_lin_pred_vs_true_sc1 9 by 7
(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper1_typeIV_ADM4, improper = T,
                                      title = paste("Model: Improper1_typeIV", title, sep = "   "))


scenario_name = "sc4"
title = "Scenario: linear trend, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model select_regions_lin_pred_vs_true_sc3 9 by 7
(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      proper2_onlyInt_ADM4, improper = F,
                                      title = paste("Model: Proper2_onlyInt", title, sep = "   "))

scenario_name = "sc6"
title = "Scenario: change point, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper2_typeIII_ADM4, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

scenario_name = "sc8"
title = "Scenario: const trend, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save
(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper1_noInt_ADM4, improper = T,
                                      title = paste("Model: Improper1_noInt", title, sep = "   "))


scenario_name = "sc10"
title = "Scenario: linear trend, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot a proper model
(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      proper1_onlyInt_ADM4, improper = F,
                                      title = paste("Model: proper1_onlyInt", title, sep = "   "))

scenario_name = "sc12"
title = "Scenario: change point, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model
(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                                      Improper1_typeIV_ADM4, improper = T,
                                      title = paste("Model: Improper1_typeIV", title, sep = "   "))






################################################################################
# Create ridgeplots for the interval scores 1, 2, and 3 years ahead for both count and rate


model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
                "Improper1_typeIII", "Improper1_typeIV",
                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
                "Improper2_typeIII", "Improper2_typeIV",
                "proper1_noInt", "proper1_onlyInt", "proper1_full", "proper1_iid",
                "proper2_noInt", "proper2_onlyInt", "proper2_full", "proper2_iid")

scenario_name = "sc2"



#########
# Based on counts
model_name = model_names[1]
load(paste("./results/model_choice/model_choice_", 
           model_name, "_", 
           scenario_name, ".RData",
           sep = ""))

tmp_ <- model_choice_for_counts

to_plot_ridge.df <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                              IS_one_year_ahead = tmp_$IS_1_year_ahead,
                              IS_two_year_ahead = tmp_$IS_2_year_ahead,
                              IS_three_year_ahead = tmp_$IS_3_year_ahead,
                              IS_tot = tmp_$total_IS)


for(model_name in model_names[2:length(model_names)]){
  ### Load in Model choice data
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_counts
  
  tmp2_ <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                      IS_one_year_ahead = tmp_$IS_1_year_ahead,
                      IS_two_year_ahead = tmp_$IS_2_year_ahead,
                      IS_three_year_ahead = tmp_$IS_3_year_ahead,
                      IS_tot = tmp_$total_IS)
  
  to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
  
}




ggplot(to_plot_ridge.df, aes(x = IS_one_year_ahead, 
                             y = model_name, 
                             fill = model_name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") + 
  xlim(12.5, 18)

ggplot(to_plot_ridge.df, aes(x = IS_two_year_ahead, 
                             y = model_name, 
                             fill = model_name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

ggplot(to_plot_ridge.df, aes(x = IS_three_year_ahead, 
                             y = model_name, 
                             fill = model_name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")


#########
# Based on rates
model_name = model_names[1]
load(paste("./results/model_choice/model_choice_", 
           model_name, "_", 
           scenario_name, ".RData",
           sep = ""))

tmp_ <- model_choice_for_rates

to_plot_ridge.df <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                               IS_one_year_ahead = tmp_$IS_1_year_ahead,
                               IS_two_year_ahead = tmp_$IS_2_year_ahead,
                               IS_three_year_ahead = tmp_$IS_3_year_ahead,
                               IS_tot = tmp_$total_IS)


for(model_name in model_names[2:length(model_names)]){
  ### Load in Model choice data
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_rates
  
  tmp2_ <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                      IS_one_year_ahead = tmp_$IS_1_year_ahead,
                      IS_two_year_ahead = tmp_$IS_2_year_ahead,
                      IS_three_year_ahead = tmp_$IS_3_year_ahead,
                      IS_tot = tmp_$total_IS)
  
  to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
  
}




ggplot(to_plot_ridge.df, aes(x = IS_one_year_ahead, 
                             y = model_name, 
                             fill = model_name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

ggplot(to_plot_ridge.df, aes(x = IS_two_year_ahead, 
                             y = model_name, 
                             fill = model_name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

ggplot(to_plot_ridge.df, aes(x = IS_three_year_ahead, 
                             y = model_name, 
                             fill = model_name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

################################################################################
# Posterior distributions









