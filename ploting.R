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

# Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)

dataset_id = 3 
dataset_id_2 = 4

################################################################################
# Splines

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
# Create ridgeplots for the interval scores 1, 2, and 3 years ahead for both count and rate


model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
                "Improper1_typeIII", "Improper1_typeIV",
                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
                "Improper2_typeIII", "Improper2_typeIV",
                "proper1_noInt", "proper1_onlyInt", "proper1_full", "proper1_iid",
                "proper2_noInt", "proper2_onlyInt", "proper2_full", "proper2_iid")


ridgeplot_counts <- function(model_names, scenario_name,
                             one_2_3_or_total, xlab,
                             xlim){
  
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
                                 IS_tot = tmp_$total_IS,
                                 MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                                 MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                                 MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                                 MSE_tot = tmp_$total_mse)
  
  
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
                        IS_tot = tmp_$total_IS,
                        MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                        MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                        MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                        MSE_tot = tmp_$total_mse)
    
    to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
    
  }
  if(one_2_3_or_total == 1){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_one_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab("Model") + 
        xlim(xlim[1], xlim[2])
      
    )
  } else if(one_2_3_or_total == 2){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_two_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.text.x = element_text(size = 10),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) + 
        xlim(xlim[1], xlim[2])
    )
    
  } else if(one_2_3_or_total == 3){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_three_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) + 
        xlim(xlim[1], xlim[2])
    )
  }
}








ridgeplot_rates <- function(model_names, scenario_name,
                             one_2_3_or_total, xlab,
                            xlim){
  
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
                                 IS_tot = tmp_$total_IS,
                                 MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                                 MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                                 MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                                 MSE_tot = tmp_$total_mse)
  
  
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
                        IS_tot = tmp_$total_IS,
                        MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                        MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                        MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                        MSE_tot = tmp_$total_mse)
    
    to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
    
  }
  
  if(one_2_3_or_total == 1){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_one_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.text.x = element_text(size = 10),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab("Model") +
        xlim(xlim[1], xlim[2])
      
    )
  } else if(one_2_3_or_total == 2){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_two_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) +
        xlim(xlim[1], xlim[2])
    )
    
  } else if(one_2_3_or_total == 3){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_three_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) +
        xlim(xlim[1], xlim[2])
    )
  }
}




ridgeplot_counts_and_rates_all_years <- function(scenario_name, xlim_counts, xlim_rates){
  
  plt1 <- ridgeplot_counts(model_names = model_names, scenario_name = scenario_name, 1, "IS count: 1 year ahead", xlim_counts)
  plt2 <- ridgeplot_counts(model_names = model_names, scenario_name = scenario_name, 2, "IS count: 2 years ahead", xlim_counts)
  plt3 <- ridgeplot_counts(model_names = model_names, scenario_name = scenario_name, 3, "IS count: 3 years ahead", xlim_counts)
  
  
  plt4 <- ridgeplot_rates(model_names = model_names, scenario_name = scenario_name, 1, "IS rate: 1 year ahead", xlim_rates)
  plt5 <- ridgeplot_rates(model_names = model_names, scenario_name = scenario_name, 2, "IS rate: 2 years ahead", xlim_rates)
  plt6 <- ridgeplot_rates(model_names = model_names, scenario_name = scenario_name, 3, "IS rate: 3 years ahead", xlim_rates)
  
  plt <- ggarrange(plt1, plt2, plt3, 
                   plt4, plt5, plt6,
                   widths = c(1.35, 1, 1),
                   ncol = 3, nrow = 2, common.legend = F)
  return(plt)
}

#Save as 15 by 15 sc1_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc1", c(8, 27), c(0.02, 0.1))

#Save as 15 by 15 sc2_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc2", c(12, 20), c(0.02, 0.12))

#Save as 15 by 15 sc3_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc3", c(10, 30), c(0.02, 0.12))


#Save as 15 by 15 sc4_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc4", c(12, 22), c(0.02, 0.13))

#Save as 15 by 15 sc5_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc5", c(10, 25), c(0.02, 0.11))

#Save as 15 by 15 sc6_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc6", c(12, 20), c(0.02, 0.13))

#Save as 15 by 15 sc7_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc7", c(8, 27), c(0.02, 0.1))

#Save as 15 by 15 sc8_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc8", c(12, 20), c(0.02, 0.12))

#Save as 15 by 15 sc9_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc9", c(8, 27), c(0.02, 0.1))

#Save as 15 by 15 sc10_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc10", c(12, 20), c(0.02, 0.12))

#Save as 15 by 15 sc11_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc11", c(8, 27), c(0.02, 0.1))

#Save as 15 by 15 sc12_IS_ridgeplot_all
ridgeplot_counts_and_rates_all_years("sc12", c(12, 20), c(0.02, 0.12))



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



### SC1
scenario_name = "sc1"
title = "Scenario: const trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_lin_pred_vs_true_sc1_Imp 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_noInt_ADM1, improper = T,
                                      title = paste("Model: Improper2_noInt", title, sep = "   "))


### Plot a proper model. Save as select_regions_lin_pred_vs_true_sc1_prop 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_onlyInt_ADM1, improper = F,
                                      title = paste("Model: proper1_onlyInt", title, sep = "   "))




### SC3

scenario_name = "sc3"
title = "Scenario: linear trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model select_regions_lin_pred_vs_true_sc3_Imp 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_typeII_ADM1, improper = T,
                                      title = paste("Model: Improper1_typeII", title, sep = "   "))


### Plot a proper model select_regions_lin_pred_vs_true_sc3_prop 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper2_noInt_ADM1, improper = F,
                                      title = paste("Model: proper2_noInt", title, sep = "   "))




### SC5

scenario_name = "sc5"
title = "Scenario: change point, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model Save as select_regions_lin_pred_vs_true_sc5_Imp 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_typeIV_ADM1, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

### Plot a proper model: select_regions_lin_pred_vs_true_sc5_prop 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_full_ADM1, improper = F,
                                      title = paste("Model: proper1_full", title, sep = "   "))



### SC7
scenario_name = "sc7"
title = "Scenario: const trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_lin_pred_vs_true_sc7_Imp 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_noInt_ADM1, improper = T,
                                      title = paste("Model: Improper2_noInt", title, sep = "   "))


### Plot a proper model. Save as select_regions_lin_pred_vs_true_sc7_prop 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_onlyInt_ADM1, improper = F,
                                      title = paste("Model: proper1_onlyInt", title, sep = "   "))




### SC9

scenario_name = "sc9"
title = "Scenario: linear trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model select_regions_lin_pred_vs_true_sc9_Imp 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_typeII_ADM1, improper = T,
                                      title = paste("Model: Improper1_typeII", title, sep = "   "))


### Plot a proper model select_regions_lin_pred_vs_true_sc9_prop 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper2_noInt_ADM1, improper = F,
                                      title = paste("Model: proper2_noInt", title, sep = "   "))




### SC11

scenario_name = "sc11"
title = "Scenario: change point, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model Save as select_regions_lin_pred_vs_true_sc11_Imp 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_typeIV_ADM1, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

### Plot a proper model: select_regions_lin_pred_vs_true_sc11_prop 9 by 7
wrapper_plt_linpredictor_vs_true_rate(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_full_ADM1, improper = F,
                                      title = paste("Model: proper1_full", title, sep = "   "))



#####
# Count






### SC1
scenario_name = "sc1"
title = "Scenario: const trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_pred_count_vs_true_sc1_Imp 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_noInt_ADM1, improper = T,
                                      title = paste("Model: Improper2_noInt", title, sep = "   "))


### Plot a proper model. Save as select_regions_pred_count_vs_true_sc1_prop 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_onlyInt_ADM1, improper = F,
                                      title = paste("Model: proper1_onlyInt", title, sep = "   "))




### SC3

scenario_name = "sc3"
title = "Scenario: linear trend, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model select_regions_pred_count_vs_true_sc3_Imp 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_typeII_ADM1, improper = T,
                                      title = paste("Model: Improper1_typeII", title, sep = "   "))


### Plot a proper model select_regions_pred_count_vs_true_sc3_prop 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper2_noInt_ADM1, improper = F,
                                      title = paste("Model: proper2_noInt", title, sep = "   "))




### SC5

scenario_name = "sc5"
title = "Scenario: change point, short range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model Save as select_regions_pred_count_vs_true_sc5_Imp 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_typeIV_ADM1, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

### Plot a proper model: select_regions_pred_count_vs_true_sc5_prop 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_full_ADM1, improper = F,
                                      title = paste("Model: proper1_full", title, sep = "   "))



### SC7
scenario_name = "sc7"
title = "Scenario: const trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model. Save as select_regions_pred_count_vs_true_sc7_Imp 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_noInt_ADM1, improper = T,
                                      title = paste("Model: Improper2_noInt", title, sep = "   "))


### Plot a proper model. Save as select_regions_pred_count_vs_true_sc7_prop 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_onlyInt_ADM1, improper = F,
                                      title = paste("Model: proper1_onlyInt", title, sep = "   "))




### SC9

scenario_name = "sc9"
title = "Scenario: linear trend, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an Improper model select_regions_pred_count_vs_true_sc9_Imp 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper1_typeII_ADM1, improper = T,
                                      title = paste("Model: Improper1_typeII", title, sep = "   "))


### Plot a proper model select_regions_pred_count_vs_true_sc9_prop 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper2_noInt_ADM1, improper = F,
                                      title = paste("Model: proper2_noInt", title, sep = "   "))




### SC11

scenario_name = "sc11"
title = "Scenario: change point, long range (ADM1)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

### Plot an improper model Save as select_regions_pred_count_vs_true_sc11_Imp 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      Improper2_typeIV_ADM1, improper = T,
                                      title = paste("Model: Improper2_typeIV", title, sep = "   "))

### Plot a proper model: select_regions_pred_count_vs_true_sc11_prop 9 by 7
wrapper_plt_predcount_vs_true_count(ADM1_grid, scenario_name, dataset_id = dataset_id,
                                      proper1_full_ADM1, improper = F,
                                      title = paste("Model: proper1_full", title, sep = "   "))




###################
# ADM4

dataset_id_2 = 4

ADM4_grid <- data.frame(
  code = c("1", "50", "100", "150", "200", "250", "300", "350", "400",
           "93", "106", "354"),
  name = c("Alb-Donau-K.", "Ansbach",
           "Munchen", "Oberhavel",
           "Celle", "Dusseldorf", 
           "Bad Kreuznach", "Stendal",
           "Wartburgkreis", "Luchow-Dannenberg",
           "Nurnberg", "Dresden"),
  row = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
  col = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
  stringsAsFactors = FALSE)


geofacet::grid_preview(ADM4_grid)

### SC2

scenario_name = "sc2"
title = "Scenario: const trend, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

# Save as timeseries_sc2 10 by 8.5
wrapper_timeseries_plt(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                       Improper1_typeIV_ADM4, improper = T,
                       title = paste("Model: Improper1_typeIV", title, sep = "   "))


### SC4

scenario_name = "sc4"
title = "Scenario: linear trend, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

# Save as timeseries_sc4 10 by 8.5
wrapper_timeseries_plt(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                       proper1_onlyInt_ADM4, improper = F,
                       title = paste("Model: proper1_onlyInt", title, sep = "   "))



### SC6

scenario_name = "sc6"
title = "Scenario: change point, short range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

# Save as timeseries_sc6 10 by 8.5
wrapper_timeseries_plt(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                       proper2_onlyInt_ADM4, improper = F,
                       title = paste("Model: proper2_onlyInt", title, sep = "   "))


### SC8

scenario_name = "sc8"
title = "Scenario: const trend, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

# Save as timeseries_sc8 10 by 8.5
wrapper_timeseries_plt(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                       Improper2_typeIII_ADM4, improper = T,
                       title = paste("Model: Improper2_typeIII", title, sep = "   "))


### SC10

scenario_name = "sc10"
title = "Scenario: linear trend, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

# Save as timeseries_sc10 10 by 8.5
wrapper_timeseries_plt(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
                       proper1_onlyInt_ADM4, improper = F,
                       title = paste("Model: proper1_onlyInt", title, sep = "   "))


### SC12

scenario_name = "sc12"
title = "Scenario: change point, long range (ADM4)"             #TeX(r'(ADM1$_{const, long}$)')

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

# Save as timeseries_sc12 10 by 8.5
wrapper_timeseries_plt(ADM4_grid, scenario_name, dataset_id = dataset_id_2,
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
# Posterior distributions


################################################################################
# Plot the fits on the maps

#####################
# Plot the continuous risk surface in one time interval (first and last)

dir = "./Data/Simulated_risk_surfaces/"
#####
# ADM1

plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc1_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                first_level_admin_map,
                                polygon_grid2 = polygon_grid2)

plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc3_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                first_level_admin_map,
                                polygon_grid2 = polygon_grid2)

plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc5_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                first_level_admin_map,
                                polygon_grid2 = polygon_grid2)


#####
# ADM 4

### SC2

# Save as cont_risk_sc2_year_1 12.5 by 7.5
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc2_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                second_level_admin_map,
                                polygon_grid2 = polygon_grid2)

# Save as cont_risk_sc2_year_12
plt_cont_risk_one_time_interval(dir = dir,
                            risk_surface_filename = "sc2_risk_surfaces.RData",
                            time_interval = 12, t_axis = t_axis,
                            second_level_admin_map,
                            polygon_grid2 = polygon_grid2)


### SC4
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc4_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                second_level_admin_map,
                                polygon_grid2 = polygon_grid2)

### SC6
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc6_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                second_level_admin_map,
                                polygon_grid2 = polygon_grid2)

### SC8
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc8_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                second_level_admin_map,
                                polygon_grid2 = polygon_grid2)

### SC10
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc10_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                second_level_admin_map,
                                polygon_grid2 = polygon_grid2)

### SC12
plt_cont_risk_one_time_interval(dir = dir,
                                risk_surface_filename = "sc12_risk_surfaces.RData",
                                time_interval = 1, t_axis = t_axis,
                                second_level_admin_map,
                                polygon_grid2 = polygon_grid2)



#####################
# Plot the discrete rate per 100 and the observed counts

plt_true_discrete_rate_four_years <- function(scenario_name, 
                              dataset_id,
                              admin_map,
                              scale_col,
                              scale){
  
  ### Load in the area-specific rate
  load(paste("Simulated_data/", scenario_name, "/",
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]
  lambda_$mu <- lambda_$lambda_it * lambda_$E_it
  
  ### Create bins
  min_mu = min(lambda_$mu); max_mu = max(lambda_$mu)
  hardcoded_bins_mu = round(seq(min_mu, max_mu, length.out = 12), 2)
  
  ### Attach the simulated values to the admin_map for the years 3, 7, 11, and 13
  tmp_map_ = admin_map
  
  #### Plot year 3
  tmp_map_$mu <- lambda_[lambda_$time_id == 3, ]$mu
  plt1 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "3")
  
  #### Plot year 7
  tmp_map_$mu <- lambda_[lambda_$time_id == 7, ]$mu
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "7")
  
  #### Plot year 11
  tmp_map_$mu <- lambda_[lambda_$time_id == 11, ]$mu
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "11")
  
  #### Plot year 13
  tmp_map_$mu <- lambda_[lambda_$time_id == 13, ]$mu
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "13")
  return(ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T,
                   legend = "bottom"))
}

plt_true_counts_four_years <- function(scenario_name, 
                                              dataset_id,
                                              admin_map,
                                              scale_col,
                                              scale){
  
  ### Load in the area-specific rate
  load(paste("Simulated_data/", scenario_name, "/",
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
  
  
  ### Create bins
  min_count = min(lambda_$sampled_counts); max_count = max(lambda_$sampled_counts)
  hardcoded_bins_count = round(seq(min_count, max_count, length.out = 12), 0)
  
  
  ### Attach the simulated values to the admin_map for the years 3, 7, 11, and 13
  tmp_map_ = admin_map
  
  #### Plot year 3
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 3, ]$sampled_counts
  plt1 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "3")
  
  
  #### Plot year 7
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 7, ]$sampled_counts
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "7")
  
  #### Plot year 11
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 11, ]$sampled_counts
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "11")
  
  #### Plot year 13
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 13, ]$sampled_counts
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "13")
  
  return(ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T,
                   legend = "bottom"))
}


scale_col = heat.colors(30, rev=TRUE)
scale = scale_col[seq(3, 30, length.out = 12)]

plt5 <- plt_true_discrete_rate_four_years(scenario_name = "sc2", 
                                          dataset_id = dataset_id_2,
                                          admin_map = second_level_admin_map,
                                          scale_col = scale_col,
                                          scale = scale)

# Save as sc2_true_rate 8.5 by 3.5
annotate_figure(plt5, 
                top = text_grob(TeX(r'(Scenario: ADM4$_{const, short}$, Simulated rate per 100 for years: )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))




plt5 <- plt_true_counts_four_years(scenario_name = "sc2", 
                                          dataset_id = dataset_id_2,
                                          admin_map = second_level_admin_map,
                                          scale_col = scale_col,
                                          scale = scale)

# Save as sc2_true_count 8.5 by 3.5
annotate_figure(plt5,
                top = text_grob(TeX(r'(Scenario: ADM4$_{const, short}$, Simulated counts for years: )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))


# AND the fitted pred. count and standard deviation ()
# Plot the best improper vs best proper four or three years


plt_fitted_rate_four_years <- function(model, 
                                       Improper = T, 
                                       admin_map,
                                       scale_col,
                                       scale,
                                       hardcoded_bins,
                                       overall_title){
  
  #Find the number of areas
  n_ADM = nrow(admin_map)
  
  # If a proper model, must sort the fitted rates
  if(!Improper){
    print("sort")
    model$summary.fitted.values$mean <- sort_proper_fitted(model$summary.fitted.values$mean, 
                                                           n_ADM, tT)
    
  }
  
  
  ### Attach the fitted values to the admin_map for the years 3, 7, 11, and 13
  tmp_map_ = admin_map
  
  
  #### Plot year 3
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (3 - 1) + 1):(n_ADM * 3)] * 100
  plt1 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "3")
  
  #### Plot year 7
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (7 - 1) + 1):(n_ADM * 7)] * 100
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "7")
  
  #### Plot year 11
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (11 - 1) + 1):(n_ADM * 11)] * 100
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "11")
  
  #### Plot year 13
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (13 - 1) + 1):(n_ADM * 13)] * 100
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "13")
  
  plt <- ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T, 
                   legend = "bottom")
  
  plt <- annotate_figure(plt,
                  top = text_grob(overall_title, 
                                  color = "black", 
                                  face = "bold", 
                                  size = 14))
  
  
  
  return(plt)
}

plt_fitted_rate_sd_four_years <- function(model, 
                                       Improper = T, 
                                       admin_map,
                                       scale_col,
                                       scale,
                                       hardcoded_bins,
                                       overall_title){
  
  #Find the number of areas
  n_ADM = nrow(admin_map)
  
  # If a proper model, must sort the fitted rates
  if(!Improper){
    print("sort")
    model$summary.fitted.values$sd <- sort_proper_fitted(model$summary.fitted.values$sd, 
                                                         n_ADM, tT)
  }
  
  
  ### Attach the fitted values to the admin_map for the years 3, 7, 11, and 13
  tmp_map_ = admin_map
  
  
  #### Plot year 3
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (3 - 1) + 1):(n_ADM * 3)] * 100
  plt1 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "3")
  
  #### Plot year 7
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (7 - 1) + 1):(n_ADM * 7)] * 100
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "7")
  
  #### Plot year 11
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (11 - 1) + 1):(n_ADM * 11)] * 100
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "11")
  
  #### Plot year 13
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (13 - 1) + 1):(n_ADM * 13)] * 100
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "13")
  
  plt <- ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T, 
                   legend = "bottom")
  
  plt <- annotate_figure(plt,
                         top = text_grob(overall_title, 
                                         color = "black", 
                                         face = "bold", 
                                         size = 14))
  
  
  
  return(plt)
}


scale_col = heat.colors(30, rev=TRUE)
scale = scale_col[seq(3, 30, length.out = 12)]




#####
#SC2
scenario_name = "sc2"
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

best_imp_model <- Improper1_typeIV_ADM4
best_prop_model <- proper2_onlyInt_ADM4  

hardcoded_bins_mean_rate = round(seq(min(best_prop_model$summary.fitted.values$mean * 100),
                                     max(best_prop_model$summary.fitted.values$mean * 100), 
                                     length.out = 12), 2)


hardcoded_bins_sd_rate = round(seq(min(best_prop_model$summary.fitted.values$sd * 100),
                                   max(best_prop_model$summary.fitted.values$sd * 100), 
                                   length.out = 12), 2)

### Improper

## Plot the mean rate 
# Save as sc2_mean_fitted_rate_Imp 8.5 by 3.5
plt_fitted_rate_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, fitted rate per 100 of Improper1_typeIV for years)'))

## Plot the sd of the rate
# Save as sc2_sd_fitted_rate_Imp 8.5 by 3.5
plt_fitted_rate_sd_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, sd of rate per 100 of Improper1_typeIV for years)'))



### Proper

## Plot the mean rate
# Save as sc2_mean_fitted_rate_prop 8.5 by 3.5
plt_fitted_rate_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, fitted rate per 100 of proper2_onlyInt for years)'))

## Plot the sd of the rate
# Save as sc2_sd_fitted_rate_prop 8.5 by 3.5
plt_fitted_rate_sd_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, sd of rate per 100 of proper2_onlyInt for years)'))




##### 
#SC4



