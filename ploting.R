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
n_ADMnew <- nrow(new_map)
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
# Posterior distributions of hyperparameters and effects

#################################
# ADM1


#################################
# ADM4

#########################
# Effects

#############
#Temporal effects
plot_fitted_temporal_trends <- function(model_on_short_range,
                                        model_on_long_range,
                                        Improper = T,
                                        n_ADM,
                                        tT,
                                        title,
                                        ylim,
                                        ylab_ = T,
                                        y_breaks = NULL){
  
  #
  years = 1:tT
  
  # Create df for plotting
  tmp <- data.frame(years = rep(years, 2), 
                    range = c(rep("short", tT),
                              rep("long", tT)), 
                    l_quant = rep(NA, 2 * tT),
                    median = rep(NA, 2 * tT),
                    u_quant = rep(NA, 2 * tT))
  
  #If Improper model
  if(Improper){
    
    # Scale the effects by the mixing param (mean[2] phi, mean[1] prec.)
    scale_short <- sqrt(model_on_short_range$summary.hyperpar$mean[2]) * 1/sqrt(model_on_short_range$summary.hyperpar$mean[1])
    scale_long <- sqrt(model_on_long_range$summary.hyperpar$mean[2]) * 1/sqrt(model_on_long_range$summary.hyperpar$mean[1])
    
    ### Insert data for short-range
    tmp$l_quant[1:tT] <- model_on_short_range$summary.random$time_id[(tT + 1):(2 * tT), 4] * scale_short
    tmp$median[1:tT] <- model_on_short_range$summary.random$time_id[(tT + 1):(2 * tT), 5] * scale_short
    tmp$u_quant[1:tT] <- model_on_short_range$summary.random$time_id[(tT + 1):(2 * tT), 6] * scale_short
    
    ### Insert data for long-range
    tmp$l_quant[(tT+1):(2*tT)] <- model_on_long_range$summary.random$time_id[(tT + 1):(2 * tT), 4] * scale_long
    tmp$median[(tT+1):(2*tT)] <- model_on_long_range$summary.random$time_id[(tT + 1):(2 * tT), 5] * scale_long
    tmp$u_quant[(tT+1):(2*tT)] <- model_on_long_range$summary.random$time_id[(tT + 1):(2 * tT), 6] * scale_long
    
    if(ylab_){
      ylab = expression(alpha[t])
    } else {
      ylab = NULL
    }
    
    
  }else{ # If proper
    fixed_effect_short <- years * model_on_short_range$summary.fixed$mean[2]
    fixed_effect_long <- years * model_on_long_range$summary.fixed$mean[2]
    
    # Insert data into df for plotting
    tmp$l_quant[1:tT] <- model_on_short_range$summary.random$time_id.copy[1:tT, 4] + fixed_effect_short
    tmp$median[1:tT] <- model_on_short_range$summary.random$time_id.copy[1:tT, 5] + fixed_effect_short
    tmp$u_quant[1:tT] <- model_on_short_range$summary.random$time_id.copy[1:tT, 6] + fixed_effect_short
    
    tmp$l_quant[(tT+1):(2*tT)] <- model_on_long_range$summary.random$time_id.copy[1:tT, 4] + fixed_effect_long
    tmp$median[(tT+1):(2*tT)] <- model_on_long_range$summary.random$time_id.copy[1:tT, 5] + fixed_effect_long
    tmp$u_quant[(tT+1):(2*tT)]<- model_on_long_range$summary.random$time_id.copy[1:tT, 6] + fixed_effect_long
    
    if(ylab_){
      ylab = expression(alpha[t] + beta*t)
    } else {
      ylab = NULL
    }
  }
  
  #Make a plot
  plt <- ggplot(data = tmp) + ggtitle(title) +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line(aes(years, median, col = range), linewidth = 1) +
    geom_line(aes(years, l_quant, col = range), linewidth = 0.6, linetype = "dashed") + 
    geom_line(aes(years, u_quant, col = range), linewidth = 0.6, linetype = "dashed") + 
    geom_vline(xintercept = 10.5, linetype = "dashed", 
               color = "darkgrey", linewidth = 0.6) +
    xlab("Year") + ylab(ylab) + 
    ylim(ylim) + guides(col = guide_legend(title="Spatial range"))
  
  if(is.null(y_breaks)){
    return(plt)
  } else {
    plt <- plt + scale_y_continuous(breaks = y_breaks)
    return(plt)
  }
  
  
  
}

#######
#Improper

# Constant trend
load("diagnostics_sc2.RData")
model_on_short_range = Improper1_noInt_ADM4

load("diagnostics_sc8.RData")
model_on_long_range = Improper1_noInt_ADM4


plt_const_imp <- plot_fitted_temporal_trends(model_on_short_range,
                                             model_on_long_range,
                                             Improper = T, 
                                             n_ADM4,
                                             tT,
                                             "Constant temporal trend",
                                             ylim = c(-0.5, 0.5))

# Linear trend
load("diagnostics_sc4.RData")
model_on_short_range = Improper1_noInt_ADM4

load("diagnostics_sc10.RData")
model_on_long_range = Improper1_noInt_ADM4


plt_lin_imp <- plot_fitted_temporal_trends(model_on_short_range,
                            model_on_long_range,
                            Improper = T, 
                            n_ADM4,
                            tT,
                            "Linear temporal trend",
                            ylim = c(-0.5, 0.5),
                            ylab = F)

# CP
load("diagnostics_sc6.RData")
model_on_short_range = Improper1_noInt_ADM4

load("diagnostics_sc12.RData")
model_on_long_range = Improper1_noInt_ADM4


plt_cp_imp <- plot_fitted_temporal_trends(model_on_short_range,
                            model_on_long_range,
                            Improper = T, 
                            n_ADM4,
                            tT,
                            "Change point temporal trend",
                            ylim = c(-0.5, 0.5),
                            ylab = F)

# Save as temporal_trends_improper 9 by 3
ggarrange(plt_const_imp, plt_lin_imp, plt_cp_imp,
          ncol = 3, nrow = 1, common.legend = T, 
          legend = "top")

#######
# Proper
# Constant trend
load("diagnostics_sc2.RData")
model_on_short_range = proper1_noInt_ADM4

load("diagnostics_sc8.RData")
model_on_long_range = proper1_noInt_ADM4


plt_const_prop <- plot_fitted_temporal_trends(model_on_short_range,
                                             model_on_long_range,
                                             Improper = F, 
                                             n_ADM4,
                                             tT,
                                             "Constant temporal trend",
                                             ylim = c(-0.5, 0.5))

# Linear trend
load("diagnostics_sc4.RData")
model_on_short_range = proper1_noInt_ADM4

load("diagnostics_sc10.RData")
model_on_long_range = proper1_noInt_ADM4


plt_lin_prop <- plot_fitted_temporal_trends(model_on_short_range,
                                           model_on_long_range,
                                           Improper = F, 
                                           n_ADM4,
                                           tT,
                                           "Linear temporal trend",
                                           ylim = c(-0.5, 0.5),
                                           ylab = F)

# CP
load("diagnostics_sc6.RData")
model_on_short_range = proper1_noInt_ADM4

load("diagnostics_sc12.RData")
model_on_long_range = proper1_noInt_ADM4


plt_cp_prop <- plot_fitted_temporal_trends(model_on_short_range,
                                          model_on_long_range,
                                          Improper = F, 
                                          n_ADM4,
                                          tT,
                                          "Change point temporal trend",
                                          ylim = c(-0.5, 0.55),
                                          ylab = F,
                                          y_breaks = c(-0.5, -0.25, 0,
                                                       0.25, 0.5))

# Save as temporal_trends_proper 9 by 3
ggarrange(plt_const_prop, plt_lin_prop, plt_cp_prop,
          ncol = 3, nrow = 1, common.legend = T, 
          legend = "top")


#############
#Spatial effects

plt_spatial_effects <- function(model,
                                Improper = T,
                                admin_map, 
                                title){
  
  n_ADM = nrow(admin_map)
  
  # Create a color scale for the heatmaps
  scale_col = heat.colors(30, rev=TRUE)
  scale = scale_col[seq(3, 30, length.out = 12)]
  
  if(Improper){
    # Get scale 
    scale_factor <- sqrt(model$summary.hyperpar$mean[4]) * 
      1/sqrt(model$summary.hyperpar$mean[3])
    
    # Get the mean effect
    mean_effect <- model$summary.random$area_id[(n_ADM + 1):(2 * n_ADM), 2] *
      scale_factor
    
    # Get sd of effect
    sd_effect <- model$summary.random$area_id[(n_ADM + 1):(2 * n_ADM), 3] *
      scale_factor
    
  } else{
    # Get the mean effect
    mean_effect <- model$summary.random$area_id[1:n_ADM, 2]
    
    # Get sd of effect
    sd_effect <- model$summary.random$area_id[1:n_ADM, 3]
  } 
  
  
  
  
  
  # Get hardcoded bins
  hardcoded_bins_effect = round(seq(min(mean_effect), max(mean_effect),
                              length.out = 8), 2)
  
  hardcoded_bins_sd = round(seq(min(sd_effect), max(sd_effect),
                                    length.out = 8), 2)
  
    
  
  
  # Actually plot
  admin_map$effect = mean_effect
  admin_map$sd = sd_effect
  
  plt_effect <- heatmap_areas(admin_map,
                              admin_map$effect,
                              scale_col = scale_col,
                              scale = scale,
                              hardcoded_bins = hardcoded_bins_effect,
                              title = TeX(r'(Mean posterior $\theta_{i}$)'))
  
  plt_sd <- heatmap_areas(admin_map,
                            admin_map$sd,
                            scale_col = scale_col,
                            scale = scale,
                            hardcoded_bins = hardcoded_bins_sd,
                            title = TeX(r'(SD posterior $\theta_{i}$)'))
  
    return(
    ggarrange(plt_effect,
            plt_sd,
            ncol = 2, nrow = 1,
            common.legend = F)
    )
}

### Improper
load("diagnostics_sc2.RData")
model_on_short_range <- Improper1_noInt_ADM4

load("diagnostics_sc8.RData")
model_on_long_range = Improper1_noInt_ADM4

plt <- plt_spatial_effects(model_on_short_range,
                            Improper = T,
                            second_level_admin_map,
                            title = "")

# Save as spatial_effect_Imp_short 5 by 3
annotate_figure(plt,  top = text_grob("Model: Improper1_noInt", 
                                      color = "black",  face = "bold", size = 14))

plt <- plt_spatial_effects(model_on_long_range,
                           Improper = T,
                           second_level_admin_map,
                           title = "")

# Save as spatial_effect_Imp_long 5 by 3
annotate_figure(plt,  top = text_grob("Model: Improper1_noInt", 
                                      color = "black",  face = "bold", size = 14))


### proper
load("diagnostics_sc2.RData")
model_on_short_range <- proper1_noInt_ADM4

load("diagnostics_sc8.RData")
model_on_long_range = proper1_noInt_ADM4

plt <- plt_spatial_effects(model_on_short_range,
                           Improper = F,
                           second_level_admin_map,
                           title = "")

# Save as spatial_effect_prop_short 5 by 3
annotate_figure(plt,  top = text_grob("Model: proper1_noInt", 
                                      color = "black",  face = "bold", size = 14))

plt <- plt_spatial_effects(model_on_long_range,
                           Improper = F,
                           second_level_admin_map,
                           title = "")

# Save as spatial_effect_prop_long 5 by 3
annotate_figure(plt,  top = text_grob("Model: proper1_noInt", 
                                      color = "black",  face = "bold", size = 14))




#############
#Space-time interaction in some way???









#########################
#Hyperparameters


##### Temporal

plt_temporal_hyperpar_imp <- function(rw1_model,
                                      rw2_model){
  
  # Transform from precision to standard deviation
  std_temporal_rw1 <-inla.tmarginal(function(x) sqrt(1/x),
                                rw1_model$marginals.hyperpar$`Precision for time_id`)
  
  std_temporal_rw2 <-inla.tmarginal(function(x) sqrt(1/x),
                                    rw2_model$marginals.hyperpar$`Precision for time_id`)
  
  # Format for plotting
  std_temporal_df <- data.frame(x_axis = c(std_temporal_rw1[, 1], std_temporal_rw2[, 1]),
                                y_axis = c(std_temporal_rw1[, 2], std_temporal_rw2[, 2]), 
                                type = c(rep("RW1", length(std_temporal_rw1[, 1])),
                                         rep("RW2", length(std_temporal_rw2[, 1]))))
 
                                
  
  phi_temporal_df <- data.frame(x_axis = c(rw1_model$marginals.hyperpar$`Phi for time_id`[, 1],
                                           rw2_model$marginals.hyperpar$`Phi for time_id`[, 1]),
                                y_axis = c(rw1_model$marginals.hyperpar$`Phi for time_id`[, 2],
                                           rw2_model$marginals.hyperpar$`Phi for time_id`[, 2]),
                                type = c(rep("RW1", length(rw1_model$marginals.hyperpar$`Phi for time_id`[, 1])),
                                         rep("RW2", length(rw2_model$marginals.hyperpar$`Phi for time_id`[, 1]))))
  
  
  
  #Extract how much of the standard deviation is captured by the structured effect
  posterior_samples_rw1 = inla.hyperpar.sample(10000, rw1_model) #Sample the posterior
  posterior_samples_rw2 = inla.hyperpar.sample(10000, rw2_model) #Sample the posterior
  
  std_explained_rw1 = rep(0, 10000); std_explained_rw2 = rep(0, 10000)
  for(i in 1:10000){
    std_explained_rw1[i] = sqrt(posterior_samples_rw1[i, 2] * 
                              as.numeric(1/posterior_samples_rw1[i, 1]))
    std_explained_rw2[i] = sqrt(posterior_samples_rw2[i, 2] * 
                                  as.numeric(1/posterior_samples_rw2[i, 1]))
  }
  
  
  sigma_sqrtPhi.df <- data.frame(x_axis = c(std_explained_rw1, std_explained_rw2),
                                 type = c(rep("RW1", 10000), rep("RW2", 10000)))
  
  
  std_temporal_plot <- ggplot(data=std_temporal_df) + 
    stat_density(aes(x = x_axis, fill = type),
                 adjust = 1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) + 
    xlab(expression(sigma)) + ylab(TeX("$f(*)$")) #expression(f(sigma)) 
  
  
  phi_temporal_plot <- ggplot(data=phi_temporal_df) + 
    stat_density(aes(x=x_axis, fill = type),
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(phi)) + ylab("") #expression(f(lambda))
  
  
  sigma_phi.plot <- ggplot(data=sigma_sqrtPhi.df) + 
    stat_density(aes(x=x_axis, fill = type),
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(TeX("$\\sigma\\sqrt{\\phi}$")) + #"$\\left(\\frac{\\phi}{\\tau}\\right)^{1/2}$"
    ylab("")

  
  plt <- ggarrange(std_temporal_plot, phi_temporal_plot, sigma_phi.plot,
                   ncol = 3, nrow = 1,
                   common.legend = T, legend = "top")
  
  return(
    annotate_figure(plt, top = text_grob("Temporal hyperparameters of improper models", 
                                       color = "black", size = 14))
    )
}

### SC2
load("diagnostics_sc2.RData")

# Save as sc2_hyperparameters_temporal_imp 9 by 3
plt_temporal_hyperpar_imp(Improper1_noInt_ADM4,
                          Improper2_typeI_ADM4)

### SC4
load("diagnostics_sc4.RData")

# Save as sc4_hyperparameters_temporal_imp 9 by 3
plt_temporal_hyperpar_imp(Improper1_noInt_ADM4,
                          Improper2_typeI_ADM4)

### SC6
load("diagnostics_sc6.RData")

# Save as sc6_hyperparameters_temporal_imp 9 by 3
plt_temporal_hyperpar_imp(Improper1_noInt_ADM4,
                          Improper2_typeI_ADM4)

### SC8
load("diagnostics_sc8.RData")

# Save as sc8_hyperparameters_temporal_imp 9 by 3
plt_temporal_hyperpar_imp(Improper1_noInt_ADM4,
                          Improper2_typeI_ADM4)

### SC10
load("diagnostics_sc10.RData")

# Save as sc10_hyperparameters_temporal_imp 9 by 3
plt_temporal_hyperpar_imp(Improper1_noInt_ADM4,
                          Improper2_typeI_ADM4)

### SC12
load("diagnostics_sc12.RData")

# Save as sc12_hyperparameters_temporal_imp 9 by 3
plt_temporal_hyperpar_imp(Improper1_noInt_ADM4,
                          Improper2_typeI_ADM4)




plt_temporal_hyperpar_prop <- function(ar1_model,
                                      ar2_model,
                                      std = T){
  
  if(std){
    std_temporal_ar1 <-inla.tmarginal(function(x) sqrt(1/x),
                                      ar1_model$marginals.hyperpar$`Precision for time_id.copy`)
    
    std_temporal_ar2 <-inla.tmarginal(function(x) sqrt(1/x),
                                      ar2_model$marginals.hyperpar$`Precision for time_id.copy`)
    xlab = expression(sigma)
    
  } else { # If problems
    std_temporal_ar1 <- ar1_model$marginals.hyperpar$`Precision for time_id.copy`
    
    std_temporal_ar2 <- ar2_model$marginals.hyperpar$`Precision for time_id.copy`
    xlab = expression(tau)
  }
  
  
  # Transform from precision to standard deviation
  
  
  
  # Format for plotting
  std_temporal_df <- data.frame(x_axis = c(std_temporal_ar1[, 1], std_temporal_ar2[, 1]),
                                y_axis = c(std_temporal_ar1[, 2], std_temporal_ar2[, 2]), 
                                type = c(rep("AR1", length(std_temporal_ar1[, 1])),
                                         rep("AR2", length(std_temporal_ar2[, 1]))))
  
  std_temporal_ar1_df <- data.frame(x_axis = std_temporal_ar1[, 1],
                                y_axis = std_temporal_ar1[, 2])
  
  std_temporal_ar2_df <- data.frame(x_axis = std_temporal_ar2[, 1],
                                    y_axis = std_temporal_ar2[, 2])
  
  
  rho_temporal_df <- data.frame(x_axis = c(ar1_model$marginals.hyperpar$`Rho for time_id.copy`[, 1],
                                           ar2_model$marginals.hyperpar$`PACF1 for time_id.copy`[, 1],
                                           ar2_model$marginals.hyperpar$`PACF2 for time_id.copy`[, 1]),
                                y_axis = c(ar1_model$marginals.hyperpar$`Rho for time_id.copy`[, 2],
                                           ar2_model$marginals.hyperpar$`PACF1 for time_id.copy`[, 2],
                                           ar2_model$marginals.hyperpar$`PACF2 for time_id.copy`[, 2]),
                                type = c(rep("AR1: rho", length(ar1_model$marginals.hyperpar$`Rho for time_id.copy`[, 1])),
                                         rep("AR2: rho 1", length(ar2_model$marginals.hyperpar$`PACF1 for time_id.copy`[, 1])),
                                         rep("AR2: rho 2", length(ar2_model$marginals.hyperpar$`PACF2 for time_id.copy`[, 1]))))
  
  
  
  std_temporal_plot <- ggplot(data=std_temporal_df) + 
    stat_density(aes(x = x_axis, fill = type),
                 adjust = 1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) + 
    xlab(xlab) + ylab(TeX("$f(*)$")) +  #expression(f(sigma)) 
    guides(fill=guide_legend(title=NULL, 
                             label.position = "top",
                             nrow = 1))
  
 
  
  rho_temporal_plot <- ggplot(data=rho_temporal_df) + 
    stat_density(aes(x=x_axis, fill = type),
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho)) + ylab("") + #expression(f(lambda))
    guides(fill=guide_legend(title=NULL, 
                             label.position = "top",
                             nrow = 1))
  rho_temporal_plot <- rho_temporal_plot + scale_fill_manual(values = c("#F8766D",
                                                                         "#00BFC4", 
                                                                         "#33cc33"),
                       labels = unname(TeX(c(" AR1: $\\rho$",
                                             " AR2: $\\rho_{1}$",
                                             " AR2: $\\rho_{2}$"))))
  
  
  plt <- ggarrange(std_temporal_plot, rho_temporal_plot,
                   ncol = 2, nrow = 1,
                   common.legend = F,
                   legend = "top")
  
  return(
    annotate_figure(plt, top = text_grob("Temporal hyperparameters of proper models", 
                                         color = "black", size = 14))
  )
}

### SC2
load("diagnostics_sc2.RData")

# Save as sc2_hyperparameters_temporal_prop 6 by 3
plt_temporal_hyperpar_prop(proper1_full_ADM4,
                           proper2_full_ADM4,
                           std = F)

### SC4
load("diagnostics_sc4.RData")

# Save as sc4_hyperparameters_temporal_prop 6 by 3
plt_temporal_hyperpar_prop(proper1_noInt_ADM4,
                           proper2_noInt_ADM4,
                           std = T)


### SC6
load("diagnostics_sc6.RData")

# Save as sc6_hyperparameters_temporal_prop 6 by 3
plt_temporal_hyperpar_prop(proper1_noInt_ADM4,
                           proper2_noInt_ADM4,
                           std = F)

### SC8
load("diagnostics_sc8.RData")

# Save as sc8_hyperparameters_temporal_prop 6 by 3
plt_temporal_hyperpar_prop(proper1_full_ADM4,
                           proper2_full_ADM4,
                           std = F)

### SC10
load("diagnostics_sc10.RData")

# Save as sc10_hyperparameters_temporal_prop 6 by 3
plt_temporal_hyperpar_prop(proper1_noInt_ADM4,
                           proper2_noInt_ADM4,
                           std = T)


### SC12
load("diagnostics_sc12.RData")

# Save as sc12_hyperparameters_temporal_prop 6 by 3
plt_temporal_hyperpar_prop(proper1_noInt_ADM4,
                           proper2_noInt_ADM4,
                           std = F)



##### Spatial

plt_spatial_hyperpar_imp <- function(model){
  
  # Transform from precision to standard deviation
  std_spatial <-inla.tmarginal(function(x) sqrt(1/x),
                                model$marginals.hyperpar$`Precision for area_id`)
  
  # Format for plotting
  std_spatial_df <- data.frame(x_axis = std_spatial[, 1],
                                y_axis = std_spatial[, 2])
  
  
  
  phi_spatial_df <- data.frame(x_axis = model$marginals.hyperpar$`Phi for area_id`[, 1], 
                                y_axis = model$marginals.hyperpar$`Phi for area_id`[, 2])
  
  #Extract how much of the standard deviation is captured by the structured effect
  posterior_samples = inla.hyperpar.sample(10000, model) #Sample the posterior
  
  #Precision for county, Phi for county col 3 and 4 respectively
  std_explained = rep(0, 10000)
  for(i in 1:10000){
    std_explained[i] = sqrt(posterior_samples[i, 4] * as.numeric(1/posterior_samples[i, 3]))
  }
  
  sigma_sqrtPhi.df <- data.frame(x_axis = std_explained)
  
  std_temporal_plot <- ggplot(data=std_spatial_df) + 
    stat_density(aes(x = x_axis), fill = "#F8766D",
                 adjust = 1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) + 
    xlab(expression(sigma)) + ylab(TeX("$f(*)$")) #expression(f(sigma)) 
  
  
  phi_temporal_plot <- ggplot(data=phi_spatial_df) + 
    stat_density(aes(x=x_axis), fill = "#F8766D",
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(phi)) + ylab("") #expression(f(lambda))
  
  
  sigma_phi.plot <- ggplot(data=sigma_sqrtPhi.df) + 
    stat_density(aes(x=x_axis), fill = "#F8766D",
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(TeX("$\\sigma\\sqrt{\\phi}$")) + #"$\\left(\\frac{\\phi}{\\tau}\\right)^{1/2}$"
    ylab("")
  
  
  plt <- ggarrange(std_temporal_plot, phi_temporal_plot, sigma_phi.plot,
                   ncol = 3, nrow = 1,
                   common.legend = T, legend = "top")
  
  return(
    annotate_figure(plt, top = text_grob("Spatial hyperparameters of improper models", 
                                         color = "black", size = 14))
  )
}

### SC2
load("diagnostics_sc2.RData")

# Save as sc2_hyperparameters_spatial_imp 9 by 3
plt_spatial_hyperpar_imp(Improper1_noInt_ADM4)

### SC4
load("diagnostics_sc4.RData")

# Save as sc4_hyperparameters_spatial_imp 9 by 3
plt_spatial_hyperpar_imp(Improper1_noInt_ADM4)

### SC6
load("diagnostics_sc6.RData")

# Save as sc6_hyperparameters_spatial_imp 9 by 3
plt_spatial_hyperpar_imp(Improper1_noInt_ADM4)

### SC8
load("diagnostics_sc8.RData")

# Save as sc8_hyperparameters_spatial_imp 9 by 3
plt_spatial_hyperpar_imp(Improper1_noInt_ADM4)

### SC10
load("diagnostics_sc10.RData")

# Save as sc10_hyperparameters_spatial_imp 9 by 3
plt_spatial_hyperpar_imp(Improper1_noInt_ADM4)

### SC12
load("diagnostics_sc12.RData")

# Save as sc12_hyperparameters_spatial_imp 9 by 3
plt_spatial_hyperpar_imp(Improper1_noInt_ADM4)



plt_spatial_hyperpar_prop <- function(model){
  
  
  # Transform from precision to standard deviation
  std_spatial <-inla.tmarginal(function(x) sqrt(1/x),
                               model$marginals.hyperpar$`Precision for area_id`)
  
  # Format for plotting
  std_spatial_df <- data.frame(x_axis = std_spatial[, 1],
                               y_axis = std_spatial[, 2])
  
  
  rho_spatial_df <- data.frame(x_axis = model$marginals.hyperpar$`Lambda for area_id`[, 1], 
                               y_axis = model$marginals.hyperpar$`Lambda for area_id`[, 2])
  
  
  std_temporal_plot <- ggplot(data=std_spatial_df) + 
    stat_density(aes(x = x_axis), fill = "#F8766D",
                 adjust = 1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) + 
    xlab(expression(sigma)) + ylab(TeX("$f(*)$")) #expression(f(sigma)) 
  
  
  rho_spatial_plot <- ggplot(data=rho_spatial_df) + 
    stat_density(aes(x=x_axis), fill = "#F8766D",
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho)) + ylab("") #expression(f(lambda))
  
  
  
  plt <- ggarrange(std_temporal_plot, rho_spatial_plot,
                   ncol = 2, nrow = 1,
                   common.legend = T, legend = "top")
  
  return(
    annotate_figure(plt, top = text_grob("Spatial hyperparameters of proper models", 
                                         color = "black", size = 14))
  )
}


### SC2
load("diagnostics_sc2.RData")

# Save as sc2_hyperparameters_spatial_prop 6 by 3
plt_spatial_hyperpar_prop(proper1_noInt_ADM4)

### SC4
load("diagnostics_sc4.RData")

# Save as sc4_hyperparameters_spatial_prop 6 by 3
plt_spatial_hyperpar_prop(proper1_noInt_ADM4)

### SC6
load("diagnostics_sc6.RData")

# Save as sc6_hyperparameters_spatial_prop 6 by 3
plt_spatial_hyperpar_prop(proper1_noInt_ADM4)

### SC8
load("diagnostics_sc8.RData")

# Save as sc8_hyperparameters_spatial_prop 6 by 3
plt_spatial_hyperpar_prop(proper1_noInt_ADM4)

### SC10
load("diagnostics_sc10.RData")

# Save as sc10_hyperparameters_spatial_prop 6 by 3
plt_spatial_hyperpar_prop(proper1_noInt_ADM4)

### SC12
load("diagnostics_sc12.RData")

# Save as sc12_hyperparameters_spatial_prop 6 by 3
plt_spatial_hyperpar_prop(proper1_noInt_ADM4)






##### Interaction

plt_interaction_hyperpar_imp <- function(typeI,
                                         typeII, 
                                         typeIII,
                                         typeIV,
                                         title = ""){
  
  # Transform from precision to standard deviation
  std_int_I <-inla.tmarginal(function(x) sqrt(1/x),
                           typeI$marginals.hyperpar$`Precision for space.time`)
  
  if(title == "Const, long"){
    std_int_II = NULL
  } else {
    std_int_II <-inla.tmarginal(function(x) sqrt(1/x),
                                typeII$marginals.hyperpar$`Precision for space.time`)
  }
  
 
  std_int_III <-inla.tmarginal(function(x) sqrt(1/x),
                             typeIII$marginals.hyperpar$`Precision for space.time`)
  
  std_int_IV <-inla.tmarginal(function(x) sqrt(1/x),
                             typeIV$marginals.hyperpar$`Precision for space.time`)
  
  if(title == "Const, long"){
    std_df <- data.frame(x_axis = c(std_int_I[, 1],
                                    std_int_III[, 1],
                                    std_int_IV[, 1]),
                         type = c(rep("type I", length(std_int_I[, 1])),
                                  rep("type III", length(std_int_III[, 1])),
                                  rep("type IV", length(std_int_IV[, 1]))))
    
  } else {
    std_df <- data.frame(x_axis = c(std_int_I[, 1],
                                    std_int_II[, 1],
                                    std_int_III[, 1],
                                    std_int_IV[, 1]),
                         type = c(rep("type I", length(std_int_I[, 1])),
                                  rep("type II", length(std_int_II[, 1])),
                                  rep("type III", length(std_int_III[, 1])),
                                  rep("type IV", length(std_int_IV[, 1]))))
    
  }
  
  if(title == "Const, long"){
    std_plot <- ggplot(data=std_df) + ggtitle(title) +
      stat_density(aes(x = x_axis, fill = type),
                   adjust = 1.5, alpha=.8, position = "identity") +
      theme_bw() +
      theme(axis.title=element_text(size=14)) +
      labs(fill = NULL) + 
      xlab(expression(sigma)) + ylab(TeX("$f(\\sigma)$")) + #expression(f(sigma)) 
      scale_fill_manual(values=c("#F8766D", "#00BFC4", "#C77CFF"))
  } else {
    std_plot <- ggplot(data=std_df) + ggtitle(title) +
      stat_density(aes(x = x_axis, fill = type),
                   adjust = 1.5, alpha=.8, position = "identity") +
      theme_bw() +
      theme(axis.title=element_text(size=14)) +
      labs(fill = NULL) + 
      xlab(expression(sigma)) + ylab(TeX("$f(\\sigma)$")) #expression(f(sigma)) 
    
  }
  
  
  
  
   return(std_plot)
}

### SC2
load("diagnostics_sc2.RData")

sc2_int_hp_imp <- plt_interaction_hyperpar_imp(Improper2_typeI_ADM4,
                             Improper1_typeII_ADM4,
                             Improper2_typeIII_ADM4,
                             Improper1_typeIV_ADM4,
                             title = "Const, short")


### SC4
load("diagnostics_sc4.RData")

sc4_int_hp_imp <- plt_interaction_hyperpar_imp(Improper2_typeI_ADM4,
                                               Improper1_typeII_ADM4,
                                               Improper2_typeIII_ADM4,
                                               Improper1_typeIV_ADM4,
                                               title = "Lin, short")

### SC6
load("diagnostics_sc6.RData")

sc6_int_hp_imp <- plt_interaction_hyperpar_imp(Improper2_typeI_ADM4,
                                               Improper1_typeII_ADM4,
                                               Improper2_typeIII_ADM4,
                                               Improper1_typeIV_ADM4,
                                               title = "CP, short")



plt <- ggarrange(sc2_int_hp_imp, sc4_int_hp_imp, sc6_int_hp_imp,
          ncol = 3, nrow = 1, 
          common.legend = T,
          legend = "top")

# Save as interactions_hyperpars_short_imp 9 by 3
annotate_figure(plt, top = text_grob("Marginal SD of interactions of improper models", 
                                     color = "black", size = 14))

### SC8
load("diagnostics_sc8.RData")

sc8_int_hp_imp <- plt_interaction_hyperpar_imp(Improper2_typeI_ADM4,
                                               Improper2_typeII_ADM4,
                                               Improper2_typeIII_ADM4,
                                               Improper1_typeIV_ADM4,
                                               title = "Const, long")


### SC10
load("diagnostics_sc10.RData")

sc10_int_hp_imp <- plt_interaction_hyperpar_imp(Improper2_typeI_ADM4,
                                               Improper1_typeII_ADM4,
                                               Improper2_typeIII_ADM4,
                                               Improper1_typeIV_ADM4,
                                               title = "Lin, long")

### SC12
load("diagnostics_sc12.RData")

sc12_int_hp_imp <- plt_interaction_hyperpar_imp(Improper2_typeI_ADM4,
                                               Improper1_typeII_ADM4,
                                               Improper2_typeIII_ADM4,
                                               Improper1_typeIV_ADM4,
                                               title = "CP, long")



plt <- ggarrange(sc10_int_hp_imp, sc12_int_hp_imp, sc8_int_hp_imp,
                 ncol = 3, nrow = 1, 
                 common.legend = T,
                 legend = "top")


# Save as interactions_hyperpars_long_imp 9 by 3
annotate_figure(plt, top = text_grob("Marginal SD of interactions of improper models", 
                                     color = "black", size = 14))



plt_interaction_hyperpar_prop <- function(prop1_onlyInt,
                                          prop2_onlyInt,
                                          prop1_full,
                                          prop2_full,
                                         title = ""){
  
  # Transform from precision to standard deviation
  std_int_1_onlyInt <-inla.tmarginal(function(x) sqrt(1/x),
                             prop1_onlyInt$marginals.hyperpar$`Precision for area_id`)
  
  std_int_2_onlyInt <-inla.tmarginal(function(x) sqrt(1/x),
                              prop2_onlyInt$marginals.hyperpar$`Precision for area_id`)
  
  std_int_1_full <-inla.tmarginal(function(x) sqrt(1/x),
                               prop1_full$marginals.hyperpar$`Precision for area_id.copy`)
  
  std_int_2_full <-inla.tmarginal(function(x) sqrt(1/x),
                              prop2_full$marginals.hyperpar$`Precision for area_id.copy`)
  
  
  std_df <- data.frame(x_axis = c(std_int_1_onlyInt[, 1],
                                  std_int_2_onlyInt[, 1],
                                  std_int_1_full[, 1],
                                  std_int_2_full[, 1]),
                       type = c(rep("Proper1_onlyInt", length(std_int_1_onlyInt[, 1])),
                                rep("Proper2_onlyInt", length(std_int_2_onlyInt[, 1])),
                                rep("Proper1_full", length(std_int_1_full[, 1])),
                                rep("Proper2_full", length(std_int_2_full[, 1]))))
  
  
  
  rho_spatial_df <- data.frame(x_axis = c(prop1_onlyInt$marginals.hyperpar$`Lambda for area_id`[, 1],
                                          prop2_onlyInt$marginals.hyperpar$`Lambda for area_id`[, 1],
                                          prop1_full$marginals.hyperpar$`Lambda for area_id.copy`[, 1], 
                                          prop2_full$marginals.hyperpar$`Lambda for area_id.copy`[, 1]),
                               type = c(rep("Proper1_onlyInt", length(prop1_onlyInt$marginals.hyperpar$`Lambda for area_id`[, 1])),
                                        rep("Proper2_onlyInt", length(prop2_onlyInt$marginals.hyperpar$`Lambda for area_id`[, 1])),
                                        rep("Proper1_full", length(prop1_full$marginals.hyperpar$`Lambda for area_id.copy`[, 1])),
                                        rep("Proper2_full", length(prop2_full$marginals.hyperpar$`Lambda for area_id.copy`[, 1]))))
  
  
  rho1_temporal_df <- data.frame(x_axis = c(prop1_onlyInt$marginals.hyperpar$`GroupRho for area_id`[, 1],
                                          prop2_onlyInt$marginals.hyperpar$`Group PACF1 for area_id`[, 1],
                                          prop1_full$marginals.hyperpar$`GroupRho for area_id.copy`[, 1], 
                                          prop2_full$marginals.hyperpar$`Group PACF1 for area_id.copy`[, 1]),
                               type = c(rep("Proper1_onlyInt", length(prop1_onlyInt$marginals.hyperpar$`GroupRho for area_id`[, 1])),
                                        rep("Proper2_onlyInt", length(prop2_onlyInt$marginals.hyperpar$`Group PACF1 for area_id`[, 1])),
                                        rep("Proper1_full", length(prop1_full$marginals.hyperpar$`GroupRho for area_id.copy`[, 1])),
                                        rep("Proper2_full", length(prop2_full$marginals.hyperpar$`Group PACF1 for area_id.copy`[, 1]))))
  
  rho2_temporal_df <- data.frame(x_axis = c(prop2_onlyInt$marginals.hyperpar$`Group PACF2 for area_id`[, 1],
                                            prop2_full$marginals.hyperpar$`Group PACF2 for area_id.copy`[, 1]),
                                 type = c(rep("Proper2_onlyInt", length(prop2_onlyInt$marginals.hyperpar$`Group PACF2 for area_id`[, 1])),
                                          rep("Proper2_full", length(prop2_full$marginals.hyperpar$`Group PACF2 for area_id.copy`[, 1]))))
  
  
  
  
  
  std_plot <- ggplot(data=std_df) + ggtitle(title) +
    stat_density(aes(x = x_axis, fill = type),
                 adjust = 1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) + 
    xlab(expression(sigma)) + ylab(TeX("$f(\\sigma)$")) #expression(f(sigma)) 
  
  
  rho_spatial_plot <- ggplot(data=rho_spatial_df) + 
    stat_density(aes(x=x_axis, fill = type), 
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho[s])) + ylab("") #expression(f(lambda))
  
  rho1_temporal_plot <- ggplot(data=rho1_temporal_df) + 
    stat_density(aes(x=x_axis, fill = type), 
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho[1])) + ylab("") #expression(f(lambda))
  
  rho2_temporal_plot <- ggplot(data=rho2_temporal_df) + 
    stat_density(aes(x=x_axis, fill = type), 
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho[2])) + ylab("") +  #expression(f(lambda))
    scale_fill_manual(values=c("#00BFC4", "#C77CFF"))
  
  return(ggarrange(std_plot, rho_spatial_plot,
                   rho1_temporal_plot, rho2_temporal_plot,
                   ncol = 2, nrow = 2,
                   common.legend = T,
                   legend = "top"))
}

### SC2
load("diagnostics_sc2.RData")

# Save as sc2_proper_interactions_hyperparams 6 by 6
plt_interaction_hyperpar_prop(proper1_onlyInt_ADM4,
                              proper2_onlyInt_ADM4,
                              proper1_full_ADM4,
                              proper2_full_ADM4)

### SC4
load("diagnostics_sc4.RData")

# Save as sc4_proper_interactions_hyperparams 6 by 6
plt_interaction_hyperpar_prop(proper1_onlyInt_ADM4,
                              proper2_onlyInt_ADM4,
                              proper1_full_ADM4,
                              proper2_full_ADM4)

### SC6
load("diagnostics_sc6.RData")

# Save as sc6_proper_interactions_hyperparams 6 by 6
plt_interaction_hyperpar_prop(proper1_onlyInt_ADM4,
                              proper2_onlyInt_ADM4,
                              proper1_full_ADM4,
                              proper2_full_ADM4)



################################################################################
# Make some figure for the distribution of the widths of the 95% CIs for the models
## In order to see how the width of the CIs change


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
                        hardcoded_bins = hardcoded_bins_mu, title = "year 3")
  
  #### Plot year 7
  tmp_map_$mu <- lambda_[lambda_$time_id == 7, ]$mu
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "year 7")
  
  #### Plot year 11
  tmp_map_$mu <- lambda_[lambda_$time_id == 11, ]$mu
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "year 11")
  
  #### Plot year 13
  tmp_map_$mu <- lambda_[lambda_$time_id == 13, ]$mu
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$mu,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_mu, title = "year 13")
  return(ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T,
                   legend = "right"))
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
                        hardcoded_bins = hardcoded_bins_count, title = "year 3")
  
  
  #### Plot year 7
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 7, ]$sampled_counts
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "year 7")
  
  #### Plot year 11
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 11, ]$sampled_counts
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "year 11")
  
  #### Plot year 13
  tmp_map_$sampled_counts <- lambda_[lambda_$time_id == 13, ]$sampled_counts
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sampled_counts,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins_count, title = "year 13")
  
  return(ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T,
                   legend = "right"))
}


scale_col = heat.colors(30, rev=TRUE)
scale = scale_col[seq(3, 30, length.out = 12)]

#temp
plt_true_discrete_rate_four_years(scenario_name = "sc13", 
                                  dataset_id = 1,
                                  admin_map = new_map,
                                  scale_col = scale_col,
                                  scale = scale)


plt_true_discrete_rate_four_years(scenario_name = "sc14", 
                                  dataset_id = 1,
                                  admin_map = new_map,
                                  scale_col = scale_col,
                                  scale = scale)


#### SC2
plt5 <- plt_true_discrete_rate_four_years(scenario_name = "sc2", 
                                          dataset_id = dataset_id_2,
                                          admin_map = second_level_admin_map,
                                          scale_col = scale_col,
                                          scale = scale)

# Save as sc2_true_rate 8.5 by 3
annotate_figure(plt5, 
                top = text_grob(TeX(r'(Scenario: ADM4$_{const, short}$, Simulated rate per 100 for )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))




plt5 <- plt_true_counts_four_years(scenario_name = "sc2", 
                                          dataset_id = dataset_id_2,
                                          admin_map = second_level_admin_map,
                                          scale_col = scale_col,
                                          scale = scale)

# Save as sc2_true_count 8.5 by 3
annotate_figure(plt5,
                top = text_grob(TeX(r'(Scenario: ADM4$_{const, short}$, Simulated counts for )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))


#### SC4
plt5 <- plt_true_discrete_rate_four_years(scenario_name = "sc4", 
                                          dataset_id = dataset_id_2,
                                          admin_map = second_level_admin_map,
                                          scale_col = scale_col,
                                          scale = scale)

# Save as sc4_true_rate 8.5 by 3
annotate_figure(plt5, 
                top = text_grob(TeX(r'(Scenario: ADM4$_{lin, short}$, Simulated rate per 100 for )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))




plt5 <- plt_true_counts_four_years(scenario_name = "sc4", 
                                   dataset_id = dataset_id_2,
                                   admin_map = second_level_admin_map,
                                   scale_col = scale_col,
                                   scale = scale)

# Save as sc4_true_count 8.5 by 3
annotate_figure(plt5,
                top = text_grob(TeX(r'(Scenario: ADM4$_{lin, short}$, Simulated counts for )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))


#### SC6
plt5 <- plt_true_discrete_rate_four_years(scenario_name = "sc6", 
                                          dataset_id = dataset_id_2,
                                          admin_map = second_level_admin_map,
                                          scale_col = scale_col,
                                          scale = scale)

# Save as sc6_true_rate 8.5 by 3
annotate_figure(plt5, 
                top = text_grob(TeX(r'(Scenario: ADM4$_{cp, short}$, Simulated rate per 100 for )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))




plt5 <- plt_true_counts_four_years(scenario_name = "sc6", 
                                   dataset_id = dataset_id_2,
                                   admin_map = second_level_admin_map,
                                   scale_col = scale_col,
                                   scale = scale)

# Save as sc6_true_count 8.5 by 3
annotate_figure(plt5,
                top = text_grob(TeX(r'(Scenario: ADM4$_{cp, short}$, Simulated counts for )'), 
                                color = "black", 
                                face = "bold", 
                                size = 14))




#####################

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
                        hardcoded_bins = hardcoded_bins, title = "year 3")
  
  #### Plot year 7
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (7 - 1) + 1):(n_ADM * 7)] * 100
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "year 7")
  
  #### Plot year 11
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (11 - 1) + 1):(n_ADM * 11)] * 100
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "year 11")
  
  #### Plot year 13
  tmp_map_$pred_rate <- model$summary.fitted.values$mean[(n_ADM * (13 - 1) + 1):(n_ADM * 13)] * 100
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$pred_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "year 13")
  
  plt <- ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T, 
                   legend = "right")
  
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
                        hardcoded_bins = hardcoded_bins, title = "year 3")
  
  #### Plot year 7
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (7 - 1) + 1):(n_ADM * 7)] * 100
  plt2 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "year 7")
  
  #### Plot year 11
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (11 - 1) + 1):(n_ADM * 11)] * 100
  plt3 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "year 11")
  
  #### Plot year 13
  tmp_map_$sd_rate <- model$summary.fitted.values$sd[(n_ADM * (13 - 1) + 1):(n_ADM * 13)] * 100
  plt4 <- heatmap_areas(map_w_values = tmp_map_, value = tmp_map_$sd_rate,
                        scale_col = scale_col, scale = scale,
                        hardcoded_bins = hardcoded_bins, title = "year 13")
  
  plt <- ggarrange(plt1, plt2, plt3, plt4, 
                   ncol = 4, nrow = 1, 
                   common.legend = T, 
                   legend = "right")
  
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
# Save as sc2_mean_fitted_rate_Imp 8.5 by 3
plt_fitted_rate_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, fitted rate per 100 of Improper1_typeIV for )'))

## Plot the sd of the rate
# Save as sc2_sd_fitted_rate_Imp 8.5 by 3
plt_fitted_rate_sd_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, sd of rate per 100 of Improper1_typeIV for )'))



### Proper

## Plot the mean rate
# Save as sc2_mean_fitted_rate_prop 8.5 by 3
plt_fitted_rate_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, fitted rate per 100 of proper2_onlyInt for )'))

## Plot the sd of the rate
# Save as sc2_sd_fitted_rate_prop 8.5 by 3
plt_fitted_rate_sd_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{const, short}$, sd of rate per 100 of proper2_onlyInt for )'))




##### 
#SC4

scenario_name = "sc4"
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
# Save as sc4_mean_fitted_rate_Imp 8.5 by 3
plt_fitted_rate_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{lin, short}$, fitted rate per 100 of Improper1_typeIV for )'))

## Plot the sd of the rate
# Save as sc4_sd_fitted_rate_Imp 8.5 by 3
plt_fitted_rate_sd_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                              scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                              overall_title = TeX(r'(Scenario: ADM4$_{lin, short}$, sd of rate per 100 of Improper1_typeIV for )'))



### Proper

## Plot the mean rate
# Save as sc4_mean_fitted_rate_prop 8.5 by 3
plt_fitted_rate_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{lin, short}$, fitted rate per 100 of proper2_onlyInt for )'))

## Plot the sd of the rate
# Save as sc4_sd_fitted_rate_prop 8.5 by 3
plt_fitted_rate_sd_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                              scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                              overall_title = TeX(r'(Scenario: ADM4$_{lin, short}$, sd of rate per 100 of proper2_onlyInt for )'))


##### 
#SC6

scenario_name = "sc6"
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
# Save as sc6_mean_fitted_rate_Imp 8.5 by 3
plt_fitted_rate_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{cp, short}$, fitted rate per 100 of Improper1_typeIV for )'))

## Plot the sd of the rate
# Save as sc6_sd_fitted_rate_Imp 8.5 by 3
plt_fitted_rate_sd_four_years(model = best_imp_model, Improper = T, admin_map = second_level_admin_map,
                              scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                              overall_title = TeX(r'(Scenario: ADM4$_{cp, short}$, sd of rate per 100 of Improper1_typeIV for )'))



### Proper

## Plot the mean rate
# Save as sc6_mean_fitted_rate_prop 8.5 by 3
plt_fitted_rate_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                           scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_mean_rate,
                           overall_title = TeX(r'(Scenario: ADM4$_{cp, short}$, fitted rate per 100 of proper2_onlyInt for )'))

## Plot the sd of the rate
# Save as sc6_sd_fitted_rate_prop 8.5 by 3
plt_fitted_rate_sd_four_years(model = best_prop_model, Improper = F, admin_map = second_level_admin_map,
                              scale_col = scale_col, scale = scale, hardcoded_bins = hardcoded_bins_sd_rate,
                              overall_title = TeX(r'(Scenario: ADM4$_{cp, short}$, sd of rate per 100 of proper2_onlyInt for )'))



