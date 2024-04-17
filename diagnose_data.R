
# Look at how the temporal trend of the data is across all regions
# i.e. plot all together or average over regions for one year
# and see if it has the temporal trend it is supposed to have


# Other diagnoses of the data

#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

# Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)

dataset_id = 3

################################################################################

#Load in data set

#Load in simulation data (done stupidly, change!)
load("./Simulated_data/sc1/sc1_data.RData")
lambda_sc1.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc1.df$mu = lambda.df$mu[, dataset_id]
lambda_sc1.df$lambda_it = lambda.df$lambda_it[, dataset_id]


load("./Simulated_data/sc2/sc2_data.RData")
lambda_sc2.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc2.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc2.df$mu = lambda.df$mu[, dataset_id]
lambda_sc2.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc3/sc3_data.RData")
lambda_sc3.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc3.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc3.df$mu = lambda.df$mu[, dataset_id]
lambda_sc3.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc4/sc4_data.RData")
lambda_sc4.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc4.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc4.df$mu = lambda.df$mu[, dataset_id]
lambda_sc4.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc5/sc5_data.RData")
lambda_sc5.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc5.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc5.df$mu = lambda.df$mu[, dataset_id]
lambda_sc5.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc6/sc6_data.RData")
lambda_sc6.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc6.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc6.df$mu = lambda.df$mu[, dataset_id]
lambda_sc6.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc7/sc7_data.RData")
lambda_sc7.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc7.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc7.df$mu = lambda.df$mu[, dataset_id]
lambda_sc7.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc8/sc8_data.RData")
lambda_sc8.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc8.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc8.df$mu = lambda.df$mu[, dataset_id]
lambda_sc8.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc9/sc9_data.RData")
lambda_sc9.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc9.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc9.df$mu = lambda.df$mu[, dataset_id]
lambda_sc9.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc10/sc10_data.RData")
lambda_sc10.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc10.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc10.df$mu = lambda.df$mu[, dataset_id]
lambda_sc10.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc11/sc11_data.RData")
lambda_sc11.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc11.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc11.df$mu = lambda.df$mu[, dataset_id]
lambda_sc11.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc12/sc12_data.RData")
lambda_sc12.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc12.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc12.df$mu = lambda.df$mu[, dataset_id]
lambda_sc12.df$lambda_it = lambda.df$lambda_it[, dataset_id]

rm(lambda.df)

################################################################################
# Produce diagnostics plots




plt1 <- plot_temporal_trend_data_one_data_set(lambda_sc1.df, "ADM1 const, short: one data set")
plt2 <- plot_temporal_trend_data_one_data_set(lambda_sc2.df, "ADM4 const, short: one data set")
plt3 <-plot_temporal_trend_data_one_data_set(lambda_sc3.df, "ADM1 lin increasing, short: one data set")
plt4 <-plot_temporal_trend_data_one_data_set(lambda_sc4.df, "ADM4 lin increasing, short: one data set")
plt5 <-plot_temporal_trend_data_one_data_set(lambda_sc5.df, "ADM1 change point, short: one data set")
plt6 <-plot_temporal_trend_data_one_data_set(lambda_sc6.df, "ADM4 change point, short: one data set")
plt7 <-plot_temporal_trend_data_one_data_set(lambda_sc7.df, "ADM1 const, long: one data set")
plt8 <-plot_temporal_trend_data_one_data_set(lambda_sc8.df, "ADM4 const, long: one data set")
plt9 <-plot_temporal_trend_data_one_data_set(lambda_sc9.df, "ADM1 lin increasing, long: one data set")
plt10 <-plot_temporal_trend_data_one_data_set(lambda_sc10.df, "ADM4 lin increasing, long: one data set")
plt11 <-plot_temporal_trend_data_one_data_set(lambda_sc11.df, "ADM1 change point, long: one data set")
plt12 <-plot_temporal_trend_data_one_data_set(lambda_sc12.df, "ADM4 change point, long: one data set")

ggarrange(plt1, plt2, plt3,
          plt4, plt5, plt6,
          plt7, plt8, plt9,
          plt10, plt11, plt12,
          ncol = 2, 
          nrow = 3)




plt1 <- plot_temporal_trend_data_all_data_set("sc1", "ADM1 const, short: all data sets")
plt2 <- plot_temporal_trend_data_all_data_set("sc2", "ADM4 const, short: all data sets")
plt3 <-plot_temporal_trend_data_all_data_set("sc3", "ADM1 lin increasing, short: all data sets")
plt4 <-plot_temporal_trend_data_all_data_set("sc4", "ADM4 lin increasing, short: all data sets")
plt5 <-plot_temporal_trend_data_all_data_set("sc5", "ADM1 change point, short: all data sets")
plt6 <-plot_temporal_trend_data_all_data_set("sc6", "ADM4 change point, short: all data sets")
plt7 <-plot_temporal_trend_data_all_data_set("sc7", "ADM1 const, long: all data sets")
plt8 <-plot_temporal_trend_data_all_data_set("sc8", "ADM4 const, long: all data sets")
plt9 <-plot_temporal_trend_data_all_data_set("sc9", "ADM1 lin increasing, long: all data sets")
plt10 <-plot_temporal_trend_data_all_data_set("sc10", "ADM4 lin increasing, long: all data sets")
plt11 <-plot_temporal_trend_data_all_data_set("sc11", "ADM1 change point, long: all data sets")
plt12 <-plot_temporal_trend_data_all_data_set("sc12", "ADM4 change point, long: all data sets")

ggarrange(plt1, plt2, plt3,
          plt4, plt5, plt6,
          plt7, plt8, plt9,
          plt10, plt11, plt12,
          ncol = 2, 
          nrow = 3)



################################################################################
# Load in risk surface for ADM1, and see how the data looks

dir = "Data/Simulated_risk_surfaces/"
columnNames_risk_surface = c("t", "polygon_id", "geometry", "x", "y", 
                "unique_id", "time_id", 
                "first_level_area_id_mapping", 
                "second_level_area_id_mapping")

columnNames_lambda = c("area_id", "time_id", "E_it", "space.time")

#####
#SC1
scenario_name = "sc1"
load(paste(dir, scenario_name, "_risk_surfaces.RData", sep = ""))
risk_surface.list_sc1 = risk_surface.list[, columnNames_risk_surface]
risk_surface.list_sc1$values = risk_surface.list$values[, 1]
rm(risk_surface.list)


plt1 <- heatmap_points(risk_surface.list_sc1,
               polygon_grid2, 
               first_level_admin_map,
               t_axis[1], paste("t: ", toString(round(t_axis[1], 1)),
                                sep = ""))

plt2 <- heatmap_points(risk_surface.list_sc1,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[2], paste("t: ", toString(round(t_axis[2], 1)),
                                        sep = ""))

plt3 <- heatmap_points(risk_surface.list_sc1,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[3], paste("t: ", toString(round(t_axis[3], 1)),
                                        sep = ""))
# Save as sc1_t_1_continuous_risk 12 by 6
ggarrange(plt1, plt2, plt3,
          ncol = 3, nrow = 1, common.legend = T)


### Load in the area-specific rate

load(paste("Simulated_data/sc1/", scenario_name, "_data.RData", sep = ""))
lambda_ <- lambda.df[, columnNames_lambda]
lambda_$lambda_it <- lambda.df$lambda_it[, 1]

tmp_map_ = first_level_admin_map
tmp_map_$lambda_it <- lambda_[lambda_$time_id == 1, ]$lambda_it

# Save as sc1_t_1_discrete_rate 5 by 5
heatmap_areas(tmp_map_, tmp_map_$lambda_it, title = "Rate for year 1")



#####
#SC7
scenario_name = "sc7"
load(paste(dir, scenario_name, "_risk_surfaces.RData", sep = ""))
risk_surface.list_sc7 = risk_surface.list[, columnNames_risk_surface]
risk_surface.list_sc7$values = risk_surface.list$values[, 1]
rm(risk_surface.list)


plt1 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[1], paste("t: ", toString(round(t_axis[1], 1)),
                                        sep = ""))

plt2 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[2], paste("t: ", toString(round(t_axis[2], 1)),
                                        sep = ""))

plt3 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[3], paste("t: ", toString(round(t_axis[3], 1)),
                                        sep = ""))
# Save as sc7_t_1_continuous_risk 12 by 6
ggarrange(plt1, plt2, plt3,
          ncol = 3, nrow = 1, common.legend = T)



scenario_name = "sc2"
load(paste(dir, scenario_name, "_risk_surfaces.RData", sep = ""))
risk_surface.list_sc7 = risk_surface.list[, columnNames_risk_surface]
risk_surface.list_sc7$values = risk_surface.list$values[, 1]
rm(risk_surface.list)


plt1 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       second_level_admin_map,
                       t_axis[1], paste("t: ", toString(round(t_axis[1], 1)),
                                        sep = ""))

plt2 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       second_level_admin_map,
                       t_axis[2], paste("t: ", toString(round(t_axis[2], 1)),
                                        sep = ""))

plt3 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       second_level_admin_map,
                       t_axis[3], paste("t: ", toString(round(t_axis[3], 1)),
                                        sep = ""))
# Save as sc7_t_1_continuous_risk 12 by 6
ggarrange(plt1, plt2, plt3,
          ncol = 3, nrow = 1, common.legend = T)


plt1 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[length(t_axis) - 2], 
                       paste("t: ", toString(round(t_axis[length(t_axis) - 2], 1)),
                                        sep = ""))

plt2 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[length(t_axis) - 1], 
                       paste("t: ", toString(round(t_axis[length(t_axis) - 1], 1)),
                                        sep = ""))

plt3 <- heatmap_points(risk_surface.list_sc7,
                       polygon_grid2, 
                       first_level_admin_map,
                       t_axis[length(t_axis)], 
                       paste("t: ", toString(round(t_axis[length(t_axis) ], 1)),
                                        sep = ""))
# Save as sc7_t_13_continuous_risk 12 by 6
ggarrange(plt1, plt2, plt3,
          ncol = 3, nrow = 1, common.legend = T)


### Load in the area-specific rate

load(paste("Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
lambda_ <- lambda.df[, columnNames_lambda]
lambda_$lambda_it <- lambda.df$lambda_it[, 1]

tmp_map_ = first_level_admin_map
tmp_map_$lambda_it <- lambda_[lambda_$time_id == 1, ]$lambda_it

# Save as sc7_t_1_discrete_rate 5 by 5
plt1 <- heatmap_areas(tmp_map_, tmp_map_$lambda_it, title = "Rate for year 1")

tmp_map_ = first_level_admin_map
tmp_map_$lambda_it <- lambda_[lambda_$time_id == 13, ]$lambda_it

# Save as sc7_t_13_discrete_rate 5 by 5
plt2 <- heatmap_areas(tmp_map_, tmp_map_$lambda_it, title = "Rate for year 13")

ggarrange(plt1, plt2, 
          ncol = 2, nrow = 1,
          common.legend = F)








