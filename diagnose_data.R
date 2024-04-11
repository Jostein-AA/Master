
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

plot_temporal_trend_data_one_data_set <- function(sim_data, title){

  # Aggregate over the years
  aggr <- aggregate(sim_data, by = list(time_id = sim_data$time_id), FUN = mean) 
  
  # Plot over the years
  return(ggplot(data = aggr[, 2:ncol(aggr)]) + ggtitle(title) + 
    geom_line(aes(x = time_id, y = mu, col = "mu"), color = "blue") + 
    geom_point(aes(x = time_id, y = sampled_counts, col = "sampled counts"), color = "black") +
    ylim(7.5, 12.5))
}


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


plot_temporal_trend_data_all_data_set <- function(scenario_name, title){
  load(paste("./Simulated_data/", scenario_name,"/", scenario_name, "_data.RData",
             sep = ""))
  
  tmp_ = lambda.df[, c("space.time", "area_id", "time_id", "E_it", "sampled_counts", "mu")]
  
  # Aggregate the sampled counts and mu over the data sets
  
  tmp_$sampled_counts <- rowMeans(tmp_$sampled_counts)
  tmp_$mu <- rowMeans(tmp_$mu)
  
  # Aggregate over the years
  aggr <- aggregate(tmp_, by = list(time_id = tmp_$time_id), FUN = mean)
  
  # Plot over the years
  return(ggplot(data = aggr[, 2:ncol(aggr)]) + ggtitle(title) + 
           geom_line(aes(x = time_id, y = mu, col = "mu"), color = "blue") + 
           geom_point(aes(x = time_id, y = sampled_counts, col = "sampled counts"), color = "black") +
           ylim(7.5, 12.5))
}

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









