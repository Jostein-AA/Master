#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("geofacet")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

## Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)


###



### state a scenario name
# scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")
# scenario_names_ADM4 = c("sc2", "sc4", "sc6", "sc8", "sc10", "sc12")

scenario_name = "sc5"


### Load in simulated data for that scenario
load(paste("./Simulated_data/", scenario_name, "/", 
           scenario_name, "_data.RData", sep = ""))

lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
lambda_$sampled_counts <- lambda.df$sampled_counts[, 2]
lambda_$lambda_it <- lambda.df$lambda_it[, 2]

#### Remove unessecary data
rm(lambda.df)

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))

model = typeIV_1_first_level

# Plot the 95% CIs over all the 13 years (also for those predicted on) for areas

### Extract the first 10 years directly
pred_to_plot <- data.frame(area_id = lambda_$area_id,
                          time_id = lambda_$time_id,
                          median = rep(NA, nrow(lambda_)),
                          quantile_0.025 = rep(NA, nrow(lambda_)),
                          quantile_0.975 = rep(NA, nrow(lambda_)),
                          sampled_counts = lambda_$sampled_counts)



#sapply(base_1_first_level$marginals.fitted.values, 
#       FUN = function(x){return(find_ul_quants_counts_single_pred(x, E_it))})



#<- find_ul_quants_counts_single_pred(lambda_marginal, E_it)

#base_1_first_level$marginals.fitted.values



pred_to_plot$median = model$summary.fitted.values$'0.5quant' 
pred_to_plot$quantile_0.025 = model$summary.fitted.values$'0.025quant' 
pred_to_plot$quantile_0.975 = model$summary.fitted.values$'0.975quant' 


### Plot
region_time_series <- function(pred_to_plot, #true_risk
                               region){
  
  tmp_ = pred_to_plot[pred_to_plot$area_id == region, ]
  
  plt <- ggplot(data = tmp_, aes(x = time_id)) + ggtitle(paste("region: ", region)) +
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = "95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) +
    xlab("Year") + ylab("") + 
    labs(col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=9),
          plot.title = element_text(size=10))
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black")) # "#00BFC4", "blue"
  return(plt)
}

region_time_series(pred_to_plot = pred_to_plot,
                   region = 9)





#grid_design()

ADM1_grid <- data.frame(
  code = c("15", "8", "6", "3", "5", "9", "4", "13", "10", "14", "16", "7", "11", "12", "2", "1"),
  name = c(" Schleswig-Holstein ", " Mecklenburg-Vorpommern ", " Hamburg", " Berlin ", " Bremen ", " Niedersachsen", " Brandenburg", " Sachsen-Anhalt", " Nordrhein-Westfalen", " Sachsen", " Th??ringen", " Hessen ", " Rheinland-Pfalz ", " Saarland", " Bayern", " Baden-W??rttemberg "),
  row = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5),
  col = c(2, 4, 3, 4, 2, 3, 5, 4, 2, 5, 4, 3, 2, 1, 4, 3),
  stringsAsFactors = FALSE
)


facet_geo(grid = "de_states_grid1")

geofacet::grid_preview(ADM1_grid)


select_regions_lin_pred_vs_true <- function(true_risk,
                                            model,
                                            regions,
                                            n, 
                                            tT){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  plt1 <- region_time_series(true_risk, model, regions[1], n, tT)
  plt2 <- region_time_series(true_risk, model, regions[2], n, tT)
  plt3 <- region_time_series(true_risk, model, regions[3], n, tT)
  plt4 <- region_time_series(true_risk, model, regions[4], n, tT)
  plt5 <- region_time_series(true_risk, model, regions[5], n, tT)
  plt6 <- region_time_series(true_risk, model, regions[6], n, tT)
  plt7 <- region_time_series(true_risk, model, regions[7], n, tT)
  plt8 <- region_time_series(true_risk, model, regions[8], n, tT)
  
  
  ggarrange(plt1, plt2, plt3, 
            plt4, plt5, plt6,
            plt7, plt8,
            ncol = 2, nrow = 4, 
            common.legend = TRUE, legend = "top")
  
  
}





################################################################################

load("improper1_noInt_fitted.RData")
load("improper1_typeI_fitted.RData")

#Load in simulation data (done stupidly, change!)
load("./Simulated_data/sc1_data_1.RData")
lambda_sc1.df <- lambda.df
lambda_sc1.df$space.time = 1:nrow(lambda_sc1.df)

load("./Simulated_data/sc2_data_1.RData")
lambda_sc2.df <- lambda.df
lambda_sc2.df$space.time = 1:nrow(lambda_sc2.df)

load("./Simulated_data/sc3_data_1.RData")
lambda_sc3.df <- lambda.df
lambda_sc3.df$space.time = 1:nrow(lambda_sc3.df)

load("./Simulated_data/sc4_data_1.RData")
lambda_sc4.df <- lambda.df
lambda_sc4.df$space.time = 1:nrow(lambda_sc4.df)

load("./Simulated_data/sc5_data_1.RData")
lambda_sc5.df <- lambda.df
lambda_sc5.df$space.time = 1:nrow(lambda_sc5.df)

load("./Simulated_data/sc6_data_1.RData")
lambda_sc6.df <- lambda.df
lambda_sc6.df$space.time = 1:nrow(lambda_sc6.df)

load("./Simulated_data/sc7_data_1.RData")
lambda_sc7.df <- lambda.df
lambda_sc7.df$space.time = 1:nrow(lambda_sc7.df)

load("./Simulated_data/sc8_data_1.RData")
lambda_sc8.df <- lambda.df
lambda_sc8.df$space.time = 1:nrow(lambda_sc8.df)

load("./Simulated_data/sc9_data_1.RData")
lambda_sc9.df <- lambda.df
lambda_sc9.df$space.time = 1:nrow(lambda_sc9.df)

load("./Simulated_data/sc10_data_1.RData")
lambda_sc10.df <- lambda.df
lambda_sc10.df$space.time = 1:nrow(lambda_sc10.df)

load("./Simulated_data/sc11_data_1.RData")
lambda_sc11.df <- lambda.df
lambda_sc11.df$space.time = 1:nrow(lambda_sc11.df)

load("./Simulated_data/sc12_data_1.RData")
lambda_sc12.df <- lambda.df
lambda_sc12.df$space.time = 1:nrow(lambda_sc12.df)



################################################################################

# Just plot the fits directly and see
plot(improper_noInt_sc1)
plot(improper_noInt_sc2)
plot(improper_noInt_sc3)
plot(improper_noInt_sc4)
plot(improper_noInt_sc5)
plot(improper_noInt_sc6)
plot(improper_noInt_sc7)
plot(improper_noInt_sc8)
plot(improper_noInt_sc9)
plot(improper_noInt_sc10)
plot(improper_noInt_sc11)
plot(improper_noInt_sc12)

plot(improper_typeI_sc1)
plot(improper_typeI_sc2)
plot(improper_typeI_sc3)
plot(improper_typeI_sc4)
plot(improper_typeI_sc5)
plot(improper_typeI_sc6)
plot(improper_typeI_sc7)
plot(improper_typeI_sc8)
plot(improper_typeI_sc9)
plot(improper_typeI_sc10)
plot(improper_typeI_sc11)
plot(improper_typeI_sc12)

## Plot the fitted linear pred. for 6 areas over time w. 95 %CIs against true values
regions = c(1, 3, 5,
            7, 9, 10,
            12, 14, 16)

#region_time_series(true_risk = lambda.df, model = improper_noInt,
#                   region = regions[6], n = nrow(germany_map_2), tT = tT)

## Plot the fitted linear predictor vs the true values for the regions = regions
select_regions_lin_pred_vs_true(true_risk = lambda_sc1.df,
                                model = improper_typeI_sc1,
                                regions = regions,
                                n = nrow(germany_map_2),
                                tT = tT)


plot(improper_noInt_sc4)

## Plot the fitted linear pred. for 6 areas over time w. 95 %CIs against true values
regions = c(1, 50, 100,
            150, 200, 250,
            300, 350, 400)

#region_time_series(true_risk = lambda.df, model = improper_noInt,
#                   region = regions[6], n = nrow(germany_map_2), tT = tT)

## Plot the fitted linear predictor vs the true values for the regions = regions
select_regions_lin_pred_vs_true(true_risk = lambda_sc4.df,
                                model = improper_noInt_sc4,
                                regions = regions,
                                n = nrow(germany_map),
                                tT = tT)

