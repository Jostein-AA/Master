#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("geofacet")
library(latex2exp)

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

## Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)
dataset_id = 3 
dataset_id_2 = 4

scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")


###
# Functions used for plotting

select_regions_lin_pred_vs_true <- function(geofacet_grid,
                                            pred_to_plot,
                                            title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_point(aes(x = time_id, y = lambda_it, col = "True rate")) + 
    geom_point((aes(x = time_id, y = count_div_pop, col = "sampled count/Eit"))) +
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) + 
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=11),
          plot.title = element_text(size=11),
          strip.text.x = element_text(size = 9),
          legend.position = c(0, 1),
          legend.justification = c("left", "top"),
          legend.box.just = "left")
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue")) # 
  return(plt)
}

select_regions_lin_pred_vs_true_improper <- function(ADM1_grid, lambda_, model, title){
  
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = model$summary.fitted.values$'0.5quant', 
                             quantile_0.025 = model$summary.fitted.values$'0.025quant', 
                             quantile_0.975 = model$summary.fitted.values$'0.975quant', 
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it,
                             count_div_pop = lambda_$sampled_counts/lambda_$E_it)
  
  select_regions_lin_pred_vs_true(ADM1_grid, pred_to_plot, title)
  
}

select_regions_lin_pred_vs_true_proper <- function(ADM1_grid, lambda_, model, title){
### NB: Must sort the proper ones
pred_to_plot <- data.frame(area_id = lambda_$area_id,
                           time_id = lambda_$time_id,
                           median = sort_proper_fitted(model$summary.fitted.values$'0.5quant', n_ADM1, tT), # model$summary.fitted.values$'0.5quant',
                           quantile_0.025 = sort_proper_fitted(model$summary.fitted.values$'0.025quant', n_ADM1, tT),# model$summary.fitted.values$'0.025quant',
                           quantile_0.975 = sort_proper_fitted(model$summary.fitted.values$'0.975quant', n_ADM1, tT), #model$summary.fitted.values$'0.975quant',
                           sampled_counts = lambda_$sampled_counts,
                           lambda_it = lambda_$lambda_it,
                           count_div_pop = lambda_$sampled_counts/lambda_$E_it)

select_regions_lin_pred_vs_true(ADM1_grid, pred_to_plot, title)
}


select_regions_lin_pred_vs_true_2 <- function(geofacet_grid,
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
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) + 
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=11),
          plot.title = element_text(size=11),
          strip.text.x = element_text(size = 9))
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue")) # 
  return(plt)
}




################################################################################
## ADM1
## Make a geofacet grid to plot onto
#ADM1_grid <- data.frame(
#  code = c("15", "8", "6", "3", "5", "9", "4", "13", "10", "14", "16", "7", "11", "12", "2", "1"),
#  name = c("Schleswig-Holstein", 
#           "Mecklenburg-Vorpommern", 
#           "Hamburg", 
#           "Berlin", 
#           "Bremen", 
#           "Niedersachsen", 
#           "Brandenburg", 
#           "Sachsen-Anhalt", 
#           "Nordrhein-Westfalen", 
#           "Sachsen", 
#           "Thuringen", 
#           "Hessen", 
#           "Rheinland-Pfalz", 
#           "Saarland", 
#           "Bayern", 
#           "Baden-Wurttemberg"),
#  row = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5),
#  col = c(2, 4, 3, 4, 2, 3, 5, 4, 2, 5, 4, 3, 2, 1, 4, 3),
#  stringsAsFactors = FALSE)
#geofacet::grid_preview(ADM1_grid)

ADM1_grid <- data.frame(
  row = c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5),
  col = c(3, 4, 4, 6, 3, 2, 6, 3, 4, 3, 2, 6, 4, 1, 4, 3),
  code = c("15", "8", "4", "6", "9", "10", "5", "13", "14", "7", "11", "3", "16", "12", "2", "1"),
  name = c("Schleswig-Holstein", "Mecklenburg-Vorpommern", "Brandenburg", "Hamburg", "Niedersachsen", "Nordrhein-Westfalen", "Bremen", "Sachsen-Anhalt", "Sachsen", "Hessen", "Rheinland-Pfalz", "Berlin", "Thuringen", "Saarland", "Bayern", "Baden-Wurttemberg"),
  stringsAsFactors = FALSE
)
geofacet::grid_preview(ADM1_grid)

################################################
### Plot the rates

scenario_name = "sc7"
title = TeX(r'(ADM1$_{const, long}$)')

### Load in simulated data for that scenario
load(paste("./Simulated_data/", scenario_name, "/", 
           scenario_name, "_data.RData", sep = ""))

lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]

#### Remove unessecary data
rm(lambda.df)

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))


### Plot a proper model
model = proper2_onlyInt_ADM1
select_regions_lin_pred_vs_true_proper(ADM1_grid, lambda_, model, title)

### Plot an improper model
model = Improper1_typeI_ADM1
select_regions_lin_pred_vs_true_improper(ADM1_grid, lambda_, model, title)


################################################
### Do the same but use the marginals to make plots on the observational level

scenario_name = "sc3"
title = TeX(r'(ADM1$_{conts, short}$)')

### Load in simulated data for that scenario
load(paste("./Simulated_data/", scenario_name, "/", 
           scenario_name, "_data.RData", sep = ""))

lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]

#### Remove unessecary data
rm(lambda.df)

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))
model = typeIV_1_first_level

ul_each <- lapply(model$marginals.fitted.values, 
                   FUN = function(x){
                     return(find_ul_quants_counts_single_pred(x, 100))
                   })

### NB: Must sort the proper models

pred_to_plot <- data.frame(area_id = lambda_$area_id,
                           time_id = lambda_$time_id,
                           median = rep(NA, nrow(lambda_)),
                           quantile_0.025 = rep(NA, nrow(lambda_)),
                           quantile_0.975 = rep(NA, nrow(lambda_)),
                           sampled_counts = lambda_$sampled_counts,
                           lambda_it = lambda_$lambda_it * lambda_$E_it,
                           count_div_pop = lambda_$sampled_counts/lambda_$E_it)

for(i in 1:nrow(pred_to_plot)){
  pred_to_plot[i, ]$median = ul_each[[i]]$median
  pred_to_plot[i, ]$quantile_0.025 = ul_each[[i]]$l
  pred_to_plot[i, ]$quantile_0.975 = ul_each[[i]]$u
}

select_regions_lin_pred_vs_true_2(ADM1_grid, pred_to_plot, title)

################################################################################

# Plot some of the ADM4

## Make a geofacet grid to plot onto
ADM4_grid <- data.frame(
  code = c("1", "50", "100", "150", "200", "250", "300", "350", "400"),
  name = c("Alb-Donau-Kreis",
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


##################################
# Rate as it is

scenario_name = "sc2"
title = TeX(r'(ADM4$_{const, short}$)')

### Load in simulated data for that scenario
load(paste("./Simulated_data/", scenario_name, "/", 
           scenario_name, "_data.RData", sep = ""))

lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]

#### Remove unessecary data
rm(lambda.df)

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))
rm(proper1_full_second_level, proper2_full_second_level)

model = typeIV_1_second_level

# Must sort the proper ones


#proper_summary_fitted_values = sort_proper_fitted(model$summary.fitted.values, n_ADM1, tT)

pred_to_plot <- data.frame(area_id = lambda_$area_id,
                           time_id = lambda_$time_id,
                           median = model$summary.fitted.values$'0.5quant',#sort_proper_fitted(model$summary.fitted.values$'0.5quant', n_ADM4, tT), # 
                           quantile_0.025 = model$summary.fitted.values$'0.025quant',#sort_proper_fitted(model$summary.fitted.values$'0.025quant', n_ADM4, tT), # 
                           quantile_0.975 = model$summary.fitted.values$'0.975quant',#sort_proper_fitted(model$summary.fitted.values$'0.975quant', n_ADM4, tT), # 
                           sampled_counts = lambda_$sampled_counts,
                           lambda_it = lambda_$lambda_it,
                           count_div_pop = lambda_$sampled_counts/lambda_$E_it)

select_regions_lin_pred_vs_true(ADM4_grid, pred_to_plot, title)



##################################
# Counts

scenario_name = "sc2"
title = TeX(r'(ADM4$_{const, short}$)')

### Load in simulated data for that scenario
load(paste("./Simulated_data/", scenario_name, "/", 
           scenario_name, "_data.RData", sep = ""))

lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]

#### Remove unessecary data
rm(lambda.df)

### Load in the models for that scenario
load(paste("diagnostics_", scenario_name, ".RData", sep = ""))
rm(proper1_full_second_level, proper2_full_second_level)

model = proper1_onlyInt_second_level

ul_each <- lapply(model$marginals.fitted.values, 
                  FUN = function(x){
                    return(find_ul_quants_counts_single_pred(x, 100))
                  })

### Must sort the proper ones
#proper_summary_fitted_values = sort_proper_fitted(model$summary.fitted.values, n_ADM1, tT)

pred_to_plot <- data.frame(area_id = lambda_$area_id,
                           time_id = lambda_$time_id,
                           median = rep(NA, nrow(lambda_)),
                           quantile_0.025 = rep(NA, nrow(lambda_)),
                           quantile_0.975 = rep(NA, nrow(lambda_)),
                           sampled_counts = lambda_$sampled_counts,
                           lambda_it = lambda_$lambda_it * lambda_$E_it,
                           count_div_pop = lambda_$sampled_counts/lambda_$E_it)

for(i in 1:nrow(pred_to_plot)){
  pred_to_plot[i, ]$median = ul_each[[i]]$median
  pred_to_plot[i, ]$quantile_0.025 = ul_each[[i]]$l
  pred_to_plot[i, ]$quantile_0.975 = ul_each[[i]]$u
}

select_regions_lin_pred_vs_true_2(ADM4_grid, pred_to_plot, title)















