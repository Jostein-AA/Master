#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("geofacet")
library(ggh4x)
library(ggridges)
library(latex2exp)
library(geoR)
library(paletteer)
library('ggsci')
library(ggstats)

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
dataset_id_new = 2



scenario_names_ADM4 <- c("sc2", "sc4", "sc6",
                         "sc8", "sc10", "sc12")



################################################################################

model_names <- c("Improper1_typeIV",
                 "Improper1_typeV")


# 
ridgeplot_mse_is_rates_all_years(model_names = model_names,
                                 "sc2", xlim_mse = c(0.35, 3.5), 
                                 xlim_is = c(2, 11))


ridgeplot_mse_is_rates_all_years(model_names = model_names,
                                 "sc4", xlim_mse = c(0.35, 3.5), 
                                 xlim_is = c(2, 11))


ridgeplot_mse_is_rates_all_years(model_names = model_names,
                                 "sc6", xlim_mse = c(0.35, 3.5), 
                                 xlim_is = c(2, 11))

ridgeplot_mse_is_rates_all_years(model_names = model_names,
                                 "sc8", xlim_mse = c(0.35, 3.5), 
                                 xlim_is = c(2, 11))


ridgeplot_mse_is_rates_all_years(model_names = model_names,
                                 "sc10", xlim_mse = c(0.35, 3.5), 
                                 xlim_is = c(2, 11))


ridgeplot_mse_is_rates_all_years(model_names = model_names,
                                 "sc12", xlim_mse = c(0.35, 3.5), 
                                 xlim_is = c(2, 11))

################################################################################



















































