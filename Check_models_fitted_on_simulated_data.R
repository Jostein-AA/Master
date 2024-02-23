#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Preliminaries.R")
source("Utilities.R")

## Necessary to define temporal domain's max
tT = 13

load("test.RData")
load("test2.RData")

load("./Data/Simulated_data/sc1_data_1.RData")
lambda_sc1.df <- lambda.df

load("./Data/Simulated_data/sc4_data_1.RData")
lambda_sc4.df <- lambda.df

load("./Data/Simulated_data/sc8_data_1.RData")
lambda_sc8.df <- lambda.df

################################################################################

plot(improper_typeI_sc1)

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

