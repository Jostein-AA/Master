#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

## Necessary to define temporal domain's max
tT = 13

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

