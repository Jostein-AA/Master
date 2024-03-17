#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

n_sim = 100

## Get the number of areas in the different maps
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)


################################################################################
# Functions for calculating MSE, IS,...

file_path <- "results/Improper1_noInt/sc1/"
load(paste(file_path, "Improper1_noInt_sc1_2.RData", sep = ""))
load("./Simulated_data/sc1/sc1_data.RData")



count_mse_one_year_one_dataset <- function(sampled_counts_one_year, 
                                           lambda_marginals_one_year,
                                           E_it){
  
  ## For each area find the expected predicted count (i.e. point prediction)
  pred_count <- 100 * as.numeric(sapply(lambda_marginals_one_year, 
                          FUN = function(x){return(mean(x[, 1]))}))
  
  ## Calculate the MSE
  mse_one_year <- mean((sampled_counts_one_year - pred_count)**2)
  
  # Return the MSE of that year
  return(mse_one_year)
}

## Testing
#one_year_margs <- marginals[(11 * n_ADM1 + 1):(12 * n_ADM1)]
#count_mse_one_year_one_dataset(lambda.df$sampled_counts[1],one_year_margs,100)


count_mean_mse_one_dataset <- function(sampled_counts, lambda_marginals, n_ADM,
                                       E_it = 100){
  
  ## Take the years predicted on
  years_pred_on <- 11:13
  
  ## Initialize memory for mse
  mse = rep(-1, length(years_pred_on) + 1)
  
  ## For each year calculate the mse that year
  for(year in years_pred_on){
    one_year_margs <- marginals[((year - 1) * n_ADM + 1):(year * n_ADM)]
    
    mse[year - years_pred_on[1] + 1] = count_mse_one_year_one_dataset(
      sampled_counts[((year - 1) * n_ADM + 1):(year * n_ADM)], one_year_margs, E_it)
  }
  
  mse[length(mse)] = mean(mse[1:(length(mse) - 1)])
  
  return(mse)
}

rate_mse_one_year_one_dataset <- function(){
  
}

rate_mean_mse_one_dataset <- function(){
  
}

## Testing
#count_mean_mse_one_dataset(lambda.df$sampled_counts[, 1],marginals,n_ADM1)

find_ul_quants_counts_single_pred <- function(lambda_marginal,
                                              E_it){
  # Function to calculate upper (u) and lower (l) quantiles for a single
  # count prediction
  
  ## Sample lambda and scale w. E_it to get an instance of Poisson
  poisson_param_sample <- E_it * inla.rmarginal(5000, lambda_marginal) #[[1]]
  
  ## Sample from Poisson
  count_sample <- sapply(poisson_param_sample, 
                         FUN = function(x){return(rpois(1, x))})
  
  ## Calculate upper and lower quantile (and median)
  u = as.numeric(quantile(count_sample, 0.975)); l = as.numeric(quantile(count_sample, 0.025))
  median = as.numeric(quantile(count_sample, 0.5))
  
  return(list(l = l, u = u, median = median))
  
}

## Test for find_ul_quants_counts_single_pred
#find_ul_quants_counts_single_pred(marginals[12 * n_ADM1 + 1], 100)

#Interval-scores
find_IS_one_obs <- function(l, u, true_value){
  IS_score = (u - l) + 
              2/0.05 * (l - true_value) * (true_value < l) + 
              2/0.05 * (true_value - u) * (true_value > u) 
  return(IS_score)
}


count_IS_one_year_one_dataset <- function(sampled_counts_one_year,
                                          lambda_marginals_one_year,
                                          E_it){
  
  ## Find the l and u for each singular instance in a year
  ul_each_one_year <- lapply(lambda_marginals_one_year, 
                        FUN = function(x){
                          return(find_ul_quants_counts_single_pred(x, 100))
                          })
  
  ## Find the IS for each singular instance in a year
  ### Initialize space for IS 
  IS_each_instance = rep(0, length(lambda_marginals_one_year))
  
  ### Calculate IS each area
  for(i in 1:length(lambda_marginals_one_year)){
    IS_each_instance[i] = find_IS_one_obs(ul_each_one_year[[i]]$l, ul_each_one_year[[i]]$u, 
                                          sampled_counts_one_year[i])
  }
  
  ## Find the average IS this year
  IS_this_year <- mean(IS_each_instance)
  
  return(IS_this_year)
}

## Test of count_IS_one_year_one_dataset
#test2 <- count_IS_one_year_one_dataset(lambda.df$sampled_counts[(10 * n_ADM1 + 1):(11 * n_ADM1), 1], marginals[(10 * n_ADM1 + 1):(11 * n_ADM1)], 100)



count_mean_IS_one_dataset <- function(sampled_counts,
                                      lambda_marginals,
                                      n_ADM,
                                      E_it){
  
  ## Take the years predicted on
  years_pred_on <- 11:13
  
  ## Initialize memory for mse
  IS = rep(-1, length(years_pred_on) + 1)
  
  ## For each year calculate the mse that year
  for(year in years_pred_on){
    one_year_margs <- lambda_marginals[((year - 1) * n_ADM + 1):(year * n_ADM)]
    
    IS[year - years_pred_on[1] + 1] = count_IS_one_year_one_dataset(
                                      sampled_counts[((year - 1) * n_ADM + 1):(year * n_ADM)], 
                                      one_year_margs, 
                                      E_it)
  }
  
  IS[length(IS)] = mean(IS[1:(length(IS) - 1)])
  
  return(IS)
  
}
## Testing
#test = count_mean_IS_one_dataset(lambda.df$sampled_counts[, 1], marginals, n_ADM1, 100)




rate_IS_one_year_one_dataset <- function(){
  
}

rate_mean_IS_one_dataset <- function(){
  
}
###########################################
# Function for loading in all the data for one model and one scenario, that calculates everything into one df
calc_model_choice <- function(model_name, scenario_name, lambda, n_sim){
  
  ## Create a file path for the predictions 
  predictions_file_path = paste("./results/", model_name, "/", 
                    scenario_name, "/", 
                    model_name, "_", scenario_name, "_", 
                    sep = "")
  
  file_ending = ".RData"
  
  ## Initialize a data frame to store results of each data set into
  model_choice <- data.frame(mse_1_year_ahead = rep(-1, n_sim), mse_2_year_ahead = rep(-1, n_sim),
                             mse_3_year_ahead = rep(-1, n_sim), IS_1_year_ahead = rep(-1, n_sim),
                             IS_2_year_ahead = rep(-1, n_sim), IS_3_year_ahead = rep(-1, n_sim),
                             total_mse = rep(-1, n_sim), total_IS = rep(-1, n_sim))
  
  ## Initialize a data frame to store the end results into
  model_choice_total <- data.frame(mse_1_year_ahead = rep(-1, 1), mse_2_year_ahead = rep(-1, 1),
                                   mse_3_year_ahead = rep(-1, 1), IS_1_year_ahead = rep(-1, 1),
                                   IS_2_year_ahead = rep(-1, 1), IS_3_year_ahead = rep(-1, 1),
                                   total_mse = rep(-1, 1), total_IS = rep(-1, 1))
  
  for(i in 1:n_sim){
    
  
    ## Load predicted linear predictors (marginals) and log-score (cpo)
    pred_filename = paste(predictions_file_path, i, 
                          file_ending, sep = "")
    load(pred_filename)
    
    ## Extract the right data
    lambda_ <- lambda[, c("area_id", "time_id", "E_it", "space.time")]
    lambda_$lambda_it <- lambda$lambda_it[, i]
    lambda_$sampled_counts <- lambda$sampled_counts[, i]
    
    ## Calculate the MSE, IS of this data set
    
    
  }
  
  ## Calculate the total MSE, IS
  
  
}





calc_model_choice("Improper1_noInt", "sc1", lambda.df, n_sim)



################################################################################
# load in data and start calculating



