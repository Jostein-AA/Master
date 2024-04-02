#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

n_sim = 100
E_it = 100

## Get the number of areas in the different maps
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)


################################################################################
# Functions for calculating MSE, IS,...

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








rate_IS_one_year_one_dataset <- function(){
  
}

rate_mean_IS_one_dataset <- function(){
  
}
###########################################




################################################################################
# load in data and start calculating

## Support function for trying to read a file
read_file <- function(file_name_){
  tryCatch(
  {
    # Try to read file_name_
    load(file_name_)
    
    return(marginals)
  },
  error = function(cond) {
    #message(paste("file does not exist:", file_name_))
    #message("Here's the original error message:")
    #message(conditionMessage(cond))
    
    # Choose a return value in case of error
    NULL
  },
  warning = function(cond) {
    #message(paste("file caused a warning:", file_name_))
    #message("Here's the original warning message:")
    #message(conditionMessage(cond))
    
    # Choose a return value in case of warning
    NULL
  },
  finally = {
    # Dont know how this stuff works
    a = 2
  }
  )
}


## Function that takes one model, one scenario, and calculates all the MSEs and IS' for ADM1
calc_mse_is_count_one_model_one_scenario_ADM1 <- function(model_name,
                                                          scenario_name,
                                                          n_sim){
  
  
  ### Initialize data frame to store the mse's and is's for each data set
  model_choice <- data.frame(mse_1_year_ahead = rep(NA, n_sim), mse_2_year_ahead = rep(NA, n_sim),
                             mse_3_year_ahead = rep(NA, n_sim), IS_1_year_ahead = rep(NA, n_sim),
                             IS_2_year_ahead = rep(NA, n_sim), IS_3_year_ahead = rep(NA, n_sim),
                             total_mse = rep(NA, n_sim), total_IS = rep(NA, n_sim))
  
  ## Find the number of files containing results
  
  predictions_dir = paste("./results/", model_name, "/",
                          scenario_name, sep = "")
  
  
  num_files <- length(list.files(predictions_dir))
  
  print(paste(model_name, scenario_name, num_files, sep = " - "))
  
  file_base_name = paste("/", model_name, "_", scenario_name, "_", sep = "")
  file_ending = ".RData"
  
  ## Iterate over data sets in scenario
  
  ### Initialize progressbar
  pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  
  #### First load in the entire simulated scenario (i.e. each of the simulated data sets of the scenario)
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData",
             sep = ""))
  lambda_df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
  
  #### Take the years predicted on
  years_pred_on <- 11:13
  
  #### Count how many data sets we open (i.e. if any are corrupted or something)
  successes = 0
  for(i in 1:n_sim){
    
    ## Should probably do a try-catch thing opening the files
    marginals = read_file(paste(predictions_dir, file_base_name, toString(i), file_ending,
                     sep = ""))
    
    ### If failure: jump to next
    if(is.null(marginals)){
      next
    }
    
    ### If not failure: count it as successfully opened!
    successes = successes + 1
    
    ### Load in the data: both counts (sampled_counts) and rate (lambda_it)
    lambda_df$sampled_counts <- lambda.df$sampled_counts[, i]
    lambda_df$lambda_it <- lambda.df$lambda[, i]
    
    #### For each year calculate the MSE and IS that year
    for(year in years_pred_on){
      ## Extract the predicted marginals for this year
      one_year_margs <- marginals[((year - 1) * n_ADM1 + 1):(year * n_ADM1)]
      
      ## Extract the sampled counts for this year
      one_year_sampled_counts = lambda_df$sampled_counts[((year - 1) * n_ADM1 + 1):(year * n_ADM1)]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      model_choice[i, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(
                                                              one_year_sampled_counts,
                                                              one_year_margs,
                                                              E_it)
      
      ## IS for count in year ahead (1, 2, or 3 years ahead)
      model_choice[i, year - years_pred_on[1] + 4] =  count_IS_one_year_one_dataset(
                                                                one_year_sampled_counts, 
                                                                one_year_margs, 
                                                                E_it)
      
    }
    
    ## Update progressbar
    setTxtProgressBar(pb, i)
  }
  
  ## Close the progressbar
  close(pb) 
  
  ### Check that number of successes equals number of files found
  if(successes != num_files){
    print("Number of successes not equal to number of files...")
    Sys.sleep(30)
  }
  
  filename = paste("./results/model_choice/", "model_choice_", model_name, "_", 
                   scenario_name, ".RData")
  
  save(model_choice,
       file = filename)
}


## Function that takes one model, one scenario, and calculates all the MSEs and IS' for ADM1
calc_mse_is_count_one_model_one_scenario_ADM4 <- function(model_name,
                                                          scenario_name,
                                                          n_sim){
  
  
  ### Initialize data frame to store the mse's and is's for each data set
  model_choice <- data.frame(mse_1_year_ahead = rep(NA, n_sim), mse_2_year_ahead = rep(NA, n_sim),
                             mse_3_year_ahead = rep(NA, n_sim), IS_1_year_ahead = rep(NA, n_sim),
                             IS_2_year_ahead = rep(NA, n_sim), IS_3_year_ahead = rep(NA, n_sim),
                             total_mse = rep(NA, n_sim), total_IS = rep(NA, n_sim))
  
  ## Find the number of files containing results
  
  predictions_dir = paste("./results/", model_name, "/",
                          scenario_name, sep = "")
  
  
  num_files <- length(list.files(predictions_dir))
  
  print(paste(model_name, scenario_name, num_files, sep = " - "))
  
  file_base_name = paste("/", model_name, "_", scenario_name, "_", sep = "")
  file_ending = ".RData"
  
  ## Iterate over data sets in scenario
  
  ### Initialize progressbar
  pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  
  #### First load in the entire simulated scenario (i.e. each of the simulated data sets of the scenario)
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData",
             sep = ""))
  lambda_df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
  
  #### Take the years predicted on
  years_pred_on <- 11:13
  
  #### Count how many data sets we open (i.e. if any are corrupted or something)
  successes = 0
  for(i in 1:n_sim){
    
    ## Should probably do a try-catch thing opening the files
    marginals = read_file(paste(predictions_dir, file_base_name, toString(i), file_ending,
                                sep = ""))
    
    ### If failure: jump to next
    if(is.null(marginals)){
      next
    }
    
    ### If not failure: count it as successfully opened!
    successes = successes + 1
    
    ### Load in the data: both counts (sampled_counts) and rate (lambda_it)
    lambda_df$sampled_counts <- lambda.df$sampled_counts[, i]
    lambda_df$lambda_it <- lambda.df$lambda[, i]
    
    #### For each year calculate the MSE and IS that year
    for(year in years_pred_on){
      ## Extract the predicted marginals for this year
      one_year_margs <- marginals[((year - years_pred_on[1]) * n_ADM4 + 1):((year - years_pred_on[1] + 1) * n_ADM4)]
      
      ## Extract the sampled counts for this year
      one_year_sampled_counts = lambda_df$sampled_counts[((year - 1) * n_ADM4 + 1):(year * n_ADM4)]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      model_choice[i, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(
        one_year_sampled_counts,
        one_year_margs,
        E_it)
      
      ## IS for count in year ahead (1, 2, or 3 years ahead)
      model_choice[i, year - years_pred_on[1] + 4] =  count_IS_one_year_one_dataset(
        one_year_sampled_counts, 
        one_year_margs, 
        E_it)
      
    }
    
    ## Update progressbar
    setTxtProgressBar(pb, i)
  }
  
  ## Close the progressbar
  close(pb) 
  
  ### Check that number of successes equals number of files found
  if(successes != num_files){
    print("Number of successes not equal to number of files...")
    Sys.sleep(30)
  }
  
  filename = paste("./results/model_choice/", "model_choice_", model_name, "_", 
                   scenario_name, ".RData")
  
  save(model_choice,
       file = filename)
}


################################################################################
model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
                "Improper1_typeIII", "Improper1_typeIV",
                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
                "Improper2_typeIII", "Improper2_typeIV",
                "proper1_noInt", "proper1_onlyInt", "proper1_full",
                "proper2_noInt", "proper2_onlyInt", "proper2_full")

scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")
scenario_names_ADM4 = c("sc2", "sc4", "sc6", "sc8", "sc10", "sc12")


### Iterate over each model for each scenario on ADM1 to calculate the model choice data frames
for(model_name in model_names){
  for(scenario_name in scenario_names_ADM1){
    ## Calculate the model choice
    calc_mse_is_count_one_model_one_scenario_ADM1(model_name = model_name,
                                       scenario_name = scenario_name,
                                       n_sim)
  }
}

### Iterate over each model for each scenario on ADM4 to calculate the model choice data frames
for(model_name in model_names){
  for(scenario_name in scenario_names_ADM4){
    calc_mse_is_count_one_model_one_scenario_ADM4(model_name = model_name,
                                                  scenario_name = scenario_name,
                                                  n_sim)
  }
}



#######
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
