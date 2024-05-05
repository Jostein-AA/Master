#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library(latex2exp)
library(tables)

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

n_sim = 100
E_it = 100

## Get the number of areas in the different maps
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)
n_ADMnew <- nrow(new_map)




scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")
scenario_names_ADM4 = c("sc2", "sc4", "sc6", "sc8", "sc10", "sc12")
scenario_names_new = c("sc13", "sc14", "sc15", "sc16", "sc17", "sc18")


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
calc_width_CI_one_scenario_ADM1 <- function(model_name,
                                                          scenario_name,
                                                          n_sim){
  
  
  ### Initialize data frame to store the mse's and is's for each data set
  width_CIs <- data.frame(width_CI_1_year_ahead = rep(NA, n_sim), width_CI_2_year_ahead = rep(NA, n_sim),
                          width_CI_3_year_ahead = rep(NA, n_sim), width_CI_avg = rep(NA, n_sim))
  
  
  
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
    
    #### For each year calculate the MSE and IS that year
    for(year in years_pred_on){
      ## Extract the predicted marginals for this year
      one_year_margs <- marginals[((year - 1) * n_ADM1 + 1):(year * n_ADM1)]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      width_CIs[i, year - years_pred_on[1] + 1] =  avg_width(one_year_margs, E_it)
      
    }
    
    #Get the average CI width for count
    width_CIs[i, 4] = mean(as.numeric(width_CIs[i, 1:3]))
    
    
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
  
  filename = paste("./results/width_CIs/", "width_CIs", model_name, "_", 
                   scenario_name, ".RData", sep = "")
  
  save(width_CIs,
       file = filename)
}


## Function that takes one model, one scenario, and calculates all the MSEs and IS' for ADM4
calc_width_CI_one_scenario_ADM4 <- function(model_name,
                                                          scenario_name,
                                                          n_sim){
  
  ### Initialize data frame to store the mse's and is's for each data set
  width_CIs <- data.frame(width_CI_1_year_ahead = rep(NA, n_sim), width_CI_2_year_ahead = rep(NA, n_sim),
                          width_CI_3_year_ahead = rep(NA, n_sim), width_CI_avg = rep(NA, n_sim))
  
  
  
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
    
    #### For each year calculate the MSE and IS that year
    for(year in years_pred_on){
      ## Extract the predicted marginals for this year
      one_year_margs <- marginals[((year - years_pred_on[1]) * n_ADM4 + 1):((year - years_pred_on[1] + 1) * n_ADM4)]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      width_CIs[i, year - years_pred_on[1] + 1] =  avg_width(one_year_margs, E_it)
      
    }
    
    #Get the average CI width for count
    width_CIs[i, 4] = mean(as.numeric(width_CIs[i, 1:3]))
    
    
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
  
  filename = paste("./results/width_CIs/", "width_CIs", model_name, "_", 
                   scenario_name, ".RData", sep = "")
  
  save(width_CIs,
       file = filename)
  
}




## Function that takes one model, one scenario, and calculates all the MSEs and IS' for ADM1
calc_width_CI_one_scenario_new <- function(model_name,
                                                         scenario_name,
                                                         n_sim){
  
  ### Initialize data frame to store the mse's and is's for each data set
  width_CIs <- data.frame(width_CI_1_year_ahead = rep(NA, n_sim), width_CI_2_year_ahead = rep(NA, n_sim),
                          width_CI_3_year_ahead = rep(NA, n_sim), width_CI_avg = rep(NA, n_sim))
  
  
  
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
    
    #### For each year calculate the MSE and IS that year
    for(year in years_pred_on){
      ## Extract the predicted marginals for this year
      one_year_margs <- marginals[((year - years_pred_on[1]) * nrow(new_map) + 1):((year - years_pred_on[1] + 1) * nrow(new_map))]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      width_CIs[i, year - years_pred_on[1] + 1] =  avg_width(one_year_margs, E_it)
      
    }
    
    #Get the average CI width for count
    width_CIs[i, 4] = mean(as.numeric(width_CIs[i, 1:3]))
    
    
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
  
  filename = paste("./results/width_CIs/", "width_CIs", model_name, "_", 
                   scenario_name, ".RData", sep = "")
  
  save(width_CIs,
       file = filename)
}






################################################################################
#model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
#                "Improper1_typeIII", "Improper1_typeIV",
#                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
#                "Improper2_typeIII", "Improper2_typeIV",
#                "proper1_noInt", "proper1_onlyInt", "proper1_full", "proper1_iid",
#                "proper2_noInt", "proper2_onlyInt", "proper2_full", "proper2_iid")

model_names = c("Improper1_typeIV", "Improper2_typeI",
                "proper1_full", "proper2_onlyInt")


### Iterate over each model for each scenario on ADM1 to calculate the model choice data frames
for(model_name in model_names){
  for(scenario_name in scenario_names_ADM1){
    ## Calculate the model choice
    calc_width_CI_one_scenario_ADM1(model_name = model_name,
                                                  scenario_name = scenario_name,
                                                  n_sim)
  }
}


### Iterate over each model for each scenario on ADM4 to calculate the model choice data frames
for(model_name in model_names){
  for(scenario_name in scenario_names_ADM4){
    calc_width_CI_one_scenario_ADM4(model_name = model_name,
                                                  scenario_name = scenario_name,
                                                  n_sim)
  }
}

### Iterate over each model for each scenario on new_map to calculate the model choice data frames
for(model_name in model_names){
  for(scenario_name in scenario_names_new){
    calc_width_CI_one_scenario_new(model_name = model_name,
                                                 scenario_name = scenario_name,
                                                 n_sim)
  }
}
