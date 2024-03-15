#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
# Functions for calculating MSE, IS,...

count_mse_one_year_one_dataset <- function(){
  
}

count_mean_mse_one_dataset <- function(){
  
}

count_IS_one_year_one_dataset <- function(){
  
}

count_mean_IS_one_dataset <- function(){
  
}

rate_mse_one_year_one_dataset <- function(){
  
}

rate_mean_mse_one_dataset <- function(){
  
}

rate_IS_one_year_one_dataset <- function(){
  
}

rate_mean_IS_one_dataset <- function(){
  
}
###########################################
# Function for loading in all the data for one model and one scenario, that calculates everything into one df
calc_model_choice <- function(model_name, scenario_name, n_sim){
  file_path = paste("./Simulated_data/", scenario_name, "/", 
                    model_name, "_", scenario_name, "_", 
                    sep = "")
  
  file_ending = ".RData"
  for(i in 1:n_sim){
    file_name = paste(file_path, i, file_ending, sep = "")
    
    load(file_name)
    
    ## Calculate the MSE, IS of this data set
    
  }
  
  ## Calculate the total MSE, IS
  
  
}
calc_model_choice("Improper1_noInt", "sc1")

################################################################################
# load in data and start calculating

load("Simulated_data/sc1/Improper1_noInt_sc1_1.RData")

