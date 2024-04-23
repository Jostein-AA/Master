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
tT = 13


#model_names = c("proper1_noInt", "proper1_onlyInt", "proper1_full", "proper1_iid",
#                "proper2_noInt", "proper2_onlyInt", "proper2_full", "proper2_iid")


#scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")

################################################################################
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

################################################################################

for(model_name in model_names){
  for(scenario_name in scenario_names_ADM1){
    print(scenario_name)
    for(i in 1:n_sim){
      print(i)
      predictions_dir = paste("./results/", model_name, "/",
                              scenario_name, sep = "")
      
      file_base_name = paste("/", model_name, "_", scenario_name, "_", sep = "")
      file_ending = ".RData"
      
      marginals = read_file(paste(predictions_dir, file_base_name, toString(i), file_ending,
                      sep = ""))
      
      if(is.null(marginals)){
        next
      }
      
      # Sort the marginals, so that I do not have to deal w. this bull
      marginals = sort_proper_fitted(marginals, n_ADM1, tT)
      
      save(marginals,
           file = paste(predictions_dir, file_base_name, toString(i), file_ending,
                  sep = ""))
    }
  }
}

model_name = "Improper1_typeII"
scenario_names = c("sc6", "sc8", "sc10", "sc12")

for(scenario_name in scenario_names){
  for(i in 1:n_sim){
    print(i)
    predictions_dir = paste("./results/", model_name, "/",
                            scenario_name, sep = "")
    
    file_base_name = paste("/", model_name, "_", scenario_name, "_", sep = "")
    file_ending = ".RData"
    
    marginals = read_file(paste(predictions_dir, file_base_name, toString(i), file_ending,
                                sep = ""))
    
    if(is.null(marginals)){
      next
    }
    
    if(length(marginals) == 13 * nrow(second_level_admin_map)){
      marginals = marginals[(n_ADM4 * 10 + 1):(n_ADM4 * 13)]
    }
    
    
    save(marginals,
         file = paste(predictions_dir, file_base_name, toString(i), file_ending,
                      sep = ""))
  }
}







