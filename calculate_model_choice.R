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
  model_choice_for_counts <- data.frame(mse_1_year_ahead = rep(NA, n_sim), mse_2_year_ahead = rep(NA, n_sim),
                                         mse_3_year_ahead = rep(NA, n_sim), total_mse = rep(NA, n_sim), 
                                         IS_1_year_ahead = rep(NA, n_sim), IS_2_year_ahead = rep(NA, n_sim), 
                                         IS_3_year_ahead = rep(NA, n_sim), total_IS = rep(NA, n_sim))
  
  model_choice_for_rates <- data.frame(mse_1_year_ahead = rep(NA, n_sim), mse_2_year_ahead = rep(NA, n_sim),
                                       mse_3_year_ahead = rep(NA, n_sim), total_mse = rep(NA, n_sim), 
                                       IS_1_year_ahead = rep(NA, n_sim), IS_2_year_ahead = rep(NA, n_sim), 
                                       IS_3_year_ahead = rep(NA, n_sim), total_IS = rep(NA, n_sim))
  
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
    
    ### For the proper models on the ADM1, the marginals must be sorted
    ### This was done by hardcoding this in for only the proper models
    # marginals = sort_proper_fitted(marginals, n_ADM1, tT)
    
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
      
      ## Extract the sampled rates for this year
      one_year_sampled_rates = lambda_df$lambda_it[((year - 1) * n_ADM1 + 1):(year * n_ADM1)]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      model_choice_for_counts[i, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(one_year_sampled_counts,
                                                                                                one_year_margs,
                                                                                                E_it)
      
      ## MSE for rate in year ahead (1, 2, or 3 years ahead)
      model_choice_for_rates[i, year - years_pred_on[1] + 1] = rate_mse_one_year_one_dataset(one_year_sampled_rates, 
                                                                                             one_year_margs)
      
      
      ## IS for count in year ahead (1, 2, or 3 years ahead)
      model_choice_for_counts[i, year - years_pred_on[1] + 5] =  count_IS_one_year_one_dataset(one_year_sampled_counts, 
                                                                                               one_year_margs, 
                                                                                               E_it)
      
      
      ## IS for rate in year ahead (1, 2, or 3 years ahead)
      model_choice_for_rates[i, year - years_pred_on[1] + 5] = rate_IS_one_year_one_dataset(one_year_sampled_rates,
                                                                                            one_year_margs)
      
    }
    
    #Get the total MSE for count
    model_choice_for_counts[i, 4] = mean(as.numeric(model_choice_for_counts[i, 1:3]))
    
    #Get the total MSE for rate
    model_choice_for_rates[i, 4] = mean(as.numeric(model_choice_for_rates[i, 1:3]))
    
    #Get the total IS for count
    model_choice_for_counts[i, 8] = mean(as.numeric(model_choice_for_counts[i, 5:7]))
    
    #Get the total IS for rate
    model_choice_for_rates[i, 8] = mean(as.numeric(model_choice_for_rates[i, 5:7]))
    
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
                   scenario_name, ".RData", sep = "")
  
  save(model_choice_for_counts,
       model_choice_for_rates,
       file = filename)
}


## Function that takes one model, one scenario, and calculates all the MSEs and IS' for ADM1
calc_mse_is_count_one_model_one_scenario_ADM4 <- function(model_name,
                                                          scenario_name,
                                                          n_sim){
  
  
  ### Initialize data frame to store the mse's and is's for each data set
  model_choice_for_counts <- data.frame(mse_1_year_ahead = rep(NA, n_sim), mse_2_year_ahead = rep(NA, n_sim),
                             mse_3_year_ahead = rep(NA, n_sim), total_mse = rep(NA, n_sim), 
                             IS_1_year_ahead = rep(NA, n_sim), IS_2_year_ahead = rep(NA, n_sim), 
                             IS_3_year_ahead = rep(NA, n_sim), total_IS = rep(NA, n_sim))
  
  model_choice_for_rates <- data.frame(mse_1_year_ahead = rep(NA, n_sim), mse_2_year_ahead = rep(NA, n_sim),
                                        mse_3_year_ahead = rep(NA, n_sim), total_mse = rep(NA, n_sim), 
                                        IS_1_year_ahead = rep(NA, n_sim), IS_2_year_ahead = rep(NA, n_sim), 
                                        IS_3_year_ahead = rep(NA, n_sim), total_IS = rep(NA, n_sim))
  
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
      
      ## Extract the sampled rates for this year
      one_year_sampled_rates = lambda_df$lambda_it[((year - 1) * n_ADM4 + 1):(year * n_ADM4)]
      
      ## MSE for count in year ahead (1, 2, or 3 years ahead)
      model_choice_for_counts[i, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(one_year_sampled_counts,
                                                                                                one_year_margs,
                                                                                                E_it)
      ## MSE for rate in year ahead (1, 2, or 3 years ahead)
      model_choice_for_rates[i, year - years_pred_on[1] + 1] = rate_mse_one_year_one_dataset(one_year_sampled_rates, 
                                                                                             one_year_margs)
      
      ## IS for count in year ahead (1, 2, or 3 years ahead)
      model_choice_for_counts[i, year - years_pred_on[1] + 5] =  count_IS_one_year_one_dataset(one_year_sampled_counts, 
                                                                                               one_year_margs, 
                                                                                               E_it)
      
      ## IS for rate in year ahead (1, 2, or 3 years ahead)
      model_choice_for_rates[i, year - years_pred_on[1] + 5] = rate_IS_one_year_one_dataset(one_year_sampled_rates,
                                                                                            one_year_margs)
      
    }
    
    #Get the total MSE for count
    model_choice_for_counts[i, 4] = mean(as.numeric(model_choice_for_counts[i, 1:3]))
    
    #Get the total MSE for rate
    model_choice_for_rates[i, 4] = mean(as.numeric(model_choice_for_rates[i, 1:3]))
    
    #Get the total IS for count
    model_choice_for_counts[i, 8] = mean(as.numeric(model_choice_for_counts[i, 5:7]))
    
    #Get the total IS for rate
    model_choice_for_rates[i, 8] = mean(as.numeric(model_choice_for_rates[i, 5:7]))
    
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
                   scenario_name, ".RData", sep = "")
  
  save(model_choice_for_counts,
       model_choice_for_rates,
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


################################################################################
# Calculate the totals and make a table for better vieweng of results

#####
# ADM1 results


library(latex2exp)
library(tables)

# The count results
for(scenario_name in scenario_names_ADM1){
  print(scenario_name)
  
  model_choice_counts.df <- data.frame(Model = rep(c("Improper1_noInt",
                                              "Improper1_typeI",
                                              "Improper1_typeII",
                                              "Improper1_typeIII",
                                              "Improper1_typeIV",
                                              "Improper2_noInt",
                                              "Improper2_typeI",
                                              "Improper2_typeII",
                                              "Improper2_typeIII",
                                              "Improper2_typeIV",
                                              "proper1_noInt",
                                              "proper1_onlyInt",
                                              "proper1_full",
                                              "proper2_noInt",
                                              "proper2_onlyInt",
                                              "proper2_full"), 8),
                                model_choice = c(rep(1, 16),rep(2, 16),
                                                 rep(3, 16),rep(4, 16),
                                                 rep(5, 16),rep(6, 16),
                                                 rep(7, 16),rep(8, 16)),
                                value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = ""))
    tmp_ <- model_choice_for_counts
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Calculate the average over each data set
    tmp2_ <- as.numeric(colMeans(tmp_))
    
    #Insert values into the model_choice_counts.df
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_counts.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_counts_", scenario_name, ".tex", sep = ""))
}

# The rate results
for(scenario_name in scenario_names_ADM1){
  print(scenario_name)
  
  model_choice_rates.df <- data.frame(Model = rep(c("Improper1_noInt",
                                                     "Improper1_typeI",
                                                     "Improper1_typeII",
                                                     "Improper1_typeIII",
                                                     "Improper1_typeIV",
                                                     "Improper2_noInt",
                                                     "Improper2_typeI",
                                                     "Improper2_typeII",
                                                     "Improper2_typeIII",
                                                     "Improper2_typeIV",
                                                     "proper1_noInt",
                                                     "proper1_onlyInt",
                                                     "proper1_full",
                                                     "proper2_noInt",
                                                     "proper2_onlyInt",
                                                     "proper2_full"), 8),
                                       model_choice = c(rep(1, 16),rep(2, 16),
                                                        rep(3, 16),rep(4, 16),
                                                        rep(5, 16),rep(6, 16),
                                                        rep(7, 16),rep(8, 16)),
                                       value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = ""))
    tmp_ <- model_choice_for_rates
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Calculate the average over each data set
    tmp2_ <- as.numeric(colMeans(tmp_))
    
    #Insert values into the model_choice_rates.df
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_rates.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_rates_", scenario_name, ".tex", sep = ""))
}


#####
# ADM4 results

model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
                "Improper1_typeIII", "Improper1_typeIV",
                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
                "Improper2_typeIII", "Improper2_typeIV",
                "proper1_noInt", "proper1_onlyInt", "proper1_full",
                "proper2_noInt", "proper2_onlyInt", "proper2_full")

# For counts
for(scenario_name in scenario_names_ADM4){
  print(scenario_name)
  
  model_choice_counts.df <- data.frame(Model = rep(c("Improper1_noInt",
                                              "Improper1_typeI",
                                              "Improper1_typeII",
                                              "Improper1_typeIII",
                                              "Improper1_typeIV",
                                              "Improper2_noInt",
                                              "Improper2_typeI",
                                              "Improper2_typeII",
                                              "Improper2_typeIII",
                                              "Improper2_typeIV",
                                              "proper1_noInt",
                                              "proper1_onlyInt",
                                              "proper1_full",
                                              "proper2_noInt",
                                              "proper2_onlyInt",
                                              "proper2_full"), 8),
                                model_choice = c(rep(1, 16),rep(2, 16),
                                                 rep(3, 16),rep(4, 16),
                                                 rep(5, 16),rep(6, 16),
                                                 rep(7, 16),rep(8, 16)),
                                value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = " "))
    tmp_ <- model_choice_for_counts
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Calculate the average over each data set
    tmp2_ <- as.numeric(colMeans(tmp_))
    
    #Insert values into the model_choice_counts.df
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_counts.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_counts_", scenario_name, ".tex", sep = ""))
}

# For rates
for(scenario_name in scenario_names_ADM4){
  print(scenario_name)
  
  model_choice_rates.df <- data.frame(Model = rep(c("Improper1_noInt",
                                                     "Improper1_typeI",
                                                     "Improper1_typeII",
                                                     "Improper1_typeIII",
                                                     "Improper1_typeIV",
                                                     "Improper2_noInt",
                                                     "Improper2_typeI",
                                                     "Improper2_typeII",
                                                     "Improper2_typeIII",
                                                     "Improper2_typeIV",
                                                     "proper1_noInt",
                                                     "proper1_onlyInt",
                                                     "proper1_full",
                                                     "proper2_noInt",
                                                     "proper2_onlyInt",
                                                     "proper2_full"), 8),
                                       model_choice = c(rep(1, 16),rep(2, 16),
                                                        rep(3, 16),rep(4, 16),
                                                        rep(5, 16),rep(6, 16),
                                                        rep(7, 16),rep(8, 16)),
                                       value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = " "))
    tmp_ <- model_choice_for_rates
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Calculate the average over each data set
    tmp2_ <- as.numeric(colMeans(tmp_))
    
    #Insert values into the model_choice_rates.df
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_rates.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_rates_", scenario_name, ".tex", sep = ""))
}

################################################################################
# Get a table for a singular dataset for each scenario, to better see maybe?

#####
# ADM1 results

# The count results
for(scenario_name in scenario_names_ADM1){
  print(scenario_name)
  
  model_choice_counts.df <- data.frame(Model = rep(c("Improper1_noInt",
                                                     "Improper1_typeI",
                                                     "Improper1_typeII",
                                                     "Improper1_typeIII",
                                                     "Improper1_typeIV",
                                                     "Improper2_noInt",
                                                     "Improper2_typeI",
                                                     "Improper2_typeII",
                                                     "Improper2_typeIII",
                                                     "Improper2_typeIV",
                                                     "proper1_noInt",
                                                     "proper1_onlyInt",
                                                     "proper1_full",
                                                     "proper2_noInt",
                                                     "proper2_onlyInt",
                                                     "proper2_full"), 8),
                                       model_choice = c(rep(1, 16),rep(2, 16),
                                                        rep(3, 16),rep(4, 16),
                                                        rep(5, 16),rep(6, 16),
                                                        rep(7, 16),rep(8, 16)),
                                       value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = ""))
    tmp_ <- model_choice_for_counts
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    # Extract first data set
    tmp2_ = tmp_[1, ]
    
    #Insert values into the model_choice_counts.df
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_counts.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_counts_singular_dataset_", scenario_name, ".tex", sep = ""))
}

# The rate results
for(scenario_name in scenario_names_ADM1){
  print(scenario_name)
  
  model_choice_rates.df <- data.frame(Model = rep(c("Improper1_noInt",
                                                    "Improper1_typeI",
                                                    "Improper1_typeII",
                                                    "Improper1_typeIII",
                                                    "Improper1_typeIV",
                                                    "Improper2_noInt",
                                                    "Improper2_typeI",
                                                    "Improper2_typeII",
                                                    "Improper2_typeIII",
                                                    "Improper2_typeIV",
                                                    "proper1_noInt",
                                                    "proper1_onlyInt",
                                                    "proper1_full",
                                                    "proper2_noInt",
                                                    "proper2_onlyInt",
                                                    "proper2_full"), 8),
                                      model_choice = c(rep(1, 16),rep(2, 16),
                                                       rep(3, 16),rep(4, 16),
                                                       rep(5, 16),rep(6, 16),
                                                       rep(7, 16),rep(8, 16)),
                                      value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = ""))
    tmp_ <- model_choice_for_rates
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Extract first data set
    tmp2_ <- tmp_[1, ]
    
    #Insert values into the model_choice_rates.df
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_rates.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_rates_singular_dataset_", scenario_name, ".tex", sep = ""))
}


#####
# ADM4 results

model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
                "Improper1_typeIII", "Improper1_typeIV",
                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
                "Improper2_typeIII", "Improper2_typeIV",
                "proper1_noInt", "proper1_onlyInt", "proper1_full",
                "proper2_noInt", "proper2_onlyInt", "proper2_full")

# For counts
for(scenario_name in scenario_names_ADM4){
  print(scenario_name)
  
  model_choice_counts.df <- data.frame(Model = rep(c("Improper1_noInt",
                                                     "Improper1_typeI",
                                                     "Improper1_typeII",
                                                     "Improper1_typeIII",
                                                     "Improper1_typeIV",
                                                     "Improper2_noInt",
                                                     "Improper2_typeI",
                                                     "Improper2_typeII",
                                                     "Improper2_typeIII",
                                                     "Improper2_typeIV",
                                                     "proper1_noInt",
                                                     "proper1_onlyInt",
                                                     "proper1_full",
                                                     "proper2_noInt",
                                                     "proper2_onlyInt",
                                                     "proper2_full"), 8),
                                       model_choice = c(rep(1, 16),rep(2, 16),
                                                        rep(3, 16),rep(4, 16),
                                                        rep(5, 16),rep(6, 16),
                                                        rep(7, 16),rep(8, 16)),
                                       value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = " "))
    tmp_ <- model_choice_for_counts
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Extract first dataset
    tmp2_ <- tmp_[1, ]
    
    #Insert values into the model_choice_counts.df
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_counts.df[model_choice_counts.df$Model == model_name & model_choice_counts.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_counts.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_counts_singular_dataset_", scenario_name, ".tex", sep = ""))
}

# For rates
for(scenario_name in scenario_names_ADM4){
  print(scenario_name)
  
  model_choice_rates.df <- data.frame(Model = rep(c("Improper1_noInt",
                                                    "Improper1_typeI",
                                                    "Improper1_typeII",
                                                    "Improper1_typeIII",
                                                    "Improper1_typeIV",
                                                    "Improper2_noInt",
                                                    "Improper2_typeI",
                                                    "Improper2_typeII",
                                                    "Improper2_typeIII",
                                                    "Improper2_typeIV",
                                                    "proper1_noInt",
                                                    "proper1_onlyInt",
                                                    "proper1_full",
                                                    "proper2_noInt",
                                                    "proper2_onlyInt",
                                                    "proper2_full"), 8),
                                      model_choice = c(rep(1, 16),rep(2, 16),
                                                       rep(3, 16),rep(4, 16),
                                                       rep(5, 16),rep(6, 16),
                                                       rep(7, 16),rep(8, 16)),
                                      value = 1:(8 * 16))
  
  
  for(model_name in model_names){
    #Load in the model_choice data
    load(paste("./results/model_choice/model_choice_", model_name, "_", scenario_name, ".RData",
               sep = " "))
    tmp_ <- model_choice_for_rates
    
    #Remove potential NAs
    tmp_ <- na.omit(tmp_)
    
    #Extract first data set
    tmp2_ <- tmp_[1, ]
    
    #Insert values into the model_choice_rates.df
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 1, ]$value = tmp2_[1]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 2, ]$value = tmp2_[2]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 3, ]$value = tmp2_[3]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 4, ]$value = tmp2_[4]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 5, ]$value = tmp2_[5]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 6, ]$value = tmp2_[6]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 7, ]$value = tmp2_[7]
    model_choice_rates.df[model_choice_rates.df$Model == model_name & model_choice_rates.df$model_choice == 8, ]$value = tmp2_[8]
    
  }
  
  #Make caption and label for latex table
  caption = "HEIHEI"
  label = paste("model-choice-", scenario_name, sep = "")
  
  #Make latex table
  latex_tabular <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, 
                               nopagebreak = "\\hline",
                               spacing = 0)~
      Heading()*Factor(model_choice, 
                       levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (total)",
                                      "IS (1)", "IS (2)", "IS (3)", "IS (total)"))*
      Heading()*value*Heading()*identity,
    data = model_choice_rates.df),
    caption = caption,
    label = label
  )
  
  #Save latex table
  cat(latex_tabular, file =paste("./results/model_choice/table_rates_singular_dataset_", scenario_name, ".tex", sep = ""))
}




