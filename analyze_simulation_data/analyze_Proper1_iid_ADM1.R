#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
# Create formulas

## Specify priors for hyperparameters of proper models
#---
### Temporal hyperparameters (prec. of AR1 and AR1's mixing param) w. corresponding priors: penalized constraint 
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), 
                 rho = list(prior = 'pc.cor1', 
                            param = c(0.5, 0.5 + 1E-2)))


### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01)))

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))

#---

## Specify precision matrices
#---
### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_first_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_first_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

#---


proper_1_iid_ADM1 <- sampled_counts ~ 1 + time_id +
                                    f(time_id.copy,
                                      model = "ar1",
                                      hyper = ar1_hyper) + 
                                    f(area_id, 
                                      model = "besagproper2",
                                      graph = Besag_prec_first_level,
                                      hyper = spatial_hyper) + 
                                    f(space.time,
                                      model = "iid", 
                                      hyper = interaction_hyper )



################################################################################

## Should probably do a try-catch something for each fit, and also save eventual failures
tryCatch_inla <- function(data,
                          data_set_id,
                          csv_tracker_filename,
                          model_name, scenario_name) {
  tryCatch(
    {
      tmp_ = inla(proper_1_iid_ADM1, 
                  data = data, 
                  family = "poisson",
                  E = E_it, #E_it
                  verbose = T,
                  control.predictor = list(compute = TRUE,
                                           link = 1),       #For predictions
                  #control.family = list(control.link = list(model = "log")),
                  control.compute = list(config = TRUE, # To see constraints later
                                         cpo = T,       # For model selection
                                         return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      
      if(tmp_$ok == FALSE){ ## INLA has crashed
        # Update tracker
        tracker.df <- read.csv(csv_tracker_filename)
        tracker.df[data_set_id, ]$error = data_set_id
        write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      }
      
      if(tmp_$mode$mode.status > 0){ ## Potentially something weird with the mode of hyperparameters
        # Update tracker
        tracker.df <- read.csv(csv_tracker_filename)
        tracker.df[data_set_id, ]$warning = data_set_id
        write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      }
      
      
      ### Save the linear predictor-marginal distribution and the CPO-values
      filename_to_save <- paste("./results/", model_name, "/", scenario_name, "/", 
                                model_name, "_", scenario_name, "_", toString(data_set_id), ".RData", 
                                sep = "")
      
      marginals = sort_proper_fitted(tmp_$marginals.fitted.values, n_ADM1, tT) 
      cpo = tmp_$cpo$cpo
      
      save(marginals, 
           cpo,
           file = filename_to_save)
    },
    error = function(cond) {
      print("!Error!")
      
      # Update tracker
      tracker.df <- read.csv(csv_tracker_filename)
      tracker.df[data_set_id, ]$error = data_set_id
      write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      
      # Choose a return value in case of error
      -1
    },
    warning = function(cond) {
      print("!warning!")
      print(cond)
      
      # Update tracker
      tracker.df <- read.csv(csv_tracker_filename)
      tracker.df[data_set_id, ]$warning = data_set_id
      write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      
      print("---")
    },
    finally = {
      print(paste("[_/(^::^)--|", model_name, scenario_name, 
                  toString(data_set_id), sep = " "))
      
      # Update tracker
      tracker.df <- read.csv(csv_tracker_filename)
      tracker.df[data_set_id, ]$analyzed = data_set_id
      write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      
    }
  )
}



################################################################################
# SC1
model_name = "proper1_iid"
scenario_name = "sc1"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
  
}

print("Do not run as is on Markov ya bastard")

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC3
model_name = "proper1_iid"
scenario_name = "sc3"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC5
model_name = "proper1_iid"
scenario_name = "sc5"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC7
model_name = "proper1_iid"
scenario_name = "sc7"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC9
model_name = "proper1_iid"
scenario_name = "sc9"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC11
model_name = "proper1_iid"
scenario_name = "sc11"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
  
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))