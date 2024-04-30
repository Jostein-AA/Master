#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
#Load libraries
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)
library(ggpubr)
library(pspline)
library(readxl)
library(OneR)
library(mgcv)
library(splines)
library(ggfortify)
library(crs)
#library(magick)
library(INLA)
#library(bigDM)
library(MASS)
#library(RColorBrewer)
#library(tmap)


library("makemyprior")

source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
## Specify precision matrices
#---
### Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)

#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))

### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_second_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

# get scaled Besag
scaled_besag_prec_second_level <- INLA::inla.scale.model(Besag_prec_second_level, 
                                                         constr = list(A = matrix(1,1,dim(Besag_prec_second_level)[1]), 
                                                                       e = 0))

INLA::inla.write.graph(
  INLA::inla.matrix2graph(RW1_prec),
  filename = "RW1_prec.graph"
)
INLA::inla.write.graph(
  INLA::inla.matrix2graph(Besag_prec_second_level),
  filename = "Besag_prec_second_level.graph"
)


RW1_prec_path <- paste0(getwd(), "/RW1_prec.graph")
Besag_prec_second_level_path <- paste0(getwd(), "/Besag_prec_second_level.graph")


#Get precision matric for type II interaction by Kronecker product
typeII_prec_second_level <- scaled_RW_prec %x% diag(nrow(second_level_admin_map))

#Get precision matric for type II interaction by Kronecker product
typeIII_prec_second_level <- diag(tT) %x% scaled_besag_prec_second_level

#Get type IV interaction precision matrix
typeIV_prec_second_level <- scaled_RW_prec %x% scaled_besag_prec_second_level

#---
# Get the sum-to-zero constraints

#Get sum-to-zero constraints for type II interaction
typeII_constraints_second_level = constraints_maker(type = "II", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)

#Get sum-to-zero constraints for type II interaction
typeIII_constraints_second_level = constraints_maker(type = "III", 
                                                     n = nrow(second_level_admin_map), 
                                                     t = tT)

#Get sum-to-zero constraints for type IV interaction
typeIV_constraints_second_level = constraints_maker(type = "IV", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)

#---
################################################################################
# Load in the separate temporal, spatial, and spatiotemporal priors,
# AND make formula

prior_data_bym2_time = readRDS("./typeV/prior_data_bym2_time_ADM4.rds")
prior_data_bym2_space = readRDS("./typeV/prior_data_bym2_space_ADM4.rds")
prior_data_interaction = readRDS("./typeV/prior_data_interaction_ADM4.rds")

# obs! merk at rekkefoelgen paa inla-formelen maa matche rekkefoelgen vi bruker log-presisjonene her!
prior_func <- function(logprec){
 return(
   eval_joint_prior(-c(logprec[1:2]), prior_data_bym2_time) +
     eval_joint_prior(-c(logprec[3:4]), prior_data_bym2_space) +
     eval_joint_prior(-c(logprec[5:8]), prior_data_interaction)
 )
}




################################################################################

## Should probably do a try-catch something for each fit, and also save eventual failures
tryCatch_inla <- function(data,
                          data_set_id,
                          csv_tracker_filename,
                          model_name, scenario_name) {
  tryCatch(
    {
      inla.setOption(inla.timeout = 1500) # Set upper-time limit to 1500 sec (25 minutes) 
      
      tmp_ = inla(, 
                  data = data, 
                  family = "poisson",
                  E = E_it, #E_it
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
      
      marginals = tmp_$marginals.fitted.values[(n_ADM4 * 10 + 1):(n_ADM4 * 13)]
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
# SC2
model_name = "Improper1_typeIV"
scenario_name = "sc2"

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
  
  ### Load in sc2 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                "space.time")]
  
  lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla <- tryCatch_inla(lambda_sc.df,
                               data_set_id,
                               csv_tracker_filename,
                               model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC4
model_name = "Improper1_typeIV"
scenario_name = "sc4"

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
  
  ### Load in simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                "space.time")]
  
  lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla <- tryCatch_inla(lambda_sc.df,
                               data_set_id,
                               csv_tracker_filename,
                               model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC6
model_name = "Improper1_typeIV"
scenario_name = "sc6"

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
  
  ### Load in simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                "space.time")]
  
  lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla <- tryCatch_inla(lambda_sc.df,
                               data_set_id,
                               csv_tracker_filename,
                               model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC8
model_name = "Improper1_typeIV"
scenario_name = "sc8"

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
  
  ### Load in simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                "space.time")]
  
  lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla <- tryCatch_inla(lambda_sc.df,
                               data_set_id,
                               csv_tracker_filename,
                               model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC10
model_name = "Improper1_typeIV"
scenario_name = "sc10"

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
  
  ### Load in simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                "space.time")]
  
  lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla <- tryCatch_inla(lambda_sc.df,
                               data_set_id,
                               csv_tracker_filename,
                               model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC12
model_name = "Improper1_typeIV"
scenario_name = "sc12"

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
  
  ### Load in simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                "space.time")]
  
  lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla <- tryCatch_inla(lambda_sc.df,
                               data_set_id,
                               csv_tracker_filename,
                               model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))













