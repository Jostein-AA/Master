#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("makemyprior")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
# Create formulas

## Specify precision matrices
#---
### Specify the RW1 precision matrix
# RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
#                            scale.model = FALSE, sparse = TRUE)
# 
# #Scale precision matrix of RW model so the geometric mean of the marginal variances is one
# scaled_RW_prec <- inla.scale.model(RW1_prec,
#                                    list(A = matrix(1, 1, dim(RW1_prec)[1]),
#                                         e = 0))
# 
# ### Make precision matrix for Besag on ADM1
# matrix4inla <- nb2mat(nb_second_level, style="B")
# mydiag = rowSums(matrix4inla)
# matrix4inla <- -matrix4inla
# diag(matrix4inla) <- mydiag
# Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse
# 
# # get scaled Besag
# scaled_besag_prec_second_level <- INLA::inla.scale.model(Besag_prec_second_level, 
#                                                         constr = list(A = matrix(1,1,dim(Besag_prec_second_level)[1]), 
#                                                                       e = 0))
# 
# 
# #Get precision matric for type II interaction by Kronecker product
# typeII_prec_second_level <- scaled_RW_prec %x% diag(nrow(second_level_admin_map))
# 
# #Get precision matric for type II interaction by Kronecker product
# typeIII_prec_second_level <- diag(tT) %x% scaled_besag_prec_second_level
# 
# #Get type IV interaction precision matrix
# typeIV_prec_second_level <- scaled_RW_prec %x% scaled_besag_prec_second_level
# 
# #---
# # Get the sum-to-zero constraints
# 
# #Get sum-to-zero constraints for type II interaction
# typeII_constraints_second_level = constraints_maker(type = "II", 
#                                                     n = nrow(second_level_admin_map), 
#                                                     t = tT)
# 
# #Get sum-to-zero constraints for type II interaction
# typeIII_constraints_second_level = constraints_maker(type = "III", 
#                                                      n = nrow(second_level_admin_map), 
#                                                      t = tT)
# 
# #Get sum-to-zero constraints for type IV interaction
# typeIV_constraints_second_level = constraints_maker(type = "IV", 
#                                                    n = nrow(second_level_admin_map), 
#                                                    t = tT)

#---
#Make formula

# typeV_formula <- sampled_counts ~ 1 + mc(time_id_iid,      # Unstructured temporal effect
#                                          model = "iid") + 
#   mc(time_id_struct,   # Structured temporal effect
#      model = "besag", 
#      graph = RW1_prec_path,
#      scale.model = TRUE,
#      constr = T) + 
#   mc(area_id_iid,      # Unstructured spatial effect
#      model = "iid") + 
#   mc(area_id_struct,   # Structured spatial effect
#      model = "besag",
#      graph = Besag_prec_new_level_path,
#      scale.model = T,
#      constr = T) + 
#   mc(space.time_I,     # Unstructured space-time interaction
#      model = "iid") + 
#   mc(space.time_II,    # Structured time interacting w. iid space
#      model = "generic0",
#      Cmatrix = typeII_prec_new_level,
#      extraconstr = typeII_constraints_new_level,
#      rankdef = nrow(new_map)) +
#   mc(space.time_III,   # iid time interacting w. structured space
#      model = "generic0",
#      Cmatrix = typeIII_prec_new_level,
#      extraconstr = typeIII_constraints_new_level,
#      rankdef = tT) +
#   mc(space.time_IV,    # Structured time interacting w. structured space
#      model = "generic0",
#      Cmatrix = typeIV_prec_new_level,
#      extraconstr = typeIV_constraints_new_level)
# 
# 
# INLA::inla.write.graph(
#   INLA::inla.matrix2graph(RW1_prec),
#   filename = "RW1_prec.graph"
# )
# INLA::inla.write.graph(
#   INLA::inla.matrix2graph(Besag_prec_new_level),
#   filename = "Besag_prec_new_level.graph"
# )
# 
# 
# RW1_prec_path <- paste0(getwd(), "/RW1_prec.graph")
# Besag_prec_new_level_path <- paste0(getwd(), "/Besag_prec_new_level.graph")


# Load in the prior
typeV_prior_ADMnew = readRDS("./typeV/typeV_prior_ADMnew.rds")


################################################################################

## Should probably do a try-catch something for each fit, and also save eventual failures
tryCatch_inla <- function(data,
                          data_set_id,
                          csv_tracker_filename,
                          model_name, scenario_name) {
  tryCatch(
    {
      inla.setOption(inla.timeout = 1000) # Set upper-time limit to 1000 sec
      
      makemyprior_inla_obj <- inference_inla(prior_obj = typeV_prior_ADMnew, 
                                             data = data, 
                                             E = E_it, 
                                             control.predictor = list(compute = TRUE,
                                                                        link = 1),       #For predictions
                                             control.compute = list(config = TRUE, # To see constraints later
                                                                    return.marginals.predictor=TRUE))
      
      #return inla object: $inla
      tmp_ <- makemyprior_inla_obj$inla
        
      ### Save the linear predictor-marginal distribution and the CPO-values
      filename_to_save <- paste("./results/", model_name, "/", scenario_name, "/", 
                                model_name, "_", scenario_name, "_", toString(data_set_id), ".RData", 
                                sep = "")
      
      marginals = tmp_$marginals.fitted.values 
      
      save(marginals, 
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
# SC13
model_name = "Improper1_typeV"
scenario_name = "sc13"

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
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", 
                                 "space.time")]
  
  lambda_$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  lambda_$area_id_iid = lambda_$area_id
  lambda_$area_id_struct = lambda_$area_id
  lambda_$time_id_iid = lambda_$time_id
  lambda_$time_id_struct = lambda_$time_id
  lambda_$space.time_I = lambda_$space.time
  lambda_$space.time_II = lambda_$space.time
  lambda_$space.time_III = lambda_$space.time
  lambda_$space.time_IV = lambda_$space.time
  
  ## Set the last three years counts to NA for the fit
  lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
  
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC14
model_name = "Improper1_typeV"
scenario_name = "sc14"

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
  lambda_sc1.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                 "space.time")]
  
  lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc1.df[lambda_sc1.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_sc1.df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC15
model_name = "Improper1_typeV"
scenario_name = "sc15"

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
  lambda_sc1.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                 "space.time")]
  
  lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc1.df[lambda_sc1.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_sc1.df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC16
model_name = "Improper1_typeV"
scenario_name = "sc16"

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
  lambda_sc1.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                 "space.time")]
  
  lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc1.df[lambda_sc1.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_sc1.df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC17
model_name = "Improper1_typeV"
scenario_name = "sc17"

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
  lambda_sc1.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                 "space.time")]
  
  lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc1.df[lambda_sc1.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_sc1.df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC18
model_name = "Improper1_typeV"
scenario_name = "sc18"

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
  lambda_sc1.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                                 "space.time")]
  
  lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_sc1.df[lambda_sc1.df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_sc1.df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))
