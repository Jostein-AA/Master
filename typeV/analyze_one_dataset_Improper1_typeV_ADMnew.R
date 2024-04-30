#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("makemyprior")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

#---
RW1_prec_path <- paste0(getwd(), "/RW1_prec.graph")
Besag_prec_new_level_path <- paste0(getwd(), "/Besag_prec_new_level.graph")


# Load in the prior
typeV_prior_ADMnew = readRDS("./typeV/typeV_prior_ADMnew.rds")

plot_prior(typeV_prior_ADMnew)
plot_tree_structure(typeV_prior_ADMnew)

### Implement on one data set
scenario_name = "sc14"

load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", 
                         "space.time")]

lambda_$sampled_counts = lambda.df$sampled_counts[, 5]
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

typeV_prior_ADMnew$data$random$time_id_iid = lambda_$time_id
typeV_prior_ADMnew$data$random$time_id_struct = lambda_$time_id
typeV_prior_ADMnew$data$random$area_id_iid = lambda_$area_id
typeV_prior_ADMnew$data$random$area_id_struct = lambda_$area_id
typeV_prior_ADMnew$data$random$space.time_I = lambda_$space.time
typeV_prior_ADMnew$data$random$space.time_II = lambda_$space.time
typeV_prior_ADMnew$data$random$space.time_III = lambda_$space.time
typeV_prior_ADMnew$data$random$space.time_IV = lambda_$space.time

typeV_prior_ADMnew$response = lambda_$sampled_counts


### Function to overwrite
#trace(FUN, edit= T)


#Issue in makemyprior:::make_inla_formula, call ->
### debug(makemyprior:::make_inla_formula)
### Then call only the inference_inla(...) below
# To see the issue. Issue is there are in prior_obj$data 
# There are no fixed or random... Meaning the inla_formula becomes 'y ~ '

makemyprior_inla_obj <- inference_inla(#formula = typeV_formula,
                                       typeV_prior_ADMnew, 
                                       #data = lambda_, 
                                       #verbose = T,
                                       E = lambda_$E_it, 
                                       control.predictor = list(compute = TRUE,
                                                                link = 1),       #For predictions
                                       control.compute = list(config = TRUE, # To see constraints later
                                                              return.marginals.predictor=TRUE))


#Lag prior for kun tid, kun for rom, kun for rom-tid


evaL_joint_prior

#plot(makemyprior_inla_obj$inla)
#makemyprior_inla_obj$inla$marginals.fitted.values






