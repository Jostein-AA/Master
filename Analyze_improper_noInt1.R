#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
# Create formulas

## Specify priors for hyperparameters of improper models
#---
### Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

### Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))
#---

## Specify precision matrices
#---
### Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)

### Make precision matrix for Besag on germany_map2
matrix4inla <- nb2mat(nb_first_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_first_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

## Make precision matrix for Besag on second-level
matrix4inla <- nb2mat(nb_second_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse
#---

## Specify base-formula on ADM1
base_formula_first_level <- sampled_counts ~ 1 + f(time_id, 
                                                   model = 'bym2',
                                                   scale.model = T, 
                                                   constr = T, 
                                                   rankdef = 1,
                                                   graph = RW1_prec,
                                                   hyper = temporal_hyper) + 
  f(area_id, 
    model = 'bym2',
    scale.model = T,
    constr = T,
    rankdef = 1,
    graph = Besag_prec_first_level,
    hyper = spatial_hyper)


## Specify base-formula on ADM4
base_formula_second_level <- sampled_counts ~ 1 + f(time_id, 
                                                    model = 'bym2',
                                                    scale.model = T, 
                                                    constr = T, 
                                                    rankdef = 1,
                                                    graph = RW1_prec,
                                                    hyper = temporal_hyper) + 
  f(area_id, 
    model = 'bym2',
    scale.model = T,
    constr = T,
    rankdef = 1,
    graph = Besag_prec_second_level,
    hyper = spatial_hyper)


################################################################################
# Try-catch system


## Should probably do a try-catch something for each fit, and also save eventual failures
tryCatch_inla <- function(data, formula) {
  tryCatch(
    {
      inla(formula, 
           data = data, 
           family = "poisson",
           E = E_it, 
           control.predictor = list(compute = TRUE),       #For predictions
           control.compute = list(config = TRUE, # To see constraints later
                                  cpo = T,       # For model selection
                                  return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      
      
    },
    error = function(cond) {
      print("!Error!")
      # Choose a return value in case of error
      NA
    },
    warning = function(cond) {
      print("!warning!")
      #message(conditionMessage(cond))
      # Choose a return value in case of warning
      NULL
    },
    finally = {
      print("[_/(^::^)--|")
    }
  )
  
}

################################################################################
# Fitting the data, storing the results, and updating the tracker

## Load in newest data set

#get_first_not_yet_analyzed
#get_csv_tracker_filename

## Set the last three years counts to NA for the fit

### For every successful analysis - must update the corresponding tracker-csv

### For an unsuccessful analysis - also update

## Store the linear-predictor density and CPO-values?






