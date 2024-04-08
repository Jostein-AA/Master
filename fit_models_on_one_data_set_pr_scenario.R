#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
# Create formulas for the Improper models

## Specify priors for hyperparameters of improper models
#---
### Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

### Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))
#---

## Specify precision matrices
#---
### Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)

### Specify the RW1 precision matrix
RW2_prec <- INLA:::inla.rw(n = tT, order = 2, 
                           scale.model = FALSE, sparse = TRUE)

### Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW1_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))

scaled_RW2_prec <- inla.scale.model(RW2_prec,
                                   list(A = matrix(1, 1, dim(RW2_prec)[1]),
                                        e = 0))

### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_first_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_first_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

### get scaled Besag on ADM1
scaled_besag_prec_first_level <- INLA::inla.scale.model(Besag_prec_first_level, 
                                                        constr = list(A = matrix(1,1,dim(Besag_prec_first_level)[1]), 
                                                                      e = 0))

### Make precision matrix for Besag on ADM4
matrix4inla <- nb2mat(nb_second_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

### get scaled Besag on ADM4
scaled_besag_prec_second_level <- INLA::inla.scale.model(Besag_prec_second_level, 
                                                        constr = list(A = matrix(1,1,dim(Besag_prec_second_level)[1]), 
                                                                      e = 0))

#---

## Specify base-formula on ADM1
### w. RW1
base_formula_1_first_level <- sampled_counts ~ 1 + f(time_id, 
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

### w. RW2
base_formula_2_first_level <- sampled_counts ~ 1 + f(time_id, 
                                                   model = 'bym2',
                                                   scale.model = T, 
                                                   constr = T, 
                                                   graph = RW2_prec,
                                                   hyper = temporal_hyper) + 
  f(area_id, 
    model = 'bym2',
    scale.model = T,
    constr = T,
    rankdef = 1,
    graph = Besag_prec_first_level,
    hyper = spatial_hyper)

## Specify base-formula on ADM4
### w. RW1
base_formula_1_second_level <- sampled_counts ~ 1 + f(time_id, 
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

### w. RW2
base_formula_2_second_level <- sampled_counts ~ 1 + f(time_id, 
                                                    model = 'bym2',
                                                    scale.model = T, 
                                                    constr = T, 
                                                    graph = RW2_prec,
                                                    hyper = temporal_hyper) + 
  f(area_id, 
    model = 'bym2',
    scale.model = T,
    constr = T,
    rankdef = 1,
    graph = Besag_prec_second_level,
    hyper = spatial_hyper)

#----
## Specify all the constraints and the precision matrices of the interactions

#Get precision matric for type II interaction by Kronecker product
### w. RW1
typeII_prec_1_first_level <- scaled_RW1_prec %x% diag(nrow(first_level_admin_map))

typeII_prec_1_second_level <- scaled_RW1_prec %x% diag(nrow(second_level_admin_map))

### w. RW2
typeII_prec_2_first_level <- scaled_RW2_prec %x% diag(nrow(first_level_admin_map))

typeII_prec_2_second_level <- scaled_RW2_prec %x% diag(nrow(second_level_admin_map))


#Get precision matric for type III interaction by Kronecker product
typeIII_prec_first_level <- diag(tT) %x% scaled_besag_prec_first_level

typeIII_prec_second_level <- diag(tT) %x% scaled_besag_prec_second_level

#Get type IV interaction precision matrix
### w. RW1
typeIV_prec_1_first_level <- scaled_RW1_prec %x% scaled_besag_prec_first_level
typeIV_prec_1_second_level <- scaled_RW1_prec %x% scaled_besag_prec_second_level

### w. RW2
typeIV_prec_2_first_level <- scaled_RW2_prec %x% scaled_besag_prec_first_level
typeIV_prec_2_second_level <- scaled_RW2_prec %x% scaled_besag_prec_second_level

#Get sum-to-zero constraints for type II interaction
### w. RW1
typeII_constraints_1_first_level = constraints_maker(type = "II", 
                                                     n = nrow(first_level_admin_map), 
                                                     t = tT)

typeII_constraints_1_second_level = constraints_maker(type = "II", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)
### w. RW2
typeII_constraints_2_first_level = constraints_maker(type = "II", 
                                                   n = nrow(first_level_admin_map), 
                                                   t = tT,
                                                   rw = "RW2",
                                                   prec_matrix = typeII_prec_2_first_level)

typeII_constraints_2_second_level = constraints_maker(type = "II", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT,
                                                    rw = "RW2",
                                                    prec_matrix = typeII_prec_2_second_level)


#Get sum-to-zero constraints for type II interaction
typeIII_constraints_first_level = constraints_maker(type = "III", 
                                                    n = nrow(first_level_admin_map), 
                                                    t = tT)

typeIII_constraints_second_level = constraints_maker(type = "III", 
                                                     n = nrow(second_level_admin_map), 
                                                     t = tT)

#Get sum-to-zero constraints for type IV interaction
### w. RW1
typeIV_constraints_1_first_level = constraints_maker(type = "IV", 
                                                   n = nrow(first_level_admin_map), 
                                                   t = tT)

typeIV_constraints_1_second_level = constraints_maker(type = "IV", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)

### w. RW2
typeIV_constraints_2_first_level = constraints_maker(type = "IV", 
                                                   n = nrow(first_level_admin_map), 
                                                   t = tT,
                                                   rw = "RW2",
                                                   prec_matrix = typeIV_prec_2_first_level)

typeIV_constraints_2_second_level = constraints_maker(type = "IV", 
                                                     n = nrow(second_level_admin_map), 
                                                     t = tT,
                                                     rw = "RW2",
                                                     prec_matrix = typeIV_prec_2_second_level)

#----
## Specify all the formulas w. interactions

# Type I formula
### w. RW1
typeI_1_formula_first_level <- update(base_formula_1_first_level,  
                                    ~. + f(space.time,
                                           model = "iid", 
                                           hyper = interaction_hyper ))

typeI_1_formula_second_level <- update(base_formula_1_second_level,  
                                     ~. + f(space.time,
                                            model = "iid", 
                                            hyper = interaction_hyper ))

### w. RW2
typeI_2_formula_first_level <- update(base_formula_2_first_level,  
                                    ~. + f(space.time,
                                           model = "iid", 
                                           hyper = interaction_hyper ))

typeI_2_formula_second_level <- update(base_formula_2_second_level,  
                                       ~. + f(space.time,
                                              model = "iid", 
                                              hyper = interaction_hyper ))


# Get typeII formula
### w. RW1
typeII_1_formula_first_level <- update(base_formula_1_first_level, 
                                     ~. + f(space.time, 
                                            model = "generic0", 
                                            Cmatrix = typeII_prec_1_first_level, 
                                            extraconstr = typeII_constraints_1_first_level, 
                                            rankdef = nrow(first_level_admin_map), 
                                            hyper = interaction_hyper))

typeII_1_formula_second_level <- update(base_formula_1_second_level, 
                                      ~. + f(space.time, 
                                             model = "generic0", 
                                             Cmatrix = typeII_prec_1_second_level, 
                                             extraconstr = typeII_constraints_1_second_level, 
                                             rankdef = nrow(second_level_admin_map), 
                                             hyper = interaction_hyper))

### w. RW2
typeII_2_formula_first_level <- update(base_formula_2_first_level, 
                                     ~. + f(space.time, 
                                            model = "generic0", 
                                            Cmatrix = typeII_prec_2_first_level, 
                                            extraconstr = typeII_constraints_2_first_level, 
                                            hyper = interaction_hyper))

typeII_2_formula_second_level <- update(base_formula_2_second_level, 
                                      ~. + f(space.time, 
                                             model = "generic0", 
                                             Cmatrix = typeII_prec_2_second_level, 
                                             extraconstr = typeII_constraints_2_second_level, 
                                             hyper = interaction_hyper))

# Get typeIII formula
### w. RW1
typeIII_1_formula_first_level <- update(base_formula_1_first_level, 
                                      ~. + f(space.time, 
                                             model = "generic0", 
                                             Cmatrix = typeIII_prec_first_level, 
                                             extraconstr = typeIII_constraints_first_level, 
                                             rankdef = tT, 
                                             hyper = interaction_hyper))

typeIII_1_formula_second_level <- update(base_formula_1_second_level, 
                                       ~. + f(space.time, 
                                              model = "generic0", 
                                              Cmatrix = typeIII_prec_second_level, 
                                              extraconstr = typeIII_constraints_second_level, 
                                              rankdef = tT, 
                                              hyper = interaction_hyper))

### w. RW2
typeIII_2_formula_first_level <- update(base_formula_2_first_level, 
                                        ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeIII_prec_first_level, 
                                               extraconstr = typeIII_constraints_first_level, 
                                               rankdef = tT, 
                                               hyper = interaction_hyper))

typeIII_2_formula_second_level <- update(base_formula_2_second_level, 
                                         ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeIII_prec_second_level, 
                                                extraconstr = typeIII_constraints_second_level, 
                                                rankdef = tT, 
                                                hyper = interaction_hyper))


#Get formula for type IV
### w. RW1
typeIV_1_formula_first_level <- update(base_formula_1_first_level,
                                     ~. + f(space.time, 
                                            model = "generic0",
                                            Cmatrix = typeIV_prec_1_first_level,
                                            extraconstr = typeIV_constraints_1_first_level,
                                            rankdef = (nrow(first_level_admin_map) + tT - 1), 
                                            hyper = interaction_hyper))

typeIV_1_formula_second_level <- update(base_formula_1_second_level,
                                      ~. + f(space.time, 
                                             model = "generic0",
                                             Cmatrix = typeIV_prec_1_second_level,
                                             extraconstr = typeIV_constraints_1_second_level,
                                             rankdef = (nrow(second_level_admin_map) + tT - 1), 
                                             hyper = interaction_hyper))

### w. RW2
typeIV_2_formula_first_level <- update(base_formula_2_first_level,
                                     ~. + f(space.time, 
                                            model = "generic0",
                                            Cmatrix = typeIV_prec_2_first_level,
                                            extraconstr = typeIV_constraints_2_first_level,
                                            hyper = interaction_hyper))

typeIV_2_formula_second_level <- update(base_formula_2_second_level,
                                       ~. + f(space.time, 
                                              model = "generic0",
                                              Cmatrix = typeIV_prec_2_second_level,
                                              extraconstr = typeIV_constraints_2_second_level,
                                              hyper = interaction_hyper))


#####
# Fit the Improper models

scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")
scenario_names_ADM4 = c("sc2", "sc4", "sc6", "sc8", "sc10", "sc12")

for(scenario_name in scenario_names_ADM1){
  print(paste("---", scenario_name, "---"))
  #Load in a singular data set
  load(paste('./Simulated_data/', scenario_name, '/', 
             scenario_name, '_data.RData',sep = ""))
  
  lambda <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda$sampled_counts = lambda.df$sampled_counts[, 2] #Just a chosen data set
  
  # Set the last three years to unkown
  lambda[lambda$time_id %in% 11:13, ]$sampled_counts = NA
  
  # Fit each model to the data set
  
  base_1_first_level <- inla(base_formula_1_first_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  base_2_first_level <- inla(base_formula_2_first_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    
  
  typeI_1_first_level <- inla(typeI_1_formula_first_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeI_2_first_level <- inla(typeI_2_formula_first_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeII_1_first_level <- inla(typeII_1_formula_first_level, 
                              data = lambda, 
                              family = "poisson",
                              E = E_it, #E_it
                              control.predictor = list(compute = TRUE, link = 1),       #For predictions
                              control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeII_2_first_level <- inla(typeII_2_formula_first_level, 
                              data = lambda, 
                              family = "poisson",
                              E = E_it, #E_it
                              control.predictor = list(compute = TRUE, link = 1),       #For predictions
                              control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    
  typeIII_1_first_level <- inla(typeIII_1_formula_first_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeIII_2_first_level <- inla(typeIII_2_formula_first_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal

  
  
  typeIV_1_first_level <- inla(typeIV_1_formula_first_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeIV_2_first_level <- inla(typeIV_2_formula_first_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  
  # Save the INLA-objects
  save(base_1_first_level,
       base_2_first_level,
       typeI_1_first_level,
       typeI_2_first_level,
       typeII_1_first_level,
       typeII_2_first_level,
       typeIII_1_first_level,
       typeIII_2_first_level,
       typeIV_1_first_level,
       typeIV_2_first_level,
       file = paste('diagnostics_', scenario_name, ".RData", sep = ""))
}

for(scenario_name in scenario_names_ADM4){
  print(paste("---", scenario_name, "---"))
  #Load in a singular data set
  load(paste('./Simulated_data/', scenario_name, '/', 
             scenario_name, '_data.RData',sep = ""))
  
  lambda <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda$sampled_counts = lambda.df$sampled_counts[, 2] #Just a chosen data set
  
  # Set the last three years to unkown
  lambda[lambda$time_id %in% 11:13, ]$sampled_counts = NA
  
  # Fit each model to the data set
  
  base_1_second_level <- inla(base_formula_1_second_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  base_2_second_level <- inla(base_formula_2_second_level, 
                             data = lambda, 
                             family = "poisson",
                             E = E_it, #E_it
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  
  typeI_1_second_level <- inla(typeI_1_formula_second_level, 
                              data = lambda, 
                              family = "poisson",
                              E = E_it, #E_it
                              control.predictor = list(compute = TRUE, link = 1),       #For predictions
                              control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeI_2_second_level <- inla(typeI_2_formula_second_level, 
                              data = lambda, 
                              family = "poisson",
                              E = E_it, #E_it
                              control.predictor = list(compute = TRUE, link = 1),       #For predictions
                              control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeII_1_second_level <- inla(typeII_1_formula_second_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeII_2_second_level <- inla(typeII_2_formula_second_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeIII_1_second_level <- inla(typeIII_1_formula_second_level, 
                                data = lambda, 
                                family = "poisson",
                                E = E_it, #E_it
                                control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeIII_2_second_level <- inla(typeIII_2_formula_second_level, 
                                data = lambda, 
                                family = "poisson",
                                E = E_it, #E_it
                                control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  
  
  typeIV_1_second_level <- inla(typeIV_1_formula_second_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  typeIV_2_second_level <- inla(typeIV_2_formula_second_level, 
                               data = lambda, 
                               family = "poisson",
                               E = E_it, #E_it
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
  
  
  # Save the INLA-objects
  save(base_1_second_level,
       base_2_second_level,
       typeI_1_second_level,
       typeI_2_second_level,
       typeII_1_second_level,
       typeII_2_second_level,
       typeIII_1_second_level,
       typeIII_2_second_level,
       typeIV_1_second_level,
       typeIV_2_second_level,
       file = paste('diagnostics_', scenario_name, ".RData", sep = ""))
}


################################################################################
# Specify the formulas for the Proper models

### Temporal hyperparameters (prec. of AR1 and AR1's mixing param) w. corresponding priors: penalized constraint 
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), 
                 rho = list(prior = 'pc.cor1', 
                            param = c(0.5, 0.5 + 1E-6))) #, mean = list(prior = 'normal', param = c(0, 1), fixed = TRUE)) 


### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01))) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) 


