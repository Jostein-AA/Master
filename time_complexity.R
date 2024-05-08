#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

# Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)
n_ADMnew <- nrow(new_map)

dataset_ids = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")
scenario_names_ADM4 = c("sc2", "sc4", "sc6", "sc8", "sc10", "sc12")
scenario_names = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11",
                   "sc2", "sc4", "sc6", "sc8", "sc10", "sc12")

# Set upper-time limit to 600 sec (10 minutes)
inla.setOption(inla.timeout = 600) 

################################################################################
# Specify the hyperpriors

### Temporal hyperparameters (Precision of iid and precision of RW1 and RW2) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

### Temporal hyperparameters (prec. of AR1 and AR1's mixing param) w. corresponding priors: penalized constraint 
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), 
                 rho = list(prior = 'pc.cor1', 
                            param = c(0.5, 0.5 + 1E-2))) #, mean = list(prior = 'normal', param = c(0, 1), fixed = TRUE)) 

### Temporal hyperparameters (prec. of AR2 and AR2's autocorrelation param) w. corresponding priors: penalized constraint 
ar2_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)),
                 pacf1 = list(prior = 'pc.cor1', 
                              param = c(0.5, 0.5 + 1E-2)),
                 pacf2 = list(prior = 'pc.cor0',
                              param = c(0.5, 0.5)))

### Group hyper
group_hyper_ar1 = list(rho = list(prior = 'pc.cor1', 
                                  param = c(0.5, 0.5 + 1E-2)))

### group hyper for proper models
group_hyper_ar2 = list(pacf1 = list(prior = 'pc.cor1', 
                                    param = c(0.5, 0.5 + 1E-2)),
                       pacf2 = list(prior = 'pc.cor0',
                                    param = c(0.5, 0.5)))

### Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper_proper = list(prec= list(prior = 'pc.prec', 
                                       param = c(1, 0.01))) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) 

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))

################################################################################



# Specify the precision matrices

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


# For interactions
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


################################################################################

# Specify all the constraints and the precision matrices of the interactions


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

################################################################################
# Specify all the formulas

## Specify base-formula on ADM1
### w. RW1
Improper1_noInt_ADM1_formula <- sampled_counts ~ 1 + f(time_id, 
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
Improper2_noInt_ADM1_formula <- sampled_counts ~ 1 + f(time_id, 
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
Improper1_noInt_ADM4_formula <- sampled_counts ~ 1 + f(time_id, 
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
Improper2_noInt_ADM4_formula <- sampled_counts ~ 1 + f(time_id, 
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


## Specify all the formulas w. interactions

# Type I formula
### w. RW1
Improper1_typeI_ADM1_formula <- update(Improper1_noInt_ADM1_formula,  
                                       ~. + f(space.time,
                                              model = "iid", 
                                              hyper = interaction_hyper ))

Improper1_typeI_ADM4_formula <- update(Improper1_noInt_ADM4_formula,  
                                       ~. + f(space.time,
                                              model = "iid", 
                                              hyper = interaction_hyper ))

### w. RW2
Improper2_typeI_ADM1_formula <- update(Improper2_noInt_ADM1_formula,  
                                       ~. + f(space.time,
                                              model = "iid", 
                                              hyper = interaction_hyper ))

Improper2_typeI_ADM4_formula <- update(Improper2_noInt_ADM4_formula,  
                                       ~. + f(space.time,
                                              model = "iid", 
                                              hyper = interaction_hyper ))


# Get typeII formula
### w. RW1
Improper1_typeII_ADM1_formula <- update(Improper1_noInt_ADM1_formula, 
                                        ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeII_prec_1_first_level, 
                                               extraconstr = typeII_constraints_1_first_level, 
                                               rankdef = nrow(first_level_admin_map), 
                                               hyper = interaction_hyper))

Improper1_typeII_ADM4_formula <- update(Improper1_noInt_ADM4_formula, 
                                        ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeII_prec_1_second_level, 
                                               extraconstr = typeII_constraints_1_second_level, 
                                               rankdef = nrow(second_level_admin_map), 
                                               hyper = interaction_hyper))

### w. RW2
Improper2_typeII_ADM1_formula <- update(Improper2_noInt_ADM1_formula, 
                                        ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeII_prec_2_first_level, 
                                               extraconstr = typeII_constraints_2_first_level, 
                                               hyper = interaction_hyper))

Improper2_typeII_ADM4_formula <- update(Improper2_noInt_ADM4_formula, 
                                        ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeII_prec_2_second_level, 
                                               extraconstr = typeII_constraints_2_second_level, 
                                               hyper = interaction_hyper))

# Get typeIII formula
### w. RW1
Improper1_typeIII_ADM1_formula <- update(Improper1_noInt_ADM1_formula, 
                                         ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeIII_prec_first_level, 
                                                extraconstr = typeIII_constraints_first_level, 
                                                rankdef = tT, 
                                                hyper = interaction_hyper))

Improper1_typeIII_ADM4_formula <- update(Improper1_noInt_ADM4_formula, 
                                         ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeIII_prec_second_level, 
                                                extraconstr = typeIII_constraints_second_level, 
                                                rankdef = tT, 
                                                hyper = interaction_hyper))

### w. RW2
Improper2_typeIII_ADM1_formula <- update(Improper2_noInt_ADM1_formula, 
                                         ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeIII_prec_first_level, 
                                                extraconstr = typeIII_constraints_first_level, 
                                                rankdef = tT, 
                                                hyper = interaction_hyper))

Improper2_typeIII_ADM4_formula <- update(Improper2_noInt_ADM4_formula, 
                                         ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeIII_prec_second_level, 
                                                extraconstr = typeIII_constraints_second_level, 
                                                rankdef = tT, 
                                                hyper = interaction_hyper))


#Get formula for type IV
### w. RW1
Improper1_typeIV_ADM1_formula <- update(Improper1_noInt_ADM1_formula,
                                        ~. + f(space.time, 
                                               model = "generic0",
                                               Cmatrix = typeIV_prec_1_first_level,
                                               extraconstr = typeIV_constraints_1_first_level,
                                               rankdef = (nrow(first_level_admin_map) + tT - 1), 
                                               hyper = interaction_hyper))

Improper1_typeIV_ADM4_formula <- update(Improper1_noInt_ADM4_formula,
                                        ~. + f(space.time, 
                                               model = "generic0",
                                               Cmatrix = typeIV_prec_1_second_level,
                                               extraconstr = typeIV_constraints_1_second_level,
                                               rankdef = (nrow(second_level_admin_map) + tT - 1), 
                                               hyper = interaction_hyper))

### w. RW2
Improper2_typeIV_ADM1_formula <- update(Improper2_noInt_ADM1_formula,
                                        ~. + f(space.time, 
                                               model = "generic0",
                                               Cmatrix = typeIV_prec_2_first_level,
                                               extraconstr = typeIV_constraints_2_first_level,
                                               hyper = interaction_hyper))

Improper2_typeIV_ADM4_formula <- update(Improper2_noInt_ADM4_formula,
                                        ~. + f(space.time, 
                                               model = "generic0",
                                               Cmatrix = typeIV_prec_2_second_level,
                                               extraconstr = typeIV_constraints_2_second_level,
                                               hyper = interaction_hyper))



proper1_noInt_ADM1_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper)

proper1_noInt_ADM4_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper)

proper1_onlyInt_ADM1_formula <- sampled_counts ~ 1 + time_id + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar1",
                         hyper = group_hyper_ar1))

proper1_onlyInt_ADM4_formula <- sampled_counts ~ 1 + time_id + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar1",
                         hyper = group_hyper_ar1))

proper1_full_ADM1_formula <- sampled_counts ~ 1 + time_id + 
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) +
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper) + 
  f(area_id.copy, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar1",
                         hyper = group_hyper_ar1)) 

proper1_full_ADM4_formula <- sampled_counts ~ 1 + time_id + 
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) +
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper) + 
  f(area_id.copy, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar1",
                         hyper = group_hyper_ar1))

proper1_iid_ADM1_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper) + 
  f(space.time,
    model = "iid", 
    hyper = interaction_hyper )

proper1_iid_ADM4_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper) + 
  f(space.time,
    model = "iid", 
    hyper = interaction_hyper )


proper2_noInt_ADM1_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar",
    order = 2,
    hyper = ar2_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper)

proper2_noInt_ADM4_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar",
    order = 2,
    hyper = ar2_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper)

proper2_onlyInt_ADM1_formula <- sampled_counts ~ 1 + time_id + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar", 
                         order = 2,
                         hyper = group_hyper_ar2))

proper2_onlyInt_ADM4_formula <- sampled_counts ~ 1 + time_id + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar", 
                         order = 2,
                         hyper = group_hyper_ar2))

proper2_full_ADM1_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar",
    order = 2,
    hyper = ar2_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper) + 
  f(area_id.copy, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar", 
                         order = 2,
                         hyper = group_hyper_ar2))


proper2_full_ADM4_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar",
    order = 2,
    hyper = ar2_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper) + 
  f(area_id.copy, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper,
    group = time_id, 
    control.group = list(model = "ar", 
                         order = 2,
                         hyper = group_hyper_ar2))

proper2_iid_ADM1_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar",
    order = 2,
    hyper = ar2_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_first_level,
    hyper = spatial_hyper_proper) + 
  f(space.time,
    model = "iid", 
    hyper = interaction_hyper )

proper2_iid_ADM4_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar",
    order = 2,
    hyper = ar2_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper_proper) + 
  f(space.time,
    model = "iid", 
    hyper = interaction_hyper )




################################################################################
# Fit the models on the first-level map

scenario_names_ADM1 = c("sc1", "sc3", "sc5", "sc7", "sc9", "sc11")


for(scenario_name in scenario_names_ADM1){
  for(dataset_id in dataset_ids){
  print(paste("---", scenario_name, "---"))
  #Load in a singular data set
  load(paste('./Simulated_data/', scenario_name, '/', 
             scenario_name, '_data.RData',sep = ""))
  
  lambda <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda$sampled_counts = lambda.df$sampled_counts[, dataset_id] #Just a chosen data set
  
  # Set the last three years to unkown
  lambda[lambda$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  
  
  # Fit each model to the data set
  time_ = Sys.time()
  tryCatch(
    {
      Improper1_noInt_ADM1 <- inla(Improper1_noInt_ADM1_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      
      
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper1_noInt_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper2_noInt_ADM1 <- inla(Improper2_noInt_ADM1_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper2_noInt_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()  
  tryCatch(
    {
      Improper1_typeI_ADM1 <- inla(Improper1_typeI_ADM1_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper1_typeI_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper2_typeI_ADM1 <- inla(Improper2_typeI_ADM1_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper2_typeI_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper1_typeII_ADM1 <- inla(Improper1_typeII_ADM1_formula, 
                                    data = lambda, 
                                    family = "poisson",
                                    E = E_it, #E_it
                                    control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                    control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper1_typeII_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper2_typeII_ADM1 <- inla(Improper2_typeII_ADM1_formula, 
                                    data = lambda, 
                                    family = "poisson",
                                    E = E_it, #E_it
                                    control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                    control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper2_typeII_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper1_typeIII_ADM1 <- inla(Improper1_typeIII_ADM1_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper1_typeIII_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper2_typeIII_ADM1 <- inla(Improper2_typeIII_ADM1_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper2_typeIII_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper1_typeIV_ADM1 <- inla(Improper1_typeIV_ADM1_formula, 
                                    data = lambda, 
                                    family = "poisson",
                                    E = E_it, #E_it
                                    control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                    control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper1_typeIV_ADM1 <- Sys.time() - time_ 
  
  time_ = Sys.time()
  tryCatch(
    {
      Improper2_typeIV_ADM1 <- inla(Improper2_typeIV_ADM1_formula, 
                                    data = lambda, 
                                    family = "poisson",
                                    E = E_it, #E_it
                                    control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                    control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_Improper2_typeIV_ADM1 <- Sys.time() - time_ 
  
  print(paste("Improper models done for", scenario_name))
  
  # Resort the data for the proper models
  ## Reorder due to change in space.time interaction
  lambda <- lambda[order(lambda$area_id, decreasing = F), ]
  rownames(lambda) <- 1:nrow(lambda)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda$area_id.copy <- lambda$area_id
  lambda$time_id.copy <- lambda$time_id
  
  print("reordering of lambda done")
  
  time_ = Sys.time()
  tryCatch(
    {
      proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                                 data = lambda, 
                                 family = "poisson",
                                 E = E_it, #E_it
                                 control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                 control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_proper1_noInt_ADM1 <- Sys.time() - time_
  
  time_ = Sys.time()
  tryCatch(
    {
      proper2_noInt_ADM1 <- inla(proper2_noInt_ADM1_formula, 
                                 data = lambda, 
                                 family = "poisson",
                                 E = E_it, #E_it
                                 control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                 control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_proper2_noInt_ADM1 <- Sys.time() - time_
  
  time_ = Sys.time()
  tryCatch(
    {
      proper1_onlyInt_ADM1 <- inla(proper1_onlyInt_ADM1_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_proper1_onlyInt_ADM1 <- Sys.time() - time_
  
  time_ = Sys.time()
  tryCatch(
    {
      proper2_onlyInt_ADM1 <- inla(proper2_onlyInt_ADM1_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_proper2_onlyInt_ADM1 <- Sys.time() - time_
  
  time_ = Sys.time()
  tryCatch(
    {
      proper1_full_ADM1 <- inla(proper1_full_ADM1_formula, 
                                data = lambda, 
                                family = "poisson",
                                E = E_it, #E_it
                                control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_proper1_full_ADM1 <- Sys.time() - time_
  
  time_ = Sys.time()
  tryCatch(
    {
      proper2_full_ADM1 <- inla(proper2_full_ADM1_formula, 
                                data = lambda, 
                                family = "poisson",
                                E = E_it, #E_it
                                control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
    },
    error = function(cond) {
      print("!Error!")
    },
    warning = function(cond) {
      print("Warning")
    },
    finally = {
      print("H")
    }
  )
  time_proper2_full_ADM1 <- Sys.time() - time_ 
  
  # Get timecomplexity into dataframe to save
  time_complexity_ADM1 <- data.frame(model_name = c("Improper1_noInt",
                                                    "Improper2_noInt", 
                                                    "Improper1_typeI",
                                                    "Improper2_typeI",
                                                    "Improper1_typeII",
                                                    "Improper2_typeII",
                                                    "Improper1_typeIII",
                                                    "Improper2_typeIII",
                                                    "Improper1_typeIV",
                                                    "Improper2_typeIV",
                                                    "proper1_noInt",
                                                    "proper2_noInt",
                                                    "proper1_onlyInt",
                                                    "proper2_onlyInt",
                                                    "proper1_full",
                                                    "proper2_full"),
                                     time = c(time_Improper1_noInt_ADM1,
                                              time_Improper2_noInt_ADM1,
                                              time_Improper1_typeI_ADM1,
                                              time_Improper2_typeI_ADM1,
                                              time_Improper1_typeII_ADM1,
                                              time_Improper2_typeII_ADM1,
                                              time_Improper1_typeIII_ADM1,
                                              time_Improper2_typeIII_ADM1,
                                              time_Improper1_typeIV_ADM1,
                                              time_Improper2_typeIV_ADM1,
                                              time_proper1_noInt_ADM1,
                                              time_proper2_noInt_ADM1,
                                              time_proper1_onlyInt_ADM1,
                                              time_proper2_onlyInt_ADM1,
                                              time_proper1_full_ADM1,
                                              time_proper2_full_ADM1))
  
  
  # Save the INLA-objects
  save(time_complexity_ADM1,
       file = paste("results/time_complexity/", scenario_name,"_", toString(dataset_id), 
                    ".RData", sep = ""))
  
  }
}

#####
# Fit the models on the second-level map
scenario_names_ADM4 = c("sc2", "sc4", "sc6", "sc8", "sc10", "sc12")


for(scenario_name in scenario_names_ADM4){
  for(dataset_id in dataset_ids){
    print(paste("---", scenario_name, "---"))
    #Load in a singular data set
    load(paste('./Simulated_data/', scenario_name, '/', 
               scenario_name, '_data.RData',sep = ""))
    
    lambda <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
    lambda$sampled_counts = lambda.df$sampled_counts[, dataset_id] #Just a chosen data set
    
    # Set the last three years to unkown
    lambda[lambda$time_id %in% 11:13, ]$sampled_counts = NA
    
    
    
    
    # Fit each model to the data set
    time_ = Sys.time()
    tryCatch(
      {
        Improper1_noInt_ADM4 <- inla(Improper1_noInt_ADM4_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
        
        
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper1_noInt_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper2_noInt_ADM4 <- inla(Improper2_noInt_ADM4_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper2_noInt_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()  
    tryCatch(
      {
        Improper1_typeI_ADM4 <- inla(Improper1_typeI_ADM4_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper1_typeI_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper2_typeI_ADM4 <- inla(Improper2_typeI_ADM4_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper2_typeI_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper1_typeII_ADM4 <- inla(Improper1_typeII_ADM4_formula, 
                                      data = lambda, 
                                      family = "poisson",
                                      E = E_it, #E_it
                                      control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                      control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper1_typeII_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper2_typeII_ADM4 <- inla(Improper2_typeII_ADM4_formula, 
                                      data = lambda, 
                                      family = "poisson",
                                      E = E_it, #E_it
                                      control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                      control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper2_typeII_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper1_typeIII_ADM4 <- inla(Improper1_typeIII_ADM4_formula, 
                                       data = lambda, 
                                       family = "poisson",
                                       E = E_it, #E_it
                                       control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                       control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper1_typeIII_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper2_typeIII_ADM4 <- inla(Improper2_typeIII_ADM4_formula, 
                                       data = lambda, 
                                       family = "poisson",
                                       E = E_it, #E_it
                                       control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                       control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper2_typeIII_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper1_typeIV_ADM4 <- inla(Improper1_typeIV_ADM4_formula, 
                                      data = lambda, 
                                      family = "poisson",
                                      E = E_it, #E_it
                                      control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                      control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper1_typeIV_ADM4 <- Sys.time() - time_ 
    
    time_ = Sys.time()
    tryCatch(
      {
        Improper2_typeIV_ADM4 <- inla(Improper2_typeIV_ADM4_formula, 
                                      data = lambda, 
                                      family = "poisson",
                                      E = E_it, #E_it
                                      control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                      control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_Improper2_typeIV_ADM4 <- Sys.time() - time_ 
    
    print(paste("Improper models done for", scenario_name))
    
    # Resort the data for the proper models
    ## Reorder due to change in space.time interaction
    lambda <- lambda[order(lambda$area_id, decreasing = F), ]
    rownames(lambda) <- 1:nrow(lambda)
    
    ## Add copies of area and time ids, INLA requires unique random effects
    lambda$area_id.copy <- lambda$area_id
    lambda$time_id.copy <- lambda$time_id
    
    print("reordering of lambda done")
    
    time_ = Sys.time()
    tryCatch(
      {
        proper1_noInt_ADM4 <- inla(proper1_noInt_ADM4_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_proper1_noInt_ADM4 <- Sys.time() - time_
    
    time_ = Sys.time()
    tryCatch(
      {
        proper2_noInt_ADM4 <- inla(proper2_noInt_ADM4_formula, 
                                   data = lambda, 
                                   family = "poisson",
                                   E = E_it, #E_it
                                   control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                   control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_proper2_noInt_ADM4 <- Sys.time() - time_
    
    time_ = Sys.time()
    tryCatch(
      {
        proper1_onlyInt_ADM4 <- inla(proper1_onlyInt_ADM4_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_proper1_onlyInt_ADM4 <- Sys.time() - time_
    
    time_ = Sys.time()
    tryCatch(
      {
        proper2_onlyInt_ADM4 <- inla(proper2_onlyInt_ADM4_formula, 
                                     data = lambda, 
                                     family = "poisson",
                                     E = E_it, #E_it
                                     control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                     control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_proper2_onlyInt_ADM4 <- Sys.time() - time_
    
    time_ = Sys.time()
    tryCatch(
      {
        proper1_full_ADM4 <- inla(proper1_full_ADM4_formula, 
                                  data = lambda, 
                                  family = "poisson",
                                  E = E_it, #E_it
                                  control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                  control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_proper1_full_ADM4 <- Sys.time() - time_
    
    time_ = Sys.time()
    tryCatch(
      {
        proper2_full_ADM4 <- inla(proper2_full_ADM4_formula, 
                                  data = lambda, 
                                  family = "poisson",
                                  E = E_it, #E_it
                                  control.predictor = list(compute = TRUE, link = 1),       #For predictions
                                  control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      },
      error = function(cond) {
        print("!Error!")
      },
      warning = function(cond) {
        print("Warning")
      },
      finally = {
        print("H")
      }
    )
    time_proper2_full_ADM4 <- Sys.time() - time_ 
    
    # Get timecomplexity into dataframe to save
    time_complexity_ADM4 <- data.frame(model_name = c("Improper1_noInt",
                                                      "Improper2_noInt", 
                                                      "Improper1_typeI",
                                                      "Improper2_typeI",
                                                      "Improper1_typeII",
                                                      "Improper2_typeII",
                                                      "Improper1_typeIII",
                                                      "Improper2_typeIII",
                                                      "Improper1_typeIV",
                                                      "Improper2_typeIV",
                                                      "proper1_noInt",
                                                      "proper2_noInt",
                                                      "proper1_onlyInt",
                                                      "proper2_onlyInt",
                                                      "proper1_full",
                                                      "proper2_full"),
                                       time = c(time_Improper1_noInt_ADM4,
                                                time_Improper2_noInt_ADM4,
                                                time_Improper1_typeI_ADM4,
                                                time_Improper2_typeI_ADM4,
                                                time_Improper1_typeII_ADM4,
                                                time_Improper2_typeII_ADM4,
                                                time_Improper1_typeIII_ADM4,
                                                time_Improper2_typeIII_ADM4,
                                                time_Improper1_typeIV_ADM4,
                                                time_Improper2_typeIV_ADM4,
                                                time_proper1_noInt_ADM4,
                                                time_proper2_noInt_ADM4,
                                                time_proper1_onlyInt_ADM4,
                                                time_proper2_onlyInt_ADM4,
                                                time_proper1_full_ADM4,
                                                time_proper2_full_ADM4))
    
    
    # Save the INLA-objects
    save(time_complexity_ADM4,
         file = paste("results/time_complexity/", scenario_name,"_", toString(dataset_id), 
                      ".RData", sep = ""))
    
  }
}



################################################################################
# Create a df containing all the times for each model
dir = "results/time_complexity/"
scenario_name = "sc2"

complete_time_complexity_df <- data.frame(model_name = rep(c("Improper1_noInt",
                                                         "Improper2_noInt", 
                                                         "Improper1_typeI",
                                                         "Improper2_typeI",
                                                         "Improper1_typeII",
                                                         "Improper2_typeII",
                                                         "Improper1_typeIII",
                                                         "Improper2_typeIII",
                                                         "Improper1_typeIV",
                                                         "Improper2_typeIV",
                                                         "proper1_noInt",
                                                         "proper2_noInt",
                                                         "proper1_onlyInt",
                                                         "proper2_onlyInt",
                                                         "proper1_full",
                                                         "proper2_full"), 12),
                                          scenario_name = c(rep("sc1", 16), rep("sc3", 16), rep("sc5", 16),
                                                            rep("sc7", 16), rep("sc9", 16), rep("sc11", 16),
                                                            rep("sc2", 16), rep("sc4", 16), rep("sc6", 16),
                                                            rep("sc8", 16), rep("sc10", 16), rep("sc12", 16)),
                                          time_1 = rep(0, 16 * 12),
                                          time_10 = rep(0, 16 * 12),
                                          time_20 = rep(0, 16 * 12),
                                          time_30 = rep(0, 16 * 12),
                                          time_40 = rep(0, 16 * 12),
                                          time_50 = rep(0, 16 * 12),
                                          time_60 = rep(0, 16 * 12),
                                          time_70 = rep(0, 16 * 12),
                                          time_80 = rep(0, 16 * 12),
                                          time_90 = rep(0, 16 * 12),
                                          time_100 = rep(0, 16 * 12),
                                          avg_time = rep(0, 16 * 12))





for(scenario_name in scenario_names){
  for(dataset_id_index in 1:length(dataset_ids)){
    dataset_id = dataset_ids[dataset_id_index]
    load(paste(dir, scenario_name, "_", dataset_id, ".RData", sep = ""))
    
    if(scenario_name %in% scenario_names_ADM1){
      for(i in 1:nrow(time_complexity_ADM1)){
        if(time_complexity_ADM1[i, ]$time >= 600){
          time_complexity_ADM1[i, ]$time = NA
        }
        complete_time_complexity_df[complete_time_complexity_df$scenario_name == scenario_name & 
          complete_time_complexity_df$model_name == time_complexity_ADM1[i, ]$model_name, 
          2 + dataset_id_index] = time_complexity_ADM1[i, ]$time
      }
    } else {
      for(i in 1:nrow(time_complexity_ADM4)){
        if(time_complexity_ADM4[i, ]$time >= 600){
          time_complexity_ADM4[i, ]$time = NA
        }
        complete_time_complexity_df[complete_time_complexity_df$scenario_name == scenario_name & 
                                      complete_time_complexity_df$model_name == time_complexity_ADM4[i, ]$model_name, 
                                    2 + dataset_id_index] = time_complexity_ADM4[i, ]$time
      }
      
      
    }
  }
}

complete_time_complexity_df[, ncol(complete_time_complexity_df)] <- rowMeans(
  complete_time_complexity_df[,(3):(ncol(complete_time_complexity_df) - 1)],
         na.rm=TRUE)


# Create plot of time-complexity

### First, format for ggplot
to_plot_time_complex <- data.frame(model_name = rep(c("Improper1_noInt",
                                                      "Improper2_noInt", 
                                                      "Improper1_typeI",
                                                      "Improper2_typeI",
                                                      "Improper1_typeII",
                                                      "Improper2_typeII",
                                                      "Improper1_typeIII",
                                                      "Improper2_typeIII",
                                                      "Improper1_typeIV",
                                                      "Improper2_typeIV",
                                                      "proper1_noInt",
                                                      "proper2_noInt",
                                                      "proper1_onlyInt",
                                                      "proper2_onlyInt",
                                                      "proper1_full",
                                                      "proper2_full"), 12 * 10),
                                   prop_improp = rep(c(rep("Improper models", 10), rep("Proper models", 6)),
                                                     12 * 10),
                                   ADM = c(rep("ADM1", 16 * 10 * 6), rep("ADM4", 16 * 10 * 6)),
                                   scenario_name = c(rep("sc1", 16 * 10), rep("sc3", 16 * 10), rep("sc5", 16 * 10),
                                                     rep("sc7", 16 * 10), rep("sc9", 16 * 10), rep("sc11", 16 * 10),
                                                     rep("sc2", 16 * 10), rep("sc4", 16 * 10), rep("sc6", 16 * 10),
                                                     rep("sc8", 16 * 10), rep("sc10", 16 * 10), rep("sc12", 16 * 10)),
                                   time = rep(NA, 12 * 16 * 10))

for(model_name in unique(to_plot_time_complex$model_name)){
  for(scenario_name in scenario_names){
    times = complete_time_complexity_df[complete_time_complexity_df$model_name == model_name & 
                                  complete_time_complexity_df$scenario_name == scenario_name, 
                                  3:(ncol(complete_time_complexity_df) - 1)]
    
    to_plot_time_complex[to_plot_time_complex$model_name == model_name & 
                           to_plot_time_complex$scenario_name == scenario_name, 5] = times
    
  }
}


### ADM 1
comp_time_ADM1_plt <- ggplot(data = to_plot_time_complex[to_plot_time_complex$ADM == "ADM1", ],
       aes(x = time, y = model_name, fill = prop_improp)) + #, 
  geom_boxplot(outlier.shape = NA) + #fill = "#00BFC4",outlier.shape = NA
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    legend.text = element_text(size = 11)
  ) + 
  labs(fill = NULL) +
  ggtitle("Computational time on ADM1 map:\n n = 16, T = 13") +
  xlab("Computation time (s)") + 
  ylab("Model")

to_plot_time_complex[to_plot_time_complex$ADM == "ADM4", ]$time = log(to_plot_time_complex[to_plot_time_complex$ADM == "ADM4", ]$time)

comp_time_ADM4_plt <-ggplot(data = to_plot_time_complex[to_plot_time_complex$ADM == "ADM4", ],
       aes(x = time, y = model_name, fill = prop_improp)) + #, 
  geom_boxplot(outlier.shape = NA) + #fill = "#00BFC4", outlier.shape = NA
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    legend.text = element_text(size = 11)
  ) + 
  labs(fill = NULL) +
  ggtitle("Computational time on ADM4 map:\n n = 402, T = 13") +
  xlab("Computation time (log(s))") + 
  ylab("Model")

# Save as comp_time_plot 8 by 4
ggarrange(comp_time_ADM1_plt, comp_time_ADM4_plt,
          ncol = 2, nrow = 1,
          common.legend = T,
          legend = "top")






























