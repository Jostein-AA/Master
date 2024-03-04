#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

## Necessary to define temporal domain's max
tT = 13
################################################################################
#Load in simulation data (done stupidly, change!)
load("./Simulated_data/sc1_data_2.RData")
lambda_sc1.df <- lambda.df
lambda_sc1.df$space.time = 1:nrow(lambda_sc1.df)

load("./Simulated_data/sc2_data_2.RData")
lambda_sc2.df <- lambda.df
lambda_sc2.df$space.time = 1:nrow(lambda_sc2.df)

load("./Simulated_data/sc3_data_2.RData")
lambda_sc3.df <- lambda.df
lambda_sc3.df$space.time = 1:nrow(lambda_sc3.df)

load("./Simulated_data/sc4_data_2.RData")
lambda_sc4.df <- lambda.df
lambda_sc4.df$space.time = 1:nrow(lambda_sc4.df)

load("./Simulated_data/sc5_data_2.RData")
lambda_sc5.df <- lambda.df
lambda_sc5.df$space.time = 1:nrow(lambda_sc5.df)

load("./Simulated_data/sc6_data_2.RData")
lambda_sc6.df <- lambda.df
lambda_sc6.df$space.time = 1:nrow(lambda_sc6.df)

load("./Simulated_data/sc7_data_2.RData")
lambda_sc7.df <- lambda.df
lambda_sc7.df$space.time = 1:nrow(lambda_sc7.df)

load("./Simulated_data/sc8_data_2.RData")
lambda_sc8.df <- lambda.df
lambda_sc8.df$space.time = 1:nrow(lambda_sc8.df)

load("./Simulated_data/sc9_data_2.RData")
lambda_sc9.df <- lambda.df
lambda_sc9.df$space.time = 1:nrow(lambda_sc9.df)

load("./Simulated_data/sc10_data_2.RData")
lambda_sc10.df <- lambda.df
lambda_sc10.df$space.time = 1:nrow(lambda_sc10.df)

load("./Simulated_data/sc11_data_2.RData")
lambda_sc11.df <- lambda.df
lambda_sc11.df$space.time = 1:nrow(lambda_sc11.df)

load("./Simulated_data/sc12_data_2.RData")
lambda_sc12.df <- lambda.df
lambda_sc12.df$space.time = 1:nrow(lambda_sc12.df)

################################################################################
# Make INLA formulas

## Specify priors for hyperparameters of improper models
### Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

### Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))



## Specify priors and precision matrices
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)


## Make precision matrix for Besag on germany_map2
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


## Specify base-formula on germany_map
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


## Specify base-formula on germany_map
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


#Get sum-to-zero constraints for type II interaction
typeII_constraints_first_level = constraints_maker(type = "II", n = nrow(first_level_admin_map), t = tT)
typeII_constraints_second_level = constraints_maker(type = "II", n = nrow(second_level_admin_map), t = tT)

#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))

#Get precision matric for type II interaction by Kronecker product
typeII_prec_first_level <- scaled_RW_prec %x% diag(nrow(first_level_admin_map))
typeII_prec_second_level <- scaled_RW_prec %x% diag(nrow(second_level_admin_map))

typeII_formula_first_level <- update(base_formula_first_level, ~. + f(space.time, 
                                              model = "generic0", 
                                              Cmatrix = typeII_prec_first_level, 
                                              extraconstr = typeII_constraints_first_level, 
                                              rankdef = nrow(first_level_admin_map), 
                                              hyper = interaction_hyper))

typeII_formula_second_level <- update(base_formula_second_level, ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeII_prec_second_level, 
                                               extraconstr = typeII_constraints_second_level, 
                                               rankdef = nrow(second_level_admin_map), 
                                               hyper = interaction_hyper))


################################################################################
# Do analysis in parallel

ptm <- Sys.time(); ptm_ = Sys.time()
improper_typeII_sc1 <- inla(typeII_formula_first_level, 
                           data = lambda_sc1.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection


print(c("improper_typeII_sc1 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()

improper_typeII_sc3 <- inla(typeII_formula_first_level, 
                           data = lambda_sc3.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection

print(paste("improper_typeII_sc3 fitted in: ", toString(Sys.time() - ptm_), sep = ""))
ptm_ = Sys.time()

improper_typeII_sc5 <- inla(typeII_formula_first_level, 
                           data = lambda_sc5.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection


print(paste("improper_typeII_sc5 fitted in: ", toString(Sys.time() - ptm_), sep = ""))
ptm_ = Sys.time()

improper_typeII_sc7 <- inla(typeII_formula_first_level, 
                           data = lambda_sc7.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection

print(paste("improper_typeII_sc7 fitted in: ", toString(Sys.time() - ptm_), sep = ""))
ptm_ = Sys.time()

improper_typeII_sc9 <- inla(typeII_formula_first_level, 
                           data = lambda_sc9.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection
print(paste("improper_typeII_sc9 fitted in: ", toString(Sys.time() - ptm_), sep = ""))
ptm_ = Sys.time()

improper_typeII_sc11 <- inla(typeII_formula_first_level, 
                            data = lambda_sc11.df, 
                            family = "poisson",
                            E = E_it, 
                            control.compute = list(config = TRUE, # To see constraints later
                                                   cpo = T)) # For model selection

print(paste("improper_typeII_sc11 fitted in: ", toString(Sys.time() - ptm_), sep = ""))
ptm_ = Sys.time()

time_improper_typeII_first_level = Sys.time() - ptm
print("----------------")
print(time_improper_typeII_first_level)
print("----------------")


ptm <- Sys.time(); ptm_ = Sys.time()
improper_typeII_sc2 <- inla(typeII_formula_second_level, 
                           data = lambda_sc2.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection


print(c("improper_typeII_sc2 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()

improper_typeII_sc4 <- inla(typeII_formula_second_level, 
                           data = lambda_sc4.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection


print(c("improper_typeII_sc4 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()

improper_typeII_sc6 <- inla(typeII_formula_second_level, 
                           data = lambda_sc6.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection


print(c("improper_typeII_sc6 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()

improper_typeII_sc8 <- inla(typeII_formula_second_level, 
                           data = lambda_sc8.df, 
                           family = "poisson",
                           E = E_it, 
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T)) # For model selection


print(c("improper_typeII_sc8 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()

improper_typeII_sc10 <- inla(typeII_formula_second_level, 
                            data = lambda_sc10.df, 
                            family = "poisson",
                            E = E_it, 
                            control.compute = list(config = TRUE, # To see constraints later
                                                   cpo = T)) # For model selection


print(c("improper_typeII_sc10 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()

improper_typeII_sc12 <- inla(typeII_formula_second_level, 
                            data = lambda_sc12.df, 
                            family = "poisson",
                            E = E_it, 
                            control.compute = list(config = TRUE, # To see constraints later
                                                   cpo = T)) # For model selection


print(c("improper_typeII_sc12 fitted in: ", Sys.time() - ptm_))
ptm_ = Sys.time()


time_improper_typeII_second_level = Sys.time() - ptm
print("----------------")
print(time_improper_typeII_second_level)
print("----------------")


save(improper_typeII_sc1,
     improper_typeII_sc2,
     improper_typeII_sc3,
     improper_typeII_sc4,
     improper_typeII_sc5,
     improper_typeII_sc6,
     improper_typeII_sc7,
     improper_typeII_sc8,
     improper_typeII_sc9,
     improper_typeII_sc10,
     improper_typeII_sc11,
     improper_typeII_sc12,
     file = "improper1_typeII_fitted.RData")