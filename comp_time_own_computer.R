#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

n_ADM4 <- nrow(second_level_admin_map)

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

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))
#---

## Specify precision matrices
#---
### Specify the RW1 precision matrix
RW2_prec <- INLA:::inla.rw(n = tT, order = 2, 
                           scale.model = FALSE, sparse = TRUE)

### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_second_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

#---

## Specify base-formula on ADM1
base_formula_second_level <- sampled_counts ~ 1 + f(time_id, 
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

#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW2_prec,
                                   list(A = matrix(1, 1, dim(RW2_prec)[1]),
                                        e = 0))

# get scaled Besag
scaled_besag_prec_second_level <- INLA::inla.scale.model(Besag_prec_second_level, 
                                                         constr = list(A = matrix(1,1,dim(Besag_prec_second_level)[1]), 
                                                                       e = 0))
#Get type IV interaction precision matrix
typeIV_prec_second_level <- scaled_RW_prec%x% scaled_besag_prec_second_level


#Get sum-to-zero constraints for type IV interaction
typeIV_constraints_second_level = constraints_maker(type = "IV", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT,
                                                    rw = "RW2",
                                                    prec_matrix = typeIV_prec_second_level)


#Get formula for type IV
typeIV_formula_second_level <- update(base_formula_second_level,
                                      ~. + f(space.time, 
                                             model = "generic0",
                                             Cmatrix = typeIV_prec_second_level,
                                             extraconstr = typeIV_constraints_second_level,
                                             hyper = interaction_hyper))



load("Simulated_data/sc2/sc2_data.RData")

lambda_sc.df <- lambda.df[, c("area_id", "time_id", "E_it", 
                              "space.time")]

lambda_sc.df$sampled_counts = lambda.df$sampled_counts[, 1]

## Set the last three years counts to NA for the fit
lambda_sc.df[lambda_sc.df$time_id %in% 11:13, ]$sampled_counts = NA

time_ = Sys.time()
tmp_ = inla(typeIV_formula_second_level, 
            data = lambda_sc.df, 
            family = "poisson",
            E = E_it, #E_it
            verbose = T,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            #control.family = list(control.link = list(model = "log")),
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE)) # Get the lin.pred.marginal


time_improper2_typeIV = Sys.time() - time_

################################################################################

ar_hyper = list(prec = list(prior = 'pc.prec', 
                            param = c(1, 0.01)),
                pacf1 = list(prior = 'pc.cor1', 
                             param = c(0.5, 0.5 + 1E-2)),
                pacf2 = list(prior = 'pc.cor0',
                             param = c(0.5, 0.5)))


### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01))) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) 


### Group hyper
group_hyper = list(pacf1 = list(prior = 'pc.cor1', 
                                param = c(0.5, 0.5 + 1E-2)),
                   pacf2 = list(prior = 'pc.cor0',
                                param = c(0.5, 0.5)))


proper_onlyInt_formula_second_level <- sampled_counts ~ 1 + time_id + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_second_level,
    hyper = spatial_hyper,
    group = time_id, 
    control.group = list(model = "ar", 
                         order = 2,
                         hyper = group_hyper))

## Reorder due to change in space.time interaction
lambda_sc.df <- lambda_sc.df[order(lambda_sc.df$area_id, decreasing = F), ]
rownames(lambda_sc.df) <- 1:nrow(lambda_sc.df)

## Add copies of area and time ids, INLA requires unique random effects
lambda_sc.df$area_id.copy <- lambda_sc.df$area_id
lambda_sc.df$time_id.copy <- lambda_sc.df$time_id

time_ = Sys.time()
tmp_ = inla(proper_onlyInt_formula_second_level, 
            data = lambda_sc.df, 
            family = "poisson",
            E = E_it, #E_it
            verbose = T,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            #control.family = list(control.link = list(model = "log")),
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
time_prop2_onlyInt = Sys.time() -  time_





