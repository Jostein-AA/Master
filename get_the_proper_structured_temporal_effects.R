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

dataset_id = 3
dataset_id_2 = 4
dataset_id_new = 2

# Set a maximum run-time of 750 sec
inla.setOption(inla.timeout = 750)

print("gulrot")

################################################################################

### Temporal hyperparameters (prec. of AR1 and AR1's mixing param) w. corresponding priors: penalized constraint 
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), 
                 rho = list(prior = 'pc.cor1', 
                            param = c(0.5, 0.5 + 1E-2))) #, mean = list(prior = 'normal', param = c(0, 1), fixed = TRUE)) 

### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper_proper = list(prec= list(prior = 'pc.prec', 
                                       param = c(1, 0.01))) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) 

################################################################################


### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_first_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_first_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

### Make precision matrix for Besag on ADM4
matrix4inla <- nb2mat(nb_second_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

### Make precision matrix for Besag on ADMnew
matrix4inla <- nb2mat(nb_new_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_new_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse


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

proper1_noInt_ADMnew_formula <- sampled_counts ~ 1 + time_id +
  f(time_id.copy,
    model = "ar1",
    hyper = ar1_hyper) + 
  f(area_id, 
    model = "besagproper2",
    graph = Besag_prec_new_level,
    hyper = spatial_hyper_proper)





################################################################################

scenario_names_ADM1 = c("sc1", "sc3", "sc5", 
                        "sc7", "sc9", "sc11")

for(scenario_name in scenario_names_ADM1){
  print(scenario_name)
  
  tmp.df <- data.frame(years = 1:13,
                       median = rep(NA, 13),
                       l = rep(NA, 13),
                       u = rep(NA, 13))
  
  load(paste('./Simulated_data/', scenario_name, '/',
             scenario_name, '_data.RData',sep = ""))
  lambda <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda$sampled_counts = lambda.df$sampled_counts[, dataset_id] #Just a chosen data set
  
  # Set the last three years to unkown
  lambda[lambda$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda <- lambda[order(lambda$area_id, decreasing = F), ]
  rownames(lambda) <- 1:nrow(lambda)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda$area_id.copy <- lambda$area_id
  lambda$time_id.copy <- lambda$time_id
  
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(1, rep(NA, 12)), time_id = 1)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                            data = lambda,
                            family = "poisson",
                            E = E_it, #E_it
                            lincomb = lc_temporal,
                            verbose = T,
                            control.predictor = list(compute = TRUE, link = 1),       #For predictions
                            control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[1, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[1, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[1, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(NA, 1, rep(NA, 11)), time_id = 2)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[2, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[2, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[2, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 2), 1, rep(NA, 10)), time_id = 3)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[3, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[3, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[3, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 3), 1, rep(NA, 9)), time_id = 4)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[4, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[4, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[4, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 4), 1, rep(NA, 8)), time_id = 5)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[5, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[5, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[5, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 5), 1, rep(NA, 7)), time_id = 6)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[6, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[6, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[6, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 6), 1, rep(NA, 6)), time_id = 7)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[7, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[7, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[7, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 7), 1, rep(NA, 5)), time_id = 8)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[8, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[8, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[8, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 8), 1, rep(NA, 4)), time_id = 9)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[9, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[9, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[9, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 9), 1, rep(NA, 3)), time_id = 10)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[10, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[10, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[10, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 10), 1, rep(NA, 2)), time_id = 11)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[11, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[11, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[11, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 11), 1, rep(NA, 1)), time_id = 12)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[12, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[12, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[12, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 12), 1), time_id = 13)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[13, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[13, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[13, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  save(tmp.df,
       file = paste('./proper_structured_temporal_effect_', scenario_name, '.RData', sep = ""))
}

# load(paste('./proper_structured_temporal_effect_', "sc5", '.RData', sep = ""))
# to_plot_sc3.df = tmp.df
# load(paste('./proper_structured_temporal_effect_', "sc11", '.RData', sep = ""))
# to_plot_sc9.df = tmp.df
# 
# to_plot.df <- data.frame(years = c(1:13, 1:13),
#                          median = c(to_plot_sc3.df$median, to_plot_sc9.df$median),
#                          l = c(to_plot_sc3.df$l, to_plot_sc9.df$l),
#                          u = c(to_plot_sc3.df$u, to_plot_sc9.df$u),
#                          type = c(rep("short", 13), rep("long", 13)))
# 
# 
# ggplot(data = to_plot.df) + 
#   geom_line(aes(x = years, y = median, col = type)) + 
#   geom_line(aes(x = years, y = l, col = type)) + 
#   geom_line(aes(x = years, y = u, col = type))


scenario_names_ADM4 = c("sc2", "sc4", "sc6",
                        "sc8", "sc10", "sc12")

for(scenario_name in scenario_names_ADM4){
  print(scenario_name)
  
  tmp.df <- data.frame(years = 1:13,
                       median = rep(NA, 13),
                       l = rep(NA, 13),
                       u = rep(NA, 13))
  
  load(paste('./Simulated_data/', scenario_name, '/',
             scenario_name, '_data.RData',sep = ""))
  lambda <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda$sampled_counts = lambda.df$sampled_counts[, dataset_id_2] #Just a chosen data set
  
  # Set the last three years to unkown
  lambda[lambda$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda <- lambda[order(lambda$area_id, decreasing = F), ]
  rownames(lambda) <- 1:nrow(lambda)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda$area_id.copy <- lambda$area_id
  lambda$time_id.copy <- lambda$time_id
  
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(1, rep(NA, 12)), time_id = 1)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[1, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[1, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[1, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(NA, 1, rep(NA, 11)), time_id = 2)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[2, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[2, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[2, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 2), 1, rep(NA, 10)), time_id = 3)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[3, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[3, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[3, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 3), 1, rep(NA, 9)), time_id = 4)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[4, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[4, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[4, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 4), 1, rep(NA, 8)), time_id = 5)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[5, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[5, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[5, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 5), 1, rep(NA, 7)), time_id = 6)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[6, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[6, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[6, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 6), 1, rep(NA, 6)), time_id = 7)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[7, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[7, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[7, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 7), 1, rep(NA, 5)), time_id = 8)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[8, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[8, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[8, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 8), 1, rep(NA, 4)), time_id = 9)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[9, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[9, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[9, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 9), 1, rep(NA, 3)), time_id = 10)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[10, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[10, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[10, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 10), 1, rep(NA, 2)), time_id = 11)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[11, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[11, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[11, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 11), 1, rep(NA, 1)), time_id = 12)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[12, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[12, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[12, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  lc_temporal <- inla.make.lincomb(time_id.copy = c(rep(NA, 12), 1), time_id = 13)
  proper1_noInt_ADM1 <- inla(proper1_noInt_ADM4_formula, 
                             data = lambda,
                             family = "poisson",
                             E = E_it, #E_it
                             lincomb = lc_temporal,
                             verbose = T,
                             control.predictor = list(compute = TRUE, link = 1),       #For predictions
                             control.compute = list(config = TRUE, return.marginals.predictor=TRUE))
  
  tmp.df[13, ]$median = proper1_noInt_ADM1$summary.lincomb.derived$'0.5quant'
  tmp.df[13, ]$l = proper1_noInt_ADM1$summary.lincomb.derived$'0.025quant'
  tmp.df[13, ]$u = proper1_noInt_ADM1$summary.lincomb.derived$'0.975quant'
  
  save(tmp.df,
       file = paste('./proper_structured_temporal_effect_', scenario_name, '.RData', sep = ""))
}


load(paste('./proper_structured_temporal_effect_', "sc2", '.RData', sep = ""))
to_plot_sc4.df = tmp.df
load(paste('./proper_structured_temporal_effect_', "sc8", '.RData', sep = ""))
to_plot_sc10.df = tmp.df

to_plot.df <- data.frame(years = c(1:13, 1:13),
                         median = c(to_plot_sc4.df$median, to_plot_sc10.df$median),
                         l = c(to_plot_sc4.df$l, to_plot_sc10.df$l),
                         u = c(to_plot_sc4.df$u, to_plot_sc10.df$u),
                         type = c(rep("short", 13), rep("long", 13)))


ggplot(data = to_plot.df) + 
  geom_line(aes(x = years, y = median, col = type)) + 
  geom_line(aes(x = years, y = l, col = type)) + 
  geom_line(aes(x = years, y = u, col = type))



#for(){
#}








# lc_temporal <- inla.make.lincombs(time_id = rep(1, 13), 
#                                   time_id.copy = c(1, rep(NA, 12)),
#                                   time_id.copy = c(NA, 1, rep(NA, 11)),
#                                   time_id.copy = c(rep(NA, 2), 1, rep(NA, 10)))
lc_temporal <- inla.make.lincomb(time_id.copy = c(1, rep(NA, 12)), time_id = 1)



#lc  1 0.001232954 0.10828 -0.2486542 0.004157137  0.2415526 0.003126038 7.247115e-05
#lc  1 0.00827994 0.1043836 -0.2198033 0.005231767  0.2537075 0.006920592 8.258116e-05
#lc  1 0.04587463 0.1018499 -0.1812256 0.04235068  0.2711566 0.02512975 5.141632e-05

# time_id = c(1, 1, 1, 1, 1, 1,
#             1, 1, 1, 1, 1, 1,
#             1)


# proper1_noInt_ADM1 <- inla(proper1_noInt_ADM1_formula, 
#                            data = lambda, 
#                            family = "poisson",
#                            E = E_it, #E_it
#                            lincomb = lc_temporal_2,
#                            verbose = T,
#                            control.predictor = list(compute = TRUE, link = 1),       #For predictions
#                            control.compute = list(config = TRUE, return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
# 
# 
# plot(proper1_noInt_ADM1)
# 
# 
# 

# 
# 
# proper1_noInt_ADM1$summary.lincomb.derived
# 
# proper1_noInt_ADM1$marginals.lincomb.derived







