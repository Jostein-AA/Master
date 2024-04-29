#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")


load("maps_and_nb.RData")
load("grids_and_mappings.RData")
load("temporal_B_spline.RData")
load("spatial_B_splines.RData")
load("B_spline_basis.RData")
load("penalization_matrices.RData")
load("scaled_tensor_prod_smooths_cov_matrices.RData")

# Define the intercept (is constant across all scenarios)
intercept = log(0.1)

# Define the number of simulated data sets per scenario
n_sim = 100

# Seed
seed = 2276543

################################################################################
# Scenario 1: first-level admin map, const. temporal trend (w. greater temporal variation),
# greater spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
Lambda_st = simulate_risk_surface(seed,
                                  Bst_20, 
                                  Bs_20,
                                  Sigma_st_13,
                                  kt,
                                  ks_20,
                                  intercept, temporal_trend = 1,
                                  sig_st = 0.05,
                                  n_sim = n_sim)


## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(new_map) * tT),
                        time_id = rep(0, nrow(new_map) * tT),
                        E_it = rep(100, nrow(new_map) * tT),
                        space.time = 1:(nrow(new_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(new_map)){
    index = (t - 1) * nrow(new_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$new_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
    print(lambda_it[index, ])
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

# Update seed
seed = seed + 1

## Sample counts
lambda.df$sampled_counts = sample_counts(seed, lambda.df)


print("Scenario_13")

save(risk_surface.list,
     file = "./Simulated_risk_surfaces/sc13.RData")

save(lambda.df,
     file = "./Simulated_data/sc13/sc13_data.RData")

################################################################################
# Scenario 14: second_level_admin_map2, increasing temporal trend (w. smaller temporal variation),
# greater spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
Lambda_st = simulate_risk_surface(seed,
                                  Bst_20, 
                                  Bs_20,
                                  Sigma_st_1415,
                                  kt,
                                  ks_20,
                                  intercept, temporal_trend = 2,
                                  sig_st = 0.05, beta1_t = 0.014, 
                                  t_axis = t_axis,
                                  n_sim = n_sim)


## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(new_map) * tT),
                        time_id = rep(0, nrow(new_map) * tT),
                        E_it = rep(100, nrow(new_map) * tT),
                        space.time = 1:(nrow(new_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(new_map)){
    index = (t - 1) * nrow(new_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$new_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

# Update seed
seed = seed + 1

## Sample counts
lambda.df$sampled_counts = sample_counts(seed, lambda.df)


print("Scenario_14")

save(risk_surface.list,
     file = "./Simulated_risk_surfaces/sc14.RData")


save(lambda.df,
     file = "./Simulated_data/sc14/sc14_data.RData")



################################################################################
# Scenario 15: second_level_admin_map2, change-point temporal trend (w. smaller temporal variation),
# greater spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
Lambda_st = simulate_risk_surface(seed,
                                  Bst_20, 
                                  Bs_20,
                                  Sigma_st_1415,
                                  kt,
                                  ks_20,
                                  intercept, temporal_trend = 3,
                                  sig_st = 0.05, beta1_t = 0.02, beta2_t = 0.015, 
                                  t_axis = t_axis,
                                  n_sim = n_sim)


## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(new_map) * tT),
                        time_id = rep(0, nrow(new_map) * tT),
                        E_it = rep(100, nrow(new_map) * tT),
                        space.time = 1:(nrow(new_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(new_map)){
    index = (t - 1) * nrow(new_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$new_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

# Update seed
seed = seed + 1

## Sample counts
lambda.df$sampled_counts = sample_counts(seed, lambda.df)


print("Scenario_15")

save(risk_surface.list,
     file = "./Simulated_risk_surfaces/sc15.RData")


save(lambda.df,
     file = "./Simulated_data/sc15/sc15_data.RData")


################################################################################
# Scenario 16: second_level_admin_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
Lambda_st = simulate_risk_surface(seed,
                                  Bst_10, 
                                  Bs_10,
                                  Sigma_st_16,
                                  kt,
                                  ks_10,
                                  intercept, temporal_trend = 1,
                                  n_sim = n_sim)


## add the risk-values to yxt_grid
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(new_map) * tT),
                        time_id = rep(0, nrow(new_map) * tT),
                        E_it = rep(100, nrow(new_map) * tT),
                        space.time = 1:(nrow(new_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(new_map)){
    index = (t - 1) * nrow(new_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$new_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

# Update seed
seed = seed + 1

## Sample counts
lambda.df$sampled_counts = sample_counts(seed, lambda.df)


print("Scenario_16")

save(risk_surface.list,
     file = "./Simulated_risk_surfaces/sc16.RData")

save(lambda.df,
     file = "./Simulated_data/sc16/sc16_data.RData")


################################################################################
# Scenario 17: second_level_admin_map2, increasing temporal trend (w. smaller temporal variation),
# smaller spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
Lambda_st = simulate_risk_surface(seed,
                                  Bst_10, 
                                  Bs_10,
                                  Sigma_st_1718,
                                  kt,
                                  ks_10,
                                  intercept, temporal_trend = 2, beta1_t = 0.014, 
                                  t_axis = t_axis,
                                  n_sim = n_sim)


## add the risk-values to yxt_grid
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(new_map) * tT),
                        time_id = rep(0, nrow(new_map) * tT),
                        E_it = rep(100, nrow(new_map) * tT),
                        space.time = 1:(nrow(new_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(new_map)){
    index = (t - 1) * nrow(new_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$new_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

# Update seed
seed = seed + 1

## Sample counts
lambda.df$sampled_counts = sample_counts(seed, lambda.df)


print("Scenario_17")

save(risk_surface.list,
     file = "./Simulated_risk_surfaces/sc17.RData")


save(lambda.df,
     file = "./Simulated_data/sc17/sc17_data.RData")



################################################################################
# Scenario 18: second_level_admin_map2, change-point temporal trend (w. smaller temporal variation),
# smaller spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
Lambda_st = simulate_risk_surface(seed,
                                  Bst_10, 
                                  Bs_10,
                                  Sigma_st_1718,
                                  kt,
                                  ks_10,
                                  intercept, temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015, 
                                  t_axis = t_axis,
                                  n_sim = n_sim)


## add the risk-values to yxt_grid
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(new_map) * tT),
                        time_id = rep(0, nrow(new_map) * tT),
                        E_it = rep(100, nrow(new_map) * tT),
                        space.time = 1:(nrow(new_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(new_map)){
    index = (t - 1) * nrow(new_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$new_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

# Update seed
seed = seed + 1

## Sample counts
lambda.df$sampled_counts = sample_counts(seed, lambda.df)


print("Scenario_18")

save(risk_surface.list,
     file = "./Simulated_risk_surfaces/sc18.RData")


save(lambda.df,
     file = "./Simulated_data/sc18/sc18_data.RData")

