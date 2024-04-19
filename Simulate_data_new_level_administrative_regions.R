#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")


load("maps_and_nb.RData")
load("grids_and_mappings.RData")
#load("temporal_B_spline.RData")
#load("spatial_B_splines.RData")
#load("B_spline_basis.RData")
#load("penalization_matrices.RData")
#load("scaled_tensor_prod_smooths_cov_matrices.RData")

# Define the intercept (is constant across all scenarios)
intercept = log(0.1)

# Define the number of simulated data sets per scenario
n_sim = 100

# Seed
seed = 134563

################################################################################
# Scenario 1: first-level admin map, const. temporal trend (w. greater temporal variation),
# greater spatial variation

# Update seed
seed = seed + 1

## Get the risk-field
#Lambda_st = simulate_risk_surface(seed,
#                                  Bst_20, 
#                                  Bs_20,
#                                  Sigma_st_12,
#                                  kt,
#                                  ks_20,
#                                  intercept, temporal_trend = 1,
#                                  n_sim = n_sim)



## Add risk-values to yxt_grid w. geoms
#risk_surface.list = yxt_geom
#risk_surface.list$values = Lambda_st

## Drop points outside Germany
#risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
#indices_risk_surface.list = nrow(risk_surface.list)
#rownames(risk_surface.list) = 1:indices_risk_surface.list

load("Simulated_risk_surfaces/sc1_risk_surfaces.RData")

yxt_geom_ = yxt_geom

## Drop points outside Germany
yxt_geom_ = yxt_geom_[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(yxt_geom_)
rownames(yxt_geom_) = 1:indices_risk_surface.list

## Add new map area-id
risk_surface.list$new_level_area_id_mapping = yxt_geom_$new_level_area_id_mapping





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


print("Scenario_1")


save(lambda.df,
     file = "./Simulated_data/sc13/sc13_data.RData")

################################################################################
# Scenario 3: second_level_admin_map2, increasing temporal trend (w. smaller temporal variation),
# greater spatial variation

# Update seed
seed = seed + 1

load("Simulated_risk_surfaces/sc3_risk_surfaces.RData")

## Add new map area-id
risk_surface.list$new_level_area_id_mapping = yxt_geom_$new_level_area_id_mapping


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


print("Scenario_3")
save(lambda.df,
     file = "./Simulated_data/sc14/sc14_data.RData")



################################################################################
# Scenario 5: second_level_admin_map2, change-point temporal trend (w. smaller temporal variation),
# greater spatial variation

# Update seed
seed = seed + 1

load("Simulated_risk_surfaces/sc5_risk_surfaces.RData")

## Add new map area-id
risk_surface.list$new_level_area_id_mapping = yxt_geom_$new_level_area_id_mapping


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


print("Scenario_5")


save(lambda.df,
     file = "./Simulated_data/sc15/sc15_data.RData")


################################################################################
# Scenario 7: second_level_admin_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

# Update seed
seed = seed + 1


load("Simulated_risk_surfaces/sc7_risk_surfaces.RData")

## Add new map area-id
risk_surface.list$new_level_area_id_mapping = yxt_geom_$new_level_area_id_mapping


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


print("Scenario_7")


save(lambda.df,
     file = "./Simulated_data/sc16/sc16_data.RData")


################################################################################
# Scenario 9: second_level_admin_map2, increasing temporal trend (w. smaller temporal variation),
# smaller spatial variation

# Update seed
seed = seed + 1


load("Simulated_risk_surfaces/sc9_risk_surfaces.RData")

## Add new map area-id
risk_surface.list$new_level_area_id_mapping = yxt_geom_$new_level_area_id_mapping


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


print("Scenario_9")


save(lambda.df,
     file = "./Simulated_data/sc17/sc17_data.RData")



################################################################################
# Scenario 11: second_level_admin_map2, change-point temporal trend (w. smaller temporal variation),
# smaller spatial variation

# Update seed
seed = seed + 1


load("Simulated_risk_surfaces/sc11_risk_surfaces.RData")

## Add new map area-id
risk_surface.list$new_level_area_id_mapping = yxt_geom_$new_level_area_id_mapping

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


print("Scenario_11")


save(lambda.df,
     file = "./Simulated_data/sc18/sc18_data.RData")

