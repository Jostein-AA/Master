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

################################################################################
# Scenario 1: first-level admin map, const. temporal trend (w. greater temporal variation),
# greater spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_20, 
                                  Bs_20,
                                  Sigma_st_12,
                                  kt,
                                  ks_20,
                                  intercept, temporal_trend = 1,
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
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_1")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc1_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc1/sc1_data.RData")

## Create a txt-file for sc7 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc1/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc1.csv", sep = ""),
          row.names = F)



################################################################################
# Scenario 3: second_level_admin_map2, increasing temporal trend (w. smaller temporal variation),
# greater spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_20, 
                                  Bs_20,
                                  Sigma_st_3456,
                                  kt,
                                  ks_20,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.014, 
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
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_3")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc3_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc3/sc3_data.RData")

## Create a txt-file for sc7 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc3/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc3.csv", sep = ""),
          row.names = F)


################################################################################
# Scenario 5: second_level_admin_map2, change-point temporal trend (w. smaller temporal variation),
# greater spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_20, 
                                  Bs_20,
                                  Sigma_st_3456,
                                  kt,
                                  ks_20,
                                  intercept, 
                                  temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015,
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
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_5")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc5_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc5/sc5_data.RData")

## Create a txt-file for sc7 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc5/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc5.csv", sep = ""),
          row.names = F)

################################################################################
# Scenario 7: second_level_admin_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_10, 
                                  Bs_10,
                                  Sigma_st_78,
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
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_7_1")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc7_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc7/sc7_data_2.RData")

## Create a txt-file for sc7 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc7/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc7.csv", sep = ""),
          row.names = F)

################################################################################
# Scenario 9: second_level_admin_map2, increasing temporal trend (w. smaller temporal variation),
# smaller spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_10, 
                                  Bs_10,
                                  Sigma_st_9101112,
                                  kt,
                                  ks_10,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.014, 
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
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_9")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc9_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc9/sc9_data.RData")

## Create a txt-file for sc9 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc9/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc9.csv", sep = ""),
          row.names = F)

################################################################################
# Scenario 11: second_level_admin_map2, change-point temporal trend (w. smaller temporal variation),
# smaller spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_10, 
                                  Bs_10,
                                  Sigma_st_9101112,
                                  kt,
                                  ks_10,
                                  intercept, 
                                  temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015,
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
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    
    # Approximately the integral of the area and time continuous risk 
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_11")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc11_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc11/sc11_data.RData")

## Create a txt-file for sc7 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc11/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc11.csv", sep = ""),
          row.names = F)

