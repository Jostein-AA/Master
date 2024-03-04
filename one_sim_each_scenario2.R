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
load("tensor_prod_smooths_cov_matrices.RData")

# Define the intercept (is constant across all scenarios)
intercept = log(0.1)

################################################################################
# Scenario 1: second_level_admin_map2, const. temporal trend (w. greater temporal variation),
# greater spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_20, 
                                  Bs_20,
                                  Sigma_st_12,
                                  kt,
                                  ks_20,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)



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
                        lambda_it = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map2_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_1_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc1_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc1_data_1.RData")

################################################################################
# Scenario 2: second_level_admin_map, const. temporal trend (w. greater temporal variation),
# greater spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_20, 
                                  Bs_20,
                                  Sigma_st_12,
                                  kt,
                                  ks_20,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)



## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT))




for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}




lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_2_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc2_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc2_data_1.RData")



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
                                  n_sim = 1)


## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list


# Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(first_level_admin_map) * tT),
                        time_id = rep(0, nrow(first_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map2_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_3_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc3_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc3_data_1.RData")


################################################################################
# Scenario 4: second_level_admin_map, increasing temporal trend (w. smaller temporal variation),
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
                                  n_sim = 1)



## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}




lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_4_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc4_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc4_data_1.RData")

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
                                  n_sim = 1)


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
                        lambda_it = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map2_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_5_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc5_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc5_data_1.RData")


################################################################################
# Scenario 6: second_level_admin_map, change-point temporal trend (w. smaller temporal variation),
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
                                  n_sim = 1)


## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}




lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_6_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc6_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc6_data_1.RData")


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
                                  n_sim = 1)


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
                        lambda_it = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map2_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_7_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc7_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc7_data_1.RData")


################################################################################
# Scenario 8: second_level_admin_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst_10, 
                                  Bs_10,
                                  Sigma_st_78,
                                  kt,
                                  ks_10,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)


## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}


lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_8_1")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc8_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc8_data_1.RData")

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
                                  n_sim = 1)


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
                        lambda_it = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map2_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_9_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc9_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc9_data_1.RData")


################################################################################
# Scenario 10: second_level_admin_map, increasing temporal trend (w. smaller temporal variation),
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
                                  n_sim = 1)
## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}




lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_10_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc10_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc10_data_1.RData")


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
                                  n_sim = 1)


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
                        lambda_it = rep(0, nrow(first_level_admin_map) * tT),
                        E_it = rep(100, nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map2_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_11_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc11_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc11_data_1.RData")

################################################################################
# Scenario 12: second_level_admin_map, change-point temporal trend (w. smaller temporal variation),
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
                                  n_sim = 1)

## Add risk-values to yxt_grid w. geoms
risk_surface.list = yxt_geom
risk_surface.list$values = Lambda_st

## Drop points outside Germany
risk_surface.list = risk_surface.list[within_germany_indices, ]

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        lambda_it = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_admin_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}




lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})


print("Scenario_12_1")


save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc12_risk_surface_1.RData")

save(lambda.df,
     file = "./Simulated_data/sc12_data_1.RData")