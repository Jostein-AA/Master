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
# Scenario 2: second_level_admin_map, const. temporal trend (w. greater temporal variation),
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

## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(second_level_admin_map) * tT),
                        time_id = rep(0, nrow(second_level_admin_map) * tT),
                        E_it = rep(100, nrow(second_level_admin_map) * tT),
                        space.time = 1:(nrow(second_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                           missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it


## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_2")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc2_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc2/sc2_data.RData")

## Create a txt-file for sc2 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc2/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc2.csv", sep = ""),
          row.names = F)

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
                                  n_sim = n_sim)



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
                        E_it = rep(100, nrow(second_level_admin_map) * tT),
                        space.time = 1:(nrow(second_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                           missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it


## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_4")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc4_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc4/sc4_data.RData")

## Create a txt-file for sc4 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc4/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc4.csv", sep = ""),
          row.names = F)

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
                                  n_sim = n_sim)


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
                        E_it = rep(100, nrow(second_level_admin_map) * tT),
                        space.time = 1:(nrow(second_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                           missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it


## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_6")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc6_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc6/sc6_data.RData")

## Create a txt-file for sc6 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc6/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc6.csv", sep = ""),
          row.names = F)

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
                                  n_sim = n_sim)


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
                        E_it = rep(100, nrow(second_level_admin_map) * tT),
                        space.time = 1:(nrow(second_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                           missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it


## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_8")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc8_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc8/sc8_data.RData")

## Create a txt-file for sc8 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc8/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc8.csv", sep = ""),
          row.names = F)

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
                                  n_sim = n_sim)
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
                        E_it = rep(100, nrow(second_level_admin_map) * tT),
                        space.time = 1:(nrow(second_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                           missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it


## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_10")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc10_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc10/sc10_data.RData")

## Create a txt-file for sc10 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc10/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc10.csv", sep = ""),
          row.names = F)

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
                                  n_sim = n_sim)

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
                        E_it = rep(100, nrow(second_level_admin_map) * tT),
                        space.time = 1:(nrow(second_level_admin_map) * tT))


lambda_it = matrix(nrow = nrow(lambda.df), 
                   ncol = dim(risk_surface.list$values)[2])


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$second_level_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
  for(i in unique(missing_mapping_second_level_admin_map.df$missing_id)){
    index = (t - 1) * nrow(second_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$missing_id == i & 
                                                                           missing_mapping_second_level_admin_map.df$time_id == t, ]$yxt_geom_unique_id,
    ]
    
    #lambda.df[index, ]$lambda_it = mean(tmp_$values)
    lambda_it[index, ] = colMeans(tmp_$values)
  }
}

## Insert the area and time specific rates
lambda.df$lambda_it = lambda_it

## Calculate the area and time specific expected number of counts
lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it


## Sample counts
lambda.df$sampled_counts = sample_counts(lambda.df)


print("Scenario_12")

save(risk_surface.list, 
     file = "./Simulated_risk_surfaces/sc12_risk_surfaces.RData")

save(lambda.df,
     file = "./Simulated_data/sc12/sc12_data.RData")

## Create a txt-file for sc12 stating how many data sets there are
## For the analysis of each model, a separate file stating which data sets they
## have analyzed will also be made
## This way a full overview of what has been analyzed is available

dir = "./Simulated_data/sc12/"
overview.df = data.frame(data_sets_id = 1:dim(lambda.df$sampled_counts)[2])

write.csv(overview.df, paste(dir, "data_sets_sc12.csv", sep = ""),
          row.names = F)