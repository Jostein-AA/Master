#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")


load("maps_and_nb.RData")
load("grids_and_mappings.RData")
load("temporal_B_spline")
load("spatial_B_splines.RData")
load("B_spline_basis.RData")
load("penalization_matrices.RData")
load("tensor_prod_smooths_cov_matrices.RData")

# Define the intercept (is constant across all scenarios)
intercept = log(0.1)

################################################################################
# Define the B-spline basis functions in 3-dimensions, with more knots for
# more spatial variation. Also define the marginal penalization matrices

## Define number of knots in 1) the spatial directions and 2) temporal direction 
n_knots_s = 20; n_knots_t = 6

## define B-spline basis degree
bdeg = 3

# Specify the knot locations and create the B-spline basis design matrices
## Find the distance between each knot in x- and y-directions
disx <- (max(x) - min(x)) / n_knots_s; disy <- (max(y) - min(y)) / n_knots_s

## Find the left (l) and right (r) end and a bitmore spline not 0 at ends!
xl <- min(x) - disx * 0.05; yl <- min(y) - disy * 0.05
xr <- max(x) + disx * 0.05; yr <- max(y) + disy * 0.05

## This is distance between each knot?
dx <- (xr - xl) / n_knots_s; dy <- (yr - yl) / n_knots_s

## Specify the knot locations in x1 and x2 direction
knotsx <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
knotsy <- seq(yl - bdeg * dy, yr + bdeg * dy, by = dy)

## Make marginal design matrices for B-splines in x and y direction 
Bx <- spline.des(knotsx, x, bdeg + 1, 0 * x)$design  ##Base marginal x
By <- spline.des(knotsy, y, bdeg + 1, 0 * y)$design  ##Base marginal y

## Perform row-wise Kronecker product to get the spatial B-spline basis
Bs = row_wise_Kronecker(Bx, By)

## Find dimensions (spatial)
kx <- dim(Bx)[2]; ky <- dim(By)[2]; ks <- dim(Bs)[2]

## Find distance 
dist <- (max(xt) - min(xt)) / n_knots_t
xtl <- min(xt) - dist * 0.05; xtr <- max(xt) + dist * 0.05
dxt <- (xtr - xtl) / n_knots_t

## Define knot locations in time
knots <- seq(xtl - bdeg * dxt, xtr + bdeg * dxt, by = dxt)

## Define temporal marginal B-spline basis design matrix
Bt <- spline.des(knots, xt, bdeg + 1)$design
kt <- ncol(Bt)

## Define the 3-dimensional B-spline basis functions
Bst <- kronecker(Bt, Bs)

# Define the penalization matrices and get the risk field
## Define the marginal penalization matrices (RW1 precision matrices)
Px <- INLA:::inla.rw(n = kx, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Py <- INLA:::inla.rw(n = ky, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Pt <- INLA:::inla.rw(n = kt, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

## Identity matrices k1 x k1, k2 x k2, and kt x kt
Ix <- diag(kx); Iy <- diag(ky); It <- diag(kt)

print("---\n!All preliminaries are done!\n----")

################################################################################
# Scenario 1: germany_map2, const. temporal trend (w. greater temporal variation),
# greater spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 35 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(100, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map2_area_id_mapping == i &
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
# Scenario 2: germany_map, const. temporal trend (w. greater temporal variation),
# greater spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 35 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(100, nrow(germany_map) * tT))




for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                                                missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id,
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
# Scenario 3: germany_map2, increasing temporal trend (w. smaller temporal variation),
# greater spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(100, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map2_area_id_mapping == i &
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
# Scenario 4: germany_map, increasing temporal trend (w. smaller temporal variation),
# greater spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(100, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                                                missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id,
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
# Scenario 5: germany_map2, change-point temporal trend (w. smaller temporal variation),
# greater spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(100, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map2_area_id_mapping == i &
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
# Scenario 6: germany_map, change-point temporal trend (w. smaller temporal variation),
# greater spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(100, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                                                missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id,
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
# Change the B-spline basis functions to have fewer knots in space to lower
# spatial variation

## Remove the old B-spline matrix to ease computations
rm(Bst)

# Define the B-spline basis functions in 3-dimensions, with more knots for
# more spatial variation. Also define the marginal penalization matrices

## Define number of knots in 1) the spatial directions and 2) temporal direction 
n_knots_s = 10; n_knots_t = 6

## define B-spline basis degree
bdeg = 3

# Specify the knot locations and create the B-spline basis design matrices
## Find the distance between each knot in x- and y-directions
disx <- (max(x) - min(x)) / n_knots_s; disy <- (max(y) - min(y)) / n_knots_s

## Find the left (l) and right (r) end and a bitmore spline not 0 at ends!
xl <- min(x) - disx * 0.05; yl <- min(y) - disy * 0.05
xr <- max(x) + disx * 0.05; yr <- max(y) + disy * 0.05

## This is distance between each knot?
dx <- (xr - xl) / n_knots_s; dy <- (yr - yl) / n_knots_s

## Specify the knot locations in x1 and x2 direction
knotsx <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
knotsy <- seq(yl - bdeg * dy, yr + bdeg * dy, by = dy)

## Make marginal design matrices for B-splines in x and y direction 
Bx <- spline.des(knotsx, x, bdeg + 1, 0 * x)$design  ##Base marginal x
By <- spline.des(knotsy, y, bdeg + 1, 0 * y)$design  ##Base marginal y

## Perform row-wise Kronecker product to get the spatial B-spline basis
Bs = row_wise_Kronecker(Bx, By)

## Find dimensions (spatial)
kx <- dim(Bx)[2]; ky <- dim(By)[2]; ks <- dim(Bs)[2]

## Find distance 
dist <- (max(xt) - min(xt)) / n_knots_t
xtl <- min(xt) - dist * 0.05; xtr <- max(xt) + dist * 0.05
dxt <- (xtr - xtl) / n_knots_t

## Define knot locations in time
knots <- seq(xtl - bdeg * dxt, xtr + bdeg * dxt, by = dxt)

## Define temporal marginal B-spline basis design matrix
Bt <- spline.des(knots, xt, bdeg + 1)$design
kt <- ncol(Bt)

## Define the 3-dimensional B-spline basis functions
Bst <- kronecker(Bt, Bs)

# Define the penalization matrices and get the risk field
## Define the marginal penalization matrices (RW1 precision matrices)
Px <- INLA:::inla.rw(n = kx, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Py <- INLA:::inla.rw(n = ky, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Pt <- INLA:::inla.rw(n = kt, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

## Identity matrices k1 x k1, k2 x k2, and kt x kt
Ix <- diag(kx); Iy <- diag(ky); It <- diag(kt)


################################################################################
# Scenario 7: germany_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 35 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(100, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map2_area_id_mapping == i &
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
# Scenario 8: germany_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 35 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(100, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
        missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                      missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id,
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
# Scenario 9: germany_map2, increasing temporal trend (w. smaller temporal variation),
# smaller spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(100, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map2_area_id_mapping == i &
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
# Scenario 10: germany_map, increasing temporal trend (w. smaller temporal variation),
# smaller spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(100, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                                                missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id,
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
# Scenario 11: germany_map2, change-point temporal trend (w. smaller temporal variation),
# smaller spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(100, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map2_area_id_mapping == i &
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
# Scenario 12: germany_map, change-point temporal trend (w. smaller temporal variation),
# smaller spatial variation

## Define parameters for simulation
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
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
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(100, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in no_problem_areas){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$germany_map_area_id_mapping == i &
                               risk_surface.list$time_id == t, ]
    
    tmp_ = drop_na(tmp_)
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$unique_id %in% 
                               missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                                                missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id,
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











