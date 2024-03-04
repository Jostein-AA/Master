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

tactical error: to remove (whats below)
###########################
## Get a polygon_grid over Germany
polygon_grid = st_as_sf(st_make_grid(germany_border, 
                                     n = c(110, 110), 
                                     square = FALSE),
                        crs = crs_$srid)

## Add unique polygon_id
polygon_grid$polygon_id = 1:nrow(polygon_grid)

## Remove polygons outside of Germany
polygon_grid2 <- sf::st_intersection(polygon_grid, 
                                     st_make_valid(germany_border))

## Get the centroids of each polygon in the grid
centroids_grid <- st_centroid(polygon_grid)

## Get the coordinates of the centroids
gpts = sf::st_coordinates(centroids_grid)

## standardize x and y
locs <- apply(gpts, 2, function(x) (x - min(x)) / (max(x) - min(x)))

## Extract coordinates lat/lon
x <- locs[, 1]
y <- locs[, 2]

## Make t_axis: spatial domain (which discretely is: 1,...,13, and cont [0, 13])
t_axis <- seq(0.0001, 13, length.out = 40); xt_ = seq(0, 13, length.out = 40)

## Scale the temporal 'coordinates'
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Transform yxt_grid to have geometry points
yxt_geom = st_as_sf(yxt_grid, 
                    coords = c("x", "y"),
                    crs = crs_$srid)

yxt_geom$x = yxt_grid$x; yxt_geom$y = yxt_grid$y
yxt_geom$unique_id = 1:nrow(yxt_geom)


## Add time-id to yxt_geom
yxt_geom$time_id = ceiling(yxt_geom$t)

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Find indices within Germany
within_germany_indices = as.numeric(rownames(sf::st_intersection(yxt_geom, 
                                                                 st_make_valid(germany_border))))


## Find mapping to each region in second_level_admin_map2
first_level_area_id_mapping <- st_join(yxt_geom, st_make_valid(first_level_admin_map),
                                       join = st_within)$ID_1

yxt_geom$first_level_area_id_mapping = first_level_area_id_mapping

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
###########################

tactical error: to remove ^^


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
                                  n_sim = 100)



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
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
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
     file = "./Simulated_risk_surfaces/sc1_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc1_data_2.RData")


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
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
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
     file = "./Simulated_risk_surfaces/sc3_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc3_data_2.RData")


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
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
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
     file = "./Simulated_risk_surfaces/sc5_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc5_data_2.RData")

################################################################################
# Scenario 7: second_level_admin_map, const. temporal trend (w. greater temporal variation),
# smaller spatial variation

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, 
                                  Bs,
                                  Sigma_st_78,
                                  kt,
                                  ks_10,
                                  intercept, temporal_trend = 1,
                                  n_sim = 4)


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
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
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
     file = "./Simulated_risk_surfaces/sc7_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc7_data_2.RData")


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
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
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
     file = "./Simulated_risk_surfaces/sc9_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc9_data_2.RData")


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
                        E_it = rep(100, nrow(first_level_admin_map) * tT),
                        space.time = 1:(nrow(first_level_admin_map) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(first_level_admin_map)){
    index = (t - 1) * nrow(first_level_admin_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[risk_surface.list$first_level_area_id_mapping == i &
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
     file = "./Simulated_risk_surfaces/sc11_risk_surface_2.RData")

save(lambda.df,
     file = "./Simulated_data/sc11_data_2.RData")

