#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Preliminaries.R")
source("Utilities.R")


## Get a polygon_grid over Germany
polygon_grid = st_as_sf(st_make_grid(germany_border, 
                                     n = c(85, 85), 
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

## Define spatial domain (which discretely is: 1,...,13, and cont [0, 13])
t_axis <- seq(0.0001, 13, length.out = 50); xt_ = seq(0, 13, length.out = 50)

## Scale the temporal 'coordinates'
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

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
# Scenario 1: germany_map2, const. temporal trend (w. greater temporal variation),
# greater spatial variation

## Define parameters for simulation
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 75 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = "./validation_plots/test/scenario_1/")



## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, st_make_valid(germany_map_2),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_1",
                                             "geometry")]


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(1E4, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_1 == i &
                                risk_surface.list2$time_id == t, ]
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 75 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
#plots_for_GIF(risk_surface.list = risk_surface.list,
#              polygons = polygon_grid2,
#              t_axis = t_axis,
#              filename_base = "./validation_plots/test/scenario_2/")

## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, 
                              st_make_valid(germany_map),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_2",
                                             "geometry")]


## Drop NA-values
risk_surface.list2 = drop_na(risk_surface.list2)


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(1E4, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in unique(germany_map$ID_2)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_2 == i &
                                risk_surface.list2$time_id == t, ]
    if(is.na(mean(tmp_$values))){
      ## Must find closest points for these areas
      tmp_ = st_join(risk_surface.list[risk_surface.list$time_id == t,], 
                     st_make_valid(germany_map[germany_map$ID_2 == i, ]),
                     join = st_nearest_feature)
      
    }
    
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = "./validation_plots/test/scenario_3/")



## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, st_make_valid(germany_map_2),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_1",
                                             "geometry")]


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(1E4, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_1 == i &
                                risk_surface.list2$time_id == t, ]
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)

## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, 
                              st_make_valid(germany_map),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_2",
                                             "geometry")]


## Drop NA-values
risk_surface.list2 = drop_na(risk_surface.list2)


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(1E4, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in unique(germany_map$ID_2)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_2 == i &
                                risk_surface.list2$time_id == t, ]
    if(is.na(mean(tmp_$values))){
      ## Must find closest points for these areas
      tmp_ = st_join(risk_surface.list[risk_surface.list$time_id == t,], 
                     st_make_valid(germany_map[germany_map$ID_2 == i, ]),
                     join = st_nearest_feature)
      
    }
    
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015,
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = "./validation_plots/test/scenario_5/")



## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, st_make_valid(germany_map_2),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_1",
                                             "geometry")]


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(1E4, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_1 == i &
                                risk_surface.list2$time_id == t, ]
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015,
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)

## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, 
                              st_make_valid(germany_map),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_2",
                                             "geometry")]


## Drop NA-values
risk_surface.list2 = drop_na(risk_surface.list2)


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(1E4, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in unique(germany_map$ID_2)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_2 == i &
                                risk_surface.list2$time_id == t, ]
    if(is.na(mean(tmp_$values))){
      ## Must find closest points for these areas
      tmp_ = st_join(risk_surface.list[risk_surface.list$time_id == t,], 
                     st_make_valid(germany_map[germany_map$ID_2 == i, ]),
                     join = st_nearest_feature)
      
    }
    
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

## Define spatial domain (which discretely is: 1,...,13, and cont [0, 13])
t_axis <- seq(0.0001, 13, length.out = 50); xt_ = seq(0, 13, length.out = 50)

## Scale the temporal 'coordinates'
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 75 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = "./validation_plots/test/scenario_7/")



## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, st_make_valid(germany_map_2),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_1",
                                             "geometry")]


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(1E4, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_1 == i &
                                risk_surface.list2$time_id == t, ]
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 75 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, temporal_trend = 1,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
#plots_for_GIF(risk_surface.list = risk_surface.list,
#              polygons = polygon_grid2,
#              t_axis = t_axis,
#              filename_base = "./validation_plots/test/scenario_2/")

## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, 
                              st_make_valid(germany_map),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_2",
                                             "geometry")]


## Drop NA-values
risk_surface.list2 = drop_na(risk_surface.list2)


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(1E4, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in unique(germany_map$ID_2)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_2 == i &
                                risk_surface.list2$time_id == t, ]
    if(is.na(mean(tmp_$values))){
      ## Must find closest points for these areas
      tmp_ = st_join(risk_surface.list[risk_surface.list$time_id == t,], 
                     st_make_valid(germany_map[germany_map$ID_2 == i, ]),
                     join = st_nearest_feature)
      
    }
    
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = "./validation_plots/test/scenario_9/")



## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, st_make_valid(germany_map_2),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_1",
                                             "geometry")]


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(1E4, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_1 == i &
                                risk_surface.list2$time_id == t, ]
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)

## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, 
                              st_make_valid(germany_map),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_2",
                                             "geometry")]


## Drop NA-values
risk_surface.list2 = drop_na(risk_surface.list2)


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(1E4, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in unique(germany_map$ID_2)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_2 == i &
                                risk_surface.list2$time_id == t, ]
    if(is.na(mean(tmp_$values))){
      ## Must find closest points for these areas
      tmp_ = st_join(risk_surface.list[risk_surface.list$time_id == t,], 
                     st_make_valid(germany_map[germany_map$ID_2 == i, ]),
                     join = st_nearest_feature)
      
    }
    
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015,
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)


## Save plots for validation of risk-surface
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = "./validation_plots/test/scenario_11/")



## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, st_make_valid(germany_map_2),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_1",
                                             "geometry")]


#Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map_2) * tT),
                        time_id = rep(0, nrow(germany_map_2) * tT),
                        lambda_it = rep(0, nrow(germany_map_2) * tT),
                        E_it = rep(1E4, nrow(germany_map_2) * tT))

for(t in 1:tT){
  print(t)
  for(i in 1:nrow(germany_map_2)){
    index = (t - 1) * nrow(germany_map_2) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_1 == i &
                                risk_surface.list2$time_id == t, ]
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
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 3, 
                                  beta1_t = 0.02, beta2_t = 0.015,
                                  t_axis = t_axis,
                                  n_sim = 1)

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y

## Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

## Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t)

## Get which points within areas of Germany
risk_surface.list2 <- st_join(risk_surface.list, 
                              st_make_valid(germany_map),
                              join = st_within)

risk_surface.list2 <- risk_surface.list2[, c("t", "values", "polygon_id", 
                                             "x", "y", "time_id", "ID_2",
                                             "geometry")]


## Drop NA-values
risk_surface.list2 = drop_na(risk_surface.list2)


## Integrate to get the lambda_it values
lambda.df <- data.frame(area_id = rep(0, nrow(germany_map) * tT),
                        time_id = rep(0, nrow(germany_map) * tT),
                        lambda_it = rep(0, nrow(germany_map) * tT),
                        E_it = rep(1E4, nrow(germany_map) * tT))


for(t in 1:tT){
  print(c("t: ", toString(t)))
  for(i in unique(germany_map$ID_2)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list2[risk_surface.list2$ID_2 == i &
                                risk_surface.list2$time_id == t, ]
    if(is.na(mean(tmp_$values))){
      ## Must find closest points for these areas
      tmp_ = st_join(risk_surface.list[risk_surface.list$time_id == t,], 
                     st_make_valid(germany_map[germany_map$ID_2 == i, ]),
                     join = st_nearest_feature)
      
    }
    
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











