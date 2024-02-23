#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Preliminaries.R")
source("Utilities.R")


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


## Find mapping to each region in germany_map2
germany_map2_area_id_mapping <- st_join(yxt_geom, st_make_valid(germany_map_2),
                                        join = st_within)$ID_1

yxt_geom$germany_map2_area_id_mapping = germany_map2_area_id_mapping

## Find mapping to each region in germany_map
germany_map_area_id_mapping = rep(NA, nrow(yxt_geom))
germany_map_area_id_mapping_tmp <- st_join(yxt_geom, st_make_valid(germany_map),
                                           join = st_within)

## Drop NA-values
germany_map_area_id_mapping_tmp = drop_na(germany_map_area_id_mapping_tmp)

## Add area IDs to yxt_geom
germany_map_area_id_mapping[germany_map_area_id_mapping_tmp$unique_id] = germany_map_area_id_mapping_tmp$ID_2
yxt_geom$germany_map_area_id_mapping = germany_map_area_id_mapping


## Find mapping for the areas w.o. associatied points for germany_map 
missing_mapping_germany_map.df <- data.frame(missing_id = rep(NA, nrow(yxt_geom)),
                                             time_id = rep(NA, nrow(yxt_geom)),
                                             yxt_geom_unique_id = yxt_geom$unique_id)


### Find the regions within Germany where there is no corresponding grid point being evaluated
missing_areas_IDs = c()
for(i in unique(germany_map$ID_2)){
  tmp_ = yxt_geom[within_germany_indices &
                    yxt_geom$germany_map_area_id_mapping == i &
                    yxt_geom$time_id == 1, ]
  tmp_ = drop_na(tmp_)
  if(is.na(mean(tmp_$polygon_id))){
    missing_areas_IDs = c(missing_areas_IDs, i)
  }
}
missing_areas_IDs

### Find the grid point closest to the area, and use polygon_id to identify area, and use
### unique_id to map for each specific time
tmp_ = yxt_geom[within_germany_indices, ]
for(missing_area in missing_areas_IDs){
  for(j in 1:length(t_axis)){
    #print(t)
    tmp2_ = tmp_[tmp_$t == t_axis[j], ]
    index = st_nearest_feature(st_make_valid(germany_map[germany_map$ID_2 == missing_area, ]),
                               tmp2_)
    
    ## Must account for increase in time
    index = index + (j - 1) * nrow(tmp2_)
    
    ## Add values
    missing_mapping_germany_map.df[index, ]$missing_id = missing_area
    missing_mapping_germany_map.df[index, ]$time_id = tmp_[index, ]$time_id 
  }
}

missing_mapping_germany_map.df = drop_na(missing_mapping_germany_map.df)

no_problem_areas = unique(germany_map$ID_2)
no_problem_areas = no_problem_areas[!(no_problem_areas %in% missing_mapping_germany_map.df$missing_id)]

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
intercept = -2.2 #Intercept
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
intercept = -2.2 #Intercept
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
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[
          risk_surface.list$unique_id %in% 
            missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                    missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id, ]
    
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
intercept = -2.2 #Intercept
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
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
intercept = -2.2 #Intercept
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
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
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[
      risk_surface.list$unique_id %in% 
        missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                         missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id, ]
    
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
intercept = -2.2 #Intercept
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
intercept = -2.2 #Intercept
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
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[
      risk_surface.list$unique_id %in% 
        missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                         missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id, ]
    
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
intercept = -2.2 #Intercept
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
intercept = -2.2 #Intercept
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
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[
      risk_surface.list$unique_id %in% 
        missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                         missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id, ]
    
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
intercept = -2.2 #Intercept
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
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
intercept = -2.2 #Intercept
tau_s = 15; tau_t = 65 #1) spatial and 2) temporal smoothness parameters

## Get the risk-field
Lambda_st = simulate_risk_surface(Bst, Px, Py, Pt,
                                  Ix, Iy, It,
                                  tau_s, tau_t,
                                  intercept, 
                                  temporal_trend = 2, beta1_t = 0.02, 
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
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[
      risk_surface.list$unique_id %in% 
        missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                         missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id, ]
    
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
intercept = -2.2 #Intercept
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
intercept = -2.2 #Intercept
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
    
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
  for(i in unique(missing_mapping_germany_map.df$missing_id)){
    index = (t - 1) * nrow(germany_map) + i
    lambda.df[index, ]$area_id = i; lambda.df[index, ]$time_id = t
    
    tmp_ = risk_surface.list[
      risk_surface.list$unique_id %in% 
        missing_mapping_germany_map.df[missing_mapping_germany_map.df$missing_id == i & 
                                         missing_mapping_germany_map.df$time_id == t, ]$yxt_geom_unique_id, ]
    
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











