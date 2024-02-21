#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Preliminaries.R")
source("Utilities.R")

################################################################################
# Define the B-spline basis functions in 3-dimensions and the marginal penalization
# matrices 

## Get a polygon_grid over Germany
polygon_grid = st_as_sf(st_make_grid(germany_border, 
                                     n = c(125, 125), 
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
t_axis <- seq(0.0001, 13, length.out = 57); xt_ = seq(0, 13, length.out = 57)

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
# Scenario 1: germany_map2, const. temporal trend, greater spatial variation


## Define parameters for simulation
intercept = -7.6 #Intercept
tau_s = 50; tau_t = 150 #1) spatial and 2) temporal smoothness parameters

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

save(Lambda_st,
     lambda.df,
     risk_surface.list,
     file = "scenario_1_1.RData")

################################################################################
# Scenario 2: germany_map, const. temporal trend, greater spatial variation















