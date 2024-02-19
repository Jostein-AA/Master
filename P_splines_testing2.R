#Clear environment
rm(list = ls())

#Load libraries
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)
library(pspline)
library(readxl)
library(OneR)
library(mgcv)
library(splines)
library(ggfortify)
library(crs)
library(magick)
library(INLA)
library(bigDM)
library(MASS)
library(RColorBrewer)
library(tmap)

#Load in data
source('Preliminaries.R')
source("Utilities.R")

#Load in structures
nb <- spdep::poly2nb(germany_map, queen = FALSE)
nb2 <- spdep::poly2nb(st_make_valid(germany_map_2), queen = FALSE)
#Get the coordinate system used
crs_ = st_crs(germany_border, parameters = TRUE)

### Define marginal basis functions
#Preliminaries for tensor product smooth using 3 P-splines

#Extract a spatial domain to find [x-min, x-max] X [y-min, y-max]
gpts = sf::st_coordinates(germany_border$geometry)

## 1) Find the x-min an x-max 2) Find y-min and y-max
x_min = min(gpts[,c("X")]); x_max = max(gpts[,c("X")])
y_min = min(gpts[,c("Y")]); y_max = max(gpts[,c("Y")])

# Define number of knots in 1) the spatial directions 2) temporal direction 
n_knots_s = 10; n_knots_t = 5

#define B-spline basis degree
bdeg = 3

## Define the spatial grid xy_grid and standardize the coordinates (important!)
x_axis = seq(x_min, x_max, length.out = 75) #75 points to evaluate on x-axis
y_axis = seq(y_min, y_max, length.out = 75) #75 points to evaluate on y-axis
#xy-grid containing 1000 points to evaluate the spatial domain on
xy_grid = expand.grid(x = x_axis, y = y_axis)        

#standardize x and y
locs <- apply(xy_grid, 2, function(x) (x - min(x)) / (max(x) - min(x)))

#Extract coordinates lat/lon
x <- locs[, 1]
y <- locs[, 2]

#Specify the knot locations and create the B-spline basis design matrices
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

#Make marginal design matrices for B-splines in x and y direction 
Bx <- spline.des(knotsx, x, bdeg + 1, 0 * x)$design  ##Base marginal x
By <- spline.des(knotsy, y, bdeg + 1, 0 * y)$design  ##Base marginal y

#Perform row-wise Kronecker product to get the spatial B-spline basis
Bs = row_wise_Kronecker(Bx, By)

#Find dimensions (spatial)
kx <- dim(Bx)[2]; ky <- dim(By)[2]; ks <- dim(Bs)[2]

# Define spatial domain (1, 13)
t_axis <- seq(0, 13, length.out = 57); xt_ = seq(0, 13, length.out = 57)

#Scale the temporal 'coordinates'
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

#Find distance 
dist <- (max(xt) - min(xt)) / n_knots_t
xtl <- min(xt) - dist * 0.05; xtr <- max(xt) + dist * 0.05
dxt <- (xtr - xtl) / n_knots_t

#Define knot locations in time
knots <- seq(xtl - bdeg * dxt, xtr + bdeg * dxt, by = dxt)

#Define temporal marginal B-spline basis design matrix
Bt <- spline.des(knots, xt, bdeg + 1)$design
kt <- ncol(Bt)

#Define the 3-dimensional P-spline basis function
Bst <- kronecker(Bt, Bs)

#Define the marginal penalization matrices (RW1 precision matrices)
Px <- INLA:::inla.rw(n = kx, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Py <- INLA:::inla.rw(n = ky, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Pt <- INLA:::inla.rw(n = kt, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)


#Identity matrices k1 x k1, k2 x k2, and kt x kt
Ix <- diag(kx); Iy <- diag(ky); It <- diag(kt)

#Combined spatial penalty matrix

tau_s = 25 # Isotropic spatial smoothness param
tau_t = 150 # Temporal smoothness parameter

Pst <- tau_s *(It %x% Py %x% Ix + It %x% Iy %x% Px) + 
        tau_t * (Pt %x% Iy %x% Ix)


#Dont know if needed, but I'll keep it
Pst <- Pst + diag(x = 0.0001, nrow = ncol(Pst))

#Set a random choosen seed

n <- 1 # number of simulations
intercept = -7.6 #Intercept


#Sampled space-time P-spline parameters
set.seed(07101984)
sampled_theta_st = inla.qsample(n = 1, Q = Pst)

#Sum-to-zero (since IGMRF approx)
sampled_theta_st = sampled_theta_st - mean(sampled_theta_st)
rownames(sampled_theta_st) <- 1:nrow(sampled_theta_st)

#Get the space-time interaction
interactions = Bst %*% sampled_theta_st[, 1]

#Make a linearly increasing temporal effect
beta_t = 0.0142
tmp_ = xt_ * beta_t; temporal_effect = rep(0, length(interactions))
for(t in 1:length(xt)){
  temporal_effect[((t - 1) * (dim(Bs)[1]) + 1):(t * dim(Bs)[1])] = tmp_[t]
}

#Get the risk-field
Lambda_st <- exp(as.vector(intercept + temporal_effect + interactions))

## make a yxt-grid
yxt_grid = expand.grid(y = y_axis, x = x_axis, t = t_axis)

# add the values to yxt_grid
yxt_grid$values = Lambda_st

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y 

##Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                                        st_make_valid(germany_border))

#Reset the indices (Necessary!)
indices_risk_surface.list = nrow(risk_surface.list)
rownames(risk_surface.list) = 1:indices_risk_surface.list

#Make grid over Germany (both polygons and centroid of the polygons)
tactical_error: Below changes a lot!
polygon_grid = st_as_sf(st_make_grid(germany_border, n = c(75, 75), square = FALSE),
                        crs = crs_$srid)
Maybe a better strategy now is to take the centroids of each polygon,
take locations, find values, add to df with polygons, and then I can plot w. ggplot

centroid_grid = st_make_grid(germany_border, n = c(75, 75), 
                             square = F, what = "centers")



plot(st_geometry(germany_border))
plot(polygon_grid, add = T)







#Code to make a GIF

for(time in unique(risk_surface.list$t)){
  if(time == risk_surface.list$t[1]){files = c(); i = 1}
  
  #Extract the risk-surface at time: t=time
  tmp_ = risk_surface.list[risk_surface.list$t == time, c("values", "geometry")]
  
  #Create a filename (and save filename) for where current time risk-surface is stored
  filename = paste("./Plots/P_spline_tests/", toString(i), ".png", sep = "")
  print(filename)
  png(filename = filename); files[i] = filename
  
  #Visualize the sampled values on Germany
  plot(st_geometry(germany_border), border = "grey", title = time) + title(main = time)
  plot(tmp_,
       add = TRUE)
  
  #Something needed to reset saving of image AND update counter i
  dev.off(); i = i + 1
}


### Make a GIF (for fun)
img = c(image_read(files[1]))
for(name in files[2:length(files)]){
  img = c(img, image_read(name))
}

image_append(image_scale(img, "x200"))

my.animation <- image_animate(image_scale(img, 
                                          "400x400"),
                              fps = 10,
                              dispose = "previous")

image_write(my.animation, "test.gif")

#State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)

# Make a identifier for area and time that can be used to map locations and times
# to correct area and time. Unique id for each area, integer time, and combined

# ids runs over areas then time i.e. area 1 time 1, area 2 time 1,..., area n time 1, area 1 time 2,...
ids <- 1:(nrow(germany_map_2) * tT)

#indices will be used to map risk_surface.list$values to correct area and time
indices = 1:nrow(risk_surface.list)

mapping.df <- data.frame(indices = indices, ids = rep(0, length(indices)))

risk_surface.list$area_id = rep(0, nrow(risk_surface.list))
risk_surface.list$time_id = rep(0, nrow(risk_surface.list))

# First add time id based on time integer value
risk_surface.list$time_id = ceiling(risk_surface.list$t) 
mapping.df$time_id = risk_surface.list$time_id

#Add area id based on which area the location is within
area_id = rep(0, nrow(risk_surface.list))
for(i in 1:nrow(germany_map_2)){
  print(i)
  tmp_ = sf::st_intersection(risk_surface.list, 
                             st_make_valid(germany_map_2$geometry[i]))
  
  indices = as.numeric(rownames(tmp_))
  area_id[indices] = germany_map_2$ID_1[i]
}
risk_surface.list$area_id = area_id
mapping.df$area_id = risk_surface.list$area_id


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
    
    tmp_ = risk_surface.list[risk_surface.list$area_id == i &
                               risk_surface.list$time_id == t, ]
    lambda.df[index, ]$lambda_it = mean(tmp_$values)
  }
}

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})

RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)

# Make precision matrix for Besag
matrix4inla <- nb2mat(nb2, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

#Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

#Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

#Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))

################################################################################
#Fit a basic model

#Specify the base formula
base_formula <- sampled_counts ~ 1 + f(time_id, 
                                       model = 'bym2',
                                       scale.model = T, 
                                       constr = T, 
                                       rankdef = 1,
                                       graph = RW1_prec,
                                       hyper = temporal_hyper) + 
  f(area_id, 
    model = 'bym2',
    scale.model = T,
    constr = T,
    rankdef = 1,
    graph = Besag_prec,
    hyper = spatial_hyper)


ptm <- Sys.time()
improper_noInt <- inla(base_formula, data = lambda.df, family = "poisson",
                       E = E_it, 
                       control.compute = list(config = TRUE, # To see constraints later
                                              cpo = T,   # For model selection
                                              waic = T)) # For model selection

time_improper_noInt = Sys.time()-ptm
print(c("Basic model fitted in: ", time_improper_noInt))

plot(improper_noInt)


################################################################################

#Get constraints for type IV interactions
typeIV_constraints <- constraints_maker(type = "IV", n = nrow(germany_map_2),
                                        t = 13)


# get scaled ICAR
scaled_ICAR_prec <- INLA::inla.scale.model(Besag_prec, 
                                           constr = list(A = matrix(1,1,dim(Besag_prec)[1]), e = 0))


#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))


#Get type IV interaction precision matrix
typeIV_prec <- scaled_RW_prec %x% scaled_ICAR_prec

#Get formula for type IV
typeIV_formula <- update(base_formula, ~. + f(space.time, 
                                              model = "generic0",
                                              Cmatrix = typeIV_prec,
                                              extraconstr = typeIV_constraints,
                                              rankdef = (nrow(germany_map_2) + 13 - 1), 
                                              hyper = interaction_hyper))



lambda.df$space.time = 1:(nrow(germany_map_2) * 13)

ptm <- Sys.time()
RW1_ICAR_IV_fit <- inla(typeIV_formula, data = lambda.df, family = "poisson",
                        E = E_it, control.compute = list(config = TRUE, 
                                                                cpo = TRUE,
                                                                waic = TRUE))
time_RW1_ICAR_IV = Sys.time() - ptm
print(c("Type IV model fitted in: ", time_RW1_ICAR_IV))

plot(RW1_ICAR_IV_fit)

















