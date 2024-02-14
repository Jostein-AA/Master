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

#Load in data
source('Preliminaries.R')
source("Utilities.R")

nb <- spdep::poly2nb(germany_map, queen = FALSE)
nb2 <- spdep::poly2nb(st_make_valid(germany_map_2), queen = FALSE)
crs_ = st_crs(germany_border, parameters = TRUE)

#Preliminaries for tensor product smooth using 3 P-splines

#Extract a spatial domain [x-min, x-max] X [y-min, y-max]
gpts = sf::st_coordinates(germany_border$geometry)

## Find the x-min an x-max
x_min = min(gpts[,c("X")]); x_max = max(gpts[,c("X")])

## Find y-min and y-max
y_min = min(gpts[,c("Y")]); y_max = max(gpts[,c("Y")])

# Define number of knots in spatial domain 
n_knots = 12

#Define number of knots in temporal domain
n_knots_t = 6

## n_knots equi-spaced knots in X
x_knots = seq(x_min, x_max, length.out = n_knots)[2:(n_knots - 1)]
x_axis = seq(x_min, x_max, length.out = 80)

## n_knots equi-spaced knots in Y
y_knots = seq(y_min, y_max, length.out = n_knots)[2:(n_knots - 1)]
y_axis = seq(y_min, y_max, length.out = 80)

## n_knots_t equi-spaced knots in time
t_knots = seq(0, 13, length.out = n_knots_t)[2:(n_knots_t - 1)]
t_axis = seq(-0.1, 13.1, length.out = 67)

# Make B-spline basis functions 

## B-spline basis functions in X
x_basis <- bs(x_axis, knots = x_knots, 
              Boundary.knots = c(x_min, x_max), degree = 3)
x_basis <- x_basis[, 1:(ncol(x_basis) - 1)]
plot(x_basis[, 1] ~ x_axis, type = "l", ylim = c(0, 1.5))
for(i in 2:ncol(x_basis)){lines(x_basis[, i] ~ x_axis)}

## B-spline basis functions in Y
y_basis <- bs(y_axis, knots = y_knots, 
              Boundary.knots = c(y_min, y_max), degree = 3)
y_basis <- y_basis[, 1:(ncol(y_basis) - 1)]
plot(y_basis[, 1] ~ y_axis, type = "l", ylim = c(0, 1.5))
for(i in 2:ncol(y_basis)){lines(y_basis[, i] ~ y_axis)}

t_basis <- bs(t_axis, knots = t_knots, 
              Boundary.knots = c(min(t_axis), max(t_axis)), degree = 3)
t_basis <- t_basis[2:(nrow(t_basis) - 1), 1:(ncol(t_basis) - 1)]
t_axis = t_axis[2:(length(t_axis) - 1)]
plot(t_basis[, 1] ~ t_axis, type = "l", ylim = c(0, 1.5))
for(i in 2:ncol(t_basis)){lines(t_basis[, i] ~ t_axis)}

################################################################################
#Plot the structures

#Plot geographical polygons and dependency structure of Germany
par(mfrow = c(1, 2))
plot(st_geometry(germany_map), border = "grey")
plot.nb(nb, st_geometry(germany_map), add = TRUE)

plot(st_geometry(germany_map_2), border = "grey")
plot.nb(nb2, st_geometry(germany_map_2), add = TRUE)
par(mfrow = c(1,1))

#Plot heatmap of population of Germany
scale_col = heat.colors(50, rev=TRUE) #Divide color gradient into 50 
scale = scale_col[c(3, 7, 15, 20, 25, 30, 35, 37, 40, 43, 45, 47, 49, 50)] #Select color scale to be more red
hardcoded_bins =  c(23800, 50000, 100000, 300000, #Hardcode bins for population 
                    600000, 1000000, 1400000, 1800000, 2100000, #Heat grade
                    2400000, 2700000, 3000000, 3300000, 3600000)
#OneR::bin(germany_map$Population, nbins = 8, method = "length") #Function to perform binning, no work

heatmap_areas(germany_map, value = germany_map$Population, scale_col, #Function that plots
              scale, hardcoded_bins, 'Population Germany') #Heatmap of map and value provided

################################################################################
#Make a two-dimensional P-spline over Germany!

### Make the combined B-matrix: B_x kronecker B_y
B_matrix <- x_basis %x% y_basis

#Draw parameters for 2-dim P-spline from prior distribution

## Find the penalization matrix
x_penalization_matrix <- INLA:::inla.rw(n = ncol(x_basis), order = 1, 
                                        scale.model = FALSE, 
                                        sparse = TRUE)
x_identity_matrix <- diag(x = 1, nrow = ncol(x_basis), ncol = ncol(x_basis))

y_penalization_matrix <- INLA:::inla.rw(n = ncol(y_basis), order = 1, 
                                        scale.model = FALSE, 
                                        sparse = TRUE)

y_identity_matrix <- diag(x = 1, nrow = ncol(y_basis), ncol = ncol(y_basis))

### The combined penalization matrix
x_penalty = 4; y_penalty = 4

S = x_penalty * (x_penalization_matrix %x% y_identity_matrix) + 
    y_penalty * (x_identity_matrix %x% y_penalization_matrix)

#Just for testing
S = S + diag(x = 0.01, nrow = nrow(S), ncol = ncol(S))

## Draw from normal distribution

#inla.qsample ???
set.seed(1234)
sampled_parameters <- inla.qsample(n = 1,
                                   Q = S)
rownames(sampled_parameters) <- 1:nrow(sampled_parameters)

#Calculate resulting field
risk_surface <- B_matrix %*% sampled_parameters

#Produce a list with the appropriate points and values

## make a xy-grid
xy_grid = expand.grid(y = y_axis, x = x_axis)

## add risk values to grid
xy_grid$values = risk_surface[, 1]

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(xy_grid, coords = c("x", "y"), crs = crs_$srid)

##Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                  st_make_valid(germany_border))


#Visualize the sampled values on Germany
plot(st_geometry(germany_border), border = "grey")
plot(risk_surface.list,
     add = TRUE)


################################################################################
#Make a 3-dim P-spline over Germany

### Make the combined B-matrix: B_t Kronecker B_x Kronecker B_y  
B_matrix <- t_basis %x% x_basis %x% y_basis 

#Draw parameters for 3-dim P-spline from prior distribution

## Find the penalization matrix
x_penalization_matrix <- INLA:::inla.rw(n = ncol(x_basis), order = 1, 
                                        scale.model = FALSE, 
                                        sparse = TRUE)
x_identity_matrix <- diag(x = 1, nrow = ncol(x_basis), ncol = ncol(x_basis))

y_penalization_matrix <- INLA:::inla.rw(n = ncol(y_basis), order = 1, 
                                        scale.model = FALSE, 
                                        sparse = TRUE)

y_identity_matrix <- diag(x = 1, nrow = ncol(y_basis), ncol = ncol(y_basis))

t_penalization_matrix <- INLA:::inla.rw(n = ncol(t_basis), order = 1, 
                                        scale.model = FALSE, 
                                        sparse = TRUE)

t_identity_matrix <- diag(x = 1, nrow = ncol(t_basis), ncol = ncol(t_basis))


### The combined penalization matrix
x_penalty = 10; y_penalty = 10; t_penalty = 200

S = x_penalty * (t_identity_matrix %x% x_penalization_matrix %x% y_identity_matrix) + 
  y_penalty * (t_identity_matrix %x% x_identity_matrix %x% y_penalization_matrix) + 
  t_penalty * (t_penalization_matrix%x% x_identity_matrix %x% y_identity_matrix)

#Just for testing
S = S + diag(x = 0.01, nrow = nrow(S), ncol = ncol(S))

## Draw from normal distribution

set.seed(6547321)
sampled_parameters <- inla.qsample(n = 1,
                                   Q = S, 
                                   mu = rep(0.3, nrow(S)))
rownames(sampled_parameters) <- 1:nrow(sampled_parameters)

#Calculate resulting field
risk_surface <- B_matrix %*% sampled_parameters

#Remove large B_matrix NB!!!!!!!!!!! Only do this if no more use of B_matrix!!!
#rm(B_matrix)

## make a yxt-grid
yxt_grid = expand.grid(y = y_axis, x = x_axis, t = t_axis)

# add the values to yxt_grid
yxt_grid$values = risk_surface[ ,1]

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

################################################################################
#Integrate to get discrete risk values for areas and times

## Is it just to sum all the points within each area at that time and divide by
## number of points within???

#rm(B_matrix) # Remove B_matrix as it is memory intensive on mah poor computah


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
                        E_it = rep(10, nrow(germany_map_2) * tT))


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

################################################################################
#Simulate values from Poisson with expected value = E_it * lambda_it as in lambda.df

lambda.df$mu = lambda.df$E_it * lambda.df$lambda_it

lambda.df$sampled_counts = apply(lambda.df, 
                                 MARGIN = 1, 
                                 FUN = function(row){return(rpois(1, row[5]))})



################################################################################
# Apply models to the sampled data set!!!













