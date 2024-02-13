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

#Load in data
source('Preliminaries.R')
source("Utilities.R")

################################################################################
#Plot the structures

#Plot geographical polygons and dependency structure of Germany
par(mfrow = c(1, 2))
nb <- spdep::poly2nb(germany_map, queen = FALSE)
plot(st_geometry(germany_map), border = "grey")
plot.nb(nb, st_geometry(germany_map), add = TRUE)

nb2 <- spdep::poly2nb(st_make_valid(germany_map_2), queen = FALSE)
plot(st_geometry(germany_map_2), border = "grey")
plot.nb(nb2, st_geometry(germany_map_2), add = TRUE)
par(mfrow = c(1,1))

crs_ = st_crs(germany_border, parameters = TRUE)

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

#Extract a spatial domain [x-min, x-max] X [y-min, y-max]
gpts = sf::st_coordinates(germany_border$geometry)

## Find the x-min an x-max
x_min = min(gpts[,c("X")]); x_max = max(gpts[,c("X")])

## Find y-min and y-max
y_min = min(gpts[,c("Y")]); y_max = max(gpts[,c("Y")])

#Find equally spaced knots (20 each)
n_knots = 15

## 20 equi-spaced knots in X
x_knots = seq(x_min, x_max, length.out = n_knots)[2:(n_knots - 1)]
x_axis = seq(x_min, x_max, length.out = 100)

## 20 equi-spaced knots in Y
y_knots = seq(y_min, y_max, length.out = n_knots)[2:(n_knots - 1)]
y_axis = seq(y_min, y_max, length.out = 100)

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
x_penalty = 0.1; y_penalty = 0.1

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

#Produce a data frame with the appropriate points and values
#x_axis, y_axis
#risk_surface.list = list()

## make a xy-grid
xy_grid = expand.grid(y = y_axis, x = x_axis)

## add risk values to grid
xy_grid$values = risk_surface[, 1]

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(xy_grid, coords = c("x", "y"), crs = crs_$srid)

##Drop points outside Germany
risk_surface.list = sf::st_intersection(risk_surface.list, 
                  st_make_valid(germany_border))


#Visualize
plot(st_geometry(germany_border), border = "grey")
plot(risk_surface.list,
     add = TRUE)
















