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
library(GpGp)

#Load in data
source('Preliminaries.R')
source("Utilities.R")


#Load in structures
nb <- spdep::poly2nb(germany_map, queen = FALSE)
nb2 <- spdep::poly2nb(st_make_valid(germany_map_2), queen = FALSE)
#Get the coordinate system used
crs_ = st_crs(germany_border, parameters = TRUE)

#Extract a spatial domain to find [x-min, x-max] X [y-min, y-max]
gpts = sf::st_coordinates(germany_border$geometry)

## 1) Find the x-min an x-max 2) Find y-min and y-max
x_min = min(gpts[,c("X")]); x_max = max(gpts[,c("X")])
y_min = min(gpts[,c("Y")]); y_max = max(gpts[,c("Y")])


## Define the spatial grid xy_grid and standardize the coordinates (important!)
x_axis = seq(x_min, x_max, length.out = 50) #50 points to evaluate on x-axis
y_axis = seq(y_min, y_max, length.out = 50) #50 points to evaluate on y-axis
#xy-grid containing 1000 points to evaluate the spatial domain on
xy_grid = expand.grid(x = x_axis, y = y_axis)

#standardize x and y
locs <- apply(xy_grid, 2, function(x) (x - min(x)) / (max(x) - min(x)))

#Extract coordinates lat/lon
x <- locs[, 1]
y <- locs[, 2]

#Find the t_axis
t_axis = seq(0, 13, length.out = 14)
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

#Make a matrix with x,y,t
xyt_matrix = matrix(data = rep(0, 3 * length(x)), 
                    ncol = 3, nrow = length(x))

xyt_matrix[, 1] = x; xyt_matrix[, 2] = y; xyt_matrix[, 3] = rep(xt[1], nrow(xyt_matrix))

for(t in 2:length(xt)){
  tmp_ = matrix(data = rep(0, 3 * length(x)), 
                ncol = 3, nrow = length(x))
  
  tmp_[, 1] = x; tmp_[, 2] = y; tmp_[, 3] = rep(xt[t], nrow(tmp_))
  xyt_matrix = rbind(xyt_matrix, tmp_)
}

sigma2 = 0.2
spatial_range = 0.13
temporal_range = 0.37
smoothness = 2.1
nugget = 0.001

covparams = c(sigma2, spatial_range, temporal_range, smoothness, nugget)

#A vector with covariance parameters in the form (, , , , nugget)

#Holy Fuck that is expensive
test = matern_spacetime(covparams, xyt_matrix)











