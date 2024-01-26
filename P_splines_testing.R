#Clear environment
rm(list = ls())

#Load libraries
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)


#read the shapefile of Germany
germany_map  <- read_sf('./Shapefiles/Germany')[, c("ID_2", "NAME_2", "geometry")]

#Make sure no NA's in the map!
germany_map <- na.omit(germany_map)

#Plot dependency structure of Germany
nb <- spdep::poly2nb(germany_map, queen = FALSE)
plot(st_geometry(germany_map), border = "grey")
plot.nb(nb, st_geometry(germany_map), add = TRUE)

#Merge all the polygons into one polygon, internal borders still somewhat present
sfu = sf::st_union(germany_map$geometry) 
plot(st_geometry(sfu), border = "grey") #Plot the polygon

#Extract the coordinates of the large polygon
gpts = sf::st_coordinates(sfu)
gmat <- matrix(c(gpts[,1], gpts[,2]), ncol = 2)

#Find the outer border (Chull) of the polygon of the entirity of Germany
cv.id <- chull(gmat) #Contains indices making the border
cv.id <- c(cv.id, cv.id[1])

#Transform border back to polygon and plot
border_coords = data.frame(coords = gmat[cv.id, ])
border_polygon <- sf::st_as_sf(border_coords, coords = c("coords.1", "coords.2"))
plot(st_geometry(border_polygon), border = "grey")
plot(st_geometry(germany_map), border = "grey", add = T)
