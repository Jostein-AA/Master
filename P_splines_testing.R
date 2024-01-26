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

#Plot dependency structure of Germany
nb <- spdep::poly2nb(germany_map, queen = FALSE)
plot(st_geometry(germany_map), border = "grey")
plot.nb(nb, st_geometry(germany_map), add = TRUE)



