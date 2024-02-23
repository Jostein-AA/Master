#Clear environment
rm(list = ls())

#Library needed to make GIFs
library(magick)

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

t_axis <- seq(0.0001, 13, length.out = 40)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc1_risk_surface_1.RData")
dir = "./Plots/scenario_1/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc1.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc3_risk_surface_1.RData")
dir = "./Plots/scenario_3/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)

gif_name = "sc3.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc5_risk_surface_1.RData")
dir = "./Plots/scenario_5/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc5.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc7_risk_surface_1.RData")
dir = "./Plots/scenario_7/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc7.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc9_risk_surface_1.RData")
dir = "./Plots/scenario_9/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc9.gif"
make_GIF(dir, gif_name)


## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc11_risk_surface_1.RData")
dir = "./Plots/scenario_11/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc11.gif"
make_GIF(dir, gif_name)



