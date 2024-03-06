#Clear environment
rm(list = ls())

#Library needed to make GIFs
library(magick)

source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")


## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc1_risk_surfaces.RData")
risk_surface.list$values <- risk_surface.list$values[, 1]

dir = "./Plots/scenario_1/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc1.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc3_risk_surfaces.RData")
risk_surface.list$values <- risk_surface.list$values[, 1]

dir = "./Plots/scenario_3/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)

gif_name = "sc3.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc5_risk_surfaces.RData")
risk_surface.list$values <- risk_surface.list$values[, 1]

dir = "./Plots/scenario_5/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc5.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc7_risk_surfaces.RData")
risk_surface.list$values <- risk_surface.list$values[, 1]

dir = "./Plots/scenario_7/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc7.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc9_risk_surfaces.RData")
risk_surface.list$values <- risk_surface.list$values[, 1]

dir = "./Plots/scenario_9/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc9.gif"
make_GIF(dir, gif_name)


## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc11_risk_surfaces.RData")
risk_surface.list$values <- risk_surface.list$values[, 1]

dir = "./Plots/scenario_11/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc11.gif"
make_GIF(dir, gif_name)



