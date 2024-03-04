#Clear environment
rm(list = ls())

#Library needed to make GIFs
library(magick)

source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")


## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc1_risk_surface_2.RData")
dir = "./Plots/scenario_1_2/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc1_2.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc3_risk_surface_2.RData")
dir = "./Plots/scenario_3_2/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)

gif_name = "sc3_2.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc5_risk_surface_2.RData")
dir = "./Plots/scenario_5_2/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc5_2.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc7_risk_surface_2.RData")
dir = "./Plots/scenario_7_2/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc7_2.gif"
make_GIF(dir, gif_name)

## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc9_risk_surface_2.RData")
dir = "./Plots/scenario_9_2/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc9_2.gif"
make_GIF(dir, gif_name)


## Save plots for validation of risk-surface
load("./Data/Simulated_risk_surfaces/sc11_risk_surface_2.RData")
dir = "./Plots/scenario_11_2/"
plots_for_GIF(risk_surface.list = risk_surface.list,
              polygons = polygon_grid2,
              t_axis = t_axis,
              filename_base = dir)


gif_name = "sc11_2.gif"
make_GIF(dir, gif_name)



