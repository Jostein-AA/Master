#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("makemyprior")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

load("typeV/one_dataset_fitted_typeV.RData")

plot(res)

res$marginals.fitted.values

res$marginals.random

res$marginals.hyperpar
