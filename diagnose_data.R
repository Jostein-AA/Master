
# Look at how the temporal trend of the data is across all regions
# i.e. plot all together or average over regions for one year
# and see if it has the temporal trend it is supposed to have


# Other diagnoses of the data

#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

# Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADM4 <- nrow(second_level_admin_map)

################################################################################




