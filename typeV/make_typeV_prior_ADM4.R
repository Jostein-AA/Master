#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
#Load libraries
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)
library(ggpubr)
library(pspline)
library(readxl)
library(OneR)
library(mgcv)
library(splines)
library(ggfortify)
library(crs)
#library(magick)
library(INLA)
#library(bigDM)
library(MASS)
#library(RColorBrewer)
#library(tmap)


library("makemyprior")

source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
## Specify precision matrices
#---
### Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)

#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))

### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_second_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_second_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

# get scaled Besag
scaled_besag_prec_second_level <- INLA::inla.scale.model(Besag_prec_second_level, 
                                                         constr = list(A = matrix(1,1,dim(Besag_prec_second_level)[1]), 
                                                                       e = 0))

INLA::inla.write.graph(
  INLA::inla.matrix2graph(RW1_prec),
  filename = "RW1_prec.graph"
)
INLA::inla.write.graph(
  INLA::inla.matrix2graph(Besag_prec_second_level),
  filename = "Besag_prec_second_level.graph"
)


RW1_prec_path <- paste0(getwd(), "/RW1_prec.graph")
Besag_prec_second_level_path <- paste0(getwd(), "/Besag_prec_second_level.graph")


#Get precision matric for type II interaction by Kronecker product
typeII_prec_second_level <- scaled_RW_prec %x% diag(nrow(second_level_admin_map))

#Get precision matric for type II interaction by Kronecker product
typeIII_prec_second_level <- diag(tT) %x% scaled_besag_prec_second_level

#Get type IV interaction precision matrix
typeIV_prec_second_level <- scaled_RW_prec %x% scaled_besag_prec_second_level

#---
# Get the sum-to-zero constraints

#Get sum-to-zero constraints for type II interaction
typeII_constraints_second_level = constraints_maker(type = "II", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)

#Get sum-to-zero constraints for type II interaction
typeIII_constraints_second_level = constraints_maker(type = "III", 
                                                     n = nrow(second_level_admin_map), 
                                                     t = tT)

#Get sum-to-zero constraints for type IV interaction
typeIV_constraints_second_level = constraints_maker(type = "IV", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)

#---


################################################################################
# Create formulas

# Load in a singular data set
load("Simulated_data/sc2/sc2_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, 1]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_II = lambda_$space.time
lambda_$space.time_III = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time


# Make a prior
source("typeV/prior_wrapper.R")

saveRDS(prior_data_bym2_time, "./typeV/prior_data_bym2_time_ADM4.rds")
saveRDS(prior_data_bym2_space, "./typeV/prior_data_bym2_space_ADM4.rds")
saveRDS(prior_data_interaction, "./typeV/prior_data_interaction_ADM4.rds")



# typeV_formula <- sampled_counts ~ 1 + f(time_id_iid,      # Unstructured temporal effect
#                                         model = "iid") + 
#                                       f(time_id_struct,   # Structured temporal effect
#                                         model = "besag", 
#                                         graph = RW1_prec_path,
#                                         scale.model = TRUE,
#                                         constr = T) + 
#                                       f(area_id_iid,      # Unstructured spatial effect
#                                         model = "iid") + 
#                                       f(area_id_struct,   # Structured spatial effect
#                                         model = "besag",
#                                         graph = Besag_prec_second_level_path,
#                                         scale.model = T,
#                                         constr = T) + 
#                                       f(space.time_I,     # Unstructured space-time interaction
#                                         model = "iid") + 
#                                       f(space.time_II,    # Structured time interacting w. iid space
#                                         model = "generic0",
#                                         Cmatrix = typeII_prec_second_level,
#                                         extraconstr = typeII_constraints_second_level,
#                                         rankdef = nrow(second_level_admin_map)) +
#                                       f(space.time_III,   # iid time interacting w. structured space
#                                         model = "generic0",
#                                         Cmatrix = typeIII_prec_second_level,
#                                         extraconstr = typeIII_constraints_second_level,
#                                         rankdef = tT) +
#                                       f(space.time_IV,    # Structured time interacting w. structured space
#                                         model = "generic0",
#                                         Cmatrix = typeIV_prec_second_level,
#                                         extraconstr = typeIV_constraints_second_level)
# 
# 
# 
# 
# # Make the prior
# raw_prior_typeV = list(
#   tree = 's1 = (time_id_iid, time_id_struct); s2 = (area_id_iid, area_id_struct); s3 = (space.time_I, space.time_II, space.time_III, space.time_IV)',
#   V = list(s1 = list(prior = "pc",
#                      param = c(1, 0.01)),
#            s2 = list(prior = "pc",
#                      param = c(1, 0.01)),
#            s3 = list(prior = "pc",
#                      param = c(1, 0.01))),
#   w = list(s1 = list(prior = "pc1",
#                      param = 0.51),
#            s2 = list(prior = "pc1",
#                      param = 0.51),
#            s3 = list(prior = "dirichlet"))
# )
# 
# 
# 
# 
# # Load in a singular data set
# load("Simulated_data/sc2/sc2_data.RData")
# 
# ### Format data for makemyprior
# lambda_ = lambda.df[, c("area_id", "time_id",
#                         "E_it", "space.time")]
# 
# 
# lambda_$sampled_counts = lambda.df$sampled_counts[, 1]
# lambda_$area_id_iid = lambda_$area_id
# lambda_$area_id_struct = lambda_$area_id
# lambda_$time_id_iid = lambda_$time_id
# lambda_$time_id_struct = lambda_$time_id
# lambda_$space.time_I = lambda_$space.time
# lambda_$space.time_II = lambda_$space.time
# lambda_$space.time_III = lambda_$space.time
# lambda_$space.time_IV = lambda_$space.time
# 
# 
# print("Starting to make the prior")
# 
# 
# typeV_prior_ADM4 <- make_prior(typeV_formula, lambda_,
#                           family = "poisson",
#                           prior = raw_prior_typeV)
# 
# ### Inspecting the prior
# #plot_prior(typeV_prior_ADM4)
# #plot_tree_structure(typeV_prior_ADM4)
# 
# 
# print("Prior made")
# 
# saveRDS(typeV_prior_ADM4, "./typeV/typeV_prior_ADM4.rds")
# 
# print("Prior saved")




