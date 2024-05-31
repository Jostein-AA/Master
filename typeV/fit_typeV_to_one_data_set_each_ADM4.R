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

dataset_id_2 = 4

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



#Get type IV interaction precision matrix
typeIV_prec_second_level <- scaled_RW_prec %x% scaled_besag_prec_second_level

#---
# Get the sum-to-zero constraints

#Get sum-to-zero constraints for type IV interaction
typeIV_constraints_second_level = constraints_maker(type = "IV", 
                                                    n = nrow(second_level_admin_map), 
                                                    t = tT)

#---
print("Precision matrices made")
################################################################################

### Load in a singular data set
# Load in a singular data set
load("Simulated_data/sc2/sc2_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, dataset_id_2]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time

prior_data_bym2_time = readRDS("./typeV/prior_data_bym2_time_ADM4.rds")
prior_data_bym2_space = readRDS("./typeV/prior_data_bym2_space_ADM4.rds")
prior_data_interaction = readRDS("./typeV/prior_data_interaction_reduced_ADM4.rds")

# obs! merk at rekkefoelgen paa inla-formelen maa matche rekkefoelgen vi bruker log-presisjonene her!
prior_func <- function(logprec){
  return(
    eval_joint_prior(-c(logprec[1:2]), prior_data_bym2_time) +
      eval_joint_prior(-c(logprec[3:4]), prior_data_bym2_space) +
      eval_joint_prior(-c(logprec[5:8]), prior_data_interaction)
  )
}


typeV_formula <- sampled_counts ~ 1 + f(time_id_iid,      # Unstructured temporal effect
                                        model = "iid") + 
  f(time_id_struct,   # Structured temporal effect
    model = "besag", 
    graph = RW1_prec_path,
    scale.model = TRUE,
    constr = T) + 
  f(area_id_iid,      # Unstructured spatial effect
    model = "iid") + 
  f(area_id_struct,   # Structured spatial effect
    model = "besag",
    graph = Besag_prec_second_level_path,
    scale.model = T,
    constr = T) + 
  f(space.time_I,     # Unstructured space-time interaction
    model = "iid") + 
  f(space.time_IV,    # Structured time interacting w. structured space
    model = "generic0",
    Cmatrix = typeIV_prec_second_level,
    extraconstr = typeIV_constraints_second_level)


# Set all to NA to see what is the prior
#lambda_$sampled_counts <- NA
#res_prior <- inla(typeV_formula, 
#                  data = lambda_, 
#                  family = "poisson",
#                  verbose = T,
#                  control.fixed= list(prec.intercept = 1),
#                  control.expert = list(jp = inla.jp.define(
#                    prior_func,
#                    prior_data_bym2_time = prior_data_bym2_time,
#                    prior_data_bym2_space = prior_data_bym2_space,
#                    prior_data_interaction = prior_data_interaction,
#                    eval_joint_prior = makemyprior::eval_joint_prior,
#                    hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
#                    calc_jac_logdet = makemyprior:::calc_jac_logdet,
#                    choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
#                    cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
#                    expit = makemyprior:::expit,
#                    get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
#                    get_indexes = makemyprior:::get_indexes,
#                    get_indexes2 = makemyprior:::get_indexes2,
#                    hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
#                    hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
#                    eval_spline_lpdf = makemyprior:::eval_spline_lpdf
#                  )))

#samps <- inla.hyperpar.sample(n = 1000, res_prior) # husk at dette naa blir paa presisjons-skala, du maa regne om til
# varianser og vekter selv!! siden vi har samples kan du transformere direkte og saa plotte sammen med prior fra makemyprior

# hente ut ferdig evaluerte marginale PC prior-fordelinger for vekter fra makemyprior:
#eval_pc_prior(x = seq(0, 1, 0.01), obj = prior_space, param = "w[area_id_iid/area_id_iid_area_id_struct]")



lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA

#lc <- inla.make.lincombs(space.time_I = rep(1, 13), space.time_IV = rep(1, 13))

res <- inla(typeV_formula, 
            data = lambda_, 
            verbose = T,
            family = "poisson",
            E = lambda_$E_it,
            #lincomb = lc,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE),
            control.expert = list(jp = inla.jp.define(
              prior_func,
              prior_data_bym2_time = prior_data_bym2_time,
              prior_data_bym2_space = prior_data_bym2_space,
              prior_data_interaction = prior_data_interaction,
              eval_joint_prior = makemyprior::eval_joint_prior,
              hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
              calc_jac_logdet = makemyprior:::calc_jac_logdet,
              choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
              cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
              expit = makemyprior:::expit,
              get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
              get_indexes = makemyprior:::get_indexes,
              get_indexes2 = makemyprior:::get_indexes2,
              hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
              hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
              eval_spline_lpdf = makemyprior:::eval_spline_lpdf
            )))

#plot(res)

save(res,
     file = "./typeV/one_dataset_fitted_typeV_sc2.RData")


# Load in a singular data set
load("Simulated_data/sc4/sc4_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, dataset_id_2]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time

lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA

res <- inla(typeV_formula, 
            data = lambda_, 
            verbose = T,
            family = "poisson",
            E = lambda_$E_it,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE),
            control.expert = list(jp = inla.jp.define(
              prior_func,
              prior_data_bym2_time = prior_data_bym2_time,
              prior_data_bym2_space = prior_data_bym2_space,
              prior_data_interaction = prior_data_interaction,
              eval_joint_prior = makemyprior::eval_joint_prior,
              hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
              calc_jac_logdet = makemyprior:::calc_jac_logdet,
              choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
              cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
              expit = makemyprior:::expit,
              get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
              get_indexes = makemyprior:::get_indexes,
              get_indexes2 = makemyprior:::get_indexes2,
              hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
              hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
              eval_spline_lpdf = makemyprior:::eval_spline_lpdf
            )))

#plot(res)

save(res,
     file = "./typeV/one_dataset_fitted_typeV_sc4.RData")




# Load in a singular data set
load("Simulated_data/sc6/sc6_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, dataset_id_2]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time

lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA

res <- inla(typeV_formula, 
            data = lambda_, 
            verbose = T,
            family = "poisson",
            E = lambda_$E_it,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE),
            control.expert = list(jp = inla.jp.define(
              prior_func,
              prior_data_bym2_time = prior_data_bym2_time,
              prior_data_bym2_space = prior_data_bym2_space,
              prior_data_interaction = prior_data_interaction,
              eval_joint_prior = makemyprior::eval_joint_prior,
              hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
              calc_jac_logdet = makemyprior:::calc_jac_logdet,
              choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
              cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
              expit = makemyprior:::expit,
              get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
              get_indexes = makemyprior:::get_indexes,
              get_indexes2 = makemyprior:::get_indexes2,
              hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
              hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
              eval_spline_lpdf = makemyprior:::eval_spline_lpdf
            )))

#plot(res)

save(res,
     file = "./typeV/one_dataset_fitted_typeV_sc6.RData")



# Load in a singular data set
load("Simulated_data/sc8/sc8_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, dataset_id_2]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time

lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA

res <- inla(typeV_formula, 
            data = lambda_, 
            verbose = T,
            family = "poisson",
            E = lambda_$E_it,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE),
            control.expert = list(jp = inla.jp.define(
              prior_func,
              prior_data_bym2_time = prior_data_bym2_time,
              prior_data_bym2_space = prior_data_bym2_space,
              prior_data_interaction = prior_data_interaction,
              eval_joint_prior = makemyprior::eval_joint_prior,
              hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
              calc_jac_logdet = makemyprior:::calc_jac_logdet,
              choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
              cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
              expit = makemyprior:::expit,
              get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
              get_indexes = makemyprior:::get_indexes,
              get_indexes2 = makemyprior:::get_indexes2,
              hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
              hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
              eval_spline_lpdf = makemyprior:::eval_spline_lpdf
            )))

#plot(res)

save(res,
     file = "./typeV/one_dataset_fitted_typeV_sc8.RData")



# Load in a singular data set
load("Simulated_data/sc10/sc10_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, dataset_id_2]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time

lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA

res <- inla(typeV_formula, 
            data = lambda_, 
            verbose = T,
            family = "poisson",
            E = lambda_$E_it,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE),
            control.expert = list(jp = inla.jp.define(
              prior_func,
              prior_data_bym2_time = prior_data_bym2_time,
              prior_data_bym2_space = prior_data_bym2_space,
              prior_data_interaction = prior_data_interaction,
              eval_joint_prior = makemyprior::eval_joint_prior,
              hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
              calc_jac_logdet = makemyprior:::calc_jac_logdet,
              choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
              cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
              expit = makemyprior:::expit,
              get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
              get_indexes = makemyprior:::get_indexes,
              get_indexes2 = makemyprior:::get_indexes2,
              hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
              hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
              eval_spline_lpdf = makemyprior:::eval_spline_lpdf
            )))

#plot(res)

save(res,
     file = "./typeV/one_dataset_fitted_typeV_sc10.RData")


# Load in a singular data set
load("Simulated_data/sc12/sc12_data.RData")

### Format data for makemyprior
lambda_ = lambda.df[, c("area_id", "time_id",
                        "E_it", "space.time")]


lambda_$sampled_counts = lambda.df$sampled_counts[, dataset_id_2]
lambda_$area_id_iid = lambda_$area_id
lambda_$area_id_struct = lambda_$area_id
lambda_$time_id_iid = lambda_$time_id
lambda_$time_id_struct = lambda_$time_id
lambda_$space.time_I = lambda_$space.time
lambda_$space.time_IV = lambda_$space.time

lambda_[lambda_$time_id %in% 11:13, ]$sampled_counts = NA

res <- inla(typeV_formula, 
            data = lambda_, 
            verbose = T,
            family = "poisson",
            E = lambda_$E_it,
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE),
            control.expert = list(jp = inla.jp.define(
              prior_func,
              prior_data_bym2_time = prior_data_bym2_time,
              prior_data_bym2_space = prior_data_bym2_space,
              prior_data_interaction = prior_data_interaction,
              eval_joint_prior = makemyprior::eval_joint_prior,
              hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf,
              calc_jac_logdet = makemyprior:::calc_jac_logdet,
              choose_prior_lpdf = makemyprior:::choose_prior_lpdf,
              cw_priors_lpdf = makemyprior:::cw_priors_lpdf,
              expit = makemyprior:::expit,
              get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter,
              get_indexes = makemyprior:::get_indexes,
              get_indexes2 = makemyprior:::get_indexes2,
              hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf,
              hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf,
              eval_spline_lpdf = makemyprior:::eval_spline_lpdf
            )))

#plot(res)

save(res,
     file = "./typeV/one_dataset_fitted_typeV_sc12.RData")

