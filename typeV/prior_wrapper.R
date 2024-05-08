####
# Print input to the prior wrapper as provided in the make_typeV... file 
# currently running

print(paste("Make the temporal BYM2 prior: ", make_temporal_BYM2_prior, sep = ""))

print(paste("Make the spatial BYM2 prior: ", make_spatial_BYM2_prior, sep = ""))

print(paste("Make Dirichlet space-time interaction prior: ", make_st_dirichlet_prior, sep = ""))

print(paste("Make reduced space-time interaction prior (PC): ", make_st_reduced_prior, sep = ""))

print(paste("Make full space-time interaction prior (PC): ", make_st_full_prior, sep = ""))

################################################################################





data_time <- lambda_[, c("time_id_iid", "time_id_struct")]
data_time$y <- 1:nrow(data_time) # er ikke saa noeye hva dataene er her, de skal uansett ikke brukes
if (make_temporal_BYM2_prior){ # denne tar mye tid
  prior_time <- make_prior(
    y ~ mc(time_id_iid, model = "iid") + 
      mc(time_id_struct, model = "besag", 
         graph = RW1_prec_path, 
         scale.model = T, 
         constr = T),
    data_time,
    family = "poisson",
    prior = list(
      tree = "s1 = (time_id_iid, time_id_struct)",
      V = list(s1 = list(prior = "pc", param = c(1, 0.01))), # default for poisson
      w = list(s1 = list(prior = "pc1", param = 0.75))
    )
  )
  
  prior_data_bym2_time <- make_eval_prior_data(prior_time)
}

# else { # denne er mye kjappere (er ikke 100% sikker paa at de blir like, men se kommentar for prior_space)
#   prior_time <- make_prior(
#     y ~ mc(time_id_iid, model = "iid") + 
#       mc(time_id_struct, model = "besag", graph = RW1_prec_path, scale.model = T, constr = T),
#     data_time[seq(1, 13*402, 402),], # hente ut bare unike indekser
#     family = "poisson",
#     prior = list(
#       tree = "s1 = (time_id_iid, time_id_struct)",
#       V = list(s1 = list(prior = "pc", param = c(1, 0.01))), # default for poisson
#       w = list(s1 = list(prior = "pc1", param = 0.75))
#     )
#   )
# }

print("Temporal part of prior made")

data_space <- lambda_[, c("area_id_iid", "area_id_struct")]
data_space$y <- 1:nrow(data_space) # er ikke saa noeye hva dataene er her, de skal uansett ikke brukes
if (make_spatial_BYM2_prior){ # denne tar mye tid
  prior_space <- make_prior(
    y ~ mc(area_id_iid, model = "iid") + 
      mc(area_id_struct, 
         model = "besag", 
         graph = Besag_prec_second_level_path, 
         scale.model = T, 
         constr = T),
    data_space,
    family = "poisson",
    prior = list(
      tree = "s1 = (area_id_iid, area_id_struct)",
      V = list(s1 = list(prior = "pc", param = c(1, 0.01))), # default for poisson
      w = list(s1 = list(prior = "pc1", param = 0.75))
    )
  )
  
  prior_data_bym2_space <- make_eval_prior_data(prior_space)
} 

# else { # denne er mye kjappere (jeg er ikke 100% sikker paa at de blir like, men denne og den med de 402*2 foerste radene blir like)
#   prior_space <- make_prior(
#     y ~ mc(area_id_iid, model = "iid") + 
#       mc(area_id_struct, model = "besag", graph = Besag_prec_second_level_path, scale.model = T, constr = T),
#     data_space[1:402,],
#     family = "poisson",
#     prior = list(
#       tree = "s1 = (area_id_iid, area_id_struct)",
#       V = list(s1 = list(prior = "pc", param = c(1, 0.01))), # default for poisson
#       w = list(s1 = list(prior = "pc1", param = 0.75))
#     )
#   )
# }

print("Spatial part of prior made")


if(make_st_dirichlet_prior){ # Make a 'stupid'-Dirichlet prior for splitting the variance between all four components
  data_interaction <- lambda_[, c("space.time_I", "space.time_II", "space.time_III", "space.time_IV")]
  data_interaction$y <- 1:nrow(data_interaction) # er ikke saa noeye hva dataene er her, de skal uansett ikke brukes
  
  
  prior_interaction <- make_prior(
      y ~ mc(space.time_I,     # Unstructured space-time interaction
             model = "iid") + 
      mc(space.time_II,    # Structured time interacting w. iid space
         model = "iid") +
      mc(space.time_III,   # iid time interacting w. structured space
         model = "iid") +
      mc(space.time_IV,    # Structured time interacting w. structured space
         model = "iid"),
    data_interaction[1:10,], # antall datapunkter ikke noeye naar vi bruker dirichlet, da dirichlet ikke bryr seg om struktur uansett
    family = "poisson",
    prior = list(
      tree = "s1 = (space.time_I, space.time_II, space.time_III, space.time_IV)",
      V = list(s1 = list(prior = "pc", param = c(1, 0.01))), # default for poisson
      w = list(s1 = list(prior = "dirichlet"))
    )
  )
  
} else if(make_st_reduced_prior){ # Make a prior for a reduced model, only having type I and IV interactions
  data_interaction <- lambda_[, c("space.time_I", "space.time_IV")]
  data_interaction$y <- 1:nrow(data_interaction) # er ikke saa noeye hva dataene er her, de skal uansett ikke brukes
  
  
  prior_interaction <- make_prior(
    y ~ mc(space.time_I,     # Unstructured space-time interaction
           model = "iid") + 
      mc(space.time_IV,    # Structured time interacting w. structured space
         model = "iid"),
    data_interaction[1:10,],
    family = "poisson",
    prior = list(
      tree = "s1 = (space.time_I, space.time_IV)",
      V = list(s1 = list(prior = "pc", param = c(1, 0.01))), 
      w = list(s1 = list(prior = "dirichlet"))
    )
  )

} else if(make_st_full_prior){ # Make a complete prior using a PC prior in first step to split between type I and rest, then Dirichlet
  prior_interaction <- make_prior(
    y ~ mc(space.time_I,     # Unstructured space-time interaction
           model = "iid") + 
      mc(space.time_II,    # Structured time interacting w. iid space
         model = "generic0",
         Cmatrix = typeII_prec_second_level,
         extraconstr = typeII_constraints_second_level,
         rankdef = nrow(second_level_admin_map)) +
      mc(space.time_III,   # iid time interacting w. structured space
         model = "generic0",
         Cmatrix = typeIII_prec_second_level,
         extraconstr = typeIII_constraints_second_level,
         rankdef = tT) +
      mc(space.time_IV,    # Structured time interacting w. structured space
         model = "generic0",
         Cmatrix = typeIV_prec_second_level,
         extraconstr = typeIV_constraints_second_level),
    data_interaction,
    family = "poisson",
    prior = list(
      tree = "s1 = (space.time_I, s2); s2 = (space.time_II, space.time_III, space.time_IV)",
      V = list(s1 = list(prior = "pc", param = c(1, 0.01))), 
      w = list(s1 = list(prior = "pc0", param = c(0.5, 0.5)),
               s2 = list(prior = "dirichlet"))
    ) 
  )
}

prior_data_interaction <- make_eval_prior_data(prior_interaction)









# if (FALSE){ # denne tar mye tid
#   prior_interaction <- make_prior(
#     y ~ mc(space.time_I,     # Unstructured space-time interaction
#            model = "iid") + 
#       mc(space.time_II,    # Structured time interacting w. iid space
#          model = "generic0",
#          Cmatrix = typeII_prec_second_level,
#          extraconstr = typeII_constraints_second_level,
#          rankdef = nrow(second_level_admin_map)) +
#       mc(space.time_III,   # iid time interacting w. structured space
#          model = "generic0",
#          Cmatrix = typeIII_prec_second_level,
#          extraconstr = typeIII_constraints_second_level,
#          rankdef = tT) +
#       mc(space.time_IV,    # Structured time interacting w. structured space
#          model = "generic0",
#          Cmatrix = typeIV_prec_second_level,
#          extraconstr = typeIV_constraints_second_level),
#     data_interaction,
#     family = "poisson",
#     prior = list(
#       tree = "s1 = (space.time_I, space.time_II, space.time_III, space.time_IV)",
#       V = list(s1 = list(prior = "pc", param = c(1.6, 0.05))), # default for poisson
#       w = list(s1 = list(prior = "dirichlet"))
#     )
#   )
#   # her deles variansen foerst "blindt" (dirichlet) til to noder som deretter deler den videre med en PC prior
#   # (jeg aaaaner ikke om dette gir mening i praksis)
#   prior_interaction_pc <- make_prior( # eksempel pC% pc prior, denne vil vC&re treeeeig
#     y ~ mc(space.time_I,     # Unstructured space-time interaction
#            model = "iid") + 
#       mc(space.time_II,    # Structured time interacting w. iid space
#          model = "generic0",
#          Cmatrix = typeII_prec_second_level,
#          extraconstr = typeII_constraints_second_level,
#          rankdef = nrow(second_level_admin_map)) +
#       mc(space.time_III,   # iid time interacting w. structured space
#          model = "generic0",
#          Cmatrix = typeIII_prec_second_level,
#          extraconstr = typeIII_constraints_second_level,
#          rankdef = tT) +
#       mc(space.time_IV,    # Structured time interacting w. structured space
#          model = "generic0",
#          Cmatrix = typeIV_prec_second_level,
#          extraconstr = typeIV_constraints_second_level),
#     data_interaction,
#     family = "poisson",
#     prior = list(
#       tree = "s1 = (s2, s3); s2 = (space.time_I, space.time_II); s3 = (space.time_III, space.time_IV))",
#       V = list(s1 = list(prior = "pc", param = c(1.6, 0.05))), # default for poisson
#       w = list(s1 = list(prior = "dirichlet"),
#                s2 = list(prior = "pcM", param = c(0.5, 0.75)),
#                s3 = list(prior = "pcM", param = c(0.5, 0.75)))
#     )
#   )
# } 



print("Space-time interaction part of prior made")


# obs! merk at rekkefoelgen paa inla-formelen maa matche rekkefoelgen vi bruker log-presisjonene her!
#prior_func <- function(logprec){
#  return(
#    eval_joint_prior(-c(logprec[1:2]), prior_data_bym2_time) +
#      eval_joint_prior(-c(logprec[3:4]), prior_data_bym2_space) +
#      eval_joint_prior(-c(logprec[5:8]), prior_data_interaction)
#  )
#}

#print("Joint prior finished")


# f1 <- 
#   sampled_counts ~ 1 + f(time_id_iid,      # Unstructured temporal effect
#                          model = "iid") + 
#   f(time_id_struct,   # Structured temporal effect
#     model = "besag", 
#     graph = RW1_prec_path,
#     scale.model = TRUE,
#     constr = T) + 
#   f(area_id_iid,      # Unstructured spatial effect
#     model = "iid") + 
#   f(area_id_struct,   # Structured spatial effect
#     model = "besag",
#     graph = Besag_prec_second_level_path,
#     scale.model = T,
#     constr = T) + 
#   f(space.time_I,     # Unstructured space-time interaction
#     model = "iid") + 
#   f(space.time_II,    # Structured time interacting w. iid space
#     model = "generic0",
#     Cmatrix = typeII_prec_second_level,
#     extraconstr = typeII_constraints_second_level,
#     rankdef = nrow(second_level_admin_map)) +
#   f(space.time_III,   # iid time interacting w. structured space
#     model = "generic0",
#     Cmatrix = typeIII_prec_second_level,
#     extraconstr = typeIII_constraints_second_level,
#     rankdef = tT) +
#   f(space.time_IV,    # Structured time interacting w. structured space
#     model = "generic0",
#     Cmatrix = typeIV_prec_second_level,
#     extraconstr = typeIV_constraints_second_level)
# 


# uten data for aa sjekke at prior er som forventet
# lambda_$sampled_counts2 <- lambda_$sampled_counts
# lambda_$sampled_counts <- NA
# res_prior <- inla(f1, data = lambda_, family = "poisson",
#             control.expert = list(jp = inla.jp.define(
#               prior_func, 
#               prior_data_bym2_time = prior_data_bym2_time,
#               prior_data_bym2_space = prior_data_bym2_space,
#               prior_data_interaction = prior_data_interaction,
#               eval_joint_prior = makemyprior::eval_joint_prior,
#               hd_prior_joint_lpdf = makemyprior:::hd_prior_joint_lpdf, 
#               calc_jac_logdet = makemyprior:::calc_jac_logdet, 
#               choose_prior_lpdf = makemyprior:::choose_prior_lpdf, 
#               cw_priors_lpdf = makemyprior:::cw_priors_lpdf, 
#               expit = makemyprior:::expit, 
#               get_dirichlet_parameter = makemyprior:::get_dirichlet_parameter, 
#               get_indexes = makemyprior:::get_indexes, 
#               get_indexes2 = makemyprior:::get_indexes2, 
#               hd_dirichlet_prior_lpdf = makemyprior:::hd_dirichlet_prior_lpdf, 
#               hd_pc_prior_lpdf = makemyprior:::hd_pc_prior_lpdf, 
#               eval_spline_lpdf = makemyprior:::eval_spline_lpdf
#             )))
# 
# samps <- inla.hyperpar.sample(n = 1000, res_prior) # husk at dette naa blir paa presisjons-skala, du maa regne om til 
# # varianser og vekter selv!! siden vi har samples kan du transformere direkte og saa plotte sammen med prior fra makemyprior
# 
# # hente ut ferdig evaluerte marginale PC prior-fordelinger for vekter fra makemyprior:
# eval_pc_prior(x = seq(0, 1, 0.01), obj = prior_space, param = "w[area_id_iid/area_id_iid_area_id_struct]")
# 
# 
# 
# 
# 
# lambda_$sampled_counts <- lambda_$sampled_counts2 # legge dataene tilbake
# # sjekk at dette faktisk stemmer, jeg tar meg ikke tid til aa vente paa resultatene hehe
# res <- inla(f1, data = lambda_, family = "poisson", # E = ...,
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

