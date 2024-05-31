#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("geofacet")
library(ggh4x)
library(ggridges)
library(latex2exp)
library(geoR)
library(paletteer)
library('ggsci')
library(ggstats)
library(tidyverse)

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

# Some useful numbers
tT = 13
n_sim = 100
E_it = 100

n_ADM4 <- nrow(second_level_admin_map)

dataset_id_2 = 4


scenario_names_ADM4 <- c("sc2", "sc4", "sc6",
                         "sc8", "sc10", "sc12")

################################################################################

load("typeV/one_dataset_fitted_typeV_sc2.RData")
load("diagnostics_sc2.RData")
rm(proper1_noInt_ADM4, proper1_full_ADM4, #proper1_onlyInt_ADM4,
   proper2_full_ADM4, proper2_noInt_ADM4, proper2_onlyInt_ADM4,
   Improper1_noInt_ADM4, Improper1_typeII_ADM4,Improper2_typeIII_ADM4,
   Improper2_typeII_ADM4, Improper2_typeI_ADM4)



plot(res)
plot(Improper1_typeIV_ADM4)
plot(proper1_onlyInt_ADM4)

proper1_onlyInt_ADM4$summary.hyperpar

#res$marginals.hyperpar$`Precision for space.time_I`

#samps <- inla.hyperpar.sample(n = 3000, res)

marg_var_IV <- inla.tmarginal(fun = function(x){1/x},
                              res$marginals.hyperpar$`Precision for space.time_IV`)

marg_var_I <- inla.tmarginal(fun = function(x){1/x},
                             res$marginals.hyperpar$`Precision for space.time_I`)


var_IV_samps <- inla.rmarginal(n = 2000, marg_var_IV)
mean(var_IV_samps)

var_I_samps <- inla.rmarginal(n = 2000, marg_var_I)
mean(var_I_samps)

var_samps <- var_I_samps + var_IV_samps

weight_st_samps_typeIV <- var_IV_samps/(var_I_samps + var_IV_samps)


var_samps.df <- data.frame(var_samps = var_samps)
variance_plot <- ggplot(var_samps.df,
                     aes(x = var_samps)) + 
                geom_density(fill = "#e9acaf", alpha = 0.8) + 
                theme_bw() + 
                xlab(TeX(r'($\sigma_{st}^2$)')) + 
                ylab(TeX(r'(Posterior density: $\sigma_{st}^2$)'))



tmp.df <- data.frame(type = c(rep("Type IV", length(var_IV_samps)),
                              rep("Type I", length(var_I_samps))),
                     var_samps <- c(var_IV_samps, var_I_samps))

split_variance_plot <- ggplot(tmp.df,
                             aes(x = var_samps)) + 
                        geom_density(aes(fill = type), 
                                     alpha = 0.8) + 
                        theme_bw() + 
                        theme(legend.position = "top") + 
                        xlab(TeX(r'($\sigma_{*}^2$)')) + 
                        ylab(TeX(r'(Posterior density: $\sigma_{*}^2$)')) + 
                        guides(fill = guide_legend(title="*"))



weight_st_samps_typeIV.df <- data.frame(weight_st_samps_typeIV = weight_st_samps_typeIV)

weight_typeIV <- ggplot(tmp.df,
                        aes(x = weight_st_samps_typeIV)) + 
                    geom_density(fill = "#69b3a2", alpha = 0.8) + 
                    theme_bw() + 
                    xlab(TeX(r'($w_{IV}$)')) + 
                    ylab(TeX(r'(Posterior density: $w_{IV}$)'))

# Save as 10 by 3.5 typeV_variance
ggarrange(split_variance_plot, variance_plot, weight_typeIV,
          ncol = 3, nrow = 1, 
          widths = c(1, 1, 1))

















