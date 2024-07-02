#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library(bigDM)
library(latex2exp)
library(tables)
library("geofacet")
library(ggh4x)
library(ggridges)
library(geoR)
library(paletteer)
library('ggsci')
library(ggstats)
library(tmap)
library(RColorBrewer)

################################################################################

### Load in data and Map
# Load in the considered lung-cancer data
data(Data_LungCancer) # Areas in numerically increasing order in terms of ID for each year?

# Load in a map of Spain containing 7907 areas
data("Carto_SpainMUN") # Areas in numerically increasing order in terms of ID?

### Remove disjointed area Llivia
# Find the ID (and other aspects of disjointed area)
problem_area = na.omit(Carto_SpainMUN[Carto_SpainMUN$name == "Llivia", ]) 

# Remove it from the map! (it has index 2454 in Carto)
Carto_SpainMUN <- Carto_SpainMUN[-as.integer(rownames(problem_area)), ]

# Remove Llivia from Data_lungCancer
Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID != problem_area$ID, ]

### Extract the data within the principality of Extremadura
# Extract for the map
Carto_SpainMUN <- Carto_SpainMUN[Carto_SpainMUN$region == "Extremadura", ]

# Find the IDs of the areas within Extremadura
IDs_Extremadura <- unique(Carto_SpainMUN$ID) # IDs either 06... or 10...

# Extract the lung cancer data within Extremadura
Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID %in% IDs_Extremadura, ]



### Get the number of years, start-year and end-year.
# Get the unique years
years = unique(Data_LungCancer$year)

# Get the number of years
tT = length(years)

# Get start-year and end-year
t.from <- min(Data_LungCancer$year)
t.to <- max(Data_LungCancer$year)

### Create IDs for INLA-effects
# Create time-ids 1,...,25 instead of 1991,...,2015 for the sake of INLA
Data_LungCancer$year_id <- Data_LungCancer$year - min(Data_LungCancer$year) + 1


### add area IDs and find number of areas
# find number of areas
n = nrow(Carto_SpainMUN)

# Add unique area_ids to each area in Carto_SpainMUN (that is more easily read for me i.e. 1,...,380)
Carto_SpainMUN$area_id = 1:n

# Add area_id to Data_LungCancer
Data_LungCancer$area_id = rep(NA, nrow(Data_LungCancer))
for(year_id in 1:tT){
  # Insert area_id for Data_LungCancer that corresponds to Carto_SpainMUN
  for(i in 1:nrow(Carto_SpainMUN)){
    Data_LungCancer[Data_LungCancer$year_id == year_id & 
                      Data_LungCancer$ID == Carto_SpainMUN[i, ]$ID, ]$area_id = Carto_SpainMUN[i, ]$area_id
  }
}


### Reset rownames and add space.time id
row.names(Carto_SpainMUN) = 1:nrow(Carto_SpainMUN)
row.names(Data_LungCancer) = 1:nrow(Data_LungCancer)
Data_LungCancer$space.time = 1:nrow(Data_LungCancer)


### Set the last five years of observations to zero - as we are to predict on them
# First store the last five years
last_five_years_obs <- Data_LungCancer[Data_LungCancer$year_id %in% 21:25, ]$obs 

#### Weird things relating to the data

for(t in 1:tT){ # How does the average population change over the years
  print(paste("year: ", t, " Average pop: ", mean(Data_LungCancer[Data_LungCancer$year_id == t, ]$pop), " Median pop: ", median(Data_LungCancer[Data_LungCancer$year_id == t, ]$pop)))
} # In the last five years, the average population drops of slightly, median pop also lower

for(t in 1:tT){ # How does the average obs change over the years
  print(paste("year: ", t, " Average obs: ", mean(Data_LungCancer[Data_LungCancer$year_id == t, ]$obs), " Median obs: ", median(Data_LungCancer[Data_LungCancer$year_id == t, ]$obs)))
}




################################################################################

###Specify priors for hyperparameters of improper models
#---
# Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

# Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

# Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))
#---

### Specify precision matrix
#---
# Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)


# Find the adjacencies
nb_spain <- spdep::poly2nb(Carto_SpainMUN, queen = T)


### Make precision matrix for Besag on Extremadura map
matrix4inla <- nb2mat(nb_spain, style = "B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_spain<- Matrix(matrix4inla, sparse = TRUE) #Make it sparse
rm(matrix4inla)
#---

### Plot the neighborhod structure (See that Besag looks all right)
plot(st_geometry(Carto_SpainMUN))  
plot(nb_spain, coords = st_geometry(Carto_SpainMUN), add = T)



### Get necessities for type IV interaction (constraints and prec. matrix)
#Get sum-to-zero constraints for type IV interaction
typeIV_constraints = constraints_maker(type = "IV", 
                                       n = n, 
                                       t = tT)


#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))

# get scaled Besag
scaled_besag_prec <- INLA::inla.scale.model(Besag_prec_spain, 
                                            constr = list(A = matrix(1,1,dim(Besag_prec_spain)[1]), 
                                                          e = 0))

#Get type IV interaction precision matrix
typeIV_prec <- scaled_RW_prec %x% scaled_besag_prec


## Specify type IV formula
typeIV_formula <- obs ~ 1 + f(year_id, 
                              model = 'bym2',
                              scale.model = T, 
                              constr = T, 
                              rankdef = 1,
                              graph = RW1_prec,
                              hyper = temporal_hyper) + 
                            f(area_id, 
                              model = 'bym2',
                              scale.model = T,
                              constr = T,
                              rankdef = 1,
                              graph = Besag_prec_spain,
                              hyper = spatial_hyper) + 
                            f(space.time, 
                              model = "generic0",
                              Cmatrix = typeIV_prec,
                              extraconstr = typeIV_constraints,
                              rankdef = (n + tT - 1), 
                              hyper = interaction_hyper)

# Then set last five years of obs to NA
Data_LungCancer[Data_LungCancer$year_id %in% 21:25, ]$obs = NA

### Fit INLA object
Improper1_typeIV_Extremadura = inla(typeIV_formula, 
                                    data = Data_LungCancer, 
                                    family = "poisson",
                                    E = pop, 
                                    verbose = T,
                                    control.predictor = list(compute = TRUE, 
                                                             link = 1),       #For predictions
                                    control.compute = list(return.marginals.predictor=TRUE)) # Get the lin.pred.marginal

#Re-insert last five years of observations 
Data_LungCancer[Data_LungCancer$year_id %in% 21:25, ]$obs = last_five_years_obs


# Quick sanity check
plot(Improper1_typeIV_Extremadura)

### Save the INLA object
save(Improper1_typeIV_Extremadura,
     file = "./case_study/Improper1_typeIV_Extremadura.RData")


################################################################################

load("case_study/Improper1_typeIV_Extremadura.RData")


### Find coverage of model



### Find number of upper quantiles equal to 2,3,... for each year, likewise for l




### Find the distances of the observations outside of the 95% pred. CIs for each year




calc_log_scores <- function(model, data,
                            n, tT, 
                            Improper = T){
  
  ### Initialize data frame to store the MSE's, IS's, and log-scores 1, 2, 3, 4, and 5 years ahead and average
  model_choice.df <- data.frame(log_1 = rep(NA, n), log_2 = rep(NA, n), log_3 = rep(NA, n),
                                log_4 = rep(NA, n), log_5 = rep(NA, n), log_avg = rep(NA, n))
  
  
  # If proper models, sort the marginals to get the correct order
  if(!Improper){
    print("Sorting")
    model$marginals.fitted.values <- sort_proper_fitted(model$marginals.fitted.values,
                                                        n, tT)
  }
  
  ### Take the years predicted on
  years_pred_on <- 21:25
  
  
  #### For each year calculate the MSE and IS that year
  for(year in years_pred_on){
    print(paste("Year: ", year, sep = ""))
    
    ## Extract the predicted marginals for this year
    one_year_margs = model$marginals.fitted.values[((year - 1) * n + 1):(year * n)]
    
    ## Extract the observed counts for this year
    one_year_counts = data$obs[((year - 1) * n + 1):(year * n)]
    
    ## Extract the populations for this year
    one_year_pop  = data$pop[((year - 1) * n + 1):(year * n)]
    
    ### Calculate the year-ahead log-scores 
    model_choice.df[, year - years_pred_on[1] + 1] =  count_log_s_one_year(counts = one_year_counts,
                                                                           marginals = one_year_margs,
                                                                           population = one_year_pop)
    
  }
  
  #Get the total MSE for count
  model_choice.df[, 6] = rowMeans(model_choice.df[, 1:5])
  
  
  return(model_choice.df)
}



tmp <- calc_log_scores(Improper1_typeIV_Extremadura,
                Data_LungCancer,
                n, tT, 
                Improper = T)




colMeans(tmp)



## Calculate model choice critera: MSE, IS
calc_model_choice <- function(model,
                              data,
                              n, tT,
                              Improper = T){
  
  ### Initialize data frame to store the MSE's, IS's, and log-scores 1, 2, 3, 4, and 5 years ahead and average
  model_choice.df <- data.frame(mse_1 = NA, mse_2 = NA, mse_3 = NA,
                                mse_4 = NA, mse_5 = NA, mse_avg = NA,
                                is_1 = NA, is_2 = NA, is_3 = NA,
                                is_4 = NA, is_5 = NA, is_avg = NA)
  
  
  # If proper models, sort the marginals to get the correct order
  if(!Improper){
    print("Sorting")
    model$marginals.fitted.values <- sort_proper_fitted(model$marginals.fitted.values,
                                                        n, tT)
  }
  
  ### Take the years predicted on
  years_pred_on <- 21:25
  
  
  #### For each year calculate the MSE and IS that year
  for(year in years_pred_on){
    print(paste("Year: ", year, sep = ""))
    
    ## Extract the predicted marginals for this year
    one_year_margs = model$marginals.fitted.values[((year - 1) * n + 1):(year * n)]
    
    ## Extract the observed counts for this year
    one_year_counts = data$obs[((year - 1) * n + 1):(year * n)]
    
    ## Extract the populations for this year
    one_year_pop  = data$pop[((year - 1) * n + 1):(year * n)]
    
    ## MSE
    model_choice.df[1, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(one_year_counts,
                                                                                      one_year_margs,
                                                                                      E_it = one_year_pop)
    
    
    ## IS for count in year ahead 
    model_choice.df[1, year - years_pred_on[1] + 7] =  count_IS_one_year_case_study(counts = one_year_counts, 
                                                                                    marginals = one_year_margs, 
                                                                                    population = one_year_pop)
    
    
    # Check the width of the 95% CI of the predicted rate. If that is increasing over the years
    # Then, the 95% CI for the predicted number of counts ought to increase as well!!! 
    
    print(paste("avg. IS: ", model_choice.df[1, year - years_pred_on[1] + 7], " avg. MSE: ", model_choice.df[1, year - years_pred_on[1] + 1],
                " avg. CI width: ", width_CI_one_year_one_dataset(one_year_margs)))
    
    
    
  }
  
  #Get the total MSE for count
  model_choice.df[1, 6] = mean(as.numeric(model_choice.df[1, 1:5]))
  
  
  #Get the total IS for count
  model_choice.df[1, 12] = mean(as.numeric(model_choice.df[1, 7:11]))
  
  return(model_choice.df)
}


tmp2 <- calc_model_choice(Improper1_typeIV_Extremadura,
                          Data_LungCancer,
                          n, tT, 
                          Improper = T)





################################################################################

### Plot the posterior predictive 95% CIs for (y_it/E_it = lambda_it) for 6 regions

#plt_predcount_vs_true_count <- function(geofacet_grid,
#                                        pred_to_plot,
#                                        title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true counts

Extremadura_grid <- data.frame(
  code = c(1, 6, 203,
           257, 309, 337),
  name = c("Acedera", "Alburquerque", "Caminomorisco",
           "Hervas", "Plasencia", "Serrejon"),
  row = c(1, 1, 1, 2, 2, 2),
  col = c(1, 2, 3, 1, 2, 3),
  stringsAsFactors = FALSE)


geofacet::grid_preview(Extremadura_grid)


pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, #lambda_$area_id
                           time_id = Data_LungCancer$year_id, # lambda_$time_id
                           median = Improper1_typeIV_Extremadura$summary.fitted.values$'0.5quant', #sort_proper_fitted(model$summary.fitted.values$'0.5quant', length(unique(lambda_$area_id)), tT) * lambda_$E_it, # model$summary.fitted.values$'0.5quant',
                           quantile_0.025 = Improper1_typeIV_Extremadura$summary.fitted.values$'0.025quant',# model$summary.fitted.values$'0.025quant',
                           quantile_0.975 = Improper1_typeIV_Extremadura$summary.fitted.values$'0.975quant',
                           y_it_div_E_it = Data_LungCancer$obs/Data_LungCancer$pop)


  
ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
    geom_point((aes(x = time_id, y = y_it_div_E_it, col = "y/E"))) +
    geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
    geom_vline(xintercept = 20.5, linetype = "longdash", 
               color = "darkgrey", linewidth = 0.6) +
    facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
    labs(title = "HEIHEI",
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=12),
          plot.title = element_text(hjust = 0.5, size=12),
          strip.text.x = element_text(size = 10),
          legend.position = "right",
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  
plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue"),
                                labels = unname(TeX(c("Posterior 95% CI",
                                                      "Posterior mean $\\hat{Y}_{it}$",
                                                      "True $Y_{it}$",
                                                      "Expected true $Y_{it}$")))) # 



### Posterior predicted counts

# Need to calculate the quantiles for the predictive distribution
## Get the upper, lower, and median quantile for the pred. counts
#ul_each <- lapply(Improper1_typeIV_Extremadura$marginals.fitted.values, 
#                  FUN = function(x){
#                    return(find_ul_quants_counts_single_pred(x, 100)) #100 aka E_it
#                  })

pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           #post_mean = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.025 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.975 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           y_it = Data_LungCancer$obs)


for(area_id in Extremadura_grid$code){
  for(t in 1:tT){
    ul <- find_ul_quants_counts_single_pred(Improper1_typeIV_Extremadura$marginals.fitted.values[[(t - 1) * n + area_id]],
                                            Data_LungCancer[Data_LungCancer$area_id == area_id & Data_LungCancer$year_id == t, ]$pop)
    
    
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$median = ul$median
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.975 = ul$u
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.025 = ul$l
  }
}




ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  geom_point((aes(x = time_id, y = y_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = "HEIHEI",
       x = "Year",
       y = "Rate",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        strip.text.x = element_text(size = 10),
        legend.position = "right",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))



### Posterior predicted counts per 100,000


pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           #post_mean = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.025 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.975 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           y_it = Data_LungCancer$obs/Data_LungCancer$pop * 1E5)


for(area_id in Extremadura_grid$code){
  for(t in 1:tT){
    ul <- find_ul_quants_counts_single_pred(Improper1_typeIV_Extremadura$marginals.fitted.values[[(t - 1) * n + area_id]],
                                            1E5)
    
    
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$median = ul$median
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.975 = ul$u
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.025 = ul$l
  }
}




ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  geom_point((aes(x = time_id, y = y_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = "HEIHEI",
       x = "Year",
       y = "Rate",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        strip.text.x = element_text(size = 10),
        legend.position = "right",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))








poppy = 1E4


pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           #post_mean = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.025 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.975 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           y_it = Data_LungCancer$obs/Data_LungCancer$pop * poppy)




for(area_id in Extremadura_grid$code){
  for(t in 1:tT){
    ul <- find_ul_quants_counts_single_pred(Improper1_typeIV_Extremadura$marginals.fitted.values[[(t - 1) * n + area_id]],
                                            Data_LungCancer[(t - 1) * n + area_id, ]$pop)
    
    
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$median = (ul$median/Data_LungCancer[(t - 1) * n + area_id, ]$pop) * poppy
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.975 = (ul$u/Data_LungCancer[(t - 1) * n + area_id, ]$pop) * poppy
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.025 = (ul$l/Data_LungCancer[(t - 1) * n + area_id, ]$pop) * poppy
  }
}




ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  geom_point((aes(x = time_id, y = y_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = "HEIHEI",
       x = "Year",
       y = "Rate",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        strip.text.x = element_text(size = 10),
        legend.position = "right",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))








Data_LungCancer[Data_LungCancer$area_id == 1, ]$pop

















################################################################################

####
load("case_study/proper2_RW1_Extremadura.RData")

if(FALSE){ # if(FALSE) added so that dont sort them unless I really want to
  proper2_RW1_Extremadura$marginals.fitted.values <- sort_proper_fitted(proper2_RW1_Extremadura$marginals.fitted.values,
                                                                        n, tT)
  
  proper2_RW1_Extremadura$summary.fitted.values$mean <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$mean,
                                                                           n, tT)
  
  proper2_RW1_Extremadura$summary.fitted.values$sd <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$sd,
                                                                         n, tT)
  
                                                                            
}

median_ <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$'0.5quant', n, tT)
l <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$'0.025quant', n, tT)
u <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$'0.975quant', n, tT)


pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = median_, 
                           quantile_0.025 = l,
                           quantile_0.975 = u,
                           y_it_div_E_it = Data_LungCancer$obs/Data_LungCancer$pop)



plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  #geom_point((aes(x = time_id, y = y_it_div_E_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = "HEIHEI",
       x = "Year",
       y = "Rate",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        strip.text.x = element_text(size = 10),
        legend.position = "right",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))

plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue"),
                                labels = unname(TeX(c("Posterior 95% CI",
                                                      "Posterior mean $\\hat{Y}_{it}$",
                                                      "True $Y_{it}$",
                                                      "Expected true $Y_{it}$")))) # 



### Posterior predicted counts

# Need to calculate the quantiles for the predictive distribution
## Get the upper, lower, and median quantile for the pred. counts
#ul_each <- lapply(Improper1_typeIV_Extremadura$marginals.fitted.values, 
#                  FUN = function(x){
#                    return(find_ul_quants_counts_single_pred(x, 100)) #100 aka E_it
#                  })

pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           #post_mean = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.025 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.975 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           y_it = Data_LungCancer$obs)


for(area_id in Extremadura_grid$code){
  for(t in 1:tT){
    ul <- find_ul_quants_counts_single_pred(Improper1_typeIV_Extremadura$marginals.fitted.values[[(t - 1) * n + area_id]],
                                            Data_LungCancer[Data_LungCancer$area_id == area_id & Data_LungCancer$year_id == t, ]$pop)
    
    
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$median = ul$median
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.975 = ul$u
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.025 = ul$l
  }
}




plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  geom_point((aes(x = time_id, y = y_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = "HEIHEI",
       x = "Year",
       y = "Rate",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        strip.text.x = element_text(size = 10),
        legend.position = "right",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))



### Posterior predicted counts per 100,000


pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           #post_mean = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.025 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.975 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           y_it = Data_LungCancer$obs/Data_LungCancer$pop * 1E5)


for(area_id in Extremadura_grid$code){
  for(t in 1:tT){
    ul <- find_ul_quants_counts_single_pred(Improper1_typeIV_Extremadura$marginals.fitted.values[[(t - 1) * n + area_id]],
                                            1E5)
    
    
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$median = ul$median
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.975 = ul$u
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.025 = ul$l
  }
}




plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  geom_point((aes(x = time_id, y = y_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = "HEIHEI",
       x = "Year",
       y = "Rate",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12),
        strip.text.x = element_text(size = 10),
        legend.position = "right",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))



################################################################################
### Plot the 10%-, 50%- and 90% empirical quantiles of the population

pop.df <- data.frame(year = t.from:t.to,
                     u = rep(NA, tT),
                     median = rep(NA, tT),
                     l = rep(NA, tT))

for(year in t.from:t.to){
  
  pop_quants_in_year <- quantile(Data_LungCancer[Data_LungCancer$year == year, ]$pop,
                                 probs = c(0.25, 0.5, 0.75))
  
  pop.df[pop.df$year == year, ]$l = pop_quants_in_year[[1]]
  pop.df[pop.df$year == year, ]$median = pop_quants_in_year[[2]]
  pop.df[pop.df$year == year, ]$u = pop_quants_in_year[[3]]
}

ggplot(data = pop.df) + 
  geom_line(aes(x = year, y = median)) + 
  geom_line(aes(x = year, y = l), linetype = "dashed") + 
  geom_line(aes(x = year, y = u), linetype = "dashed") + 
  theme_bw()










