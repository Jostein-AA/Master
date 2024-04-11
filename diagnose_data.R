
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

dataset_id = 3

################################################################################

#Load in data set

#Load in simulation data (done stupidly, change!)
load("./Simulated_data/sc1/sc1_data.RData")
lambda_sc1.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc1.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc1.df$mu = lambda.df$mu[, dataset_id]
lambda_sc1.df$lambda_it = lambda.df$lambda_it[, dataset_id]


load("./Simulated_data/sc2/sc2_data.RData")
lambda_sc2.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc2.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc2.df$mu = lambda.df$mu[, dataset_id]
lambda_sc2.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc3/sc3_data.RData")
lambda_sc3.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc3.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc3.df$mu = lambda.df$mu[, dataset_id]
lambda_sc3.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc4/sc4_data.RData")
lambda_sc4.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc4.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc4.df$mu = lambda.df$mu[, dataset_id]
lambda_sc4.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc5/sc5_data.RData")
lambda_sc5.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc5.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc5.df$mu = lambda.df$mu[, dataset_id]
lambda_sc5.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc6/sc6_data.RData")
lambda_sc6.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc6.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc6.df$mu = lambda.df$mu[, dataset_id]
lambda_sc6.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc7/sc7_data.RData")
lambda_sc7.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc7.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc7.df$mu = lambda.df$mu[, dataset_id]
lambda_sc7.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc8/sc8_data.RData")
lambda_sc8.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc8.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc8.df$mu = lambda.df$mu[, dataset_id]
lambda_sc8.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc9/sc9_data.RData")
lambda_sc9.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc9.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc9.df$mu = lambda.df$mu[, dataset_id]
lambda_sc9.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc10/sc10_data.RData")
lambda_sc10.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc10.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc10.df$mu = lambda.df$mu[, dataset_id]
lambda_sc10.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc11/sc11_data.RData")
lambda_sc11.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc11.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc11.df$mu = lambda.df$mu[, dataset_id]
lambda_sc11.df$lambda_it = lambda.df$lambda_it[, dataset_id]

load("./Simulated_data/sc12/sc12_data.RData")
lambda_sc12.df <- lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
lambda_sc12.df$sampled_counts = lambda.df$sampled_counts[, dataset_id]
lambda_sc12.df$mu = lambda.df$mu[, dataset_id]
lambda_sc12.df$lambda_it = lambda.df$lambda_it[, dataset_id]

rm(lambda.df)

################################################################################
# Produce diagnostics plots

plot_temporal_trend_data_one_data_set <- function(sim_data, title){

  # Aggregate over the years
  aggr <- aggregate(sim_data, by = list(time_id = sim_data$time_id), FUN = mean) 
  
  # Plot over the years
  return(ggplot(data = aggr[, 2:ncol(aggr)]) + ggtitle(title) + 
    geom_line(aes(x = time_id, y = mu, col = "mu"), color = "blue") + 
    geom_point(aes(x = time_id, y = sampled_counts, col = "sampled counts"), color = "black") +
    ylim(7.5, 12.5))
}


plt1 <- plot_temporal_trend_data_one_data_set(lambda_sc1.df, "ADM1 const, short: one data set")
plt2 <- plot_temporal_trend_data_one_data_set(lambda_sc2.df, "ADM4 const, short: one data set")
plt3 <-plot_temporal_trend_data_one_data_set(lambda_sc3.df, "ADM1 lin increasing, short: one data set")
plt4 <-plot_temporal_trend_data_one_data_set(lambda_sc4.df, "ADM4 lin increasing, short: one data set")
plt5 <-plot_temporal_trend_data_one_data_set(lambda_sc5.df, "ADM1 change point, short: one data set")
plt6 <-plot_temporal_trend_data_one_data_set(lambda_sc6.df, "ADM4 change point, short: one data set")
plt7 <-plot_temporal_trend_data_one_data_set(lambda_sc7.df, "ADM1 const, long: one data set")
plt8 <-plot_temporal_trend_data_one_data_set(lambda_sc8.df, "ADM4 const, long: one data set")
plt9 <-plot_temporal_trend_data_one_data_set(lambda_sc9.df, "ADM1 lin increasing, long: one data set")
plt10 <-plot_temporal_trend_data_one_data_set(lambda_sc10.df, "ADM4 lin increasing, long: one data set")
plt11 <-plot_temporal_trend_data_one_data_set(lambda_sc11.df, "ADM1 change point, long: one data set")
plt12 <-plot_temporal_trend_data_one_data_set(lambda_sc12.df, "ADM4 change point, long: one data set")

ggarrange(plt1, plt2, plt3,
          plt4, plt5, plt6,
          plt7, plt8, plt9,
          plt10, plt11, plt12,
          ncol = 2, 
          nrow = 3)


plot_temporal_trend_data_all_data_set <- function(scenario_name, title){
  load(paste("./Simulated_data/", scenario_name,"/", scenario_name, "_data.RData",
             sep = ""))
  
  tmp_ = lambda.df[, c("space.time", "area_id", "time_id", "E_it", "sampled_counts", "mu")]
  
  # Aggregate the sampled counts and mu over the data sets
  
  tmp_$sampled_counts <- rowMeans(tmp_$sampled_counts)
  tmp_$mu <- rowMeans(tmp_$mu)
  
  # Aggregate over the years
  aggr <- aggregate(tmp_, by = list(time_id = tmp_$time_id), FUN = mean)
  
  # Plot over the years
  return(ggplot(data = aggr[, 2:ncol(aggr)]) + ggtitle(title) + 
           geom_line(aes(x = time_id, y = mu, col = "mu"), color = "blue") + 
           geom_point(aes(x = time_id, y = sampled_counts, col = "sampled counts"), color = "black") +
           ylim(7.5, 12.5))
}

plt1 <- plot_temporal_trend_data_all_data_set("sc1", "ADM1 const, short: all data sets")
plt2 <- plot_temporal_trend_data_all_data_set("sc2", "ADM4 const, short: all data sets")
plt3 <-plot_temporal_trend_data_all_data_set("sc3", "ADM1 lin increasing, short: all data sets")
plt4 <-plot_temporal_trend_data_all_data_set("sc4", "ADM4 lin increasing, short: all data sets")
plt5 <-plot_temporal_trend_data_all_data_set("sc5", "ADM1 change point, short: all data sets")
plt6 <-plot_temporal_trend_data_all_data_set("sc6", "ADM4 change point, short: all data sets")
plt7 <-plot_temporal_trend_data_all_data_set("sc7", "ADM1 const, long: all data sets")
plt8 <-plot_temporal_trend_data_all_data_set("sc8", "ADM4 const, long: all data sets")
plt9 <-plot_temporal_trend_data_all_data_set("sc9", "ADM1 lin increasing, long: all data sets")
plt10 <-plot_temporal_trend_data_all_data_set("sc10", "ADM4 lin increasing, long: all data sets")
plt11 <-plot_temporal_trend_data_all_data_set("sc11", "ADM1 change point, long: all data sets")
plt12 <-plot_temporal_trend_data_all_data_set("sc12", "ADM4 change point, long: all data sets")

ggarrange(plt1, plt2, plt3,
          plt4, plt5, plt6,
          plt7, plt8, plt9,
          plt10, plt11, plt12,
          ncol = 2, 
          nrow = 3)



################################################################################

E_it2 = 500


lambda_sc1.df$sampled_counts2 = sapply(lambda_sc1.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc2.df$sampled_counts2 = sapply(lambda_sc2.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc3.df$sampled_counts2 = sapply(lambda_sc3.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc4.df$sampled_counts2 = sapply(lambda_sc4.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc5.df$sampled_counts2 = sapply(lambda_sc5.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc6.df$sampled_counts2 = sapply(lambda_sc6.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc7.df$sampled_counts2 = sapply(lambda_sc7.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc8.df$sampled_counts2 = sapply(lambda_sc8.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc9.df$sampled_counts2 = sapply(lambda_sc9.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc10.df$sampled_counts2 = sapply(lambda_sc10.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc11.df$sampled_counts2 = sapply(lambda_sc11.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})
lambda_sc12.df$sampled_counts2 = sapply(lambda_sc12.df$lambda_it * E_it2, FUN = function(x){return(rpois(1, x))})




plot_alternate <- function(sim_data, title){
  
  # Aggregate over the years
  aggr <- aggregate(sim_data, by = list(time_id = sim_data$time_id), FUN = mean)
  
  aggr$mu = aggr$mu * E_it2/100
  
  # Plot over the years
  return(ggplot(data = aggr[, 2:ncol(aggr)]) + ggtitle(title) + 
           geom_line(aes(x = time_id, y = mu, col = "mu"), color = "blue") + 
           geom_point(aes(x = time_id, y = sampled_counts2, col = "sampled counts"), color = "black") +
           ylim(35, 65))
}


plt1 <- plot_alternate(lambda_sc1.df, "ADM1 const, short: E_it = 500")
plt2 <- plot_alternate(lambda_sc2.df, "ADM4 const, short: E_it = 500")
plt3 <-plot_alternate(lambda_sc3.df, "ADM1 lin increasing, short: E_it = 500")
plt4 <-plot_alternate(lambda_sc4.df, "ADM4 lin increasing, short: E_it = 500")
plt5 <-plot_alternate(lambda_sc5.df, "ADM1 change point, short: E_it = 500")
plt6 <-plot_alternate(lambda_sc6.df, "ADM4 change point, short: E_it = 500")
plt7 <-plot_alternate(lambda_sc7.df, "ADM1 const, long: E_it = 500")
plt8 <-plot_alternate(lambda_sc8.df, "ADM4 const, long: E_it = 500")
plt9 <-plot_alternate(lambda_sc9.df, "ADM1 lin increasing, long: E_it = 500")
plt10 <-plot_alternate(lambda_sc10.df, "ADM4 lin increasing, long: E_it = 500")
plt11 <-plot_alternate(lambda_sc11.df, "ADM1 change point, long: E_it = 500")
plt12 <-plot_alternate(lambda_sc12.df, "ADM4 change point, long: E_it = 500")

ggarrange(plt1, plt2, plt3,
          plt4, plt5, plt6,
          plt7, plt8, plt9,
          plt10, plt11, plt12,
          ncol = 2, 
          nrow = 3)

################################################################################
# Check how a couple of fitted models fair

### Temporal hyperparameters (Precision of iid and precision of RW1 and RW2) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

### Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))

### Temporal hyperparameters (prec. of AR1 and AR1's mixing param) w. corresponding priors: penalized constraint 
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), 
                 rho = list(prior = 'pc.cor1', 
                            param = c(0.5, 0.5 + 1E-2)))

### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper_proper = list(prec= list(prior = 'pc.prec', 
                                       param = c(1, 0.01)))

### Group hyper
group_hyper_ar1 = list(rho = list(prior = 'pc.cor1', 
                                  param = c(0.5, 0.5 + 1E-2)))




### Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)

### Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW1_prec <- inla.scale.model(RW1_prec,
                                    list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                         e = 0))


### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_first_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_first_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

### get scaled Besag on ADM1
scaled_besag_prec_first_level <- INLA::inla.scale.model(Besag_prec_first_level, 
                                                        constr = list(A = matrix(1,1,dim(Besag_prec_first_level)[1]), 
                                                                      e = 0))

### Get type IV interaction precision matrix
typeIV_prec_1_first_level <- scaled_RW1_prec %x% scaled_besag_prec_first_level

### Get type IV constraints
typeIV_constraints_1_first_level = constraints_maker(type = "IV", 
                                                     n = nrow(first_level_admin_map), 
                                                     t = tT)


Improper1_noInt_ADM1_formula <- sampled_counts2 ~ 1 + f(time_id, 
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
    graph = Besag_prec_first_level,
    hyper = spatial_hyper)

Improper1_typeIV_ADM1_formula <- update(Improper1_noInt_ADM1_formula,
                                        ~. + f(space.time, 
                                               model = "generic0",
                                               Cmatrix = typeIV_prec_1_first_level,
                                               extraconstr = typeIV_constraints_1_first_level,
                                               rankdef = (nrow(first_level_admin_map) + tT - 1), 
                                               hyper = interaction_hyper))





proper1_onlyInt_ADM1_formula <- sampled_counts2 ~ 1 + time_id + 
                                          f(area_id, 
                                            model = "besagproper2",
                                            graph = Besag_prec_first_level,
                                            hyper = spatial_hyper_proper,
                                            group = time_id, 
                                            control.group = list(model = "ar1",
                                                                 hyper = group_hyper_ar1))






lambda_sc3.df$E_it2 = E_it2

tmp_ = inla(Improper1_typeIV_ADM1_formula, 
            data = lambda_sc3.df, 
            family = "poisson",
            E = E_it2, #E_it
            control.predictor = list(compute = TRUE,
                                     link = 1),       #For predictions
            control.compute = list(config = TRUE, # To see constraints later
                                   cpo = T,       # For model selection
                                   return.marginals.predictor=TRUE))

library("geofacet")
select_regions_lin_pred_vs_true <- function(geofacet_grid,
                                            pred_to_plot,
                                            title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_point(aes(x = time_id, y = lambda_it, col = "True rate")) + 
    geom_point((aes(x = time_id, y = count_div_pop, col = "sampled count/Eit"))) +
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) + 
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=11),
          plot.title = element_text(size=11),
          strip.text.x = element_text(size = 9))
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue")) # 
  return(plt)
}

select_regions_lin_pred_vs_true_2 <- function(geofacet_grid,
                                              pred_to_plot,
                                              title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true counts
  
  plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
    geom_point((aes(x = time_id, y = lambda_it, col = "True rate per 100"))) +
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) + 
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=11),
          plot.title = element_text(size=11),
          strip.text.x = element_text(size = 9))
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue")) # 
  return(plt)
}


ADM1_grid <- data.frame(
  code = c("15", "8", "6", "3", "5", "9", "4", "13", "10", "14", "16", "7", "11", "12", "2", "1"),
  name = c(" Schleswig-Holstein ", 
           " Mecklenburg-Vorpommern ", 
           " Hamburg", 
           " Berlin ", 
           " Bremen ", 
           " Niedersachsen", 
           " Brandenburg", 
           " Sachsen-Anhalt", 
           " Nordrhein-Westfalen", 
           " Sachsen", 
           " Thuringen", 
           " Hessen ", 
           " Rheinland-Pfalz ", 
           " Saarland", 
           " Bayern", 
           " Baden-Wurttemberg "),
  row = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5),
  col = c(2, 4, 3, 4, 2, 3, 5, 4, 2, 5, 4, 3, 2, 1, 4, 3),
  stringsAsFactors = FALSE)



pred_to_plot <- data.frame(area_id = lambda_sc3.df$area_id,
                           time_id = lambda_sc3.df$time_id,
                           median =  tmp_$summary.fitted.values$'0.5quant', # sort_proper_fitted(tmp_$summary.fitted.values$'0.5quant', n_ADM1, tT),
                           quantile_0.025 = tmp_$summary.fitted.values$'0.025quant', # sort_proper_fitted(tmp_$summary.fitted.values$'0.025quant', n_ADM1, tT),
                           quantile_0.975 =  tmp_$summary.fitted.values$'0.975quant', # sort_proper_fitted(tmp_$summary.fitted.values$'0.975quant', n_ADM1, tT),
                           sampled_counts = lambda_sc3.df$sampled_counts2,
                           lambda_it = lambda_sc3.df$lambda_it,
                           count_div_pop = lambda_sc3.df$sampled_counts2/lambda_sc3.df$E_it2)

select_regions_lin_pred_vs_true(ADM1_grid, pred_to_plot, "")


ul_each <- lapply(tmp_$marginals.fitted.values, 
                  FUN = function(x){
                    return(find_ul_quants_counts_single_pred(x, E_it2))
                  })

### NB: Must sort the proper models

pred_to_plot <- data.frame(area_id = lambda_sc3.df$area_id,
                           time_id = lambda_sc3.df$time_id,
                           median = rep(NA, nrow(lambda_sc3.df)),
                           quantile_0.025 = rep(NA, nrow(lambda_sc3.df)),
                           quantile_0.975 = rep(NA, nrow(lambda_sc3.df)),
                           sampled_counts = lambda_sc3.df$sampled_counts2,
                           lambda_it = lambda_sc3.df$lambda_it * lambda_sc3.df$E_it2,
                           count_div_pop = lambda_sc3.df$sampled_counts2/lambda_sc3.df$E_it2)

for(i in 1:nrow(pred_to_plot)){
  pred_to_plot[i, ]$median = ul_each[[i]]$median
  pred_to_plot[i, ]$quantile_0.025 = ul_each[[i]]$l
  pred_to_plot[i, ]$quantile_0.975 = ul_each[[i]]$u
}

select_regions_lin_pred_vs_true_2(ADM1_grid, pred_to_plot, "")














