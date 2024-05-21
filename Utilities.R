################################################################################
# Functions for calculating MSE, IS,...

count_mse_one_year_one_dataset <- function(sampled_counts_one_year, 
                                           lambda_marginals_one_year,
                                           E_it = 100){
  
  
  ## For each area find the expected predicted count (i.e. point prediction)
  pred_count <- E_it * as.numeric(sapply(lambda_marginals_one_year, 
                                        FUN = function(x){return(mean(inla.rmarginal(2000, 
                                                                                     x)))}))
  
  
  
  ## Calculate the MSE
  mse_one_year <- mean((sampled_counts_one_year - pred_count)**2)
  
  
  
  
  # Return the MSE of that year
  return(mse_one_year)
}

## Testing
#one_year_margs <- marginals[(11 * n_ADM1 + 1):(12 * n_ADM1)]
#count_mse_one_year_one_dataset(lambda.df$sampled_counts[1],one_year_margs,100)




rate_mse_one_year_one_dataset <- function(sampled_rates_one_year, 
                                          lambda_marginals_one_year,
                                          E_it = 100){
  
  pred_rate <- as.numeric(sapply(lambda_marginals_one_year,
                                 FUN = function(x){
                                   return(inla.emarginal(function(x){x}, marginal = x))
                                   })
                          )
  
  ## For each area find the expected predicted count (i.e. point prediction)
  #pred_rate <- as.numeric(sapply(lambda_marginals_one_year, 
  #                                FUN = function(x){return(mean(inla.rmarginal(2000, 
  #                                                                             x))
  #                                                         )}))
  
  ## Calculate the MSE of rate per 100
  mse_one_year <- mean((sampled_rates_one_year * E_it - pred_rate * E_it)**2)
  
  # Return the MSE of that year
  return(mse_one_year)
}



## Testing
#count_mean_mse_one_dataset(lambda.df$sampled_counts[, 1],marginals,n_ADM1)

find_ul_quants_counts_single_pred <- function(lambda_marginal,
                                              E_it){
  # Function to calculate upper (u) and lower (l) quantiles for a single
  # count prediction
  
  ## Sample lambda and scale w. E_it to get an instance of Poisson
  poisson_param_sample <- E_it * inla.rmarginal(5000, lambda_marginal) #[[1]]
  
  ## Sample from Poisson
  count_sample <- sapply(poisson_param_sample, 
                         FUN = function(x){return(rpois(1, x))})
  
  ## Calculate upper and lower quantile (and median)
  u = as.numeric(quantile(count_sample, 0.975)); l = as.numeric(quantile(count_sample, 0.025))
  median = as.numeric(quantile(count_sample, 0.5))
  mean <- as.numeric(mean(count_sample))
  
  return(list(l = l, u = u, median = median,
              post_mean = mean))
  
}

## Test for find_ul_quants_counts_single_pred
#find_ul_quants_counts_single_pred(marginals[12 * n_ADM1 + 1], 100)

#Interval-scores
find_IS_one_obs <- function(l, u, true_value){
  IS_score = (u - l) + 
    2/0.05 * (l - true_value) * (true_value < l) + 
    2/0.05 * (true_value - u) * (true_value > u) 
  return(IS_score)
}



count_IS_one_year_case_study <- function(counts,
                                         marginals,
                                         population){
  
  IS_each_instance = rep(0, length(marginals))
  
  # Iterate over each marginal
  for(i in 1:length(marginals)){
    ul_each_one_year <- find_ul_quants_counts_single_pred(marginals[[i]], 
                                                          population[i])
    
    
    IS_each_instance[i] = find_IS_one_obs(ul_each_one_year$l, ul_each_one_year$u, 
                                          counts[i])
  }
  
  ## Find the average IS this year
  IS_this_year <- mean(IS_each_instance)
  
  return(IS_this_year)
}



count_IS_one_year_one_dataset <- function(sampled_counts_one_year,
                                          lambda_marginals_one_year,
                                          E_it = 100){
  
  ## Find the l and u for each singular instance in a year
  ul_each_one_year <- lapply(lambda_marginals_one_year, 
                             FUN = function(x){
                               return(find_ul_quants_counts_single_pred(x, E_it))
                             })
  
  ## Find the IS for each singular instance in a year
  ### Initialize space for IS 
  IS_each_instance = rep(0, length(lambda_marginals_one_year))
  
  ### Calculate IS each area
  for(i in 1:length(lambda_marginals_one_year)){
    IS_each_instance[i] = find_IS_one_obs(ul_each_one_year[[i]]$l, ul_each_one_year[[i]]$u, 
                                          sampled_counts_one_year[i])
  }
  
  ## Find the average IS this year
  IS_this_year <- mean(IS_each_instance)
  
  return(IS_this_year)
}

## Test of count_IS_one_year_one_dataset
#test2 <- count_IS_one_year_one_dataset(lambda.df$sampled_counts[(10 * n_ADM1 + 1):(11 * n_ADM1), 1], marginals[(10 * n_ADM1 + 1):(11 * n_ADM1)], 100)


find_ul_quants_rate_single_pred <- function(lambda_marginal){
  # Function to calculate upper (u) and lower (l) quantiles for a single
  # rate prediction
  
  ## Sample lambda and scale w. E_it to get an instance of Poisson
  rate_sample <- inla.rmarginal(5000, lambda_marginal) #[[1]]
  
  ## Calculate upper and lower quantile (and median)
  u = as.numeric(quantile(rate_sample, 0.975)); l = as.numeric(quantile(rate_sample, 0.025))
  median = as.numeric(quantile(rate_sample, 0.5))
  
  return(list(l = l, u = u, median = median))
  
}

width_CI_one_year_one_dataset <- function(lambda_marginals_one_year,
                                          E_it = 100){
  
  
  ## Find the l and u for each singular instance in a year
  ul_each_one_year <- lapply(lambda_marginals_one_year, 
                             FUN = function(x){
                               return(find_ul_quants_rate_single_pred(x))
                             })
  
  width_each_instance = rep(0, length(lambda_marginals_one_year))
  
  ### Calculate IS each area of rate per 100
  for(i in 1:length(lambda_marginals_one_year)){
    width_each_instance[i] = (ul_each_one_year[[i]]$u - ul_each_one_year[[i]]$l) * E_it
  }
  
  ## Find the average IS this year
  avg_width <- mean(width_each_instance)
  
  return(avg_width)
  
}





rate_IS_one_year_one_dataset <- function(sampled_rates_one_year,
                                         lambda_marginals_one_year,
                                         E_it = 100){
  
  ## Find the l and u for each singular instance in a year
  ul_each_one_year <- lapply(lambda_marginals_one_year, 
                             FUN = function(x){
                               return(find_ul_quants_rate_single_pred(x))
                             })
  
  ## Find the IS for each singular instance in a year
  ### Initialize space for IS 
  IS_each_instance = rep(0, length(lambda_marginals_one_year))
  
  ### Calculate IS each area of rate per 100
  for(i in 1:length(lambda_marginals_one_year)){
    IS_each_instance[i] = find_IS_one_obs(ul_each_one_year[[i]]$l * E_it, ul_each_one_year[[i]]$u * E_it, 
                                          sampled_rates_one_year[i] * E_it)
  }
  
  ## Find the average IS this year
  IS_this_year <- mean(IS_each_instance)
  
  return(IS_this_year)
}


count_log_s_one_year <- function(counts,
                                 marginals,
                                 population){
  

  log_s_each = rep(0, length(marginals))
  
  # Iterate over each marginal
  for(i in 1:length(marginals)){
    #scaled_marg <- inla.tmarginal(fun = function(x){x * population[i]},
    #                              marginal = marginals[[i]])
    
    #This dont include the Poisson noise...
    
    log_s_each[i] <- inla.dmarginal(counts[i], 
                                    scaled_marg) 
  }
  
  
  
  ## Find the average IS this year
  avg_log_s <- mean(log_s_each)
  
  return(avg_log_s)
  
  
}

get_pred_SD <- function(rate_marginals,
                        population,
                        n, 
                        year){
  one_year_marg <- rate_marginals[(n * (year - 1) + 1):(n * year)]
  one_year_pop <- population[(n * (year - 1) + 1):(n * year)]
  
  #ul_each_one_year <- find_ul_quants_counts_single_pred(marginals[[i]], 
  #                                                      population[i])
  tmp_sd <- rep(NA, n)
  
  
  
  for(i in 1:n){
    poisson_param_sample <- 1E5 * inla.rmarginal(5000, one_year_marg[[i]]) #[[1]]
    
    ## Sample from Poisson
    count_sample <- sapply(poisson_param_sample, 
                           FUN = function(x){return(rpois(1, x))})
    
    tmp_sd[i] <- sd(count_sample)
  }
  
  return(tmp_sd)
}

calc_width_CI_and_count_obs_outside <- function(marginals,
                                                observations, 
                                                population,
                                                n, tT = 25, start_of_predictions = 21){
  ### Initialize data frame to store the avg width of the 95% CI 1,...,5 years ahead
  ### and the number of obs that fall outside the 95% CIs
  tmp.df <- data.frame(width_CI_1_year_ahead = NA, width_CI_2_year_ahead = NA,
                       width_CI_3_year_ahead = NA, width_CI_4_year_ahead = NA,
                       width_CI_5_year_ahead = NA, width_CI_avg = NA,
                       obs_outside_1_year = NA, obs_outside_2_year = NA,
                       obs_outside_3_year = NA, obs_outside_4_year = NA,
                       obs_outside_5_year = NA, obs_outside_tot = NA)
  
  
  #### Take the years predicted on
  years_pred_on <- start_of_predictions:tT
  
  
  #### For each year calculate the width of the avg width of the CI that year
  #### and the number of obs outside the 95% CIs
  for(year in years_pred_on){
    print(year)
    
    ## Extract the predicted marginals for this year
    one_year_margs <- marginals[((year - 1) * n + 1):(year * n)]
    
    ### Extract the observed counts for this year
    one_year_obs <- observations[((year - 1) * n + 1):(year * n)]
    
    ### Extrect the population for this year
    one_year_pop <- population[((year - 1) * n + 1):(year * n)]
    
    ### Calculate the 95% CIs for each area, check if obs outside record
    one_year_widths <- rep(NA, n)
    one_year_count_misses = 0
    for(i in 1:n){
      
      # Get samples of the poisson parameters
      poisson_param_sample <- one_year_pop[i] * inla.rmarginal(5000, one_year_margs[[i]]) #[[1]]
      
      ## Sample from Poisson
      count_sample <- sapply(poisson_param_sample,
                             FUN = function(x){return(rpois(1, x))})
      
      # Calculate upper and lower quantiles to get 95% CI
      u = as.numeric(quantile(count_sample, 0.975)); l = as.numeric(quantile(count_sample, 0.025))
      
      # Get the width of the 95% CI
      one_year_widths[i] = u - l
      
      # Check if observation is outside of 95% CI, if so record it as a miss
      if(one_year_obs[i] < l |
         one_year_obs[i] > u){
        one_year_count_misses = one_year_count_misses + 1
      }
    }
    
    tmp.df[1, year - years_pred_on[1] + 1] =  mean(one_year_widths)
    tmp.df[1, year - years_pred_on[1] + 7] =  one_year_count_misses
    
    
  }
  tmp.df[1, 6] = mean(as.numeric(tmp.df[1, 1:5]))
  tmp.df[1, 12] = sum(as.numeric(tmp.df[1, 7:11]))
  
  return(tmp.df)
}


################################################################################
# CSV-tracking files and functions

initialize_csv_tracker <- function(model_name, 
                                   scenario){
  #Function that creates a data frame in a csv that tracks which files have
  # been analyzed, potential errors and warnings!
  filename = paste("./Simulated_data/", scenario, "/", 
                   model_name, "_", scenario, "_tracker.csv", sep = "")
  
  df_ = data.frame(analyzed = rep(NA, 100), 
                   error = rep(NA, 100),
                   warning = rep(NA, 100))
  
  write.csv(df_, file = filename, row.names = F)
}

get_csv_tracker_filename <- function(model_name, scenario){
  return(paste("./Simulated_data/", scenario, "/", 
               model_name, "_", scenario, "_tracker.csv", sep = ""))
}

get_first_not_yet_analyzed <- function(model_name, scenario){
  df_ = read.csv(get_csv_tracker_filename(model_name = model_name, 
                                          scenario = scenario))
  
  
  for(index in 1:nrow(df_)){
    if(is.na(df_$analyzed[index])){
      if(is.na(df_$error[index])){
        break
      }
    }
  }
  
  return(index)
}



################################################################################
#Plotting functions

################
# Time series of specific areas

timeseries_plt <- function(geofacet_grid,
                          pred_to_plot,
                          title){
  
  
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  if(nrow(geofacet_grid) == 16){
    theme <- theme(axis.title=element_text(size=15),
                   plot.title = element_text(hjust = 0.5, size=15),
                   strip.text.x = element_text(size = 15),
                   legend.position = c(0.15, 1),
                   legend.justification = c("right", "top"),
                   legend.box.just = "right",
                   legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"),
                   legend.text = element_text(size=15),
                   axis.text=element_text(size=13))
  } else{
    theme <- theme(axis.title=element_text(size=15),
                   plot.title = element_text(hjust = 0.5, size=15),
                   strip.text.x = element_text(size = 15),
                   legend.background = element_rect(linetype = 1,
                                                    linewidth = 1,
                                                    colour = "grey"),
                   legend.position = "bottom",
                   legend.text = element_text(size=13),
                   axis.text=element_text(size=13))
  }
  
  plt <- ggplot(data = pred_to_plot, aes(time_id)) + 
    geom_ribbon(aes(x = time_id, ymin = pred_count_quantile_0.025, 
                    ymax = pred_count_quantile_0.975, 
                    col = "A"), 
                fill = "#CC33CC", alpha = 0.15) +
    geom_ribbon(aes(x = time_id, ymin = pred_rate_quantile_0.025, 
                    ymax = pred_rate_quantile_0.975, 
                    col = "B"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_line(aes(x = time_id, y = pred_rate_median, 
                  col = "C")) + 
    geom_point(aes(x = time_id, 
                   y = lambda_it, 
                   col = "D")) + 
    geom_point((aes(x = time_id, 
                    y = sampled_counts, 
                    col = "E"))) + # TeX(r'($Y_{it}/E_{it}$)')
    geom_vline(xintercept = 10.5, linetype = "longdash", 
               color = "darkgrey", linewidth = 0.6) +
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate per 100",
         col = NULL) +
    theme_bw() + 
    theme
  
  plt <- plt + scale_color_manual(values = c("#CC33CC", 
                                             "#F8766D",
                                             "black",
                                             "#00BFC4", 
                                             "blue"),
                                  labels = unname(TeX(c(" Posterior 95% CI: $Y_{it}$",
                                                        " Posterior 95% CI: $\\lambda_{it}E_{it}$",
                                                        " Posterior median $\\lambda_{it}E_{it}$",
                                                        " True: $\\lambda_{it}E_{it}$",
                                                        " True: $Y_{it}$")))) # 
  
  plt <- plt + scale_x_continuous(breaks = c(1, 4, 7, 10, 13))
  
  return(plt)
}












wrapper_timeseries_plt <- function(ADM_grid, 
                                   scenario_name,
                                   dataset_id,
                                   model,
                                   improper = T,
                                   title){
  
  ### Load in simulated data for that scenario
  load(paste("./Simulated_data/", scenario_name, "/", 
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
  lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]
  
  #### Remove unessecary data
  rm(lambda.df)
  
  ## If proper sort the marginals
  print("All right love")
  
  if(!improper){
    model$marginals.fitted.values <- sort_proper_fitted(model$marginals.fitted.values,
                                                        length(unique(lambda_$area_id)),
                                                        tT)
    
    model$summary.fitted.values$'0.5quant' <- sort_proper_fitted(model$summary.fitted.values$'0.5quant',
                                                      length(unique(lambda_$area_id)),
                                                      tT)
    
    model$summary.fitted.values$'0.025quant' <- sort_proper_fitted(model$summary.fitted.values$'0.025quant',
                                                                 length(unique(lambda_$area_id)),
                                                                 tT)
    
    model$summary.fitted.values$'0.975quant' <- sort_proper_fitted(model$summary.fitted.values$'0.975quant',
                                                                 length(unique(lambda_$area_id)),
                                                                 tT)
    
  } 
  
  ## Get the upper, lower, and median quantile for the pred. counts
  ul_each <- lapply(model$marginals.fitted.values, 
                    FUN = function(x){
                      return(find_ul_quants_counts_single_pred(x, 100)) #100 aka E_it
                    })
  
  print("HELL")
  
  
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             pred_rate_median = model$summary.fitted.values$'0.5quant' * lambda_$E_it, 
                             pred_rate_quantile_0.025 = model$summary.fitted.values$'0.025quant' * lambda_$E_it, 
                             pred_rate_quantile_0.975 = model$summary.fitted.values$'0.975quant' * lambda_$E_it, 
                             pred_count_mean = rep(NA, nrow(lambda_)),
                             pred_count_quantile_0.025 = rep(NA, nrow(lambda_)),
                             pred_count_quantile_0.975 = rep(NA, nrow(lambda_)),
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it * lambda_$E_it)
  
  print("Cheeky bugger")
  
  
  for(i in 1:nrow(pred_to_plot)){
    pred_to_plot[i, ]$pred_count_mean = ul_each[[i]]$post_mean
    pred_to_plot[i, ]$pred_count_quantile_0.025 = ul_each[[i]]$l
    pred_to_plot[i, ]$pred_count_quantile_0.975 = ul_each[[i]]$u
  }
  
  timeseries_plt(ADM_grid, pred_to_plot, title)
}






plt_linpredictor_vs_true_rate <- function(geofacet_grid,
                                          pred_to_plot,
                                          title){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  if(nrow(geofacet_grid) == 16){
    theme <- theme(axis.title=element_text(size=12),
                   plot.title = element_text(hjust = 0.5, size=12),
                   strip.text.x = element_text(size = 10),
                   legend.position = c(0.15, 1),
                   legend.justification = c("right", "top"),
                   legend.box.just = "right",
                   legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  } else{
    theme <- theme(axis.title=element_text(size=12),
                   plot.title = element_text(hjust = 0.5, size=12),
                   strip.text.x = element_text(size = 10),
                   legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  }
  
  plt <- ggplot(data = pred_to_plot, aes(time_id, median)) + 
    geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = "Posterior 95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_point(aes(x = time_id, y = lambda_it, col = "True rate")) + 
    geom_point((aes(x = time_id, y = sampled_counts, col = "sampled count"))) + # TeX(r'($Y_{it}/E_{it}$)')
    geom_line(aes(x = time_id, y = median, col = "Posterior median risk")) + 
    geom_vline(xintercept = 10.5, linetype = "longdash", 
               color = "darkgrey", linewidth = 0.6) +
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate per 100",
         col = NULL) +
    theme_bw() + 
    theme
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue"),
                                  labels = unname(TeX(c("Posterior 95% CI of rate per 100",
                                                        "Posterior median rate per 100",
                                                        "Simulated rate per 100: $\\lambda_{it}\\cdot E_{it}$",
                                                        "Simulated count $Y_{it}$")))) # 
  
  plt <- plt + scale_x_continuous(breaks = c(1, 4, 7, 10, 13))
  
  return(plt)
}

plt_linpredictor_vs_true_rate_improper <- function(ADM_grid, lambda_, model, title){
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = model$summary.fitted.values$'0.5quant' * lambda_$E_it, 
                             quantile_0.025 = model$summary.fitted.values$'0.025quant' * lambda_$E_it, 
                             quantile_0.975 = model$summary.fitted.values$'0.975quant' * lambda_$E_it, 
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it * lambda_$E_it)
  
  plt_linpredictor_vs_true_rate(ADM_grid, pred_to_plot, title)
  
}

plt_linpredictor_vs_true_rate_proper <- function(ADM_grid, lambda_, model, title){
  ### NB: Must sort the proper ones
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = sort_proper_fitted(model$summary.fitted.values$'0.5quant',
                                                         length(unique(lambda_$area_id)),
                                                         tT) * lambda_$E_it, # model$summary.fitted.values$'0.5quant',
                             quantile_0.025 = sort_proper_fitted(model$summary.fitted.values$'0.025quant',
                                                                 length(unique(lambda_$area_id)),
                                                                 tT) * lambda_$E_it,# model$summary.fitted.values$'0.025quant',
                             quantile_0.975 = sort_proper_fitted(model$summary.fitted.values$'0.975quant',
                                                                 length(unique(lambda_$area_id)),
                                                                 tT) * lambda_$E_it, #model$summary.fitted.values$'0.975quant',
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it * lambda_$E_it)
  
  plt_linpredictor_vs_true_rate(ADM_grid, pred_to_plot, title)
}


wrapper_plt_linpredictor_vs_true_rate <- function(ADM_grid, 
                                                  scenario_name,
                                                  dataset_id,
                                                  model,
                                                  improper = T,
                                                  title){
  
  ### Load in simulated data for that scenario
  load(paste("./Simulated_data/", scenario_name, "/", 
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
  lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]
  
  #### Remove unessecary data
  rm(lambda.df)
  
  ## Actually plot it
  if(improper){
    plt_linpredictor_vs_true_rate_improper(ADM_grid, lambda_, model, title)
  } else{
    plt_linpredictor_vs_true_rate_proper(ADM_grid, lambda_, model, title)
  }
}




plt_predcount_vs_true_count <- function(geofacet_grid,
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
    geom_line(aes(x = time_id, y = post_mean, col = "Posterior mean risk")) +
    geom_vline(xintercept = 10.5, linetype = "longdash", 
               color = "darkgrey", linewidth = 0.6) +
    facet_geo(~ area_id, grid = geofacet_grid, label = "name") + 
    labs(title = title,
         x = "Year",
         y = "Rate",
         col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=12),
          plot.title = element_text(hjust = 0.5, size=12),
          strip.text.x = element_text(size = 10),
          legend.position = c(0.15, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))
  
  plt <- plt + scale_color_manual(values = c("#F8766D", "black", "#00BFC4", "blue"),
                                  labels = unname(TeX(c("Posterior 95% CI",
                                                        "Posterior mean $\\hat{Y}_{it}$",
                                                        "True $Y_{it}$",
                                                        "Expected true $Y_{it}$")))) # 
  return(plt)
}

wrapper_plt_predcount_vs_true_count <- function(ADM1_grid, 
                                                scenario_name,
                                                dataset_id,
                                                model,
                                                improper = T,
                                                title){
  
  ### Load in simulated data for that scenario
  load(paste("./Simulated_data/", scenario_name, "/", 
             scenario_name, "_data.RData", sep = ""))
  
  lambda_ <- lambda.df[, c("area_id", "time_id", "E_it", "space.time")]
  lambda_$sampled_counts <- lambda.df$sampled_counts[, dataset_id]
  lambda_$lambda_it <- lambda.df$lambda_it[, dataset_id]
  
  #### Remove unessecary data
  rm(lambda.df)
  
  ## If proper sort the marginals
  if(!improper){
    model$marginals.fitted.values <- sort_proper_fitted(model$marginals.fitted.values,
                                                        length(unique(lambda_$area_id)),
                                                        tT)
  } 
  
  ## Get the upper, lower, and median quantile for the pred. counts
  ul_each <- lapply(model$marginals.fitted.values, 
                    FUN = function(x){
                      return(find_ul_quants_counts_single_pred(x, 100)) #100 aka E_it
                    })
  
  pred_to_plot <- data.frame(area_id = lambda_$area_id,
                             time_id = lambda_$time_id,
                             median = rep(NA, nrow(lambda_)),
                             post_mean = rep(NA, nrow(lambda_)),
                             quantile_0.025 = rep(NA, nrow(lambda_)),
                             quantile_0.975 = rep(NA, nrow(lambda_)),
                             sampled_counts = lambda_$sampled_counts,
                             lambda_it = lambda_$lambda_it * lambda_$E_it,
                             count_div_pop = lambda_$sampled_counts/lambda_$E_it)
  
  for(i in 1:nrow(pred_to_plot)){
    pred_to_plot[i, ]$median = ul_each[[i]]$median
    pred_to_plot[i, ]$post_mean = ul_each[[i]]$post_mean
    pred_to_plot[i, ]$quantile_0.025 = ul_each[[i]]$l
    pred_to_plot[i, ]$quantile_0.975 = ul_each[[i]]$u
  }
  
  plt_predcount_vs_true_count(ADM1_grid, pred_to_plot, title)
}


################
# Every other plotting function
ridgeplot_counts <- function(model_names, scenario_name,
                             one_2_3_or_total, xlab,
                             xlim){
  
  model_name = model_names[1]
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_counts
  
  to_plot_ridge.df <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                                 IS_one_year_ahead = tmp_$IS_1_year_ahead,
                                 IS_two_year_ahead = tmp_$IS_2_year_ahead,
                                 IS_three_year_ahead = tmp_$IS_3_year_ahead,
                                 IS_tot = tmp_$total_IS,
                                 MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                                 MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                                 MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                                 MSE_tot = tmp_$total_mse)
  
  
  for(model_name in model_names[2:length(model_names)]){
    ### Load in Model choice data
    load(paste("./results/model_choice/model_choice_", 
               model_name, "_", 
               scenario_name, ".RData",
               sep = ""))
    
    tmp_ <- model_choice_for_counts
    
    tmp2_ <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                        IS_one_year_ahead = tmp_$IS_1_year_ahead,
                        IS_two_year_ahead = tmp_$IS_2_year_ahead,
                        IS_three_year_ahead = tmp_$IS_3_year_ahead,
                        IS_tot = tmp_$total_IS,
                        MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                        MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                        MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                        MSE_tot = tmp_$total_mse)
    
    to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
    
  }
  if(one_2_3_or_total == 1){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_one_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab("Model") + 
        xlim(xlim[1], xlim[2])
      
    )
  } else if(one_2_3_or_total == 2){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_two_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.text.x = element_text(size = 10),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) + 
        xlim(xlim[1], xlim[2])
    )
    
  } else if(one_2_3_or_total == 3){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_three_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) + 
        xlim(xlim[1], xlim[2])
    )
  }
}




ridgeplot_rates <- function(model_names, scenario_name,
                            one_2_3_or_total, xlab,
                            xlim){
  
  model_name = model_names[1]
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_rates
  
  to_plot_ridge.df <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                                 IS_one_year_ahead = tmp_$IS_1_year_ahead,
                                 IS_two_year_ahead = tmp_$IS_2_year_ahead,
                                 IS_three_year_ahead = tmp_$IS_3_year_ahead,
                                 IS_tot = tmp_$total_IS,
                                 MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                                 MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                                 MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                                 MSE_tot = tmp_$total_mse)
  
  
  for(model_name in model_names[2:length(model_names)]){
    ### Load in Model choice data
    load(paste("./results/model_choice/model_choice_", 
               model_name, "_", 
               scenario_name, ".RData",
               sep = ""))
    
    tmp_ <- model_choice_for_rates
    
    tmp2_ <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                        IS_one_year_ahead = tmp_$IS_1_year_ahead,
                        IS_two_year_ahead = tmp_$IS_2_year_ahead,
                        IS_three_year_ahead = tmp_$IS_3_year_ahead,
                        IS_tot = tmp_$total_IS,
                        MSE_one_year_ahead = tmp_$mse_1_year_ahead,
                        MSE_two_year_ahead = tmp_$mse_2_year_ahead,
                        MSE_three_year_ahead = tmp_$mse_3_year_ahead,
                        MSE_tot = tmp_$total_mse)
    
    to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
    
  }
  
  if(one_2_3_or_total == 1){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_one_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.text.x = element_text(size = 10),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab("Model") +
        xlim(xlim[1], xlim[2])
      
    )
  } else if(one_2_3_or_total == 2){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_two_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) +
        xlim(xlim[1], xlim[2])
    )
    
  } else if(one_2_3_or_total == 3){
    return(
      ggplot(to_plot_ridge.df, aes(x = IS_three_year_ahead, 
                                   y = model_name, 
                                   fill = model_name)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none",
              axis.text.y=element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
              axis.title=element_text(size=13)) + 
        xlab(xlab) + 
        ylab(NULL) +
        xlim(xlim[1], xlim[2])
    )
  }
}




ridgeplot_counts_and_rates_all_years <- function(scenario_name, xlim_counts, xlim_rates){
  
  plt1 <- ridgeplot_counts(model_names = model_names, scenario_name = scenario_name, 1, "IS count: 1 year ahead", xlim_counts)
  plt2 <- ridgeplot_counts(model_names = model_names, scenario_name = scenario_name, 2, "IS count: 2 years ahead", xlim_counts)
  plt3 <- ridgeplot_counts(model_names = model_names, scenario_name = scenario_name, 3, "IS count: 3 years ahead", xlim_counts)
  
  
  plt4 <- ridgeplot_rates(model_names = model_names, scenario_name = scenario_name, 1, "IS rate: 1 year ahead", xlim_rates)
  plt5 <- ridgeplot_rates(model_names = model_names, scenario_name = scenario_name, 2, "IS rate: 2 years ahead", xlim_rates)
  plt6 <- ridgeplot_rates(model_names = model_names, scenario_name = scenario_name, 3, "IS rate: 3 years ahead", xlim_rates)
  
  plt <- ggarrange(plt1, plt2, plt3, 
                   plt4, plt5, plt6,
                   widths = c(1.35, 1, 1),
                   ncol = 3, nrow = 2, common.legend = F)
  return(plt)
}



# Make ridgeplot w. MSE and IS for rates
ridgeplot_mse_is_rates <- function(to_plot.df, 
                                   value_to_plot, 
                                   one_2_3_or_total, 
                                   IS_or_MSE = "IS", 
                                   xlab, 
                                   xlim){
  
  
  to_plot.df$value_to_plot = value_to_plot
  
  if(one_2_3_or_total == 1 & IS_or_MSE == "IS"){
    ylab = "Model"
    axis.text.y = element_text(size = 15)
    title = NULL
  } else if(one_2_3_or_total == 2 & IS_or_MSE == "IS"){
    ylab = NULL
    axis.text.y = element_blank()
    title = NULL
  } else if(one_2_3_or_total == 3 & IS_or_MSE == "IS"){
    ylab = NULL
    axis.text.y = element_blank()
    title = NULL
  } else if(one_2_3_or_total == 4 & IS_or_MSE == "IS"){
    ylab = NULL
    axis.text.y = element_blank()
    title = NULL
    
  } else if(one_2_3_or_total == 1 & IS_or_MSE == "MSE"){
    ylab = "Model"
    axis.text.y = element_text(size = 15)
    title = "1 year ahead"
  } else if(one_2_3_or_total == 2 & IS_or_MSE == "MSE"){
    ylab = NULL
    axis.text.y = element_blank()
    title = "2 years ahead"
  } else if(one_2_3_or_total == 3 & IS_or_MSE == "MSE"){
    ylab = NULL
    axis.text.y = element_blank()
    title = "3 years ahead"
  } else if(one_2_3_or_total == 4 & IS_or_MSE == "MSE"){
    ylab = NULL
    axis.text.y = element_blank()
    title = "Average all years"
  }
  
  return(ggplot(data = to_plot.df, 
                aes(x = value_to_plot, 
                    y = model_name, 
                    fill = model_name)) + 
           ggtitle(title) +
           geom_density_ridges() +
           theme_ridges() + 
           theme(legend.position = "none",
                 axis.text.y = axis.text.y,
                 axis.title.y = element_text(hjust = 0.5, vjust = 0.5,
                                             size = 15),
                 axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
                 axis.text.x = element_text(size = 15),
                 axis.title = element_text(size=15),
                 plot.title = element_text(size = 15, hjust = 0.5,
                                           face = "plain")) + 
           xlab(xlab) + 
           ylab(ylab) +
           xlim(xlim[1], xlim[2]))
}

ridgeplot_mse_is_rates_all_years <- function(model_names,
                                             scenario_name, 
                                             xlim_mse, 
                                             xlim_is,
                                             presentation = F){
  
  model_name = model_names[1]
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_rates
  
  to_plot_ridge.df <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                                 IS_1_year_ahead = tmp_$IS_1_year_ahead,
                                 IS_2_year_ahead = tmp_$IS_2_year_ahead,
                                 IS_3_year_ahead = tmp_$IS_3_year_ahead,
                                 IS_tot = tmp_$total_IS,
                                 mse_1_year_ahead = tmp_$mse_1_year_ahead,
                                 mse_2_year_ahead = tmp_$mse_2_year_ahead,
                                 mse_3_year_ahead = tmp_$mse_3_year_ahead,
                                 MSE_tot = tmp_$total_mse)
  
  
  
  
  for(model_name in model_names[2:length(model_names)]){
    ### Load in Model choice data
    load(paste("./results/model_choice/model_choice_", 
               model_name, "_", 
               scenario_name, ".RData",
               sep = ""))
    
    tmp_ <- model_choice_for_rates
    
    tmp2_ <- data.frame(model_name = rep(model_name, nrow(tmp_)),
                        IS_1_year_ahead = tmp_$IS_1_year_ahead,
                        IS_2_year_ahead = tmp_$IS_2_year_ahead,
                        IS_3_year_ahead = tmp_$IS_3_year_ahead,
                        IS_tot = tmp_$total_IS,
                        mse_1_year_ahead = tmp_$mse_1_year_ahead,
                        mse_2_year_ahead = tmp_$mse_2_year_ahead,
                        mse_3_year_ahead = tmp_$mse_3_year_ahead,
                        MSE_tot = tmp_$total_mse)
    
    
    to_plot_ridge.df = rbind(to_plot_ridge.df, tmp2_)
  }
  
  widths = c(1.4, 1, 1, 1)
  
  if(presentation){
    widths = c(1.5, 1, 1)
    
    to_plot_ridge.df[to_plot_ridge.df$model_name == "Improper1_noInt", ]$model_name = "Improper1 noInt"
    to_plot_ridge.df[to_plot_ridge.df$model_name == "Improper1_typeIV", ]$model_name = "Improper1 typeIV"
    to_plot_ridge.df[to_plot_ridge.df$model_name == "proper2_onlyInt", ]$model_name = "Proper2 onlyInt"
    to_plot_ridge.df[to_plot_ridge.df$model_name == "proper2_propInt_Improp_temporal", ]$model_name = "Proper2 & RW1"
  }
  
  plt1 <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                 value_to_plot = to_plot_ridge.df$mse_1_year_ahead,
                                 one_2_3_or_total = 1, 
                                 IS_or_MSE = "MSE",
                                 xlab = TeX(r'(MSE$\left(\widehat{\lambda_{11}E_{11}}\right)$)'), xlim_mse)
  plt2 <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                 value_to_plot = to_plot_ridge.df$mse_2_year_ahead,
                                 one_2_3_or_total = 2, 
                                 IS_or_MSE = "MSE",
                                 xlab = TeX(r'(MSE$\left(\widehat{\lambda_{12}E_{12}}\right)$)'), xlim_mse)
  
  plt3 <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                 value_to_plot = to_plot_ridge.df$mse_3_year_ahead,
                                 one_2_3_or_total = 3, 
                                 IS_or_MSE = "MSE",
                                 xlab = TeX(r'(MSE$\left(\widehat{\lambda_{13}E_{13}}\right)$)'), xlim_mse)
  
  plt_mse_tot <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                        value_to_plot = to_plot_ridge.df$MSE_tot,
                                        one_2_3_or_total = 4, 
                                        IS_or_MSE = "MSE",
                                        xlab = TeX(r'(MSE$\left(\widehat{\lambda_{11,12,13}E_{11,12,13}}\right)$)'), xlim_mse)
  
  
  plt4 <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                 value_to_plot = to_plot_ridge.df$IS_1_year_ahead,
                                 one_2_3_or_total = 1, 
                                 IS_or_MSE = "IS",
                                 xlab = TeX(r'(IS$\left(\widehat{\lambda_{11}E_{11}}\right)$)'), xlim_is)
  
  plt5 <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                 value_to_plot = to_plot_ridge.df$IS_2_year_ahead,
                                 one_2_3_or_total = 2,
                                 IS_or_MSE = "IS",
                                 xlab = TeX(r'(IS$\left(\widehat{\lambda_{12}E_{12}}\right)$)'), xlim_is)
  
  
  plt6 <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                 value_to_plot = to_plot_ridge.df$IS_3_year_ahead,
                                 one_2_3_or_total = 3,
                                 IS_or_MSE = "IS",
                                 xlab = TeX(r'(IS$\left(\widehat{\lambda_{13}E_{13}}\right)$)'), xlim_is)
  
  plt_is_tot <- ridgeplot_mse_is_rates(to_plot.df = to_plot_ridge.df,
                                       value_to_plot = to_plot_ridge.df$IS_tot,
                                       one_2_3_or_total = 4, 
                                       IS_or_MSE = "IS",
                                       xlab = TeX(r'(IS$\left(\widehat{\lambda_{11,12,13}E_{11,12,13}}\right)$)'), xlim_is)
  
  
  plt <- ggarrange(plt1, plt2, plt3, plt_mse_tot,
                   plt4, plt5, plt6, plt_is_tot,
                   widths = widths,
                   ncol = 4, nrow = 2, common.legend = F)
  
  
  
  return(plt)
}


plt_cont_risk_one_time_interval <- function(dir, risk_surface_filename,
                                            time_interval, t_axis,
                                            admin_map,
                                            polygon_grid2,
                                            dataset_id){
  
  ## Load in a specific scenario
  filename = paste(dir, risk_surface_filename, sep = "")
  load(filename)
  col_names <- colnames(risk_surface.list)[colnames(risk_surface.list) != "values"]
  tmp_ = risk_surface.list[, col_names]
  tmp_$values = risk_surface.list$values[, dataset_id]
  rm(risk_surface.list)
  
  ## get the times
  t_axis_indices = which(time_interval - 1 < t_axis & 
                           t_axis <= time_interval)
  
  
  if(length(t_axis_indices) == 3){
    ### Plot the continuous true risk surface for times within time_interval
    plt_1 <- heatmap_points(tmp_,
                            polygon_grid2,
                            admin_map,
                            t_axis[t_axis_indices[1]],
                            title = NULL)
    
    plt_2 <- heatmap_points(tmp_,
                            polygon_grid2,
                            admin_map,
                            t_axis[t_axis_indices[2]],
                            title = NULL)
    
    plt_3 <- heatmap_points(tmp_,
                            polygon_grid2,
                            admin_map,
                            t_axis[t_axis_indices[3]],
                            title = NULL)
    
    
    
    ggarrange(plt_1, plt_2, plt_3,
              ncol = 3, nrow = 1,
              common.legend = T, legend = "right")
    
  } else {
    print("Unsuited")
  }
}



plot_temporal_trend_data_one_data_set <- function(scenario_name, data_set_id, title,
                                                  xlab, ylab){
  #Function to plot the trend of the simulated data, one data set
  
  # Load in the simulated data
  load(paste("./Simulated_data/", scenario_name,"/", scenario_name, "_data.RData",
             sep = ""))
  
  # Extract only the simulated data of data_set_id
  sim_data = lambda.df[, c("space.time", "area_id", "time_id", "E_it")]
  sim_data$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  sim_data$mu = lambda.df$mu[, data_set_id]
  
  # Aggregate over the years
  aggr <- aggregate(sim_data, by = list(time_id = sim_data$time_id), FUN = mean)
  aggr <- aggr[, 2:ncol(aggr)]
  
  aggr1 <- aggr[, c("time_id", "E_it")]
  aggr1$aggregated_value = aggr$sampled_counts
  aggr1$type = rep("Aggregated sampled counts", nrow(aggr))
  aggr2 <- aggr[, c("time_id", "E_it")]
  aggr2$aggregated_value = aggr$mu
  aggr2$type = rep("Aggregated expected count", nrow(aggr))
  
  #Combine the two
  aggr3 <- rbind(aggr1, aggr2)
  
  # Plot over the years
  return(
    ggplot(data = aggr3) + 
      ggtitle(title) + 
      geom_line(aes(x = time_id, y = aggregated_value, col = type, linetype = type)) + #linetype = type
      scale_linetype_manual(values = c("longdash", "solid")) +
      scale_color_manual(values = c("darkred", "blue")) +
      ylim(8, 14) + 
      theme_bw() + 
      theme(axis.title=element_text(size=11),
            plot.title = element_text(size=11),
            strip.text.x = element_text(size = 9)) +
      labs(col = NULL) + 
      guides(linetype = "none") +
      xlab(xlab) + 
      ylab(ylab) 
  )
}


plot_temporal_trend_data_all_data_set <- function(scenario_name, title){
  # Function to plot the simulated average trend of the data across all the simulations for
  # that scenario
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


plot_fitted_temporal_trends <- function(model_on_first_level_20_knots,
                                        model_on_first_level_10_knots,
                                        model_on_second_level_20_knots,
                                        model_on_second_level_10_knots,
                                        tT){
  time = 1:tT
  
  #Scale
  scaling_factor1 <- sqrt(model_on_first_level_20_knots$summary.hyperpar$mean[2]) * 1/sqrt(model_on_first_level_20_knots$summary.hyperpar$mean[1])
  
  tmp1 <- data.frame(time = time)
  tmp1$lower_quant <- model_on_first_level_20_knots$summary.random$time_id[(tT + 1):(2 * tT), 4] * scaling_factor1
  tmp1$median <- model_on_first_level_20_knots$summary.random$time_id[(tT + 1):(2 * tT), 5] * scaling_factor1
  tmp1$upper_quant <- model_on_first_level_20_knots$summary.random$time_id[(tT + 1):(2 * tT), 6] * scaling_factor1
  
  plt1 <- ggplot(data = tmp1, aes(time, median)) + ggtitle("First-level 20 knots") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line() + 
    geom_line(data = tmp1, aes(time, lower_quant), linetype = "dashed") + 
    geom_line(data = tmp1, aes(time, upper_quant), linetype = "dashed") + 
    xlab("t") + ylab(expression(alpha[t]))
  
  scaling_factor2 <- sqrt(model_on_first_level_10_knots$summary.hyperpar$mean[2]) * 1/sqrt(model_on_first_level_10_knots$summary.hyperpar$mean[1])
  
  tmp2 <- data.frame(time = time)
  tmp2$lower_quant <- model_on_first_level_10_knots$summary.random$time_id[(tT + 1):(2 * tT), 4] * scaling_factor2
  tmp2$median <- model_on_first_level_10_knots$summary.random$time_id[(tT + 1):(2 * tT), 5] * scaling_factor2
  tmp2$upper_quant <- model_on_first_level_10_knots$summary.random$time_id[(tT + 1):(2 * tT), 6] * scaling_factor2
  
  plt2 <- ggplot(data = tmp2, aes(time, median)) + ggtitle("First-level 10 knots") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line() + 
    geom_line(data = tmp2, aes(time, lower_quant), linetype = "dashed") + 
    geom_line(data = tmp2, aes(time, upper_quant), linetype = "dashed") + 
    xlab("t") + ylab(expression(alpha[t]))
  
  scaling_factor3 <- sqrt(model_on_second_level_20_knots$summary.hyperpar$mean[2]) * 1/sqrt(model_on_second_level_20_knots$summary.hyperpar$mean[1])
  
  tmp3 <- data.frame(time = time)
  tmp3$lower_quant <- model_on_second_level_20_knots$summary.random$time_id[(tT + 1):(2 * tT), 4] * scaling_factor3
  tmp3$median <- model_on_second_level_20_knots$summary.random$time_id[(tT + 1):(2 * tT), 5] * scaling_factor3
  tmp3$upper_quant <- model_on_second_level_20_knots$summary.random$time_id[(tT + 1):(2 * tT), 6] * scaling_factor3
  
  plt3 <- ggplot(data = tmp3, aes(time, median)) + ggtitle("Second-level 20 knots") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line() + 
    geom_line(data = tmp3, aes(time, lower_quant), linetype = "dashed") + 
    geom_line(data = tmp3, aes(time, upper_quant), linetype = "dashed") + 
    xlab("t") + ylab(expression(alpha[t]))
  
  scaling_factor4 <- sqrt(model_on_second_level_10_knots$summary.hyperpar$mean[2]) * 1/sqrt(model_on_second_level_10_knots$summary.hyperpar$mean[1])
  
  tmp4 <- data.frame(time = time)
  tmp4$lower_quant <- model_on_second_level_10_knots$summary.random$time_id[(tT + 1):(2 * tT), 4] * scaling_factor4
  tmp4$median <- model_on_second_level_10_knots$summary.random$time_id[(tT + 1):(2 * tT), 5] * scaling_factor4
  tmp4$upper_quant <- model_on_second_level_10_knots$summary.random$time_id[(tT + 1):(2 * tT), 6] * scaling_factor4
  
  plt4 <- ggplot(data = tmp4, aes(time, median)) + ggtitle("Second-level 10 knots") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line() + 
    geom_line(data = tmp4, aes(time, lower_quant), linetype = "dashed") + 
    geom_line(data = tmp4, aes(time, upper_quant), linetype = "dashed") + 
    xlab("t") + ylab(expression(alpha[t]))
  
  ggarrange(plt1, plt2, 
            plt3, plt4,
            ncol = 2, nrow = 2)
  
}



heatmap_areas <- function(map_w_values,
                          value,
                          scale_col = NULL,
                          scale = NULL,
                          hardcoded_bins = NULL,
                          title = NULL){
  
  map_w_values$to_plot = value
  
  if(is.null(scale_col)){
    scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  }
  if(is.null(scale)){
    scale = scale_col[c(3, 8, 12, 15, 19, 23, 26, 30)]
  }
  if(is.null(hardcoded_bins)){
    ggplot(data = map_w_values) +  
      geom_sf(aes(fill = to_plot), 
              alpha = 1,
              color="black") + ggtitle(title) + 
      theme(plot.title = element_text(size = 19, hjust = 0.5),
            axis.title.x = element_blank(), #Remove axis and background grid
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            plot.margin =  unit(c(0, 0, 0, 0), "inches"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.text = element_text(size = 15),
            panel.spacing = unit(1, 'lines')) +
      guides(fill=guide_legend(title=NULL,
                               reverse = TRUE,
                               label.position = "right",
                               drop = F)) + #Remove colorbar title
      binned_scale( #Scaling the color
        aesthetics = "fill",
        scale_name = "gradientn",
        palette = function(x) c(scale),
        labels = function(x){x},
        guide = "colorscale") 
    #+ scale_fill_manual(drop = F)
  } else {
    ggplot(data = map_w_values) +  
      geom_sf(aes(fill = to_plot), 
              alpha = 1,
              color="black") + ggtitle(title) + 
      theme(plot.title = element_text(size = 19, hjust = 0.5,
                                      vjust = -0.1),
            axis.title.x = element_blank(), #Remove axis and background grid
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            plot.margin =  unit(c(0, 0, 0, 0), "inches"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.text = element_text(size = 15),
            panel.spacing = unit(1, 'lines')) +
      guides(fill=guide_legend(title=NULL,
                               reverse = TRUE, 
                               label.position = "right",
                               ncol = 1,
                               drop = F)) + #Remove colorbar title
      binned_scale( #Scaling the color
        aesthetics = "fill",
        scale_name = "gradientn",
        palette = function(x) c(scale),
        labels = function(x){x},
        breaks = hardcoded_bins,
        guide = "colorscale") 
    #+ scale_fill_manual(drop = F)
    }
}

heatmap_points <- function(risk_surface.list, 
                           polygons,
                           admin_map,
                           t,
                           title,
                           legends.title = NULL){
  
  ## join risk_surface.list to polygon_grid2 based on polygon_id in order to plot
  tmp_ = data.frame(values = risk_surface.list$values, 
                    t = risk_surface.list$t,
                    polygon_id = risk_surface.list$polygon_id)
  
  tmp2_ = data.frame(polygon_id = polygons$polygon_id,
                     geometry = polygons$x)
  
  tmp3_ = merge(tmp_, tmp2_)
  
  tmp_ = st_set_geometry(tmp3_[, c("t", "values", "polygon_id")], 
                         tmp3_$geometry)
  
  
  # Make it so that each heatmap is plotted on similar color scale 
  scale_col = heat.colors(50, rev=TRUE)          #Divide color gradient into 30 
  scale = scale_col[seq(3, 50, length.out = 15)] #Select color scale to be more red
  risk.min = min(tmp_$values); risk.max = max(tmp_$values) 
  hardcoded_bins =  round(seq(risk.min, risk.max, length.out = 15), 4)
    
  #Extract the values for time = t
  tmp2_ = tmp_[tmp_$t == t, ]
    
  p <- ggplot(data = tmp2_) +  
    
    geom_sf(aes(fill = values), 
            alpha = 1,
            color = NA) + ggtitle(title) + #"lightgrey"
    theme(plot.title = element_text(size = 15, hjust = 0.5), 
          axis.title.x = element_blank(), #Remove axis and background grid
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.margin =  unit(c(0, 0, 0, 0), "inches"),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
          panel.spacing = unit(1, 'lines')) +
    guides(fill=guide_legend(title = legends.title, 
                             reverse = TRUE, 
                             label.position = "right")) + #Remove colorbar title
    binned_scale( #Scaling the color
      aesthetics = "fill",
      scale_name = "gradientn",
      palette = function(x) c(scale),
      labels = function(x){x},
      breaks = hardcoded_bins,
      guide = "colorscale") + 
    geom_sf(data = admin_map,
            aes(), 
            alpha = 0,
            color="black")
    
  return(p)
}



region_time_series <- function(true_risk,
                               model,
                               region,
                               n,
                               tT){
  years <- 1:tT
  values.df <- data.frame(years = years, 
                          true_risk = true_risk$lambda_it[seq(region, n * tT, by = n)],
                          count_div_eit = true_risk$sampled_counts[seq(region, n * tT, by = n)] / 
                                                true_risk$E_it[seq(region, n * tT, by = n)],
                          fitted_rate = model$summary.fitted.values[seq(region,n*tT,by=n), 4],
                          lower_quant = model$summary.fitted.values[seq(region,n*tT,by=n), 3],
                          upper_quant = model$summary.fitted.values[seq(region,n*tT,by=n), 5])
  
  plt <- ggplot(data = values.df, aes(x = years)) + ggtitle(paste("region: ", region)) +
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "#F8766D", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Posterior median risk")) +
    geom_point(aes(x = years, y = true_risk, col = "True risk")) + 
    geom_point((aes(x = years, y = count_div_eit, col = "sampled count/Eit"))) +
    xlab("Time") + ylab("") + 
    labs(col = NULL) +
    theme_bw() + 
    theme(axis.title=element_text(size=9),
          plot.title = element_text(size=10))
  
  plt <- plt + scale_color_manual(values=c("#F8766D", "black", "#00BFC4", "blue"))
  return(plt)
}

select_regions_lin_pred_vs_true <- function(true_risk,
                                            model,
                                            regions,
                                            n, 
                                            tT){
  # Function that plots for select regions the fitted linear predictor of
  # the provided model along w. corresponding 95% CI's
  # against the true risk
  
  plt1 <- region_time_series(true_risk, model, regions[1], n, tT)
  plt2 <- region_time_series(true_risk, model, regions[2], n, tT)
  plt3 <- region_time_series(true_risk, model, regions[3], n, tT)
  plt4 <- region_time_series(true_risk, model, regions[4], n, tT)
  plt5 <- region_time_series(true_risk, model, regions[5], n, tT)
  plt6 <- region_time_series(true_risk, model, regions[6], n, tT)
  plt7 <- region_time_series(true_risk, model, regions[7], n, tT)
  plt8 <- region_time_series(true_risk, model, regions[8], n, tT)
  
  
  ggarrange(plt1, plt2, plt3, 
            plt4, plt5, plt6,
            plt7, plt8,
            ncol = 2, nrow = 4, 
            common.legend = TRUE, legend = "top")
  
  
}





plots_for_GIF <- function(risk_surface.list, 
                          polygons, 
                          t_axis,
                          filename_base){
  ## join risk_surface.list to polygon_grid2 based on polygon_id in order to plot
  tmp_ = data.frame(values = risk_surface.list$values, 
                    t = risk_surface.list$t,
                    polygon_id = risk_surface.list$polygon_id)
  
  tmp2_ = data.frame(polygon_id = polygons$polygon_id,
                     geometry = polygons$x)
  
  tmp3_ = merge(tmp_, tmp2_)
  
  tmp_ = st_set_geometry(tmp3_[, c("t", "values", "polygon_id")], 
                         tmp3_$geometry)
  
  
  # Make it so that each heatmap is plotted on similar color scale 
  scale_col = heat.colors(30, rev=TRUE)          #Divide color gradient into 30 
  scale = scale_col[seq(3, 30, length.out = 15)] #Select color scale to be more red
  risk.min = min(tmp_$values); risk.max = max(tmp_$values) 
  hardcoded_bins =  round(seq(risk.min, risk.max, length.out = 15), 6)
  files = c(); i = 1
  for(t in t_axis){
    #if(t == t_axis[1]){files = c(); i = 1}
    
    #Extract the values for time = t
    tmp2_ = tmp_[tmp_$t == t, ]
    
    #Create a filename (and save filename) for where current time risk-surface is stored
    filename = paste(filename_base, toString(i), ".png", sep = "")
    print(filename); files[i] = filename
    
    p <- ggplot(data = tmp2_) +  
      geom_sf(aes(fill = values), 
              alpha = 1,
              color = NA) + ggtitle(round(t, 1)) + 
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_blank(), #Remove axis and background grid
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            plot.margin =  unit(c(0, 0, 0, 0), "inches"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            panel.spacing = unit(1, 'lines')) +
      guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
      binned_scale( #Scaling the color
        aesthetics = "fill",
        scale_name = "gradientn",
        palette = function(x) c(scale),
        labels = function(x){x},
        breaks = hardcoded_bins,
        guide = "colorscale")
    
    ## Save the plot to filename
    ggsave(filename = filename, plot=p, 
           width=4,height=4,units="in",scale=1)
    
    #update counter i
    i = i + 1
  }
  return(files)
}


################################################################################
#GIF functions

make_GIF <- function(dir, gif_name){
  files <- list.files(path = dir, pattern = "*.png")
  for(i in 1:length(files)){
    if(paste(toString(i), ".png", sep = "") %in% files){
      if(i == 1){
        img = c(image_read(paste(dir, 
                                 paste(toString(i), ".png", sep = ""), 
                                 sep = "")))
      } else{
        img = c(img, image_read(paste(dir, 
                                      paste(toString(i), ".png", sep = ""), 
                                      sep = "")))
      }
    } else{
      print(paste("No such file:", paste(toString(i), ".png", sep = "")))
    }
  }
  
  image_append(image_scale(img, "x200"))
  
  my.animation <- image_animate(image_scale(img, 
                                            "400x400"),
                                fps = 1,
                                dispose = "previous")
  print("Writing")
  image_write(my.animation, gif_name)
  
}


################################################################################
# Functions used for sampling

#Define row-vise Kronecker product
row_wise_Kronecker <- function(X1, X2){
  #Function that returns the row-wise Kronecker product!
  one_k1 <- matrix(1, 1, ncol(X1)); one_k2 <- matrix(1, 1, ncol(X2))
  return((X2 %x% one_k1) * (one_k2 %x% X1))
}


simulate_risk_surface <- function(seed,
                                  Bst, 
                                  Bs = NULL,
                                  Sigma_st, 
                                  kt,
                                  ks, 
                                  intercept, temporal_trend = 1, sig_st = 0.05,
                                  t_axis = NULL, beta1_t = NULL, beta2_t = NULL,
                                  n_sim = 1){
  #function that simulates true risk-surface
  
  # Set seed
  set.seed(seed)

  ## Draw the tensor product smooth parameters: Each row is one sample 
  parameters <- mvrnorm(n = n_sim, rep(0, kt * ks), sig_st * Sigma_st) 
  
  ## Calculate interactions: Each column now represents a different data set
  interactions <- Bst %*% t(parameters)
  
  ## Make the temporal trend
  if(temporal_trend == 1){
    temporal_effect = rep(0, dim(parameters)[2])
  } else if(temporal_trend == 2){
    tmp_ = t_axis * beta1_t; temporal_effect = rep(0, dim(parameters)[2])
    for(t in 1:length(t_axis)){
      temporal_effect[((t - 1) * (dim(Bs)[1]) + 1):(t * dim(Bs)[1])] = tmp_[t]
      }
    } else if(temporal_trend == 3){
      tmp_ = ifelse(t_axis < 7, 
                    beta1_t * t_axis, 
                    beta1_t * 7 - beta2_t * (t_axis - 7))
      temporal_effect = rep(0, dim(parameters)[2])
      for(t in 1:length(t_axis)){
        temporal_effect[((t - 1) * (dim(Bs)[1]) + 1):(t * dim(Bs)[1])] = tmp_[t]
      }
    }
  
  ## Get the sampled risk-fields
  #Lambda_st <- exp(as.vector(intercept + temporal_effect + interactions))
  Lambda_st <- exp(intercept + temporal_effect + interactions)
  
  return(Lambda_st)
}


sample_counts <- function(seed,
                          lambda.df){
  ## Draw from a Poisson distribution the simulated counts
  
  ## Extract matrix w. n_sim columns each corresponding to a data set to be simulated
  Mu = lambda.df$mu
  
  ## Initialize a matrix to save the sampled counts to
  sampled_counts = matrix(nrow = dim(Mu)[1], ncol = dim(Mu)[2])
  
  ## Set seed
  set.seed(seed)
  
  ## For each col (data set) sample from a poisson distribution
  for(col_id in 1:dim(Mu)[2]){
    sampled_counts[, col_id] = sapply(Mu[, col_id], FUN = function(x){return(rpois(1, x))})
  }
  
  return(sampled_counts)
}

################################################################################
#Functions used for model fits

### Function needed to sort the proper-models back again
sort_proper_fitted <- function(proper_fitted, n, tT){
  sorted_proper_fitted <- proper_fitted
  for(t in 1:tT){
    sorted_proper_fitted[((t-1)*n + 1):(t*n)] =  proper_fitted[seq(t, n*tT, by = tT)]
  }
  return(sorted_proper_fitted)
}

## Constraint maker function
constraints_maker <- function(type = NULL, n = NULL, t = NULL,
                              rw = "RW1", prec_matrix = NULL){
  #Type: specifies what interaction type and hence what type of constraint is desired
  #n: specifies number of areas
  #t: specifies number of time points
  #rw: specifies if temporal random effects follows a RW1 or RW2
  #prec_matrix: Precision matrix, used to define constraints using eigenvectors (only for RW2)
  
  if(rw == "RW1"){ #Define constraints for RW(1)
    if(type == "II"){
      #For a type II interaction, there is a RW(1) over each interaction (assuming RW1)
      #Hence for each county, constrain the RW to sum-to-zero
      A <- matrix(0, nrow = n, ncol = n * t)
      for (i in 1:(n - 1)) {
        A[i, which((1:(n * t))%%n == i)] <- 1
      }
      A[n, which((1:(n * t))%%n == 0)] <- 1
      
    } else if(type == "III"){
      #For a type III interaction, there is a indep. ICAR at each time point
      #Need the ICAR at each time point to sum-to-zero
      A <- matrix(0, nrow = t, ncol = n * t)
      for (i in 1:t) {
        # The ICAR at each time point needs to sum to 0
        A[i, ((i - 1) * n + 1):(i * n)] <- 1
      }
      
    } else if(type == "IV"){
      #For a type IV interaction, we have to do both sum-to-zero
      #over each RW on each county, and for each time point sum-to-zero
      #over each ICAR
      time_constr <- matrix(0, nrow = n, ncol = n * t)
      for (i in 1:(n - 1)) {
        time_constr[i, which((1:(n * t))%%n == i)] <- 1
      }
      time_constr[n, which((1:(n * t))%%n == 0)] <- 1
      
      space_constr <- matrix(0, nrow = t-1, ncol = n * t)
      for (i in 1:(t-1)) { 
        space_constr[i, ((i - 1) * n + 1):(i * n)] <- 1
      }
      
      A <- rbind(time_constr, space_constr)
    }
  } else{ #Define constraints with RW2
    if(type == "II"){
      #For a type II interaction, there is a RW(2) over each interaction (assuming RW2)
      eigens <- eigen(prec_matrix)
      
      #Extract 2n last eigenvectors corresponding to eigenvalues=0
      A <- t(eigens$vectors[ ,(nrow(eigens$vectors) - 2 * n + 1):nrow(eigens$vectors)])
      
    } else if(type == "III"){
      #For a type III interaction, there is a indep. ICAR at each time point
      #Need the ICAR at each time point to sum-to-zero
      A <- matrix(0, nrow = t, ncol = n * t)
      for (i in 1:t) {
        # The ICAR at each time point needs to sum to 0
        A[i, ((i - 1) * n + 1):(i * n)] <- 1
      }
      
    } else if(type == "IV"){
      #Calculate eigenvectors
      eigens <- eigen(prec_matrix)
      
      #Extract 2n + T - 2 last eigenvectors corresponding to eigenvalues=0
      A <- t(eigens$vectors[ ,(nrow(eigens$vectors) - 2 * n - t + 3):nrow(eigens$vectors)])
    }
  }
  
  #Get constraints in INLA format
  constr.st <- list(A = A, e = rep(0, dim(A)[1]))
  return(constr.st)
}