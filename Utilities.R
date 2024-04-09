################################################################################
# Functions for calculating MSE, IS,...

count_mse_one_year_one_dataset <- function(sampled_counts_one_year, 
                                           lambda_marginals_one_year,
                                           E_it){
  
  ## For each area find the expected predicted count (i.e. point prediction)
  pred_count <- 100 * as.numeric(sapply(lambda_marginals_one_year, 
                                        FUN = function(x){return(mean(x[, 1]))}))
  
  ## Calculate the MSE
  mse_one_year <- mean((sampled_counts_one_year - pred_count)**2)
  
  # Return the MSE of that year
  return(mse_one_year)
}

## Testing
#one_year_margs <- marginals[(11 * n_ADM1 + 1):(12 * n_ADM1)]
#count_mse_one_year_one_dataset(lambda.df$sampled_counts[1],one_year_margs,100)




rate_mse_one_year_one_dataset <- function(sampled_rates_one_year, 
                                          lambda_marginals_one_year){
  
  ## For each area find the expected predicted count (i.e. point prediction)
  pred_count <- as.numeric(sapply(lambda_marginals_one_year, 
                                  FUN = function(x){return(mean(x[, 1]))}))
  
  ## Calculate the MSE
  mse_one_year <- mean((sampled_rates_one_year - pred_count)**2)
  
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
  
  return(list(l = l, u = u, median = median))
  
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


count_IS_one_year_one_dataset <- function(sampled_counts_one_year,
                                          lambda_marginals_one_year,
                                          E_it){
  
  ## Find the l and u for each singular instance in a year
  ul_each_one_year <- lapply(lambda_marginals_one_year, 
                             FUN = function(x){
                               return(find_ul_quants_counts_single_pred(x, 100))
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





rate_IS_one_year_one_dataset <- function(sampled_rates_one_year,
                                         lambda_marginals_one_year){
  
  ## Find the l and u for each singular instance in a year
  ul_each_one_year <- lapply(lambda_marginals_one_year, 
                             FUN = function(x){
                               return(find_ul_quants_rate_single_pred(x))
                             })
  
  ## Find the IS for each singular instance in a year
  ### Initialize space for IS 
  IS_each_instance = rep(0, length(lambda_marginals_one_year))
  
  ### Calculate IS each area
  for(i in 1:length(lambda_marginals_one_year)){
    IS_each_instance[i] = find_IS_one_obs(ul_each_one_year[[i]]$l, ul_each_one_year[[i]]$u, 
                                          sampled_rates_one_year[i])
  }
  
  ## Find the average IS this year
  IS_this_year <- mean(IS_each_instance)
  
  return(IS_this_year)
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
      theme(plot.title = element_text(size = 15, hjust = 0.5),
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
        guide = "colorscale")
  } else {
    ggplot(data = map_w_values) +  
      geom_sf(aes(fill = to_plot), 
              alpha = 1,
              color="black") + ggtitle(title) + 
      theme(plot.title = element_text(size = 15, hjust = 0.5),
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
    }
}

heatmap_points <- function(risk_surface.list, 
                           polygons,
                           admin_map,
                           t,
                           title){
  
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
    guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
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