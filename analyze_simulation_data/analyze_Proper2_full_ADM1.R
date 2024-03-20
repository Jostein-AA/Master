#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

################################################################################
# Create formulas


#inla.pc.ar.lambda = function(a = 0.5, b = 0.8, p = 4)
#{
#  inla.pc.ar.solve.lambda = function(pred.err.factors, nseq = 1000L) {
## pred.err.factor = E(1-rho^2)
#    pred.err = function(lambda) {
#      return(0.5 * lambda * sqrt(pi) * exp(lambda^2/4 + log_erfc(lambda/2)))
#    }
## find lower and upper limit of lambda for the pred.err.factors given
#    lambda.min = lambda.max = 1
#    val.max = max(pred.err.factors)
#    val.min = min(pred.err.factors)
#    while(pred.err(lambda.min) > val.min) {
#      lambda.min = lambda.min / 2.0
#    }
#   while(pred.err(lambda.max) < val.max) {
#      lambda.max = lambda.max * 2.0
#    }
#    lambda = lambda.min * exp(seq(0, log(lambda.max/lambda.min), length = nseq))
#    fun = splinefun(pred.err(lambda), lambda, method = "monoH.FC")
#    lambdas = fun(pred.err.factors)
#    return (lambdas)
#  }
#  pred.err.factors = 1.0 - (1.0-a)*b^(0:(p-1))
#  lambda = inla.pc.ar.solve.lambda(pred.err.factors)
#  return (lambda)
#}

#test = inla.pc.ar.lambda()

#inla.pc.ar.test1 = function(p=1, p.est = p+1, n = 100)
#{
#  pacf = c(runif(p), rep(0, p.est-p))
#  phi = inla.ar.pacf2phi(pacf)
#  sd.x = sqrt(1/prod(1-pacf^2))
#  x = arima.sim(n = n, model = list(ar = phi))/sd.x
#  lambda = c(inla.pc.ar.lambda(p = p.est, b = 0.5), rep(1, 10))
#  initial = c(inla.models()$latent$ar$hyper$theta2$to.theta(pacf), rep(0, 10))
#  r = (inla(
#    x ~ -1 + f(time, model = "ar", order = p.est,
#               hyper = list(
#                 prec = list(param = c(3, 0.01), initial = 0),
#                 pacf1 = list(param = c(lambda[1], 0), initial = initial[1]),
#                 pacf2 = list(param = c(lambda[2], 0), initial = initial[2]),
#                 pacf3 = list(param = c(lambda[3], 0), initial = initial[3]),
#                 pacf4 = list(param = c(lambda[4], 0), initial = initial[4]),
#                 pacf5 = list(param = c(lambda[5], 0), initial = initial[5]),
#                 pacf6 = list(param = c(lambda[6], 0), initial = initial[6]),
#                 pacf7 = list(param = c(lambda[7], 0), initial = initial[7]),
#                 pacf8 = list(param = c(lambda[8], 0), initial = initial[8]),
#                 pacf9 = list(param = c(lambda[9], 0), initial = initial[9]),
#                 pacf10 = list(param = c(lambda[10], 0), initial = initial[10]))),
#    data = data.frame(x=x, time = 1:length(x)),
#    control.family = list(hyper = list(prec = list(initial = 12,
#                                                   fixed=TRUE)))))
#  result = cbind(est = r$summary.hyperpar$mean[-1], true = pacf)
#  inside = c()
#  for(i in 1:p.est) {
#    int = inla.hpdmarginal(0.95, r$marginals.hyperpar[[i+1]])
#    inside[i] = (int[1] < pacf[i] && pacf[i] < int[2])
#  }
#  result = cbind(result, coverage = inside)
#  print(round(result, digits=3))
#  result = (cbind(acf.est = inla.ar.pacf2acf(r$summary.hyperpar$mean[-1], lag.max = 10),
#                  acf.emp = c(acf(x, lag.max=10, plot=FALSE)$acf),
#                  acf.true = inla.ar.pacf2acf(pacf, lag.max=10)))
#  print(round(result, digits=3))
# return (invisible())
#}










## Specify priors for hyperparameters of proper models
#---
### Temporal hyperparameters (prec. of AR1 and AR1's mixing param) w. corresponding priors: penalized constraint 
ar_hyper = list(prec = list(prior = 'pc.prec', 
                            param = c(1, 0.01))) #, mean = list(prior = 'normal', param = c(0, 1), fixed = TRUE)) 


### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01))) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) 

#---

## Specify precision matrices
#---
### Make precision matrix for Besag on ADM1
matrix4inla <- nb2mat(nb_first_level, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_first_level <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse


# get scaled Besag
#scaled_besag_prec_first_level <- INLA::inla.scale.model(Besag_prec_first_level, 
#                                                        constr = list(A = matrix(1,1,dim(Besag_prec_first_level)[1]), 
#                                                                      e = 0))
#---

#Define a model with intercept, fixed temporal effect and spatiotemporal interaction by ar2 and properbesag
proper_full_formula_first_level <- sampled_counts ~ 1 + time_id +
                                                      f(time_id.copy,
                                                        model = "ar",
                                                        order = 2,
                                                        hyper = ar_hyper) + 
                                                      f(area_id, 
                                                        model = "besagproper2",
                                                        graph = Besag_prec_first_level,
                                                        hyper = spatial_hyper) + 
                                                      f(area_id.copy, 
                                                        model = "besagproper2",
                                                        graph = Besag_prec_first_level,
                                                        hyper = spatial_hyper,
                                                        group = time_id, 
                                                        control.group = list(model = "ar", 
                                                                             order = 2))

################################################################################

## Should probably do a try-catch something for each fit, and also save eventual failures
tryCatch_inla <- function(data,
                          data_set_id,
                          csv_tracker_filename,
                          model_name, scenario_name) {
  tryCatch(
    {
      ## Set an upper-time limit for inla before a timeout
      inla.setOption(inla.timeout = 480) # Set to 480 sec (8 minutes!)
      
      tmp_ = inla(proper_full_formula_first_level, 
                  data = data, 
                  family = "poisson",
                  E = E_it, #E_it
                  control.predictor = list(compute = TRUE,
                                           link = 1),       #For predictions
                  #control.family = list(control.link = list(model = "log")),
                  control.compute = list(config = TRUE, # To see constraints later
                                         cpo = T,       # For model selection
                                         return.marginals.predictor=TRUE)) # Get the lin.pred.marginal
      
      if(tmp_$ok == FALSE){ ## INLA has crashed
        # Update tracker
        tracker.df <- read.csv(csv_tracker_filename)
        tracker.df[data_set_id, ]$error = data_set_id
        write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      }
      
      if(tmp_$mode$mode.status > 0){ ## Potentially something weird with the mode of hyperparameters
        # Update tracker
        tracker.df <- read.csv(csv_tracker_filename)
        tracker.df[data_set_id, ]$warning = data_set_id
        write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      }
      
      
      ### Save the linear predictor-marginal distribution and the CPO-values
      filename_to_save <- paste("./results/", model_name, "/", scenario_name, "/", 
                                model_name, "_", scenario_name, "_", toString(data_set_id), ".RData", 
                                sep = "")
      
      marginals = tmp_$marginals.fitted.values 
      cpo = tmp_$cpo$cpo
      
      save(marginals, 
           cpo,
           file = filename_to_save)
    },
    error = function(cond) {
      print("!Error!")
      
      # Update tracker
      tracker.df <- read.csv(csv_tracker_filename)
      tracker.df[data_set_id, ]$error = data_set_id
      write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      
      # Choose a return value in case of error
      -1
    },
    warning = function(cond) {
      print("!warning!")
      print(cond)
      
      # Update tracker
      tracker.df <- read.csv(csv_tracker_filename)
      tracker.df[data_set_id, ]$warning = data_set_id
      write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      
      print("---")
    },
    finally = {
      print(paste("[_/(^::^)--|", model_name, scenario_name, 
                  toString(data_set_id), sep = " "))
      
      # Update tracker
      tracker.df <- read.csv(csv_tracker_filename)
      tracker.df[data_set_id, ]$analyzed = data_set_id
      write.csv(tracker.df, file = csv_tracker_filename, row.names = F)
      
    }
  )
}



################################################################################
# SC1
model_name = "proper2_full"
scenario_name = "sc1"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
  
}


tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC3
model_name = "proper2_full"
scenario_name = "sc3"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC5
model_name = "proper2_full"
scenario_name = "sc5"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC7
model_name = "proper2_full"
scenario_name = "sc7"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC9
model_name = "proper2_full"
scenario_name = "sc9"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))

################################################################################
# SC11
model_name = "proper2_full"
scenario_name = "sc11"

## Get the tracker-filename
csv_tracker_filename = get_csv_tracker_filename(model_name, scenario_name)


not_finished = T
while(not_finished){
  ## Load in newest data set
  ### Start by evaluating the tracker
  data_set_id = get_first_not_yet_analyzed(model_name, scenario_name)
  tracker.df = read.csv(csv_tracker_filename)
  
  ## If all the data sets have been analyzed for this scenario, move on!
  if(data_set_id == nrow(tracker.df)){
    if(is.na(tracker.df[nrow(tracker.df), ]$analyzed)){
      print("last")
    } else{
      not_finished = F
      print("Finished")
      break
    }
  }
  
  ### Load in sc1 simulated data
  load(paste("./Simulated_data/", scenario_name, "/", scenario_name, "_data.RData", sep = ""))
  lambda_df <- lambda.df[, c("area_id", "time_id", "E_it", 
                             "space.time")]
  
  lambda_df$sampled_counts = lambda.df$sampled_counts[, data_set_id]
  
  ## Set the last three years counts to NA for the fit
  lambda_df[lambda_df$time_id %in% 11:13, ]$sampled_counts = NA
  
  ## Reorder due to change in space.time interaction
  lambda_df <- lambda_df[order(lambda_df$area_id, decreasing = F), ]
  rownames(lambda_df) <- 1:nrow(lambda_df)
  
  ## Add copies of area and time ids, INLA requires unique random effects
  lambda_df$area_id.copy <- lambda_df$area_id
  lambda_df$time_id.copy <- lambda_df$time_id
  
  
  ## Do tryCatch
  fitted_inla_sc1 <- tryCatch_inla(lambda_df,
                                   data_set_id,
                                   csv_tracker_filename,
                                   model_name, scenario_name)
  
}

tracker.df = read.csv(csv_tracker_filename)
print(paste("Number of errors: ", sum(!is.na(tracker.df$error))))