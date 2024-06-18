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


#potential errors seem to occur here (below) for some reason unbeknownst to me:::


### Find coverage of model



### Find number of upper quantiles equal to 2,3,... for each year, likewise for l




### Find the distances of the observations outside of the 95% pred. CIs for each year






















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













