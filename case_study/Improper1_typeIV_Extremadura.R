#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library(bigDM)
################################################################################

# Load in the considered lung-cancer data
data(Data_LungCancer)

# Load in a map of Spain containing 7907 areas
data("Carto_SpainMUN")

# Remove disjointed area
map_Spain <- Carto_SpainMUN[-2454, ]
row.names(map_Spain) = 1:nrow(map_Spain)
problem_area = Carto_SpainMUN[2454, ]

Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID != problem_area$ID, ]

# Extract the areas within the principality of Extremadura
map_Spain <- map_Spain[map_Spain$region == "Extremadura", ]
IDs_Extremadura <- unique(map_Spain$ID)

# Extract the data within Extremadura
Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID %in% IDs_Extremadura, ]

# Get the years and areas
years = unique(Data_LungCancer$year)
tT = length(years)
t.from <- min(Data_LungCancer$year)
t.to <- max(Data_LungCancer$year)

# Create time-ids 1,...,25 instead of 1991,...,2015 for the sake of INLA
Data_LungCancer$year_id <- Data_LungCancer$year - min(Data_LungCancer$year) + 1

# For sake of INLA, transform ID from chr to num
Data_LungCancer$ID <- as.numeric(Data_LungCancer$ID) 

areas = unique(Data_LungCancer$ID)
n = length(areas)

# Ad a area-id starting at 1 and ending at 380 to both map_Spain and Data_LungCancer
map_Spain$area_id = 1:nrow(map_Spain)
Data_LungCancer$area_id = rep(NA, nrow(Data_LungCancer))
for(year_id in 1:tT){
  Data_LungCancer[Data_LungCancer$year_id == year_id, ]$area_id = map_Spain$area_id
}

Data_LungCancer$space.time <- 1:nrow(Data_LungCancer)

################################################################################
## Specify priors for hyperparameters of improper models
#---
### Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec',  param = c(1, 0.01)), 
                      phi = list(prior = 'pc',  param = c(0.5, 0.5))) 

### Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', param = c(1, 0.01)), 
                     phi = list(prior = 'pc', param = c(0.5, 0.5)))

### Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec", param=c(1,0.01)))
#---

## Specify precision matrix
#---
### Specify the RW1 precision matrix
RW1_prec <- INLA:::inla.rw(n = tT, order = 1, 
                           scale.model = FALSE, sparse = TRUE)


# Find the adjacencies
nb_spain <- spdep::poly2nb(map_Spain, queen = T)


### Make precision matrix for Besag on Extremadura map
matrix4inla <- nb2mat(nb_spain, style = "B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_spain<- Matrix(matrix4inla, sparse = TRUE) #Make it sparse
rm(matrix4inla)
#---

## Specify base-formula
base_formula <- obs ~ 1 + 
  f(year_id, 
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
    hyper = spatial_hyper) 


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

#Get formula for type IV
typeIV_formula <- update(base_formula,
                      ~. + f(space.time, 
                             model = "generic0",
                             Cmatrix = typeIV_prec,
                             extraconstr = typeIV_constraints,
                             rankdef = (n + tT - 1), 
                             hyper = interaction_hyper))

Data_LungCancer[Data_LungCancer$year_id %in% 21:25, ]$obs = NA

Improper1_typeIV_Extremadura = inla(typeIV_formula, 
                               data = Data_LungCancer, 
                               family = "poisson",
                               E = pop, 
                               verbose = T,
                               control.predictor = list(compute = TRUE, link = 1),       #For predictions
                               control.compute = list(return.marginals.predictor=TRUE)) # Get the lin.pred.marginal



save(Improper1_typeIV_Extremadura,
     file = "./case_study/Improper1_typeIV_Extremadura.RData")














