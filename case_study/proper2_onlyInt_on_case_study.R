#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

# library("geofacet")
# library(ggh4x)
# library(ggridges)
# library(latex2exp)
# library(geoR)
# library(paletteer)q()
# library('ggsci')
# library(ggstats)
library(bigDM)


#library(tmap)
#library(RColorBrewer)

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

str(Data_LungCancer)
str(map_Spain)

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

# Ad a area-id starting at 1 and ending at 7906 to both map_Spain and Data_LungCancer
map_Spain$area_id = 1:nrow(map_Spain)
Data_LungCancer$area_id = rep(NA, nrow(Data_LungCancer))
for(year_id in 1:tT){
  Data_LungCancer[Data_LungCancer$year_id == year_id, ]$area_id = map_Spain$area_id
}

# Calculate expected number of counts per 100,000
Data_LungCancer$exp_per_1E5 <- (Data_LungCancer$exp * 1E5)/Data_LungCancer$pop
################################################################################
# Create precision matrix required, and specify the proper2_onlyInt formula

## Specify precision matrix
#---
nb_spain <- spdep::poly2nb(map_Spain, queen = T)


### Make precision matrix for Besag on Carto_SpainMUN
matrix4inla <- nb2mat(nb_spain, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec_spain<- Matrix(matrix4inla, sparse = TRUE) #Make it sparse
rm(matrix4inla)
#---

## Specify priors for hyperparameters of proper models
#---
### Temporal hyperparameters (prec. of AR2 and AR2's autocorrelation param) w. corresponding priors: penalized constraint 
ar_hyper = list(prec = list(prior = 'pc.prec', 
                            param = c(1, 0.01)),
                pacf1 = list(prior = 'pc.cor1', 
                             param = c(0.5, 0.5 + 1E-2)),
                pacf2 = list(prior = 'pc.cor0',
                             param = c(0.5, 0.5)))


### Spatial hyperparameters (Leroux prec. and Leroux mixing param) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01))) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) #, lambda = list(prior = 'gaussian', param = c(0, 0.45)) 


### Group hyper
group_hyper = list(pacf1 = list(prior = 'pc.cor1', 
                                param = c(0.5, 0.5 + 1E-2)),
                   pacf2 = list(prior = 'pc.cor0',
                                param = c(0.5, 0.5)))
#---





#Define a model with intercept, fixed temporal effect and spatiotemporal interaction by ar1 and properbesag
proper2_onlyInt_formula <- obs ~ 1 + year_id + 
                            f(area_id, 
                              model = "besagproper2",
                              graph = Besag_prec_spain,
                              hyper = spatial_hyper,
                              group = year_id, 
                              control.group = list(model = "ar", 
                                                   order = 2,
                                                   hyper = group_hyper))



Data_LungCancer[Data_LungCancer$year_id %in% 21:25, ]$obs = NA

## Reorder due to change in space.time interaction
Data_LungCancer <- Data_LungCancer[order(Data_LungCancer$area_id, decreasing = F), ]
rownames(Data_LungCancer) <- 1:nrow(Data_LungCancer)


proper2_onlyInt_Spain = inla(proper2_onlyInt_formula, 
                              data = Data_LungCancer, 
                              family = "poisson",
                              E = pop, 
                              verbose = T,
                              control.predictor = list(compute = TRUE, link = 1),       #For predictions
                              control.compute = list(return.marginals.predictor=TRUE)) # Get the lin.pred.marginal



save(proper2_onlyInt_Spain,
     file = "./case_study/proper2_onlyInt_Spain.RData")












