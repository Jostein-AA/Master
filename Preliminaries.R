#Clear environment
rm(list = ls())

#Load libraries
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)
library(ggpubr)
library(pspline)
library(readxl)
library(OneR)
library(mgcv)
library(splines)
library(ggfortify)
library(crs)
library(INLA)
library(MASS)

source("Utilities.R")

################################################################################
# Preliminaries w. respect to Maps

## Load the raw map-data
second_level_admin_map  <- read_sf('./Shapefiles/Germany')[,c('NAME_2', 'TYPE_2', 'ID_2', 'geometry')]
first_level_admin_map <- read_sf('./Shapefiles/Germany_first_level_unpacked/')[,c('NAME_1', 'TYPE_1', 'ID_1', 'geometry')]
germany_border <- read_sf('./Shapefiles/Germany_border/')[, c('geometry')]
Population_data <- read_excel('./Data/germany_population.xlsx')

new_map <- read_sf('./Shapefiles/DEU_adm/DEU_2/')[, c("ID_2", "NAME_2", "geometry")]

## Make sure no NA's in the map!
second_level_admin_map <- na.omit(second_level_admin_map)

## Remove water body from second_level_admin_map
second_level_admin_map <- second_level_admin_map[second_level_admin_map$TYPE_2 != 'Water body', ]

## Omit population data that is not from Germany
Population_data <- Population_data[grepl('DE', Population_data$Code), ]

## Create TYPE_2 column for population data
NAME_2 = rep('', nrow(Population_data)); TYPE_2 = rep('', nrow(Population_data))
for(i in 1:nrow(Population_data)){
  if(grepl(',', Population_data$Name[i])){
    temp = strsplit(Population_data$Name[i], split = ', ')
    NAME_2[i] = temp[[1]][1]; TYPE_2[i] = temp[[1]][2]
  } else{NAME_2[i] =Population_data$Name[i]; TYPE_2[i] = ''}
}
Population_data$NAME_2 = NAME_2; Population_data$TYPE_2 = TYPE_2

## Merge population data to map of Germany
second_level_admin_map$Population <- rep(0, nrow(second_level_admin_map))
for(i in 1:nrow(Population_data)){
  #If name is in Germany map names
  if(Population_data$NAME_2[i] %in% second_level_admin_map$NAME_2){
    if(Population_data$TYPE_2[i] == ''){ #If landskreis, statdkreis, etc not specified merge directly
      second_level_admin_map$Population[
        second_level_admin_map$NAME_2 == Population_data$NAME_2[i]] <- Population_data$Pop[i]
    } else{ #If landskreis, statdkreis, etc specified, merge accordingly
      second_level_admin_map$Population[
        second_level_admin_map$NAME_2 == Population_data$NAME_2[i] &
          second_level_admin_map$TYPE_2 == Population_data$TYPE_2[i]] <- Population_data$Pop[i]
    }
  }
}
## Special cases
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Dillingen an der Donau'] = Population_data$Pop[Population_data$NAME_2 == 'Dillingen a.d. Donau']
second_level_admin_map$Population[98] = 115436
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Neumarkt in der Oberpfalz'] = Population_data$Pop[Population_data$NAME_2 == 'Neumarkt i. d. OPf.']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Neustadt an der Aisch-Bad Windsheim'] = Population_data$Pop[Population_data$NAME_2 == 'Neustadt a. d. Aisch-Bad Windsheim']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Neustadt an der Waldnaab'] = Population_data$Pop[Population_data$NAME_2 == 'Neustadt a. d. Waldnaab']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Pfaffenhofen an der Ilm'] = Population_data$Pop[Population_data$NAME_2 == 'Pfaffenhofen a. d. Ilm']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Weiden in der Oberpfalz'] = Population_data$Pop[Population_data$NAME_2 == 'Weiden i. d. Opf']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Wunsiedel im Fichtelgebirge'] = Population_data$Pop[Population_data$NAME_2 == 'Wunsiedel i. Fichtelgebirge']
second_level_admin_map$Population[193] = Population_data$Pop[Population_data$NAME_2 == 'Landkreis Rostock']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Friesland'] = Population_data$Pop[Population_data$NAME_2 == 'Friesland (DE)']
second_level_admin_map$Population[223] = Population_data$Pop[Population_data$NAME_2 == 'Oldenburg (Oldenburg)']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Osterode am Harz'] = 23993
second_level_admin_map$Population[382] = 42000
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Ilm-Kreis'] = Population_data$Pop[Population_data$NAME_2 == 'Ilm-Kreis (NUTS 2021)']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Saalfeld-Rudolstadt'] = Population_data$Pop[Population_data$NAME_2 == 'Saalfeld-Rudolstadt (NUTS 2021)']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Schmalkalden-Meiningen'] = Population_data$Pop[Population_data$NAME_2 == 'Schmalkalden-Meiningen (NUTS 2021)']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Sonneberg'] = Population_data$Pop[Population_data$NAME_2 == 'Sonneberg (NUTS 2021)']
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Suhl'] = 37000
second_level_admin_map$Population[second_level_admin_map$NAME_2 == 'Wartburgkreis'] = Population_data$Pop[Population_data$NAME_2 == 'Wartburgkreis (NUTS 2021)']

## Reset the ID as some rows have been removed
second_level_admin_map$ID_2 = 1:nrow(second_level_admin_map)

## Load in structures
nb_second_level <- spdep::poly2nb(second_level_admin_map, queen = FALSE)
nb_first_level <- spdep::poly2nb(st_make_valid(first_level_admin_map), queen = FALSE)
nb_new_level <- spdep::poly2nb(st_make_valid(new_map), queen = FALSE)

## Get the coordinate system used
crs_ = st_crs(germany_border, parameters = TRUE)

#Remove things we do not wish to keep going forward
rm(list = c('temp', 'i', 'NAME_2', 'TYPE_2', 'Population_data'))

save(crs_, nb_first_level, nb_second_level, nb_new_level,
     germany_border,
     first_level_admin_map,
     second_level_admin_map,
     new_map,
     file = "maps_and_nb.RData")

print("Maps done")
################################################################################
# Preliminaries w. respect to creating points where the risk is evaluated at

## Get a polygon_grid over Germany
polygon_grid = st_as_sf(st_make_grid(germany_border, 
                                     n = c(110, 110), 
                                     square = FALSE),
                        crs = crs_$srid)

## Plot the entire polygon_grid
plot(st_geometry(germany_border), border = "black")
plot(st_geometry(polygon_grid), add = T)

## Add unique polygon_id
polygon_grid$polygon_id = 1:nrow(polygon_grid)

## Remove polygons outside of Germany
polygon_grid2 <- sf::st_intersection(polygon_grid, 
                                     st_make_valid(germany_border))

## Plot the polygon_grid2
plot(st_geometry(polygon_grid2), border = "black")


## Get the centroids of each polygon in the grid
centroids_grid <- st_centroid(polygon_grid)

## Get the coordinates of the centroids
gpts = sf::st_coordinates(centroids_grid)

## standardize x and y
locs <- apply(gpts, 2, function(x) (x - min(x)) / (max(x) - min(x)))

## Extract coordinates lat/lon
x <- locs[, 1]
y <- locs[, 2]

## Make t_axis: spatial domain (which discretely is: 1,...,13, and cont [0, 13])
t_axis <- seq(0.0001, 13, length.out = 40); xt_ = seq(0, 13, length.out = 40)

## Scale the temporal 'coordinates'
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Transform yxt_grid to have geometry points
yxt_geom = st_as_sf(yxt_grid, 
                    coords = c("x", "y"),
                    crs = crs_$srid)

yxt_geom$x = yxt_grid$x; yxt_geom$y = yxt_grid$y
yxt_geom$unique_id = 1:nrow(yxt_geom)


## Add time-id to yxt_geom
yxt_geom$time_id = ceiling(yxt_geom$t)

## State first and last time point
t0 = round(t_axis[1], digits = 0); tT = round(t_axis[length(t_axis)], digits = 0)


## Find indices within Germany
within_germany_indices = as.numeric(rownames(sf::st_intersection(yxt_geom, 
                                                                 st_make_valid(germany_border))))


## Find mapping to each region in second_level_admin_map2
first_level_area_id_mapping <- st_join(yxt_geom, st_make_valid(first_level_admin_map),
                                        join = st_within)$ID_1

yxt_geom$first_level_area_id_mapping = first_level_area_id_mapping

## Find mapping to each region in new map
new_level_area_id_mapping <- st_join(yxt_geom, st_make_valid(new_map),
                                       join = st_within)$ID_2

yxt_geom$new_level_area_id_mapping = new_level_area_id_mapping

## Find mapping to each region in second_level_admin_map
second_level_area_id_mapping = rep(NA, nrow(yxt_geom))
second_level_area_id_mapping_tmp <- st_join(yxt_geom, st_make_valid(second_level_admin_map),
                                           join = st_within)

## Drop NA-values
second_level_area_id_mapping_tmp = drop_na(second_level_area_id_mapping_tmp)

## Add area IDs to yxt_geom
second_level_area_id_mapping[second_level_area_id_mapping_tmp$unique_id] = second_level_area_id_mapping_tmp$ID_2
yxt_geom$second_level_area_id_mapping = second_level_area_id_mapping


## Find mapping for the areas w.o. associatied points for second_level_admin_map 
missing_mapping_second_level_admin_map.df <- data.frame(missing_id = rep(NA, nrow(yxt_geom)),
                                             time_id = rep(NA, nrow(yxt_geom)),
                                             yxt_geom_unique_id = yxt_geom$unique_id)


### Find the regions within Germany where there is no corresponding grid point being evaluated
missing_areas_IDs = c()
for(i in unique(second_level_admin_map$ID_2)){
  tmp_ = yxt_geom[within_germany_indices &
                    yxt_geom$second_level_area_id_mapping == i &
                    yxt_geom$time_id == 1, ]
  tmp_ = drop_na(tmp_)
  if(is.na(mean(tmp_$polygon_id))){
    missing_areas_IDs = c(missing_areas_IDs, i)
  }
}

### Find the grid point closest to the area, and use polygon_id to identify area, and use
### unique_id to map for each specific time
tmp_ = yxt_geom[within_germany_indices, ]
for(missing_area in missing_areas_IDs){
  for(j in 1:length(t_axis)){
    #print(t)
    tmp2_ = tmp_[tmp_$t == t_axis[j], ]
    index_tmp2_ = st_nearest_feature(st_make_valid(second_level_admin_map[second_level_admin_map$ID_2 == missing_area, ]),
                                     tmp2_)
    
    
    ## Add values
    missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$yxt_geom_unique_id == tmp2_[index_tmp2_, ]$unique_id, ]$missing_id = missing_area
    missing_mapping_second_level_admin_map.df[missing_mapping_second_level_admin_map.df$yxt_geom_unique_id == tmp2_[index_tmp2_, ]$unique_id,  ]$time_id = tmp2_[index_tmp2_, ]$time_id 
  }
}

## Drop NAs
missing_mapping_second_level_admin_map.df = drop_na(missing_mapping_second_level_admin_map.df)

## Store which areas are unproblematic
no_problem_areas = unique(second_level_admin_map$ID_2)
no_problem_areas = no_problem_areas[!(no_problem_areas %in% missing_mapping_second_level_admin_map.df$missing_id)]


save(no_problem_areas,
     yxt_geom,
     within_germany_indices,
     missing_mapping_second_level_admin_map.df,
     t0, tT, t_axis, xt, xt_,
     polygon_grid, polygon_grid2,
     file = "grids_and_mappings.RData")

print("Grid and mapping done")

################################################################################
# Preliminaries related to the B-spline basis functions, P-splines and tensor prod. smooths

## define B-spline basis degree
bdeg = 3

#########################################
## Make the temporal B-spline basis functions

## Specify number of intervals in time (knots - 1)
n_intervals_t = 5

## Find distance between knots
dist <- (max(xt) - min(xt)) / n_intervals_t
xtl <- min(xt) - dist * 0.05; xtr <- max(xt) + dist * 0.05
dxt <- (xtr - xtl) / n_intervals_t

## Define knot locations in time
knots <- seq(xtl - bdeg * dxt, xtr + bdeg * dxt, by = dxt)

## Define temporal marginal B-spline basis design matrix
Bt <- spline.des(knots, xt, bdeg + 1)$design
kt <- ncol(Bt)



save(Bt, kt,
     file = "temporal_B_spline.RData")

#########################################
## 20 knots in latitude and longitude directions

## Define number of intervals in each of the spatial directions (knots - 1) 
n_intervals_s = 19 

# Specify the knot locations and create the B-spline basis design matrices
## Find the distance between each knot in x- and y-directions
disx <- (max(x) - min(x)) / n_intervals_s; disy <- (max(y) - min(y)) / n_intervals_s

## Find the left (l) and right (r) end and a bitmore spline not 0 at ends!
xl <- min(x) - disx * 0.05; yl <- min(y) - disy * 0.05
xr <- max(x) + disx * 0.05; yr <- max(y) + disy * 0.05

## This is distance between each knot?
dx <- (xr - xl) / n_intervals_s; dy <- (yr - yl) / n_intervals_s

## Specify the knot locations in x1 and x2 direction
knotsx <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
knotsy <- seq(yl - bdeg * dy, yr + bdeg * dy, by = dy)

## Make marginal design matrices for B-splines in x and y direction 
Bx_20 <- spline.des(knotsx, x, bdeg + 1, 0 * x)$design  ##Base marginal x
By_20 <- spline.des(knotsy, y, bdeg + 1, 0 * y)$design  ##Base marginal y

## Perform row-wise Kronecker product to get the spatial B-spline basis
Bs_20 = row_wise_Kronecker(Bx_20, By_20)

## Find dimensions (spatial)
kx_20 <- dim(Bx_20)[2]; ky_20 <- dim(By_20)[2]; ks_20 <- dim(Bs_20)[2]

#########################################
## 10 knots in latitude and longitude directions

## Define number of intervals in each of the spatial directions (knots - 1)  
n_intervals_s = 9 

# Specify the knot locations and create the B-spline basis design matrices
## Find the distance between each knot in x- and y-directions
disx <- (max(x) - min(x)) / n_intervals_s; disy <- (max(y) - min(y)) / n_intervals_s

## Find the left (l) and right (r) end and a bitmore spline not 0 at ends!
xl <- min(x) - disx * 0.05; yl <- min(y) - disy * 0.05
xr <- max(x) + disx * 0.05; yr <- max(y) + disy * 0.05

## This is distance between each knot?
dx <- (xr - xl) / n_intervals_s; dy <- (yr - yl) / n_intervals_s

## Specify the knot locations in x1 and x2 direction
knotsx <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
knotsy <- seq(yl - bdeg * dy, yr + bdeg * dy, by = dy)

## Make marginal design matrices for B-splines in x and y direction 
Bx_10 <- spline.des(knotsx, x, bdeg + 1, 0 * x)$design  ##Base marginal x
By_10 <- spline.des(knotsy, y, bdeg + 1, 0 * y)$design  ##Base marginal y

## Perform row-wise Kronecker product to get the spatial B-spline basis
Bs_10 = row_wise_Kronecker(Bx_10, By_10)

## Find dimensions (spatial)
kx_10 <- dim(Bx_10)[2]; ky_10 <- dim(By_10)[2]; ks_10 <- dim(Bs_10)[2]

save(Bs_20, kx_20, ky_20, ks_20, Bx_20, By_20,
     Bs_10, kx_10, ky_10, ks_10, Bx_10, By_10,
     file = "spatial_B_splines.RData")

print("Marginal B-splines done done")
#########################################

## Define the 3-dimensional B-spline basis functions
Bst_20 <- kronecker(Bt, Bs_20)
Bst_10 <- kronecker(Bt, Bs_10)

save(Bst_20, Bst_10,
     file = "B_spline_basis.RData")

print("3-dim. B-spline matrix done")

# Define the penalization matrices and get the risk field
## Define the marginal penalization matrices (RW1 precision matrices)
Pt <- INLA:::inla.rw(n = kt, order = 1, 
                        scale.model = FALSE, 
                        sparse = TRUE)

Px_20 <- INLA:::inla.rw(n = kx_20, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Py_20 <- INLA:::inla.rw(n = ky_20, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Px_10 <- INLA:::inla.rw(n = kx_10, order = 1, 
                        scale.model = FALSE, 
                        sparse = TRUE)

Py_10 <- INLA:::inla.rw(n = ky_10, order = 1, 
                        scale.model = FALSE, 
                        sparse = TRUE)


## Identity matrices k1 x k1, k2 x k2, and kt x kt
It <- diag(kt)
Ix_20 <- diag(kx_20); Iy_20 <- diag(ky_20)
Ix_10 <- diag(kx_10); Iy_10 <- diag(ky_10)

## Make the covariance matrix using generalized inverse

### Scenario 1 and 2
tau_s = 20; tau_t = 40 #1) spatial and 2) temporal smoothness parameters
Pst_sc_12 <- tau_s *(It %x% Py_20 %x% Ix_20 + It %x% Iy_20 %x% Px_20) + 
  tau_t * (Pt %x% Iy_20 %x% Ix_20)

### Scenario 3, 4, 5, and 6
tau_s = 20; tau_t = 80 #1) spatial and 2) temporal smoothness parameters
Pst_sc_3456 <- tau_s *(It %x% Py_20 %x% Ix_20 + It %x% Iy_20 %x% Px_20) + 
  tau_t * (Pt %x% Iy_20 %x% Ix_20)

### Scenario 7 and 8
tau_s = 20; tau_t = 40 #1) spatial and 2) temporal smoothness parameters
Pst_sc_78 <- tau_s *(It %x% Py_10 %x% Ix_10 + It %x% Iy_10 %x% Px_10) + 
  tau_t * (Pt %x% Iy_10 %x% Ix_10)

### Scenario 9, 10, 11, and 12
tau_s = 20; tau_t = 80 #1) spatial and 2) temporal smoothness parameters
Pst_sc_9101112 <- tau_s *(It %x% Py_10 %x% Ix_10 + It %x% Iy_10 %x% Px_10) + 
                  tau_t * (Pt %x% Iy_10 %x% Ix_10)

### Scenario 13
tau_s = 10; tau_t = 10 #1) spatial and 2) temporal smoothness parameters
Pst_sc_13 <- tau_s *(It %x% Py_20 %x% Ix_20 + It %x% Iy_20 %x% Px_20) + 
  tau_t * (Pt %x% Iy_20 %x% Ix_20)

### Scenario 14, 15
tau_s = 10; tau_t = 20 #1) spatial and 2) temporal smoothness parameters
Pst_sc_1415 <- tau_s *(It %x% Py_20 %x% Ix_20 + It %x% Iy_20 %x% Px_20) + 
  tau_t * (Pt %x% Iy_20 %x% Ix_20)

### Scenario 16
tau_s = 10; tau_t = 10 #1) spatial and 2) temporal smoothness parameters
Pst_sc_16 <- tau_s *(It %x% Py_10 %x% Ix_10 + It %x% Iy_10 %x% Px_10) + 
  tau_t * (Pt %x% Iy_10 %x% Ix_10)

### Scenario 17, 18
tau_s = 10; tau_t = 20 #1) spatial and 2) temporal smoothness parameters
Pst_sc_1718 <- tau_s *(It %x% Py_10 %x% Ix_10 + It %x% Iy_10 %x% Px_10) + 
  tau_t * (Pt %x% Iy_10 %x% Ix_10)

save(Pst_sc_12,
     Pst_sc_3456,
     Pst_sc_78,
     Pst_sc_9101112,
     Pst_sc_13,
     Pst_sc_1415,
     Pst_sc_16,
     Pst_sc_1718,
     file = "penalization_matrices.RData")

print("Penalization matrices made done")

## Scale the penalization matrices
scaled_Pst_sc_12 <- INLA::inla.scale.model(Pst_sc_12, 
                                           constr = list(A = matrix(1,1,dim(Pst_sc_12)[1]), 
                                                         e = 0))

scaled_Pst_sc_3456 <- INLA::inla.scale.model(Pst_sc_3456, 
                                           constr = list(A = matrix(1,1,dim(Pst_sc_3456)[1]), 
                                                         e = 0))

scaled_Pst_sc_78 <- INLA::inla.scale.model(Pst_sc_78, 
                                           constr = list(A = matrix(1,1,dim(Pst_sc_78)[1]), 
                                                         e = 0))

scaled_Pst_sc_9101112 <- INLA::inla.scale.model(Pst_sc_9101112, 
                                             constr = list(A = matrix(1,1,dim(Pst_sc_9101112)[1]), 
                                                           e = 0))

scaled_Pst_sc_13 <- INLA::inla.scale.model(Pst_sc_13, 
                                           constr = list(A = matrix(1,1,dim(Pst_sc_13)[1]), 
                                                         e = 0))

scaled_Pst_sc_1415 <- INLA::inla.scale.model(Pst_sc_1415, 
                                             constr = list(A = matrix(1,1,dim(Pst_sc_1415)[1]), 
                                                           e = 0))

scaled_Pst_sc_16 <- INLA::inla.scale.model(Pst_sc_16, 
                                           constr = list(A = matrix(1,1,dim(Pst_sc_16)[1]), 
                                                         e = 0))

scaled_Pst_sc_1718 <- INLA::inla.scale.model(Pst_sc_1718, 
                                                constr = list(A = matrix(1,1,dim(Pst_sc_1718)[1]), 
                                                              e = 0))


save(scaled_Pst_sc_12,
     scaled_Pst_sc_3456,
     scaled_Pst_sc_78,
     scaled_Pst_sc_9101112,
     scaled_Pst_sc_13,
     scaled_Pst_sc_1415,
     scaled_Pst_sc_16,
     scaled_Pst_sc_1718,
     file = "scaled_penalization_matrices.RData")


### Take the generalized inverse and save
Sigma_st_12      <- ginv(as.matrix(scaled_Pst_sc_12))
print("Sigma_st_12")
Sigma_st_3456    <- ginv(as.matrix(scaled_Pst_sc_3456))
print("Sigma_st_3456")
Sigma_st_78      <- ginv(as.matrix(scaled_Pst_sc_78))
print("Sigma_st_78")
Sigma_st_9101112 <- ginv(as.matrix(scaled_Pst_sc_9101112))
print("Sigma_st_9101112")

Sigma_st_13      <- ginv(as.matrix(scaled_Pst_sc_13))
print("Sigma_st_13")
Sigma_st_1415    <- ginv(as.matrix(scaled_Pst_sc_1415))
print("Sigma_st_1415")
Sigma_st_16      <- ginv(as.matrix(scaled_Pst_sc_16))
print("Sigma_st_16")
Sigma_st_1718 <- ginv(as.matrix(scaled_Pst_sc_1718))
print("Sigma_st_1718")

save(Sigma_st_12,
     Sigma_st_3456,
     Sigma_st_78,
     Sigma_st_9101112,
     Sigma_st_13,
     Sigma_st_1415,
     Sigma_st_16,
     Sigma_st_1718,
     file = "scaled_tensor_prod_smooths_cov_matrices.RData")


print("Done")






