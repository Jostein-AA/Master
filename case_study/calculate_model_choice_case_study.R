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

load("case_study/imp_typeIV_test.RData")
load("case_study/proper2_onlyInt_Spain.RData")

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

################################################################################



