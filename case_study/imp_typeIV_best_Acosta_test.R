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
row.names(Data_LungCancer) = 1:nrow(Data_LungCancer)

str(Data_LungCancer)
str(map_Spain)

# Get the years and areas
years = unique(Data_LungCancer$year)
tT = length(years)
t.from <- min(Data_LungCancer$year)
t.to <- max(Data_LungCancer$year)

# Create time-ids 1,...,25 instead of 1991,...,2015 for the sake of INLA
Data_LungCancer$year_id <- Data_LungCancer$year - min(Data_LungCancer$year) + 1

# Remove Loads of data
Data_LungCancer <- Data_LungCancer[Data_LungCancer$year_id %in% 1:7, ]
tT = 7

# For sake of INLA, transform ID from chr to num
#Data_LungCancer$ID <- as.numeric(Data_LungCancer$ID) 

#areas = unique(Data_LungCancer$ID)
#n     = length(areas)

# Create provinces IDs
map_Spain$ID.prov <- substr(map_Spain$ID, 1, 2)

# Ad a area-id starting at 1 and ending at 7906 to both map_Spain and Data_LungCancer
map_Spain$area_id = 1:nrow(map_Spain)
Data_LungCancer$area_id = rep(NA, nrow(Data_LungCancer))
for(year_id in 1:tT){
  Data_LungCancer[Data_LungCancer$year_id == year_id, ]$area_id = map_Spain$area_id
}



# Remove several regions of Spain for this test scenairo
map_Spain <- map_Spain[map_Spain$region %in% c("Extremadura", 
                                               "Andalucia", 
                                               "Murcia"), 
                       ]

# ggplot(data = map_Spain) + 
#   geom_sf(alpha = 1,
#           color="black") +  
#   theme(plot.title = element_text(size = 15,  hjust = 0.5),
#         axis.title.x = element_blank(), #Remove axis and background grid
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_blank(),
#         plot.margin =  unit(c(0, 0, 0, 0), "inches"),
#         legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
#         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
#         panel.spacing = unit(1, 'lines')) +
#   guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) 


#Extract area_ids of remaining areas, and remove areas from Data_LungCancer
remaining_area_ids <- unique(map_Spain$area_id)
Data_LungCancer <- Data_LungCancer[Data_LungCancer$area_id %in% remaining_area_ids, ]


################################################################################




# To get predictions 2 years ahead!
Data_LungCancer[Data_LungCancer$year_id %in% 6:7, ]$obs = NA

#reset row_names
rownames(Data_LungCancer) = 1:nrow(Data_LungCancer)
rownames(map_Spain)       = 1:nrow(map_Spain)

imp_typeIV <- STCAR_INLA(carto    = map_Spain, 
                         data     = Data_LungCancer,
                         ID.area  = "ID", 
                         ID.year  = "year", 
                         O        = "obs", 
                         E        = "pop", 
                         ID.group = "ID.prov",
                         spatial  = "BYM2", 
                         temporal = "rw1", 
                         interaction="TypeIV",
                         model    = "partition", 
                         k        = 2, 
                         scale.model = T,
                         PCpriors = T,
                         merge.strategy = "original",
                         strategy = "gaussian",
                         compute.fitted.values = T) #control.compute=list(return.marginals.predictor=TRUE)

save(imp_typeIV,
     file = "./case_study/imp_typeIV_test.RData")
