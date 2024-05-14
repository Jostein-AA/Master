#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library("geofacet")
library(ggh4x)
library(ggridges)
library(latex2exp)
library(geoR)
library(paletteer)
library('ggsci')
library(ggstats)
library(bigDM)
library(tmap)
library(RColorBrewer)

load("case_study/imp_typeIV.RData")

load("case_study/proper2_onlyInt_Spain.RData")

load("case_study/proper1Int_impEffects_Spain.RData")
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
# plot objects themselves



plot(imp_typeIV)

plot(proper2_onlyInt_Spain)

plot(proper1Int_imp1Effects_Spain)

################################################################################
# Plot the fit on the map

# Get the same paleta
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(seq(min(min(proper2_onlyInt_Spain$summary.fitted.values$mean  * 1E5),
                    min(imp_typeIV$summary.fitted.values$mean  * 1E5)),
                max(max(proper2_onlyInt_Spain$summary.fitted.values$mean  * 1E5),
                    max(imp_typeIV$summary.fitted.values$mean  * 1E5)), 
                length.out = 8),
            Inf)

# Estimated rate proper2_onlyInt_Spain
est_r_prop2_OI <- matrix(proper2_onlyInt_Spain$summary.fitted.values$mean * 1E5, 
                         nrow = n, 
                         ncol = tT, 
                         byrow = F)

colnames(est_r_prop2_OI) = paste("Year", seq(t.from, t.to), sep = ".")

carto_prop2_OI <- cbind(map_Spain, est_r_prop2_OI)


Map.risks <- tm_shape(carto_prop2_OI) +
  tm_polygons(col=paste("Year",round(seq(t.from,t.to,length.out=9)),sep= "."),
              palette=paleta, title="proper2_onlyInt:\n Posterior mean of\n rate per 100,000", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, midpoint=0, interval.closure="left") +
  tm_grid(n.x=5, n.y=5, alpha=0.2, labels.format=list(scientific=T),
          labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="", main.title.position="center", panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01),
            panel.labels=as.character(round(seq(t.from,t.to,length.out=9)))) +
  tm_facets(nrow=3, ncol=3)

print(Map.risks)

# Estimated rate proper2_onlyInt_Spain
est_r_prop1_impEff <- matrix(proper1Int_imp1Effects_Spain$summary.fitted.values$mean * 1E5, 
                         nrow = n, 
                         ncol = tT, 
                         byrow = F)

colnames(est_r_prop1_impEff) = paste("Year", seq(t.from, t.to), sep = ".")

carto_prop1_impEff <- cbind(map_Spain, est_r_prop1_impEff)


Map.risks <- tm_shape(carto_prop1_impEff) +
  tm_polygons(col=paste("Year",round(seq(t.from,t.to,length.out=9)),sep= "."),
              palette=paleta, title="proper1_impEff:\n Posterior mean of\n rate per 100,000", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, midpoint=0, interval.closure="left") +
  tm_grid(n.x=5, n.y=5, alpha=0.2, labels.format=list(scientific=T),
          labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="", main.title.position="center", panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01),
            panel.labels=as.character(round(seq(t.from,t.to,length.out=9)))) +
  tm_facets(nrow=3, ncol=3)

print(Map.risks)



# Estimated rate Improper1\_typeIV
est_r_imp_typeIV <- matrix(imp_typeIV$summary.fitted.values$mean * 1E5, 
                         nrow = n, 
                         ncol = tT, 
                         byrow = F)

colnames(est_r_imp_typeIV) = paste("Year", seq(t.from, t.to), sep = ".")

carto_imp_typeIV <- cbind(map_Spain, est_r_imp_typeIV)


Map.risks <- tm_shape(carto_imp_typeIV) +
  tm_polygons(col=paste("Year",round(seq(t.from,t.to,length.out=9)),sep= "."),
              palette=paleta, title="Improper1_typeIV:\n Posterior mean of\n rate per 100,000", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, midpoint=0, interval.closure="left") +
  tm_grid(n.x=5, n.y=5, alpha=0.2, labels.format=list(scientific=T),
          labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="", main.title.position="center", panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01),
            panel.labels=as.character(round(seq(t.from,t.to,length.out=9)))) +
  tm_facets(nrow=3, ncol=3)

print(Map.risks)



















