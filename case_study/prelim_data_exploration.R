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

################################################################################

# Load in the considered lung-cancer data
data(Data_LungCancer)

# Load in a map of Spain containing 7907 areas
data("Carto_SpainMUN")

str(Data_LungCancer)
str(Carto_SpainMUN)

# Get the years and areas
years = unique(Data_LungCancer$year)
tT = length(years)
t.from <- min(Data_LungCancer$year)
t.to <- max(Data_LungCancer$year)

areas = unique(Data_LungCancer$ID)
n = length(areas)

# Calculate expected number of counts per 100,000
Data_LungCancer$exp_per_1E5 <- (Data_LungCancer$exp * 1E5)/Data_LungCancer$pop
################################################################################
map_Spain = Carto_SpainMUN
map_Spain$disjoint = "no"
map_Spain[2454, ]$disjoint = "yes"

ggplot(data = map_Spain) + 
  geom_sf(aes(fill = disjoint), 
          alpha = 1,
          color="black") +  
  theme(plot.title = element_text(size = 15,  hjust = 0.5),
        axis.title.x = element_blank(), #Remove axis and background grid
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin =  unit(c(0, 0, 0, 0), "inches"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.spacing = unit(1, 'lines')) +
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + 
  scale_fill_manual(values = c("no" = "lightgrey", "yes" = "red"))


################################################################################
# Plot the expected counts in each and the observed number of counts for select years

#####
# Expected nnumber of counts

expected <- matrix(Data_LungCancer$exp_per_1E5, nrow = n, ncol = tT, byrow = F)
colnames(expected) = paste("Year", seq(t.from, t.to), sep = ".")

carto <- cbind(Carto_SpainMUN, expected)


paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(min(Data_LungCancer$exp_per_1E5), 45, 70, 90,
            110, 140, 170, 200, Inf)

Map.risks <- tm_shape(carto) +
  tm_polygons(col=paste("Year",round(seq(t.from,t.to,length.out=9)),sep= "."),
              palette=paleta, title="Expected count\n per 100,000", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, midpoint=0, interval.closure="left") +
  tm_grid(n.x=5, n.y=5, alpha=0.2, labels.format=list(scientific=T),
          labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="", main.title.position="center", panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01),
            panel.labels=as.character(round(seq(t.from,t.to,length.out=9)))) +
  tm_facets(nrow=3, ncol=3)

print(Map.risks)

#####
# Observed number of counts


observed <- matrix(Data_LungCancer$obs, nrow = n, ncol = tT, byrow = F)
colnames(observed) = paste("Year", seq(t.from, t.to), sep = ".")

carto <- cbind(Carto_SpainMUN, observed)


paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(min(Data_LungCancer$obs), 10, 20, 40,
            60, 90, 150, 200, Inf)

Map.risks <- tm_shape(carto) +
  tm_polygons(col=paste("Year",round(seq(t.from,t.to,length.out=9)),sep= "."),
              palette=paleta, title="Observed counts", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, midpoint=0, interval.closure="left") +
  tm_grid(n.x=5, n.y=5, alpha=0.2, labels.format=list(scientific=T),
          labels.inside.frame=F, labels.col="white") +
  tm_layout(main.title="", main.title.position="center", panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01),
            panel.labels=as.character(round(seq(t.from,t.to,length.out=9)))) +
  tm_facets(nrow=3, ncol=3)

print(Map.risks)









