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


map  <- read_sf('./Shapefiles/DEU_adm/DEU_2/')

ggplot(data = map) + 
  geom_sf(aes(), 
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
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right"))


#temp
scale_col = heat.colors(30, rev=TRUE)
scale = scale_col[seq(3, 30, length.out = 12)]

plt_true_discrete_rate_four_years(scenario_name = "sc13", 
                                  dataset_id = 1,
                                  admin_map = map,
                                  scale_col = scale_col,
                                  scale = scale)


plt_true_discrete_rate_four_years(scenario_name = "sc14", 
                                  dataset_id = 1,
                                  admin_map = map,
                                  scale_col = scale_col,
                                  scale = scale)

plt_true_discrete_rate_four_years(scenario_name = "sc15", 
                                  dataset_id = 1,
                                  admin_map = map,
                                  scale_col = scale_col,
                                  scale = scale)

plt_true_discrete_rate_four_years(scenario_name = "sc16", 
                                  dataset_id = 1,
                                  admin_map = map,
                                  scale_col = scale_col,
                                  scale = scale)

plt_true_discrete_rate_four_years(scenario_name = "sc17", 
                                  dataset_id = 1,
                                  admin_map = map,
                                  scale_col = scale_col,
                                  scale = scale)

plt_true_discrete_rate_four_years(scenario_name = "sc18", 
                                  dataset_id = 1,
                                  admin_map = map,
                                  scale_col = scale_col,
                                  scale = scale)



