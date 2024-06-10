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
library(tidyverse)
library(bigDM)
library(tables)
library(tmap)
library(RColorBrewer)

load("maps_and_nb.RData")
load("grids_and_mappings.RData")

# Some useful numbers
tT = 13
n_sim = 100
E_it = 100
n_ADM1 <- nrow(first_level_admin_map)
n_ADMnew <- nrow(new_map)
n_ADM4 <- nrow(second_level_admin_map)

dataset_id = 3 
dataset_id_2 = 4
dataset_id_new = 2


scenario_names_ADM4 <- c("sc2", "sc4", "sc6",
                         "sc8", "sc10", "sc12")

scenario_names_ADM1 <- c("sc1", "sc3", "sc5",
                         "sc7", "sc9", "sc11")


############################################################################################################################################

plt_overall_results_all_scenarios <- function(model_names, 
                                              scenario_names,
                                              xlab,
                                              title){
  
  scenario_name = scenario_names[1]
  model_name = model_names[1]
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_rates
  
  values.df <- data.frame(model_name = model_name,
                           scenario_name = scenario_name,
                           IS_1_year_ahead = mean(tmp_$IS_1_year_ahead, na.rm = TRUE),
                           IS_2_year_ahead = mean(tmp_$IS_2_year_ahead, na.rm = TRUE),
                           IS_3_year_ahead = mean(tmp_$IS_3_year_ahead, na.rm = TRUE),
                           total_IS = mean(tmp_$total_IS, na.rm = TRUE))
  
  for(scenario_name in scenario_names){
    for(model_name in model_names){
      if(scenario_name == scenario_names[1] & model_name == model_names[1]){
        # Skip this one as it is all ready done
        next
      }
      
      load(paste("./results/model_choice/model_choice_", 
                 model_name, "_", 
                 scenario_name, ".RData",
                 sep = ""))
      
      tmp_ <- model_choice_for_rates
      
      tmp2_ <- data.frame(model_name = model_name,
                          scenario_name = scenario_name,
                          IS_1_year_ahead = mean(tmp_$IS_1_year_ahead, na.rm = TRUE),
                          IS_2_year_ahead = mean(tmp_$IS_2_year_ahead, na.rm = TRUE),
                          IS_3_year_ahead = mean(tmp_$IS_3_year_ahead, na.rm = TRUE),
                          total_IS = mean(tmp_$total_IS, na.rm = TRUE))
      
      
      values.df = rbind(values.df, tmp2_)
    }
  }
  
  # Average over scenarios
  to_plot.df <- data.frame(model_name = rep(model_names, 4),
                           div = c(rep("1 year ahead", length(model_names)),
                                   rep("2 years ahead", length(model_names)),
                                   rep("3 years ahead", length(model_names)),
                                   rep("Overall", length(model_names))),
                           values = 1:(length(model_names) * 4))
                           
                           
                           
  
  for(model_name in model_names){
    to_plot.df[to_plot.df$model_name == model_name & 
                 to_plot.df$div == "1 year ahead", ]$values <- mean(as.numeric(values.df[values.df$model_name == model_name, ]$IS_1_year_ahead))
    to_plot.df[to_plot.df$model_name == model_name & 
                 to_plot.df$div == "2 years ahead", ]$values <- mean(as.numeric(values.df[values.df$model_name == model_name, ]$IS_2_year_ahead))
    to_plot.df[to_plot.df$model_name == model_name & 
                 to_plot.df$div == "3 years ahead", ]$values <- mean(as.numeric(values.df[values.df$model_name == model_name, ]$IS_3_year_ahead))
    to_plot.df[to_plot.df$model_name == model_name & 
                 to_plot.df$div == "Overall", ]$values <- mean(as.numeric(values.df[values.df$model_name == model_name, ]$total_IS))
  }
  
  to_plot.df[to_plot.df$model_name == "Improper1_noInt", ]$model_name = "No interaction"
  to_plot.df[to_plot.df$model_name == "Improper1_typeI", ]$model_name = "Knorr-Held I"
  to_plot.df[to_plot.df$model_name == "Improper1_typeIV", ]$model_name = "Knorr-Held IV"
  to_plot.df[to_plot.df$model_name == "proper2_propInt_Improp_temporal", ]$model_name = "proper interaction"
  
  #level_order = scenario_names
  
  to_plot.df <- to_plot.df %>% 
                mutate(model_name = fct_reorder(.f = model_name,
                                                .x = values,
                                                .fun = mean,
                                                .na_rm = TRUE,
                                                .desc = TRUE))
  print(to_plot.df$model_name)
  
  
  return(
    # to_plot.df %>% 
    # mutate(model_name = fct_reorder(.f = model_name,
    #                                 .x = values,
    #                                 .fun = mean,
    #                                 .na_rm = TRUE,
    #                                 .desc = TRUE)) %>%
    ggplot(data = to_plot.df,
           aes(x = div, 
               y = values,
               fill = model_name)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 18.5),
          axis.title = element_text(size = 18.5),
          legend.title = element_text(size = 18.5),
          legend.text = element_text(size = 16),
          legend.position = "top",
          plot.title = element_text(vjust = -7, size = 18.5)) + 
    ggtitle(title) + guides(fill=guide_legend(title="Model:")) + 
    xlab(NULL) + 
    ylab("Average IS"))
}





model_names = c("Improper1_noInt",
                "Improper1_typeI", 
                "Improper1_typeIV",
                "proper2_propInt_Improp_temporal")


# Save as AAA_res_poster 11 by 4
plt_overall_results_all_scenarios(model_names,
                                  scenario_names_ADM4,
                                  xlab = "",
                                  title = "n = 402.")

plt_ADM1 <- plt_overall_results_all_scenarios(model_names,
                                              scenario_names_ADM1,
                                              xlab = "",
                                              title = "n = 16")

ggarrange(plt_ADM1, plt_ADM4,
          ncol = 2, nrow = 1,
          common.legend = T, legend = "right")


plt_overall_results_all_scenarios_2 <- function(model_names, 
                                              scenario_names,
                                              xlab,
                                              title){
  
  scenario_name = scenario_names[1]
  model_name = model_names[1]
  load(paste("./results/model_choice/model_choice_", 
             model_name, "_", 
             scenario_name, ".RData",
             sep = ""))
  
  tmp_ <- model_choice_for_rates
  
  values.df <- data.frame(model_name = model_name,
                          scenario_name = scenario_name,
                          IS_1_year_ahead = mean(tmp_$IS_1_year_ahead, na.rm = TRUE),
                          IS_2_year_ahead = mean(tmp_$IS_2_year_ahead, na.rm = TRUE),
                          IS_3_year_ahead = mean(tmp_$IS_3_year_ahead, na.rm = TRUE),
                          total_IS = mean(tmp_$total_IS, na.rm = TRUE))
  
  for(scenario_name in scenario_names){
    for(model_name in model_names){
      if(scenario_name == scenario_names[1] & model_name == model_names[1]){
        # Skip this one as it is all ready done
        next
      }
      
      load(paste("./results/model_choice/model_choice_", 
                 model_name, "_", 
                 scenario_name, ".RData",
                 sep = ""))
      
      tmp_ <- model_choice_for_rates
      
      tmp2_ <- data.frame(model_name = model_name,
                          scenario_name = scenario_name,
                          IS_1_year_ahead = mean(tmp_$IS_1_year_ahead, na.rm = TRUE),
                          IS_2_year_ahead = mean(tmp_$IS_2_year_ahead, na.rm = TRUE),
                          IS_3_year_ahead = mean(tmp_$IS_3_year_ahead, na.rm = TRUE),
                          total_IS = mean(tmp_$total_IS, na.rm = TRUE))
      
      
      values.df = rbind(values.df, tmp2_)
    }
  }
  
  #to_plot.df <- data.frame(model_name = rep(model_names, 4),
  #                         div = c(rep("1 year ahead", length(model_names)),
  #                                 rep("2 years ahead", length(model_names)),
  #                                 rep("3 years ahead", length(model_names)),
  #                                 rep("Overall", length(model_names))),
  #                         values = 1:(length(model_names) * 4))
  
  to_plot.df <- data.frame(model_name = rep(model_names, 4 * 600),
                          div = c(rep("1 year ahead", length(model_names) * 600),
                                  rep("2 years ahead", length(model_names) * 600),
                                  rep("3 years ahead", length(model_names) * 600),
                                  rep("Overall", length(model_names) * 600)),
                          values = 1:(length(model_names) * 4 * 600))
  
  
  
  
  for(model_name in model_names){
    
    to_plot.df[to_plot.df$model_name == model_name &
                 to_plot.df$div == "1 year ahead", ]$values <- values.df[values.df$model_name == model_name, ]$IS_1_year_ahead
    to_plot.df[to_plot.df$model_name == model_name &
                 to_plot.df$div == "2 years ahead", ]$values <- values.df[values.df$model_name == model_name, ]$IS_2_year_ahead
    to_plot.df[to_plot.df$model_name == model_name &
                 to_plot.df$div == "3 years ahead", ]$values <- values.df[values.df$model_name == model_name, ]$IS_3_year_ahead
    to_plot.df[to_plot.df$model_name == model_name &
                 to_plot.df$div == "Overall", ]$values <- values.df[values.df$model_name == model_name, ]$total_IS
  
    
    }
  
  to_plot.df[to_plot.df$model_name == "Improper1_noInt", ]$model_name = "No interaction"
  to_plot.df[to_plot.df$model_name == "Improper1_typeI", ]$model_name = "Knorr-Held I"
  to_plot.df[to_plot.df$model_name == "Improper1_typeIV", ]$model_name = "Knorr-Held IV"
  to_plot.df[to_plot.df$model_name == "proper2_propInt_Improp_temporal", ]$model_name = "proper interaction"
  
  return(to_plot.df <- to_plot.df %>% 
           mutate(model_name = fct_reorder(.f = model_name,
                                           .x = values,
                                           .fun = mean,
                                           .na_rm = TRUE,
                                           .desc = TRUE)) %>%
    ggplot(aes(x = div, 
               y = values,
               fill = model_name)) + 
      geom_boxplot() + #position="dodge", stat="identity"
      #geom_stripped_cols(odd = "#FFFFFF00", even = "#00000001", alpha = .08) +
      theme_bw() + 
      theme(axis.text = element_text(size = 18.5),
            axis.title = element_text(size = 18.5),
            legend.title = element_text(size = 18.5),
            legend.text = element_text(size = 16),
            legend.position = "top",
            plot.title = element_text(vjust = -7, size = 20.5, face="bold")) + 
      ggtitle(title) + guides(fill=guide_legend(title="Model:")) + 
      xlab(NULL) + 
      ylab("Average IS"))
}

# Save as AAA_res_poster 11 by 4
plt_overall_results_all_scenarios_2(model_names,
                                  scenario_names_ADM4,
                                  xlab = "",
                                  title = "n = 402.")



#######################################################################################################
# Plot the temporal trends

# For presentation


temporal_pattern.df <- data.frame(year = rep(1:13, 3),
                                  value = c(rep(exp(log(0.1)), 13) * 100,
                                            exp(log(0.1) + 0.014 * 1:13) * 100,
                                            c(exp(log(0.1) + 0.02 * 1:7) * 100, 
                                              exp(log(0.1) + (0.02 * 7) - 0.015 * (8-7):(13-7)) * 100)),
                                  trend = c(rep("Constant", 13),
                                            rep("Linear", 13),
                                            rep("Changing", 13)))


# Save as temporal_trends_poster 5 by 5
ggplot(data = temporal_pattern.df) + 
  geom_line(aes(x = year, y = value, color = trend)) + 
  ylim(9, 12.5) + 
  ylab("Rate per 100") + 
  theme_bw() + 
  theme(legend.box.background = element_rect(color="grey", size=1),
        axis.title=element_text(size=19.5),
        strip.text.x = element_text(size=19.5),
        axis.text=element_text(size=19.5),
        legend.position = c(0.25, .85),
        legend.title = element_text(size=19.5),
        legend.text = element_text(size=19.5),
        plot.title = element_text(size=19.5)) +
  scale_x_continuous(name = "Year", 
                    breaks = c(1, 4, 7, 10, 13)) + 
  guides(color = guide_legend(title = NULL)) + 
  ggtitle("Temporal trends")
  


###########################################################################################################
# Plot spatial range

# Extract two risk-surfaces for single time-point

load("Data/Simulated_risk_surfaces/sc2_risk_surfaces.RData")
risk_surface.list_sc2_t1 <- risk_surface.list[risk_surface.list$t == t_axis[1], ]
risk_surface.list_sc2_t1$values = risk_surface.list[risk_surface.list$t == t_axis[1], ]$values[, 1]
risk_surface.list_sc2_t1$values = log(risk_surface.list_sc2_t1$values) - log(0.1)
rm(risk_surface.list)

load("Data/Simulated_risk_surfaces/sc8_risk_surfaces.RData")
risk_surface.list_sc8_t1 <- risk_surface.list[risk_surface.list$t == t_axis[1], ]
risk_surface.list_sc8_t1$values = risk_surface.list[risk_surface.list$t == t_axis[1], ]$values[, 1]
risk_surface.list_sc8_t1$values = log(risk_surface.list_sc8_t1$values) - log(0.1)
rm(risk_surface.list)

coords = matrix(0, nrow = nrow(risk_surface.list_sc8_t1), ncol = 2)
coords[, 1] = risk_surface.list_sc8_t1$x; coords[, 2] = risk_surface.list_sc8_t1$y

#Standardize coords
coords <- apply(coords, 2, function(x) (x - min(x)) / (max(x) - min(x)))

risk_surface.list_sc2_t1$coords = coords
risk_surface.list_sc8_t1$coords = coords

# Sample a set of indices for the risk_surfaces
indices_to_use <- sample(1:nrow(risk_surface.list_sc8_t1), 1000, replace = F)
#indices_to_use = 1:nrow(risk_surface.list_sc8_t1)

ML_sc2_t1 <- likfit(#geodata = risk_surface.list_sc2_t1[indices_to_use, ],
                    coords = risk_surface.list_sc2_t1[indices_to_use, ]$coords,
                    data = risk_surface.list_sc2_t1[indices_to_use, ]$values,
                    ini.cov.pars = c(.1, 1),
                    messages = TRUE,
                    fix.nugget = TRUE)

ML_sc2_t1_est_sigma <- ML_sc2_t1$parameters.summary$values[3]
ML_sc2_t1_est_range <- ML_sc2_t1$parameters.summary$values[4]



temp_cov_sc2 = cov.spatial(seq(0, 10, length.out = 100), 
                           cov.model = "exponential",#Calculate the covariance based on the distances
                           cov.pars = c(1, 
                                        ML_sc2_t1_est_range))


ML_sc8_t1 <- likfit(#geodata = risk_surface.list_sc8_t1[indices_to_use, ],
                    coords = risk_surface.list_sc8_t1[indices_to_use, ]$coords,
                    data = risk_surface.list_sc8_t1[indices_to_use, ]$values,
                    ini.cov.pars = c(0.1, 1),
                    messages = TRUE,
                    fix.nugget = TRUE)


ML_sc8_t1_est_sigma <- ML_sc8_t1$parameters.summary$values[3]
ML_sc8_t1_est_range <- ML_sc8_t1$parameters.summary$values[4]


temp_cov_sc8 = cov.spatial(seq(0, 10, length.out = 100), 
                           cov.model = 'exponential', #Calculate the correlation based on the distances
                           cov.pars = c(1, 
                                        ML_sc8_t1_est_range))

spatial_correlation.df = data.frame(distances = c(seq(0, 10, length.out = 100),
                                                  seq(0, 10, length.out = 100)),
                                    corrs = c(temp_cov_sc2,
                                              temp_cov_sc8),
                                    type = c(rep("Short", 100),
                                             rep("Long", 100)))

# For presentation: spatial_range_poster 5 by 5
ggplot(data = spatial_correlation.df, aes(distances, corrs)) + 
  geom_line(aes(col = type)) + 
  guides(col=guide_legend(title="Spatial range:", 
                          label.position = "right")) +
  theme_bw() + 
  theme(legend.position="top",
        axis.title=element_text(size=19.5),
        legend.text = element_text(size=19.5),
        legend.title = element_text(size=19.5),
        axis.text=element_text(size=17.5)) +
  ylab(TeX(r'(Corr$[X+h,\;X]$)')) + #"Correlation"
  xlab(TeX(r'(h)')) #Distance


###########################################################################################
# Plot for the case study



load("case_study/proper2_RW1_Spain_del_Extremadura.RData")
load("case_study/imp_typeIV_Spain_del_Extremadura.RData")

# Load in the considered lung-cancer data
data(Data_LungCancer)

# Load in a map of Spain containing 7907 areas
data("Carto_SpainMUN")

# Remove disjointed area
map_Spain <- Carto_SpainMUN[-2454, ]
row.names(map_Spain) = 1:nrow(map_Spain)
problem_area = Carto_SpainMUN[2454, ]

Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID != problem_area$ID, ]

# Extract the areas NOT within the principality of Extremadura
map_Spain <- map_Spain[map_Spain$region != "Extremadura", ]
IDs_NOT_within_Extremadura <- unique(map_Spain$ID)

# Extract the data NOT within Extremadura
Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID %in% IDs_NOT_within_Extremadura, ]

# Get the years and areas
years = unique(Data_LungCancer$year)
tT = length(years)
t.from <- min(Data_LungCancer$year)
t.to <- max(Data_LungCancer$year)

# Create time-ids 1,...,25 instead of 1991,...,2015 for the sake of INLA
Data_LungCancer$year_id <- Data_LungCancer$year - min(Data_LungCancer$year) + 1

areas = unique(Data_LungCancer$ID)
n = length(areas)

# Create provinces IDs
map_Spain$ID.prov <- substr(map_Spain$ID, 1, 2)

# Ad a area-id starting at 1 and ending at 7906 to both map_Spain and Data_LungCancer
map_Spain$area_id = 1:nrow(map_Spain)
Data_LungCancer$area_id = rep(NA, nrow(Data_LungCancer))
for(year_id in 1:tT){
  Data_LungCancer[Data_LungCancer$year_id == year_id, ]$area_id = map_Spain$area_id
}

# Just for the sake of simplicity, sort the marginals of the proper models now
if(FALSE){ # if(FALSE) added so that dont sort them unless I really want to. Need sorting only once
  proper2_RW1_Spain_del_Extremadura$marginals.fitted.values <- sort_proper_fitted(proper2_RW1_Spain_del_Extremadura$marginals.fitted.values,
                                                                                  n, tT)
  
  proper2_RW1_Spain_del_Extremadura$summary.fitted.values$mean <- sort_proper_fitted(proper2_RW1_Spain_del_Extremadura$summary.fitted.values$mean,
                                                                                     n, tT)
  
  proper2_RW1_Spain_del_Extremadura$summary.fitted.values$sd <- sort_proper_fitted(proper2_RW1_Spain_del_Extremadura$summary.fitted.values$sd,
                                                                                   n, tT)
  
}

# Predicted SD of counts per 100,000
sd_2015 <- get_pred_SD(proper2_RW1_Spain_del_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 25)

#pred_SD_prop2_RW1 <- matrix(c(sd_2011, sd_2013, sd_2015), nrow = n, ncol = 3, byrow = F)
#colnames(pred_SD_prop2_RW1) = paste("Year", c(2011, 2013, 2015), sep = ".")


# Predicted SD of counts per 100,000
sd_IV_2015 <- get_pred_SD(imp_typeIV$marginals.fitted.values, Data_LungCancer$pop, n, 25)

#pred_SD_impIV <- matrix(c(sd_IV_2011, sd_IV_2013, sd_IV_2015), nrow = n, ncol = 3, byrow = F)
#colnames(pred_SD_impIV) = paste("Year", c(2011, 2013, 2015), sep = ".")

pred_SD <- matrix(c(sd_2015, sd_IV_2015), nrow = n, ncol = 2, byrow = F)
colnames(pred_SD) = c("Proper", "Improper")


# Create a common color-palet
paleta <- brewer.pal(8,"RdYlGn")[8:1]
pred_count_qs <- quantile(pred_SD, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
values <- c(round(min(pred_SD), 1), round(pred_count_qs, 1), Inf)


# Plot the predicted SD of counts per 100,000
carto_prop2_RW1 <- cbind(map_Spain, pred_SD)

Map.risks_prop2_RW1 <- tm_shape(carto_prop2_RW1) +
  tm_polygons(col=c("Proper", "Improper"),
              palette=paleta, 
              title = "", #"SD of the posterior\npredicted counts\nper 100,000" 
              legend.show=T, 
              border.col="transparent",
              legend.reverse=T, 
              style="fixed", 
              breaks=values, 
              midpoint=0, 
              interval.closure="left") +
  tm_grid(n.x=5, 
          n.y=5, 
          alpha=0.2, 
          labels.format=list(scientific=T),
          labels.inside.frame=F, 
          labels.col="white") +
  tm_layout(main.title="SD of the posterior predicted counts per 100,000", 
            main.title.position="center",
            main.title.size = 1.6,
            bg.color = "white", # Background color, white
            outer.bg.color = "white", 
            panel.label.size=1.5,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.text.size = 1.3, 
            legend.outside.size=0.175, #0.15
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0, 0, 0, 0), #0.01, 0.01, 0.01, 0.01
            between.margin = 0, #0.01
            panel.labels=as.character(c("Proper interaction: 2015", "Divide-and-conquer: 2015"))) +
  tm_facets(nrow=1, ncol=2)


print(Map.risks_prop2_RW1)

# Save Map
tmap_save(Map.risks_prop2_RW1,
          filename = "Plots/case_study_poster.pdf",
          width = 10,
          height = 4)











# Predicted number of counts per 100,000
pred_count <- matrix(c(proper2_RW1_Spain_del_Extremadura$summary.fitted.values$mean[(n * (25 - 1) + 1):(25 * n)] * 1E5, #
                       imp_typeIV$summary.fitted.values$mean[(n * (25 - 1) + 1):(25 * n)] * 1E5), 
                     nrow = n, ncol = 2, byrow = F)
colnames(pred_count) = c("Proper", "Improper")  #paste("Year", seq(t.from, t.to), sep = ".")


# Create a common color-palet
paleta <- brewer.pal(8,"RdYlGn")[8:1]
pred_count_qs <- quantile(pred_count, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
values <- c(min(pred_count), pred_count_qs, Inf)



# Plot the predicted counts per 100,000 for the proper2_RW1
carto_prop2_RW1 <- cbind(map_Spain, pred_count)


# title="Predicted count\n per 100,000\n proper2_RW1",


Map.risks_prop2_RW1 <- tm_shape(carto_prop2_RW1) +
  tm_polygons(col=c("Proper", "Improper"),
              palette=paleta, 
              title = "", #"Mean posterior\npredicted counts\nper 100,000" 
              legend.show=T, 
              border.col="transparent",
              legend.reverse=T, 
              style="fixed", 
              breaks=values, 
              midpoint=0, 
              interval.closure="left") +
  tm_grid(n.x=5, 
          n.y=5, 
          alpha=0.2, 
          labels.format=list(scientific=T),
          labels.inside.frame=F, 
          labels.col="white") +
  tm_layout(main.title="Mean posterior predicted counts per 100,000", 
            main.title.position="center",
            main.title.size = 1.6,
            bg.color = "white", # Background color, white
            outer.bg.color = "white", 
            panel.label.size=1.5,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.text.size = 1.3, 
            legend.outside.size=0.175, #0.15
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0, 0, 0, 0), #0.01, 0.01, 0.01, 0.01
            between.margin = 0, #0.01
            panel.labels=as.character(c("Proper interaction: 2015", "Divide-and-conquer: 2015"))) +
  tm_facets(nrow=1, ncol=2)

print(Map.risks_prop2_RW1)

tmap_save(Map.risks_prop2_RW1,
          filename = "Plots/case_study_poster_2.pdf",
          width = 10,
          height = 4)


################################################################################

























