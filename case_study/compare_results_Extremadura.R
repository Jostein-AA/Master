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

load("case_study/proper2_RW1_Extremadura.RData")
load("case_study/proper2_impEff_Extremadura.RData")
load("case_study/Improper1_typeIV_Extremadura.RData")

################################################################################
plot(proper2_RW1_Extremadura)

plot(proper2_impEff_Extremadura)

plot(Improper1_typeIV_Extremadura)

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



################################################################################
# Calculate model choice

## Calculate model choice critera: MSE, IS, log-score
calc_model_choice <- function(model,
                              data,
                              n, tT,
                              Improper = T){
  
  ### Initialize data frame to store the MSE's, IS's, and log-scores 1, 2, 3, 4, and 5 years ahead and average
  model_choice.df <- data.frame(mse_1 = NA, mse_2 = NA, mse_3 = NA,
                                mse_4 = NA, mse_5 = NA, mse_avg = NA,
                                is_1 = NA, is_2 = NA, is_3 = NA,
                                is_4 = NA, is_5 = NA, is_avg = NA,
                                log_s_1 = NA, log_s_2 = NA, log_s_3 = NA,
                                log_s_4 = NA, log_s_5 = NA, log_s_avg = NA)
  
  
  # If proper models, sort the marginals to get the correct order
  if(!Improper){
    model$marginals.fitted.values <- sort_proper_fitted(model$marginals.fitted.values,
                                                        n, tT)
  }
  
  ### Take the years predicted on
  years_pred_on <- 21:25
  
    
  #### For each year calculate the MSE and IS that year
  for(year in years_pred_on){
    print(paste("Year: ", year, sep = ""))
    
    ## Extract the predicted marginals for this year
    one_year_margs <- model$marginals.fitted.values[((year - 1) * n + 1):(year * n)]
    
    ## Extract the observed counts for this year
    one_year_counts = data$obs[((year - 1) * n + 1):(year * n)]
    
    ## Extract the populations for this year
    one_year_pop <- data$pop[((year - 1) * n + 1):(year * n)]
    
    
    ## MSE
    model_choice.df[1, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(one_year_counts,
                                                                                      one_year_margs,
                                                                                      E_it = one_year_pop)
    
    
    ## IS for count in year ahead 
    model_choice.df[1, year - years_pred_on[1] + 7] =  count_IS_one_year_case_study(one_year_counts, 
                                                                                     one_year_margs, 
                                                                                    population = one_year_pop)
    
    
    
    
    ## Log-score
    #model_choice.df[1, year - years_pred_on[1] + 13] = count_log_s_one_year(one_year_counts,
    #                                                                        one_year_margs,
    #                                                                        population = one_year_pop)
    
    
  }
    
    #Get the total MSE for count
  model_choice.df[1, 6] = mean(as.numeric(model_choice.df[1, 1:5]))
    
    
    #Get the total IS for count
  model_choice.df[1, 12] = mean(as.numeric(model_choice.df[1, 7:11]))
    
  
  
  #filename = paste("./results/model_choice/", "model_choice_", model_name, "_", 
  #                 scenario_name, ".RData", sep = "")
  
  #save(,
  #     file = filename)
  return(model_choice.df)
}


tmp <- calc_model_choice(proper2_RW1_Extremadura,
                  Data_LungCancer,
                  n, tT, 
                  Improper = F)


tmp2 <- calc_model_choice(Improper1_typeIV_Extremadura,
                  Data_LungCancer,
                  n, tT, 
                  Improper = T)

tmp3 <- calc_model_choice(proper2_impEff_Extremadura,
                  Data_LungCancer,
                  n, tT, 
                  Improper = F)

tmp.df <- data.frame(model = c(rep("proper2_RW1", 12),
                               rep("proper2_impEff", 12),
                               rep("Improper1_typeIV", 12)),
                     model_choice = rep(1:12, 3),
                     value = 1:(12*3))
for(i in 1:12){
  tmp.df[tmp.df$model == "proper2_RW1", ]$value[i] <- tmp[1, i]
  tmp.df[tmp.df$model == "proper2_impEff", ]$value[i] <- tmp2[1, i]
  tmp.df[tmp.df$model == "Improper1_typeIV", ]$value[i] <- tmp3[1, i]
}

### Create a table
#Make caption and label for latex table
caption = "HEIHEI"
label = "model choice Extremadura"

#Make latex table
latex_tabular <- latexTable(tabular(
  Heading("Model")*RowFactor(model, 
                             nopagebreak = "\\hline",
                             spacing = 0)~
    Heading()*Factor(model_choice, 
                     levelnames = c("MSE (1)", "MSE (2)", "MSE (3)", "MSE (4)", "MSE (5)", "MSE (total)",
                                    "IS (1)", "IS (2)", "IS (3)", "IS (4)", "IS (5)", "IS (total)"))*
    Heading()*value*Heading()*identity,
  data = tmp.df),
  caption = caption,
  label = label
)

#Save latex table
cat(latex_tabular, file = "./case_study/Extremadura_model_choice.tex")

  


################################################################################

# Just for the sake of simplicity, sort the marginals of the proper models now
if(FALSE){ # if(FALSE) added so that dont sort them unless I really want to
  proper2_RW1_Extremadura$marginals.fitted.values <- sort_proper_fitted(proper2_RW1_Extremadura$marginals.fitted.values,
                                                                        n, tT)
  
  proper2_RW1_Extremadura$summary.fitted.values$mean <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$mean,
                                                                           n, tT)
  
  proper2_RW1_Extremadura$summary.fitted.values$sd <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$sd,
                                                                         n, tT)
  
  proper2_impEff_Extremadura$marginals.fitted.values <- sort_proper_fitted(proper2_impEff_Extremadura$marginals.fitted.values,
                                                                           n, tT)
  
  proper2_impEff_Extremadura$summary.fitted.values$mean <- sort_proper_fitted(proper2_impEff_Extremadura$summary.fitted.values$mean,
                                                                           n, tT)
  
  proper2_impEff_Extremadura$summary.fitted.values$sd <- sort_proper_fitted(proper2_impEff_Extremadura$summary.fitted.values$sd,
                                                                         n, tT)
}

#########
# Mean pred count

# Predicted number of counts per 100,000
pred_count_prop2_RW1 <- matrix(proper2_RW1_Extremadura$summary.fitted.values$mean * 1E5, nrow = n, ncol = tT, byrow = F)
colnames(pred_count_prop2_RW1) = paste("Year", seq(t.from, t.to), sep = ".")


# Predicted number of counts per 100,000
pred_count_impIV <- matrix(Improper1_typeIV_Extremadura$summary.fitted.values$mean * 1E5, 
                           nrow = n, ncol = tT, byrow = F)
colnames(pred_count_impIV) = paste("Year", seq(t.from, t.to), sep = ".")


# Create a common color-palet
paleta <- brewer.pal(8,"RdYlGn")[8:1]
pred_count_qs <- quantile(pred_count_impIV, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
values <- c(min(pred_count_impIV), pred_count_qs, Inf)



# Plot the predicted counts per 100,000 for the proper2_RW1
carto_prop2_RW1 <- cbind(map_Spain, pred_count_prop2_RW1)

# paleta <- brewer.pal(8,"RdYlGn")[8:1]
# pred_count_qs <- quantile(pred_count_prop2_RW1, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
# values <- c(min(pred_count_prop2_RW1), pred_count_qs, Inf)


Map.risks_prop2_RW1 <- tm_shape(carto_prop2_RW1) +
  tm_polygons(col=paste("Year", c(2011, 2013, 2015),sep= "."),
              palette=paleta, 
              title="Predicted count\n per 100,000\n proper2_RW1", 
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
  tm_layout(main.title="", 
            main.title.position="center",
            bg.color = "white", # Background color, white
            outer.bg.color = "white", 
            panel.label.size=1.5,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.outside.size=0.15, 
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0.01, 0.01, 0.01, 0.01),
            between.margin = 0.01,
            panel.labels=as.character(c(2011, 2013, 2015))) +
  tm_facets(nrow=1, ncol=3)


print(Map.risks_prop2_RW1)

# Save Map
tmap_save(Map.risks_prop2_RW1,
          filename = "Plots/Extremadura_pred_counts_prop2_RW1.pdf",
          width = 12,
          height = 4)





# Plot the predicted counts Improper1_typeIV
carto_impIV <- cbind(map_Spain, pred_count_impIV)

Map.risks_impIV <- tm_shape(carto_impIV) +
  tm_polygons(col=paste("Year", c(2011, 2013, 2015),sep= "."),
              palette=paleta, 
              title="Predicted count\n per 100,000\n Improper1_typeIV", 
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
  tm_layout(main.title="", 
            main.title.position="center",
            bg.color = "white", # Background color, white
            outer.bg.color = "white", 
            panel.label.size=1.5,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.outside.size=0.15, 
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0.01, 0.01, 0.01, 0.01),
            between.margin = 0.01,
            panel.labels=as.character(c(2011, 2013, 2015))) +
  tm_facets(nrow=1, ncol=3)


print(Map.risks_impIV)


# Save Map
tmap_save(Map.risks_impIV,
          filename = "Plots/Extremadura_pred_counts_Improper1_typeIV.pdf",
          width = 12, #12 by 4
          height = 4)


########
# SD count




# Predicted SD of counts per 100,000
sd_2011 <- get_pred_SD(proper2_RW1_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 21)
sd_2013 <- get_pred_SD(proper2_RW1_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 23)
sd_2015 <- get_pred_SD(proper2_RW1_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 25)

pred_SD_prop2_RW1 <- matrix(c(sd_2011, sd_2013, sd_2015), nrow = n, ncol = 3, byrow = F)
colnames(pred_SD_prop2_RW1) = paste("Year", c(2011, 2013, 2015), sep = ".")


# Predicted SD of counts per 100,000
sd_IV_2011 <- get_pred_SD(Improper1_typeIV_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 21)
sd_IV_2013 <- get_pred_SD(Improper1_typeIV_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 23)
sd_IV_2015 <- get_pred_SD(Improper1_typeIV_Extremadura$marginals.fitted.values, Data_LungCancer$pop, n, 25)

pred_SD_impIV <- matrix(c(sd_IV_2011, sd_IV_2013, sd_IV_2015), nrow = n, ncol = 3, byrow = F)
colnames(pred_SD_impIV) = paste("Year", c(2011, 2013, 2015), sep = ".")


# Create a common color-palet
paleta <- brewer.pal(8,"RdYlGn")[8:1]
pred_count_qs <- quantile(pred_SD_impIV, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
values <- c(min(pred_SD_impIV), pred_count_qs, Inf)



# Plot the predicted SD of counts per 100,000 for the proper2_RW1
carto_prop2_RW1 <- cbind(map_Spain, pred_SD_prop2_RW1)

# paleta <- brewer.pal(8,"RdYlGn")[8:1]
# pred_count_qs <- quantile(pred_count_prop2_RW1, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
# values <- c(min(pred_count_prop2_RW1), pred_count_qs, Inf)


Map.risks_prop2_RW1 <- tm_shape(carto_prop2_RW1) +
  tm_polygons(col=paste("Year", c(2011, 2013, 2015),sep= "."),
              palette=paleta, 
              title="Posterior SD of the \npredicted counts\nper 100,000\nproper2_RW1", 
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
  tm_layout(main.title="", 
            main.title.position="center",
            bg.color = "white", # Background color, white
            outer.bg.color = "white", 
            panel.label.size=1.5,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.outside.size=0.15, 
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0.01, 0.01, 0.01, 0.01),
            between.margin = 0.01,
            panel.labels=as.character(c(2011, 2013, 2015))) +
  tm_facets(nrow=1, ncol=3)


print(Map.risks_prop2_RW1)

# Save Map
tmap_save(Map.risks_prop2_RW1,
          filename = "Plots/Extremadura_SD_counts_prop2_RW1.pdf",
          width = 12,
          height = 4)





# Plot the predicted counts Improper1_typeIV
carto_impIV <- cbind(map_Spain, pred_SD_impIV)

Map.risks_impIV <- tm_shape(carto_impIV) +
  tm_polygons(col=paste("Year", c(2011, 2013, 2015),sep= "."),
              palette=paleta, 
              title="Posterior SD of the \npredicted counts\nper 100,000\n Improper1_typeIV", 
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
  tm_layout(main.title="", 
            main.title.position="center",
            bg.color = "white", # Background color, white
            outer.bg.color = "white", 
            panel.label.size=1.5,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.outside.size=0.15, 
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0.01, 0.01, 0.01, 0.01),
            between.margin = 0.01,
            panel.labels=as.character(c(2011, 2013, 2015))) +
  tm_facets(nrow=1, ncol=3)


print(Map.risks_impIV)


# Save Map
tmap_save(Map.risks_impIV,
          filename = "Plots/Extremadura_SD_counts_Improper1_typeIV.pdf",
          width = 12, #12 by 4
          height = 4)

################################################################################
# Calculate widths of 95% CIs for the predicted counts, and also count number of
# Times observations land outside the 95% CI

widths_and_misses_prop2_RW1 <- calc_width_CI_and_count_obs_outside(proper2_RW1_Extremadura$marginals.fitted.values,
                                                                    Data_LungCancer$obs,
                                                                    Data_LungCancer$pop,
                                                                    n)  


widths_and_misses_prop2_impEff <- calc_width_CI_and_count_obs_outside(proper2_impEff_Extremadura$marginals.fitted.values,
                                                                   Data_LungCancer$obs,
                                                                   Data_LungCancer$pop,
                                                                   n)

widths_and_misses_impIV <- calc_width_CI_and_count_obs_outside(Improper1_typeIV_Extremadura$marginals.fitted.values,
                                                                      Data_LungCancer$obs,
                                                                      Data_LungCancer$pop,
                                                                      n)










