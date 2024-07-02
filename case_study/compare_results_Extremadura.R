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
### Load in data and Map
# Load in the considered lung-cancer data
data(Data_LungCancer) # Areas in numerically increasing order in terms of ID for each year?

# Load in a map of Spain containing 7907 areas
data("Carto_SpainMUN") # Areas in numerically increasing order in terms of ID?

### Remove disjointed area Llivia
# Find the ID (and other aspects of disjointed area)
problem_area = na.omit(Carto_SpainMUN[Carto_SpainMUN$name == "Llivia", ]) 

# Remove it from the map! (it has index 2454 in Carto)
Carto_SpainMUN <- Carto_SpainMUN[-as.integer(rownames(problem_area)), ]

# Remove Llivia from Data_lungCancer
Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID != problem_area$ID, ]

### Extract the data within the principality of Extremadura
# Extract for the map
Carto_SpainMUN <- Carto_SpainMUN[Carto_SpainMUN$region == "Extremadura", ]

# Find the IDs of the areas within Extremadura
IDs_Extremadura <- unique(Carto_SpainMUN$ID) # IDs either 06... or 10...

# Extract the lung cancer data within Extremadura
Data_LungCancer <- Data_LungCancer[Data_LungCancer$ID %in% IDs_Extremadura, ]



### Get the number of years, start-year and end-year.
# Get the unique years
years = unique(Data_LungCancer$year)

# Get the number of years
tT = length(years)

# Get start-year and end-year
t.from <- min(Data_LungCancer$year)
t.to <- max(Data_LungCancer$year)

### Create IDs for INLA-effects
# Create time-ids 1,...,25 instead of 1991,...,2015 for the sake of INLA
Data_LungCancer$year_id <- Data_LungCancer$year - min(Data_LungCancer$year) + 1


### add area IDs and find number of areas
# find number of areas
n = nrow(Carto_SpainMUN)

# Add unique area_ids to each area in Carto_SpainMUN (that is more easily read for me i.e. 1,...,380)
Carto_SpainMUN$area_id = 1:n

# Add area_id to Data_LungCancer
Data_LungCancer$area_id = rep(NA, nrow(Data_LungCancer))
for(year_id in 1:tT){
  # Insert area_id for Data_LungCancer that corresponds to Carto_SpainMUN
  for(i in 1:nrow(Carto_SpainMUN)){
    Data_LungCancer[Data_LungCancer$year_id == year_id & 
                      Data_LungCancer$ID == Carto_SpainMUN[i, ]$ID, ]$area_id = Carto_SpainMUN[i, ]$area_id
  }
}


### Reset rownames and add space.time id
row.names(Carto_SpainMUN) = 1:nrow(Carto_SpainMUN)
row.names(Data_LungCancer) = 1:nrow(Data_LungCancer)
Data_LungCancer$space.time = 1:nrow(Data_LungCancer)


### Plot the 0.25, 0.5, adn 0.75 quantiles of the population for each year
pop.df <- data.frame(year = t.from:t.to,
                     u = rep(NA, tT),
                     median = rep(NA, tT),
                     l = rep(NA, tT))

for(year in t.from:t.to){
  
  pop_quants_in_year <- quantile(Data_LungCancer[Data_LungCancer$year == year, ]$pop,
                                 probs = c(0.25, 0.5, 0.75))
  
  pop.df[pop.df$year == year, ]$l = pop_quants_in_year[[1]]
  pop.df[pop.df$year == year, ]$median = pop_quants_in_year[[2]]
  pop.df[pop.df$year == year, ]$u = pop_quants_in_year[[3]]
}

# Save as Extremadura_population 5 by 3
ggplot(data = pop.df) + 
  geom_line(aes(x = year, y = median)) + 
  geom_ribbon(aes(x = year, ymin = l, ymax = u), alpha = 0.3) + 
  theme_bw() + xlab("Year") + ylab("Population") + ylim(0, 1250) + 
  theme(axis.title = element_text(size = 14, face = "plain"),
        axis.text = element_text(size = 14)) + 
  scale_x_continuous(breaks = c(1991, 1995, 2000, 2005, 2010, 2015))
  




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
                                is_4 = NA, is_5 = NA, is_avg = NA)
  
  
  # If proper models, sort the marginals to get the correct order
  if(!Improper){
    print("Sorting")
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
    
    print(length(unique(data$year[((year - 1) * n + 1):(year * n)])))
    print(data$year[(year - 1) * n + 1])
    
    
    ## MSE
    model_choice.df[1, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset(one_year_counts,
                                                                                      one_year_margs,
                                                                                      E_it = one_year_pop)
    
    
    ## IS for count in year ahead 
    model_choice.df[1, year - years_pred_on[1] + 7] =  count_IS_one_year_case_study(counts = one_year_counts, 
                                                                                    marginals = one_year_margs, 
                                                                                    population = one_year_pop)
    
    
    
    
    
    
    # Check the width of the 95% CI of the predicted rate. If that is increasing over the years
    # Then, the 95% CI for the predicted number of counts ought to increase as well!!! 
    
    
    
    
    
    
    print(width_CI_one_year_case_study(marginals = one_year_margs, population = one_year_pop))
    
    #print(model_choice.df[1, year - years_pred_on[1] + 7])
    
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

## Calculate model choice critera: MSE, IS, log-score
calc_model_choice_scaled <- function(model,
                              data,
                              n, tT, poppy = 1E4,
                              Improper = T){
  
  ### Initialize data frame to store the MSE's, IS's, and log-scores 1, 2, 3, 4, and 5 years ahead and average
  model_choice.df <- data.frame(mse_1 = NA, mse_2 = NA, mse_3 = NA,
                                mse_4 = NA, mse_5 = NA, mse_avg = NA,
                                is_1 = NA, is_2 = NA, is_3 = NA,
                                is_4 = NA, is_5 = NA, is_avg = NA)
  
  
  # If proper models, sort the marginals to get the correct order
  if(!Improper){
    print("Sorting")
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
    
    print(length(unique(data$year[((year - 1) * n + 1):(year * n)])))
    print(data$year[(year - 1) * n + 1])
    
    
    ## MSE
    model_choice.df[1, year - years_pred_on[1] + 1] =  count_mse_one_year_one_dataset_scaled(one_year_counts,
                                                                                      one_year_margs,
                                                                                      E_it = one_year_pop,
                                                                                      poppy,
                                                                                      version_1 = T)
    
    
    ## IS for count in year ahead 
    model_choice.df[1, year - years_pred_on[1] + 7] =  count_IS_one_year_case_study_scaled(counts = one_year_counts, 
                                                                                    marginals = one_year_margs, 
                                                                                    population = one_year_pop,
                                                                                    poppy,
                                                                                    version_1 = T)
    
    
  }
  
  #Get the total MSE for count
  model_choice.df[1, 6] = mean(as.numeric(model_choice.df[1, 1:5]))
  
  
  #Get the total IS for count
  model_choice.df[1, 12] = mean(as.numeric(model_choice.df[1, 7:11]))
  
  return(model_choice.df)
}



tmp <- calc_model_choice_scaled(proper2_RW1_Extremadura,
                         Data_LungCancer,
                         n, tT, poppy = 1E4,
                         Improper = F)


tmp2 <- calc_model_choice_scaled(Improper1_typeIV_Extremadura,
                          Data_LungCancer,
                          n, tT, poppy = 1E4,
                          Improper = T)

tmp3 <- calc_model_choice_scaled(proper2_impEff_Extremadura,
                          Data_LungCancer,
                          n, tT, poppy = 1E4,
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
cat(latex_tabular, file = "./case_study/Extremadura_model_choice_scaled.tex")

################################################################################

# Just for the sake of simplicity, sort the marginals of the proper models now
if(FALSE){ # if(FALSE) added so that dont sort them unless I really want to
  proper2_RW1_Extremadura$marginals.fitted.values <- sort_proper_fitted(proper2_RW1_Extremadura$marginals.fitted.values,
                                                                        n, tT)
  
  proper2_RW1_Extremadura$summary.fitted.values$mean <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$mean,
                                                                           n, tT)
  
  proper2_RW1_Extremadura$summary.fitted.values$'0.5quant' <- sort_proper_fitted(proper2_RW1_Extremadura$summary.fitted.values$'0.5quant',
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
pred_counts <- matrix(c(proper2_RW1_Extremadura$summary.fitted.values$'0.5quant' * 1E5,
                       Improper1_typeIV_Extremadura$summary.fitted.values$'0.5quant' * 1E5), nrow = n, ncol = 2 * tT, byrow = F)

colnames(pred_counts) = c(paste("proper2_RW1: Year", seq(t.from, t.to), sep = "."),
                          paste("Improper1_typeIV: Year", seq(t.from, t.to), sep = "."))   #paste("Year", seq(t.from, t.to), sep = ".")


# Predicted number of counts per 100,000
pred_count_impIV <- matrix(Improper1_typeIV_Extremadura$summary.fitted.values$'0.5quant' * 1E5, 
                           nrow = n, ncol = tT, byrow = F)
colnames(pred_count_impIV) = paste("Year", seq(t.from, t.to), sep = ".")


# Create a common color-palet
paleta <- brewer.pal(8,"RdYlGn")[8:1]
pred_count_qs <- quantile(pred_counts, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
values <- c(min(pred_counts), pred_count_qs, Inf)



# Plot the predicted counts per 100,000
carto <- cbind(Carto_SpainMUN, pred_counts)

# paleta <- brewer.pal(8,"RdYlGn")[8:1]
# pred_count_qs <- quantile(pred_count_prop2_RW1, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
# values <- c(min(pred_count_prop2_RW1), pred_count_qs, Inf)

columns_ <- c(paste("proper2_RW1..Year", c(2011, 2013, 2015),sep= "."),
              paste("Improper1_typeIV..Year", c(2011, 2013, 2015),sep= "."))

Map.risks <- tm_shape(carto) +
  tm_polygons(col = columns_,
              palette=paleta, 
              title="Posterior median\npredicted rate\nper 100,000", 
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
            panel.label.size=1.55,
            legend.outside=T, 
            legend.outside.position="right", 
            legend.frame=F,
            legend.outside.size = 0.2, #0.15
            legend.title.size = 1.5,
            legend.text.size = 1.25,
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0.001, 0.001, 0.001, 0.001),
            between.margin = 0.001,
            panel.labels=as.character(c(paste("proper2_RW1:", c(2011, 2013, 2015), sep = " "),
                                        paste("Improper1_typeIV:", c(2011, 2013, 2015), sep = " ")))
            ) +
  tm_facets(nrow=2, ncol=3)


print(Map.risks)

# Save Map
tmap_save(Map.risks,
          filename = "Plots/Extremadura_pred_count.pdf",
          width = 12,
          height = 6)

### Difference plot
diff <- (pred_count_prop2_RW1 - pred_count_impIV)/(pred_count_impIV) 
quantiles.diff <- quantile(diff, probs = c(0.025, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.975))

carto_diff <- cbind(map_Spain, diff)

paleta.diff <- brewer.pal(9,"RdBu")[9:1]
values.diff <- c(-Inf, quantiles.diff, Inf)
  
  
Map.diff <- tm_shape(carto_diff) +
  tm_polygons(col=paste("Year", c(2011, 2013, 2015),sep= "."),
              palette=paleta.diff, 
              title="Relative difference\nin posterior mean\npredicted counts", #(Improper1_typeIV\n - proper2_RW1)/\nproper2_RW1
              legend.show=T, 
              border.col="transparent",
              legend.reverse=T, 
              style="fixed", 
              breaks = values.diff,
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


print(Map.diff)

tmap_save(Map.diff,
          filename = "Plots/Extremadura_pred_counts_relative_diff.pdf",
          width = 12, #12 by 4
          height = 4)

########
# SD count

# Predicted SD of counts per 100,000
sd_2011 <-  get_pred_SD(proper2_RW1_Extremadura$marginals.fitted.values, 1E5, n, 21)
sd_2013 <- get_pred_SD(proper2_RW1_Extremadura$marginals.fitted.values, 1E5, n, 23)
sd_2015 <- get_pred_SD(proper2_RW1_Extremadura$marginals.fitted.values, 1E5, n, 25)
# Predicted SD of counts per 100,000
sd_IV_2011 <- get_pred_SD(Improper1_typeIV_Extremadura$marginals.fitted.values, 1E5, n, 21)
sd_IV_2013 <- get_pred_SD(Improper1_typeIV_Extremadura$marginals.fitted.values, 1E5, n, 23)
sd_IV_2015 <- get_pred_SD(Improper1_typeIV_Extremadura$marginals.fitted.values, 1E5, n, 25)


pred_SD <- matrix(c(sd_2011, sd_2013, sd_2015,
                    sd_IV_2011, sd_IV_2013, sd_IV_2015), nrow = n, ncol = 2 * 3, byrow = F)

colnames(pred_SD) = c(paste("proper2_RW1: Year", c(2011, 2013, 2015), sep = "."),
                                paste("Improper1_typeIV: Year", c(2011, 2013, 2015), sep = "."))



# Create a common color-palet
paleta <- brewer.pal(8,"RdYlGn")[8:1]
pred_count_qs <- quantile(pred_SD, probs = c(0.1, 0.25, 0.45, 0.6, 0.75, 0.875, 0.975))
values <- round(c(min(pred_SD) - 0.1, pred_count_qs, Inf), 1)



# Plot the predicted SD of counts per 100,000 for the proper2_RW1
carto_SD <- cbind(Carto_SpainMUN, pred_SD)

# paleta <- brewer.pal(8,"RdYlGn")[8:1]
# pred_count_qs <- quantile(pred_count_prop2_RW1, probs = c(0.15, 0.3, 0.45, 0.6, 0.75, 0.875, 0.975))
# values <- c(min(pred_count_prop2_RW1), pred_count_qs, Inf)


Map.risks_SD <- tm_shape(carto_SD) +
  tm_polygons(col = columns_,
              palette=paleta, 
              title="Posterior SD of\npredicted rate\nper 100,000", 
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
            legend.outside.size = 0.2, #0.15
            legend.title.size = 1.5,
            legend.text.size = 1.25,
            outer.margins=c(0.01,0.01,0.02,0.01),
            inner.margins = c(0.001, 0.001, 0.001, 0.001),
            between.margin = 0.001,
            panel.labels=as.character(c(paste("proper2_RW1:", c(2011, 2013, 2015), sep = " "),
                                        paste("Improper1_typeIV:", c(2011, 2013, 2015), sep = " ")))
            ) +
  tm_facets(nrow=2, ncol=3)


print(Map.risks_SD)

# Save Map
tmap_save(Map.risks_SD,
          filename = "Plots/Extremadura_SD.pdf",
          width = 12,
          height = 6)





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


################################################################################

Extremadura_grid <- data.frame(
  code = c(15,  203, 356),
  name = c("Badajoz", "Caminomorisco", "Trujillo"),
  row = c(1, 1, 1),
  col = c(1, 2, 3),
  stringsAsFactors = FALSE)


geofacet::grid_preview(Extremadura_grid)


pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, #lambda_$area_id
                           time_id = Data_LungCancer$year_id, # lambda_$time_id
                           median = Improper1_typeIV_Extremadura$summary.fitted.values$'0.5quant' * 1E5, #sort_proper_fitted(model$summary.fitted.values$'0.5quant', length(unique(lambda_$area_id)), tT) * lambda_$E_it, # model$summary.fitted.values$'0.5quant',
                           quantile_0.025 = Improper1_typeIV_Extremadura$summary.fitted.values$'0.025quant' * 1E5,# model$summary.fitted.values$'0.025quant',
                           quantile_0.975 = Improper1_typeIV_Extremadura$summary.fitted.values$'0.975quant' * 1E5,
                           y_it_div_E_it = Data_LungCancer$obs/Data_LungCancer$pop)


# Save as Extremadura_timeseries 10 by 3.5
ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI for the rate per 100,000"), 
              fill = "#F8766D", alpha = 0.6) +
  #geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  #geom_point((aes(x = time_id, y = y_it_div_E_it, col = "y/E"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median rate per 100,000")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(title = NULL,
       x = "Year",
       y = "Rate per 100,000",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=13),
        plot.title = element_text(hjust = 0.5, size=13),
        strip.text.x = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.position = "top",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))




### Posterior predicted counts



pred_to_plot <- data.frame(area_id = Data_LungCancer$area_id, 
                           time_id = Data_LungCancer$year_id, 
                           median = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           #post_mean = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.025 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           quantile_0.975 = rep(NA, nrow(Improper1_typeIV_Extremadura$summary.fitted.values)),
                           y_it = Data_LungCancer$obs)


for(area_id in Extremadura_grid$code){
  for(t in 1:tT){
    ul <- find_ul_quants_counts_single_pred(Improper1_typeIV_Extremadura$marginals.fitted.values[[(t - 1) * n + area_id]],
                                            Data_LungCancer[Data_LungCancer$area_id == area_id & Data_LungCancer$year_id == t, ]$pop)
    
    
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$median = ul$median
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.975 = ul$u
    pred_to_plot[pred_to_plot$area_id == area_id &
                   pred_to_plot$time_id == t, ]$quantile_0.025 = ul$l
  }
}



# Save as Extremadura_timeseries_counts 10 by 3.5
ggplot(data = pred_to_plot, aes(time_id, median)) + 
  geom_ribbon(aes(x = time_id, ymin = quantile_0.025, ymax = quantile_0.975, col = " Posterior 95% CI for the count"), 
              fill = "#F8766D", alpha = 0.6) +
  # geom_point(aes(x = time_id, y = sampled_counts, col = "True count")) + 
  geom_point((aes(x = time_id, y = y_it, col = "Observed y"))) +
  geom_line(aes(x = time_id, y = median, col = "Posterior median count")) +
  geom_vline(xintercept = 20.5, linetype = "longdash", 
             color = "darkgrey", linewidth = 0.6) +
  facet_geo(~ area_id, grid = Extremadura_grid, label = "name") + 
  labs(x = "Year",
       y = "y",
       col = NULL) +
  theme_bw() + 
  theme(axis.title=element_text(size=13),
        plot.title = element_text(hjust = 0.5, size=13),
        strip.text.x = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.position = "top",
        legend.background = element_rect(linetype = 1, linewidth = 1, colour = "grey"))













