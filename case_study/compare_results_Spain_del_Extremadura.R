#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source("libraries.R")
source("Utilities.R")

library(bigDM)

load("case_study/proper2_RW1_Spain_del_Extremadura.RData")
load("case_study/proper2_impEff_Spain_del_Extremadura.RData")
load("case_study/imp_typeIV_Spain_del_Extremadura.RData")

################################################################################
plot(proper2_RW1_Spain_del_Extremadura)

plot(proper2_impEff_Spain_del_Extremadura)

plot(imp_typeIV)

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


tmp <- calc_model_choice(proper2_RW1_Spain_del_Extremadura,
                         Data_LungCancer,
                         n, tT, 
                         Improper = F)


tmp2 <- calc_model_choice(imp_typeIV,
                          Data_LungCancer,
                          n, tT, 
                          Improper = T)

tmp3 <- calc_model_choice(proper2_impEff_Spain_del_Extremadura,
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
cat(latex_tabular, file = "./case_study/Spain_del_Extremadura_model_choice.tex")

################################################################################

# Plot things...
