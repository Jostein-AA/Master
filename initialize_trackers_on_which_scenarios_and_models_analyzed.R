#Clear environment
rm(list = ls())

#Libraries
source("libraries.R")
source("Utilities.R")

## Need to initialize csv's containing an empty data frame.
## Each data frame corresponds to one model and one scenario
## When the model has analyzed one data set in that scenario
## The data set index is added to the data frame and that is written as the
## new csv! If an error occurs, that is also written into the data frame.

#get_csv_tracker_filename(model_name = "Improper1_noInt", 
#                         scenario = "sc1")

#get_first_not_yet_analyzed(model_name = "Improper1_noInt", 
#                           scenario = "sc1")

model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
                "Improper1_typeIII", "Improper1_typeIV",
                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
                "Improper2_typeIII", "Improper2_typeIV",
                "proper1_noInt", "proper1_onlyInt", "proper1_full", "proper1_iid",
                "proper2_noInt", "proper2_onlyInt", "proper2_full", "proper2_iid")

scenario_names = c("sc1", "sc2", "sc3", "sc4",
                   "sc5", "sc6", "sc7", "sc8",
                   "sc9", "sc10", "sc11", "sc12")

for(model_name in model_names){
  for(scenario_name in scenario_names){
    print(get_csv_tracker_filename(model_name = model_name,
                                   scenario = scenario_name))
    
    initialize_csv_tracker(model_name = model_name,
                           scenario = scenario_name)
  }
}















