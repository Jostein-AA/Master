
source("Utilities.R")

### These are here to see names of models and scenarios
#---
#model_names = c("Improper1_noInt", "Improper1_typeI", "Improper1_typeII",
#                "Improper1_typeIII", "Improper1_typeIV",
#                "Improper2_noInt", "Improper2_typeI", "Improper2_typeII",
#                "Improper2_typeIII", "Improper2_typeIV",
#                "proper1_noInt", "proper1_onlyInt", "proper1_full",
#                "proper2_noInt", "proper2_onlyInt", "proper2_full")

#scenario_names = c("sc1", "sc2", "sc3", "sc4",
#                   "sc5", "sc6", "sc7", "sc8",
#                   "sc9", "sc10", "sc11", "sc12")
#---


## Specify model_name and scenario name, and function does the rest

initialize_csv_tracker(model_name = "proper1_noInt",
                       scenario = "sc1")

