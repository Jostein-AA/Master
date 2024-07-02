This readme was written after the deadline of the thesis, and serves only as an overview of the code. 
# Master

## Simulation study

### Preliminaries and simulating data
- In "Preliminaries.R", the necessary code to generate the B-spline basis functions and tensor product smooth penalization matrices is contained.
- In "Simulate_data_first_level_administrative_regions.R", the code for simulating data over the ADM1 map is contained. Running the "Preliminaries.R" is necessary to do beforehand.
- In "Simulate_data_second_level_administrative_regions.R", the code for simulating data over the ADM4 map is contained.
The function "simulate_risk_surface" is within the "Utilities.R" file.
In the "Validation_risk_surfaces.R" you can generate GIFs of a simulated risk-surface, for choosen scenarios. GIFs of the simulated data for the thesis are "sc1": ADM1 const, short; "sc3": ADM1 lin, short; "sc5": ADM1 cp, short; "sc7": ADM1 const, long; "sc9": ADM1 lin, long; "sc11": ADM1 cp, long; "sc2": ADM4 const, short; "sc4": ADM4 lin, short; "sc6": ADM4 cp, short; "sc8": ADM4 const, long; "sc10": ADM4 lin, long; "sc12": ADM4 cp, long. 

### Analysis
- Prerquiste before analyzing is initializing trackers of each model, this is done in "initialize_trackers_on_which_scenarios_and_models_analyzed.R".
- The files containing the code to analyze the simulation data are in the folder "analyze_simulation_data". There are two files corresponding to each model, one for the analysis over the ADM1 map, and the other for the analysis over the ADM4 map. Additionally, there are files wrapping over these files for the analysis of only one model. If you are to run the analysis, I would recommend doing this on a high performance computer as it is very time consuming.
- Model chocie critera are calculated in "calculate_model_choice"
- In "fit_models_on_one_data_set_pr_scenario", each model is fitted to one data set per considered scenario, and the entire INLA object is stored

- The plots are made in the file "ploting".

# Case study

- The code related to the case study is found in the folder "case_study"
### 
- The files containing the code to carry out the analysis are for the Extremadura case in files with the ending "Extremadura.R" and for the analysis over the remainder of Spain in files with the ending "Spain_del_Extremadura.R". The Div_and_conquer model is here called Improper1_typeIV.
- Plots and model choice critera over Extremadura are within "compare_results_Extremadura.R"
- Plots and model choice critera over Spain excluding Extremadura are within "compare_results_Spain_del_Extremadura.R"


# Reparametrized space-time interaction
- The code for the reparameterized space-time interaction are found within the folder "typeV"
###
- The file "prior_wrapper.R" contains the code to create the marginal priors on unspecified data.
- The file "make_typeV_prior_ADM4.R" contains the code calling "prior_wrapper.R" with the specifics of our considered analysis.
- Plots are made in "typeV/inference_typeV_vs_typeIV.R".
- Model choice are calculated using the "calculate_model_choice.R" file in the main folder.


