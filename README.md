# Master

## Simulation study

### Preliminaries and simulating data
- In "Preliminaries.R", the necessary code to generate the B-spline basis functions and tensor product smooth penalization matrices is contained.
- In "Simulate_data_first_level_administrative_regions.R", the code for simulating data over the ADM1 map is contained. Running the "Preliminaries.R" is necessary to do beforehand.
- In "Simulate_data_second_level_administrative_regions.R", the code for simulating data over the ADM4 map is contained.
The function "simulate_risk_surface" is within the "Utilities.R" file.
In the "Validation_risk_surfaces.R" you can generate GIFs of a simulated risk-surface, for choosen scenarios. GIFs of the simulated data for the thesis are "sc1": ADM1 const, short; "sc3": ADM1 lin, short; "sc5": ADM1 cp, short; "sc7": ADM1 const, long; "sc9": ADM1 lin, long; "sc11": ADM1 cp, long; "sc2": ADM4 const, short; "sc4": ADM4 lin, short; "sc6": ADM4 cp, short; "sc8": ADM4 const, long; "sc10": ADM4 lin, long; "sc12": ADM4 cp, long. 

### Analysis
- The files containing the code to analyze the simulation data are in the folder "analyze_simulation_data".
- 

