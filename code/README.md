# Code

**last updated: 2024-10-07**

## Set up
- `00a-set_up.R`: libraries, functions, operators, objects  
- `00b-getting_data_from_EDI.R`: a script to download a zipped version of the data used in this analysis from the Environmental Data Initiative. **Note:** this code only has to be run _once_ to use the code; if the code is already downloaded in the correct place, then start with `00a-set_up.R`.  

## Kelp recovery
- `01a-kelp_recovery.R`: sources `00-set_up.R`, analysis of kelp removal and recovery deltas  
- `01b-kelp_BACIPS.R`: sources `01-kelp_recovery.R`, kelp Before-After-Control-Impact-Paired-Series (BACIPS) analysis  

## Community dynamics
- `02a-community_recovery.R`: sources `01a-kelp_recovery.R`, community biomass analysis  
- `02b-community_composition.R`: sources `02a-community_recovery.R`, community composition analysis  

## Relationship between algae/epilithic invertebrates and kelp
- `03a-algae_epi_kelp.R`: sources `02a-community_recovery.R`, delta algae and delta epi vs delta kelp analysis

## Additional resources
- `castorani`: code from Castorani et al. 2021  
- `Thiault`: script for BACIPS analysis (reproduced in `01b-kelp_BACIPS.R`)