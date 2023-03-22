# Code

**last updated: 2023-03-22**

    ├── code/                        # all analysis code
    │   └── 00-set_up.R                                                   
    │   └── 01a-kelp_recovery.R                             
    │   └── 01b-kelp_BACIPS.R              
    │   └── 01c-kelp_recovery_control.R     
    │   └── 02a-community_recovery.R                        
    |   └── 02b-community_composition.R 
    |   └── 03a-algae_epi_kelp.R   
    |   └── README.md

## Set up
- `00-set_up.R`: libraries, functions, operators, objects  

## Kelp recovery
- `01a-kelp_recovery.R`: sources `00-set_up.R`, analysis of kelp removal and recovery deltas  
- `01b-kelp_BACIPS.R`: sources `01-kelp_recovery.R`, kelp Before-After-Control-Impact-Paired-Series (BACIPS) analysis  
- `01c-kelp_recovery_control.R`: sources `01-kelp_recovery.R`, examines how giant kelp biomass in control plots influences time to recovery metric  

## Community dynamics
- `02a-community_recovery.R`: sources `01a-kelp_recovery.R`, community biomass analysis  
- `02b-community_composition.R`: sources `02a-community_recovery.R`, community composition analysis  

## Relationship between algae/epilithic invertebrates and kelp
- `03a-algae_epi_kelp.R`: sources `02a-community_recovery.R`, delta algae and delta epi vs delta kelp analysis