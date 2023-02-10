# Kelp recovery

This is a repository for analyzing kelp and kelp-associated community recovery dynamics from the long term kelp removal experiment in the Santa Barbara Coastal Long Term Ecological Research (SBC LTER) network. All data are from the [SBC LTER Data Catalog](https://sbclter.msi.ucsb.edu/data/catalog/).

## Repository structure:

    .
    ├── code                                       # all analysis code
    │   └── 00-set_up.R                                     # all set up libraries, functions, objects  
    │   └── 01a-kelp_recovery.R                             # analysis, figures, tables for kelp recovery
    │   └── 01b-kelp_BACIPS                                 # analysis, figures, tables for kelp Before-AFter-Control-Impact-Paired-Series (BACIPS) analysis
    │   └── 02a-community_recovery.R                        # analysis, figures, tables for community biomass analysis
    |   └── 02b-community_composition.R                     # analysis, figures for community composition analysis
    |   └── README.md
    |
    ├── data/                                      # raw data files
    |   └── algae                                          
    |   └── all-species-biomass
    |   └── allometrics
    |   └── benthics
    |   └── detritus
    |   └── fish
    |   └── inverts
    |   └── irradiance
    |   └── kelp-fronds
    |   └── percent-cover
    |   └── quads-swaths
    |   └── substrate
    |   └── traits
    |   └── transect-depths
    |   └── understory
    |   └── urchins
    |   └── README.md  
    |
    ├── figures/                                  # folder containing the code for App #1, a single-file-app 
    |   └── ms-figures                                       # app code
    |    
    ├── tables/ms-tables/                                     # folder containing the code for App #2, a two-file-app 
    ├── .gitignore        
    ├── README.md
    └── SBC-LTE.Rproj






