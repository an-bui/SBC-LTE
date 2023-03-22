# Kelp recovery

**last updated: 2023-03-22**  

This is a repository for analyzing kelp and kelp-associated community recovery dynamics from the long term kelp removal experiment in the Santa Barbara Coastal Long Term Ecological Research (SBC LTER) network. All data are from the [SBC LTER Data Catalog](https://sbclter.msi.ucsb.edu/data/catalog/).

## Repository structure:

    .
    ├── code/                        # all analysis code
    │   └── 00-set_up.R                                     
    │   └── 01a-kelp_recovery.R                             
    │   └── 01b-kelp_BACIPS.R
    │   └── 01c-kelp_recovery_control.R
    │   └── 02a-community_recovery.R                        
    |   └── 02b-community_composition.R 
    |   └── 03a-algae_epi_kelp.R   
    |   └── README.md
    |
    ├── data/                        # raw data files
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
    ├── figures/                      # figures generated from code 
    |   └── ms-figures                                      
    |    
    ├── tables/ms-tables/             # tables generated from code 
    |  
    ├── .gitignore        
    ├── README.md
    └── SBC-LTE.Rproj

Please visit each directory for detailed descriptions of items.




