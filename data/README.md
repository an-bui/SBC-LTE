## Data

All files were pulled from [SBC LTER Data Catalog](https://sbclter.msi.ucsb.edu/data/catalog/).

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

| Directory | Files | Notes | 
| -------- | ----- | ----------- | 
| [algae](https://doi.org/10.6073/pasta/47db4ee01f516b0a47b7c585fd552645) | LTE_Algae_Biomass_at_section_20220422.csv, LTE_Algae_Biomass_at_transect_20220314.csv | combined universal point count (UPC) and density converted to biomass |
| [all-species-biomass](https://doi.org/10.6073/pasta/47db4ee01f516b0a47b7c585fd552645) | LTE_All_Species_Biomass_at_transect_20210209.csv, LTE_All_Species_Biomass_at_transect_20220202.csv, LTE_All_Species_Biomass_at_transect_20220208.csv, LTE_All_Species_Biomass_at_transect_20220314.csv | Biomass taken from UPC (cover) or density (quad swath), converted using coefficients for estimating biomass from body size or percent cover for [macroalgae, invertebrates, and fish](https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-sbc&identifier=127) |
| [allometrics](https://doi.org/10.6073/pasta/edff194a827cf11ce80e2ce07d14bf2f) | Allometrics_All_Years_20131009.csv | [morphometric measurements](https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Long_Term_Kelp_Removal/Long_Term_Experiment_Protocol_Understory_Kelp_Allometrics.pdf) of _Pterygophora californica_ and _Laminaria farlowii_


| [Cover of sessile organisms, UPC](https://doi.org/10.6073/pasta/9ef0a3d317f6553e1600a0e5af016e43) | LTE_Cover_All_Years_20200605.csv, LTE_Substrate_All_Years_20200605.csv | 80 points along transect, species percent cover is determined as the fraction of points a species intercepts x 100, kelps (_Macrocystis_, _Pterygophora_, _Eisenia_, _Laminaria_) only measured using holdfasts, includes all sessile organisms encountered | 
| [Abundance and size of Giant Kelp](https://doi.org/10.6073/pasta/5bf131bc3b03ec9f59dc885629065824) | LTE_Kelp_All_Years_20200605.csv | Density measured along whole transect |
| [Fish abundance](https://doi.org/10.6073/pasta/ecf2e269db7a4807bcaa765422d8186c) | LTE_All_Fish_All_Years_20200605.csv | Density measured along whole transect |
| [Invertebrate and algal density](https://doi.org/10.6073/pasta/731d8515e67243716ccb4ee7a28b8843) | LTE_Quad_Swath_All_Years_20200605.csv | Density measured in 6 quadrats along transect and/or in 20 x 1m swaths for select species (list in methods [here](https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Long_Term_Kelp_Removal/Long%20Term%20Experiment%20Protocol%20-%20Density%20of%20algae%20and%20invertebrates_5-30-20.pdf)) |
| [Biomass of kelp forest species](https://doi.org/10.6073/pasta/47db4ee01f516b0a47b7c585fd552645) | LTE_All_Species_Biomass_at_transect_20200605.csv | Biomass taken from UPC (cover) or density (quad swath), converted using coefficients for estimating biomass from body size or percent cover for [macroalgae, invertebrates, and fish](https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-sbc&identifier=127) |
| [Detritus Biomass](https://doi.org/10.6073/pasta/25ae07a87d5764c8eca62b88d695dd50) | LTE_Detritus_Biomass_All_Years_20200605.csv | Measured in 6 quadrats along transect, detritus is collected and weighed in lab |
| [Urchin size frequency distribution](https://doi.org/10.6073/pasta/5a1e9ef03aa47bd2225c0bb98a02a63b) | LTE_Urchin_All_Years_20200605.csv | At least 50 individuals of _Strongylocentrotus franciscanus_ and _S. purpuratus_ measured along transects |
| [Understory kelp allometrics](https://doi.org/10.6073/pasta/53f4dea6d9cc028760859d386be6169c) | Allometrics_All_Years_20131009.csv | Taken for _Laminaria_ and _Pterygophora_ between 2008-2012 to predict biomass from density data of different size classes |
| [Transect depth data](https://doi.org/10.6073/pasta/5b9116a15e1b2b47177ac835b6652596) | LTE_Transect_Depths_coors.csv | taken once by dive computer |
| [Taxon-specific seasonal net primary production (NPP) for macroalgae](https://doi.org/10.6073/pasta/d338c48ec580c052a59aec02c847c2bc) | Understory_NPP_All_Year_season_20200605.csv, Understory_Summer_Bmass...20191122.csv | [methods -_-](https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Long_Term_Kelp_Removal/SBC_LTER_protocol_Reed_LTE_NPP_macroalgae_20200821.pdf) |
| [Hourly photon irradiance at the surface and seafloor](https://doi.org/10.6073/pasta/803abbcd7fb33bbfa9eff08521a397e8) | Hourly_Irrandiance_All_Year_20200320.csv | PAR sensors |
