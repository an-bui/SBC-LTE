# Foundation species recovery yields inconsistent recovery of associated community: a long-term experiment

**last updated: 2024-11-14**  

This is a repository for analyzing kelp and kelp-associated community recovery dynamics from the long term kelp removal experiment in the Santa Barbara Coastal Long Term Ecological Research (SBC LTER) network. The analysis accompanies Bui et al. 2024, "Foundation species recovery yields inconsistent recovery of associated community: a long-term experiment" (in preparation for submission to _Ecology Letters_).

## Data citation

This analysis relies on data from: Reed, D. and R. Miller. 2024. SBC LTER: Reef: Long-term experiment: biomass of kelp forest species, ongoing since 2008 ver 11. Environmental Data Initiative. https://doi.org/10.6073/pasta/1d4a114b80d8d3ceb5adc668d5fbe497 (Accessed 2024-10-07).

## Repository structure:
```
.
├── README.md
├── SBC-LTE.Rproj
├── code
│   ├── 00a-set_up.R
│   ├── 00b-getting_data_from_EDI.R
│   ├── 01a-kelp_recovery.R
│   ├── 01b-kelp_BACIPS.R
│   ├── 02a-community_recovery.R
│   ├── 02b-community_composition.R
│   ├── 03a-algae_epi_kelp.R
│   ├── README.md
│   ├── old-code
│   └── resources
│       ├── Thiault
│       │   └── mee312655-sup-0001-appendixs1.r
│       └── castorani
│           ├── LTE_guild_data.csv
│           ├── lte_analysis_2022-02-15.R
│           └── sbc_canopy_project-5_sites_rev3_2021-07-01.r
├── data
│   ├── README.md
│   └── all-species-biomass
│       ├── LTE_All_Species_Biomass_at_transect_20240501.csv
│       └── biomass.RDS
├── figures
│   ├── fg-recovery
│   ├── ms-figures
│   └── old
├── tables
│   ├── ms-tables
│   └── old
```

