# Foundation species recovery yields inconsistent recovery of associated community: a long-term experiment

**last updated: 2026-10-22**  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14165444.svg)](https://doi.org/10.5281/zenodo.14165444)

This is a repository for analyzing kelp and kelp-associated community recovery dynamics from the long term kelp removal experiment in the Santa Barbara Coastal Long Term Ecological Research (SBC LTER) network. The analysis accompanies Bui et al. 2026, "Foundation species recovery yields inconsistent recovery of associated community: A long-term experiment" (in preparation for submission).  

## Data citation

This analysis relies on data from: Reed, D. and R. Miller. 2026. SBC LTER: Reef: Long-term experiment: biomass of kelp forest species ver 16. Environmental Data Initiative. https://doi.org/10.6073/pasta/9c2001b2540f2f064c98d40cc2795362 (Accessed 2026-10-22).

Please see the README in the [`data`](https://github.com/an-bui/SBC-LTE/tree/submission/data) directory for information on downloading the data from EDI.

## Repository structure:

Please see the READMEs in the [`code`](https://github.com/an-bui/SBC-LTE/tree/submission/code) and [`data`](https://github.com/an-bui/SBC-LTE/tree/submission/data) directories for more information.  

Files in the `figures` and `tables` directories are not listed in full; please see those directories for all files included in the manuscript.

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
│   └── resources
│       ├── Thiault
│       │   └── mee312655-sup-0001-appendixs1.r
│       └── castorani
│           └── LTE_guild_data.csv
├── data
│   ├── README.md
│   └── all-species-biomass
│       └── knb-lter-sbc.119.13
│           ├── LTE_All_Species_Biomass_at_transect_20241115.csv
│           ├── knb-lter-sbc.119.13.report.xml
│           ├── knb-lter-sbc.119.13.txt
│           ├── knb-lter-sbc.119.13.xml
│           └── manifest.txt
├── figure
│   ├── icons
│   └── ms-figures 
└── tables             
    └── ms-tables
```

## License

This work is licensed under a
[Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).  

[![CC BY 4.0](https://i.creativecommons.org/l/by/4.0/88x31.png)](http://creativecommons.org/licenses/by/4.0/)

