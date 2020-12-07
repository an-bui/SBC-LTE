# libraries
library(tidyverse)
library(here)
library(janitor)

# data
biomass <- read_csv(here::here("data", "LTE_All_Species_Biomass_at_transect_20200605.csv")) %>% 
  clean_names() %>% 
  mutate(scientific_name = str_replace(scientific_name, " ", "_"), 
         treatment = str_replace(treatment, " ", "_")) %>% 
  unite("sample_ID", site, treatment, year, remove = FALSE) 
