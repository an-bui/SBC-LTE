
# 1. libraries ------------------------------------------------------------

library(tidyverse)
library(here)
library(janitor)
library(lubridate)

# 2. data -----------------------------------------------------------------

############################################
# a. biomass
############################################

biomass <- read_csv(here::here("data", "LTE_All_Species_Biomass_at_transect_20200605.csv")) %>% 
  clean_names() %>% 
  # ANOB is incorrectly coded as having "SESSILE" mobility
  mutate(mobility = ifelse(sp_code == "ANOB", "MOBILE", mobility)) %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # create "functional groups" of group and mobility
  unite("group_mobility", group, mobility, remove = FALSE, sep = "_") %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "group_mobility", "growth_morph"), str_to_lower)


# 3. operators ------------------------------------------------------------

# not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))


# 4. start and end dates --------------------------------------------------

############################################
# a. Naples (NAPL)
############################################

napl_start_dates <- c("NAPL_CONTROL_2008-01-10", "NAPL_ANNUAL_2008-01-10", "NAPL_CONTINUAL_2010-04-27")

napl_start_date <- as_date("2008-01-10")

napl_after_date <- as_date("2016-02-09")

############################################
# b. Mohawk (MOHK)
############################################

mohk_start_dates <- c("MOHK_ANNUAL_2008-01-17", "MOHK_CONTROL_2008-01-17", "MOHK_CONTINUAL_2010-05-05")

mohk_start_date <- as_date("2008-01-17")

mohk_after_date <- as_date("2017-02-13")

############################################
# c. Arroyo Quemado (AQUE)
############################################

aque_start_dates <- c("AQUE_CONTROL_2008-01-30", "AQUE_ANNUAL_2008-01-30", "AQUE_CONTINUAL_2010-04-26")

aque_start_date <- as_date("2008-01-30")

aque_after_date <- as_date("2017-03-02")

############################################
# d. Carpinteria (CARP)
############################################

carp_start_dates <- c("CARP_CONTROL_2008-02-12", "CARP_ANNUAL_2008-02-12", "CARP_CONTINUAL_2010-04-23")

carp_start_date <- as_date("2008-02-12")

carp_after_date <- as_date("2017-02-15")

############################################
# e. Isla Vista (IVEE)
############################################

ivee_start_dates <- c("IVEE_CONTROL_2011-10-26", "IVEE_ANNUAL_2011-10-26")

ivee_start_date <- as_date("2011-10-26")

ivee_after_date <- as_date("2016-02-18")


# 5. useful vectors and data frames ---------------------------------------

############################################
# a. species
############################################

# data frame of each species including scientific name, common name, and group
spp_names <- biomass %>% 
  select(group, group_mobility, sp_code, scientific_name, common_name) %>% 
  unique() 

# vector of each species code (length: 221)
spp_codes <- biomass %>% 
  pull(sp_code) %>% 
  unique()

# vector of algae species
algae_spp <- spp_names %>% 
  filter(group == "ALGAE") %>% 
  pull(sp_code)

# vector of fish species
fish_spp <- spp_names %>% 
  filter(group == "FISH") %>% 
  pull(sp_code)

# vector of invert species
invert_spp <- spp_names %>% 
  filter(group == "INVERT") %>% 
  pull(sp_code)

# vector of group_mobility
group_mobility <- biomass %>% 
  pull(group_mobility) %>% 
  unique()

# data frame of mobile invert species
invert_mobile <- biomass %>% 
  filter(group_mobility == "invert_mobile") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

# data frame of sessile invert species
invert_sessile <- biomass %>% 
  filter(group_mobility == "invert_sessile") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

# data frame of mobile fish species
fish_mobile <- biomass %>% 
  filter(group_mobility == "fish_mobile") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

# vector of group_mobility
group_mobility <- biomass %>% 
  pull(group_mobility) %>% 
  unique()

# data frame of mobile invert species
invert_mobile <- biomass %>% 
  filter(group_mobility == "invert_mobile") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

# data frame of sessile invert species
invert_sessile <- biomass %>% 
  filter(group_mobility == "invert_sessile") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

# data frame of mobile fish species
fish_mobile <- biomass %>% 
  filter(group_mobility == "fish_mobile") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

############################################
# b. sites
############################################

# vector of sites
LTE_sites <- biomass %>% 
  pull(site) %>% 
  unique()

