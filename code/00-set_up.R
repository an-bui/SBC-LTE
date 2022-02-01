
# 1. libraries ------------------------------------------------------------

# general cleaning, visualization, etc.
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(patchwork)
library(calecopal)
library(ggrepel)
library(plotly)
library(fuzzyjoin)
library(gt)
library(rlang)

# analysis
library(vegan)
library(vegclust)
library(ecotraj)
library(FD)
library(BiodiversityR)


# 2. start and end dates --------------------------------------------------

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

# 3. data -----------------------------------------------------------------

############################################
# a. biomass
############################################

biomass <- read_csv(here::here("data", "LTE_All_Species_Biomass_at_transect_20210209.csv")) %>% 
  clean_names() %>% 
  # ANOB is incorrectly coded as having "SESSILE" mobility
  mutate(mobility = replace(mobility, sp_code == "ANOB", "MOBILE")) %>% 
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA)) %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # create "functional groups" of group and mobility
  unite("group_mobility", group, mobility, remove = FALSE, sep = "_") %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "group_mobility", "growth_morph"), str_to_lower) %>% 
  # make a new column designating "start", "during" and "after" removal
  mutate(exp_dates = case_when(
    site == "NAPL" & sample_ID %in% napl_start_dates ~ "start",
    site == "NAPL" & date > napl_after_date ~ "after",
    site == "MOHK" & sample_ID %in% mohk_start_dates ~ "start",
    site == "MOHK" & date > mohk_after_date ~ "after",
    site == "AQUE" & sample_ID %in% aque_start_dates ~ "start",
    site == "AQUE" & date > aque_after_date ~ "after",
    site == "CARP" & sample_ID %in% carp_start_dates ~ "start",
    site == "CARP" & date > carp_after_date ~ "after",
    site == "IVEE" & sample_ID %in% ivee_start_dates ~ "start",
    site == "IVEE" & date > ivee_after_date ~ "after",
    TRUE ~ "during"
  ),
  exp_dates = factor(exp_dates, levels = c("start", "during", "after"))) 

############################################
# b. percent cover
############################################

percov <- read_csv(here::here("data", "LTE_Cover_All_Years_20200605.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph"), str_to_lower) %>% 
  # make a new column designating "start", "during" and "after" removal
  mutate(exp_dates = case_when(
    site == "NAPL" & sample_ID %in% napl_start_dates ~ "start",
    site == "NAPL" & date > napl_after_date ~ "after",
    site == "MOHK" & sample_ID %in% mohk_start_dates ~ "start",
    site == "MOHK" & date > mohk_after_date ~ "after",
    site == "AQUE" & sample_ID %in% aque_start_dates ~ "start",
    site == "AQUE" & date > aque_after_date ~ "after",
    site == "CARP" & sample_ID %in% carp_start_dates ~ "start",
    site == "CARP" & date > carp_after_date ~ "after",
    site == "IVEE" & sample_ID %in% ivee_start_dates ~ "start",
    site == "IVEE" & date > ivee_after_date ~ "after",
    TRUE ~ "during"
  ),
  exp_dates = factor(exp_dates, levels = c("start", "during", "after"))) 

############################################
# c. traits
############################################

traits <- read_csv(here::here("data/traits", "00-coarse_traits.csv"))

############################################
# d. LTER annual monitoring
############################################

percov_annual <- read_csv(here::here("data/benthics", "Annual_Cover_All_Years_20210108.csv"))

biomass_annual <- read_csv(here::here("data/benthics", "Annual_All_Species_Biomass_at_transect_20210108.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph"), str_to_lower)

############################################
# e. LTE kelp biomass
############################################

kelp_fronds <- read_csv(here::here("data", "LTE_Kelp_All_Years_20210209.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph"), str_to_lower) %>% 
  # make a new column designating "start", "during" and "after" removal
  mutate(exp_dates = case_when(
    site == "NAPL" & sample_ID %in% napl_start_dates ~ "start",
    site == "NAPL" & date > napl_after_date ~ "after",
    site == "MOHK" & sample_ID %in% mohk_start_dates ~ "start",
    site == "MOHK" & date > mohk_after_date ~ "after",
    site == "AQUE" & sample_ID %in% aque_start_dates ~ "start",
    site == "AQUE" & date > aque_after_date ~ "after",
    site == "CARP" & sample_ID %in% carp_start_dates ~ "start",
    site == "CARP" & date > carp_after_date ~ "after",
    site == "IVEE" & sample_ID %in% ivee_start_dates ~ "start",
    site == "IVEE" & date > ivee_after_date ~ "after",
    TRUE ~ "during"
  ),
  exp_dates = factor(exp_dates, levels = c("start", "during", "after"))) %>% 
  select(1:14) %>% 
  filter(fronds > -1) %>% 
  mutate(season = case_when(
    month %in% c(12, 1, 2) ~ "winter",
    month %in% c(3, 4, 5) ~ "spring",
    month %in% c (6, 7, 8) ~ "summer",
    month %in% c(9, 10, 11) ~ "fall"
  ),
  season = fct_relevel(season, "winter", "spring", "summer", "fall"))

# 4. operators ------------------------------------------------------------

# not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))


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
  filter(group == "algae") %>% 
  pull(sp_code)

# vector of fish species
fish_spp <- spp_names %>% 
  filter(group == "fish") %>% 
  pull(sp_code)

# vector of invert species
invert_spp <- spp_names %>% 
  filter(group == "invert") %>% 
  pull(sp_code)

# vector of groups
all_groups <- biomass %>% 
  pull(group) %>% 
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

# data frame of algae species
algae_spp <- biomass %>% 
  filter(group == "algae") %>% 
  select(sp_code, scientific_name, common_name) %>% 
  unique()

# most abundant algae
common_algae <- c("PH", "PTCA", # Pterygophora californica 
                  "DL", # Desmarestia ligulata
                  "R", # Rhodymenia californica 
                  "CC", # Chondracanthus corymbiferus 
                  "POLA", # Polyneura latissima 
                  "CYOS", # Stephanocystis osmundacea 
                  "FTHR", # Pterosiphonia dendroidea 
                  "CO", # Corallina officinalis var. chilensis 
                  "LX", # Osmundea spectabilis
                  "GS", # Gracilaria spp. 
                  "BR", # Halymenia spp.
                  "BO", # Bossiella orbigniana 
                  "FB", # Ectocarpaceae spp. 
                  "BF", # Cryptopleura ruprechtiana 
                  "LAFA", # Laminaria farlowii 
                  "CF", # Callophyllis rhynchocarpa 
                  "DP" # Dictyota spp. 
)

common_algae_df <- algae_spp %>% 
  filter(sp_code %in% common_algae)
# not found: 
# Pseudolithophyllum neofarlowii 
# Polysiphonia spp. 

############################################
# b. sites
############################################

# vector of sites
LTE_sites <- biomass %>% 
  pull(site) %>% 
  unique()

LTER_sites <- biomass_annual %>% 
  select(site) %>% 
  mutate(site = fct_relevel(site, c("BULL", "AQUE", "AHND", "NAPL", "IVEE", "GOLB", "ABUR", "MOHK", "CARP", "SCDI", "SCTW"))) %>% 
  unique() %>% 
  pull(site)

site_quality <- tribble(
  ~site, ~quality,
  "MOHK", "high",
  "IVEE", "high",
  "AQUE", "medium",
  "NAPL", "medium",
  "CARP", "low"
) %>% 
  mutate(quality = fct_relevel(quality, c("low", "medium", "high")))

############################################
# c. sample IDs
############################################

# NAPL after
napl_after_sampleIDs <- percov %>% 
  filter(site == "NAPL" & exp_dates == "after") %>% 
  pull(sample_ID)

# MOHK after
mohk_after_sampleIDs <- percov %>% 
  filter(site == "MOHK" & exp_dates == "after") %>% 
  pull(sample_ID)

# all after
all_after_sampleIDs <- percov %>% 
  filter(exp_dates == "after") %>% 
  pull(sample_ID)


# 6.  plot aesthetics -----------------------------------------------------

AQUE_col <- "#FDD989"
CARP_col <- "#516238"
IVEE_col <- "#4CA2B0"
MOHK_col <- "#985E5C"
NAPL_col <- "#D98A63"

annual_col <- "#54662C"
continual_col <- "#009BB0"
control_col <- "#000000"
kelp_col <- "#6D5A18"
under_col <- "#CC7540"

annual_shape <- 19
continual_shape <- 17
control_shape <- 21




