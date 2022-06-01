
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
library(multcompView)
library(ggeffects)
library(ggnewscale)

# analysis
library(vegan)
library(vegclust)
library(ecotraj)
library(FD)
library(BiodiversityR)
library(minpack.lm) # Fitting non-linear models
library(nls2) # Fitting non-linear models
library(AICcmodavg) # calculate second order AIC (AICc)
library(MuMIn)
library(boot)
library(lme4)
library(nlme)
library(DHARMa)
library(performance)
library(emmeans)
library(gtsummary)


# 2. start and end dates --------------------------------------------------

# ⟞ a. Arroyo Quemado (AQUE) ---------------------------------------------

aque_start_dates <- c("AQUE_CONTROL_2008-01-30", 
                      "AQUE_ANNUAL_2008-01-30", 
                      "AQUE_CONTINUAL_2010-04-26")

aque_start_date <- as_date("2008-01-30")

aque_after_date <- as_date("2017-03-02")

aque_after_date_annual <- as_date("2018-05-10")

aque_after_date_continual <- as_date("2017-08-16")


# ⟞ b. Naples (NAPL) -----------------------------------------------------

napl_start_dates <- c("NAPL_CONTROL_2008-01-10", 
                      "NAPL_ANNUAL_2008-01-10", 
                      "NAPL_CONTINUAL_2010-04-27")

napl_start_date <- as_date("2008-01-10")

napl_after_date <- as_date("2016-02-19") # wrong in methods? 

napl_after_date_annual <- as_date("2017-05-16")

napl_after_date_continual <- as_date("2016-08-14")


# ⟞ c. Isla Vista (IVEE) --------------------------------------------------

ivee_start_dates <- c("IVEE_CONTROL_2011-10-26", 
                      "IVEE_ANNUAL_2011-10-26")

ivee_start_date <- as_date("2011-10-26")

ivee_after_date <- as_date("2016-02-18")

ivee_after_date_annual <- as_date("2017-05-15")

# ⟞ d. Mohawk (MOHK) ------------------------------------------------------

mohk_start_dates <- c("MOHK_ANNUAL_2008-01-17", 
                      "MOHK_CONTROL_2008-01-17", 
                      "MOHK_CONTINUAL_2010-05-05")

mohk_start_date <- as_date("2008-01-17")

mohk_after_date <- as_date("2017-02-13")

mohk_after_date_annual <- as_date("2018-05-15")

mohk_after_date_continual <- as_date("2017-08-11")

# ⟞ e. Carpinteria (CARP)  ------------------------------------------------

carp_start_dates <- c("CARP_CONTROL_2008-02-12", 
                      "CARP_ANNUAL_2008-02-12", 
                      "CARP_CONTINUAL_2010-04-23")

carp_start_date <- as_date("2008-02-12")

carp_after_date <- as_date("2017-02-15")

carp_after_date_annual <- as_date("2018-05-22")

carp_after_date_continual <- as_date("2017-08-10")

# 3. useful functions for wrangling ---------------------------------------

# create a column for "after" experimental removal
after_dates_column <- function(df) {
  df %>% 
    mutate(after_dates = case_when(
      site == "aque" & treatment == "annual" ~ aque_after_date_annual,
      site == "aque" & treatment == "continual" ~ aque_after_date_continual,
      site == "aque" & treatment == "control" ~ aque_after_date_annual,
      site == "napl" & treatment == "annual" ~ napl_after_date_annual,
      site == "napl" & treatment == "continual" ~ napl_after_date_continual,
      site == "napl" & treatment == "control" ~ napl_after_date_annual,
      site == "ivee" & treatment == "annual" ~ ivee_after_date_annual,
      site == "mohk" & treatment == "annual" ~ mohk_after_date_annual,
      site == "mohk" & treatment == "continual" ~ mohk_after_date_continual,
      site == "mohk" & treatment == "control" ~ mohk_after_date_annual,
      site == "carp" & treatment == "annual" ~ carp_after_date_annual,
      site == "carp" & treatment == "continual" ~ carp_after_date_continual,
      site == "carp" & treatment == "control" ~ carp_after_date_annual
    ))
}

# make a new column for during and after and set factor levels
exp_dates_column <- function(df) {
  df %>% 
    mutate(exp_dates = case_when(
      # after for annual removal:
      site == "aque" & treatment == "annual" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "annual" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "annual" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "annual" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "annual" & date > carp_after_date_annual ~ "after",
      # after for continual removal:
      site == "aque" & treatment == "continual" & date > aque_after_date_continual ~ "after",
      site == "napl" & treatment == "continual" & date > napl_after_date_continual ~ "after",
      site == "mohk" & treatment == "continual" & date > mohk_after_date_continual ~ "after",
      site == "carp" & treatment == "continual" & date > carp_after_date_continual ~ "after",
      # after for control:
      site == "aque" & treatment == "control" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "control" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "control" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "control" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "control" & date > carp_after_date_annual ~ "after",
      # everything else is "during" the experiment
      TRUE ~ "during"
    ),
    exp_dates = fct_relevel(exp_dates, c("during", "after")))  
}

# annual removal: make a new column for during and after and set factor levels
exp_dates_column_annual <- function(df) {
  df %>% 
    mutate(exp_dates = case_when(
      # after for annual removal:
      site == "aque" & treatment == "annual" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "annual" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "annual" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "annual" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "annual" & date > carp_after_date_annual ~ "after",
      # after for control:
      site == "aque" & treatment == "control" & date > aque_after_date_annual ~ "after",
      site == "napl" & treatment == "control" & date > napl_after_date_annual ~ "after",
      site == "ivee" & treatment == "control" & date > ivee_after_date_annual ~ "after",
      site == "mohk" & treatment == "control" & date > mohk_after_date_annual ~ "after",
      site == "carp" & treatment == "control" & date > carp_after_date_annual ~ "after",
      # everything else is "during" the experiment
      TRUE ~ "during"
    ),
    exp_dates = fct_relevel(exp_dates, c("during", "after")))  
}

# continual removal: make a new column for during and after and set factor levels
exp_dates_column_continual <- function(df) {
  df %>% 
    mutate(exp_dates = case_when(
      # after for continual removal:
      site == "aque" & treatment == "continual" & date > aque_after_date_continual ~ "after",
      site == "napl" & treatment == "continual" & date > napl_after_date_continual ~ "after",
      site == "mohk" & treatment == "continual" & date > mohk_after_date_continual ~ "after",
      site == "carp" & treatment == "continual" & date > carp_after_date_continual ~ "after",
      # after for control:
      site == "aque" & treatment == "control" & date > aque_after_date_continual ~ "after",
      site == "napl" & treatment == "control" & date > napl_after_date_continual ~ "after",
      site == "ivee" & treatment == "control" & date > ivee_after_date_continual ~ "after",
      site == "mohk" & treatment == "control" & date > mohk_after_date_continual ~ "after",
      site == "carp" & treatment == "control" & date > carp_after_date_continual ~ "after",
      # everything else is "during" the experiment
      TRUE ~ "during"
    ),
    exp_dates = fct_relevel(exp_dates, c("during", "after")))  
}

# create a new column for season and set factor levels
season_column <- function(df) {
  df %>% 
    mutate(season = case_when(
      month %in% c(12, 1, 2, 3) ~ "winter",
      month %in% c(4, 5) ~ "spring",
      month %in% c (6, 7, 8) ~ "summer",
      month %in% c(9, 10, 11) ~ "fall"
    ),
    season = fct_relevel(season, "spring", "summer", "fall", "winter")) 
}

# create a new column for new groupings
new_group_column <- function(df) {
  df %>% 
    mutate(new_group = case_when(
      group == "algae" ~ "algae",
      group == "fish" ~ "fish",
      group == "invert" & taxon_family != "Pholadidae" & mobility == "sessile" ~ "epi_inverts",
      taxon_family == "Pholadidae" ~ "endo_inverts",
      sp_code %in% c("OPSP", "SPL", "LIGL", "SFL", "MECR") ~ "herb_inverts",
      sp_code %in% c("COCA", "AML", "PGL", "PAIN", "PHL", "LOGR", "DIL", "KEKE", "PBL") ~ "carn_inverts"
    ))
}

# function to calculate standard error
se <- function(x,...){
  sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
}

# 4. data -----------------------------------------------------------------

# ⟞ a. Max's guild data ---------------------------------------------------

guilds <- read_csv(here::here("code/castorani", "LTE_guild_data.csv")) %>% 
  mutate(sp.code = replace_na(sp.code, "Nandersoniana")) %>% 
  rename("new_group" = biomass.guild)


# ⟞ b. LTE all species biomass --------------------------------------------

biomass <- read_csv(here::here("data", "LTE_All_Species_Biomass_at_transect_20220314.csv")) %>% 
  clean_names() %>% 
  # ANOB is incorrectly coded as having "SESSILE" mobility
  mutate(mobility = replace(mobility, sp_code == "ANOB", "MOBILE")) %>% 
  # replace NA sp_code with Nandersoniana
  mutate(sp_code = case_when(
    scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
    TRUE ~ sp_code
  )) %>% 
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA)) %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "treatment", "site"), str_to_lower) %>% 
  # # make a new column for after dates
  # after_dates_column() %>% 
  # make a new column for during and after and set factor levels
  exp_dates_column() %>% 
  # create a new column for season and set factor levels
  season_column() %>% 
  # new group column %>% 
  left_join(., guilds, by = c("sp_code" = "sp.code")) %>% 
  # take out all the first dates 
  filter(!(sample_ID %in% c(aque_start_dates, napl_start_dates, ivee_start_dates, mohk_start_dates, carp_start_dates)))

# ⟞ c. LTE algae biomass at section ---------------------------------------

biomass_section <- read_csv(here::here("data", "LTE_Algae_Biomass_at_section_20220422.csv")) %>% 
  clean_names() %>% 
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA)) %>% 
  # create a sample_ID for each sampling date at each treatment at each site at each section
  unite("sample_ID", site, treatment, date, quad, side, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "treatment", "site"), str_to_lower) %>% 
  # make a new column for during and after and set factor levels
  exp_dates_column() %>% 
  # create a new column for season and set factor levels
  season_column() %>% 
  # take out all the first dates 
  filter(!(sample_ID %in% c(aque_start_dates, napl_start_dates, ivee_start_dates, mohk_start_dates, carp_start_dates)))

# ⟞ d. LTE percent cover --------------------------------------------------

percov <- read_csv(here::here("data", "LTE_Cover_All_Years_20200605.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "site", "treatment"), str_to_lower) %>% 
  # make a new column for after dates
  after_dates_column() %>% 
  # make a new column for during and after and set factor levels
  exp_dates_column() %>% 
  # create a new column for season and set factor levels
  season_column()


# ⟞ e. annual benthics percent cover --------------------------------------

percov_annual <- read_csv(here::here("data/benthics", "Annual_Cover_All_Years_20210108.csv"))


# ⟞ f. annual benthics biomass --------------------------------------------

biomass_annual <- read_csv(here::here("data/benthics", "Annual_All_Species_Biomass_at_transect_20210108.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "site"), str_to_lower)

# ⟞ g. LTE kelp fronds ----------------------------------------------------

kelp_fronds <- read_csv(here::here("data", "LTE_Kelp_All_Years_20220202.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "treatment", "site"), str_to_lower) %>% 
  # # make a new column for after dates
  # after_dates_column() %>% 
  # make a new column for during and after and set factor levels
  exp_dates_column() %>% 
  # create a new column for season and set factor levels
  season_column() %>% 
  # replace -99999 values with NA
  mutate(fronds = replace(fronds, fronds < 0, NA))

# 5. operators ------------------------------------------------------------

# not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# 6. useful vectors and data frames ---------------------------------------


# ⟞ a. species -----------------------------------------------------------

# data frame of each species including scientific name, common name, and group
spp_names <- biomass %>% 
  dplyr::select(group, sp_code, scientific_name, common_name) %>% 
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

# from rank abundance curves
algae_selection <- c("PTCA", "DL", "CYOS", "CC", "R", "EC", "POLA", "RAT", "CO", "EGME", "BO", "AU")


# ⟞ b. sites --------------------------------------------------------------

# vector of sites
LTE_sites <- c("aque", "napl", "ivee", "mohk", "carp")

LTER_sites <- biomass_annual %>% 
  dplyr::select(site) %>% 
  mutate(site = fct_relevel(site, c("bull", "aque", "ahnd", "napl", "ivee", "golb", "abur", "mohk", "carp", "scdi", "sctw"))) %>% 
  unique() %>% 
  pull(site)

site_quality <- tribble(
  ~site, ~quality,
  "mohk", "high",
  "ivee", "high",
  "aque", "medium",
  "napl", "medium",
  "carp", "low"
) %>% 
  mutate(quality = fct_relevel(quality, c("low", "medium", "high")))

# ⟞ c. sample IDs ---------------------------------------------------------

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


# ⟞ d. today's date -------------------------------------------------------

todays_date <- Sys.Date()


# 7.  plot aesthetics -----------------------------------------------------

# ⟞ a. site colors and shapes --------------------------------------------

aque_col <- "#FDD989"
carp_col <- "#516238"
ivee_col <- "#4CA2B0"
mohk_col <- "#985E5C"
napl_col <- "#D98A63"

color_palette_site <- c("aque" = aque_col, 
                        "napl" = napl_col, 
                        "mohk" = mohk_col, 
                        "carp" = carp_col)

aque_shape <- 21
napl_shape <- 22
mohk_shape <- 23
carp_shape <- 24

shape_palette_site <- c("aque" = aque_shape, 
                        "napl" = napl_shape, 
                        "mohk" = mohk_shape, 
                        "carp" = carp_shape)

annual_col <- "#54662C"
continual_col <- "#009BB0"
control_col <- "grey"
kelp_col <- "#6D5A18"
under_col <- "#CC7540"

# ⟞ b. treatment colors and shapes ----------------------------------------

color_palette <- c("annual" = annual_col, 
                   "continual" = continual_col, 
                   "control" = control_col)

annual_shape <- 19
continual_shape <- 17
control_shape <- 21

shape_palette <- c("annual" = annual_shape, 
                   "continual" = continual_shape, 
                   "control" = control_shape)

# ⟞ c. site full names ----------------------------------------------------

aque_full <- "Arroyo Quemado (AQUE)"
napl_full <- "Naples (NAPL)"
ivee_full <- "Isla Vista (IVEE)"
mohk_full <- "Mohawk (MOHK)"
carp_full <- "Carpinteria (CARP)"

sites_full_levels <- c("Arroyo Quemado (AQUE)", 
                       "Naples (NAPL)", 
                       "Isla Vista (IVEE)", 
                       "Mohawk (MOHK)", 
                       "Carpinteria (CARP)")

sites_full <- setNames(c("Arroyo Quemado (AQUE)",
                         "Naples (NAPL)",
                         "Isla Vista (IVEE)",
                         "Mohawk (MOHK)",
                         "Carpinteria (CARP)"),
                       c("aque",
                         "napl",
                         "ivee",
                         "mohk",
                         "carp"))

sites_continual_full <- setNames(c("Arroyo Quemado (AQUE)",
                         "Naples (NAPL)",
                         "Mohawk (MOHK)",
                         "Carpinteria (CARP)"),
                       c("aque",
                         "napl",
                         "mohk",
                         "carp"))








