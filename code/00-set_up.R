
##########################################################################-
# 1. libraries ------------------------------------------------------------
##########################################################################-

# general cleaning, visualization, etc.
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(patchwork)
library(calecopal)
library(PNWColors)
library(ggrepel)
library(plotly)
library(fuzzyjoin)
library(gt)
library(rlang)
library(multcompView)
library(ggeffects)
library(ggnewscale)
library(scales)
library(ggpubr)

# analysis
library(vegan)
library(vegclust)
library(ecotraj)
library(FD)
# library(BiodiversityR)
library(minpack.lm) # Fitting non-linear models
library(nls2) # Fitting non-linear models
library(AICcmodavg) # calculate second order AIC (AICc)
library(MuMIn)
library(boot)
library(lmerTest) # also loads `lme4`
library(glmmTMB)
library(nlme)
library(DHARMa)
library(performance)
library(emmeans)
library(gtsummary)

##########################################################################-
# 2. start and end dates --------------------------------------------------
##########################################################################-

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

##########################################################################-
# 3. useful wranging functions --------------------------------------------
##########################################################################-

# create a column for "after" experimental removal
after_dates_column <- function(df) {
  df %>% 
    mutate(after_dates = case_when(
      site == "aque" & treatment == "annual" ~ aque_after_date_annual,
      site == "aque" & treatment == "continual" ~ aque_after_date_continual,
      site == "aque" & treatment == "control" ~ aque_after_date_continual,
      site == "napl" & treatment == "annual" ~ napl_after_date_annual,
      site == "napl" & treatment == "continual" ~ napl_after_date_continual,
      site == "napl" & treatment == "control" ~ napl_after_date_continual,
      site == "ivee" & treatment == "annual" ~ ivee_after_date_annual,
      site == "mohk" & treatment == "annual" ~ mohk_after_date_annual,
      site == "mohk" & treatment == "continual" ~ mohk_after_date_continual,
      site == "mohk" & treatment == "control" ~ mohk_after_date_continual,
      site == "carp" & treatment == "annual" ~ carp_after_date_annual,
      site == "carp" & treatment == "continual" ~ carp_after_date_continual,
      site == "carp" & treatment == "control" ~ carp_after_date_continual
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

# time since start and time since end
# expects data frame with the following columns:
# site, year, month, date, control, annual, delta_annual
time_since_columns_annual <- function(df) {
  df %>% 
    # create a column for quarter
    mutate(quarter = case_when(
      month <= 3 ~ "Q1",
      month <= 6 ~ "Q2",
      month <= 9 ~ "Q3",
      TRUE ~ "Q4"
    )) %>% 
    # calculate time since start of experiment
    mutate(time_yrs = case_when(
      # AQUE, NAPL, MOHK, CARP: control and annual started in 2008
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q1" ~ year + 0.125 - 2008,
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q2" ~ year + 0.375 - 2008,
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q3" ~ year + 0.625 - 2008,
      site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q4" ~ year + 0.875 - 2008, 
      # IVEE control and annual started in 2011
      site == "ivee" & quarter == "Q1" ~ year + 0.125 - 2011,
      site == "ivee" & quarter == "Q2" ~ year + 0.375 - 2011,
      site == "ivee" & quarter == "Q3" ~ year + 0.625 - 2011,
      site == "ivee" & quarter == "Q4" ~ year + 0.875 - 2011
    )) %>% 
    group_by(site) %>% 
    mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
    ungroup() %>% 
    # calculate time since end of experiment
    group_by(site, exp_dates) %>% 
    # if "after", then simple: the time in years - the minimum time in years
    # if "during", then more complex: take the max time in years and add 0.25, then subtract the time in years
    mutate(time_since_end = case_when(
      exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
      exp_dates == "after" ~ time_yrs - min(time_yrs)
    )) %>% 
    ungroup()
}

time_since_columns_continual <- function(df) {
  df %>% 
    # create a column for quarter
    mutate(quarter = case_when(
      month <= 3 ~ "Q1",
      month >= 4 & month <= 6 ~ "Q2",
      month >= 7 & month <= 9 ~ "Q3",
      month >= 10 & month <= 12 ~ "Q4"
    )) %>% 
    # mutate(quarter = case_when(
    #   month <= 3 ~ "Q1",
    #   month <= 6 ~ "Q2",
    #   month <= 9 ~ "Q3",
    #   TRUE ~ "Q4"
    # )) %>% 
    # calculate time since start of experiment
    mutate(time_yrs = case_when(
      # AQUE, NAPL, MOHK, CARP: continual started in 2010
      quarter == "Q1" ~ year + 0.125 - 2010,
      quarter == "Q2" ~ year + 0.375 - 2010,
      quarter == "Q3" ~ year + 0.625 - 2010,
      quarter == "Q4" ~ year + 0.875 - 2010
    )) %>% 
    group_by(site) %>% 
    mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
    ungroup() %>% 
    # calculate time since end of experiment
    group_by(site, exp_dates) %>% 
    mutate(time_since_end = case_when(
      exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
      exp_dates == "after" ~ time_yrs - min(time_yrs)
    )) %>% 
    mutate(test_min_time_yrs = min(time_yrs)) %>% 
    ungroup()
}

# kelp year
# expects data frame that already has quarter from time_since_columns functions
kelp_year_column <- function(df) {
  df %>% 
  # create a new column for "kelp year"
  mutate(quarter = fct_relevel(quarter, "Q2", "Q3", "Q4", "Q1")) %>% 
    mutate(kelp_year = case_when(
      quarter %in% c("Q2", "Q3", "Q4") ~ year,
      quarter == "Q1" ~ year - 1
    )) %>% 
    mutate(kelp_year = paste("kelp_", kelp_year, "-", kelp_year + 1, sep = "")) 
}

# comparison column for annual removal
# first 2 or 3 years of start, last 2 or 3 years during experimental removal, most recent 2 or 3 years
comparison_column_annual <- function(df) {
  df %>% 
  mutate(comp_2yrs = case_when(
    site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010") ~ "start",
    site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during", 
    site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
    site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
    kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022") ~ "after"
  )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    # create a column for the points to compare for "3 year interval"
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010", "kelp_2010-2011") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2019-2020", "kelp_2020-2021", "kelp_2021-2022") ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    # add in full names of sites
    left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
    rename(site_full = value) %>% 
    mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Isla Vista", "Mohawk", "Carpinteria"))
}

comparison_column_continual <- function(df) {
  df %>% 
    # create a column for the points to compare for "2 year interval"
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022") ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    # create a column for the points to compare for "3 year interval"
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012", "kelp_2012-2013") ~ "start",
      site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during", 
      site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
      site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
      kelp_year %in% c("kelp_2019-2020", "kelp_2020-2021", "kelp_2021-2022") ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    # create a new sample ID that is site, date, quarter
    unite("sample_ID", site, date, quarter, remove = FALSE)
}

##########################################################################-
# 4. data -----------------------------------------------------------------
##########################################################################-

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

##########################################################################-
# 5. operators ------------------------------------------------------------
##########################################################################-

# not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

##########################################################################-
# 6. useful vectors and data frames ---------------------------------------
##########################################################################-

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
algae_spp_names <- biomass %>% 
  filter(group == "algae") %>% 
  dplyr::select(sp_code, scientific_name, common_name, taxon_phylum, taxon_class, taxon_order, taxon_family) %>% 
  unique() %>% 
  # fill in missing orders with phyla
  mutate(tax_group = case_when(
    taxon_order == -99999  ~ paste("Unidentified ", taxon_phylum, sep = ""),
    taxon_order = TRUE ~ taxon_order
  ))

algae_orders <- unique(algae_spp_names$tax_group)

epi_spp_names <- biomass %>% 
  filter(new_group == "epilithic.sessile.invert") %>% 
  dplyr::select(sp_code, scientific_name, common_name, taxon_phylum, taxon_class, taxon_order, taxon_family) %>% 
  unique() %>% 
  # fill in missing classes with phyla
  mutate(tax_group = case_when(
    taxon_class == -99999  ~ paste("Unidentified ", taxon_phylum, sep = ""),
    taxon_class = TRUE ~ taxon_class
  ))

epi_classes <- unique(epi_spp_names$tax_group)

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

site_quality <- tribble(
  ~site, ~quality,
  "mohk", "high",
  "ivee", "high",
  "aque", "medium",
  "napl", "medium",
  "carp", "low"
) %>% 
  mutate(quality = fct_relevel(quality, c("low", "medium", "high")))


# 7.  plot aesthetics -----------------------------------------------------

# ⟞ a. site colors and shapes --------------------------------------------

aque_col <- '#0f85a0'
napl_col <- '#dd4124'
ivee_col <- '#ed8b00'
mohk_col <- '#edd746'
carp_col <- '#00496f'

color_palette_site <- c("aque" = aque_col, 
                        "napl" = napl_col, 
                        "mohk" = mohk_col, 
                        "carp" = carp_col)

color_palette_site_annual <- c("aque" = aque_col, 
                               "napl" = napl_col, 
                               "ivee" = ivee_col,
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


# ⟞ d. start-during-after colors ------------------------------------------

start_col <- "#FFBF00"
during_col <- "#54662C"
after_col <- "#114C54"


# ⟞ e. site full names ----------------------------------------------------

aque_full <- "Arroyo Quemado"
napl_full <- "Naples"
ivee_full <- "Isla Vista"
mohk_full <- "Mohawk"
carp_full <- "Carpinteria"

sites_full_levels <- c("Arroyo Quemado", 
                       "Naples", 
                       "Isla Vista", 
                       "Mohawk", 
                       "Carpinteria")

sites_full <- setNames(c("Arroyo Quemado",
                         "Naples",
                         "Isla Vista",
                         "Mohawk",
                         "Carpinteria"),
                       c("aque",
                         "napl",
                         "ivee",
                         "mohk",
                         "carp"))

sites_continual_full <- setNames(c("Arroyo Quemado",
                         "Naples",
                         "Mohawk",
                         "Carpinteria"),
                       c("aque",
                         "napl",
                         "mohk",
                         "carp"))








