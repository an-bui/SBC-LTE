
##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "00-set_up.R"))

##########################################################################-
# 1. data frames ----------------------------------------------------------
##########################################################################-

# ⟞ a. delta_annual -------------------------------------------------------

delta_annual <- biomass %>% 
  filter(sp_code == "MAPY" & treatment %in% c("control", "annual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
  # fill in values for sampling dates where control and annual were surveyed on different days
  mutate(control = case_when(
    site == "aque" & date == "2008-03-07" ~ 106.82000,
    site == "napl" & date == "2012-09-26" ~ 143.22088,
    TRUE ~ control
  )) %>% 
  mutate(delta_annual = annual - control) %>%  
  # missing dates are from AQUE 2008-03-05, NAPL 2012-09-25, NAPL 2008-10-10
  drop_na(delta_annual) %>% 
  mutate(exp_dates = case_when(
    # after for continual removal:
    site == "aque" & date >= aque_after_date_annual ~ "after",
    site == "napl" & date >= napl_after_date_annual ~ "after",
    site == "ivee" & date >= ivee_after_date_annual ~ "after", 
    site == "mohk" & date >= mohk_after_date_annual ~ "after",
    site == "carp" & date >= carp_after_date_annual ~ "after",
    # everything else is "during" the experiment
    TRUE ~ "during"
  ),
  exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
  time_since_columns_annual() %>% 
  kelp_year_column() %>% 
  comparison_column_annual() 

# ⟞ b. delta_continual ----------------------------------------------------

delta_continual <- biomass %>% 
  filter(sp_code == "MAPY" & treatment %in% c("control", "continual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
  mutate(delta_continual = continual - control) %>%  
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for continual removal:
    site == "aque" & date >= aque_after_date_continual ~ "after",
    site == "napl" & date >= napl_after_date_continual ~ "after",
    site == "mohk" & date >= mohk_after_date_continual ~ "after",
    site == "carp" & date >= carp_after_date_continual ~ "after",
    # everything else is "during" the experiment
    TRUE ~ "during"
  ),
  exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual() %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria"))

