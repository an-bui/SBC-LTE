##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

##########################################################################-
# 1. data frames ----------------------------------------------------------
##########################################################################-

# ⟞ a. continual removal communities --------------------------------------

# community data frames
comm_df <- biomass %>% 
  filter(sp_code != "MAPY") %>% 
  new_group_column() %>% 
  dplyr::select(site, year, month, treatment, date, new_group, sp_code, dry_gm2) %>% 
  pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
  dplyr::select(-annual) %>% 
  mutate(delta_continual = continual - control) %>% 
  drop_na(delta_continual) %>% 
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
  rename(site_full = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  full_join(., site_quality, by = "site") %>% 
  dplyr::select(-delta_continual) %>% 
  pivot_longer(cols = c(control, continual), names_to = "treatment", values_to = "dry_gm2") %>% 
  # create sample_ID
  unite("sample_ID", site, treatment, date, sep = "_", remove = FALSE) %>% 
  filter(site != "ivee")

# metadata for continual removal plots
comm_meta_continual <- comm_df %>% 
  select(sample_ID:exp_dates, quarter:treatment, -sp_code, -new_group) %>% 
  filter(treatment == "continual") %>% 
  drop_na(comp_2yrs) %>% 
  # weird point?
  filter(sample_ID != "mohk_control_2020-08-12") %>% 
  unique()

# metadata for control plots
comm_meta_control <- comm_df %>% 
  select(sample_ID:exp_dates, quarter:treatment, -sp_code, -new_group) %>% 
  filter(treatment == "control") %>% 
  drop_na(comp_2yrs) %>% 
  # weird point?
  filter(sample_ID != "mohk_control_2020-08-12") %>% 
  unique()

# ⟞ b. algae --------------------------------------------------------------

comm_mat_continual_algae <- comm_df %>% 
  # algae only
  filter(new_group == "algae") %>% 
  # only include sampling from removal plots
  filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
  # select columns of interest
  select(sample_ID, sp_code, dry_gm2) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  # make the sample_ID column row names
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)

comm_mat_control_algae <- comm_df %>% 
  filter(new_group == "algae") %>% 
  # only include sampling from control plots
  filter(sample_ID %in% comm_meta_control$sample_ID) %>% 
  select(sample_ID, sp_code, dry_gm2) %>% 
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)

# ⟞ c. epilithic ----------------------------------------------------------

comm_mat_continual_epi_inverts <- comm_df %>% 
  filter(new_group == "epi_inverts") %>% 
  filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
  select(sample_ID, sp_code, dry_gm2) %>% 
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)

comm_mat_control_epi_inverts <- comm_df %>% 
  filter(new_group == "epi_inverts") %>% 
  filter(sample_ID %in% comm_meta_control$sample_ID) %>% 
  select(sample_ID, sp_code, dry_gm2) %>% 
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)











