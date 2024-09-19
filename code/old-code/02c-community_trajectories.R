##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "02b-community_composition.R"))

##########################################################################-
# 1. algae trajectory -----------------------------------------------------
##########################################################################-

pre_df <- biomass %>% 
  new_group_column() %>% 
  # select columns of interest
  dplyr::select(site, year, month, treatment, date, new_group, sp_code, dry_gm2) %>% 
  exp_dates_column_continual() %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual() %>% 
  # only including continual + control sampling dates
  # filtered from kelp delta data frame created in upstream script
  filter(sample_ID %in% (delta_continual$sample_ID)) %>% 
  filter(treatment %in% c("control", "continual")) %>% 
  full_join(., site_quality, by = "site") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename(site_full = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>%
  # create new sample ID with treatment
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # algae only
  filter(new_group == "algae" & sp_code != "MAPY")

comm_df_all <- pre_df %>% 
  # select columns of interest
  select(sample_ID, sp_code, dry_gm2) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  # make the sample_ID column row names
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)

comm_meta_all <- pre_df %>% 
  select(sample_ID, site, date, year, month, treatment, exp_dates, quarter, time_yrs, time_since_start, time_since_end, kelp_year, quality, site_full) %>% 
  unique()

# bray distance object: algae_pt_dist
algae_pt_dist_all <- vegdist(comm_df_all, method = "bray")
algae_pt_pcoa <- pcoa(algae_pt_dist_all, correction = "none")

algae_pt_pcoa_vectors <- as_tibble(algae_pt_pcoa$vectors, rownames = "sample_ID") %>% 
  left_join(., comm_meta_all, by = "sample_ID") %>% 
  clean_names() %>% 
  group_by(kelp_year, site, treatment) %>% 
  summarize(axis_1 = mean(axis_1),
            axis_2 = mean(axis_2)) 

aque_contin <- algae_pt_pcoa_vectors %>% filter(site == "aque" & treatment == "continual")
napl_contin <- algae_pt_pcoa_vectors %>% filter(site == "napl" & treatment == "continual")
mohk_contin <- algae_pt_pcoa_vectors %>% filter(site == "mohk" & treatment == "continual")
carp_contin <- algae_pt_pcoa_vectors %>% filter(site == "carp" & treatment == "continual")

aque_contr <- algae_pt_pcoa_vectors %>% filter(site == "aque" & treatment == "control")
napl_contr <- algae_pt_pcoa_vectors %>% filter(site == "napl" & treatment == "control")
mohk_contr <- algae_pt_pcoa_vectors %>% filter(site == "mohk" & treatment == "control")
carp_contr <- algae_pt_pcoa_vectors %>% filter(site == "carp" & treatment == "control")

algae_contin_segments <- ggplot() +
  geom_segment(data = aque_contin, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = aque_col) +
  annotate("text", x = -0.1, y = -0.15, label = "Arroyo Quemado", color = aque_col) +
  geom_segment(data = napl_contin, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = napl_col) +
  annotate("text", x = 0, y = -0.3, label = "Naples", color = napl_col) +
  geom_segment(data = mohk_contin, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = mohk_col) +
  annotate("text", x = 0.2, y = 0.4, label = "Mohawk", color = mohk_col) +
  geom_segment(data = carp_contin, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = carp_col) +
  annotate("text", x = -0.5, y = 0.1, label = "Carpinteria", color = carp_col) +
  labs(title = "Continual")
algae_contin_segments

algae_contr_segments <- ggplot() +
  geom_segment(data = aque_contr, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = aque_col) +
  geom_segment(data = napl_contr, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = napl_col) +
  geom_segment(data = mohk_contr, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = mohk_col) +
  geom_segment(data = carp_contr, aes(x = axis_1, y = axis_2, 
                                       xend = c(tail(axis_1, n = -1), NA), yend = c(tail(axis_2, n = -1), NA)),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = carp_col) +
  labs(title = "Control")
algae_contr_segments

algae_contin_segments + algae_contr_segments

segmentDistances(algae_pt_dist_all, comm_meta_all$site, comm_meta_all$sample_ID) 
trajectoryDistances(algae_pt_dist_all, comm_meta_all$site, comm_meta_all$sample_ID) 
trajectoryConvergence(algae_pt_dist_all, comm_meta_all$site, comm_meta_all$sample_ID) 

