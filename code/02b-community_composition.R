##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

##########################################################################-
# 1. data frames and wranging functions -----------------------------------
##########################################################################-

# ⟞ a. continual removal communities --------------------------------------


comm_df <- biomass %>% 
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
  # only include 2 year sampling sites
  drop_na(comp_2yrs)

# metadata for all plots
comm_meta <- comm_df %>% 
  select(sample_ID, site, date, year, month, treatment, exp_dates, quarter, time_yrs, time_since_start, time_since_end, kelp_year, comp_2yrs, comp_3yrs, quality, site_full) %>% 
  unique()

# metadata for continual removal plots
comm_meta_continual <- comm_meta %>% 
  filter(treatment == "continual")

# metadata for control plots
comm_meta_control <- comm_meta %>% 
  filter(treatment == "control") %>% 
  # weird point?
  filter(sample_ID != "mohk_control_2020-08-12") 

# ⟞ b. widening function --------------------------------------------------

widen <- function(df) {
  df %>% 
    # select columns of interest
    select(sample_ID, sp_code, dry_gm2) %>% 
    # get into wide format for community analysis 
    pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
    # make sample_ID column row names
    column_to_rownames("sample_ID") %>% 
    replace(is.na(.), 0)
}

# ⟞ c. algae --------------------------------------------------------------

comm_mat_algae <- comm_df %>% 
  # algae only
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  # select columns of interest
  select(sample_ID, sp_code, dry_gm2) %>% 
  # get into wide format for community analysis
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  # make the sample_ID column row names
  column_to_rownames("sample_ID") %>% 
  replace(is.na(.), 0)

comm_meta_algae <- comm_meta %>% 
  filter(sample_ID %in% rownames(comm_mat_algae))

comm_mat_continual_algae <- comm_df %>% 
  filter(new_group == "algae") %>% 
  # only include sampling from removal plots
  filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
  widen()

comm_mat_control_algae <- comm_df %>% 
  filter(new_group == "algae") %>% 
  # only include sampling from control plots
  filter(sample_ID %in% comm_meta_control$sample_ID) %>% 
  widen()

# ⟞ d. epilithic ----------------------------------------------------------

comm_mat_epi <- comm_df %>% 
  # eplithic inverts only
  filter(new_group == "epi_inverts") %>% 
  widen()

comm_mat_continual_epi <- comm_df %>% 
  filter(new_group == "epi_inverts") %>% 
  filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
  widen()

comm_mat_control_epi <- comm_df %>% 
  filter(new_group == "epi_inverts") %>% 
  filter(sample_ID %in% comm_meta_control$sample_ID) %>% 
  widen()

##########################################################################-
# 2. ordination -----------------------------------------------------------
##########################################################################-

# ⟞ a. algae --------------------------------------------------------------

# ⟞ ⟞ i. period*treatment -------------------------------------------------

algae_pt_dist <- vegdist(comm_mat_algae, "bray")
algae_pt_nmds <- metaMDS(comm_mat_algae, "bray")
stressplot(algae_pt_nmds)

plot(algae_pt_nmds, dispaly)

algae_pt_plotdf <- scores(algae_pt_nmds, display = "sites") %>%
  as_tibble(rownames = "sample_ID") %>% 
  left_join(., comm_meta, by = "sample_ID")

simper_algae <- simper(comm_mat_algae, comm_meta_algae$treatment)
summary(simper_algae)
simper_algae_df <- simper_algae$control_continual %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average)
# pull top species from simper analysis
simper_algae_spp <- simper_algae_df %>% 
  # take top 10 species
  head(10) %>% 
  pull(species)

algae_pt_species <- scores(algae_pt_nmds, display = "species", tidy = TRUE) %>% 
  as_tibble(rownames = "sp_code") %>% 
  filter(sp_code %in% simper_algae_spp) %>% 
  left_join(., spp_names, by = "sp_code")

algae_pt_perma <- adonis2(comm_mat_algae ~ treatment*comp_2yrs, data = comm_meta)
algae_pt_perma

algae_betadisper <- betadisper(algae_pt_dist, comm_meta$treatment)
anova(algae_betadisper)

algae_pt_continual_plot <- algae_pt_plotdf %>% 
  filter(treatment == "continual") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  # make sure plot is square (ratio comes from axis limits)
  coord_fixed(ratio = 3.75/2.25) +
  scale_x_continuous(limits = c(-1.75, 2)) +
  scale_y_continuous(limits = c(-1.25, 1), breaks = seq(-1, 1, by = 1)) +
  # x and y axes at 0
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  # site points
  geom_point(aes(shape = site_full, fill = comp_2yrs), size = 4, alpha = 0.9) +
  # ellipse
  stat_ellipse(aes(color = comp_2yrs), size = 1) +
  # colors and linetypes
  scale_color_manual(values = c(start_col, during_col, after_col)) +
  scale_fill_manual(values = c(start_col, during_col, after_col), guide = "none") +
  # scale_linetype_manual(values = c("continual" = 1)) +
  scale_shape_manual(values = c(aque_shape, napl_shape, mohk_shape, carp_shape)) +
  # arrows for species from SIMPER
  geom_text_repel(data = algae_pt_species, 
                  aes(x = NMDS1, y = NMDS2, 
                      label = stringr::str_wrap(scientific_name, 10, width = 40)),
                  color = "#C70000", lineheight = 0.8) +
  geom_segment(data = algae_pt_species,
               aes(x = 0, y = 0,
                   xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.5, "cm")), 
               color = "#C70000", size = 1) +
  # plot aesthetics
  theme_classic() +
  # scale_x_continuous(limits = c(-1.75, 2)) +
  # scale_y_continuous(limits = c(-1.25, 1), breaks = seq(-1, 1, by = 1)) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.position = c(-1.2, 1.2)) +
  labs(shape = "Site",
       color = "Time period",
       fill = "Time period") +
  # ellipse labels. clown shit
  # annotate("text", x = -1.2, y = 1.05, label = "Start of", size = 10, col = start_col) +
  # annotate("text", x = -1.2, y = 0.9, label = "removal", size = 10, col = start_col) +
  # annotate("text", x = 1.7, y = 0.85, label = "End of", size = 10, col = during_col) +
  # annotate("text", x = 1.7, y = 0.7, label = "removal", size = 10, col = during_col) +
  # annotate("text", x = -1.4, y = 1.05, label = "Recovery", size = 10, col = after_col) +
  # annotate("text", x = -1.4, y = 0.9, label = "period", size = 10, col = after_col) +
  # stress annotation
  annotate("text", x = -1.4, y = -1.25, label = "Stress = 0.2", size = 5)
algae_pt_continual_plot

algae_pt_control_plot <- algae_pt_plotdf %>% 
  filter(treatment == "control") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  # make sure plot is square (ratio comes from axis limits)
  # coord_fixed(ratio = 3.75/2.25) +
  # scale_x_continuous(limits = c(-1.75, 2)) +
  # scale_y_continuous(limits = c(-1.25, 1), breaks = seq(-1, 1, by = 1)) +
  # x and y axes at 0
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  # site points
  geom_point(aes(shape = site_full, fill = comp_2yrs), size = 4, alpha = 0.9) +
  # ellipse
  stat_ellipse(aes(color = comp_2yrs), size = 1) +
  # colors and linetypes
  scale_color_manual(values = c(start_col, during_col, after_col)) +
  scale_fill_manual(values = c(start_col, during_col, after_col), guide = "none") +
  # scale_linetype_manual(values = c("continual" = 1)) +
  scale_shape_manual(values = c(aque_shape, napl_shape, mohk_shape, carp_shape)) +
  # arrows for species from SIMPER
  geom_text_repel(data = algae_pt_species, 
                  aes(x = NMDS1, y = NMDS2, 
                      label = stringr::str_wrap(scientific_name, 10, width = 40)),
                  color = "#C70000", lineheight = 0.8) +
  geom_segment(data = algae_pt_species,
               aes(x = 0, y = 0,
                   xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.5, "cm")), 
               color = "#C70000", size = 1) +
  # plot aesthetics
  theme_classic() +
  # scale_x_continuous(limits = c(-1.75, 2)) +
  # scale_y_continuous(limits = c(-1.25, 1), breaks = seq(-1, 1, by = 1)) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.position = c(-1.2, 1.2)) +
  labs(shape = "Site",
       color = "Time period",
       fill = "Time period") +
  # ellipse labels. clown shit
  # annotate("text", x = -1.2, y = 1.05, label = "Start of", size = 10, col = start_col) +
  # annotate("text", x = -1.2, y = 0.9, label = "removal", size = 10, col = start_col) +
  # annotate("text", x = 1.7, y = 0.85, label = "End of", size = 10, col = during_col) +
  # annotate("text", x = 1.7, y = 0.7, label = "removal", size = 10, col = during_col) +
  # annotate("text", x = -1.4, y = 1.05, label = "Recovery", size = 10, col = after_col) +
  # annotate("text", x = -1.4, y = 0.9, label = "period", size = 10, col = after_col) +
  # stress annotation
  annotate("text", x = -1.4, y = -1.25, label = "Stress = 0.2", size = 5)
algae_pt_control_plot

