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

# ⟞ ⟞ ⟞  1. analysis ------------------------------------------------------

# Bray-Curtis
algae_pt_bray <- metaMDS(comm_mat_algae, "bray")

# stress plot
stressplot(algae_pt_bray)

# preliminary plot
plot(algae_pt_bray)

# SIMPER analysis
simper_algae <- simper(comm_mat_algae, comm_meta_algae$treatment)
summary(simper_algae)

# permanova
algae_pt_perma_2yrs <- adonis2(comm_mat_algae ~ treatment*comp_2yrs, data = comm_meta)
algae_pt_perma_2yrs

algae_pt_perma_3yrs <- adonis2(comm_mat_algae ~ treatment*comp_3yrs, data = comm_meta)
algae_pt_perma_3yrs # same as 2 yrs

# beta dispersion
algae_pt_dist <- vegdist(comm_mat_algae, "bray")
algae_betadisper <- betadisper(algae_pt_dist, comm_meta$comp_2yrs)
anova(algae_betadisper)

# ⟞ ⟞ ⟞  2. plotting ------------------------------------------------------

# points into data frame for plotting
algae_pt_bray_plotdf <- scores(algae_pt_bray, display = "sites") %>% 
  as_tibble(rownames = "sample_ID") %>% 
  # join with metadata
  left_join(., comm_meta, by = "sample_ID")

# create data frame from simper analysis
simper_algae_df <- simper_algae$control_continual %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average)

# pull top species from simper analysis
simper_algae_spp <- simper_algae_df %>% 
  # take top 10 species
  head(10) %>% 
  pull(species)

# get species points
algae_pt_bray_species <- scores(algae_pt_bray, display = "species", tidy = TRUE) %>% 
  as_tibble(rownames = "sp_code") %>% 
  # keep species from simper analysis only
  filter(sp_code %in% simper_algae_spp) %>% 
  left_join(., spp_names, by = "sp_code")

# continual removal plots only
algae_pt_bray_continual_plot <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "continual", algae_pt_bray_species
  ) +
  # axis limits
  scale_x_continuous(limits = c(-1.5, 1.6), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.6, 1.4), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) +
  # scale_x_continuous(limits = c(-1.75, 1.4), breaks = seq(-1, 1, by = 1)) +
  # scale_y_continuous(limits = c(-1.7, 1.45), breaks = seq(-1, 1, by = 1)) +
  # ordination_theme() +
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "(a) Understory algae") +
  # ellipse labels. clown shit
  # annotate("text", x = -1, y = 0.5, label = "Start of", size = 10, col = start_col) +
  # annotate("text", x = -1, y = 0.4, label = "removal", size = 10, col = start_col) +
  # annotate("text", x = 1.4, y = 0.75, label = "End of", size = 10, col = during_col) +
  # annotate("text", x = 1.4, y = 0.65, label = "removal", size = 10, col = during_col) +
  # annotate("text", x = 0.85, y = -0.9, label = "Recovery", size = 10, col = after_col) +
  # annotate("text", x = 0.85, y = -1, label = "period", size = 10, col = after_col) +
  # stress annotation
  annotate("text", x = -1.15, y = -1.51, label = "Stress = 0.2", size = 2) +
  theme(legend.position = "none")
algae_pt_bray_continual_plot


# control plots only
algae_pt_bray_control_plot <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "control", algae_pt_bray_species
) +
  # plot aesthetics
  # ordination_theme() +
  scale_x_continuous(limits = c(-1.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.85, 1.5), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) +
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "(a) Understory algae") +
  # ellipse labels. clown shit
  # annotate("text", x = -1.2, y = 1.05, label = "Start of", size = 10, col = start_col) +
  # annotate("text", x = -1.2, y = 0.9, label = "removal", size = 10, col = start_col) +
  # annotate("text", x = 1.7, y = 0.85, label = "End of", size = 10, col = during_col) +
  # annotate("text", x = 1.7, y = 0.7, label = "removal", size = 10, col = during_col) +
  # annotate("text", x = -1.4, y = 1.05, label = "Recovery", size = 10, col = after_col) +
  # annotate("text", x = -1.4, y = 0.9, label = "period", size = 10, col = after_col) +
  # stress annotation
  annotate("text", x = -0.9, y = -1.75, label = "Stress = 0.2", size = 2) +
  theme(legend.position = "none")
algae_pt_bray_control_plot

# both treatments together
algae_pt_bray_both_plot <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "both", algae_pt_bray_species
) +
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
algae_pt_bray_both_plot

# ⟞ b. epi inverts --------------------------------------------------------

# ⟞ ⟞ i. period*treatment -------------------------------------------------

# ⟞ ⟞ ⟞  1. analysis ------------------------------------------------------

# Bray-Curtis
epi_pt_bray <- metaMDS(comm_mat_epi, "bray")

# stress plot
stressplot(epi_pt_bray)

# preliminary plot
plot(epi_pt_bray)

# SIMPER analysis
simper_epi <- simper(comm_mat_epi, comm_meta$treatment)
summary(simper_epi)

# permanova
epi_pt_perma_2yrs <- adonis2(comm_mat_epi ~ treatment*comp_2yrs, data = comm_meta)
epi_pt_perma_2yrs

epi_pt_perma_3yrs <- adonis2(comm_mat_epi ~ treatment*comp_3yrs, data = comm_meta)
epi_pt_perma_3yrs # same as 2 years

# beta dispersion
epi_pt_dist <- vegdist(comm_mat_epi, "bray")
epi_betadisper <- betadisper(epi_pt_dist, comm_meta$comp_2yrs)
anova(epi_betadisper)
# not significantly different dispersions

# ⟞ ⟞ ⟞  2. plotting ------------------------------------------------------

# points into data frame for plotting
epi_pt_bray_plotdf <- scores(epi_pt_bray, display = "sites") %>% 
  as_tibble(rownames = "sample_ID") %>% 
  # join with metadata
  left_join(., comm_meta, by = "sample_ID")

# create data frame from simper analysis
simper_epi_df <- simper_epi$control_continual %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average)

# pull top species from simper analysis
simper_epi_spp <- simper_epi_df %>% 
  # take top 10 species
  head(10) %>% 
  pull(species)

# get species points
epi_pt_bray_species <- scores(epi_pt_bray, display = "species", tidy = TRUE) %>% 
  as_tibble(rownames = "sp_code") %>% 
  # keep species from simper analysis only
  filter(sp_code %in% simper_epi_spp) %>% 
  left_join(., spp_names, by = "sp_code")

# continual removal plots only
epi_pt_bray_continual_plot <- nmds_plot_fxn(
  epi_pt_bray_plotdf, "continual", epi_pt_bray_species
) +
  coord_fixed() +
  # scale_x_continuous(limits = c(-1.7, 1.4), breaks = seq(-1, 1, by = 1)) + # length = 3.1
  # scale_y_continuous(limits = c(-1.3, 1), breaks = seq(-1, 1, by = 1)) +
  scale_x_continuous(limits = c(-1.75, 1.4), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.7, 1.45), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) + # length = 3.15
  annotate("text", x = -1.4, y = -1.6, label = "Stress = 0.2", size = 2) +
  theme(legend.position = "none") +
  labs(title = "(b) Epilithic invertebrates") +
  theme(legend.position = "right")
epi_pt_bray_continual_plot


# control plots only
epi_pt_bray_control_plot <- nmds_plot_fxn(
  epi_pt_bray_plotdf, "control", epi_pt_bray_species
) +
  scale_x_continuous(limits = c(-1.75, 1.4), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.7, 1.45), breaks = seq(-1, 1, by = 1), expand = c(0, 0)) +
  annotate("text", x = -1.2, y = -1.6, label = "Stress = 0.2", size = 2) +
  labs(title = "(b) Epilithic invertebrates") 
epi_pt_bray_control_plot


# both treatments together
epi_pt_bray_both_plot <- nmds_plot_fxn(
  epi_pt_bray_plotdf, "both", epi_pt_bray_species
) +
  scale_x_continuous(limits = c(-1.75, 1.4), breaks = seq(-1, 1, by = 1)) +
  scale_y_continuous(limits = c(-1.7, 1.45), breaks = seq(-1, 1, by = 1)) +
  annotate("text", x = -1, y = -1.7, label = "Stress = 0.2", size = 2)
epi_pt_bray_both_plot

##########################################################################-
# 3. manuscript tables ----------------------------------------------------
##########################################################################-

anova_2yrs_tables <- rbind(anova_summary_fxn(algae_pt_perma_2yrs), anova_summary_fxn(epi_pt_perma_2yrs)) %>% 
  rename_with(., .fn = ~paste(., "_2yrs", sep = "", .cols = everything(cols)))
anova_2yrs_tables 

anova_3yrs_tables <- rbind(anova_summary_fxn(algae_pt_perma_3yrs), anova_summary_fxn(epi_pt_perma_3yrs)) %>% 
  rename_with(., .fn = ~paste(., "_3yrs", sep = "", .cols = everything(cols)))
anova_3yrs_tables

anova_together_tables <- cbind(anova_2yrs_tables, anova_3yrs_tables) %>% 
  # take out unwanted columns
  select(!c("model_2yrs1", "model_3yrs1", 
            "SumOfSqs_2yrs1", "SumOfSqs_3yrs1",
            "R2_2yrs1", "R2_3yrs1")) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  # group labels
  tab_row_group(
    label = "Epilithic invertebrates", rows = 6:10
  ) %>% 
  tab_row_group(
    label = "Algae", rows = 1:5
  ) %>% 
  # 2 and 3 year comparisons
  tab_spanner(
    label = "2 year comparison",
    columns = c(variables_2yrs1, Df_2yrs1, F_2yrs1, p_2yrs1)
  ) %>% 
  tab_spanner(
    label = "3 year comparison",
    columns = c(variables_3yrs1, Df_3yrs1, F_3yrs1, p_3yrs1)
  ) %>% 
  # change column names
  cols_label(
    variables_2yrs1 = "Source of variation",
    Df_2yrs1 = "df",
    F_2yrs1 = "pseudo-F",
    p_2yrs1 = "p", 
    variables_3yrs1 = "Source of variation",
    Df_3yrs1 = "df",
    F_3yrs1 = "pseudo-F",
    p_3yrs1 = "p"
  ) %>% 
  # increase spacing between cells
  tab_style(
    style = "padding-left:15px;padding-right:15px;",
    locations = cells_body()
  ) %>% 
  # align columns
  cols_align(columns = everything(),
             align = "center") %>% 
  # bold p < 0.05
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = p_2yrs1,
      rows = p_2yrs1 < 0.05
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = p_3yrs1,
      rows = p_3yrs1 < 0.05
    )
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>% 
  sub_missing(
    columns = everything(),
    missing_text = "--"
  ) %>% 
  tab_options(table.font.names = "Times New Roman") 
anova_together_tables

# gtsave(anova_together_tables,
#        here::here("tables", "ms-tables", paste("tbl-S4_", today(), ".png", sep = "")),
#        vwidth = 1500, vheight = 1000)


##########################################################################-
# 4. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. continual removal --------------------------------------------------

comm_comp_together <- algae_pt_bray_continual_plot + epi_pt_bray_continual_plot 

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)

# ⟞ b. control ------------------------------------------------------------

comm_comp_control <- algae_pt_bray_control_plot + epi_pt_bray_control_plot 

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S7_", today(), ".jpg", sep = "")),
#        plot = comm_comp_control,
#        height = 10, width = 16, units = "cm",
#        dpi = 300)



