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
  # only include 3 year sampling sites
  drop_na(comp_3yrs)

# metadata for all plots
comm_meta <- comm_df %>% 
  select(sample_ID, site, date, year, month, treatment, exp_dates, quarter, time_yrs, time_since_start, time_since_end, kelp_year, comp_1yr, comp_2yrs, comp_3yrs, quality, site_full) %>% 
  unique()

# metadata for continual removal plots
comm_meta_continual <- comm_meta %>% 
  filter(treatment == "continual")

# metadata for control plots
comm_meta_control <- comm_meta %>% 
  filter(treatment == "control")
# weird point?
# filter(sample_ID != "mohk_control_2020-08-12") 

# making a nested data frame
comm_df_nested <- comm_df %>% 
  # take out giant kelp
  filter(sp_code != "MAPY") %>% 
  # create nested columns
  nest(data = everything(), .by = new_group) %>% 
  # get data into community matrix format
  mutate(comm_mat = map(data, 
                        ~ # select columns of interest
                          select(.x, sample_ID, sp_code, dry_gm2) %>% 
                          # get into wide format for community analysis
                          pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
                          # make the sample_ID column row names
                          column_to_rownames("sample_ID") %>% 
                          replace(is.na(.), 0))) %>% 
  # get metadata
  mutate(comm_meta = map(data, 
                         ~ select(.x, sample_ID, site, date, year, month, treatment,
                                  exp_dates, quarter, time_yrs, time_since_start,
                                  time_since_end, kelp_year, comp_2yrs, comp_3yrs,
                                  quality, site_full) %>% 
                           unique()
  )) %>% 
  # find species that do occur (exclude species that never occur)
  mutate(keepspp = map(comm_mat, 
                       ~ colSums(.x, na.rm = FALSE) %>% 
                         enframe() %>% 
                         filter(value > 0) %>% 
                         pull(name))) %>% 
  # find surveys where there was actually something there
  mutate(keepsurveys = map(comm_mat,
                           ~ rowSums(.x, na.rm = FALSE) %>% 
                             enframe() %>% 
                             filter(value > 0) %>% 
                             pull(name))) %>% 
  # only retain species that do occur
  mutate(comm_mat_step1 = map2(comm_mat, keepspp,
                               ~ select(.x, .y))) %>% 
  # only retain surveys where there was something there
  mutate(comm_mat_filtered = map2(comm_mat_step1, keepsurveys,
                                  ~ rownames_to_column(.x, "surveys") %>% 
                                    filter(surveys %in% .y) %>% 
                                    column_to_rownames("surveys"))) %>% 
  # filter metadata
  mutate(comm_meta_filtered = map2(comm_meta, comm_mat_filtered, 
                                   ~ filter(.x, sample_ID %in% rownames(.y))))


comm_analyses <- comm_df_nested %>% 
  select(new_group, comm_mat_filtered, comm_meta_filtered) %>% 
  # bray-curtis: will get a warning for fish and endolithic inverts
  mutate(nmds_bray = map(comm_mat_filtered, 
                         ~ metaMDS(.x, "bray"))) %>% 
  # mutate(nmds_bray_plot = map(nmds_bray,
  #                             ~ plot(.x))) %>% 
  # # jaccard
  # mutate(nmds_jacc = map(comm_mat_filtered,
  #                         ~ metaMDS(.x, "jacc"))) %>% 
  # mutate(nmds_jacc_plot = map(nmds_jacc,
  #                             ~ plot(.x))) %>% 
  mutate(simper = map2(comm_mat_filtered, comm_meta_filtered,
                       ~ simper(.x, .y$treatment, ordered = TRUE)))

# putting in separate df because it takes a while
simper_analysis <- comm_analyses %>% 
  # create a vector of species identified as important in SIMPER analyses
  mutate(simper_spp = map(simper, 
                          ~ .x$control_continual$cusum %>% 
                            enframe() %>% 
                            head(10) %>% 
                            pull(name))) %>% 
  mutate(simper_spp_info = map(simper_spp,
                               ~ spp_names %>% 
                                 filter(sp_code %in% .x)))


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
  filter(new_group == "epilithic.sessile.invert") %>% 
  widen()

comm_mat_continual_epi <- comm_df %>% 
  filter(new_group == "epilithic.sessile.invert") %>% 
  filter(sample_ID %in% comm_meta_continual$sample_ID) %>% 
  widen()

comm_mat_control_epi <- comm_df %>% 
  filter(new_group == "epilithic.sessile.invert") %>% 
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

# permanova
comp_1yr_meta <- comm_meta %>% 
  drop_na(comp_1yr)
comp_1yr_sampleID <- comp_1yr_meta %>% 
  pull(sample_ID)
comp_1yr_algae <- comm_mat_algae[comp_1yr_sampleID, ]

algae_pt_perma_1yr <- adonis2(comp_1yr_algae ~ treatment*comp_1yr, 
                              data = comp_1yr_meta,
                              strata = comp_1yr_meta$site)
algae_pt_perma_1yr

comp_2yrs_meta <- comm_meta %>% 
  drop_na(comp_2yrs)
comp_2yrs_sampleID <- comp_2yrs_meta %>% 
  pull(sample_ID)
comp_2yrs_algae <- comm_mat_algae[comp_2yrs_sampleID, ]
algae_pt_perma_2yrs <- adonis2(comp_2yrs_algae ~ treatment*comp_2yrs, 
                               data = comp_2yrs_meta,
                               strata = comp_2yrs_meta$site)
algae_pt_perma_2yrs

comp_3yrs_meta <- comm_meta %>% 
  drop_na(comp_3yrs) %>% 
  unite("combo", comp_3yrs, treatment, sep = "-", remove = FALSE)
comp_3yrs_sampleID <- comp_3yrs_meta %>% 
  pull(sample_ID)
comp_3yrs_algae <- comm_mat_algae[comp_3yrs_sampleID, ]
algae_pt_perma_3yrs <- adonis2(comp_3yrs_algae ~ treatment*comp_3yrs, 
                               data = comp_3yrs_meta,
                               strata = comp_3yrs_meta$site)
algae_pt_perma_3yrs 

# beta dispersion
algae_pt_dist <- vegdist(comp_3yrs_algae, "bray")
algae_betadisper <- betadisper(algae_pt_dist, comp_3yrs_meta$combo)
permutest(algae_betadisper, pairwise = TRUE)
anova(algae_betadisper)
# not significantly different dispersions

# ⟞ ⟞ ⟞ 2. SIMPER ---------------------------------------------------------

# three year comparison
comp_3yrs_sampleID_continual <- comp_3yrs_meta %>% 
  filter(treatment == "continual") %>% 
  pull(sample_ID)

# comparing continual removal plots across time periods
simper_algae_3yrs <- simper(
  # subset algae matrix to only include continual removal plots
  comm = comp_3yrs_algae[comp_3yrs_sampleID_continual, ], 
  # only 2 year comparisons
  group = comp_3yrs_meta %>% 
    filter(treatment == "continual") %>% 
    pull(comp_3yrs)
)

algae_SA_comp3yrs <- simper_algae_3yrs$start_after %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average) %>% 
  head(10) %>% 
  mutate(comparison = "start-after")
# PTCA, CYOS, CO, CC, EC, DL, R, POLA, EGME, RAT
algae_SD_comp3yrs <- simper_algae_3yrs$start_during %>% 
  as_tibble() %>% 
  arrange(-average) %>% 
  head(10) %>% 
  mutate(comparison = "start-during")
# CYOS, PTCA, CC, EGME, DL, R, EC, SAMU, RAT, POLA
algae_DA_comp3yrs <- simper_algae_3yrs$during_after %>% 
  as_tibble() %>% 
  arrange(-average) %>% 
  head(10) %>% 
  mutate(comparison = "during-after")
# PTCA, CYOS, CO, EC, EGME, CC, R, SAMU, DL, RAT

algae_comp3yrs <- rbind(algae_SA_comp3yrs, algae_SD_comp3yrs, algae_DA_comp3yrs)

# comparing treatment plots across time periods
simper_algae_treatment <- simper(
  comm = comp_3yrs_algae, 
  group = comp_3yrs_meta$treatment
)

algae_simper_treatment <- simper_algae_treatment$control_continual %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average) %>% 
  head(10)
# PTCA, CYOS, DL, CC, EC, CO, R, POLA, RAT, EGME

# ⟞ ⟞ ⟞ 3. plotting -------------------------------------------------------

# points into data frame for plotting
algae_pt_bray_plotdf <- scores(algae_pt_bray, display = "sites") %>% 
  as_tibble(rownames = "sample_ID") %>% 
  # join with metadata
  left_join(., comm_meta, by = "sample_ID") %>% 
  # taking out outlier for visualization
  filter(sample_ID != "mohk_control_2020-08-12")

# pull top species from simper analysis
simper_algae_spp <- algae_comp3yrs %>% 
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
  scale_x_continuous(limits = c(-1.45, 1.6)) +
  scale_y_continuous(limits = c(-1.625, 1.425)) +
  labs(shape = "Time period",
       color = "Time period", 
       fill = "Time period",
       title = "(a) Removal") +
  theme(legend.position = "none", # c(0.2, 0.8), 
        # panel.grid = element_blank(),
        legend.background = element_blank())
algae_pt_bray_continual_plot 

algae_pt_bray_continual_plot_arrows <- algae_pt_bray_continual_plot +
  geom_text_repel(data = algae_pt_bray_species,
                  aes(x = NMDS1, y = NMDS2,
                      label = stringr::str_wrap(scientific_name, 4, width = 40)),
                  color = "#C70000", lineheight = 0.8, max.overlaps = 100, size = 1.5) +
  geom_segment(data = algae_pt_bray_species,
               aes(x = 0, y = 0,
                   xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "#C70000", linewidth = 0.5)

algae_pt_bray_continual_plot_arrows


# control plots only
algae_pt_bray_control_plot <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "control", algae_pt_bray_species
) +
  # axis limits
  scale_x_continuous(limits = c(-1.45, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.55)) +
  # labels
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "(b) Reference") +
  # theme
  theme(legend.position = "none",
        aspect.ratio = 1
        #panel.grid = element_blank()
        ) 
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

# permanova
comp_1yr_meta <- comm_meta %>% 
  drop_na(comp_1yr)
comp_1yr_sampleID <- comp_1yr_meta %>% 
  pull(sample_ID)
comp_1yr_epi <- comm_mat_epi[comp_1yr_sampleID, ]
epi_pt_perma_1yr <- adonis2(comp_1yr_epi ~ treatment*comp_1yr, 
                            data = comp_1yr_meta,
                            strata = comp_1yr_meta$site)
epi_pt_perma_1yr

comp_2yrs_meta <- comm_meta %>% 
  drop_na(comp_2yrs)
comp_2yrs_sampleID <- comp_2yrs_meta %>% 
  pull(sample_ID)
comp_2yrs_epi <- comm_mat_epi[comp_2yrs_sampleID, ]
epi_pt_perma_2yrs <- adonis2(comp_2yrs_epi ~ treatment*comp_2yrs, 
                             data = comp_2yrs_meta,
                             strata = comp_2yrs_meta$site)
epi_pt_perma_2yrs

comp_3yrs_meta <- comm_meta %>% 
  drop_na(comp_3yrs) %>% 
  unite("combo", comp_3yrs, treatment, sep = "-", remove = FALSE)
comp_3yrs_sampleID <- comp_3yrs_meta %>% 
  pull(sample_ID)
comp_3yrs_epi <- comm_mat_epi[comp_3yrs_sampleID, ]
epi_pt_perma_3yrs <- adonis2(comm_mat_epi ~ treatment*comp_3yrs,
                             data = comp_3yrs_meta,
                             strata = comp_3yrs_meta$site)
epi_pt_perma_3yrs # same as 2 years

# beta dispersion
epi_pt_dist <- vegdist(comm_mat_epi, "bray")
epi_betadisper <- betadisper(epi_pt_dist, comp_3yrs_meta$combo)
anova(epi_betadisper)
permutest(epi_betadisper, pairwise = TRUE)
TukeyHSD(epi_betadisper, which = "group", ordered = FALSE,
         conf.level = 0.95)
# significantly different dispersions

# ⟞ ⟞ ⟞ 2. SIMPER ---------------------------------------------------------

# comparing continual removal plots across time periods
simper_epi_3yrs <- simper(
  # subset epi matrix to only include continual removal plots
  comm = comp_3yrs_epi[comp_3yrs_sampleID_continual, ], 
  # only 2 year comparisons
  group = comp_3yrs_meta %>% 
    filter(treatment == "continual") %>% 
    pull(comp_3yrs)
)

epi_SA_comp3yrs <- simper_epi_3yrs$start_after %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average) %>% 
  head(10) %>% 
  mutate(comparison = "start-after")
# MUCA, TEAU, ANSP, CRGI, URLO, STMO, DIOR, ES, DC, PRUB
epi_SD_comp3yrs <- simper_epi_3yrs$start_during %>% 
  as_tibble() %>% 
  arrange(-average) %>% 
  head(10) %>% 
  mutate(comparison = "start-during")
# TEAU, MUCA, ANSP, URLO, STMO, BAEL, CRGI, DIOR, BA, PRUB
epi_DA_comp3yrs <- simper_epi_3yrs$during_after %>% 
  as_tibble() %>% 
  arrange(-average) %>% 
  head(10) %>% 
  mutate(comparison = "during-after")
# MUCA, CRGI, TEAU, ANSP, URLO, DIOR, ES, DC, BAEL, TC

epi_comp3yrs <- rbind(epi_SA_comp3yrs, epi_SD_comp3yrs, epi_DA_comp3yrs)

# comparing treatment plots across time periods
simper_epi_treatment <- simper(
  comm = comp_3yrs_epi, 
  group = comp_3yrs_meta$treatment
)

epi_simper_treatment <- simper_epi_treatment$control_continual %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average) %>% 
  head(10)
# MUCA, ANSP, TEAU, CRGI, STMO, URLO, CUSP, DIOR, ES, BAEL

# ⟞ ⟞ ⟞ 3. plotting -------------------------------------------------------

# points into data frame for plotting
epi_pt_bray_plotdf <- scores(epi_pt_bray, display = "sites") %>% 
  as_tibble(rownames = "sample_ID") %>% 
  # join with metadata
  left_join(., comm_meta, by = "sample_ID")

# pull top species from simper analysis
simper_epi_spp <- epi_comp3yrs %>% 
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
  scale_x_continuous(limits = c(-1.6, 1.2)) +
  scale_y_continuous(limits = c(-1.6, 1.2)) +
  theme(legend.position = "none"
        # panel.grid = element_blank()
        ) +
  labs(title = "(c) Removal") 
epi_pt_bray_continual_plot

epi_pt_bray_continual_plot_arrows <- epi_pt_bray_continual_plot +
  # arrows
  geom_text_repel(data = epi_pt_bray_species,
                  aes(x = NMDS1, y = NMDS2,
                      label = stringr::str_wrap(scientific_name, 4, width = 40)),
                  color = "#C70000", lineheight = 0.8, max.overlaps = 100, size = 1.5) +
  geom_segment(data = epi_pt_bray_species,
               aes(x = 0, y = 0,
                   xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "#C70000", linewidth = 0.5) 
epi_pt_bray_continual_plot_arrows

# control plots only
epi_pt_bray_control_plot <- nmds_plot_fxn(
  epi_pt_bray_plotdf, "control", epi_pt_bray_species
) +
  scale_x_continuous(limits = c(-1.3, 1.3), breaks = seq(-1, 1, by = 1)) +
  scale_y_continuous(limits = c(-1.25, 1.35), breaks = seq(-1, 1, by = 1)) +
  labs(title = "(d) Reference") +
  theme(legend.position = "none"
        # panel.grid = element_blank()
        ) 
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
# 3. individual species ---------------------------------------------------
##########################################################################-

# ⟞ a. algae species ------------------------------------------------------

algae_comp3yrs_spp <- algae_comp3yrs %>% 
  select(species) %>% 
  unique() %>% 
  pull()

algae_biomass_time_df <- biomass %>% 
  filter(treatment == "continual") %>% 
  # select columns of interest 
  dplyr::select(site, year, month, treatment, date, new_group, sp_code, dry_gm2) %>% 
  exp_dates_column_continual() %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual() %>% 
  filter(sp_code %in% algae_comp3yrs_spp) %>% 
  drop_na(dry_gm2) %>% 
  left_join(., algae_spp_names, by = "sp_code") %>% 
  mutate(scientific_name = case_match(
    scientific_name, 
    "Chondracanthus corymbiferus; Chondracanthus exasperatus" ~ "Chondracanthus spp.", 
    .default = scientific_name)
  )

algae_biomass_time_plot <- ggplot(data = algae_biomass_time_df,
                                  aes(x = time_since_end, y = dry_gm2)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(col = site), linewidth = 2) +
  geom_point(aes(shape = site, col = site), fill = "#FFFFFF", size = 1) +
  # continual
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full), guide = "none") +
  scale_color_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  facet_wrap(~scientific_name, scales = "free_y") +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        # plot.margin = margin(0, 0, 0, 0),
        # legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Time since end of removal (years)", 
       y = expression(Biomass~"("~dry~g/m^{"2"}~")"), 
       color = "Site")

algae_biomass_time_plot

# ⟞ b. epi invert species -------------------------------------------------

epi_comp3yrs_spp <- epi_comp3yrs %>% 
  select(species) %>% 
  unique() %>% 
  pull()

epi_biomass_time_df <- biomass %>% 
  filter(treatment == "continual") %>% 
  # select columns of interest 
  dplyr::select(site, year, month, treatment, date, new_group, sp_code, dry_gm2) %>% 
  exp_dates_column_continual() %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual() %>% 
  filter(sp_code %in% epi_comp3yrs_spp) %>% 
  drop_na(dry_gm2) %>% 
  left_join(., epi_spp_names, by = "sp_code") %>% 
  mutate(scientific_name = case_match(
    scientific_name, 
    "Unidentified Demospongiae spp." ~ "Demospongiae spp.", 
    .default = scientific_name)
  )

epi_biomass_time_plot <- ggplot(data = epi_biomass_time_df,
                                aes(x = time_since_end, y = dry_gm2)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(col = site), linewidth = 2) +
  geom_point(aes(shape = site, col = site), fill = "#FFFFFF", size = 1) +
  # continual
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full), guide = "none") +
  scale_color_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  facet_wrap(~scientific_name, scales = "free_y") +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        # plot.margin = margin(0, 0, 0, 0),
        # legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Time since end of removal (years)", 
       y = expression(Biomass~"("~dry~g/m^{"2"}~")"), 
       color = "Site")

epi_biomass_time_plot


##########################################################################-
# 4. manuscript tables ----------------------------------------------------
##########################################################################-

anova_1yr_tables <- rbind(anova_summary_fxn(algae_pt_perma_1yr), anova_summary_fxn(epi_pt_perma_1yr)) %>% 
  rename_with(., .fn = ~paste(., "_1yr", sep = "", .cols = everything(cols)))
anova_1yr_tables 

anova_2yrs_tables <- rbind(anova_summary_fxn(algae_pt_perma_2yrs), anova_summary_fxn(epi_pt_perma_2yrs)) %>% 
  rename_with(., .fn = ~paste(., "_2yrs", sep = "", .cols = everything(cols)))
anova_2yrs_tables 

anova_3yrs_tables <- rbind(anova_summary_fxn(algae_pt_perma_3yrs), anova_summary_fxn(epi_pt_perma_3yrs)) %>% 
  rename_with(., .fn = ~paste(., "_3yrs", sep = "", .cols = everything(cols)))
anova_3yrs_tables

anova_together_tables <- cbind(anova_1yr_tables, anova_2yrs_tables, anova_3yrs_tables) %>% 
  # take out unwanted columns
  select(!c("model_1yr1", "model_2yrs1", "model_3yrs1", 
            "SumOfSqs_1yr1", "SumOfSqs_2yrs1", "SumOfSqs_3yrs1",
            "R2_1yr1", "R2_2yrs1", "R2_3yrs1")) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  # group labels
  tab_row_group(
    label = "Epilithic invertebrates", rows = 6:10
  ) %>% 
  tab_row_group(
    label = "Understory macroalgae", rows = 1:5
  ) %>% 
  # 1, 2, and 3 year comparisons
  tab_spanner(
    label = "1 year comparison",
    columns = c(variables_1yr1, Df_1yr1, F_1yr1, p_1yr1)
  ) %>% 
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
    variables_1yr1 = "Source of variation",
    Df_1yr1 = "df",
    F_1yr1 = "pseudo-F",
    p_1yr1 = "p", 
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
      columns = p_1yr1,
      rows = p_1yr1 < 0.05
    )
  ) %>% 
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
#        here::here("tables", "ms-tables", paste("tbl-S4_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)


##########################################################################-
# 5. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. continual removal --------------------------------------------------

comm_comp_together <- algae_pt_bray_continual_plot + epi_pt_bray_continual_plot 

comm_comp_together_arrows <- algae_pt_bray_continual_plot_arrows + epi_pt_bray_continual_plot_arrows

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S8_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together_arrows,
#        height = 8, width = 16, units = "cm",
#        dpi = 400)

# ⟞ b. control ------------------------------------------------------------

comm_comp_control <- algae_pt_bray_control_plot + epi_pt_bray_control_plot 

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S7_", today(), ".jpg", sep = "")),
#        plot = comm_comp_control,
#        height = 6, width = 16, units = "cm",
#        dpi = 300)


# ⟞ c. raw species biomass ------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S9_", today(), ".jpg", sep = "")),
#        plot = algae_biomass_time_plot,
#        height = 18, width = 24, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S10_", today(), ".jpg", sep = "")),
#        plot = epi_biomass_time_plot,
#        height = 18, width = 24, units = "cm",
#        dpi = 300)

# ⟞ d. reference and removal plots together -------------------------------

# putting algae plots together
algae_plots <- plot_grid(algae_pt_bray_continual_plot, algae_pt_bray_control_plot, 
                         nrow = 2)

# putting invert plots together
epi_plots <- plot_grid(epi_pt_bray_continual_plot, epi_pt_bray_control_plot, 
                       nrow = 2)

# group plots with labels
algae_labelled <- plot_grid(algae_title, algae_plots, 
                            rel_heights = c(1, 12), 
                            ncol = 1)


epi_labelled <- plot_grid(epi_title, epi_plots, 
                          rel_heights = c(1, 12), 
                          ncol = 1)

# getting the legend as separate object
plot_legend <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "continual", algae_pt_bray_species
) +
  labs(shape = "Time period",
       color = "Time period", 
       fill = "Time period") 
legend <- cowplot::get_legend(plot_legend)

# putting plots together
all_plots <- plot_grid(algae_labelled, epi_labelled, ncol = 2,
                                rel_widths = c(1, 1))

# putting plots with legend
final_plot <- plot_grid(all_plots, legend, ncol = 2, rel_widths = c(1, 0.3))

# saving
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_", today(), ".jpg", sep = "")),
#        plot = final_plot,
#        height = 13, width = 15, units = "cm",
#        dpi = 300)

