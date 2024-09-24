
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. set up -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. source -------------------------------------------------------------

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

# ⟞ b. functions ----------------------------------------------------------

# This is a function to generate the NMDS plots in section 3a.

nmds_plot_fxn <- function(df, time_period) {
  df %>% 
    filter(comp_3yrs == time_period) %>% 
    ggplot(aes(x = NMDS1, y = NMDS2,
               color = treatment,
               fill = treatment,
               shape = treatment,
               linetype = treatment)) +
    coord_fixed(ratio = 1) +
    geom_vline(xintercept = 0, color = "grey", lty = 2) +
    geom_hline(yintercept = 0, color = "grey", lty = 2) +
    geom_point(size = 1, alpha = 0.9) +
    stat_ellipse() +
    scale_linetype_manual(values = c("continual" = 1, 
                                     "control" = 2)) +
    scale_color_manual(values = c("continual" = removal_col, 
                                  "control" = reference_col)) +
    theme_bw() +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 7),
          legend.text = element_text(size = 7), 
          legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(size = 8),
          plot.title.position = "plot",
          legend.key.size = unit(0.5, units = "cm"),
          aspect.ratio = 1)
}

# This is a function to generate the PERMANOVA summary tables in section 4a.

anova_summary_fxn <- function(adonis2.obj, name) {
 adonis2.obj %>% 
    # turn adonis2 result into data frame
    as.data.frame() %>% 
    # make rownames "variables"
    rownames_to_column("variables") %>% 
    # rename Pr(>F) column into something intelligible
    rename(p = `Pr(>F)`) %>% 
    # round values to 2 decimal points
    mutate(across(SumOfSqs:p, ~ round(.x, digits = 3))) %>% 
    # replace comp_.yrs with time period
    mutate(variables = str_replace(variables, "comp_.yrs", "time period")) %>% 
    mutate(variables = str_replace(variables, "comp_.yr", "time period")) %>% 
    # make object name column
    mutate(model = name) 
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------- 1. data frames and wrangling functions -----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section contains code for creating two data frames:
# `comm_df`: a data frame with all species biomass and metadata
# `comm_meta`: a data frame with just the metadata for each sampling point

# ⟞ a. species biomass with metadata --------------------------------------

comm_df <- biomass %>% 
  filter(treatment %in% c("control", "continual")) %>% 
  # select columns of interest 
  dplyr::select(site, year, month, treatment, 
                date, new_group, sp_code, dry_gm2) %>% 
  unite("sample_ID_short", site, date, remove = FALSE) %>% 
  # filtered from kelp delta data frame created in upstream script
  filter(sample_ID_short %in% (delta_continual$sample_ID_short)) %>% 
  # add column for experiment during/after
  exp_dates_column_continual() %>% 
  # add column for time since end of the experiment
  time_since_columns_continual() %>% 
  # add column for kelp year
  kelp_year_column() %>% 
  # add column for 1 year, 2 years, 3 years comparison
  comparison_column_continual_new() %>% 
  # join with site quality data frame and data frame of full names of site
  full_join(., site_quality, 
            by = "site") %>% 
  left_join(., enframe(sites_full), 
            by = c("site" = "name")) %>% 
  rename(site_full = value) %>% 
  mutate(site_full = fct_relevel(
    site_full, 
    "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>%
  # create new sample ID with treatment
  unite("sample_ID", site, treatment, date, remove = FALSE) %>% 
  # only include 3 year sampling sites
  drop_na(comp_3yrs)

# ⟞ b. metadata only ------------------------------------------------------

# metadata for all plots
comm_meta <- comm_df %>% 
  select(sample_ID, site, date, year, month, 
         treatment, exp_dates, quarter, time_yrs, 
         time_since_start, time_since_end, kelp_year, 
         comp_1yrs, comp_2yrs, comp_3yrs, quality, site_full) %>% 
  unique()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 2. community analyses -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section uses nested data frames to clean and analyse the data in 
# `comm_df` created in the previous section. 

# ⟞ a. community matrix and metadata --------------------------------------

# This section creates a nested data frame called `comm_df_nested`. This data
# frame includes nested data frames for understory algae and sessile inverts.

# The column `comm_mat_filtered` contains the community matrix. Rows are 
# sampling points and columns are species, with each cell representing biomass. 
# The community matrices are filtered to only include surveys where there were 
# any species detected (300/300 surveys) and only include species that were 
# detected at least once (thus excludes HPAC, MT, SC for sessile invertebrates, 
# and EA, FR, IR, NEO, SELO for understory algae).

# The column `comm_meta_filtered` contains the metadata for each
# sampling point.

# making a nested data frame
comm_df_nested <- comm_df %>% 
  # take out giant kelp
  filter(sp_code != "MAPY") %>% 
  # create nested columns
  nest(data = everything(), .by = new_group) %>% 
  # only include understory algae and sessile inverts
  filter(new_group %in% c("algae", "epilithic.sessile.invert")) %>% 
  # get data into community matrix format
  mutate(comm_mat = map(
    data, 
    ~ # select columns of interest
      select(.x, sample_ID, sp_code, dry_gm2) %>% 
      # get into wide format for community analysis
      pivot_wider(names_from = sp_code, 
                  values_from = dry_gm2) %>% 
      # make the sample_ID column row names
      column_to_rownames("sample_ID") %>% 
      replace(is.na(.), 0))) %>% 
  # get metadata
  mutate(comm_meta = map(
    data, 
    ~ select(.x, sample_ID, site, date, year, month, 
             treatment, exp_dates, quarter, time_yrs, 
             time_since_start, time_since_end, kelp_year, 
             comp_1yrs, comp_2yrs, comp_3yrs,
             quality, site_full) %>% 
      unique()
  )) %>% 
  # find species that do occur (exclude species that never occur)
  mutate(keepspp = map(
    comm_mat, 
    ~ colSums(.x, na.rm = FALSE) %>% 
      enframe() %>% 
      filter(value > 0) %>% 
      pull(name))) %>% 
  # find surveys where there was actually something there
  mutate(keepsurveys = map(
    comm_mat,
    ~ rowSums(.x, na.rm = FALSE) %>% 
      enframe() %>% 
      filter(value > 0) %>% 
      pull(name))) %>% 
  # only retain species that do occur
  mutate(comm_mat_step1 = map2(
    comm_mat, keepspp,
    ~ select(.x, any_of(.y)))) %>% 
  # only retain surveys where there was something there
  mutate(comm_mat_filtered = map2(
    comm_mat_step1, keepsurveys,
    ~ rownames_to_column(.x, "surveys") %>% 
      filter(surveys %in% .y) %>% 
      column_to_rownames("surveys"))) %>% 
  # filter metadata
  mutate(comm_meta_filtered = map2(
    comm_meta, comm_mat_filtered, 
    ~ filter(.x, sample_ID %in% rownames(.y))))

# ⟞ b. non-metric multidimensional scaling (NMDS) -------------------------

# This section generates `comm_nmds`, another nested data frame that includes
# ordinations using Bray-Curtis dissimilarity and modified Gower dissimilarity
# for understory algae and sessile inverts. This data frame also includes
# stress plots to visualize ordination stress and biplots as a preview of what
# the ordination should look like.

set.seed(1)

comm_nmds <- comm_df_nested %>% 
  select(new_group, comm_mat_filtered, comm_meta_filtered) %>% 
  # bray-curtis
  mutate(nmds_bray = map(
    comm_mat_filtered, 
    ~ metaMDS(.x, "bray"))) %>% 
  # modified Gower
  mutate(nmds_altgower = map(
    comm_mat_filtered, 
    ~ metaMDS(.x, "altGower"))) %>% 
  # bray stress plot
  mutate(bray_stressplot = map(
    nmds_bray,
    ~ stressplot(.x))) %>% 
  # bray biplot
  mutate(bray_plot = map(
    nmds_bray,
    ~ plot(.x))) %>% 
  # modified Gower stress plot
  mutate(altgower_stressplot = map(
    nmds_altgower,
    ~ stressplot(.x))) %>% 
  # modified Gower biplot
  mutate(altgower_plot = map(
    nmds_altgower,
    ~ plot(.x)))

# ⟞ c. analyses -----------------------------------------------------------

# This section creates `comm_permanova`, a nested data frame that includes
# subsets of the community matrices and metadata data frames corresponding to 
# the 1 year, 2 year, and 3 year comparisons. There are columns displaying the
# output of the permutational analysis of variance (PERMANOVA) for each 
# comparison with treatment (reference or removal) and time period (start,
# during, and after) as predictors, with site as grouping variable. As in
# section b, these analyses are repeated using Bray-Curtis dissimilarity and
# modified Gower dissimilarity.

set.seed(1)

comm_permanova <- comm_nmds %>% 
  
  # filter metadata data frame to only include comparison years of interest
  mutate(comp_1yr_meta = map(
    comm_meta_filtered,
    ~ drop_na(.x, comp_1yrs))) %>% 
  # pull sample_ID column to get names of surveys to keep
  mutate(comp_1yr_sampleID = map(
    comp_1yr_meta,
    ~ pull(.x, sample_ID))) %>% 
  # filter community matrix to only include surveys of interest
  mutate(comp_1yr_mat = map2(
    comm_mat_filtered, comp_1yr_sampleID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
  # same for the two year comparison
  mutate(comp_2yrs_meta = map(
    comm_meta_filtered,
    ~ drop_na(.x, comp_2yrs))) %>% 
  mutate(comp_2yrs_sampleID = map(
    comp_2yrs_meta,
    ~ pull(.x, sample_ID))) %>% 
  mutate(comp_2yrs_mat = map2(
    comm_mat_filtered, comp_2yrs_sampleID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
  # same for the 3 year comparison
  mutate(comp_3yrs_meta = map(
    comm_meta_filtered,
    ~ drop_na(.x, comp_3yrs))) %>% 
  mutate(comp_3yrs_sampleID = map(
    comp_3yrs_meta,
    ~ pull(.x, sample_ID))) %>% 
  mutate(comp_3yrs_mat = map2(
    comm_mat_filtered, comp_3yrs_sampleID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
  
  # 1 year comparison with bray-curtis
  mutate(bray_comp_1yr_permanova = map2(
    comp_1yr_mat, comp_1yr_meta,
    ~ adonis2(.x ~ treatment*comp_1yrs,
              data = .y,
              strata = .y$site,
              method = "bray"))) %>% 
  # 2 year comparison with bray-curtis
  mutate(bray_comp_2yrs_permanova = map2(
    comp_2yrs_mat, comp_2yrs_meta,
    ~ adonis2(.x ~ treatment*comp_2yrs,
              data = .y,
              strata = .y$site,
              method = "bray"))) %>% 
  # 3 year comparison with bray-curtis
  mutate(bray_comp_3yrs_permanova = map2(
    comp_3yrs_mat, comp_3yrs_meta,
    ~ adonis2(.x ~ treatment*comp_3yrs,
              data = .y,
              strata = .y$site,
              method = "bray"))) %>% 
  # 1 year comparison with modified Gower
  mutate(altgower_comp_1yr_permanova = map2(
    comp_1yr_mat, comp_1yr_meta,
    ~ adonis2(.x ~ treatment*comp_1yrs,
              data = .y,
              strata = .y$site,
              method = "altGower"))) %>% 
  # 2 year comparison with modified Gower
  mutate(altgower_comp_2yrs_permanova = map2(
    comp_2yrs_mat, comp_2yrs_meta,
    ~ adonis2(.x ~ treatment*comp_2yrs,
              data = .y,
              strata = .y$site,
              method = "altGower"))) %>% 
  # 3 year comparison with modified Gower
  mutate(altgower_comp_3yrs_permanova = map2(
    comp_3yrs_mat, comp_3yrs_meta,
    ~ adonis2(.x ~ treatment*comp_3yrs,
              data = .y,
              strata = .y$site,
              method = "altGower")))

algae_permanova_results <- pluck(comm_permanova, 24, 2)

epi_permanova_results <- pluck(comm_permanova, 24, 1)

# PERMANOVA results indicated that there was no interaction between treatment 
# and time period for understory algae though there was a significant 
# interaction between these two terms for sessile inverts.

# Thus, in this next section, I compared beta dispersion for treatment and time
# period independently for understory algae, and in combination for sessile
# invertebrates, following guidance by Marti Anderson in this R-sig-eco forum
# post from 2010-09-02: 
# https://stat.ethz.ch/pipermail/r-sig-ecology/2010-September/001524.html

set.seed(1)

comm_betadisper <- comm_df_nested %>% 
  select(new_group, comm_mat_filtered, comm_meta_filtered) %>% 
  mutate(comm_meta_filtered = map(
    comm_meta_filtered,
    ~ unite(., "combo", comp_3yrs, treatment, sep = "-", remove = FALSE) 
  )) %>% 
  mutate(dist_bray = map(
    comm_mat_filtered,
    ~ vegdist(., method = "bray")
  )) %>% 
  mutate(dist_altgower = map(
    comm_mat_filtered,
    ~ vegdist(., method = "altGower")
  )) 

set.seed(1)

# calculating beta dispersion between treatments for algae
algae_betadisper_treatment <- betadisper(
  comm_betadisper[[5]][[2]], 
  comm_betadisper[[3]][[2]]$treatment) 

permutest(algae_betadisper_treatment)
# difference in dispersion between reference and removal (so differences
# in community composition could be due to changes in location and dispersion)

set.seed(1)

# calculating beta dispersion between time periods for algae
algae_betadisper_time <- betadisper(
  comm_betadisper[[5]][[2]], 
  comm_betadisper[[3]][[2]]$comp_3yrs) 

permutest(algae_betadisper_time)
# no difference in dispersion through time (so differences in community 
# composition are due to changes in location, not dispersion)

set.seed(1)

# calculating beta dispersion between treatments and time periods for inverts
epi_betadisper_combo <- betadisper(
  comm_betadisper[[5]][[1]], 
  comm_betadisper[[3]][[1]]$combo)

permutest(epi_betadisper_combo)
# there is a difference in dispersions

TukeyHSD(epi_betadisper_combo)
# pairwise comparisons of differences between dispersions that are meaningful:
# during control and after control: reference plots differed in dispersion
# in the last 3 years of the experiment and in the recovery period
# start control and start continual: reference and removal plots differed in
# dispersion at the start of the experiment

# ⟞ d. control pairwise comparisons ---------------------------------------

# This section compares composition in reference plots between time periods:
# the "start" of the removal, "during" the removal, and "after" the removal 
# (in the kelp recovery period). The purpose is to show that reference plot
# community composition has also changed through time.

# The `control_pairwise` data frame is a nested data frame that includes 
# PERMANOVA output that is put into summary tables to compare start-during, 
# start-after, and during-after compositions for understory algae and sessile
# inverts. The code to create the tables is in section 4b.

set.seed(1)

control_pairwise <- comm_permanova %>% 
  select(new_group, comp_3yrs_meta, comp_3yrs_mat) %>% 
  # create subsetted metadata data frames 
  # to only include time periods being compared
  mutate(control_start_during_meta = map(
    comp_3yrs_meta,
    ~ filter(.x, 
             comp_3yrs %in% c("start", "during") & treatment == "control"))) %>% 
  mutate(control_during_after_meta = map(
    comp_3yrs_meta,
    ~ filter(.x, 
             comp_3yrs %in% c("during", "after") & treatment == "control"))) %>% 
  mutate(control_start_after_meta = map(
    comp_3yrs_meta,
    ~ filter(.x, 
             comp_3yrs %in% c("start", "after") & treatment == "control"))) %>% 
  
  # extract sample_IDs corresponding to surveys to compare
  mutate(control_start_during_sample_ID = map(
    control_start_during_meta,
    ~ pull(.x, sample_ID)
  )) %>% 
  mutate(control_during_after_sample_ID = map(
    control_during_after_meta,
    ~ pull(.x, sample_ID)
  )) %>% 
  mutate(control_start_after_sample_ID = map(
    control_start_after_meta,
    ~ pull(.x, sample_ID)
  )) %>% 
  
  # create subsetted community matrices that only include surveys to compare
  mutate(control_start_during_mat = map2(
    comp_3yrs_mat, control_start_during_sample_ID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
  mutate(control_during_after_mat = map2(
    comp_3yrs_mat, control_during_after_sample_ID,
    ~ filter(.x, row.names(.x) %in% .y))) %>%   
  mutate(control_start_after_mat = map2(
    comp_3yrs_mat, control_start_after_sample_ID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
  
  # compare composition using PERMANOVA 
  mutate(start_during_altgower = map2(
    control_start_during_mat, control_start_during_meta,
    ~ adonis2(.x ~ comp_3yrs,
              data = .y,
              strata = .y$site,
              method = "altGower") 
  )) %>% 
  mutate(during_after_altgower = map2(
    control_during_after_mat, control_during_after_meta,
    ~ adonis2(.x ~ comp_3yrs,
              data = .y,
              strata = .y$site,
              method = "altGower") 
  )) %>% 
  mutate(start_after_altgower = map2(
    control_start_after_mat, control_start_after_meta,
    ~ adonis2(.x ~ comp_3yrs,
              data = .y,
              strata = .y$site,
              method = "altGower") 
  )) %>% 

  # generate tables from PERMANOVA output
  mutate(start_during_altgower_table = map(
    start_during_altgower,
    ~ as_tibble(.x) %>% 
      mutate(Components = c("Time period", "Residual", "Total")) %>% 
      relocate(Components, .before = Df) %>% 
      mutate(comparison = "start-during") %>% 
      mutate(across(SumOfSqs:`F`, ~ round(., digits = 2)))
  )) %>% 
  mutate(start_after_altgower_table = map(
    start_after_altgower,
    ~ as_tibble(.x) %>% 
      mutate(Components = c("Time period", "Residual", "Total")) %>% 
      relocate(Components, .before = Df) %>% 
      mutate(comparison = "start-after") %>% 
      mutate(across(SumOfSqs:`F`, ~ round(., digits = 2)))
  )) %>% 
  mutate(during_after_altgower_table = map(
    during_after_altgower,
    ~ as_tibble(.x) %>% 
      mutate(Components = c("Time period", "Residual", "Total")) %>% 
      relocate(Components, .before = Df) %>% 
      mutate(comparison = "during-after") %>% 
      mutate(across(SumOfSqs:`F`, ~ round(., digits = 2)))
  )) %>% 
  
  # combine PERMANOVA tables for downstream table making
  mutate(full_table = pmap(
    list(start_during_altgower_table,
         start_after_altgower_table,
         during_after_altgower_table),
    bind_cols
  )) %>% 
  mutate(full_table = map2(
    full_table, new_group,
    ~ mutate(.x, group = .y) %>% 
      relocate(group, .before = "Components...1")
  ))


# ⟞ e. SIMPER analysis ----------------------------------------------------

# This section of code calculates similarity percentages (SIMPER) to show what
# species contribute to differences between time periods within the removal
# plots only. SIMPER is based on Bray-Curtis dissimilarity, not modified Gower; 
# thus, we did this analysis to determine which species to highlight in
# individual supplementary figures rather than to determine exactly how each 
# species contributes to differences between communities at each time period.

comm_simper <- comm_df_nested %>%
  select(new_group, comm_mat_filtered, comm_meta_filtered) %>% 
  # only include surveys in the continual treatment
  mutate(removal_3yrs_meta = map(
    comm_meta_filtered,
    ~ filter(.x, treatment == "continual") %>% 
      drop_na(comp_3yrs)
  )) %>% 
  mutate(removal_3yrs_IDs = map(
    removal_3yrs_meta,
    ~ pull(.x, sample_ID)
  )) %>% 
  mutate(removal_3yrs_mat = map2(
    comm_mat_filtered, removal_3yrs_IDs,
    ~ filter(.x, rownames(.x) %in% .y)
  )) %>% 
  
  # SIMPER analysis comparing time periods
  mutate(simper_time = map2(
    removal_3yrs_mat, removal_3yrs_meta,
    ~ simper(comm = .x, 
             group = .y$comp_3yrs)
  )) %>% 
  
  # extract the top 10 species contributing to differences between time periods
  # repeating this step for all three pairwise comparisons
  mutate(simper_start_during_spp = map(
    simper_time,
    ~ .x$start_during$cusum %>% 
      enframe() %>% 
      head(10) %>% 
      pull(name)
  )) %>% 
  mutate(simper_start_after_spp = map(
    simper_time,
    ~ .x$start_after$cusum %>% 
      enframe() %>% 
      head(10) %>% 
      pull(name)
  )) %>% 
  mutate(simper_during_after_spp = map(
    simper_time,
    ~ .x$during_after$cusum %>% 
      enframe() %>% 
      head(10) %>% 
      pull(name)
  )) %>% 
  
  # generating a list of all species codes
  mutate(simper_list = pmap(
    list(simper_start_during_spp, 
         simper_start_after_spp, 
         simper_during_after_spp),
    c 
  )) %>% 
  
  # getting rid of duplicate species codes
  mutate(simper_no_dupes = map(
    simper_list,
    ~ unique(.x)
  ))

# extracting species codes for downstream visualizations

algae_simper_spp <- pluck(comm_simper, 12, 2)

epi_simper_spp <- pluck(comm_simper, 12, 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 3. visualizations ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. ordinations --------------------------------------------------------

# This section includes code to create the NMDS plots in figure 3a-d.

comm_visualizations <- comm_nmds %>% 
  select(new_group, comm_meta_filtered, nmds_altgower) %>% 
  # extracting scores (NMDS1 and NMDS2 coordinates), joining with metadata
  mutate(nmds_scores = map2(
    nmds_altgower, comm_meta_filtered,
    ~ scores(.x, display = "sites") %>% 
      as_tibble(rownames = "sample_ID") %>% 
      left_join(., .y, by = "sample_ID")
  ))

# extracting scores for downstream visualizations

algae_scores <- pluck(comm_visualizations, 4, 2)

epi_scores <- pluck(comm_visualizations, 4, 1)

# plotting start of experiment reference and removal plots

# making axis limits consistent for easy comparison between panels
algae_axis_limits <- list(
  scale_x_continuous(limits = c(-0.55, 0.55), 
                     breaks = c(-0.5, 0, 0.5)),
  scale_y_continuous(limits = c(-0.55, 0.55), 
                     breaks = c(-0.5, 0, 0.5))
)

algae_start_plot <- nmds_plot_fxn(algae_scores,
                                     "start") +
  algae_axis_limits +
  annotate("text", 
           x = -0.35, y = 0.5, 
           label = "Reference", 
           color = reference_col, 
           size = 2) +
  annotate("text", 
           x = 0.4, y = 0.1, 
           label = "Removal", 
           color = removal_col, 
           size = 2) +
  labs(title = "(a) Start of removal") 

algae_during_plot <- nmds_plot_fxn(algae_scores,
                                     "during") +
  algae_axis_limits +
  labs(title = "(b) End of removal") 

algae_after_plot <- nmds_plot_fxn(algae_scores,
                                      "after") +
  algae_axis_limits +
  labs(title = "(c) Recovery period") 

epi_axis_limits <- list(
  scale_x_continuous(limits = c(-0.35, 0.25)),
  scale_y_continuous(limits = c(-0.3, 0.3))
)

epi_start_plot <- nmds_plot_fxn(epi_scores,
                                     "start") +
  epi_axis_limits +
  annotate("text", 
           x = 0.15, y = 0.25, 
           label = "Reference", 
           color = reference_col, 
           size = 2) +
  annotate("text", 
           x = -0.2, y = -0.2, 
           label = "Removal", 
           color = removal_col, 
           size = 2) +
  labs(title = "(d) Start of removal") 

epi_during_plot <- nmds_plot_fxn(epi_scores,
                                      "during") +
  epi_axis_limits +
  labs(title = "(e) End of removal") 

epi_after_plot <- nmds_plot_fxn(epi_scores,
                                     "after") +
  epi_axis_limits +
  labs(title = "(f) Recovery period") 

# ⟞ b. individual species through time ------------------------------------

# This section includes code to create Supplemental Material figures ______.
# The figures rely on the `biomass` data frame, which is generated in the 
# 00-set_up.R script.

spp_biomass_df <- biomass %>% 
  filter(treatment == "continual") %>% 
  select(site, year, month, treatment, date, 
         new_group, sp_code, scientific_name, dry_gm2) %>% 
  exp_dates_column_continual() %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual_new() %>% 
  # only include species that came out of SIMPER analysis
  filter(sp_code %in% c(algae_simper_spp, epi_simper_spp)) %>% 
  # recode scientific names for clarity
  mutate(scientific_name = case_match(
    scientific_name, 
    "Chondracanthus corymbiferus; Chondracanthus exasperatus" ~ "Chondracanthus spp.", 
    "Unidentified Demospongiae spp." ~ "Demospongiae spp.", 
    .default = scientific_name)
  ) %>% 
  drop_na(dry_gm2)

algae_biomass_time_plot <- spp_biomass_df %>% 
  filter(new_group == "algae") %>% 
  ggplot(aes(x = time_since_end, 
             y = dry_gm2)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(col = site), linewidth = 2) +
  geom_point(aes(shape = site, 
                 col = site), fill = "#FFFFFF", size = 1) +
  scale_shape_manual(values = shape_palette_site, 
                     labels = c("aque" = aque_full, 
                                "napl" = napl_full, 
                                "mohk" = mohk_full, 
                                carp = carp_full), 
                     guide = "none") +
  scale_color_manual(values = color_palette_site, 
                     labels = c("aque" = aque_full, 
                                "napl" = napl_full, 
                                "mohk" = mohk_full, 
                                "carp" = carp_full)) +
  facet_wrap(~scientific_name, scales = "free_y") +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Time since end of removal (years)", 
       y = "Biomass (dry g/m\U00B2)",
       color = "Site")

epi_biomass_time_plot <- spp_biomass_df %>% 
  filter(new_group == "epilithic.sessile.invert") %>% 
  ggplot(aes(x = time_since_end, 
             y = dry_gm2)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(col = site), linewidth = 2) +
  geom_point(aes(shape = site, col = site), fill = "#FFFFFF", size = 1) +
  # continual
  scale_shape_manual(values = shape_palette_site, 
                     labels = c("aque" = aque_full, 
                                "napl" = napl_full, 
                                "mohk" = mohk_full, 
                                carp = carp_full), 
                     guide = "none") +
  scale_color_manual(values = color_palette_site, 
                     labels = c("aque" = aque_full, 
                                "napl" = napl_full, 
                                "mohk" = mohk_full, 
                                carp = carp_full)) +
  facet_wrap(~scientific_name, scales = "free_y") +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Time since end of removal (years)", 
       y = "Biomass (dry g/m\U00B2)",
       color = "Site")

# ⟞ c. saving outputs -----------------------------------------------------

# ⟞ ⟞ i. ordinations ------------------------------------------------------

# This section of code includes plot arranging using functions from `cowplot`.

# putting algae plots together
algae_altgower_plots <- plot_grid(algae_start_plot, 
                                  algae_during_plot, 
                                  algae_after_plot,
                                  ncol = 1)

epi_altgower_plots <- plot_grid(epi_start_plot, 
                                epi_during_plot, 
                                epi_after_plot,
                                ncol = 1)

# group plots with labels
algae_altgower_labelled <- plot_grid(algae_title, algae_altgower_plots, 
                                     rel_heights = c(1, 9), 
                                     ncol = 1)


epi_altgower_labelled <- plot_grid(epi_title, epi_altgower_plots, 
                                   rel_heights = c(1, 9), 
                                   ncol = 1)

# putting plots together
all_altgower_plots <- plot_grid(algae_altgower_labelled, 
                                epi_altgower_labelled, 
                                ncol = 2,
                                rel_widths = c(1, 1))

# saving
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_altgower_v2_", today(), ".jpg", sep = "")),
#        plot = all_altgower_plots,
#        height = 16, width = 12, units = "cm",
#        dpi = 300)


# ⟞ ⟞ ii. individual biomass timeseries -----------------------------------

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 4. tables -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create summary tables for PERMANOVA outputs.

# ⟞ a. PERMANOVA summary tables -------------------------------------------

# This section creates Table _____ and Table _______.

comm_permanova_tables <- comm_permanova %>% 
  select(new_group, 
         altgower_comp_1yr_permanova, 
         altgower_comp_2yrs_permanova, 
         altgower_comp_3yrs_permanova) %>% 
  mutate(perm1yr_table = map(
    altgower_comp_1yr_permanova, 
    ~ anova_summary_fxn(.x, "1yr")
  )) %>% 
  mutate(perm2yr_table = map(
    altgower_comp_2yrs_permanova, 
    ~ anova_summary_fxn(.x, "2yr")
  )) %>% 
  mutate(perm3yr_table = map(
    altgower_comp_3yrs_permanova, 
    ~ anova_summary_fxn(.x, "3yr")
  )) %>% 
  mutate(full_table = pmap(
    list(perm1yr_table, perm2yr_table, perm3yr_table),
    bind_cols
  )) %>% 
  mutate(full_table = map2(
    full_table, new_group,
    ~ mutate(.x, group = .y) %>% 
      relocate(group, .before = "variables...1")
  ))

algae_permanova <- pluck(comm_permanova_tables, 8, 2)

epi_permanova <- pluck(comm_permanova_tables, 8, 1)

# create table using gt
all_permanova <- bind_rows(algae_permanova, epi_permanova) %>% 
  mutate(group = case_when(
    group == "algae" ~ "Understory macroalgae",
    group == "epilithic.sessile.invert" ~ "Sessile invertebrates"
  )) %>% 
  select(!contains("model")) %>% 
  gt(
    groupname_col = "group",
    row_group_as_column = TRUE
  ) %>% 
  sub_missing(
    columns = everything(),
    missing_text = "--"
  ) %>% 
  tab_spanner(label = "1 year",
              columns = c("variables...2", 
                          "Df...3", 
                          "SumOfSqs...4", 
                          "R2...5", 
                          "F...6", 
                          "p...7")) %>% 
  tab_spanner(label = "2 year",
              columns = c("variables...9", 
                          "Df...10", 
                          "SumOfSqs...11", 
                          "R2...12", 
                          "F...13", 
                          "p...14")) %>% 
  tab_spanner(label = "3 year",
              columns = c("variables...16", 
                          "Df...17",
                          "SumOfSqs...18",
                          "R2...19",
                          "F...20",
                          "p...21")) %>% 
  # change column names
  cols_label(
    `variables...2` = "Variables",
    `variables...9` = "Variables",
    `variables...16` = "Variables",
    `Df...3` = "df",
    `Df...10` = "df",
    `Df...17` = "df",
    `SumOfSqs...4` = "Sum of squares",
    `SumOfSqs...11` = "Sum of squares",
    `SumOfSqs...18` = "Sum of squares",
    `R2...5` = "R\U00B2",
    `R2...12` = "R\U00B2",
    `R2...19` = "R\U00B2",
    `F...6` = "F",
    `F...13` = "F",
    `F...20` = "F",
    `p...7` = "Pr(>F)",
    `p...14` = "Pr(>F)",
    `p...21` = "Pr(>F)"
  ) 

permanova_12_table <- all_permanova %>% 
  cols_hide(columns = c("variables...16", 
                        "Df...17",
                        "SumOfSqs...18",
                        "R2...19",
                        "F...20",
                        "p...21"))

permanova_3_table <- all_permanova %>% 
  cols_hide(!(columns = c("variables...16", 
                          "Df...17",
                          "SumOfSqs...18",
                          "R2...19",
                          "F...20",
                          "p...21")))


# ⟞ b. control pairwise tables --------------------------------------------

# extract algae tables
algae_pairwise <- control_pairwise %>% 
  pluck(19, 2)

# extract sessile invert tables
epi_pairwise <- control_pairwise %>% 
  pluck(19, 1)

# create table using gt
all_pairwise <- bind_rows(algae_pairwise, epi_pairwise) %>% 
  mutate(group = case_when(
    group == "algae" ~ "Understory algae",
    group == "epilithic.sessile.invert" ~ "Sessile invertebrates"
  )) %>% 
  select(!contains("comparison")) %>% 
  gt(
    groupname_col = "group",
    row_group_as_column = TRUE
  ) %>% 
  sub_missing(
    columns = everything(),
    missing_text = "--"
  ) %>% 
  tab_spanner(label = "start-during",
              columns = c("Components...2", 
                          "Df...3", 
                          "SumOfSqs...4", 
                          "R2...5", 
                          "F...6", 
                          "Pr(>F)...7")) %>% 
  tab_spanner(label = "start-after",
              columns = c("Components...9", 
                          "Df...10", 
                          "SumOfSqs...11", 
                          "R2...12", 
                          "F...13", 
                          "Pr(>F)...14")) %>% 
  tab_spanner(label = "during-after",
              columns = c("Components...16", 
                          "Df...17",
                          "SumOfSqs...18",
                          "R2...19",
                          "F...20",
                          "Pr(>F)...21")) %>% 
  # change column names
  cols_label(
    `Components...2` = "Components",
    `Components...9` = "Components",
    `Components...16` = "Components",
    `Df...3` = "df",
    `Df...10` = "df",
    `Df...17` = "df",
    `SumOfSqs...4` = "Sum of squares",
    `SumOfSqs...11` = "Sum of squares",
    `SumOfSqs...18` = "Sum of squares",
    `R2...5` = "R\U00B2",
    `R2...12` = "R\U00B2",
    `R2...19` = "R\U00B2",
    `F...6` = "F",
    `F...13` = "F",
    `F...20` = "F",
    `Pr(>F)...7` = "Pr(>F)",
    `Pr(>F)...14` = "Pr(>F)",
    `Pr(>F)...21` = "Pr(>F)"
  ) 

# ⟞ c. saving outputs -----------------------------------------------------

# ⟞ ⟞ i. PERMANOVA tables -------------------------------------------------

# gtsave(all_permanova,
#        here::here("tables",
#                   "ms-tables",
#                   paste0("tbl-S4_altgower_", today(), ".docx")),
#        vwidth = 1500, vheight = 1000)
# 
# gtsave(permanova_12_table,
#        here::here("tables",
#                   "ms-tables",
#                   paste0("tbl-S4_altgower_12comp_", today(), ".docx")),
#        vwidth = 1500, vheight = 1000)
# 
# gtsave(permanova_3_table,
#        here::here("tables",
#                   "ms-tables",
#                   paste0("tbl-S4_altgower_3comp_", today(), ".docx")),
#        vwidth = 1500, vheight = 1000)

# ⟞ ⟞ ii. control pairwise PERMANOVA table --------------------------------

# save table
# gtsave(all_pairwise,
#        here::here("tables",
#                   "ms-tables",
#                   paste0("tbl-S6_", today(), ".docx")),
#        vwidth = 1500, vheight = 1000)
