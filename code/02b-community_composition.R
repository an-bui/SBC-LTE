# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. set up -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

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

# metadata for continual removal plots
# comm_meta_continual <- comm_meta %>% 
#   filter(treatment == "continual")

# metadata for control plots
# comm_meta_control <- comm_meta %>% 
#   filter(treatment == "control")

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
# during, and after) as predictors, with site as grouping variable. Similarly
# to section b, these analyses are repeated using Bray-Curtis dissimilarity
# and modified Gower dissimilarity.

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

# PERMANOVA results indicated that there was no interaction between treatment 
# and time period for understory algae though there was a significant 
# interaction between these two terms for sessile inverts.

# Thus, in this next section, I compared beta dispersion for treatment and time
# period independently for understory algae, and in combination for sessile
# invertebrates, following guidance by Marti Anderson in this R-sig-eco forum
# post from 2010-09-02: 
# https://stat.ethz.ch/pipermail/r-sig-ecology/2010-September/001524.html

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

algae_betadisper_treatment <- betadisper(
  comm_betadisper[[5]][[2]], comm_betadisper[[3]][[2]]$treatment) 
permutest(algae_betadisper_treatment)
# difference in dispersion between reference and removal (so differences
# in community composition could be due to changes in location and dispersion)

algae_betadisper_time <- betadisper(
  comm_betadisper[[5]][[2]], comm_betadisper[[3]][[2]]$comp_3yrs) 
permutest(algae_betadisper_time)
# no difference in dispersion through time (so differences in community 
# composition are due to changes in location, not dispersion)

epi_betadisper_combo <- betadisper(
  comm_betadisper[[5]][[1]], comm_betadisper[[3]][[1]]$combo
)
permutest(epi_betadisper_combo)
# difference in dispersions
TukeyHSD(epi_betadisper_combo)
# pairwise comparisons of differences between dispersions that are meaningful:
# during control and after control: reference plots differed in dispersion
# in the last 3 years of the experiment and in the recovery period
# start control and start continual: reference and removal plots differed in
# dispersion at the start of the experiment

# ⟞ d. control pairwise comparisons ---------------------------------------

# This section 

control_pairwise <- comm_analyses %>% 
  select(new_group, comp_3yrs_meta, comp_3yrs_mat) %>% 
  mutate(control_start_during_meta = map(
    comp_3yrs_meta,
    ~ filter(.x, comp_3yrs %in% c("start", "during") & treatment == "control"))) %>% 
  mutate(control_during_after_meta = map(
    comp_3yrs_meta,
    ~ filter(.x, comp_3yrs %in% c("during", "after") & treatment == "control"))) %>% 
  mutate(control_start_after_meta = map(
    comp_3yrs_meta,
    ~ filter(.x, comp_3yrs %in% c("start", "after") & treatment == "control"))) %>% 
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
  mutate(control_start_during_mat = map2(
    comp_3yrs_mat, control_start_during_sample_ID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
  mutate(control_during_after_mat = map2(
    comp_3yrs_mat, control_during_after_sample_ID,
    ~ filter(.x, row.names(.x) %in% .y))) %>%   
  mutate(control_start_after_mat = map2(
    comp_3yrs_mat, control_start_after_sample_ID,
    ~ filter(.x, row.names(.x) %in% .y))) %>% 
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
  mutate(start_during_altgower_table = map(
    start_during_altgower,
    ~ as_tibble(.x) %>% 
      mutate(Components = c("Time period", "Residual", "Total")) %>% 
      relocate(Components, .before = Df) %>% 
      mutate(comparison = "start-during") %>% 
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
  mutate(start_after_altgower_table = map(
    start_after_altgower,
    ~ as_tibble(.x) %>% 
      mutate(Components = c("Time period", "Residual", "Total")) %>% 
      relocate(Components, .before = Df) %>% 
      mutate(comparison = "start-after") %>% 
      mutate(across(SumOfSqs:`F`, ~ round(., digits = 2)))
  ))

algae_pairwise <- bind_cols(
  control_pairwise[[16]][[2]], control_pairwise[[17]][[2]], control_pairwise[[18]][[2]]
  ) %>% 
  mutate(group = "understory algae") %>% 
  relocate(group, .before = "Components...1")

inverts_pairwise <- bind_cols(
  control_pairwise[[16]][[1]], control_pairwise[[17]][[1]], control_pairwise[[18]][[1]]
) %>% 
  mutate(group = "sessile inverts") %>% 
  relocate(group, .before = "Components...1")

all_pairwise <- bind_rows(algae_pairwise, inverts_pairwise) %>% 
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
  tab_spanner(label = "during-after",
              columns = c("Components...9", 
                          "Df...10", 
                          "SumOfSqs...11", 
                          "R2...12", 
                          "F...13", 
                          "Pr(>F)...14")) %>% 
  tab_spanner(label = "start-after",
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

all_pairwise

# gtsave(all_pairwise,
#        here::here("tables", "ms-tables", paste("tbl-S6_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)

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

# ⟞ ⟞ i. Bray-Curtis ------------------------------------------------------

# ⟞ ⟞ ⟞  1. analysis ------------------------------------------------------

# Bray-Curtis
algae_pt_bray <- metaMDS(comm_mat_algae, "bray")

# stress plot
stressplot(algae_pt_bray)

# preliminary plot
plot(algae_pt_bray)

# permanova
comp_1yr_meta <- comm_meta %>% 
  drop_na(comp_1yrs)
comp_1yr_sampleID <- comp_1yr_meta %>% 
  pull(sample_ID)
comp_1yr_algae <- comm_mat_algae[comp_1yr_sampleID, ]

algae_pt_perma_1yr <- adonis2(comp_1yr_algae ~ treatment*comp_1yrs, 
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
  unite("combo", comp_3yrs, treatment, sep = "-", remove = FALSE) %>% 
  as.data.frame() %>% 
  mutate(site = as_factor(site))
comp_3yrs_sampleID <- comp_3yrs_meta %>% 
  pull(sample_ID)
comp_3yrs_algae <- comm_mat_algae[comp_3yrs_sampleID, ]
algae_pt_perma_3yrs <- adonis2(comp_3yrs_algae ~ comp_3yrs*treatment, 
                               data = comp_3yrs_meta,
                               strata = comp_3yrs_meta$site)
algae_pt_perma_3yrs 

algae_pairwise <- pairwise.adonis2(comp_3yrs_algae ~ comp_3yrs*treatment, 
                 data = comp_3yrs_meta,
                 strata = comp_3yrs_meta$site)

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "during"))
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_algae[test_ID, ]
algae_start_during_pairwise_bray <- adonis2(test_mat ~ comp_3yrs*treatment,
        data = test_meta,
        strata = test_meta$site) %>% 
  as_tibble() %>% 
  mutate(names = c("comp_3yrs", "treatment", "comp_3yrs:treatment", "Residual", "Total")) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  relocate(names, .before = Df) %>% 
  flextable()

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("after", "during"))
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_algae[test_ID, ]
algae_after_during_pairwise_bray <- adonis2(test_mat ~ comp_3yrs*treatment,
                                            data = test_meta,
                                            strata = test_meta$site) %>% 
  as_tibble() %>% 
  mutate(names = c("comp_3yrs", "treatment", "comp_3yrs:treatment", "Residual", "Total")) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  relocate(names, .before = Df) %>% 
  flextable()

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "after"))
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_algae[test_ID, ]
algae_start_after_pairwise_bray <- adonis2(test_mat ~ comp_3yrs*treatment,
                                            data = test_meta,
                                            strata = test_meta$site) %>% 
  as_tibble() %>% 
  mutate(names = c("comp_3yrs", "treatment", "comp_3yrs:treatment", "Residual", "Total")) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  relocate(names, .before = Df) %>% 
  flextable()

# beta dispersion
algae_pt_dist <- vegdist(comp_3yrs_algae)
# testing with a x b factors following:
# https://stat.ethz.ch/pipermail/r-sig-ecology/2010-September/001524.html
algae_betadisper <- betadisper(algae_pt_dist, comp_3yrs_meta$combo)
algae_permu <- permutest(algae_betadisper, pairwise = TRUE)
algae_pval_bray <- algae_permu$pairwise$permuted %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("pairs") %>% 
  separate_wider_delim(cols = pairs, 
                       names = c("first", "second", "third", "fourth"), 
                       delim = "-") %>% 
  unite("pair1", first, second, sep = "-") %>% 
  unite("pair2", third, fourth, sep = "-") %>% 
  pivot_wider(names_from = pair2, values_from = value) %>% 
  flextable()
algae_pval_bray
TukeyHSD(algae_betadisper)
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
# PTCA, CYOS, CO, CC, DL, EC, R, POLA, EGME, RAT
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

algae_simper_treatment <- simper_algae_treatment$continual_control %>% 
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
  left_join(., comm_meta, by = "sample_ID") # %>% 
  # taking out outlier for visualization
  # filter(sample_ID != "mohk_control_2020-08-12") # nothing there except Gelidium robustum

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
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.6, 1.5)) +
  labs(shape = "Time period",
       color = "Time period", 
       fill = "Time period",
       title = "(a) Removal") +
  theme(legend.position = "none", # c(0.2, 0.8), 
        panel.grid = element_blank(),
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
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.6, 1.5)) +
  # labels
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "(b) Reference") +
  # theme
  theme(legend.position = "none",
        panel.grid = element_blank()) 
algae_pt_bray_control_plot

# control plot with outlier: mohk_control_2020-08-12
algae_pt_bray_control_outlier_plot <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "control", algae_pt_bray_species
) +
  # labels
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "Understory macroalgae | Reference plot") +
  # theme
  theme(legend.position = "none",
        panel.grid = element_blank()) 
algae_pt_bray_control_outlier_plot

# both treatments together
algae_pt_bray_both_plot <- nmds_plot_fxn(
  algae_pt_bray_plotdf, "both", algae_pt_bray_species
) +
  labs(shape = "Time period",
       color = "Time period",
       fill = "Time period",
       linetype = "Treatment",
       alpha = "Treatment") +
  # ellipse labels. clown shit
  # annotate("text", x = -1.2, y = 1.05, label = "Start of", size = 10, col = start_col) +
  # annotate("text", x = -1.2, y = 0.9, label = "removal", size = 10, col = start_col) +
  # annotate("text", x = 1.7, y = 0.85, label = "End of", size = 10, col = during_col) +
  # annotate("text", x = 1.7, y = 0.7, label = "removal", size = 10, col = during_col) +
  # annotate("text", x = -1.4, y = 1.05, label = "Recovery", size = 10, col = after_col) +
  # annotate("text", x = -1.4, y = 0.9, label = "period", size = 10, col = after_col) +
  # stress annotation
  annotate("text", x = -1.2, y = -1.25, label = "Stress = 0.22", size = 4)
algae_pt_bray_both_plot


# ⟞ ⟞ ii. modified Gower --------------------------------------------------

# ⟞ ⟞ ⟞ 1. analysis -------------------------------------------------------

algae_pt_altgower <- metaMDS(comm_mat_algae, "altGower")

# stress plot
stressplot(algae_pt_altgower)

# preliminary plot
plot(algae_pt_altgower)

# permanova
# comp_1yr_meta <- comm_meta %>% 
#   drop_na(comp_1yrs)
# comp_1yr_sampleID <- comp_1yr_meta %>% 
#   pull(sample_ID)
# comp_1yr_algae <- comm_mat_algae[comp_1yr_sampleID, ]

algae_pt_perma_1yr_altgower <- adonis2(comp_1yr_algae ~ treatment*comp_1yrs, 
                              data = comp_1yr_meta,
                              strata = comp_1yr_meta$site,
                              method = "altGower")
algae_pt_perma_1yr_altgower

algae_pt_perma_2yrs_altgower <- adonis2(comp_2yrs_algae ~ treatment*comp_2yrs, 
                               data = comp_2yrs_meta,
                               strata = comp_2yrs_meta$site,
                               method = "altGower")
algae_pt_perma_2yrs_altgower

algae_pt_perma_3yrs_altgower <- adonis2(comp_3yrs_algae ~ treatment*comp_3yrs, 
                               data = comp_3yrs_meta,
                               strata = comp_3yrs_meta$site,
                               method = "altGower")
algae_pt_perma_3yrs_altgower 

algae_pairwise_altgower <- pairwise.adonis2(comp_3yrs_algae ~ comp_3yrs*treatment, 
                                            data = comp_3yrs_meta,
                                            strata = comp_3yrs_meta$site,
                                            method = "altGower")

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "during") & treatment == "control")
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_algae[test_ID, ]
algae_start_during_pairwise_altgower <- adonis2(test_mat ~ comp_3yrs,
        data = test_meta,
        strata = test_meta$site,
        method = "altGower") %>% 
  as_tibble() %>% 
  mutate(names = c("comp_3yrs", "Residual", "Total")) %>% 
  relocate(names, .before = Df) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  flextable()
# when comparing start and during, reference communities are different

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("after", "during") & treatment == "control")
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_algae[test_ID, ]
algae_after_during_pairwise_altgower <- adonis2(test_mat ~ comp_3yrs,
                                                data = test_meta,
                                                strata = test_meta$site,
                                                method = "altGower") %>% 
  as_tibble() %>% 
  mutate(names = c("comp_3yrs", "Residual", "Total")) %>% 
  relocate(names, .before = Df) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  flextable()
# when comparing after and during, reference communities are not different

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "after") & treatment == "control")
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_algae[test_ID, ]
algae_start_after_pairwise_altgower <- adonis2(test_mat ~ comp_3yrs,
                                                data = test_meta,
                                                strata = test_meta$site,
                                                method = "altGower") %>% 
  as_tibble() %>% 
  mutate(names = c("comp_3yrs", "Residual", "Total")) %>% 
  relocate(names, .before = Df) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  flextable()
# when comparing start and after, reference communities are not different

# beta dispersion
algae_pt_dist_altgower <- vegdist(comp_3yrs_algae, method = "altGower")
algae_betadisper_altgower <- betadisper(algae_pt_dist_altgower, comp_3yrs_meta$combo)
permutest(algae_betadisper_altgower, pairwise = TRUE)
TukeyHSD(algae_betadisper_altgower)
# significantly different dispersions
# control dispersions different between during and after

# ⟞ ⟞ ⟞ 2. SIMPER ---------------------------------------------------------

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
# PTCA, CYOS, CO, CC, DL, EC, R, POLA, EGME, RAT
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

algae_simper_treatment <- simper_algae_treatment$continual_control %>% 
  as_tibble() %>% 
  # arrange by greatest to least contribution to dissimilarity
  arrange(-average) %>% 
  head(10)
# PTCA, CYOS, DL, CC, EC, CO, R, POLA, RAT, EGME

# ⟞ ⟞ ⟞ 3. plotting -------------------------------------------------------

# points into data frame for plotting
algae_pt_altgower_plotdf <- scores(algae_pt_altgower, display = "sites") %>% 
  as_tibble(rownames = "sample_ID") %>% 
  # join with metadata
  left_join(., comm_meta, by = "sample_ID") 

# pull top species from simper analysis
simper_algae_spp <- algae_comp3yrs %>% 
  pull(species)

# get species points
algae_pt_altgower_species <- scores(algae_pt_altgower, display = "species", tidy = TRUE) %>% 
  as_tibble(rownames = "sp_code") %>% 
  # keep species from simper analysis only
  filter(sp_code %in% simper_algae_spp) %>% 
  left_join(., spp_names, by = "sp_code")

# continual removal plots only
algae_pt_altgower_continual_plot <- nmds_plot_fxn(
  algae_pt_altgower_plotdf, "continual", algae_pt_altgower_species
) +
  # axis limits
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  labs(shape = "Time period",
       color = "Time period", 
       fill = "Time period",
       title = "(a) Removal") +
  theme(legend.position = "none", # c(0.2, 0.8), 
        panel.grid = element_blank(),
        legend.background = element_blank())
algae_pt_altgower_continual_plot 

algae_pt_altgower_continual_plot_arrows <- algae_pt_altgower_continual_plot +
  geom_text_repel(data = algae_pt_altgower_species,
                  aes(x = NMDS1, y = NMDS2,
                      label = stringr::str_wrap(scientific_name, 4, width = 40)),
                  color = "#C70000", lineheight = 0.8, max.overlaps = 100, size = 1.5) +
  geom_segment(data = algae_pt_altgower_species,
               aes(x = 0, y = 0,
                   xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "#C70000", linewidth = 0.5) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) 

algae_pt_altgower_continual_plot_arrows


# control plots only
algae_pt_altgower_control_plot <- nmds_plot_fxn(
  algae_pt_altgower_plotdf, "control", algae_pt_altgower_species
) +
  # axis limits
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  # labels
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "(b) Reference") +
  # theme
  theme(legend.position = "none",
        panel.grid = element_blank()) 
algae_pt_altgower_control_plot

# both treatments together
algae_pt_altgower_both_plot <- nmds_plot_fxn(
  algae_pt_altgower_plotdf, "both", algae_pt_altgower_species
) +
  labs(shape = "Time period",
       color = "Time period",
       fill = "Time period",
       linetype = "Treatment",
       alpha = "Treatment") 
  # stress annotation
  # annotate("text", x = -1.2, y = -1.25, label = "Stress = 0.22", size = 4)
algae_pt_altgower_both_plot

algae_pt_altgower_start_plot <- algae_pt_altgower_plotdf %>% 
  filter(comp_3yrs == "start") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment, linetype = treatment)) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  geom_point(size = 1, alpha = 0.9) +
  stat_ellipse() +
  scale_linetype_manual(values = c("continual" = 1, "control" = 2)) +
  scale_color_manual(values = c("continual" = removal_col, "control" = reference_col)) +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  annotate("text", x = -0.35, y = 0.5, label = "Reference", color = reference_col, size = 2) +
  annotate("text", x = 0.4, y = 0.1, label = "Removal", color = removal_col, size = 2) +
  labs(title = "(a) Start of removal") +
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
algae_pt_altgower_start_plot

algae_pt_altgower_during_plot <- algae_pt_altgower_plotdf %>% 
  filter(comp_3yrs == "during") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment, linetype = treatment)) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  geom_point(size = 1, alpha = 0.9) +
  stat_ellipse() +
  scale_linetype_manual(values = c("continual" = 1, "control" = 2)) +
  scale_color_manual(values = c("continual" = removal_col, "control" = reference_col)) +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  labs(title = "(b) End of removal") +
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
algae_pt_altgower_during_plot

algae_pt_altgower_after_plot <- algae_pt_altgower_plotdf %>% 
  filter(comp_3yrs == "after") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment, linetype = treatment)) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  geom_point(size = 1, alpha = 0.9) +
  stat_ellipse() +
  scale_linetype_manual(values = c("continual" = 1, "control" = 2)) +
  scale_color_manual(values = c("continual" = removal_col, "control" = reference_col)) +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5)) +
  labs(title = "(c) Recovery period") +
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
algae_pt_altgower_after_plot


# ⟞ b. epi inverts --------------------------------------------------------

# ⟞ ⟞ i. Bray-Curtis ------------------------------------------------------

# ⟞ ⟞ ⟞  1. analysis ------------------------------------------------------

# Bray-Curtis
epi_pt_bray <- metaMDS(comm_mat_epi, "bray")

# stress plot
stressplot(epi_pt_bray)

# preliminary plot
plot(epi_pt_bray)

# permanova
comp_1yr_meta <- comm_meta %>% 
  drop_na(comp_1yrs)
comp_1yr_sampleID <- comp_1yr_meta %>% 
  pull(sample_ID)
comp_1yr_epi <- comm_mat_epi[comp_1yr_sampleID, ]
epi_pt_perma_1yr <- adonis2(comp_1yr_epi ~ treatment*comp_1yrs, 
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
  unite("combo", comp_3yrs, treatment, sep = "-", remove = FALSE) %>% 
  as.data.frame() %>% 
  mutate(site = as_factor(site),
         treatment = as_factor(treatment))
comp_3yrs_sampleID <- comp_3yrs_meta %>% 
  pull(sample_ID)
comp_3yrs_epi <- comm_mat_epi[comp_3yrs_sampleID, ]
epi_pt_perma_3yrs <- adonis2(comm_mat_epi ~ comp_3yrs*treatment,
                             data = comp_3yrs_meta,
                             strata = comp_3yrs_meta$site)
epi_pt_perma_3yrs # same as 2 years

epi_pairwise <- pairwise.adonis2(comm_mat_epi ~ comp_3yrs*treatment,
                                 data = comp_3yrs_meta,
                                 strata = comp_3yrs_meta$site)
epi_pairwise
# all different from each other

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "during"))
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_epi[test_ID, ]
adonis2(test_mat ~ comp_3yrs*treatment,
        data = test_meta,
        strata = test_meta$site)

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "after"))
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_epi[test_ID, ]
adonis2(test_mat ~ comp_3yrs*treatment,
        data = test_meta,
        strata = test_meta$site)

test_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("after", "during"))
test_ID <- test_meta %>% pull(sample_ID)
test_mat <- comm_mat_epi[test_ID, ]
adonis2(test_mat ~ comp_3yrs*treatment,
        data = test_meta,
        strata = test_meta$site)

# beta dispersion
epi_pt_dist <- vegdist(comm_mat_epi)
epi_betadisper <- betadisper(epi_pt_dist, comp_3yrs_meta$combo)
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

epi_simper_treatment <- simper_epi_treatment$continual_control %>% 
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
  scale_x_continuous(limits = c(-1.7, 1.2)) +
  scale_y_continuous(limits = c(-1.6, 1.3)) +
  theme(legend.position = "none",
        panel.grid = element_blank()
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
  scale_x_continuous(limits = c(-1.4, 1.4), breaks = seq(-1, 1, by = 1)) +
  scale_y_continuous(limits = c(-1.35, 1.45), breaks = seq(-1, 1, by = 1)) +
  labs(title = "(d) Reference") +
  theme(legend.position = "none",
        panel.grid = element_blank()) 
epi_pt_bray_control_plot

# both treatments together
epi_pt_bray_both_plot <- nmds_plot_fxn(
  epi_pt_bray_plotdf, "both", epi_pt_bray_species
) +
  scale_x_continuous(limits = c(-1.75, 1.4), breaks = seq(-1, 1, by = 1)) +
  scale_y_continuous(limits = c(-1.7, 1.45), breaks = seq(-1, 1, by = 1)) +
  annotate("text", x = -1, y = -1.7, label = "Stress = 0.2", size = 2)
epi_pt_bray_both_plot

# ⟞ ⟞ ii. modified Gower --------------------------------------------------


# ⟞ ⟞ ⟞ 1. analysis -------------------------------------------------------

epi_pt_altgower <- metaMDS(comm_mat_epi, "altGower")

# stress plot
stressplot(epi_pt_altgower)

# preliminary plot
plot(epi_pt_altgower)

# permanova
# comp_1yr_meta <- comm_meta %>% 
#   drop_na(comp_1yrs)
# comp_1yr_sampleID <- comp_1yr_meta %>% 
#   pull(sample_ID)
# comp_1yr_epi <- comm_mat_epi[comp_1yr_sampleID, ]

epi_pt_perma_1yr_altgower <- adonis2(comp_1yr_epi ~ treatment*comp_1yrs, 
                                       data = comp_1yr_meta,
                                       strata = comp_1yr_meta$site,
                                       method = "altGower")
epi_pt_perma_1yr_altgower

epi_pt_perma_2yrs_altgower <- adonis2(comp_2yrs_epi ~ treatment*comp_2yrs, 
                                        data = comp_2yrs_meta,
                                        strata = comp_2yrs_meta$site,
                                        method = "altGower")
epi_pt_perma_2yrs_altgower

epi_pt_perma_3yrs_altgower <- adonis2(comp_3yrs_epi ~ comp_3yrs*treatment, 
                                        data = comp_3yrs_meta,
                                        strata = comp_3yrs_meta$site,
                                        method = "altGower")
epi_pt_perma_3yrs_altgower 

epi_pairwise_altgower <- pairwise.adonis2(comp_3yrs_epi ~ comp_3yrs*treatment, 
                                            data = comp_3yrs_meta,
                                            strata = comp_3yrs_meta$site,
                                            method = "altGower")

epi_start_during_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "during") & treatment == "control")
epi_start_during_ID <- epi_start_during_meta %>% pull(sample_ID)
epi_start_during_mat <- comm_mat_epi[epi_start_during_ID, ]
epi_start_during_pairwise_altgower <- adonis2(
    epi_start_during_mat ~ comp_3yrs,
    data = epi_start_during_meta,
    strata = epi_start_during_meta$site,
    method = "altGower") %>% 
  as_tibble() %>% 
  mutate(names = c("Time period", "Residual", "Total")) %>% 
  relocate(names, .before = Df) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  mutate(group = "inverts", 
         comparison = "start-during") 
# when comparing start and during, composition is different

epi_during_after_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("after", "during") & treatment == "control")
epi_during_after_ID <- epi_during_after_meta %>% pull(sample_ID)
epi_during_after_mat <- comm_mat_epi[epi_during_after_ID, ]
epi_during_after_pairwise_altgower <- adonis2(
  epi_during_after_mat ~ comp_3yrs,
    data = epi_during_after_meta,
    strata = epi_during_after_meta$site,
    method = "altGower") %>% 
  as_tibble() %>% 
  mutate(names = c("Time period", "Residual", "Total")) %>% 
  relocate(names, .before = Df) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  mutate(group = "inverts", 
         comparison = "during-after") 
# when comparing after and during, composition is different

start_after_meta <- comp_3yrs_meta %>% 
  filter(comp_3yrs %in% c("start", "after") & treatment == "control")
start_after_ID <- test_meta %>% pull(sample_ID)
start_after_mat <- comm_mat_epi[test_ID, ]
epi_start_after_pairwise_altgower <- adonis2(start_after_mat ~ comp_3yrs,
                                             data = start_after_meta,
                                             strata = start_after_meta$site,
                                             method = "altGower") %>% 
  as_tibble() %>% 
  mutate(names = c("Time period", "Residual", "Total")) %>% 
  relocate(names, .before = Df) %>% 
  mutate(across(SumOfSqs:`F`, ~ round(., digits = 2))) %>% 
  mutate(group = "inverts", 
       comparison = "start-after") 
# when comparing start and after, composition is different

table <- bind_rows(epi_start_during_pairwise_altgower, 
                   epi_after_during_pairwise_altgower,
                   epi_start_after_pairwise_altgower) %>% 
  relocate(group, .before = names) %>% 
  relocate(comparison, .after = group) %>% 
  flextable()

table

# beta dispersion
epi_pt_dist_altgower <- vegdist(comp_3yrs_epi, method = "altGower")
epi_betadisper_altgower <- betadisper(epi_pt_dist_altgower, comp_3yrs_meta$combo)
permutest(epi_betadisper_altgower, pairwise = TRUE)
TukeyHSD(epi_betadisper_altgower)
# significantly different dispersions
# control dispersions different between during and after, start and after
# start control and start removal dispersions different

# ⟞ ⟞ ⟞ 2. SIMPER ---------------------------------------------------------


# ⟞ ⟞ ⟞ 3. plotting -------------------------------------------------------

# points into data frame for plotting
epi_pt_altgower_plotdf <- scores(epi_pt_altgower, display = "sites") %>% 
  as_tibble(rownames = "sample_ID") %>% 
  # join with metadata
  left_join(., comm_meta, by = "sample_ID") 

# pull top species from simper analysis
simper_epi_spp <- epi_comp3yrs %>% 
  pull(species)

# get species points
epi_pt_altgower_species <- scores(epi_pt_altgower, display = "species", tidy = TRUE) %>% 
  as_tibble(rownames = "sp_code") %>% 
  # keep species from simper analysis only
  filter(sp_code %in% simper_epi_spp) %>% 
  left_join(., spp_names, by = "sp_code")

# continual removal plots only
epi_pt_altgower_continual_plot <- nmds_plot_fxn(
  epi_pt_altgower_plotdf, "continual", epi_pt_altgower_species
) +
  # axis limits
  scale_x_continuous(limits = c(-0.35, 0.25)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  labs(shape = "Time period",
       color = "Time period", 
       fill = "Time period",
       title = "(c) Removal") +
  theme(legend.position = "none", # c(0.2, 0.8), 
        panel.grid = element_blank(),
        legend.background = element_blank())
epi_pt_altgower_continual_plot 

epi_pt_altgower_continual_plot_arrows <- epi_pt_altgower_continual_plot +
  geom_text_repel(data = epi_pt_altgower_species,
                  aes(x = NMDS1, y = NMDS2,
                      label = stringr::str_wrap(scientific_name, 4, width = 40)),
                  color = "#C70000", lineheight = 0.8, max.overlaps = 100, size = 1.5) +
  geom_segment(data = epi_pt_altgower_species,
               aes(x = 0, y = 0,
                   xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "#C70000", linewidth = 0.5) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) 

epi_pt_altgower_continual_plot_arrows


# control plots only
epi_pt_altgower_control_plot <- nmds_plot_fxn(
  epi_pt_altgower_plotdf, "control", epi_pt_altgower_species
) +
  # axis limits
  scale_x_continuous(limits = c(-0.35, 0.25)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  # labels
  labs(shape = "Site",
       color = "Time period", fill = "Time period",
       title = "(d) Reference") +
  # theme
  theme(legend.position = "none",
        panel.grid = element_blank()) 
epi_pt_altgower_control_plot

# both treatments together
epi_pt_altgower_both_plot <- nmds_plot_fxn(
  epi_pt_altgower_plotdf, "both", epi_pt_altgower_species
) +
  labs(shape = "Time period",
       color = "Time period",
       fill = "Time period",
       linetype = "Treatment",
       alpha = "Treatment") 
# stress annotation
# annotate("text", x = -1.2, y = -1.25, label = "Stress = 0.22", size = 4)
epi_pt_altgower_both_plot

epi_pt_altgower_start_plot <- epi_pt_altgower_plotdf %>% 
  filter(comp_3yrs == "start") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment, linetype = treatment)) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  geom_point(size = 1, alpha = 0.9) +
  stat_ellipse() +
  scale_linetype_manual(values = c("continual" = 1, "control" = 2)) +
  scale_color_manual(values = c("continual" = removal_col, "control" = reference_col)) +
  scale_x_continuous(limits = c(-0.35, 0.25)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  annotate("text", x = 0.15, y = 0.2, label = "Reference", color = reference_col, size = 2) +
  annotate("text", x = -0.2, y = -0.2, label = "Removal", color = removal_col, size = 2) +
  labs(title = "(d) Start of removal") +
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
epi_pt_altgower_start_plot

epi_pt_altgower_during_plot <- epi_pt_altgower_plotdf %>% 
  filter(comp_3yrs == "during") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment, linetype = treatment)) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  geom_point(size = 1, alpha = 0.9) +
  stat_ellipse() +
  scale_linetype_manual(values = c("continual" = 1, "control" = 2)) +
  scale_color_manual(values = c("continual" = removal_col, "control" = reference_col)) +
  scale_x_continuous(limits = c(-0.35, 0.25)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  labs(title = "(e) End of removal") +
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
epi_pt_altgower_during_plot

epi_pt_altgower_after_plot <- epi_pt_altgower_plotdf %>% 
  filter(comp_3yrs == "after") %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment, linetype = treatment)) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 0, color = "grey", lty = 2) +
  geom_hline(yintercept = 0, color = "grey", lty = 2) +
  geom_point(size = 1, alpha = 0.9) +
  stat_ellipse() +
  scale_linetype_manual(values = c("continual" = 1, "control" = 2)) +
  scale_color_manual(values = c("continual" = removal_col, "control" = reference_col)) +
  scale_x_continuous(limits = c(-0.35, 0.25)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  labs(title = "(f) Recovery period") +
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
epi_pt_altgower_after_plot

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

# ⟞ a. Bray-Curtis --------------------------------------------------------

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
  mutate(group = case_when(
    str_detect(model_1yr1, "algae") ~ "Understory macroalgae",
    str_detect(model_1yr1, "epi") ~ "Sessile invertebrates"
  )) %>% 
  # take out unwanted columns
  select(!c("model_1yr1", "model_2yrs1", "model_3yrs1", 
            "SumOfSqs_1yr1", "SumOfSqs_2yrs1", "SumOfSqs_3yrs1",
            "R2_1yr1", "R2_2yrs1", "R2_3yrs1")) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  # group labels
  # tab_row_group(
  #   label = "Sessile invertebrates", rows = 6:10
  # ) %>% 
  # tab_row_group(
  #   label = "Understory macroalgae", rows = 1:5
  # ) %>% 
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
    p_1yr1 = "p-value", 
    variables_2yrs1 = "Source of variation",
    Df_2yrs1 = "df",
    F_2yrs1 = "pseudo-F",
    p_2yrs1 = "p-value", 
    variables_3yrs1 = "Source of variation",
    Df_3yrs1 = "Degrees of freedom",
    F_3yrs1 = "pseudo-F",
    p_3yrs1 = "p-value"
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
#        here::here("tables", "ms-tables", paste("tbl-S4_bray_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)


# ⟞ b. modified Gower -----------------------------------------------------


anova_1yr_altgower_tables <- rbind(anova_summary_fxn(algae_pt_perma_1yr_altgower), anova_summary_fxn(epi_pt_perma_1yr_altgower)) %>% 
  rename_with(., .fn = ~paste(., "_1yr_altgower", sep = "", .cols = everything(cols))) 

anova_1yr_altgower_tables 

anova_2yrs_altgower_tables <- rbind(anova_summary_fxn(algae_pt_perma_2yrs_altgower), anova_summary_fxn(epi_pt_perma_2yrs_altgower)) %>% 
  rename_with(., .fn = ~paste(., "_2yrs_altgower", sep = "", .cols = everything(cols)))
anova_2yrs_altgower_tables 

anova_3yrs_altgower_tables <- rbind(anova_summary_fxn(algae_pt_perma_3yrs_altgower), anova_summary_fxn(epi_pt_perma_3yrs_altgower)) %>% 
  rename_with(., .fn = ~paste(., "_3yrs_altgower", sep = "", .cols = everything(cols)))
anova_3yrs_altgower_tables

anova_12_tables_altgower <- cbind(anova_1yr_altgower_tables, anova_2yrs_altgower_tables) %>% 
  mutate(group = case_when(
    str_detect(model_1yr_altgower1, "algae") ~ "Understory macroalgae",
    str_detect(model_1yr_altgower1, "epi") ~ "Sessile invertebrates"
  )) %>% 
  # take out unwanted columns
  select(!c("model_1yr_altgower1", "model_2yrs_altgower1", 
            "SumOfSqs_1yr_altgower1", "SumOfSqs_2yrs_altgower1",
            "R2_1yr_altgower1", "R2_2yrs_altgower1")) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  # group labels
  # tab_row_group(
  #   label = "Sessile invertebrates", rows = 6:10
  # ) %>% 
  # tab_row_group(
  #   label = "Understory macroalgae", rows = 1:5
  # ) %>% 
  # 1, 2, and 3 year comparisons
  tab_spanner(
    label = "1 year comparison",
    columns = c(variables_1yr_altgower1, Df_1yr_altgower1, F_1yr_altgower1, p_1yr_altgower1)
  ) %>% 
  tab_spanner(
    label = "2 year comparison",
    columns = c(variables_2yrs_altgower1, Df_2yrs_altgower1, F_2yrs_altgower1, p_2yrs_altgower1)
  ) %>% 
  # change column names
  cols_label(
    variables_1yr_altgower1 = "Source of variation",
    Df_1yr_altgower1 = "df",
    F_1yr_altgower1 = "pseudo-F",
    p_1yr_altgower1 = "p-value", 
    variables_2yrs_altgower1 = "Source of variation",
    Df_2yrs_altgower1 = "df",
    F_2yrs_altgower1 = "pseudo-F",
    p_2yrs_altgower1 = "p-value"
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
      columns = p_1yr_altgower1,
      rows = p_1yr_altgower1 < 0.05
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = p_2yrs_altgower1,
      rows = p_2yrs_altgower1 < 0.05
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
anova_12_tables_altgower

anova_3_tables_altgower <- anova_3yrs_altgower_tables %>% 
  mutate(group = case_when(
    str_detect(model_3yrs_altgower1, "algae") ~ "Understory macroalgae",
    str_detect(model_3yrs_altgower1, "epi") ~ "Sessile invertebrates"
  )) %>% 
  # take out unwanted columns
  select(!c("model_3yrs_altgower1", 
            "SumOfSqs_3yrs_altgower1",
            "R2_3yrs_altgower1")) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  tab_spanner(
    label = "3 year comparison",
    columns = c(variables_3yrs_altgower1, Df_3yrs_altgower1, F_3yrs_altgower1, p_3yrs_altgower1)
  ) %>% 
  # change column names
  cols_label(
    variables_3yrs_altgower1 = "Source of variation",
    Df_3yrs_altgower1 = "Degrees of freedom",
    F_3yrs_altgower1 = "pseudo-F",
    p_3yrs_altgower1 = "p-value"
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
      columns = p_3yrs_altgower1,
      rows = p_3yrs_altgower1 < 0.05
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
anova_3_tables_altgower

anova_together_tables_altgower <- cbind(anova_1yr_altgower_tables, anova_2yrs_altgower_tables, anova_3yrs_altgower_tables) %>% 
  mutate(group = case_when(
    str_detect(model_1yr_altgower1, "algae") ~ "Understory macroalgae",
    str_detect(model_1yr_altgower1, "epi") ~ "Sessile invertebrates"
  )) %>% 
  # take out unwanted columns
  select(!c("model_1yr_altgower1", "model_2yrs_altgower1", "model_3yrs_altgower1", 
            "SumOfSqs_1yr_altgower1", "SumOfSqs_2yrs_altgower1", "SumOfSqs_3yrs_altgower1",
            "R2_1yr_altgower1", "R2_2yrs_altgower1", "R2_3yrs_altgower1")) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  # group labels
  # tab_row_group(
  #   label = "Sessile invertebrates", rows = 6:10
  # ) %>% 
  # tab_row_group(
  #   label = "Understory macroalgae", rows = 1:5
  # ) %>% 
  # 1, 2, and 3 year comparisons
  tab_spanner(
    label = "1 year comparison",
    columns = c(variables_1yr_altgower1, Df_1yr_altgower1, F_1yr_altgower1, p_1yr_altgower1)
  ) %>% 
  tab_spanner(
    label = "2 year comparison",
    columns = c(variables_2yrs_altgower1, Df_2yrs_altgower1, F_2yrs_altgower1, p_2yrs_altgower1)
  ) %>% 
  tab_spanner(
    label = "3 year comparison",
    columns = c(variables_3yrs_altgower1, Df_3yrs_altgower1, F_3yrs_altgower1, p_3yrs_altgower1)
  ) %>% 
  # change column names
  cols_label(
    variables_1yr_altgower1 = "Source of variation",
    Df_1yr_altgower1 = "df",
    F_1yr_altgower1 = "pseudo-F",
    p_1yr_altgower1 = "p-value", 
    variables_2yrs_altgower1 = "Source of variation",
    Df_2yrs_altgower1 = "df",
    F_2yrs_altgower1 = "pseudo-F",
    p_2yrs_altgower1 = "p-value", 
    variables_3yrs_altgower1 = "Source of variation",
    Df_3yrs_altgower1 = "Degrees of freedom",
    F_3yrs_altgower1 = "pseudo-F",
    p_3yrs_altgower1 = "p-value"
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
      columns = p_1yr_altgower1,
      rows = p_1yr_altgower1 < 0.05
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = p_2yrs_altgower1,
      rows = p_2yrs_altgower1 < 0.05
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = p_3yrs_altgower1,
      rows = p_3yrs_altgower1 < 0.05
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
anova_together_tables_altgower

# gtsave(anova_together_tables_altgower,
#        here::here("tables", "ms-tables", paste("tbl-S4_altgower_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)

# gtsave(anova_12_tables_altgower,
#        here::here("tables", "ms-tables", paste("tbl-S4_altgower_12comp_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)

# gtsave(anova_3_tables_altgower,
#        here::here("tables", "ms-tables", paste("tbl-S4_altgower_3comp_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)

##########################################################################-
# 5. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. continual removal --------------------------------------------------

comm_comp_together_bray <- algae_pt_bray_continual_plot + epi_pt_bray_continual_plot 

comm_comp_together_arrows_bray <- algae_pt_bray_continual_plot_arrows + epi_pt_bray_continual_plot_arrows

comm_comp_together_altgower <- algae_pt_altgower_continual_plot + epi_pt_altgower_continual_plot 

comm_comp_together_arrows_altgower <- algae_pt_altgower_continual_plot_arrows + epi_pt_altgower_continual_plot_arrows

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_bray_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together_bray,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S8_bray_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together_arrows_bray,
#        height = 8, width = 16, units = "cm",
#        dpi = 400)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_altgower_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together_altgower,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S8_altgower_", today(), ".jpg", sep = "")),
#        plot = comm_comp_together_arrows_altgower,
#        height = 8, width = 16, units = "cm",
#        dpi = 400)

# ⟞ b. control ------------------------------------------------------------

comm_comp_control <- algae_pt_bray_control_plot + epi_pt_bray_control_plot 

algae_pt_bray_control_outlier_plot

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S7_", today(), ".jpg", sep = "")),
#        plot = comm_comp_control,
#        height = 6, width = 16, units = "cm",
#        dpi = 300)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S3_", today(), ".jpg", sep = "")),
#        plot = algae_pt_bray_control_outlier_plot,
#        height = 10, width = 10, units = "cm",
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

# ⟞ ⟞ i. Bray-Curtis ------------------------------------------------------

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
#                   paste("fig-3_bray_", today(), ".jpg", sep = "")),
#        plot = final_plot,
#        height = 13, width = 15, units = "cm",
#        dpi = 300)


# ⟞ ⟞ ii. modified Gower --------------------------------------------------

# putting algae plots together
algae_altgower_plots <- plot_grid(algae_pt_altgower_continual_plot, algae_pt_altgower_control_plot, 
                         nrow = 2)

# putting invert plots together
epi_altgower_plots <- plot_grid(epi_pt_altgower_continual_plot, epi_pt_altgower_control_plot, 
                       nrow = 2)

# group plots with labels
algae_altgower_labelled <- plot_grid(algae_title, algae_altgower_plots, 
                            rel_heights = c(1, 12), 
                            ncol = 1)


epi_altgower_labelled <- plot_grid(epi_title, epi_altgower_plots, 
                          rel_heights = c(1, 12), 
                          ncol = 1)

# getting the legend as separate object
plot_legend_altgower <- nmds_plot_fxn(
  algae_pt_altgower_plotdf, "continual", algae_pt_altgower_species
) +
  labs(shape = "Time period",
       color = "Time period", 
       fill = "Time period") 
legend <- cowplot::get_legend(plot_legend_altgower)

# putting plots together
all_altgower_plots <- plot_grid(algae_altgower_labelled, epi_altgower_labelled, ncol = 2,
                       rel_widths = c(1, 1))

# putting plots with legend
final_altgower_plot <- plot_grid(all_altgower_plots, legend, ncol = 2, rel_widths = c(1, 0.3))

# saving
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_altgower_", today(), ".jpg", sep = "")),
#        plot = final_altgower_plot,
#        height = 13, width = 15, units = "cm",
#        dpi = 300)

# ⟞ ⟞ iii. modified Gower v2 ----------------------------------------------

# putting algae plots together
algae_altgower_plots <- plot_grid(algae_pt_altgower_start_plot, 
                                  algae_pt_altgower_during_plot, 
                                  algae_pt_altgower_after_plot,
                                  ncol = 1)

epi_altgower_plots <- plot_grid(epi_pt_altgower_start_plot, 
                                  epi_pt_altgower_during_plot, 
                                  epi_pt_altgower_after_plot,
                                ncol = 1)

# group plots with labels
algae_altgower_labelled <- plot_grid(algae_title, algae_altgower_plots, 
                                     rel_heights = c(1, 9), 
                                     ncol = 1)


epi_altgower_labelled <- plot_grid(epi_title, epi_altgower_plots, 
                                   rel_heights = c(1, 9), 
                                   ncol = 1)

# putting plots together
all_altgower_plots <- plot_grid(algae_altgower_labelled, epi_altgower_labelled, ncol = 2,
                                rel_widths = c(1, 1))

# saving
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_altgower_v2_", today(), ".jpg", sep = "")),
#        plot = all_altgower_plots,
#        height = 16, width = 12, units = "cm",
#        dpi = 300)


