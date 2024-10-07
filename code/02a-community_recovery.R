
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. set up -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only have to run this once per session
source(here::here("code", "01a-kelp_recovery.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 1. data frame preparation -----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. functions ----------------------------------------------------------

# This section includes code to prepare data frames for down stream analysis.

group_biomass_wrangle <- function(df) {
  df %>% 
    dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
    group_by(site, year, month, treatment, date) %>% 
    summarize(total_dry = sum(dry_gm2, na.rm = TRUE)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = treatment, values_from = total_dry) %>% 
    mutate(delta_continual = continual - control) %>% 
    select(site, year, month, date, control, continual, delta_continual) %>% 
    # take out years where continual removal hadn't happened yet
    drop_na(delta_continual) %>% 
    select(!delta_continual) %>% 
      mutate(exp_dates = case_when(
        site == "aque" & date >= aque_start_date_continual & date < aque_after_date_continual ~ "during",
        site == "aque" & date >= aque_after_date_continual ~ "after",
        site == "napl" & date >= napl_start_date_continual & date < napl_after_date_continual ~ "during",
        site == "napl" & date >= napl_after_date_continual ~ "after",
        site == "mohk" & date >= mohk_start_date_continual & date < mohk_after_date_continual ~ "during",
        site == "mohk" & date >= mohk_after_date_continual ~ "after",
        site == "carp" & date >= carp_start_date_continual & date < carp_after_date_continual ~ "during",
        site == "carp" & date >= carp_after_date_continual ~ "after"
      ),
      exp_dates = fct_relevel(exp_dates, "during", "after")) %>%
    time_since_columns_continual() %>% 
    kelp_year_column() %>% 
    comparison_column_continual() %>% 
    # make it longer
    pivot_longer(cols = c(control, continual)) %>% 
    # rename columns
    rename(treatment = name, biomass = value) %>% 
    # change treatment names
    mutate(treatment = case_match(
      treatment, 
      "control" ~ "reference", 
      "continual" ~ "removal"),
      treatment = as_factor(treatment)) %>% 
    # create a new sample ID that is site, year, quarter, treatment
    unite("sample_ID", site, date, quarter, treatment, remove = FALSE) %>% 
    # join only with kelp biomass (long format)
    left_join(., select(.data = continual_long, sample_ID, kelp_biomass), 
              by = "sample_ID")
}

delta_biomass_wrangle <- function(df) {
  df %>% 
    dplyr::select(site, year, month, treatment, date, biomass) %>% 
    pivot_wider(names_from = treatment, values_from = biomass) %>% 
    mutate(delta_continual = removal - reference) %>% 
    select(site, year, month, date, reference, removal, delta_continual) %>% 
    # take out years where continual removal hadn't happened yet
    drop_na(delta_continual) %>% 
      mutate(exp_dates = case_when(
        site == "aque" & date >= aque_start_date_continual & date < aque_after_date_continual ~ "during",
        site == "aque" & date >= aque_after_date_continual ~ "after",
        site == "napl" & date >= napl_start_date_continual & date < napl_after_date_continual ~ "during",
        site == "napl" & date >= napl_after_date_continual ~ "after",
        site == "mohk" & date >= mohk_start_date_continual & date < mohk_after_date_continual ~ "during",
        site == "mohk" & date >= mohk_after_date_continual ~ "after",
        site == "carp" & date >= carp_start_date_continual & date < carp_after_date_continual ~ "during",
        site == "carp" & date >= carp_after_date_continual ~ "after"
      ),
      exp_dates = fct_relevel(exp_dates, "during", "after")) %>%
    time_since_columns_continual() %>% 
    kelp_year_column() %>% 
    comparison_column_continual() %>% 
    left_join(., site_quality, by = "site") %>% 
    left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
    rename("site_full" = value) %>% 
    mutate(site_full = fct_relevel(
      site_full, 
      "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
    mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))
}

# ⟞ b. understory macroalgae ----------------------------------------------

algae_continual_long <- biomass %>% 
  # exclude giant kelp
  filter(new_group == "algae" & sp_code != "MAPY") %>% 
  group_biomass_wrangle() %>% 
  mutate(group = "understory algae")

# ⟞ c. sessile invertebrates ----------------------------------------------

epi_continual_long <- biomass %>% 
  filter(new_group == "epilithic.sessile.invert") %>% 
  # take out endos
  filter(taxon_family != "Pholadidae") %>% 
  group_biomass_wrangle() %>% 
  mutate(group = "sessile inverts")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 2. linear models ----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to construct, diagnose, and extract predictions
# from linear models describing the relationship between time since the end of 
# the removal experiment with understory algae and sessile invertebrate biomass. 
# I used a nested data frame to fit the same models to subsets of the full 
# dataset corresponding to each group. 

# The `models` data frame includes the fit models during the kelp removal and
# after the kelp removal, residuals for each, calculations for R2, and model
# predictions. I chose  to visualize the residuals outside of the nested data 
# frame for easier handling. One of the original models only included site as a 
# random effect, but the experimental design necessitated using both site and 
# year as random effects. Thus, the final analysis presented here includes two 
# random effects in the models.

# ⟞ a. model fitting ------------------------------------------------------

models <- bind_rows(algae_continual_long, epi_continual_long) %>% 
  # create a nested data frame
  nest(data = everything(), .by = group) %>% 
  
  # fit model for "during" experimental removal
  mutate(fit_model_during = map(
    data,
    ~ glmmTMB(biomass ~ time_since_end*treatment + (1|site) + (1|year),
              data = .x %>% filter(exp_dates == "during"),
              na.action = na.pass,
              family = ziGamma(link = "log"),
              ziformula = ~1)
  )) %>% 
  # fit model for "after" experimental removal
  mutate(fit_model_after = map(
    data,
    ~ glmmTMB(biomass ~ time_since_end*treatment + (1|site) + (1|year),
              data = .x %>% filter(exp_dates == "after"),
              na.action = na.pass,
              family = ziGamma(link = "log"),
              ziformula = ~1)
  )) %>% 
  
  # simulate residuals for during and after models using DHARMa
  mutate(during_residuals = map(
    fit_model_during,
    ~ simulateResiduals(.x)
  )) %>% 
  mutate(after_residuals = map(
    fit_model_after,
    ~ simulateResiduals(.x)
  )) %>%  
  
  # calculate R2 for during and after models using MuMIn
  mutate(r2_during = map(
    fit_model_during,
    ~ r.squaredGLMM(.x)
  )) %>% 
  mutate(r2_after = map(
    fit_model_after,
    ~ r.squaredGLMM(.x)
  )) %>% 
  
  # generate predictions from models using ggeffects
  mutate(during_predictions = map(
    fit_model_during,
    ~ ggpredict(.x,
                terms = c("time_since_end[-7.5:0 by = 0.25]", "treatment"),
                type = "fixed")
  )) %>% 
  mutate(after_predictions = map(
    fit_model_after,
    ~ ggpredict(.x,
                terms = c("time_since_end[0:6.75 by = 0.25]", "treatment"),
                type = "fixed")
  )) %>% 
  
  # calculate delta biomass
  mutate(delta_biomass = map(
    data,
    ~ delta_biomass_wrangle(.x)
  )) %>% 
  mutate(delta_during_predictions = map(
    during_predictions,
    ~ as.data.frame(.x) %>% 
      select(x, group, predicted) %>% 
      pivot_wider(names_from = group, values_from = predicted) %>% 
      mutate(delta = removal - reference) %>% 
      mutate(exp_dates = "during")
  )) %>% 
  mutate(delta_after_predictions = map(
    after_predictions,
    ~ as.data.frame(.x) %>% 
      select(x, group, predicted) %>% 
      pivot_wider(names_from = group, values_from = predicted) %>% 
      mutate(delta = removal - reference) %>% 
      mutate(exp_dates = "after")
  ))

# This throws warnings from r.squaredGLMM().

# ⟞ b. model diagnostics --------------------------------------------------

# understory algae during
plot(pluck(models, 5, 1))

# sessile invertebrates during
plot(pluck(models, 5, 2))

# understory algae after
plot(pluck(models, 6, 1))

# sessile invertebrates after
plot(pluck(models, 6, 2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 3. model visualizations -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. model predictions --------------------------------------------------

# This section includes code to create panels b and d from Figure 2. These 
# figures depict modeled predictions for understory macroalgae and sessile
# invertebrate biomass as a function of time since the end of the removal and
# treatment (reference or removal). The data underlying these predictions (i.e.
# the observed biomass at each sampling point) is plotted underneath the
# predictions.

# The figures in this section and the following section rely on aesthetics 
# that are saved in objects in the `00-set_up.R` script because they are also
# used in the `01a-kelp_recovery.R` script.

# ⟞ ⟞ i. understory macroalgae --------------------------------------------

overall_algae_predictions <- ggplot() +
  model_predictions_background +
  
  # raw data 
  geom_point(data = algae_continual_long, 
             aes(x = time_since_end, 
                 y = biomass, 
                 color = treatment), 
             shape = 21,
             alpha = 0.15,
             size = 0.75) +
  
  # model predictions
  geom_line(data = models[[9]][[1]], 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +
  geom_ribbon(data = models[[9]][[1]], 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.05) +
  geom_line(data = models[[10]][[1]], 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +
  geom_ribbon(data = models[[10]][[1]], 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.05) +
  # theming
  model_predictions_theme + 
  model_predictions_aesthetics + 
  coord_cartesian(ylim = c(30, 800)) +
  labs(title = "(c)",
       y = "") 

# ⟞ ⟞ ii. sessile invertebrates -------------------------------------------

overall_epi_predictions <- ggplot() +
  model_predictions_background +
  
  # raw data 
  geom_point(data = epi_continual_long, 
             aes(x = time_since_end, 
                 y = biomass, 
                 color = treatment), 
             shape = 21,
             alpha = 0.15,
             size = 0.75) +
  
  # model predictions
  geom_line(data = models[[9]][[2]], 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +
  geom_ribbon(data = models[[9]][[2]], 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.05) +
  geom_line(data = models[[10]][[2]], 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +
  geom_ribbon(data = models[[10]][[2]], 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.05) +
  # theming
  model_predictions_theme + 
  model_predictions_aesthetics + 
  coord_cartesian(ylim = c(5, 155)) +
  labs(title = "(e)",
       y = "") 


# ⟞ b. delta biomass  -----------------------------------------------------

# This section includes code to create panels b and d from Figure 2. These are
# figures depicting the "deltas", or difference in predicted biomass between
# reference and removal plots throughout the experimental removal and recovery
# period. The data for the lines is taken from the `models` object created in
# section 3a. 

# ⟞ ⟞ i. understory macroalgae --------------------------------------------

delta_algae_predictions <- ggplot() +
  model_predictions_background +
  geom_point(data = models[[11]][[1]],
             aes(x = time_since_end, 
                 y = delta_continual), 
             shape = 2, 
             alpha = 0.15,
             size = 0.75) +
  
  # delta biomass
  geom_line(data = models[[12]][[1]], 
            aes(x = x, 
                y = delta), 
            linewidth = 1) +
  geom_line(data = models[[13]][[1]], 
            aes(x = x, 
                y = delta), 
            linewidth = 1) +
  
  delta_aesthetics +
  model_predictions_theme +
  labs(title = "(d)",
       y = "")

# ⟞ ⟞ ii. sessile invertebrates -------------------------------------------

delta_epi_predictions <- ggplot() +
  model_predictions_background +
  geom_point(data = models[[11]][[2]],
             aes(x = time_since_end, 
                 y = delta_continual), 
             shape = 2, 
             alpha = 0.15,
             size = 0.75) +
  
  # delta biomass
  geom_line(data = models[[12]][[2]], 
            aes(x = x, 
                y = delta), 
            linewidth = 1) +
  geom_line(data = models[[13]][[2]], 
            aes(x = x, 
                y = delta), 
            linewidth = 1) +
  
  delta_aesthetics +
  model_predictions_theme +
  labs(title = "(f)",
       y = "")

# ⟞ c. saving outputs -----------------------------------------------------

# This section combines the `overall_kelp` and `delta_kelp_predictions` objects
# from the `01a-kelp_recovery.R` script to create Figure 2 panels a-f. 

# The legend is extracted into a separate object called `legend`, then all the
# plots are arranged into columns.

# The kelp plots are generated in the `01a-kelp_recovery.R` script and combined
# with the other plots in this section.

# ⟞ ⟞ i. extracting legend ------------------------------------------------

obj <- ggplot() +
  model_predictions_background +
  
  # model predictions
  geom_line(data = pluck(models, 9, 1), 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +

  # theming
  model_predictions_theme + 
  model_predictions_aesthetics + 
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))

legend <- get_plot_component(obj, "guide-box-right", return_all = TRUE)

# ⟞ ⟞ ii. plot arrangement ------------------------------------------------

kelp_column <- plot_grid(kelp_title, 
                         overall_kelp_predictions, 
                         delta_kelp_predictions, 
                         nrow = 3, 
                         rel_heights = c(1, 13, 13)) 

algae_column <- plot_grid(algae_title, 
                          overall_algae_predictions, 
                          delta_algae_predictions, 
                          nrow = 3, 
                          rel_heights = c(1, 13, 13))

epi_column <- plot_grid(epi_title, 
                        overall_epi_predictions, 
                        delta_epi_predictions, 
                        nrow = 3, 
                        rel_heights = c(1, 13, 13))

legend_column <- plot_grid(legend, NULL, nrow = 2)

all_columns_with_legend <- plot_grid(kelp_column, 
                                     algae_column, 
                                     epi_column, 
                                     legend_column,
                                     nrow = 1,
                                     rel_widths = c(4, 4, 4, 1.5))

# ⟞ ⟞ iii. saving ---------------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste0("fig-2_new-model_v1_", today(), ".jpg")),
#        plot = all_columns_with_legend,
#        height = 16, width = 24, units = "cm",
#        dpi = 400)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------- 4. timeseries plots --------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to generate timeseries of understory algae and
# sessile invertebrate biomass in reference and removal plots with separate
# panels for each site. These correspond to Supplemental Material figures ____.

# The `raw_biomass_plot_theme()` used below is generated in the 00-set_up.R 
# script because it is also used in the 01a-kelp_recovery.R script.

# Running the code for both these figures will throw an error about removing
# rows outside the scale range; this is due to the "Removal" and "Recovery
# labels that are drawn from a column in the actual data frame, which contains
# one row of labels with the rest of the rows being "NA".

# ⟞ a. understory macroalgae ----------------------------------------------

algae_biomass_timeseries <- algae_continual_long %>% 
  mutate(strip = case_when(
    site == "aque" ~ "(a) Arroyo Quemado",
    site == "napl" ~ "(b) Naples",
    site == "mohk" ~ "(c) Mohawk",
    site == "carp" ~ "(d) Carpinteria"
  )) %>% 
  mutate(removal_annotation = case_when(
    sample_ID == "aque_2010-06-15_Q2_reference" ~ "Removal"
  ),
  recovery_annotation = case_when(
    sample_ID == "aque_2010-06-15_Q2_reference" ~ "Recovery"
  ),
  annotation_y = case_when(
    sample_ID == "aque_2010-06-15_Q2_reference" ~ 360
  )) %>% 
  ggplot(aes(x = time_since_end,
             y = biomass,
             color = treatment,
             group = treatment, 
             linetype = treatment)) +
  model_predictions_background +
  geom_line(alpha = 0.9,
            linewidth = 0.5) +
  geom_text(aes(x = -6.75, 
                y = annotation_y, 
                label = removal_annotation), 
            size = 2, color = "black") +
  geom_text(aes(x = 5.5, 
                y = annotation_y, 
                label = recovery_annotation), 
            size = 2, color = "black") +
  model_predictions_aesthetics +
  raw_biomass_plot_theme +
  labs(y = "Understory macroalgae biomass (dry g/m\U00B2)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.95),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  facet_wrap(~strip, scales = "free_y", nrow = 4)

# ⟞ b. sessile invertebrates ----------------------------------------------

epi_biomass_timeseries <- epi_continual_long %>% 
  mutate(strip = case_when(
    site == "aque" ~ "(a) Arroyo Quemado",
    site == "napl" ~ "(b) Naples",
    site == "mohk" ~ "(c) Mohawk",
    site == "carp" ~ "(d) Carpinteria"
  )) %>% 
  mutate(removal_annotation = case_when(
    sample_ID == "aque_2010-06-15_Q2_reference" ~ "Removal"
  ),
  recovery_annotation = case_when(
    sample_ID == "aque_2010-06-15_Q2_reference" ~ "Recovery"
  ),
  annotation_y = case_when(
    sample_ID == "aque_2010-06-15_Q2_reference" ~ 60
  )) %>% 
  ggplot(aes(x = time_since_end,
             y = biomass,
             color = treatment,
             group = treatment,
             linetype = treatment)) +
  model_predictions_background +
  geom_line(alpha = 0.9,
            linewidth = 0.5) +
  geom_text(aes(x = -6.75, 
                y = annotation_y, 
                label = removal_annotation), 
            size = 2, color = "black") +
  geom_text(aes(x = 5.5, 
                y = annotation_y, 
                label = recovery_annotation), 
            size = 2, color = "black") +
  model_predictions_aesthetics +
  raw_biomass_plot_theme +
  labs(y = "Sessile invertebrate biomass (dry g/m\U00B2)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.91, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  facet_wrap(~strip, scales = "free_y", nrow = 4)

# ⟞ c. saving outputs -----------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_", today(), ".jpg", sep = "")),
#        plot = algae_biomass_timeseries,
#        height = 12, width = 10, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S5_", today(), ".jpg", sep = "")),
#        plot = epi_biomass_timeseries,
#        height = 12, width = 10, units = "cm",
#        dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 5. tables --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to extract the model summaries for each model. It
# relies on the `model_summary_fxn()` from the `00-set_up.R` script.

# ⟞ a. wrangling ----------------------------------------------------------

model_summaries <- models %>% 
  select(group, fit_model_during, fit_model_after) %>% 
  mutate(during_summary = map(
    fit_model_during,
    ~ model_summary_fxn(.x)
  )) %>% 
  mutate(after_summary = map(
    fit_model_after,
    ~ model_summary_fxn(.x) 
  )) %>% 
  # put all columns together
  mutate(merged_tables = pmap(
    list(group, during_summary, after_summary),
    bind_cols
  )) %>% 
  mutate(merged_tables = map(
    merged_tables,
    ~ rename(.x, 
             "group" = "...1",
             "term_removal" = "term...2",
             "estimate_removal" = "estimate...3",
             "p.value_removal" = "p.value...4",
             "ci_interval_removal" = "ci_interval...5",
             "signif_removal" = "signif...6",
             "term_recovery" = "term...7",
             "estimate_recovery" = "estimate...8",
             "p.value_recovery" = "p.value...9",
             "ci_interval_recovery" = "ci_interval...10",
             "signif_recovery" = "signif...11")
  ))

# ⟞ b. table creation -----------------------------------------------------

model_summary_table <- bind_rows(
  kelp_model_summaries_combined,
  pluck(model_summaries, 6, 1),
  pluck(model_summaries, 6, 2)
  ) %>% 
  mutate(group = case_when(
    group == "kelp" ~ "Giant kelp",
    group == "understory algae" ~ "Understory macroalgae",
    group == "sessile inverts" ~ "Sessile invertebrates"
  )) %>% 
  
  # turn object into a flextable, select columns to display
  flextable(col_keys = c("group",
                         "term_removal",
                         "estimate_removal",
                         "p.value_removal",
                         "ci_interval_removal",
                         "term_recovery",
                         "estimate_recovery",
                         "p.value_recovery",
                         "ci_interval_recovery")) %>%
  # add a header row and center align
  add_header_row(top = TRUE, 
                 values = c("", "Kelp removal", "Recovery"), 
                 colwidths = c(1, 4, 4)) %>% 
  align(i = 1, 
        j = NULL, 
        align = "center", 
        part = "header") %>% 
  
  # change the column names
  set_header_labels("group" = " ",
                    "term_removal" = "Term",
                    "estimate_removal" = "Estimate",
                    "p.value_removal" = "p-value",
                    "ci_interval_removal" = "95% CI",
                    "term_recovery" = "Term",
                    "estimate_recovery" = "Estimate",
                    "p.value_recovery" = "p-value",
                    "ci_interval_recovery" = "95% CI") %>% 

  
  # bold p values if they are significant
  style(i = ~ signif_removal == "yes",
        j = "p.value_removal",
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  style(i = ~ signif_recovery == "yes",
        j = "p.value_recovery",
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  
  # add a footnote for 95% CI
  footnote(
    i = 2, 
    j = c(5, 9),
    ref_symbols = "1",
    value = as_paragraph("Confidence interval"),
    part = "header"
  ) %>% 
  
  # merge group cells to create a grouping column
  merge_v(j = ~ group) %>% 
  valign(j = ~ group,
        i = NULL,
        valign = "top") %>% 
  
  # final formatting
  autofit %>% 
  fit_to_width(10) %>% 
  font(fontname = "Times New Roman",
       part = "all")

# ⟞ c. saving output ------------------------------------------------------

# model_summary_table %>%
#   save_as_docx(path = here::here(
#     "tables",
#     "ms-tables",
#     paste("tbl-S1_", today(), ".docx", sep = "")
#     ))
 