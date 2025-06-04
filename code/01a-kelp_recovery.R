
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. set up -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only have to run this once per session
source(here::here("code", "00a-set_up.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 1. data frame preparation -----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create data frames that are used in downstream
# analysis in this script and joined with other data frames in other scripts. 
# Thus, this script is the "source" for downstream scripts.

# ⟞ a. delta kelp biomass -------------------------------------------------

# This data frame includes the calculation of "delta" biomass: the difference
# in kelp biomass between removal (continual) and reference (control) plots,
# along with the metadata for each sampling event.

# The biomass dataset includes data from an experimental removal where kelp was
# removed once a year (as opposed to quarterly). The annual removal experiment
# started before the continual removal experiment, so the dataset includes 
# observations from before the focal dataset for this project. Thus, there is a
# step in the wrangling process in which observations before the continual
# removal experiment started are dropped.

# This code also relies on wrangling functions that are created in the 
# `00-set_up.R` script.

delta_continual <- biomass %>% 
  filter(sp_code == "MAPY" & treatment %in% c("control", "continual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2, exp_dates:kelp_year) %>% 
  pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
  # calculate delta
  mutate(delta_continual = continual - control) %>% 
  left_join(., site_quality, by = "site") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(
    site_full, 
    "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  mutate(site = fct_relevel(
    site, 
    "aque", "napl", "mohk", "carp")) %>% 
  # take out the missing survey from Naples
  drop_na(delta_continual) %>% 
  unite("sample_ID_short", site, date, sep = "_", remove = FALSE)

# ⟞ b. kelp biomass in long format ----------------------------------------

# kelp biomass in long format
continual_long <- delta_continual %>% 
  select(!delta_continual) %>% 
  pivot_longer(cols = c(control, continual)) %>% 
  rename(kelp_biomass = value, treatment = name) %>% 
  mutate(treatment = case_match(
    treatment, 
    "control" ~ "reference", 
    "continual" ~ "removal"),
    treatment = as_factor(treatment)) %>% 
  unite("sample_ID", site, date, quarter, treatment, remove = FALSE) %>% 
  # take out the missing NAPL survey
  filter(sample_ID != "napl_2014-11-14_Q4_removal") %>% 
  # calculate variation by observation
  group_by(site, exp_dates, treatment) %>% 
  mutate(mean = mean(kelp_biomass),
         variation = 
           ((kelp_biomass - mean(kelp_biomass))/mean(kelp_biomass))^2) %>% 
  ungroup()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 2. linear models ----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to model the effect of time since the end of the
# experiment and treatment (reference or removal) on kelp biomass with random
# effects of site and year using a zero-inflated Gamma distribution with a 
# log link. 

# Before settling on this distribution, we tried many different models, which 
# are archived in the GitHub repo for this project in the `main` branch in the
# `old-code` folder.

kelp_models <- continual_long %>% 
  nest(data = everything(), .by = exp_dates) %>% 
  # fit the model
  mutate(kelp_model = map(
    data,
    ~ glmmTMB(
        kelp_biomass ~ time_since_end*treatment + (1|site) + (1|year),
        data = .x,
        family = ziGamma(link = "log"),
        ziformula = ~1
    )
  )) %>% 
  # simulate residuals
  mutate(residuals = map(
    kelp_model,
    ~ simulateResiduals(.x, quantreg = TRUE)
  )) %>% 
  # calculate R2
  mutate(r2 = map(
    kelp_model,
    ~ r.squaredGLMM(.x)
  )) %>% 
  # generate model predictions
  mutate(predictions = case_when(
    exp_dates == "during" ~ map(
      kelp_model,
      ~ ggpredict(.x,
                  terms = c("time_since_end[-7.25:0 by = 0.25]", "treatment"),
                  type = "fixed")
    ),
    exp_dates == "after" ~ map(
      kelp_model,
      ~ ggpredict(.x,
                  terms = c("time_since_end[0:6.75 by = 0.25]", "treatment"),
                  type = "fixed")
    )
  )) %>% 
  # generate deltas (removal - reference biomass)
  mutate(delta_predictions = map2(
    predictions, exp_dates,
    ~ as.data.frame(.x) %>% 
      select(x, group, predicted) %>% 
      pivot_wider(names_from = group, values_from = predicted) %>% 
      mutate(delta = removal - reference) %>% 
      mutate(exp_dates = .y)
  ))

# this throws a warning from `r.squaredGLMM`

# residuals for model during experimental removal
plot(pluck(kelp_models, 4, 1))

# residuals for model after experimental removal (during recovery period)
plot(pluck(kelp_models, 4, 2))

# plotting residuals against predictors during experimental removal
plotResiduals(pluck(kelp_models, 3, 1), 
              form = pluck(kelp_models, 2, 1)$time_since_end)

# plotting residuals against predictors after experimental removal
plotResiduals(pluck(kelp_models, 3, 2), 
              form = pluck(kelp_models, 2, 2)$time_since_end)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 3. model visualizations -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. model predictions --------------------------------------------------

# This section includes code to create panels a and b from Figure 2. These 
# figures depict modeled predictions for giant kelp biomass as a function of 
# time since the end of the removal and treatment (reference or removal). The 
# data underlying these predictions (i.e. the observed biomass at each sampling 
# point) is plotted underneath the predictions.

# The  figures in this section and the following section rely on aesthetics 
# that are saved in objects in the `00-set_up.R` script. These plots are then
# combined with plots created in the `02a-community_recovery.R` script to 
# create the final figure.

overall_kelp_predictions <- ggplot() +
  model_predictions_background +
  
  # raw data 
  geom_point(data = continual_long, 
             aes(x = time_since_end, 
                 y = kelp_biomass, 
                 color = treatment), 
             shape = 21,
             alpha = 0.7,
             size = 0.75) +
  
  # model predictions
  geom_ribbon(data = pluck(kelp_models, 6, 1), 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.2,
              fill = "grey") +
  geom_line(data = pluck(kelp_models, 6, 1), 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +
  geom_ribbon(data = pluck(kelp_models, 6, 2), 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.2,
              fill = "grey") +
  geom_line(data = pluck(kelp_models, 6, 2), 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +

  # theming
  model_predictions_theme +
  model_predictions_aesthetics +
  transparent_theme +
  coord_cartesian(ylim = c(-10, 1650))

overall_kelp_predictions

recovery_kelp_only <- ggplot() +
  model_predictions_background +
  
  # raw data 
  geom_point(data = continual_long |> filter(exp_dates == "after"), 
             aes(x = time_since_end, 
                 y = kelp_biomass, 
                 color = treatment), 
             shape = 21,
             alpha = 0.7,
             size = 0.75) +
  
  # model predictions
  # 
  geom_ribbon(data = pluck(kelp_models, 6, 2), 
              aes(x = x, 
                  ymax = conf.high, 
                  ymin = conf.low, 
                  group = group), 
              alpha = 0.2,
              fill = "grey") +
  geom_line(data = pluck(kelp_models, 6, 2), 
            aes(x = x, 
                y = predicted, 
                color = group, 
                linetype = group), 
            linewidth = 1) +

  # theming
  model_predictions_theme +
  scale_color_manual(values = c(reference = reference_col, 
                                removal = removal_col),
                     labels = c(reference = "Reference", 
                                removal = "Removal")) +
  scale_linetype_manual(values = c(reference = "22", 
                                   removal = "solid"),
                        labels = c(reference = "Reference", 
                                   removal = "Removal")) +
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1) +
  labs(x = "Time since end of removal (years)", 
       y = "Biomass (dry g/m\U00B2)", 
       linetype = "Treatment",
       color = "Treatment",
       shape = "Treatment",
       size = "Treatment",
       group = "Treatment") +
  transparent_theme +
  coord_cartesian(ylim = c(-10, 1650))

ggsave(
  plot = overall_kelp_predictions,
  filename = here("figures", "talk-figures",
                  paste0("kelp-biomass_", today(), ".png")),
  width = 12,
  height = 6,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)

ggsave(
  plot = recovery_kelp_only,
  filename = here("figures", "talk-figures",
                  paste0("kelp-biomass-recovery_", today(), ".png")),
  width = 8,
  height = 8,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)

# ⟞ b. delta biomass  -----------------------------------------------------

# This section includes code to create panel b from Figure 2. This figure 
# depicts the "deltas", or difference in predicted biomass between reference 
# and removal plots throughout the experimental removal and recovery period. 
# The data for the lines is taken from the `kelp_models` object created in
# section 2.

delta_kelp_predictions <- ggplot() +
  model_predictions_background +
  geom_point(data = delta_continual,
             aes(x = time_since_end, 
                 y = delta_continual), 
             shape = 2, 
             alpha = 0.7,
             size = 0.75,
             color = "grey") +
  
  # delta biomass
  geom_line(data = pluck(kelp_models, 7, 1), 
            aes(x = x, 
                y = delta), 
            linewidth = 1,
            color = "white") +
  geom_line(data = pluck(kelp_models, 7, 2), 
            aes(x = x, 
                y = delta), 
            linewidth = 1,
            color = "white") +
  
  delta_aesthetics +
  model_predictions_theme +
  transparent_theme + 
  scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1800, 1000))

delta_kelp_predictions

ggsave(
  plot = delta_kelp_predictions,
  filename = here("figures", "talk-figures",
                  paste0("delta-kelp-biomass_", today(), ".png")),
  width = 12,
  height = 6,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)


# ⟞ c. intercepts ---------------------------------------------------------

intercepts <- kelp_models |> 
  mutate(intercepts = map(
    kelp_model,
    ~ ggpredict(.x, terms = c("treatment"))
  )) |> 
  select(exp_dates, intercepts) |> 
  unnest(cols = c(intercepts)) |> 
  rename(treatment = x) 

kelp_intercepts <- ggplot() +
  geom_point(data = continual_long |> filter(treatment == "removal"),
             aes(x = exp_dates,
                 y = kelp_biomass,
                 color = treatment), 
             shape = 21,
             position = position_jitter(
               width = 0.2,
               height = 0,
               seed = 666
             )) +
  geom_pointrange(data = intercepts |> filter(treatment == "removal"),
                  aes(x = exp_dates,
                      y = predicted,
                      ymin = conf.low,
                      ymax = conf.high),
                  color = "white") +
  scale_color_manual(values = c(reference = reference_col, 
                                removal = removal_col)) +
  scale_x_discrete(labels = c("during" = "During removal",
                              "after" = "After removal")) +
  labs(y = "Giant kelp biomass (dry g/m\U00B2)") +
  model_predictions_theme +
  transparent_theme +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave(
  plot = kelp_intercepts,
  filename = here("figures", 
                  "talk-figures",
                  paste0("intercepts-kelp-biomass_", today(), ".png")),
  width = 6,
  height = 10,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)


# ⟞ d. example plot -------------------------------------------------------

example_plot <- ggplot() +
model_predictions_background +
  
  # raw data 
  geom_point(data = continual_long, 
             aes(x = time_since_end, 
                 y = kelp_biomass, 
                 color = treatment), 
             alpha = 0) +
  # theming
  model_predictions_theme +
  model_predictions_aesthetics +
  transparent_theme +
  coord_cartesian(ylim = c(-10, 1650))

example_plot

ggsave(
  plot = example_plot,
  filename = here("figures", 
                  "talk-figures",
                  paste0("example-model-predictions_", today(), ".png")),
  width = 12,
  height = 6,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------- 4. timeseries plot ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section contains code to create timeseries of kelp biomass for each
# site, which is Figure 1.

kelp_biomass_timeseries <- continual_long %>% 
  # create a column for the strip label
  mutate(strip = case_when(
    site == "aque" ~ paste("(a) ", site_full, sep = ""),
    site == "napl" ~ paste("(b) ", site_full, sep = ""),
    site == "mohk" ~ paste("(c) ", site_full, sep = ""),
    site == "carp" ~ paste("(d) ", site_full, sep = "")
  )) %>% 
  ggplot(aes(x = time_since_end,
             y = kelp_biomass,
             color = treatment,
             group = treatment,
             linetype = treatment)) +
  model_predictions_background +
  geom_line(alpha = 0.9,
            linewidth = 0.5) +
  model_predictions_aesthetics +
  raw_biomass_plot_theme +
  labs(y = "Giant kelp biomass (dry g/m\U00B2)") +
  transparent_theme +
  theme(legend.position = "none",
        legend.position.inside = c(0.9, 0.95),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.text = element_text(color = "white")) +
  facet_wrap(~strip, scales = "free_y", nrow = 4)

kelp_biomass_timeseries

ggsave(here::here("figures", "talk-figures",
                  paste("fig-S1_", today(), ".png", sep = "")),
       plot = kelp_biomass_timeseries,
       height = 12, width = 10, units = "cm",
       dpi = 300,
       bg = "transparent")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------- 5. coefficient of variation ----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to calculate the coefficient of variation (CV), 
# which represents temporal variation in giant kelp biomass. 

# We then asked if temporal variation differed between removal and reference 
# plots in the recovery period (section b). If so, this would indicate that 
# giant kelp did not recover its temporal variation in recovery.

# We also tested whether temporal variation differed in the reference plots
# between removal and recovery periods (section c). This would indicate that 
# temporal variation changed in reference plots.

# We found that temporal variation did not differ between removal and reference
# plots in the recovery period, suggesting that giant kelp recovered its 
# temporal variation. We also found that temporal variation did not differ in
# reference plots between removal and recovery periods, suggesting that
# temporal variation in the reference plots remained consistent between the two
# study periods.

# ⟞ a. wrangling ----------------------------------------------------------

variation_site <- continual_long %>% 
  select(exp_dates, treatment, site, kelp_biomass) %>% 
  group_by(exp_dates, treatment, site) %>% 
  # calculating CV
  summarize(mean = mean(kelp_biomass, na.rm = TRUE),
            variation = sd(kelp_biomass, na.rm = TRUE)/
              mean(kelp_biomass, na.rm = TRUE)) %>% 
  ungroup()

# ⟞ b. variation in recovery period ---------------------------------------

# testing differences in variation using t-test
t.test(variation ~ treatment, 
       data = variation_site %>% filter(exp_dates == "after"))

# panel A
recovery_variation <- ggplot(variation_site %>% filter(exp_dates == "after"),
                             aes(x = treatment,
                                 y = variation,
                                 color = treatment)) +
  geom_point(position = position_jitter(width = 0.1, seed = 666),
             alpha = 0.8, shape = 21) +
  stat_summary(fun.data = mean_se, 
               geom = "pointrange") +
  scale_color_manual(values = c(reference = reference_col, 
                                removal = removal_col)) +
  scale_x_discrete(labels = c(reference = "Reference", 
                              removal = "Removal")) +
  labs(x = "Treatment",
       y = "Coefficient of variation",
       title = "a) Comparison of treatments in \nrecovery period") +
  theme_bw() +
  theme(legend.position = "none") +
  transparent_theme 


# ⟞ c. variation in reference plots ---------------------------------------

# testing differences in variation using t-test
t.test(variation ~ exp_dates, 
       data = variation_site %>% filter(treatment == "reference"))

reference_variation <- ggplot(variation_site %>% filter(treatment == "reference"),
                             aes(x = exp_dates,
                                 y = variation)) +
  geom_point(position = position_jitter(width = 0.1, seed = 666),
             alpha = 0.8, shape = 21,
             color = reference_col) +
  stat_summary(fun.data = mean_se, 
               geom = "pointrange",
               color = reference_col) +
  scale_x_discrete(labels = c(during = "Removal period", 
                              after = "Recovery period")) +
  labs(x = "Time period",
       y = "Coefficient of variation",
       title = "b) Comparison of reference plots \nbetween time periods") +
  theme_bw() +
  transparent_theme

# ⟞ d. saving outputs -----------------------------------------------------

variation_plots <- recovery_variation + reference_variation & 
  theme(plot.background = element_rect(fill='transparent',
                                       color = "transparent"),
        legend.background = element_rect(fill = 'transparent')
  )

ggsave(here::here("figures", "talk-figures",
                  paste("cov_plot_", today(), ".png", sep = "")),
       plot = variation_plots,
       height = 7, width = 16, units = "cm",
       dpi = 300,
       bg = "transparent")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------- 6. means plots ----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create conditional means plots to visualize
# predicted kelp biomass in reference and removal plots in the recovery period.
# There are two panels: one panel showing the difference between reference and
# removal plots when time since end = 0 (i.e. at the beginning of the recovery
# period), and another panel showing time since end = 4 (in the most recent
# year of the recovery period).

predicted_kelp_after_0 <- pluck(kelp_models, 3, 2) %>% 
  ggpredict(terms = c("time_since_end [0]", "treatment"), 
            type = "fixed")

predicted_kelp_after_both <- pluck(kelp_models, 3, 2) %>% 
  ggpredict(terms = c("time_since_end [4]", "treatment"), 
            type = "fixed") %>% 
  bind_rows(predicted_kelp_after_0) %>% 
  rename("treatment" = group, 
         "kelp_biomass" = predicted,
         "time_since_end" = x)

means <- continual_long %>% 
  filter(time_since_end == 0 | time_since_end == 4) %>% 
  ggplot(aes(x = treatment, 
             y = kelp_biomass, 
             color = treatment)) +
  geom_point(position = position_jitter(width = 0.1, seed = 1),
             shape = 21, alpha = 0.8, size = 1) +
  geom_pointrange(data = predicted_kelp_after_both,
                  aes(ymin = conf.low, 
                      ymax = conf.high)) +
  scale_color_manual(values = c(reference = reference_col, 
                                removal = removal_col)) +
  scale_x_discrete(labels = c(reference = "Reference", 
                              removal = "Removal")) +
  labs(x = "Treatment",
       y = "Giant kelp biomass (dry g/m\U00B2)") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.text = element_text(hjust = 0, size = 10, color = "white"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  transparent_theme + 
  facet_wrap(~time_since_end, 
             labeller = labeller(
               time_since_end = c("0" = "(a) Time since end = 0", 
                                  "4" = "(b) Time since end = 4")))

# ggsave(here::here("figures", "talk-figures",
#                   paste("kelp_means_plot_", today(), ".png", sep = "")),
#        plot = means,
#        bg = "transparent",
#        height = 7, width = 14, units = "cm",
#        dpi = 200)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 7. tables --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create summary tables for models from section
# 2. These tables are then joined with tables created in the 
# `02a-community_recovery.R` script to create table _____ in the manuscript. It
# relies on the `model_summary_fxn()` from the `00-set_up.R` script.

kelp_model_summaries <- kelp_models %>% 
  select(exp_dates, kelp_model) %>% 
  mutate(model_summary = map(
    kelp_model,
    ~ model_summary_fxn(.x) %>% 
      mutate(group = "kelp") %>% 
      relocate(group, .before = term)
  )) 

kelp_model_summaries_combined <- bind_cols(
  pluck(kelp_model_summaries, 3, 1),
  pluck(kelp_model_summaries, 3, 2)
) %>% 
  select(!group...7) %>% 
  rename("group" = "group...1",
         "term_removal" = "term...2",
         "estimate_removal" = "estimate...3",
         "p.value_removal" = "p.value...4",
         "ci_interval_removal" = "ci_interval...5",
         "signif_removal" = "signif...6",
         "term_recovery" = "term...8",
         "estimate_recovery" = "estimate...9",
         "p.value_recovery" = "p.value...10",
         "ci_interval_recovery" = "ci_interval...11",
         "signif_recovery" = "signif...12")

