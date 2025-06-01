
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. set up -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 1. linear models ----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create a nested dataframe to wrangle data and 
# construct linear models examining how delta kelp (the difference in biomass
# between the removal and reference plots) predicts delta biomass for 
# understory algae and sessile invertebrates in the recovery period. The 
# purpose of this analysis is to explore how change in kelp influences change
# in the biomass of the associated community.

# Similarly to model construction in previous sections, I constructed models
# and simulated residuals within the nested data frame, but extracted the 
# residuals to visualize outside of the nested data frame. The nested data 
# frame relies on the `models` object created in the`02a-community_recovery.R` 
# script and the `delta_continual` object created in the `01a-kelp_recovery.R` 
# script.

# ⟞ a. model fitting ------------------------------------------------------

delta_biomass <- bind_rows(
  # "delta_biomass" column in the `models` object
  pluck(models, 11, 1) %>% mutate(group = "algae"),
  pluck(models, 11, 2) %>% mutate(group = "sessile inverts")
) %>% 
  # nest the data frame and select columns of interest
  nest(data = everything(), .by = group) %>% 
  mutate(data = map(
    data, 
    ~ select(.x,
             sample_ID_short, reference, removal, delta_continual) 
  )) %>% 
  # join data frames with `delta_continual` (giant kelp deltas)
  mutate(data = map(
    data,
    ~ left_join(.x, delta_continual, by = "sample_ID_short") %>% 
      rename(delta_group = delta_continual.x,
             delta_kelp = delta_continual.y) %>% 
      filter(exp_dates == "after")
  )) %>% 
  # construct linear model
  mutate(model = map(
    data,
    ~ lmer(
      delta_group ~ delta_kelp + (1|site) + (1|year),
      data = .x
    )
  )) %>% 
  # simulate residuals
  mutate(residuals = map(
    model,
    ~ simulateResiduals(.x)
  )) 

# ⟞ b. model diagnostics --------------------------------------------------

# for understory algae
plot(pluck(delta_biomass, 4, 1))

# for sessile invertebrates
plot(pluck(delta_biomass, 4, 2))

# ⟞ c. R2 values ----------------------------------------------------------

# for understory algae
r.squaredGLMM(pluck(delta_biomass, 3, 1))

# for sessile inverts
r.squaredGLMM(pluck(delta_biomass, 3, 2))

# understory algae summary
summary(pluck(delta_biomass, 3, 1))

# sessile invert summary
summary(pluck(delta_biomass, 3, 2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 2. model visualizations -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create panels a and b from Figure 4. The model
# visualizations for understory algae and sessile invertebrates are created in
# separate sections. Because delta kelp did not have a significant effect on
# delta sessile inverts, the visualization for sessile inverts does not show
# model predictions.

# ⟞ a. model predictions --------------------------------------------------

algae_predictions <- ggpredict(
  pluck(delta_biomass, 3, 1),
  terms = c("delta_kelp"),
  type = "fixed"
)

# ⟞ b. visualizations -----------------------------------------------------

# ⟞ ⟞ i. understory algae -------------------------------------------------

algae_vs_kelp_plot <- pluck(delta_biomass, 2, 1) %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = delta_kelp, y = delta_group)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_point(size = 1, shape = 21, alpha = 0.4, color = under_col) + 
  geom_ribbon(data = algae_predictions, 
              aes(x = x, 
                  y = predicted, 
                  ymin = conf.low, 
                  ymax = conf.high),
              alpha = 0.1,
              fill = "grey") +
  geom_line(data = algae_predictions, 
            aes(x = x, y = predicted), 
            linewidth = 1,
            color = under_col) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), 
                     minor_breaks = seq(-2000, 2000, by = 500)) +
  scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)") +
  model_predictions_theme +
  transparent_theme

example_plot <- pluck(delta_biomass, 2, 1) %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = delta_kelp, y = delta_group)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_point(size = 1, shape = 21, alpha = 0, color = under_col) + 
  geom_ribbon(data = algae_predictions, 
              aes(x = x, 
                  y = predicted, 
                  ymin = conf.low, 
                  ymax = conf.high),
              alpha = 0,
              fill = "grey") +
  geom_line(data = algae_predictions, 
            aes(x = x, y = predicted), 
            linewidth = 1,
            color = under_col,
            alpha = 0) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), 
                     minor_breaks = seq(-2000, 2000, by = 500)) +
  scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)") +
  model_predictions_theme +
  transparent_theme

ggsave(
  filename = here("figures", "talk-figures",
                  paste0("algae-vs-kelp_", today(), ".png")),
  plot = algae_vs_kelp_plot,
  bg = "transparent",
  width = 10,
  height = 10,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = here("figures", "talk-figures",
                  paste0("example-algae-kelp-plot_", today(), ".png")),
  plot = example_plot,
  bg = "transparent",
  width = 10,
  height = 10,
  units = "cm",
  dpi = 300
)

# ⟞ ⟞ ii. sessile invertebrates -------------------------------------------

epi_vs_kelp_plot <- pluck(delta_biomass, 2, 2) %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = delta_kelp, y = delta_group)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_point(size = 1, shape = 21, alpha = 0.4, color = "white") + 
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 sessile invertebrate biomass\n(removal - reference, dry g/m\U00B2)") +
  model_predictions_theme +
  transparent_theme

# ⟞ c. saving outputs -----------------------------------------------------

group_vs_kelp <- plot_grid(algae_vs_kelp_plot, epi_vs_kelp_plot, ncol = 2)

ggsave(here::here("figures", "talk-figures",
                  paste("fig-4_", today(), ".png", sep = "")),
       plot = group_vs_kelp,
       bg = "transparent",
       height = 8, width = 16, units = "cm",
       dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 3. tables --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. wrangling ----------------------------------------------------------

group_vs_kelp_table_fxn <- function(model) {
  tidy(model, conf.int = TRUE) %>% 
    filter(effect == "fixed") %>% 
    select(term, estimate, p.value, conf.low, conf.high) %>% 
    # create a new column that indicates whether an effect is significant
    mutate(signif = case_when(
      p.value <= 0.05 ~ "yes",
      TRUE ~ "no"
    )) %>% 
    # create a p-value column that converts very small values to < 0.001
    # and rounds all other values to relevant digits
    mutate(p.value = case_when(
      between(p.value, 0, 0.001) ~ "<0.001",
      between(p.value, 0.001, 0.01) ~ as.character(round(p.value, digits = 3)),
      between(p.value, 0.01, 1) ~ as.character(round(p.value, digits = 2))
    )) %>%
    # round other numeric values to two digits
    mutate(across(where(is.numeric), ~ round(., digits = 2))) %>%
    # create a confidence interval column
    unite(ci_interval, conf.low, conf.high, sep = ", ") %>% 
    mutate(term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "delta_kelp" ~ "\U0394 giant kelp"
    ))
  
}

algae_table <- group_vs_kelp_table_fxn(pluck(delta_biomass, 3, 1)) %>% 
  mutate(group = "Understory macroalgae") %>% 
  relocate(group, .before = term)

epi_table <- group_vs_kelp_table_fxn(pluck(delta_biomass, 3, 2)) %>% 
  mutate(group = "Sessile invertebrates") %>% 
  relocate(group, .before = term)

# ⟞ b. table creation -----------------------------------------------------

group_vs_kelp_table <- bind_rows(
  algae_table,
  epi_table) %>% 
  
  # change object into a flextable, select columns to display
  flextable(col_keys = c("group",
                         "term",
                         "estimate",
                         "p.value",
                         "ci_interval")) %>% 
  
  # change column names
  set_header_labels("group" = " ",
                    "term" = "Term",
                    "estimate" = "Estimate",
                    "p.value" = "p-value",
                    "ci_interval" = "95% CI") %>% 
  # add footnote for 95% CI
  footnote(
    i = 1, 
    j = 5,
    ref_symbols = "1",
    value = as_paragraph("Confidence interval"),
    part = "header"
  ) %>% 
  
  # bold p-values if significant
  style(i = ~ signif == "yes",
        j = "p.value",
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  
  # merge group cells to create a grouping column
  merge_v(j = ~ group) %>% 
  valign(j = ~ group,
         i = NULL,
         valign = "top") %>% 
  
  # final formatting
  autofit %>% 
  fit_to_width(5) %>% 
  font(fontname = "Times New Roman",
       part = "all")

# ⟞ c. saving output ------------------------------------------------------

# group_vs_kelp_table %>%
#   save_as_docx(path = here::here("tables",
#                                  "ms-tables",
#                                  paste0("tbl-S5_", today(), ".docx")))
