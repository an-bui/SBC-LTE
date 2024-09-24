
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
# residuals to visualize outside of the nested data frame. This process 
# revealed two outliers for understory algae: 1) Naples on 2023-05-18 likely 
# because of high Pterygophora californica biomass in the reference plot, and 
# 2) Mohawk on 2019-11-19 likely because of high giant kelp biomass in the 
# removal plot. I excluded these observations from the data set, reran the 
# model, and checked the residuals to verify that these outliers improved model
# fit. In the manuscript, we present the model without outliers in the main 
# text, and provide a visualization of the model with both outliers in the 
# supplemental material.

# The nested data frame relies on the `models` object created in the 
# `02a-community_recovery.R` script and the `delta_continual` object created in
# the `01a-kelp_recovery.R` script.

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

# nothing from naples 2014-11-14 in continual removal plot
# nothing from mohk 2010-06-14 in continual removal plot

# ⟞ b. model diagnostics --------------------------------------------------

# for understory algae
plot(pluck(delta_biomass, 4, 1))

# for sessile invertebrates
plot(pluck(delta_biomass, 4, 2))

# ⟞ c. reworking understory algae model -----------------------------------

# check for outliers using `performance::check_outliers`
check_outliers(pluck(delta_biomass, 3, 1), c("cook"))
# case 96: napl_2023-05-18_Q2

algae_no_outlier_v1 <- lmer(
  delta_group ~ delta_kelp + (1|site) + (1|year),
  data = pluck(delta_biomass, 2, 1) %>% 
    filter(sample_ID_short %!in% c("napl_2023-05-18"))
)

plot(simulateResiduals(algae_no_outlier_v1)) 

check_outliers(algae_no_outlier_v1, c("cook"))
# case 56: mohk_2019_11-19

algae_no_outlier_v2 <- lmer(
  delta_group ~ delta_kelp + (1|site) + (1|year),
  data = pluck(delta_biomass, 2, 1) %>% 
    filter(sample_ID_short %!in% c("napl_2023-05-18", "mohk_2019-11-19"))
)

plot(simulateResiduals(algae_no_outlier_v2)) 
# residuals look much better

check_outliers(algae_no_outlier_v2, c("cook"))
# no outliers

# ⟞ d. R2 values ----------------------------------------------------------

# for understory algae
r.squaredGLMM(algae_no_outlier_v2)

# for sessile inverts
r.squaredGLMM(pluck(delta_biomass, 3, 2))

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
  algae_no_outlier_v2,
  terms = c("delta_kelp"),
  type = "fixed"
)

# ⟞ b. visualizations -----------------------------------------------------

# ⟞ ⟞ i. understory algae -------------------------------------------------

algae_vs_kelp_plot <- pluck(delta_biomass, 2, 1) %>% 
  filter(exp_dates == "after" & 
           sample_ID_short %!in% c("napl_2023-05-18", "mohk_2019-11-19")) %>% 
  ggplot(aes(x = delta_kelp, y = delta_group)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_point(size = 1, shape = 21, alpha = 0.4, color = under_col) + 
  geom_ribbon(data = algae_predictions, 
              aes(x = x, 
                  y = predicted, 
                  ymin = conf.low, 
                  ymax = conf.high), 
              alpha = 0.1) +
  geom_line(data = algae_predictions, 
            aes(x = x, y = predicted), 
            linewidth = 1,
            color = under_col) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), 
                     minor_breaks = seq(-2000, 2000, by = 500)) +
  scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)", 
       title = "(a) Understory macroalgae") +
  annotate("text", x = -1100, y = -200,
           label = "conditional R\U00B2 = 0.51\nmarginal R\U00B2 = 0.42\np < 0.001",
           size = 1.5) +
  model_predictions_theme

# ⟞ ⟞ ii. sessile invertebrates -------------------------------------------

epi_vs_kelp_plot <- pluck(delta_biomass, 2, 2) %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = delta_kelp, y = delta_group)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_point(size = 1, shape = 21, alpha = 0.4) + 
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 sessile invertebrate biomass\n(removal - reference, dry g/m\U00B2)", 
       title = "(b) Sessile invertebrates") +
  model_predictions_theme

# ⟞ c. saving outputs -----------------------------------------------------

group_vs_kelp <- plot_grid(algae_vs_kelp_plot, epi_vs_kelp_plot, ncol = 2)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-4_", today(), ".jpg", sep = "")),
#        plot = group_vs_kelp,
#        height = 8, width = 14, units = "cm",
#        dpi = 300)

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

algae_table <- group_vs_kelp_table_fxn(algae_no_outlier_v2) %>% 
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







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------- OLD CODE BELOW HERE --------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# lm_delta_algae_kelp_after_m3 <- lmer(
#   delta_continual_algae ~ delta_continual + (1|site), 
#   data = delta_algae_continual %>% filter(exp_dates == "after")
# )
# lm_algae_kelp_after_m1 <- lmer(
#   continual_algae ~ continual + (1|year) + (1|site),
#   data = delta_algae_continual %>% filter(exp_dates == "after") %>% filter(!(sample_ID %in% c("mohk_2019-11-19_Q4", "mohk_2021-05-13_Q2", "napl_2023-05-18_Q2")))
# )
# lm_algae_kelp_after_m2 <- lmer(
#   continual_algae ~ continual*time_since_end + (1|site),
#   data = delta_algae_continual %>% filter(exp_dates == "after")
# )

# diagnostics
# check_model(lm_delta_algae_kelp_after_m1)
# simulateResiduals(lm_delta_algae_kelp_after_m1, plot = TRUE)

simulateResiduals(lm_delta_algae_kelp_after_m2, plot = TRUE)

outlier_check <- check_outliers(lm_delta_algae_kelp_after_m2_outlier) %>% plot()
simulateResiduals(lm_delta_algae_kelp_after_m2_outlier, plot = TRUE)
# check_model(lm_delta_algae_kelp_after_m3)

# check_model(lm_algae_kelp_after_m1)
# plot(simulateResiduals(lm_algae_kelp_after_m1))

# Rsquared
# r.squaredGLMM(lm_delta_algae_kelp_after_m1)
r.squaredGLMM(lm_delta_algae_kelp_after_m2)
# r.squaredGLMM(lm_delta_algae_kelp_after_m3)

# r.squaredGLMM(lm_algae_kelp_after_m1)
# r.squaredGLMM(lm_algae_kelp_after_m2)

# summaries
# summary(lm_delta_algae_kelp_after_m1)
summary(lm_delta_algae_kelp_after_m2)
# summary(lm_delta_algae_kelp_after_m3)

lm_delta_algae_kelp_after_summary <- lm_delta_algae_kelp_after_m2 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Slope**"
  ) 

# AIC comparison
# AICc(lm_delta_algae_kelp_after_m1, lm_delta_algae_kelp_after_m2, lm_delta_algae_kelp_after_m3) %>% 
#   arrange(AICc)
# AICc(lm_algae_kelp_after_m1, lm_algae_kelp_after_m2)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# raw biomass
# predicted_algae_kelp_after <- ggpredict(lm_algae_kelp_after_m1, terms = ~ continual, type = "fixed")

# deltas
predicted_delta_algae_vs_kelp <- ggpredict(lm_delta_algae_kelp_after_m2, terms = ~ delta_continual, type = "fixed")

predicted_delta_algae_vs_kelp_outlier <- ggpredict(lm_delta_algae_kelp_after_m2_outlier, terms = ~ delta_continual, type = "fixed")

# ⟞ c. figures ------------------------------------------------------------

# ⟞ ⟞ i. raw biomass ------------------------------------------------------

# algae_kelp_after_plot <- ggplot(data = delta_algae_continual %>% 
#                                   filter(exp_dates == "after"), 
#                                 aes(x = continual, y = continual_algae)) +
#   geom_point(shape = 21) +
#   geom_line(data = predicted_algae_kelp_after, aes(x = x, y = predicted), lty = 1) +
#   geom_ribbon(data = predicted_algae_kelp_after, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2)
# 
# algae_kelp_after_plot

# ⟞ ⟞ ii. deltas ----------------------------------------------------------

delta_algae_vs_kelp_lm <- delta_algae_continual %>% 
  filter(exp_dates == "after" & sample_ID != "napl_2023-05-18_Q2") %>% 
  # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
  ggplot(aes(x = delta_continual, y = delta_continual_algae)) +
  geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.5) +
  geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.5) +
  geom_point(size = 1, shape = 5, alpha = 0.4, color = under_col) + 
  geom_ribbon(data = predicted_delta_algae_vs_kelp, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_line(data = predicted_delta_algae_vs_kelp, aes(x = x, y = predicted), 
            linewidth = 1,
            color = under_col) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
  scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)", 
       title = "(a) Understory macroalgae") +
  annotate("text", x = -1100, y = -200,
           label = "conditional R\U00B2 = 0.44\nmarginal R\U00B2 = 0.26\np < 0.001",
           size = 1.5) +
  theme_bw() + 
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        plot.title = element_text(size = 8),
        plot.title.position = "plot",
        panel.grid = element_blank()) 
delta_algae_vs_kelp_lm

delta_algae_vs_kelp_lm_outlier <- delta_algae_continual %>% 
  filter(exp_dates == "after") %>% 
  # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
  ggplot(aes(x = delta_continual, y = delta_continual_algae)) +
  geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.5) +
  geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.5) +
  geom_point(size = 1, shape = 5, alpha = 0.4, color = under_col) + 
  geom_ribbon(data = predicted_delta_algae_vs_kelp_outlier, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_line(data = predicted_delta_algae_vs_kelp_outlier, aes(x = x, y = predicted), 
            linewidth = 1,
            color = under_col) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
  scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)", 
       title = "(a) Understory macroalgae") +
  theme_bw() + 
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        plot.title = element_text(size = 8),
        plot.title.position = "plot",
        panel.grid = element_blank()) 
delta_algae_vs_kelp_lm_outlier

# delta_algae_vs_kelp_pearson <- delta_algae_continual %>% 
#   mutate(exp_dates = case_when(
#     exp_dates == "during" ~ "During removal",
#     exp_dates == "after" ~ "Recovery period"
#   )) %>% 
#   ggplot(aes(x = delta_continual, y = delta_continual_algae, linetype = exp_dates, color = exp_dates, fill = exp_dates)) +
#   geom_hline(aes(yintercept = 0), lty = 2) +
#   geom_vline(aes(xintercept = 0), lty = 2) +
#   geom_point(aes(shape = exp_dates, fill = exp_dates), size = 5, shape = 21, color = "#000000") +
#   stat_cor(method = "pearson", cor.coef.name = "rho", inherit.aes = TRUE) +
#   geom_smooth(aes(color = exp_dates), method = "lm", color = "black", linewidth = 3, se = FALSE) +
#   scale_linetype_manual(values = c("During removal" = 2, "Recovery period" = 1)) +
#   scale_color_manual(values = c("During removal" = "grey", "Recovery period" = under_col)) +
#   scale_fill_manual(values = c("During removal" = "#FFFFFF", "Recovery period" = under_col)) +
#   scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
#   labs(x = "\U0394 kelp biomass (treatment - control)",
#        y = "\U0394 understory macroalgae biomass (treatment - control)") +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 18),
#         plot.title = element_text(size = 18),
#         axis.text = element_text(size = 16),
#         legend.position = c(0.83, 0.93),
#         legend.text = element_text(size = 18),
#         legend.title = element_blank())
# delta_algae_vs_kelp_pearson

##########################################################################-
# 2. epilithic inverts ----------------------------------------------------
##########################################################################-

# ⟞ a. correlation --------------------------------------------------------

delta_epi_after <- delta_epi_continual %>% 
  filter(exp_dates == "after")

cor.test(delta_epi_after$delta_continual_epi, delta_epi_after$delta_continual,
         method = "pearson")

# ⟞ b. linear models ------------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# lm_delta_epi_kelp_after_m1 <- lm(
#   delta_continual_epi ~ delta_continual, 
#   data = delta_epi_continual %>% filter(exp_dates == "after"),
#   na.action = na.omit
# )
lm_delta_epi_kelp_after_m2 <- lmer(
  delta_continual_epi ~ delta_continual + (1|year) + (1|site), 
  data = delta_epi_continual %>% filter(exp_dates == "after")
)
# lm_delta_epi_kelp_after_m3 <- lmer(
#   delta_continual_epi ~ delta_continual + (1|site), 
#   data = delta_epi_continual %>% filter(exp_dates == "after")
# )

# lm_epi_kelp_after_m1 <- lmer(
#   continual_epi ~ continual + (1|year) + (1|site),
#   data = delta_epi_continual %>% filter(exp_dates == "after")
# )

# lm_epi_kelp_after_m2 <- lmer(
#   continual_epi ~ continual*time_since_end + (1|site),
#   data = delta_epi_continual %>% filter(exp_dates == "after")
# )

# diagnostics
# simulateResiduals(lm_delta_epi_kelp_after_m1, plot = T)
# check_model(lm_delta_epi_kelp_after_m1)

simulateResiduals(lm_delta_epi_kelp_after_m2, plot = T)

# simulateResiduals(lm_delta_epi_kelp_after_m3, plot = T)
# check_model(lm_delta_epi_kelp_after_m3)

# check_model(lm_epi_kelp_after_m1)
# check_model(lm_epi_kelp_after_m2)

# Rsquared
# r.squaredGLMM(lm_delta_epi_kelp_after_m1)
r.squaredGLMM(lm_delta_epi_kelp_after_m2)
# r.squaredGLMM(lm_delta_epi_kelp_after_m3)

# r.squaredGLMM(lm_epi_kelp_after_m1)
# r.squaredGLMM(lm_epi_kelp_after_m2)

# summaries
# summary(lm_delta_epi_kelp_after_m1)
summary(lm_delta_epi_kelp_after_m2)
# summary(lm_delta_epi_kelp_after_m3)

lm_delta_epi_kelp_after_summary <- lm_delta_epi_kelp_after_m2 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Slope**"
  ) 
lm_delta_epi_kelp_after_summary

# AICc
# AICc(lm_delta_epi_kelp_after_m1, lm_delta_epi_kelp_after_m2, lm_delta_epi_kelp_after_m3) %>% 
#   arrange(AICc)

AICc(lm_epi_kelp_after_m1, lm_epi_kelp_after_m2) # same?

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# raw biomass
# predicted_epi_vs_kelp <- ggpredict(lm_epi_kelp_after_m1, terms = ~ continual, type = "fixed")

# deltas
predicted_delta_epi_vs_kelp <- ggpredict(lm_delta_epi_kelp_after_m2, terms = ~ delta_continual, type = "fixed")

# ⟞ c. figures ------------------------------------------------------------

# ⟞ ⟞ i. raw biomass ------------------------------------------------------

# delta_epi_continual %>% 
#   filter(exp_dates == "during") %>% 
#   ggplot(aes(x = continual, y = continual_epi)) +
#   geom_point() +
#   geom_line(data = predicted_epi_vs_kelp, aes(x = x, y = predicted))

# ⟞ ⟞ ii. deltas ----------------------------------------------------------

delta_epi_vs_kelp_lm <- delta_epi_continual %>% 
  filter(exp_dates == "after") %>% 
  # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
  ggplot(aes(x = delta_continual, y = delta_continual_epi)) +
  geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.5) +
  geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.5) +
  geom_point(size = 1, shape = 5, alpha = 0.4) + 
  # geom_ribbon(data = predicted_delta_epi_vs_kelp, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  # geom_line(data = predicted_delta_epi_vs_kelp, aes(x = x, y = predicted), size = 2) +
  # scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
  labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
       y = "\U0394 sessile invertebrate biomass\n(removal - reference, dry g/m\U00B2)", 
       title = "(b) Sessile invertebrates") +
  theme_bw() + 
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        plot.title = element_text(size = 8),
        plot.title.position = "plot",
        panel.grid = element_blank()) 
delta_epi_vs_kelp_lm

# epi_vs_kelp_pearson <- delta_epi_continual %>% 
#   mutate(exp_dates = case_when(
#     exp_dates == "during" ~ "During removal",
#     exp_dates == "after" ~ "Post-removal"
#   )) %>% 
#   ggplot(aes(x = delta_continual, y = delta_continual_epi, linetype = exp_dates, color = exp_dates, fill = exp_dates)) +
#   stat_cor(aes(color = exp_dates), method = "pearson", cor.coef.name = "rho") +
#   geom_hline(aes(yintercept = 0), lty = 2) +
#   geom_vline(aes(xintercept = 0), lty = 2) +
#   geom_point(size = 5, shape = 25, color = "#000000") + 
#   geom_smooth(aes(linetype = exp_dates), method = "lm", se = FALSE, color = "black", size = 3) +
#   scale_linetype_manual(values = c("During removal" = 2, "Post-removal" = 1)) +
#   scale_fill_manual(values = c("During removal" = "#FFFFFF", "Post-removal" = "#54662C")) +
#   scale_color_manual(values = c("During removal" = "grey", "Post-removal" = "#54662C")) +
#   scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
#   labs(x = "\U0394 kelp biomass (treatment - control)",
#        y = "\U0394 epilithic invertebrate biomass (treatment - control)") +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 18),
#         plot.title = element_text(size = 18),
#         axis.text = element_text(size = 16),
#         legend.position = c(0.83, 0.93),
#         legend.text = element_text(size = 18),
#         legend.title = element_blank())
# epi_vs_kelp_pearson

##########################################################################-
# 3. manuscript tables ----------------------------------------------------
##########################################################################-

lm_vs_kelp_summary_tables <- tbl_stack(
  tbls = list(lm_delta_algae_kelp_after_summary, lm_delta_epi_kelp_after_summary),
  group_header = c("Understory macroalgae", "Sessile invertebrates"),
  quiet = TRUE) %>% 
  as_flex_table() %>% 
  font(fontname = "Times New Roman", part = "all")

# lm_vs_kelp_summary_tables %>%
#   save_as_docx(path = here::here("tables", "ms-tables", paste("tbl-S5_", today(), ".docx", sep = "")))


##########################################################################-
# 4. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. correlation --------------------------------------------------------

algae_vs_kelp_spearman

# ⟞ b. linear models ------------------------------------------------------

group_vs_kelp <- plot_grid(delta_algae_vs_kelp_lm, delta_epi_vs_kelp_lm, ncol = 2)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-4_", today(), ".jpg", sep = "")),
#        plot = group_vs_kelp,
#        height = 6, width = 12, units = "cm",
#        dpi = 300)


# ⟞ c. outlier check ------------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S12_", today(), ".jpg", sep = "")),
#        plot = outlier_check,
#        height = 8, width = 14, units = "cm",
#        dpi = 300)
