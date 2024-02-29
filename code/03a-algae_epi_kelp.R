##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

##########################################################################-
# 1. algae ----------------------------------------------------------------
##########################################################################-

# ⟞ a. correlation --------------------------------------------------------

delta_algae_after <- delta_algae_continual %>% 
  filter(exp_dates == "after")

cor.test(delta_algae_after$delta_continual_algae, delta_algae_after$delta_continual,
         method = "pearson")

# ⟞ b. linear models ------------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# models
# lm_delta_algae_kelp_after_m1 <- lm(
#   delta_continual_algae ~ delta_continual, 
#   data = delta_algae_continual %>% filter(exp_dates == "after"),
#   na.action = na.omit
# )
lm_delta_algae_kelp_after_m2_outlier <- lmer(
  delta_continual_algae ~ delta_continual + (1|site) + (1|year),
  data = delta_algae_continual %>% filter(exp_dates == "after"))

check_outliers(lm_delta_algae_kelp_after_m2_outlier, c("cook"))
# outlier: napl_2023-05-18_Q2 because of a lot of PTCA biomass in control plot

lm_delta_algae_kelp_after_m2 <- lmer(
  delta_continual_algae ~ delta_continual + (1|site) + (1|year),
  data = delta_algae_continual %>% filter(exp_dates == "after" & sample_ID != "napl_2023-05-18_Q2"))

ggpredict(lm_delta_algae_kelp_after_m2, terms = c("delta_continual")) %>% plot(show_data = TRUE)

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
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_vline(aes(xintercept = 0), lty = 2) +
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
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_vline(aes(xintercept = 0), lty = 2) +
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
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_vline(aes(xintercept = 0), lty = 2) +
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
