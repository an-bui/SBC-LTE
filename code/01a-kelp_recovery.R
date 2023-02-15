
##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "00-set_up.R"))

##########################################################################-
# 1. data frames ----------------------------------------------------------
##########################################################################-

# ⟞ a. annual removal deltas ----------------------------------------------

delta_annual <- biomass %>% 
  filter(sp_code == "MAPY" & treatment %in% c("control", "annual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
  # fill in values for sampling dates where control and annual were surveyed on different days
  mutate(control = case_when(
    site == "aque" & date == "2008-03-07" ~ 106.82000,
    site == "napl" & date == "2012-09-26" ~ 143.22088,
    TRUE ~ control
  )) %>% 
  mutate(delta_annual = annual - control) %>%  
  # missing dates are from AQUE 2008-03-05, NAPL 2012-09-25, NAPL 2008-10-10
  drop_na(delta_annual) %>% 
  mutate(exp_dates = case_when(
    # after for continual removal:
    site == "aque" & date >= aque_after_date_annual ~ "after",
    site == "napl" & date >= napl_after_date_annual ~ "after",
    site == "ivee" & date >= ivee_after_date_annual ~ "after", 
    site == "mohk" & date >= mohk_after_date_annual ~ "after",
    site == "carp" & date >= carp_after_date_annual ~ "after",
    # everything else is "during" the experiment
    TRUE ~ "during"
  ),
  exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
  time_since_columns_annual() %>% 
  kelp_year_column() %>% 
  comparison_column_annual() 

# ⟞ b. continual removal deltas -------------------------------------------

delta_continual <- biomass %>% 
  filter(sp_code == "MAPY" & treatment %in% c("control", "continual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
  mutate(delta_continual = continual - control) %>%  
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for continual removal:
    site == "aque" & date >= aque_after_date_continual ~ "after",
    site == "napl" & date >= napl_after_date_continual ~ "after",
    site == "mohk" & date >= mohk_after_date_continual ~ "after",
    site == "carp" & date >= carp_after_date_continual ~ "after",
    # everything else is "during" the experiment
    TRUE ~ "during"
  ),
  exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
  time_since_columns_continual() %>% 
  kelp_year_column() %>% 
  comparison_column_continual() %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria"))

##########################################################################-
# 2. timeseries plots -----------------------------------------------------
##########################################################################-

# ⟞ a. delta annual plot --------------------------------------------------

delta_annual_plot <- ggplot(delta_annual, aes(x = time_since_end, y = delta_annual)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(col = exp_dates)) +
  geom_smooth(method = "lm", aes(col = exp_dates)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  labs(x = "Time since end of experiment", y = "\U0394 biomass (treatment - control)",
       title = "\U0394 annual") +
  facet_wrap(~site_full, ncol = 1, scales = "free_y")

delta_annual_plot

# ⟞ b. delta continual plot -----------------------------------------------

delta_continual_plot <- ggplot(delta_continual, aes(x = time_since_end, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(col = exp_dates)) +
  geom_smooth(method = "lm", aes(col = exp_dates)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  labs(x = "Time since end of experiment", y = "\U0394 biomass (treatment - control)",
       title = "\U0394 continual") +
  facet_wrap(~site_full, ncol = 1, scales = "free_y")

delta_continual_plot

# ⟞ c. raw biomass through time plot --------------------------------------

delta_continual_sites_raw <- delta_continual %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("A. ", site_full, sep = ""),
    site == "napl" ~ paste("B. ", site_full),
    site == "mohk" ~ paste("C. ", site_full),
    site == "carp" ~ paste("D. ", site_full)
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(x = time_since_end, y = control, col = site), alpha = 0.5, size = 2.5) +
  # control
  geom_point(aes(x = time_since_end, y = control, shape = site), size = 1.5, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual, col = site), size = 2.5) +
  geom_point(aes(x = time_since_end, y = continual, shape = site, col = site), size = 1.5, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16), 
        legend.position = "none", 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 20, hjust = 0),
        panel.grid.minor = element_line(color = "white")) +
  labs(x = "Time since end of experiment (years)", 
       y = expression(Giant~kelp~biomass~(dry~g/m^{"2"}))) +
  facet_wrap(~strip, scales = "free_y")

delta_continual_sites_raw

##########################################################################-
# 3. linear models --------------------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_kelp_during_lmer <- lmer(
  delta_continual ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass)
lm_kelp_during_lme_ar1 <- lme(
  delta_continual ~ time_since_end, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass,
  correlation = corAR1())

# diagnostics
plot(simulateResiduals(lm_kelp_during_lmer))
check_model(lm_kelp_during_lmer)
check_model(lm_kelp_during_lme_ar1)
plot(ACF(lm_kelp_during_lme_ar1))

# Rsquared
MuMIn::r.squaredGLMM(lm_kelp_during_lmer)
MuMIn::r.squaredGLMM(lm_kelp_during_lme_ar1)

# summaries
summary(lm_kelp_during_lmer)
summary(lm_kelp_during_lme_ar1)
lm_kelp_during_summary <- lm_kelp_during_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_kelp_during_summary

# AIC comparison
AICc(lm_kelp_during_lmer, lm_kelp_during_lme_ar1)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_kelp_during_overall <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "fixed")
predicted_kelp_during_aque <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
predicted_kelp_during_napl <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
predicted_kelp_during_mohk <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
predicted_kelp_during_carp <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_kelp_recovery_lmer <- lmer(
  delta_continual ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass)
lm_kelp_recovery_lme_ar1 <- nlme::lme(
  delta_continual ~ time_since_end, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass,
  correlation = corAR1())

# diagnostics
plot(simulateResiduals(lm_kelp_recovery_lmer)) # convergence problems?
check_model(lm_kelp_recovery_lmer)
check_model(lm_kelp_recovery_lme_ar1_m1)
plot(ACF(lm_kelp_recovery_lme_ar1_m1))

# Rsquared
MuMIn::r.squaredGLMM(lm_kelp_recovery_lmer)
MuMIn::r.squaredGLMM(lm_kelp_recovery_lme_ar1)

# summary
summary(lm_kelp_recovery_lmer)
summary(lm_kelp_recovery_lme_ar1)
lm_kelp_recovery_summary <- lm_kelp_recovery_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_kelp_recovery_summary

# AIC comparisons
AICc(lm_kelp_recovery_lmer, lm_kelp_recovery_lme_ar1_m1)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# all sites
predicted_kelp_after_overall <- ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end", type = "fixed")
# predicted line crosses 0 at 3.97
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.97:4.0 by = 0.01]", type = "fixed")
# lower bound crosses 0 at 2.6
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [2.60:2.63 by = 0.001]", type = "fixed")
# upper bound crosses 0 at 6.4
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [6.44:6.45 by = 0.001]", type = "fixed")

# aque: 3.8 years
predicted_kelp_after_aque <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.7:3.8 by = 0.001]", type = "random", condition = c(site = "aque"))

# napl: 3.4 years
predicted_kelp_after_napl <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.3:3.5 by = 0.001]", type = "random", condition = c(site = "napl"))

# mohk: 5.4 years
predicted_kelp_after_mohk <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [5.37:5.41 by = 0.01]", type = "random", condition = c(site = "mohk"))

# carp: 3.3 years
predicted_kelp_after_carp <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.2:3.4 by = 0.01]", type = "random", condition = c(site = "carp"))

# ⟞ c. figures -------------------------------------------------------------

# ⟞ ⟞ i. overall model predictions -----------------------------------------

overall_ms <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual, 
             aes(x = time_since_end, y = delta_continual, fill = site, shape = site), size = 4, alpha = 0.9) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  # new_scale("color") + 
  # overall
  geom_line(data = predicted_kelp_after_overall, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  geom_line(data = predicted_kelp_during_overall, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(-1500, 1000, by = 1000), limits = c(-1800, 900)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        # plot.margin = margin(0, 0, 0, 0),
        legend.position = c(0.12, 0.885)) +
  labs(x = "Time since end of removal (years)", 
       y = expression("\U0394"~giant~kelp~biomass~"(treatment - control, "~dry~g/m^{"2"}~")"), 
       fill = "Site",
       shape = "Site")

overall_ms

# ⟞ ⟞ ii. site level predictions -------------------------------------------

aque <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "aque"), aes(x = time_since_end, y = delta_continual), shape = aque_shape, fill = aque_col, size = 3, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_aque, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_aque, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1000, by = 1000), limits = c(-1800, 1000)) +
  # geom_text(aes(x = -6.6, y = 600), label = "aque", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = aque_full)

napl <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "napl"), aes(x = time_since_end, y = delta_continual), shape = napl_shape, fill = napl_col, size = 3, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_napl, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_napl, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1500, by = 1000), limits = c(-1800, 1500)) +
  # geom_text(aes(x = -6.8, y = 700), label = "napl", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = napl_full)

mohk <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "mohk"), aes(x = time_since_end, y = delta_continual), shape = mohk_shape, fill = mohk_col, size = 3, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_mohk, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_mohk, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1500, by = 1000), limits = c(-1800, 1500)) +
  # geom_text(aes(x = -6.6, y = 900), label = "mohk", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = mohk_full)

carp <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "carp"), aes(x = time_since_end, y = delta_continual), shape = carp_shape, fill = carp_col, size = 3, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_carp, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_carp, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1500, by = 1000), limits = c(-1800, 1500)) +
  # geom_text(aes(x = -6.8, y = 1000), label = "carp", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = carp_full)

plots_together_sites <- (aque + napl) / (mohk + carp) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30))
plots_together_sites


##########################################################################-
# 4. recovery time vs biomass ---------------------------------------------
##########################################################################-

# ⟞ a. data frame ---------------------------------------------------------

# Not sure how to do this reproducibly. But pulled time to recovery from above predicted values with prediction intervals. Then calculated mean kelp biomass in the kelp year of recovery from `delta_continual`. For MOHK, calculated biomass using the most recent kelp year (2021-2022).

mean_kelp_all_sites <- delta_continual %>% 
  filter(exp_dates == "after") %>% 
  group_by(site) %>% 
  # filter(site == "aque" & exp_dates == "after") %>% 
  summarize(mean_control = mean(control),
            sd_control = sd(control),
            se_control = se(control),
            mean_continual = mean(continual),
            sd_continual = sd(continual),
            se_continual = se(continual),
            mean_delta = mean(delta_continual),
            sd_delta = sd(delta_continual),
            se_delta = se(delta_continual)) %>% 
  ungroup()

rec_time <- tribble(
  ~site, ~time_to_recovery, ~pred_int_low, ~pred_int_high,
  "aque",       3.80,          -674.27,         676.81,
  "napl",       3.42,          -672.08,         675.42,
  "mohk",       5.4,           -688.62,         689.42,
  "carp",       3.33,          -672.28,         674.48
) %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  left_join(., mean_kelp_all_sites, by = "site") 

# ⟞ b. figure -------------------------------------------------------------

rec_time_plot <- ggplot(rec_time, aes(x = mean_control, y = time_to_recovery, shape = site, fill = site)) +
  geom_errorbar(aes(xmin = mean_control - se_control, xmax = mean_control + se_control)) +
  geom_point(size = 10) +
  geom_text_repel(aes(label = site_full), seed = 666,
                  size = 7, point.padding = 35, nudge_y = 0.1) +
  scale_shape_manual(values = shape_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  theme_bw() +
  scale_y_continuous(breaks = c(seq(3, 6, by = 1)), limits = c(3, 6)) +
  labs(x = expression(Mean~giant~kelp~biomass~"("~dry~g/m^{"2"}~")"),
       y = "Predicted time to recovery (years)") +
  
  theme_bw() + 
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.position = "none")

rec_time_plot


##########################################################################-
# 5. manuscript tables ----------------------------------------------------
##########################################################################-

lm_kelp_tables <- tbl_merge(tbls = list(lm_kelp_during_summary, lm_kelp_recovery_summary), 
                            tab_spanner = c("**Removal**", "**Recovery**")) 

# lm_kelp_tables %>%
#   as_gt() %>%
#   gtsave(here::here("tables", "ms-tables", paste("lm_kelp_tables_", today(), ".png", sep = "")),
#        vwidth = 1500, vheight = 1000)

##########################################################################-
# 6. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. delta kelp through time -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_", today(), ".jpg", sep = "")),
#        plot = overall_ms,
#        height = 8, width = 14, dpi = 150)

# ⟞ b. raw kelp biomass through time --------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S2_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_raw,
#        height = 8, width = 16, dpi = 150)

# ⟞ c. recovery time vs biomass -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_", today(), ".jpg", sep = "")),
#        plot = rec_time_plot,
#        height = 8, width = 10, dpi = 150)

# ⟞ c. recovery time vs biomass -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S1_", today(), ".jpg", sep = "")),
#        plot = plots_together_sites,
#        height = 8, width = 16, dpi = 150)







