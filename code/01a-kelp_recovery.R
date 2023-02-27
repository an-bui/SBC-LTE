
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
    # after removal ended:
    site == "aque" & date >= aque_after_date_annual ~ "after",
    site == "napl" & date >= napl_after_date_annual ~ "after",
    site == "ivee" & date >= ivee_after_date_annual ~ "after", 
    site == "mohk" & date >= mohk_after_date_annual ~ "after",
    site == "carp" & date >= carp_after_date_annual ~ "after",
    # everything else is "during" removal
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
    # after removal ended
    site == "aque" & date >= aque_after_date_continual ~ "after",
    site == "napl" & date >= napl_after_date_continual ~ "after",
    site == "mohk" & date >= mohk_after_date_continual ~ "after",
    site == "carp" & date >= carp_after_date_continual ~ "after",
    # everything else is "during" removal
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
    site == "aque" ~ paste("(a) ", site_full, sep = ""),
    site == "napl" ~ paste("(b) ", site_full),
    site == "mohk" ~ paste("(c) ", site_full),
    site == "carp" ~ paste("(d) ", site_full)
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(x = time_since_end, y = control, col = site), alpha = 0.5, size = 2) +
  # control
  geom_point(aes(x = time_since_end, y = control, shape = site), size = 1, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual, col = site), linewidth = 2) +
  geom_point(aes(x = time_since_end, y = continual, shape = site, col = site), size = 1, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7), 
        legend.position = "none", 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 8, hjust = 0),
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
# normal model
lm_kelp_during_lmer <- lmerTest::lmer(
  delta_continual ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass)
# normal model with season
lm_kelp_during_season <- lmerTest::lmer(
  delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass)
# continuous AR1
lm_kelp_during_lme_car1 <- nlme::lme(
  delta_continual ~ time_since_end, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass,
  correlation = corCAR1())
# continuous AR1 with season
lm_kelp_during_lme_car1_season <- nlme::lme(
  delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass,
  correlation = corCAR1())
# ARMA 4
lm_kelp_during_lme_ar4 <- nlme::lme(
  delta_continual ~ time_since_end, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass,
  correlation = corARMA(p = 4, q = 0)) 
# ARMA 4 with season
lm_kelp_during_lme_ar4_season <- nlme::lme(
  delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass,
  correlation = corARMA(p = 4, q = 0)) 
# GLS CAR1
lm_kelp_during_gls_car1 <- nlme::gls(
  delta_continual ~ time_since_end,
  data = delta_continual %>% filter(exp_dates == "during"), 
  correlation = corCAR1(form = ~ 1|site)
)

# diagnostics
# normal model
plot(DHARMa::simulateResiduals(lm_kelp_during_lmer))
performance::check_model(lm_kelp_during_lmer)
performance::check_autocorrelation(lm_kelp_during_lmer) # Durbin-Watson-Test

# normal model with season
plot(DHARMa::simulateResiduals(lm_kelp_during_season))
performance::check_model(lm_kelp_during_season)

# continuous AR1
plot(fitted(lm_kelp_during_lme_car1), resid(lm_kelp_during_lme_car1))
plot(density(resid(lm_kelp_during_lme_car1)))
performance::check_model(lm_kelp_during_lme_car1)

# continuous AR1 with season
resid_plot_fxn(lm_kelp_during_lme_car1_season)
plot(density(resid(lm_kelp_during_lme_car1_season)))
performance::check_model(lm_kelp_during_lme_car1_season)

# ARMA 4
resid_plot_fxn(lm_kelp_during_lme_ar4)

# ARMA 4 with season
resid_plot_fxn(lm_kelp_during_lme_ar4_season)
check_model(lm_kelp_during_lme_ar4_season)

# GLS CAR1
performance::check_model(lm_kelp_during_gls_car1)

# model checks
# normal model
performance::check_convergence(lm_kelp_during_lmer)
performance::check_normality(lm_kelp_during_lmer) # Shapiro test
performance::check_homogeneity(lm_kelp_during_lmer) # Bartlett test
performance::check_heteroskedasticity(lm_kelp_during_lmer) # Breusch-Pagan test

# plot ACF/PACF
# normal model
acf(resid(lm_kelp_during_lmer))
pacf(resid(lm_kelp_during_lmer))

# continuous AR1
acf(resid(lm_kelp_during_lme_car1))
pacf(resid(lm_kelp_during_lme_car1))

# ARMA 4
acf(resid(lm_kelp_during_lme_ar4))
pacf(resid(lm_kelp_during_lme_ar4))

# GLS CAR1
acf(resid(lm_kelp_during_gls_car1))
pacf(resid(lm_kelp_during_gls_car1))

# Rsquared
MuMIn::r.squaredGLMM(lm_kelp_during_lmer)
MuMIn::r.squaredGLMM(lm_kelp_during_lme_car1)
MuMIn::r.squaredGLMM(lm_kelp_during_lme_ar4)

# summaries
summary(lm_kelp_during_lmer) # significant slope
summary(lm_kelp_during_lme_car1) # significant slope
summary(lm_kelp_during_gls_car1) # non significant slope
summary(lm_kelp_during_lme_ar4)
summary(lm_kelp_during_lme_ar4_season)
lm_kelp_during_summary <- lm_kelp_during_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_kelp_during_summary

# AIC comparison
MuMIn::AICc(lm_kelp_during_lmer, lm_kelp_during_season,
            lm_kelp_during_lme_car1, lm_kelp_during_lme_car1_season, 
            lm_kelp_during_lme_ar4, lm_kelp_during_lme_ar4_season,
            lm_kelp_during_gls_car1)
# ARMA 4 best model?

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_kelp_during_overall <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "fixed")
predicted_kelp_during_aque <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
predicted_kelp_during_napl <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
predicted_kelp_during_mohk <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
predicted_kelp_during_carp <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
# normal model
lm_kelp_recovery_lmer <- lmerTest::lmer(
  delta_continual ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass)
# continuous AR1
lm_kelp_recovery_lme_ar1 <- nlme::lme(
  delta_continual ~ time_since_end, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass,
  correlation = corAR1())
# ARMA 2
lm_kelp_recovery_lme_ar2 <- nlme::lme(
  delta_continual ~ time_since_end, random = ~1|site,
  data = delta_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass,
  correlation = corARMA(p = 2, q = 0))
# GLS AR1
lm_kelp_recovery_gls_ar1 <- nlme::gls(
  delta_continual ~ time_since_end, 
  data = delta_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass,
  correlation = corAR1(form = ~1|site))

# check for autocorrelation
performance::check_autocorrelation(lm_kelp_recovery_lmer)

# diagnostics
# normal model
plot(DHARMa::simulateResiduals(lm_kelp_recovery_lmer))
performance::check_model(lm_kelp_recovery_lmer)

# continuous AR1
performance::check_model(lm_kelp_recovery_lme_ar1)

# GLS AR1
qqnorm(lm_kelp_recovery_gls_ar1)
plot(fitted(lm_kelp_recovery_gls_ar1), residuals(lm_kelp_recovery_gls_ar1))

# model checks
# normal model
check_convergence(lm_kelp_recovery_lmer)
check_normality(lm_kelp_recovery_lmer)
check_heteroscedasticity(lm_kelp_recovery_lmer)

# plot ACF
# normal model
acf(residuals(lm_kelp_recovery_lmer))
pacf(residuals(lm_kelp_recovery_lmer))

# ARMA 2
acf(residuals(lm_kelp_recovery_lme_ar2))
pacf(residuals(lm_kelp_recovery_lme_ar2))

# Rsquared
MuMIn::r.squaredGLMM(lm_kelp_recovery_lmer)
MuMIn::r.squaredGLMM(lm_kelp_recovery_lme_ar1)

# summary
summary(lm_kelp_recovery_lmer)
summary(lm_kelp_recovery_lme_ar1)
summary(lm_kelp_recovery_gls_ar1)
lm_kelp_recovery_summary <- lm_kelp_recovery_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_kelp_recovery_summary

# AIC comparisons
MuMIn::AICc(lm_kelp_recovery_lmer, lm_kelp_recovery_lme_ar1, lm_kelp_recovery_lme_ar2, lm_kelp_recovery_gls_ar1)
# GLS AR1 best model?

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
# predicted time to recovery: 3.8 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.7:3.8 by = 0.001]", type = "random", condition = c(site = "aque"))
# PI: -3.5 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-3.5:-3.44 by = 0.001]", type = "random", condition = c(site = "aque"))
# PI: 12.0 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [12:12.1 by = 0.001]", type = "random", condition = c(site = "aque"))

# napl: 3.4 years
predicted_kelp_after_napl <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
# predicted time to recovery: 3.4 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.3:3.5 by = 0.001]", type = "random", condition = c(site = "napl"))
# PI: -4.0 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-4:-3.9 by = 0.001]", type = "random", condition = c(site = "napl"))
# PI: 11.5 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [11.5:11.6 by = 0.001]", type = "random", condition = c(site = "napl"))

# mohk: 5.4 years
predicted_kelp_after_mohk <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
# predicted time to recovery: 5.4 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [5.37:5.41 by = 0.01]", type = "random", condition = c(site = "mohk"))
# PI: -1.6 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-1.63:-1.5 by = 0.01]", type = "random", condition = c(site = "mohk"))
# PI: 14.4 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [14.37:14.5 by = 0.01]", type = "random", condition = c(site = "mohk"))

# carp: 3.3 years
predicted_kelp_after_carp <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))
# predicted time to recovery: 3.3 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.2:3.4 by = 0.01]", type = "random", condition = c(site = "carp"))
# PI: -4.1 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-4.25:-4 by = 0.01]", type = "random", condition = c(site = "carp"))
# PI: 11.4 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [11.3:11.4 by = 0.01]", type = "random", condition = c(site = "carp"))

# ⟞ c. figures -------------------------------------------------------------

# ⟞ ⟞ i. overall model predictions -----------------------------------------

overall_ms <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual, 
             aes(x = time_since_end, y = delta_continual, fill = site, shape = site), size = 2, alpha = 0.9) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  # new_scale("color") + 
  # overall
  geom_line(data = predicted_kelp_after_overall, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  geom_line(data = predicted_kelp_during_overall, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(-1500, 1000, by = 1000), limits = c(-1800, 900)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        # plot.margin = margin(0, 0, 0, 0),
        legend.position = c(0.1, 0.858),
        legend.key.size = unit(0.3, units = "cm")) +
  labs(x = "Time since end of removal (years)", 
       y = expression("\U0394"~giant~kelp~biomass~"(treatment - control, "~dry~g/m^{"2"}~")"), 
       fill = "Site",
       shape = "Site")

overall_ms

# ⟞ ⟞ ii. site level predictions -------------------------------------------

aque <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "aque"), aes(x = time_since_end, y = delta_continual), shape = aque_shape, fill = aque_col, size = 2, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_aque, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_aque, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1000, by = 1000), limits = c(-1800, 1000)) +
  # geom_text(aes(x = -6.6, y = 600), label = "aque", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = paste("(a) ", aque_full, sep = ""))

napl <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "napl"), aes(x = time_since_end, y = delta_continual), shape = napl_shape, fill = napl_col, size = 2, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_napl, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_napl, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1500, by = 1000), limits = c(-1800, 1500)) +
  # geom_text(aes(x = -6.8, y = 700), label = "napl", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = paste("(b) ", napl_full, sep = ""))

mohk <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "mohk"), aes(x = time_since_end, y = delta_continual), shape = mohk_shape, fill = mohk_col, size = 2, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_mohk, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_mohk, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1500, by = 1000), limits = c(-1800, 1500)) +
  # geom_text(aes(x = -6.6, y = 900), label = "mohk", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = paste("(c) ", mohk_full, sep = ""))

carp <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_continual %>% filter(site == "carp"), aes(x = time_since_end, y = delta_continual), shape = carp_shape, fill = carp_col, size = 2, alpha = 0.9) +
  # during
  geom_line(data = predicted_kelp_during_carp, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_during_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # after
  geom_line(data = predicted_kelp_after_carp, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_kelp_after_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-1500, 1500, by = 1000), limits = c(-1800, 1500)) +
  # geom_text(aes(x = -6.8, y = 1000), label = "carp", size = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  labs(x = "Time since end of removal", 
       y = "\U0394 giant kelp biomass",
       title = paste("(d) ", carp_full, sep = ""))

plots_together_sites <- (aque + napl) / (mohk + carp)
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
  "aque",       3.80,          -3.5,         12,
  "napl",       3.42,          -4,           11.5,
  "mohk",       5.4,           -1.6,         14.4,
  "carp",       3.33,          -4.1,         11.4
) %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  left_join(., mean_kelp_all_sites, by = "site") 

# ⟞ b. figure -------------------------------------------------------------

rec_time_plot <- ggplot(rec_time, aes(x = mean_control, y = time_to_recovery, shape = site, fill = site)) +
  geom_errorbar(aes(xmin = mean_control - se_control, xmax = mean_control + se_control), width = 1) +
  geom_errorbar(aes(ymin = pred_int_low, ymax = pred_int_high)) +
  geom_point(size = 4) +
  # geom_text_repel(aes(label = site_full), seed = 666,
  #                 size = 7, point.padding = 35, nudge_y = 0.1) +
  scale_shape_manual(values = shape_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  theme_bw() +
  # scale_y_continuous(breaks = c(seq(3, 6, by = 1)), limits = c(3, 6)) +
  labs(x = expression(Mean~giant~kelp~biomass~"("~dry~g/m^{"2"}~")"),
       y = "Predicted time to recovery (years)") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7),
        legend.position = "none")

rec_time_plot

##########################################################################-
# 5. control vs removal kelp biomass plot ---------------------------------
##########################################################################-

# arrows
delta_continual %>% 
  ggplot(aes(x = control, y = continual)) +
  scale_color_manual(values = color_palette_site) +
  geom_abline(slope = 1) +
  geom_segment(
    data = delta_continual %>% filter(site == "aque", exp_dates == "after"),
    aes(x = control, y = continual,
        xend = c(tail(control, n = -1), NA),
        yend = c(tail(continual, n = -1), NA),
        color = site),
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  geom_segment(
    data = delta_continual %>% filter(site == "carp", exp_dates == "after"),
    aes(x = control, y = continual,
        xend = c(tail(control, n = -1), NA),
        yend = c(tail(continual, n = -1), NA),
        color = site),
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  geom_segment(
    data = delta_continual %>% filter(site == "mohk", exp_dates == "after"),
    aes(x = control, y = continual,
        xend = c(tail(control, n = -1), NA),
        yend = c(tail(continual, n = -1), NA),
        color = site),
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  geom_segment(
    data = delta_continual %>% filter(site == "napl", exp_dates == "after"),
    aes(x = control, y = continual,
        xend = c(tail(control, n = -1), NA),
        yend = c(tail(continual, n = -1), NA),
        color = site),
    arrow = arrow(length = unit(0.3, "cm"))
  )

delta_vs_biomass <- delta_continual %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = control, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(fill = site, shape = site), size = 2) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  geom_smooth(method = "lm", color = "#000000", se = FALSE, linewidth = 1) +
  labs(x = expression(Biomass~"("~dry~g/m^{"2"}~")"), 
       y = "\U0394 biomass (treatment - control)",
       fill = "Site", shape = "Site") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))
  
delta_vs_biomass

##########################################################################-
# 6. manuscript tables ----------------------------------------------------
##########################################################################-

lm_kelp_tables <- tbl_merge(tbls = list(lm_kelp_during_summary, lm_kelp_recovery_summary), 
                            tab_spanner = c("**Removal**", "**Recovery**")) 

# lm_kelp_tables %>%
#   as_gt() %>%
#   gtsave(here::here("tables", "ms-tables", paste("lm_kelp_tables_", today(), ".png", sep = "")),
#        vwidth = 1500, vheight = 1000)

##########################################################################-
# 7. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. delta kelp through time -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_", today(), ".jpg", sep = "")),
#        plot = overall_ms,
#        height = 8, width = 14, units = "cm",
#        dpi = 400)

# ⟞ b. raw kelp biomass through time --------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S2_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_raw,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)

# ⟞ c. recovery time vs biomass -------------------------------------------

s4_panels <- rec_time_plot + delta_vs_biomass +
  plot_layout(widths = c(2, 2.1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12))

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_", today(), ".jpg", sep = "")),
#        plot = rec_time_plot,
#        height = 8, width = 12, units = "cm",
#        dpi = 300)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_panels_", today(), ".jpg", sep = "")),
#        plot = s4_panels,
#        height = 8, width = 16, units = "cm",
#        dpi = 400)

# ⟞ d. raw biomass by site ------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S1_", today(), ".jpg", sep = "")),
#        plot = plots_together_sites,
#        height = 10, width = 16, units = "cm",
#        dpi = 300)







