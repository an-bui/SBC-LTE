##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "00-set_up.R"))

##########################################################################-
# 1. data frames ----------------------------------------------------------
##########################################################################-

# ⟞ a. annual removal deltas ----------------------------------------------

# delta_annual <- biomass %>% 
#   filter(sp_code == "MAPY" & treatment %in% c("control", "annual")) %>% 
#   dplyr::select(-sp_code) %>% 
#   dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
#   pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
#   # fill in values for sampling dates where control and annual were surveyed on different days
#   mutate(control = case_when(
#     site == "aque" & date == "2008-03-07" ~ 106.82000,
#     site == "napl" & date == "2012-09-26" ~ 143.22088,
#     TRUE ~ control
#   )) %>% 
#   mutate(delta_annual = annual - control) %>%  
#   # missing dates are from AQUE 2008-03-05, NAPL 2012-09-25, NAPL 2008-10-10
#   drop_na(delta_annual) %>% 
#   mutate(exp_dates = case_when(
#     # after removal ended:
#     site == "aque" & date >= aque_after_date_annual ~ "after",
#     site == "napl" & date >= napl_after_date_annual ~ "after",
#     site == "ivee" & date >= ivee_after_date_annual ~ "after", 
#     site == "mohk" & date >= mohk_after_date_annual ~ "after",
#     site == "carp" & date >= carp_after_date_annual ~ "after",
#     # everything else is "during" removal
#     TRUE ~ "during"
#   ),
#   exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
#   time_since_columns_annual() %>% 
#   kelp_year_column() %>% 
#   comparison_column_annual() 

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
  comparison_column_continual_new() %>% 
  left_join(., site_quality, by = "site") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))

# kelp biomass in long format
continual_long <- delta_continual %>% 
  select(!delta_continual) %>% 
  pivot_longer(cols = c(control, continual)) %>% 
  rename(kelp_biomass = value, treatment = name) %>% 
  mutate(treatment = case_match(treatment, "control" ~ "reference", "continual" ~ "removal")) %>% 
  unite("sample_ID", site, date, quarter, treatment, remove = FALSE) %>% 
  # calculate variation by observation
  group_by(exp_dates, treatment) %>% 
  mutate(mean = mean(kelp_biomass),
         variation = ((kelp_biomass - mean(kelp_biomass))/mean(kelp_biomass))^2) %>% 
  ungroup()

# ⟞ c. density ------------------------------------------------------------

density <- biomass %>% 
  filter(sp_code == "MAPY" & treatment %in% c("control", "continual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, density) %>% 
  pivot_wider(names_from = treatment, values_from = density) %>% 
  mutate(delta_continual_density = continual - control) %>%  
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual_density) %>% 
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
  left_join(., site_quality, by = "site") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))

# ⟞ d. fronds -------------------------------------------------------------

fronds_clean <- fronds %>% 
  filter(treatment %in% c("control", "continual")) %>% 
  dplyr::select(-sp_code) %>% 
  dplyr::select(site, year, month, treatment, date, fronds) %>% 
  pivot_wider(names_from = treatment, values_from = fronds) %>% 
  mutate(delta_continual_fronds = continual - control) %>%  
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual_fronds) %>% 
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
  left_join(., site_quality, by = "site") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))

##########################################################################-
# 2. timeseries plots -----------------------------------------------------
##########################################################################-

# ⟞ a. delta annual plot --------------------------------------------------

# delta_annual_plot <- ggplot(delta_annual, aes(x = time_since_end, y = delta_annual)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(aes(col = exp_dates)) +
#   geom_smooth(method = "lm", aes(col = exp_dates)) +
#   theme_bw() + 
#   theme(legend.position = "none",
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 18),
#         axis.text = element_text(size = 10),
#         strip.text = element_text(size = 10)) +
#   labs(x = "Time since end of experiment", y = "\U0394 biomass (removal - reference)",
#        title = "\U0394 annual") +
#   facet_wrap(~site_full, ncol = 1, scales = "free_y")

# delta_annual_plot

# ⟞ b. delta continual plot -----------------------------------------------

# delta_continual_plot <- ggplot(delta_continual, aes(x = time_since_end, y = delta_continual)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(aes(col = exp_dates)) +
#   geom_smooth(method = "lm", aes(col = exp_dates)) +
#   theme_bw() + 
#   theme(legend.position = "none",
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 18),
#         axis.text = element_text(size = 10),
#         strip.text = element_text(size = 10)) +
#   labs(x = "Time since end of experiment", y = "\U0394 biomass (removal - reference)",
#        title = "\U0394 continual") +
#   facet_wrap(~site_full, ncol = 1, scales = "free_y")
# 
# delta_continual_plot

# ⟞ c. raw biomass through time plot --------------------------------------

continual_sites_raw <- delta_continual %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("(a) ", site_full, sep = ""),
    site == "napl" ~ paste("(b) ", site_full, sep = ""),
    site == "mohk" ~ paste("(c) ", site_full, sep = ""),
    site == "carp" ~ paste("(d) ", site_full, sep = "")
  )) %>% 
  mutate(removal_annotation = case_when(
    sample_ID == "aque_2010-06-15_Q2" ~ "Removal"
    # sample_ID == "napl_2010-06-11_Q2" ~ "Removal",
    # sample_ID == "mohk_2010-07-23_Q3" ~ "Removal",
    # sample_ID == "carp_2010-06-09_Q2" ~ "Removal"
  ),
  recovery_annotation = case_when(
    sample_ID == "aque_2010-06-15_Q2" ~ "Recovery"
    # sample_ID == "napl_2010-06-11_Q2" ~ "Recovery",
    # sample_ID == "mohk_2010-07-23_Q3" ~ "Recovery",
    # sample_ID == "carp_2010-06-09_Q2" ~ "Recovery"
  ),
  annotation_y = case_when(
    sample_ID == "aque_2010-06-15_Q2" ~ 840
    # sample_ID == "napl_2010-06-11_Q2" ~ 1010,
    # sample_ID == "mohk_2010-07-23_Q3" ~ 1700,
    # sample_ID == "carp_2010-06-09_Q2" ~ 650
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey") +
  annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
            fill = "grey", alpha = 0.3) +
  geom_text(aes(x = -6.75, y = annotation_y, label = removal_annotation), size = 2) +
  geom_text(aes(x = 5.5, y = annotation_y, label = recovery_annotation), size = 2) +
  # control
  geom_line(aes(x = time_since_end, y = control), 
            alpha = 0.9, 
            linewidth = 0.5,
            linetype = "B3",
            color = reference_col) +
  # geom_point(aes(x = time_since_end, y = control), shape = 21, color = "#FFFFFF", stroke = 0.5, fill = reference_col, size = 0.75) +
  # continual
  geom_line(aes(x = time_since_end, y = continual), 
            alpha = 0.9, 
            linewidth = 0.5,
            linetype = 1,
            color = removal_col) +
  # geom_point(aes(x = time_since_end, y = continual), shape = 21, color = "#FFFFFF", stroke = 0.5, fill = removal_col, size = 0.75) +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
  raw_biomass_plot_theme() +
  labs(x = "Time since end of experiment (years)", 
       y = "Giant kelp biomass (dry g/m\U00B2)") +
  facet_wrap(~strip, scales = "free_y", nrow = 4)

continual_sites_raw


# ⟞ d. summer biomass -----------------------------------------------------

summer_biomass <- continual_long %>% 
  filter(quarter == "Q3") %>% 
  ggplot(aes(x = time_since_end, y = kelp_biomass)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(color = site, alpha = treatment), linewidth = 2) +
  scale_color_manual(values = color_palette_site) +
  geom_point(color = "#FFFFFF", size = 1) +
  scale_alpha_manual(values = c(reference = 0.4, removal = 1)) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  facet_wrap(~site_full) +
  labs(x = "Time since end (years)",
       y = "Giant kelp biomass (dry g/m\U00B2)",
       title = "Summer kelp biomass")
 
summer_biomass


# ⟞ e. density through time -----------------------------------------------
  
density_timeseries <- density %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("(a) ", site_full, sep = ""),
    site == "napl" ~ paste("(b) ", site_full, sep = ""),
    site == "mohk" ~ paste("(c) ", site_full, sep = ""),
    site == "carp" ~ paste("(d) ", site_full, sep = "")
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  # control
  geom_line(aes(x = time_since_end, y = control, col = site), alpha = 0.3, linewidth = 2) +
  geom_point(aes(x = time_since_end, y = control, shape = site), size = 1, alpha = 0.3, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual, col = site), linewidth = 2) +
  geom_point(aes(x = time_since_end, y = continual, shape = site, col = site), size = 1, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(breaks = seq(-8, 8, by = 1), minor_breaks = NULL) +
  raw_biomass_plot_theme() +
  theme(plot.title = element_text(size = 10)) +
  labs(x = "Time since end of experiment (years)", 
       y = "Giant kelp density",
       title = "Density") +
  facet_wrap(~strip, scales = "free_y")

density_timeseries


# ⟞ f. fronds through time ------------------------------------------------

fronds_timeseries <- fronds_clean %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("(a) ", site_full, sep = ""),
    site == "napl" ~ paste("(b) ", site_full, sep = ""),
    site == "mohk" ~ paste("(c) ", site_full, sep = ""),
    site == "carp" ~ paste("(d) ", site_full, sep = "")
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  # control
  geom_line(aes(x = time_since_end, y = control, col = site), alpha = 0.3, linewidth = 2) +
  geom_point(aes(x = time_since_end, y = control, shape = site), size = 1, alpha = 0.3, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual, col = site), linewidth = 2) +
  geom_point(aes(x = time_since_end, y = continual, shape = site, col = site), size = 1, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(breaks = seq(-8, 8, by = 1), minor_breaks = NULL) +
  raw_biomass_plot_theme() +
  theme(plot.title = element_text(size = 10)) +
  labs(x = "Time since end of experiment (years)", 
       y = "Giant kelp fronds (number/m\U00B2)",
       title = "Fronds") +
  facet_wrap(~strip, scales = "free_y")

fronds_timeseries

##########################################################################-
# 3. linear models --------------------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# delta_hist <- delta_continual %>% 
#   filter(exp_dates == "during") %>% 
#   ggplot(aes(x = delta_continual)) +
#   geom_histogram(bins = 10)
# delta_hist

# model
# normal model
# lm_kelp_during_lmer <- lmerTest::lmer(
#   delta_continual ~ time_since_end + (1|site),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass)
# normal model with zero inflated negative binomial structure
# lm_kelp_during_zinb <- glmmTMB(
#   delta_continual ~ time_since_end + (1|site),
#   data = delta_continual %>% filter(exp_dates == "during"),
#   ziformula = ~time_since_end, 
#   family = nbinom2,
#   na.action = na.pass
# )
# normal model using glmmTMB with tweedie structure
# lm_kelp_during_tweedie <- glmmTMB(
#   delta_continual ~ time_since_end + (1|site),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   family = tweedie(),
#   na.action = na.pass, 
#   REML = TRUE)
# normal model with season
# lm_kelp_during_season <- lmerTest::lmer(
#   delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass)
# normal model with season using glmmTM with tweedie structure
# lm_kelp_during_season_tweedie <- glmmTMB(
#   delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   family = tweedie(),
#   na.action = na.pass)
# normal model with season as random effect (failed to converge)
# lm_kelp_during_season_re <- lmerTest::lmer(
#   delta_continual ~ time_since_end + (1|site) + (1|quarter),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass)
# continuous AR1
# lm_kelp_during_lme_car1 <- nlme::lme(
#   delta_continual ~ time_since_end, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corCAR1())
# continuous AR1 with season
# lm_kelp_during_lme_car1_season <- nlme::lme(
#   delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corCAR1())
# continuous AR1 with season as random effect
# lm_kelp_during_lme_car1_season_re <- nlme::lme(
#   delta_continual ~ time_since_end, random = list(site = ~1, quarter = ~1),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corCAR1())
# AR1
# lm_kelp_during_lme_ar1 <- nlme::lme(
#   delta_continual ~ time_since_end, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corAR1())
# AR1 with season
# lm_kelp_during_lme_ar1_season <- nlme::lme(
#   delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corAR1())
# AR1 with season as random effect
# lm_kelp_during_lme_ar1_season_re <- nlme::lme(
#   delta_continual ~ time_since_end, random = list(site = ~1, quarter = ~1),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corAR1())
# ARMA 4
# lm_kelp_during_lme_ar4 <- nlme::lme(
#   delta_continual ~ time_since_end, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corARMA(p = 4, q = 0)) 
# with season no autoregressive
# lm_kelp_during_glmmtmb_season <- glmmTMB(
#   delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   family = tweedie(),
#   na.action = na.pass) 

# ARMA 4 with season
# lm_kelp_during_lme_ar4_season <- nlme::lme(
#   delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corARMA(p = 4, q = 0)) 
# with season
# lm_kelp_during_glmmTMB_season <- glmmTMB(
#   delta_continual ~ quarter + ar1(time_since_end + 0 |site),
#   data = delta_continual %>% filter(exp_dates == "during") %>% mutate(time_since_end = as_factor(time_since_end)), 
#   family = tweedie(),
#   na.action = na.pass) 
# ARMA 4 with season as random effect
# lm_kelp_during_lme_ar4_season_re <- nlme::lme(
#   delta_continual ~ time_since_end, random = list(site = ~1, quarter = ~1),
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   correlation = corARMA(p = 4, q = 0)) 
# GLS CAR1
# lm_kelp_during_gls_car1 <- nlme::gls(
#   delta_continual ~ time_since_end,
#   data = delta_continual %>% filter(exp_dates == "during"), 
#   correlation = corCAR1(form = ~ 1|site)
# )
# linear transformation: multiply by -1, add 8
# transform <- delta_continual %>% 
#   mutate(delta_continual_tf = delta_continual*-1 + 8) %>% 
#   mutate(logratio = log(continual)/log(control)) %>% 
#   filter(exp_dates == "during")
# ggplot(transform, aes(x = logratio)) +
#   geom_histogram()
# 
# ggplot(transform, aes(x = control)) +
#   geom_histogram(bins = 10)
# transform_fxn <- function(x) x*-1 + 8
# backtransform <- function(x) (x-8)*-1
# taking out aque, mohk, or carp makes residuals look normal
# lm_kelp_during_tf_gamma <- glmmTMB::glmmTMB(
#   transform_fxn(delta_continual) ~ time_since_end + (1|site),
#   data = transform,
#   na.action = na.pass,
#   family = Gamma(link = "log"))
# simulateResiduals(lm_kelp_during_tf_gamma, plot = TRUE)
# summary(lm_kelp_during_tf_gamma)
# tf_pred <- ggpredict(lm_kelp_during_tf_gamma, terms = ~ time_since_end) %>% 
#   mutate(transform = backtransform(predicted),
#          tf_conf.low = backtransform(conf.low),
#          tf_conf.high = backtransform(conf.high))
# ggplot(delta_continual, aes(x = time_since_end, y = delta_continual)) +
#   geom_point() +
#   geom_line(data = tf_pred, aes(x = x, y = transform), color = "blue") +
#   geom_ribbon(data = tf_pred, aes(x = x, y = transform, ymin = tf_conf.low, ymax = tf_conf.high), alpha = 0.2)
# lm_kelp_during_tf_tweedie <- glmmTMB::glmmTMB(
#   delta_continual_tf ~ time_since_end + (1|site),
#   data = transform, 
#   na.action = na.pass,
#   family = tweedie(link = "log")) 
# lm_kelp_during_tf_season <- glmmTMB::glmmTMB(
#   delta_continual_tf ~ time_since_end + quarter + time_since_end*quarter + (1|site),
#   data = transform, 
#   na.action = na.pass,
#   family = Gamma(link = "inverse"),
#   start = list(psi = c(-1, 1))) 

# raw kelp biomass model
## site random effect
lm_kelp_during_zigamma_01 <- glmmTMB(
  kelp_biomass ~ time_since_end*treatment + (1|site),
  data = continual_long %>% filter(exp_dates == "during"), 
  family = ziGamma(link = "log"),
  ziformula = ~1)

## site and year random effect
lm_kelp_during_zigamma_02 <- glmmTMB(
  kelp_biomass ~ time_since_end*treatment + (1|site) + (1|year),
  data = continual_long %>% filter(exp_dates == "during"), 
  family = ziGamma(link = "log"),
  ziformula = ~1)

# df <- delta_continual %>% 
#   filter(exp_dates == "during") %>% 
#   cbind(., residuals(lm_kelp_during_lmer)) %>% 
#   dplyr::rename(resid_m1 = 'residuals(lm_kelp_during_lmer)') %>% 
#   cbind(., residuals(lm_kelp_during_season)) %>% 
#   dplyr::rename(resid_m2 = 'residuals(lm_kelp_during_season)') %>% 
#   cbind(., residuals(lm_kelp_during_lme_ar4_season)) %>% 
#   dplyr::rename(resid_m3 = 'residuals(lm_kelp_during_lme_ar4_season)')
# 
# ggplot(df, aes(x = time_since_end, y = resid_m3)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE)

# diagnostics
# normal model
# DHARMa::simulateResiduals(lm_kelp_during_lmer, plot = T)
# performance::check_model(lm_kelp_during_lmer)
# performance::check_autocorrelation(lm_kelp_during_lmer) # Durbin-Watson-Test

# normal model with season
# simulateResiduals(lm_kelp_during_season, plot = T)
# check_model(lm_kelp_during_season)

# normal model with season with random effect
# simulateResiduals(lm_kelp_during_season_re, plot = T)
# check_model(lm_kelp_during_season_re)

# continuous AR1
# resid_plot_fxn(lm_kelp_during_lme_car1)
# plot(density(resid(lm_kelp_during_lme_car1)))
# check_model(lm_kelp_during_lme_car1)

# continuous AR1 with season
# resid_plot_fxn(lm_kelp_during_lme_car1_season)
# plot(density(resid(lm_kelp_during_lme_car1_season)))
# check_model(lm_kelp_during_lme_car1_season)

# continuous AR1 with season as random effect
# resid_plot_fxn(lm_kelp_during_lme_car1_season_re)
# plot(density(resid(lm_kelp_during_lme_car1_season_re)))
# check_model(lm_kelp_during_lme_car1_season_re)

# AR1
# resid_plot_fxn(lm_kelp_during_lme_ar1)
# check_model(lm_kelp_during_lme_ar1)

# AR1 with season
# resid_plot_fxn(lm_kelp_during_lme_ar1_season)

# AR1 with season as random effect
# resid_plot_fxn(lm_kelp_during_lme_ar1_season_re)

# ARMA 4
# resid_plot_fxn(lm_kelp_during_lme_ar4)
# qqnorm(residuals(lm_kelp_during_lme_ar4))
# qqline(residuals(lm_kelp_during_lme_ar4))
# ks.test(residuals(lm_kelp_during_lme_ar4), "pnorm")

# ARMA 4 with season
# resid_plot_fxn(lm_kelp_during_lme_ar4_season)
# check_model(lm_kelp_during_lme_ar4_season)

# ARMA 4 with season as random effect
# resid_plot_fxn(lm_kelp_during_lme_ar4_season_re)
# check_model(lm_kelp_during_lme_ar4_season_re)

# GLS CAR1
# check_model(lm_kelp_during_gls_car1)

# model checks
# normal model
# testUniformity(lm_kelp_during_lmer)
# performance::check_convergence(lm_kelp_during_lmer)
# performance::check_normality(lm_kelp_during_lmer) # Shapiro test
# performance::check_homogeneity(lm_kelp_during_lmer) # Bartlett test
# performance::check_heteroskedasticity(lm_kelp_during_lmer) # Breusch-Pagan test

# plot ACF/PACF
# normal model
# acf(resid(lm_kelp_during_lmer))
# pacf(resid(lm_kelp_during_lmer))

# continuous AR1
# acf(resid(lm_kelp_during_lme_car1))
# pacf(resid(lm_kelp_during_lme_car1))

# ARMA 4
# acf(resid(lm_kelp_during_lme_ar4))
# pacf(resid(lm_kelp_during_lme_ar4))

# GLS CAR1
# acf(resid(lm_kelp_during_gls_car1))
# pacf(resid(lm_kelp_during_gls_car1))

# raw kelp biomass model
# DHARMa::simulateResiduals(lm_kelp_during_zigamma_01, plot = T)
# performance::check_model(lm_kelp_during_zigamma_01)

DHARMa::simulateResiduals(lm_kelp_during_zigamma_02, plot = T)
performance::check_model(lm_kelp_during_zigamma_02)

# Rsquared
# MuMIn::r.squaredGLMM(lm_kelp_during_lmer)
# r.squaredGLMM(lm_kelp_during_lme_car1)
# r.squaredGLMM(lm_kelp_during_lme_ar4)
# r.squaredGLMM(lm_kelp_during_zigamma_01)
r.squaredGLMM(lm_kelp_during_zigamma_02)

# summaries
# summary(lm_kelp_during_lmer)
# summary(lm_kelp_during_zigamma_01)
summary(lm_kelp_during_zigamma_02)

lm_kelp_during_zigamma_summary <- lm_kelp_during_zigamma_02 %>% 
  tbl_regression(intercept = TRUE) %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Estimate**"
  ) %>% 
  modify_column_indent(
    columns = label, 
    rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))

# filter out zero-inflated component
lm_kelp_during_zigamma_summary$table_body <- lm_kelp_during_zigamma_summary$table_body %>% 
  filter(component != "zi")
# change labels
lm_kelp_during_zigamma_summary$table_body$label <- c(
  `(Intercept)` = "(Intercept)",
  time_since_end = "Time since end",
  treatmentremoval = "Treatment (removal)",
  `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
)

# final table 
lm_kelp_during_zigamma_summary


# lm_kelp_during_summary <- lm_kelp_during_lmer %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**",
#     df = "**df**"
#   ) 
# lm_kelp_during_summary

# AIC comparison
# MuMIn::AICc(lm_kelp_during_lmer, 
#             lm_kelp_during_season, lm_kelp_during_season_re,
#             lm_kelp_during_lme_car1, lm_kelp_during_lme_car1_season,
#             lm_kelp_during_lme_car1_season_re, lm_kelp_during_lme_ar1,
#             lm_kelp_during_lme_ar1_season, lm_kelp_during_lme_ar1_season_re,
#             lm_kelp_during_lme_ar4, lm_kelp_during_lme_ar4_season,
#             lm_kelp_during_lme_ar4_season_re, lm_kelp_during_gls_car1, 
#             lm_kelp_during_tweedie, lm_kelp_during_season_tweedie
#             ) %>% 
#   arrange(AICc)

# MuMIn::AICc(lm_kelp_during_zigamma_01, lm_kelp_during_zigamma_02)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# predicted_kelp_during_overall <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "fixed")
# predicted_kelp_during_aque <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
# predicted_kelp_during_napl <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
# predicted_kelp_during_mohk <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
# predicted_kelp_during_carp <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))

# raw kelp biomass
predicted_kelp_during_raw <- ggpredict(lm_kelp_during_zigamma_02,
                                      terms = c("time_since_end [-7.25:0 by = 0.25]", "treatment"), type = "fixed")


# ⟞ ⟞ iii. means ----------------------------------------------------------

# estimates predicted mean values of response for each level in predictor
# using modelbased
during_response_means <- estimate_means(lm_kelp_during_zigamma_02, 
                               transform = "response") %>% 
  mutate(exp_dates = "during")
# using emmeans
during_response_emmeans <- emmeans(lm_kelp_during_zigamma_02, 
                          ~ treatment | time_since_end, type = "response") 

# calculates differences between levels of categorical predictor
contrasts <- estimate_contrasts(lm_kelp_during_zigamma_02, 
                                transform = "response") 
# when transforming to response scale, gives ratio of means (weird?)
contrasts

# estimates slopes of numeric predictors for levels of categorical predictors
# in summary object: estimate = difference in slope of time since end (removal - reference)
during_slopes <- estimate_slopes(lm_kelp_during_zigamma_02, 
                                 trend = "time_since_end",
                                 at = "treatment") %>% 
  as_tibble() %>% 
  mutate(exp_dates = "Experimental removal")
during_emtrends <- emmeans::emtrends(lm_kelp_during_zigamma_02, "treatment", var = "time_since_end")

estimate_relation(lm_kelp_during_zigamma_02) %>% plot()

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
# normal model
# lm_kelp_recovery_lmer <- lmerTest::lmer(
#   delta_continual ~ time_since_end + (1|site),
#   data = delta_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass)
# continuous AR1
# lm_kelp_recovery_lme_ar1 <- nlme::lme(
#   delta_continual ~ time_since_end, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass,
#   correlation = corAR1())
# ARMA 2
# lm_kelp_recovery_lme_ar2 <- nlme::lme(
#   delta_continual ~ time_since_end, random = ~1|site,
#   data = delta_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass,
#   correlation = corARMA(p = 2, q = 0))
# GLS AR1
# lm_kelp_recovery_gls_ar1 <- nlme::gls(
#   delta_continual ~ time_since_end, 
#   data = delta_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass,
#   correlation = corAR1(form = ~1|site))

# raw kelp biomass model
## site random effect
lm_kelp_recovery_zigamma_01 <- glmmTMB(
  kelp_biomass ~ time_since_end*treatment + (1|site),
  data = continual_long %>% filter(exp_dates == "after"), 
  family = ziGamma(link = "log"), 
  ziformula = ~1)

## site and year random effect
lm_kelp_recovery_zigamma_02 <- glmmTMB(
  kelp_biomass ~ time_since_end*treatment + (1|site) + (1|year),
  data = continual_long %>% filter(exp_dates == "after"), 
  family = ziGamma(link = "log"), 
  ziformula = ~1)


# check for autocorrelation
# performance::check_autocorrelation(lm_kelp_recovery_lmer)
performance::check_autocorrelation(lm_kelp_recovery_zigamma_01)
performance::check_autocorrelation(lm_kelp_recovery_zigamma_02)

# diagnostics
# normal model
# plot(DHARMa::simulateResiduals(lm_kelp_recovery_lmer))
# performance::check_model(lm_kelp_recovery_lmer)

# continuous AR1
# performance::check_model(lm_kelp_recovery_lme_ar1)

# GLS AR1
# qqnorm(lm_kelp_recovery_gls_ar1)
# plot(fitted(lm_kelp_recovery_gls_ar1), residuals(lm_kelp_recovery_gls_ar1))

# model checks
# normal model
# check_convergence(lm_kelp_recovery_lmer)
# check_normality(lm_kelp_recovery_lmer)
# check_heteroscedasticity(lm_kelp_recovery_lmer)

# plot ACF
# normal model
# acf(residuals(lm_kelp_recovery_lmer))
# pacf(residuals(lm_kelp_recovery_lmer))

# ARMA 2
# acf(residuals(lm_kelp_recovery_lme_ar2))
# pacf(residuals(lm_kelp_recovery_lme_ar2))

# raw kelp biomass model
plot(simulateResiduals(lm_kelp_recovery_zigamma_01))
performance::check_model(lm_kelp_recovery_zigamma_01)

plot(simulateResiduals(lm_kelp_recovery_zigamma_02))
performance::check_model(lm_kelp_recovery_zigamma_02)

# Rsquared
# MuMIn::r.squaredGLMM(lm_kelp_recovery_lmer)
# MuMIn::r.squaredGLMM(lm_kelp_recovery_lme_ar1)
MuMIn::r.squaredGLMM(lm_kelp_recovery_zigamma_01)
MuMIn::r.squaredGLMM(lm_kelp_recovery_zigamma_02)

# summary
# summary(lm_kelp_recovery_lmer)
# summary(lm_kelp_recovery_lme_ar1)
# summary(lm_kelp_recovery_gls_ar1)
# lm_kelp_recovery_summary <- lm_kelp_recovery_lmer %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**",
#     df = "**df**"
#   ) 
# lm_kelp_recovery_summary

summary(lm_kelp_recovery_zigamma_02)
parameters::model_parameters(lm_kelp_recovery_zigamma_02, exponentiate = TRUE)

lm_kelp_recovery_zigamma_summary <- lm_kelp_recovery_zigamma_02 %>% 
  tbl_regression(intercept = TRUE) %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Estimate**"
  ) %>% 
  modify_column_indent(
    columns = label, 
    rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))

# filter out zero-inflated component
lm_kelp_recovery_zigamma_summary$table_body <- lm_kelp_recovery_zigamma_summary$table_body %>% 
  filter(component != "zi")
# change labels
lm_kelp_recovery_zigamma_summary$table_body$label <- c(
  `(Intercept)` = "(Intercept)",
  time_since_end = "Time since end",
  treatmentremoval = "Treatment (removal)",
  `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
)

# final table 
lm_kelp_recovery_zigamma_summary


# AIC comparisons
# MuMIn::AICc(lm_kelp_recovery_lmer, lm_kelp_recovery_lme_ar1, lm_kelp_recovery_lme_ar2, lm_kelp_recovery_gls_ar1) %>% 
#   arrange(AICc)
# gls best model when carp is not taken out

# MuMIn::AICc(lm_kelp_recovery_zigamma_01, lm_kelp_recovery_zigamma_02)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# # all sites
# predicted_kelp_after_overall <- ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end", type = "fixed")
# # predicted line crosses 0 at 4.6
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [4:5 by = 0.01]", type = "fixed")
# # upper bound crosses 0 at 3.1
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3:4 by = 0.01]", type = "fixed")
# # lower bound crosses 0 at 6.8
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [6.7:7 by = 0.001]", type = "fixed")
# 
# # aque: 4.2 years
# predicted_kelp_after_aque <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
# # predicted time to recovery: 4.2 years
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [4:5 by = 0.001]", type = "random", condition = c(site = "aque")) 
# 
# # napl: 4.1 years
# predicted_kelp_after_napl <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
# # predicted time to recovery: 4.1 years
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [4:5 by = 0.001]", type = "random", condition = c(site = "napl"))
# 
# # mohk: 5.9 years
# predicted_kelp_after_mohk <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
# # predicted time to recovery: 5.9 years
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [5.5:6.5 by = 0.01]", type = "random", condition = c(site = "mohk"))
# 
# # carp: 3.9 years
# predicted_kelp_after_carp <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))
# # predicted time to recovery: 3.9 years
# ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.5:4 by = 0.01]", type = "random", condition = c(site = "carp"))

# raw kelp biomass
predicted_kelp_after_raw <- ggpredict(lm_kelp_recovery_zigamma_02,
                                      terms = c("time_since_end [0:6.75 by = 0.25]", "treatment"), type = "fixed")
# kelp in treatment plot matches kelp in reference plot in about 4 years

# aque
predicted_kelp_after_aque <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "aque"))
# predicted time to recovery: 4 years
ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.001]", "treatment"), type = "random", condition = c(site = "aque")) %>% plot()

# napl
predicted_kelp_after_napl <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "napl"))
# predicted time to recovery: 4 years
ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.001]", "treatment"), type = "random", condition = c(site = "napl")) %>% plot()

# mohk
predicted_kelp_after_mohk <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "mohk")) 
# predicted time to recovery: 4 years
ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.01]", "treatment"), type = "random", condition = c(site = "mohk")) %>% plot()

# carp
predicted_kelp_after_carp <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "carp"))
# predicted time to recovery: 4 years
ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.01]", "treatment"), type = "random", condition = c(site = "carp")) %>% plot()

# ⟞ ⟞ iii. means ----------------------------------------------------------

# estimates predicted mean values of response for each level in predictor
recovery_response_means <- estimate_means(lm_kelp_recovery_zigamma_02) %>% 
  mutate(exp_dates = "after")
recovery_response_emmeans <- emmeans(lm_kelp_recovery_zigamma_02, ~ treatment | time_since_end, type = "response") %>% plot()

# calculates differences between levels of categorical predictor
recovery_contrasts <- modelbased::estimate_contrasts(lm_kelp_recovery_zigamma_02, transform = "response") 
# when transforming to response scale, gives ratio of means (weird?)
recovery_contrasts

# estimates slopes of numeric predictors for levels of categorical predictors
# in summary object: estimate = difference in slope of time since end removal - reference
recovery_slopes <- modelbased::estimate_slopes(lm_kelp_recovery_zigamma_02, 
                                               trend = "time_since_end",
                                               at = "treatment") %>% 
  as_tibble() %>% 
  mutate(exp_dates = "Recovery") %>% 
  rbind(during_slopes) %>% 
  mutate(exp_dates = fct_relevel(exp_dates, "Experimental removal", "Recovery"))

recovery_emtrends <- emmeans::emtrends(lm_kelp_recovery_zigamma_02, "treatment", var = "time_since_end") 

# ggsave(filename = here::here("figures", "ms-figures", paste("fig-S11_", today(), ".jpg", sep = "")),
#        effplot,
#        width = 6, height = 4, dpi = 150)

# ⟞ c. figures -------------------------------------------------------------

# ⟞ ⟞ i. overall model predictions -----------------------------------------

# overall_ms <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_continual, 
#              aes(x = time_since_end, y = delta_continual, fill = site, shape = site), size = 2, alpha = 0.9) +
#   scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
#   scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
#   # new_scale("color") + 
#   # overall
#   geom_line(data = predicted_kelp_after_overall, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_after_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   geom_line(data = predicted_kelp_during_overall, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_during_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(breaks = seq(-1500, 1000, by = 1000), limits = c(-1800, 900)) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         # plot.margin = margin(0, 0, 0, 0),
#         legend.position = c(0.1, 0.858),
#         legend.key.size = unit(0.3, units = "cm")) +
#   labs(x = "Time since end of removal (years)", 
#        y = expression("\U0394"~giant~kelp~biomass~"(removal - reference, "~dry~g/m^{"2"}~")"), 
#        fill = "Site",
#        shape = "Site")


# ⟞ ⟞ new model -----------------------------------------------------------

overall_kelp <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  
  # raw data points
  geom_point(data = continual_long, aes(x = time_since_end, y = kelp_biomass, color = treatment),
             alpha = 0.15,
             size = 0.75,
             shape = 21) +
  
  # prediction lines
  geom_line(data = predicted_kelp_during_raw, aes(x = x, y = predicted, color = group), linewidth = 1) +
  geom_line(data = predicted_kelp_after_raw, aes(x = x, y = predicted, color = group), linewidth = 1) +
  
  # confidence intervals
  geom_ribbon(data = predicted_kelp_during_raw, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
  geom_ribbon(data = predicted_kelp_after_raw, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
  
  # colors and shapes
  scale_color_manual(values = c(reference = reference_col, removal = removal_col),
                     labels = c("Reference", "Removal")) +
  scale_linetype_manual(values = c(reference = 2, removal = 1),
                     labels = c("Reference", "Removal")) +
  
  # removal/recovery labels
  annotate(geom = "text", x = -6.75, y = 1925, label = "Removal", size = 3) +
  annotate(geom = "text", x = 5.5, y = 1925, label = "Recovery", size = 3) +
  
  theme_bw() + 
  scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
  coord_cartesian(ylim = c(-10, 2000)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.background = element_blank(),
        # plot.margin = margin(0, 0, 0, 0),
        # legend.position = c(0.85, 0.9),
        legend.position = "none",
        legend.key.size = unit(0.5, units = "cm"),
        legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
        legend.spacing.y = unit(0.1, units = "cm"),
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
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
       title = "(a)")

overall_kelp

# napl_raw <- ggplot(data = continual_long %>% filter(site == "napl"),
#        aes(x = time_since_end, y = kelp_biomass, color = treatment)) +
#   # x at 0 and y at 0 lines
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   
#   # raw data points
#   geom_point(shape = 1, size = 1) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~exp_dates)
# 
# napl_delta <- ggplot(data = delta_continual %>% filter(site == "napl"),
#        aes(x = time_since_end, y = delta_continual)) +
#   geom_point()

overall_kelp_removal <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
  
  # raw data points
  geom_point(data = continual_long %>% filter(treatment == "removal"), aes(x = time_since_end, y = kelp_biomass), shape = 1, size = 1, alpha = 0.4, color = removal_col) +
  
  # prediction lines
  geom_line(data = predicted_kelp_during_raw %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
  geom_line(data = predicted_kelp_after_raw %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
  
  # confidence intervals
  geom_ribbon(data = predicted_kelp_during_raw %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  geom_ribbon(data = predicted_kelp_after_raw %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  scale_y_continuous(limits = c(-10, 2000)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        # plot.margin = margin(0, 0, 0, 0),
        legend.position = c(0.88, 0.73),
        legend.key.size = unit(0.5, units = "cm"),
        legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
        legend.spacing.y = unit(0.1, units = "cm"),
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1) +
  labs(x = "Time since end of removal (years)", 
       y = "Giant kelp biomass (dry g/m\U00B2)",
       title = "(a) Removal")

overall_kelp_removal

overall_kelp_reference <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
  
  # raw data points
  geom_point(data = continual_long %>% filter(treatment == "reference"), aes(x = time_since_end, y = kelp_biomass), shape = 1, size = 1, alpha = 0.4, color = reference_col) +
  
  # prediction lines
  geom_line(data = predicted_kelp_during_raw %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
  geom_line(data = predicted_kelp_after_raw %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
  
  # confidence intervals
  geom_ribbon(data = predicted_kelp_during_raw %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  geom_ribbon(data = predicted_kelp_after_raw %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  scale_y_continuous(limits = c(-10, 2000)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        # plot.margin = margin(0, 0, 0, 0),
        legend.position = c(0.88, 0.73),
        legend.key.size = unit(0.5, units = "cm"),
        legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
        legend.spacing.y = unit(0.1, units = "cm"),
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1) +
  labs(x = "Time since end of removal (years)", 
       y = "Giant kelp biomass (dry g/m\U00B2)",
       title = "(b) Reference")

overall_kelp_reference

# ⟞ ⟞ ii. site level predictions -------------------------------------------

# aque <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_continual %>% filter(site == "aque"), aes(x = time_since_end, y = delta_continual), shape = aque_shape, fill = aque_col, size = 2, alpha = 0.9) +
#   # during
#   geom_line(data = predicted_kelp_during_aque, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_during_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   # after
#   geom_line(data = predicted_kelp_after_aque, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_after_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1250, 1100), expand = c(0, 0)) +
#   # geom_text(aes(x = -6.6, y = 600), label = "aque", size = 8) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 8),
#         plot.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         plot.title.position = "plot") +
#   labs(x = "Time since end of removal", 
#        y = "\U0394 giant kelp biomass",
#        title = paste("(a) ", aque_full, sep = ""))
# 
# napl <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_continual %>% filter(site == "napl"), aes(x = time_since_end, y = delta_continual), shape = napl_shape, fill = napl_col, size = 2, alpha = 0.9) +
#   # during
#   geom_line(data = predicted_kelp_during_napl, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_during_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   # after
#   geom_line(data = predicted_kelp_after_napl, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_after_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1250, 1100), expand = c(0, 0)) +
#   # geom_text(aes(x = -6.8, y = 700), label = "napl", size = 8) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 8),
#         plot.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         plot.title.position = "plot") +
#   labs(x = "Time since end of removal", 
#        y = "\U0394 giant kelp biomass",
#        title = paste("(b) ", napl_full, sep = ""))
# 
# mohk <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_continual %>% filter(site == "mohk"), aes(x = time_since_end, y = delta_continual), shape = mohk_shape, fill = mohk_col, size = 2, alpha = 0.9) +
#   # during
#   geom_line(data = predicted_kelp_during_mohk, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_during_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   # after
#   geom_line(data = predicted_kelp_after_mohk, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_after_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1800, 1100), expand = c(0, 0)) +
#   # geom_text(aes(x = -6.6, y = 900), label = "mohk", size = 8) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 8),
#         plot.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         plot.title.position = "plot") +
#   labs(x = "Time since end of removal", 
#        y = "\U0394 giant kelp biomass",
#        title = paste("(c) ", mohk_full, sep = ""))
# 
# carp <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_continual %>% filter(site == "carp"), aes(x = time_since_end, y = delta_continual), shape = carp_shape, fill = carp_col, size = 2, alpha = 0.9) +
#   # during
#   geom_line(data = predicted_kelp_during_carp, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_during_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   # after
#   geom_line(data = predicted_kelp_after_carp, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_kelp_after_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1250, 1100), expand = c(0, 0)) +
#   # geom_text(aes(x = -6.8, y = 1000), label = "carp", size = 8) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         plot.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         plot.title.position = "plot") +
#   labs(x = "Time since end of removal", 
#        y = "\U0394 giant kelp biomass",
#        title = paste("(d) ", carp_full, sep = ""))
# 
# plots_together_sites <- (aque + napl) / (mohk + carp)
# plots_together_sites


# ⟞ ⟞ iii. delta from predictions -----------------------------------------

# data frame of predictions
delta_predictions_during <- predicted_kelp_during_raw %>% 
  as.data.frame() %>% 
  select(x, group, predicted) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta = removal - reference) %>% 
  mutate(exp_dates = "during")

delta_predictions_after <- predicted_kelp_after_raw %>% 
  as.data.frame() %>% 
  select(x, group, predicted) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta = removal - reference) %>% 
  mutate(exp_dates = "after")

overall_predictions <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
  annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(data = delta_continual,
             aes(x = time_since_end, y = delta_continual), 
             shape = 2, 
             alpha = 0.15,
             size = 0.75) +

  # overall
  geom_line(data = delta_predictions_during, aes(x = x, y = delta), linewidth = 1) +
  geom_line(data = delta_predictions_after, aes(x = x, y = delta), linewidth = 1) +

  scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1800, 1000)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        # plot.margin = margin(0, 0, 0, 0),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
  labs(x = "Time since end of removal (years)", 
       y = "\U0394 biomass\n(removal \U2212 reference, dry g/m\U00B2)",
       fill = "Site",
       shape = "Site",
       title = "(b) Removal \U2212 reference")

overall_predictions

fig1_new <- overall_kelp / overall_predictions
# fig1_new_v2 <- overall_kelp_removal / overall_kelp_reference / overall_predictions


# ⟞ ⟞ iv. effect plots ----------------------------------------------------

kelp_means_df <- rbind(during_response_means, recovery_response_means) %>% 
  as_tibble() %>% 
  mutate(exp_dates = fct_relevel(exp_dates, "during", "after"))

# mean predicted responses: what are the differences in mean giant kelp biomass between treatments during each time period?
kelp_means_plot <- continual_long %>% 
  # filter(exp_dates == "during") %>% 
  ggplot(aes(x = treatment, y = kelp_biomass, color = treatment)) +
  facet_wrap(~exp_dates, 
             labeller = labeller(exp_dates = c("during" = "Experimental removal",
                                               "after" = "Recovery"))) +
  geom_point(position = position_jitter(width = 0.1, seed = 666),
             alpha = 0.3, shape = 21) +
  geom_pointrange(data = kelp_means_df, 
                  aes(x = treatment, y = Mean, ymin = CI_low, ymax = CI_high),
                  size = 0.5, linewidth = 0.5) +
  scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
  scale_fill_manual(values = c(reference = reference_col, removal = removal_col)) +
  scale_x_discrete(labels = c("reference" = "Reference", "removal" = "Removal")) +
  scale_y_continuous(limits = c(-5, 2000), expand = c(0.01, 0.01)) +
  labs(x = "Treatment", 
       y = "Giant kelp biomass (dry g/m\U00B2)") +
  theme_classic() +
  theme(legend.position = "none",
        strip.placement = "outer",
        strip.text = element_text(hjust = 0.5),
        strip.background = element_blank()) 

kelp_means_plot

kelp_means_nodata_plot <- ggplot(data = kelp_means_df,
                            aes(x = treatment, y = kelp_biomass, color = treatment)) +
  facet_wrap(~exp_dates, 
             labeller = labeller(exp_dates = c("during" = "Experimental removal",
                                               "after" = "Recovery"))) +
  geom_pointrange(aes(x = treatment, y = Mean, ymin = CI_low, ymax = CI_high),
                  size = 0.5, linewidth = 0.5) +
  scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
  scale_fill_manual(values = c(reference = reference_col, removal = removal_col)) +
  scale_x_discrete(labels = c("reference" = "Reference", "removal" = "Removal")) +
  scale_y_continuous(limits = c(-5, 550), 
                     breaks = seq(from = 0, to = 550, by = 125),
                     expand = c(0.01, 0.01)) +
  labs(x = "Treatment", 
       y = "Predicted giant kelp biomass \n(dry g/m\U00B2)") +
  theme_classic() +
  theme(legend.position = "none",
        strip.placement = "outer",
        strip.text = element_text(hjust = 0.5),
        strip.background = element_blank()) 
kelp_means_nodata_plot

# mean predicted slopes: what are the differences in rates of giant kelp biomass change?
effplot <- ggplot(data = recovery_slopes, aes(x = treatment, y = Coefficient)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high)) +
  labs(x = "Treatment", 
       y = "Effect of time since end (years)") +
  theme_bw() +
  raw_biomass_plot_theme() +
  facet_wrap(~exp_dates)
effplot

# ⟞ ⟞ v. conditional means ------------------------------------------------

predicted_kelp_after_0_raw <- ggpredict(lm_kelp_recovery_zigamma_02,
                                      terms = c("time_since_end [0]", "treatment"), type = "fixed")

predicted_kelp_after_both_raw <- ggpredict(lm_kelp_recovery_zigamma_02,
                                        terms = c("time_since_end [4]", "treatment"), 
                                        type = "fixed") %>% 
  rbind(predicted_kelp_after_0_raw) %>% 
  rename("treatment" = group, 
         "kelp_biomass" = predicted,
         "time_since_end" = x)

means <- continual_long %>% 
  filter(time_since_end == 0 | time_since_end == 4) %>% 
  ggplot(aes(x = treatment, y = kelp_biomass, color = treatment)) +
  geom_point(position = position_jitter(width = 0.1, seed = 1),
             shape = 21, alpha = 0.8, size = 1) +
  geom_pointrange(data = predicted_kelp_after_both_raw,
                  aes(ymin = conf.low, ymax = conf.high)) +
  scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
  scale_x_discrete(labels = c(reference = "Reference", removal = "Removal")) +
  labs(x = "Treatment",
       y = "Giant kelp biomass (dry g/m\U00B2)") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  #geom_pointrange(data = predicted_kelp_after_0, ) +
  facet_wrap(~time_since_end, labeller = labeller(time_since_end = c("0" = "(a) Time since end = 0", "4" = "(b) Time since end = 4")))
means


##########################################################################-
# 4. variation plots ------------------------------------------------------
##########################################################################-

# variation_during <- aov(variation ~ treatment, data = continual_long %>% filter(exp_dates == "during"))
# summary(variation_during)
# var_during_df <- ggpredict(variation_during, terms = c("treatment")) %>% 
#   as_tibble() %>% 
#   mutate(exp_dates = "during")
# variation_after <- aov(variation ~ treatment, data = continual_long %>% filter(exp_dates == "after"))
# summary(variation_after)
# var_after_df <- ggpredict(variation_after, terms = c("treatment")) %>% 
#   as_tibble() %>% 
#   mutate(exp_dates = "after")
# 
# var_df <- rbind(var_during_df, var_after_df) %>% 
#   rename(treatment = x) %>% 
#   mutate(exp_dates = fct_relevel(exp_dates, "during", "after"))
# 
# ggplot(data = var_df, aes(x = treatment, y = predicted)) +
#   geom_point(data = continual_long, aes(x = treatment, y = variation),
#              alpha = 0.2, shape = 21,
#              position = position_jitter(width = 0.2)) +
#   geom_pointrange(aes(x = treatment, y = predicted, ymin = conf.low, ymax = conf.high)) +
#   facet_wrap(~exp_dates) +
#   labs(x = "Treatment", 
#        y = "Variation") +
#   theme_bw() +
#   theme(panel.grid = element_blank())
# 
# ggplot(data = var_df, aes(x = treatment, y = predicted)) +
#   geom_pointrange(aes(x = treatment, y = predicted, ymin = conf.low, ymax = conf.high)) +
#   facet_wrap(~exp_dates) +
#   labs(x = "Treatment", 
#        y = "Variation") +
#   theme_bw() +
#   theme(panel.grid = element_blank())

variation_site <- continual_long %>% 
  select(exp_dates, treatment, site, kelp_biomass) %>% 
  group_by(exp_dates, treatment, site) %>% 
  summarize(mean = mean(kelp_biomass),
            variation = sd(kelp_biomass)/mean(kelp_biomass)) %>% 
  ungroup()

variation_plot <- ggplot(variation_site, aes(x = treatment, y = variation, color = treatment)) +
  geom_point(position = position_jitter(width = 0.1, seed = 666),
             alpha = 0.8, shape = 21) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
  scale_x_discrete(labels = c(reference = "Reference", removal = "Removal")) +
  facet_wrap(~exp_dates, labeller = labeller(exp_dates = c("during" = "(a) Experimental removal", "after" = "(b) Recovery period"))) +
  labs(x = "Treatment",
       y = "Coefficient of variation") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") 

variation_plot

##########################################################################-
# 5. manuscript tables ----------------------------------------------------
##########################################################################-

# lm_kelp_tables <- tbl_merge(tbls = list(lm_kelp_during_summary, lm_kelp_recovery_summary), 
#                             tab_spanner = c("**Removal**", "**Recovery**")) 
# this table is compiled with others in the `02a-community_recovery.R` script

lm_kelp_zigamma_tables <- tbl_merge(tbls = list(lm_kelp_during_zigamma_summary, lm_kelp_recovery_zigamma_summary),
                                    tab_spanner = c("**Experimental removal**", "**Recovery**"))

##########################################################################-
# 6. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. delta kelp through time -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_", today(), ".jpg", sep = "")),
#        plot = overall_ms,
#        height = 8, width = 14, units = "cm",
#        dpi = 400)

# ⟞ b. raw kelp biomass through time --------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S1_", today(), ".jpg", sep = "")),
#        plot = continual_sites_raw,
#        height = 12, width = 8, units = "cm",
#        dpi = 300)

# ⟞ c. recovery time vs biomass -------------------------------------------

# s4_panels <- rec_time_plot + delta_vs_biomass +
#   plot_layout(widths = c(2, 2.1)) +
#   plot_annotation(tag_levels = "a", tag_suffix = ")") &
#   theme(plot.tag = element_text(size = 12))
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_", today(), ".jpg", sep = "")),
#        plot = rec_time_plot,
#        height = 8, width = 12, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_panels_", today(), ".jpg", sep = "")),
#        plot = s4_panels,
#        height = 8, width = 16, units = "cm",
#        dpi = 400)

# ⟞ d. predictions by site ------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S2_", today(), ".jpg", sep = "")),
#        plot = plots_together_sites,
#        height = 10, width = 16, units = "cm",
#        dpi = 300)


# ⟞ e. new model ----------------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_new-model_", today(), ".jpg", sep = "")),
#        plot = fig1_new,
#        height = 17, width = 13, units = "cm",
#        dpi = 300)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_new-model_v2_", today(), ".jpg", sep = "")),
#        plot = fig1_new_v2,
#        height = 17, width = 13, units = "cm",
#        dpi = 300)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_new-model_removal", today(), ".jpg", sep = "")),
#        plot = overall_kelp_removal,
#        height = 8, width = 14, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-1_new-model_reference", today(), ".jpg", sep = "")),
#        plot = overall_kelp_reference,
#        height = 8, width = 14, units = "cm",
#        dpi = 300)

# Ecology max figure size: 18 x 24

# ⟞ f. means plots --------------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("kelp_means_plot_", today(), ".jpg", sep = "")),
#        plot = kelp_means_plot,
#        height = 7, width = 14, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("kelp_means_nodata_plot_", today(), ".jpg", sep = "")),
#        plot = kelp_means_nodata_plot,
#        height = 7, width = 14, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("kelp-means-recovery-plot_", today(), ".jpg", sep = "")),
#        plot = means,
#        height = 7, width = 14, units = "cm",
#        dpi = 300)

# ⟞ g. kelp density and fronds --------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("kelp-density-timeseries_", today(), ".jpg", sep = "")),
#        plot = density_timeseries,
#        height = 9, width = 16, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("kelp-frond-timeseries_", today(), ".jpg", sep = "")),
#        plot = fronds_timeseries,
#        height = 9, width = 16, units = "cm",
#        dpi = 300)

# ⟞ h. variation plots --------------------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("cov_plot_", today(), ".jpg", sep = "")),
#        plot = variation_plot,
#        height = 7, width = 14, units = "cm",
#        dpi = 300)


