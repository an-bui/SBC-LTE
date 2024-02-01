
##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "01a-kelp_recovery.R"))


##########################################################################-
# 1. data frames ----------------------------------------------------------
##########################################################################-

# ⟞ a. algae -------------------------------------------------------------

# total biomass
algae_biomass <- biomass %>% 
  filter(new_group == "algae" & sp_code != "MAPY")

# delta biomass
delta_algae_biomass <- algae_biomass %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  group_by(site, year, month, treatment, date) %>% 
  summarize(total_dry = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = treatment, values_from = total_dry) %>% 
  mutate(delta_annual = annual - control,
         delta_continual = continual - control) 

# joined with kelp deltas and raw biomass
delta_algae_continual <- delta_algae_biomass %>% 
  dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for annual removal:
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
  # create a new sample ID that is site, year, quarter
  unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
  # create new column that is the column to join with 
  # unite("join_ID", site, date, remove = FALSE) %>% 
  rename("control_algae" = control,
         "continual_algae" = continual,
         "delta_continual_algae" = delta_continual) %>% 
  full_join(., delta_continual %>% dplyr::select(sample_ID, continual, control, delta_continual), by = "sample_ID") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))

# long format
algae_continual_long <- delta_algae_biomass %>% 
  dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  select(!delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for annual removal:
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
  # make it longer
  pivot_longer(cols = c(control, continual)) %>% 
  # rename columns
  rename(treatment = name, algae_biomass = value) %>% 
  # change treatment names
  mutate(treatment = case_match(treatment, "control" ~ "reference", "continual" ~ "removal")) %>% 
  # create a new sample ID that is site, year, quarter, treatment
  unite("sample_ID", site, date, quarter, treatment, remove = FALSE) %>% 
  # join only with kelp biomass (long format)
  left_join(., select(.data = continual_long, sample_ID, kelp_biomass), by = "sample_ID")

# ⟞ b. epilithic invertebrates -------------------------------------------

# total biomass
epi_biomass <- biomass %>% 
  filter(new_group == "epilithic.sessile.invert") %>% 
  # # take out endos
  filter(taxon_family != "Pholadidae")

# delta biomass
delta_epi_biomass <- epi_biomass %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  group_by(site, year, month, treatment, date) %>% 
  summarize(total_dry = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = treatment, values_from = total_dry) %>% 
  mutate(delta_annual = annual - control,
         delta_continual = continual - control) 

# join with kelp deltas and raw biomass
delta_epi_continual <- delta_epi_biomass %>% 
  dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for annual removal:
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
  # create a new sample ID that is site, year, quarter
  unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
  # create new column that is the column to join with 
  # unite("join_ID", site, date, remove = FALSE) %>% 
  rename("control_epi" = control,
         "continual_epi" = continual,
         "delta_continual_epi" = delta_continual) %>% 
  full_join(., delta_continual %>% dplyr::select(sample_ID, continual, control, delta_continual), by = "sample_ID") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))

# long format
epi_continual_long <- delta_epi_biomass %>% 
  dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
  # take out years where continual removal hadn't happened yet
  drop_na(delta_continual) %>% 
  select(!delta_continual) %>% 
  mutate(exp_dates = case_when(
    # after for annual removal:
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
  # make it longer
  pivot_longer(cols = c(control, continual)) %>% 
  # rename columns
  rename(treatment = name, epi_biomass = value) %>% 
  # change treatment names
  mutate(treatment = case_match(treatment, "control" ~ "reference", "continual" ~ "removal")) %>% 
  # create a new sample ID that is site, year, quarter, treatment
  unite("sample_ID", site, date, quarter, treatment, remove = FALSE) %>% 
  # join only with kelp biomass (long format)
  left_join(., select(.data = continual_long, sample_ID, kelp_biomass), by = "sample_ID")

# ⟞ c. endolithic invertebrates ------------------------------------------

# total biomass
# endo_biomass <- biomass %>% 
#   filter(taxon_family == "Pholadidae")

# delta biomass
# delta_endo_biomass <- endo_biomass %>% 
#   dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
#   group_by(site, year, month, treatment, date) %>% 
#   summarize(total_dry = sum(dry_gm2)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = treatment, values_from = total_dry) %>% 
#   mutate(delta_annual = annual - control,
#          delta_continual = continual - control) 

# joining with kelp deltas and biomass
# delta_endo_continual <- delta_endo_biomass %>% 
#   dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
#   # take out years where continual removal hadn't happened yet
#   drop_na(delta_continual) %>% 
#   mutate(exp_dates = case_when(
#     # after for annual removal:
#     site == "aque" & date >= aque_after_date_continual ~ "after",
#     site == "napl" & date >= napl_after_date_continual ~ "after",
#     site == "mohk" & date >= mohk_after_date_continual ~ "after",
#     site == "carp" & date >= carp_after_date_continual ~ "after",
#     # everything else is "during" the experiment
#     TRUE ~ "during"
#   ),
#   exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
#   time_since_columns_continual() %>% 
#   kelp_year_column() %>% 
#   comparison_column_continual() %>% 
#   # create a new sample ID that is site, year, quarter
#   unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
#   # create new column that is the column to join with 
#   # unite("join_ID", site, date, remove = FALSE) %>% 
#   rename("control_endo" = control,
#          "continual_endo" = continual,
#          "delta_continual_endo" = delta_continual) %>% 
#   full_join(., delta_continual %>% dplyr::select(sample_ID, continual, control, delta_continual), by = "sample_ID") %>% 
#   left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
#   rename("site_full" = value) %>% 
#   mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
#   mutate(site = fct_relevel(site, "aque", "napl", "mohk", "carp"))

##########################################################################-
# 2. timeseries plots -----------------------------------------------------
##########################################################################-

# ⟞ a. raw biomass --------------------------------------------------------

# ⟞ ⟞ i. algae -----------------------------------------------------------

delta_continual_sites_algae_raw <- delta_algae_continual %>% 
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
  geom_line(aes(x = time_since_end, y = control_algae, col = site), alpha = 0.5, linewidth = 2) +
  geom_point(aes(x = time_since_end, y = control_algae, shape = site), size = 1, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual_algae, col = site), linewidth = 2) +
  geom_point(aes(x = time_since_end, y = continual_algae, shape = site, col = site), size = 1, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  raw_biomass_plot_theme() +
  labs(x = "Time since end of experiment (years)", 
       y = expression(Understory~macroalgae~biomass~(dry~g/m^{"2"}))) +
  facet_wrap(~strip, scales = "free_y")
delta_continual_sites_algae_raw

# ⟞ ⟞ ii. epi inverts ----------------------------------------------------

delta_continual_sites_epi_raw <- delta_epi_continual %>% 
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
  geom_line(aes(x = time_since_end, y = control_epi, col = site), alpha = 0.5, linewidth = 2) +
  geom_point(aes(x = time_since_end, y = control_epi, shape = site), size = 1, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual_epi, col = site), linewidth = 2) +
  geom_point(aes(x = time_since_end, y = continual_epi, shape = site, col = site), size = 1, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  raw_biomass_plot_theme() +
  labs(x = "Time since end of experiment (years)", 
       y = expression(Epilithic~invertebrate~biomass~(dry~g/m^{"2"}))) +
  facet_wrap(~strip, scales = "free_y")
delta_continual_sites_epi_raw

# ⟞ ⟞ iii. endo inverts --------------------------------------------------

# delta_continual_sites_endo_raw <- delta_endo_continual %>% 
#   mutate(strip = case_when(
#     site == "aque" ~ paste("(a) ", site_full, sep = ""),
#     site == "napl" ~ paste("(b) ", site_full, sep = ""),
#     site == "mohk" ~ paste("(c) ", site_full, sep = ""),
#     site == "carp" ~ paste("(d) ", site_full, sep = "")
#   )) %>% 
#   ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_line(aes(x = time_since_end, y = control_endo, col = site), alpha = 0.5, linewidth = 2) +
#   # control
#   geom_point(aes(x = time_since_end, y = control_endo, shape = site), size = 1, alpha = 0.5, fill = "#FFFFFF") +
#   # continual
#   geom_line(aes(x = time_since_end, y = continual_endo, col = site), linewidth = 2) +
#   geom_point(aes(x = time_since_end, y = continual_endo, shape = site, col = site), size = 1, fill = "#FFFFFF") +
#   scale_shape_manual(values = shape_palette_site) +
#   scale_color_manual(values = color_palette_site) +
#   scale_fill_manual(values = color_palette_site) +
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   raw_biomass_plot_theme() +
#   labs(x = "Time since end of experiment (years)", 
#        y = expression(Endolithic~invertebrate~biomass~(dry~g/m^{"2"}))) +
#   facet_wrap(~strip, scales = "free_y")
# 
# delta_continual_sites_endo_raw



##########################################################################-
# 3. algae linear model ---------------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
# lm_algae_during_lmer <- lmer(
#   delta_continual_algae ~ time_since_end + (1|site), 
#   data = delta_algae_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass)

# lm_raw_algae_during_zigamma_01 <- glmmTMB(
#   algae_biomass ~ time_since_end*treatment + (1|site), 
#   data = algae_continual_long %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1)

lm_raw_algae_during_zigamma_02 <- glmmTMB(
  algae_biomass ~ time_since_end*treatment + (1|site) + (1|year), 
  data = algae_continual_long %>% filter(exp_dates == "during"), 
  na.action = na.pass,
  family = ziGamma(link = "log"),
  ziformula = ~1)


# diagnostics
# plot(simulateResiduals(lm_algae_during_lmer))
# check_model(lm_algae_during_lmer)

# plot(simulateResiduals(lm_raw_algae_during_zigamma_01))
# check_convergence(lm_raw_algae_during_zigamma_01)

plot(simulateResiduals(lm_raw_algae_during_zigamma_02))

# Rsquared
# r.squaredGLMM(lm_algae_during_lmer)
# r.squaredGLMM(lm_raw_algae_during_zigamma_01)
r.squaredGLMM(lm_raw_algae_during_zigamma_02)

# summary
# summary(lm_algae_during_lmer)
# lm_algae_during_summary <- lm_algae_during_lmer %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**",
#     df = "**df**"
#   ) 
# lm_algae_during_summary

lm_raw_algae_during_zigamma_summary <- lm_raw_algae_during_zigamma_02 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Beta**"
  ) %>% 
  modify_column_indent(
    columns = label, 
    rows = variable %in% c("treatment", "time_since_end", "time_since_end:treatment"))

# filter out zero-inflated component
lm_raw_algae_during_zigamma_summary$table_body <- lm_raw_algae_during_zigamma_summary$table_body %>% 
  filter(component != "zi")
# change labels
lm_raw_algae_during_zigamma_summary$table_body$label <- c(
  time_since_end = "Time since end",
  treatmentremoval = "Treatment (removal)",
  `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
)

# final table 
lm_raw_algae_during_zigamma_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# predicted_algae_during <- ggpredict(lm_algae_during_lmer, terms = ~ time_since_end, type = "fixed")

predicted_raw_algae_during <- ggpredict(lm_raw_algae_during_zigamma_02, terms = c("time_since_end", "treatment"), type = "fixed")

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
# lm_algae_recovery_lmer <- lmer(
#   delta_continual_algae ~ time_since_end + (1|site), 
#   data = delta_algae_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass)

# lm_raw_algae_recovery_zigamma_01 <- glmmTMB(
#   algae_biomass ~ time_since_end*treatment + (1|site),
#   data = algae_continual_long %>% filter(exp_dates == "after"),
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1)

lm_raw_algae_recovery_zigamma_02 <- glmmTMB(
  algae_biomass ~ time_since_end*treatment + (1|site) + (1|year), 
  data = algae_continual_long %>% filter(exp_dates == "after"), 
  na.action = na.pass,
  family = ziGamma(link = "log"),
  ziformula = ~1)

# diagnostics
# plot(simulateResiduals(lm_algae_recovery_lmer))
# check_model(lm_algae_recovery_lmer)

# plot(simulateResiduals(lm_raw_algae_recovery_zigamma_01))

plot(simulateResiduals(lm_raw_algae_recovery_zigamma_02))

# Rsquared
# r.squaredGLMM(lm_algae_recovery_lmer)
# r.squaredGLMM(lm_raw_algae_recovery_zigamma_01)
r.squaredGLMM(lm_raw_algae_recovery_zigamma_02)

# summary
# summary(lm_algae_recovery_lmer)
# lm_algae_recovery_summary <- lm_algae_recovery_lmer %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**",
#     df = "**df**"
#   ) 
# lm_algae_recovery_summary

lm_raw_algae_recovery_zigamma_summary <- lm_raw_algae_recovery_zigamma_02 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Beta**"
  ) %>% 
  modify_column_indent(
    columns = label, 
    rows = variable %in% c("treatment", "time_since_end", "time_since_end:treatment"))

# filter out zero-inflated component
lm_raw_algae_recovery_zigamma_summary$table_body <- lm_raw_algae_recovery_zigamma_summary$table_body %>% 
  filter(component != "zi")
# change labels
lm_raw_algae_recovery_zigamma_summary$table_body$label <- c(
  time_since_end = "Time since end",
  treatmentremoval = "Treatment (removal)",
  `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
)

# final table 
lm_raw_algae_recovery_zigamma_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# predicted_algae_recovery <- ggpredict(lm_algae_recovery_lmer, terms = ~time_since_end, type = "fixed")

predicted_raw_algae_recovery <- ggpredict(lm_raw_algae_recovery_zigamma_02, terms = c("time_since_end[0:6.75, by = 0.25]", "treatment"), type = "fixed")

# ⟞ c. figure ------------------------------------------------------------

# algae_time <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_algae_continual, 
#              aes(x = time_since_end, y = delta_continual_algae, fill = site, shape = site), size = 2, alpha = 0.9) +
#   scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
#   scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
#   # new_scale("color") + 
#   # overall
#   geom_line(data = predicted_algae_recovery, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_algae_recovery, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   geom_line(data = predicted_algae_during, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_algae_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
#   delta_timeseries_theme("algae") +
#   labs(x = "Time since end of removal (years)",
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(b)",
#        fill = "Site", shape = "Site")
# 
# algae_time


# ⟞ ⟞ new model -----------------------------------------------------------

raw_algae_time <- ggplot() +
  # reference lines
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # raw data
  geom_point(data = algae_continual_long, 
             aes(x = time_since_end, y = algae_biomass, color = treatment, shape = treatment, size = treatment), alpha = 0.2) +
  
  # model predictions
  geom_line(data = predicted_raw_algae_recovery, aes(x = x, y = predicted, lty = group, color = group), linewidth = 1) +
  geom_ribbon(data = predicted_raw_algae_recovery, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
  geom_line(data = predicted_raw_algae_during, aes(x = x, y = predicted, lty = group, color = group), linewidth = 1) +
  geom_ribbon(data = predicted_raw_algae_during, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
  
  # colors and shapes
  scale_color_manual(values = c(reference = reference_col, removal = removal_col),
                     labels = c("Reference", "Removal")) +
  # scale_fill_manual(values = c(reference = "#6D5A1800", removal = "#CC754066"),
  #                    labels = c("Reference", "Removal")) +
  scale_linetype_manual(values = c(reference = 2, removal = 1),
                        labels = c("Reference", "Removal")) +
  scale_shape_manual(values = c(reference = 1, removal = 16),
                     labels = c("Reference", "Removal")) +
  scale_size_manual(values = c(reference = 1, removal = 1.3),
                    labels = c("Reference", "Removal")) +
  
  # theming
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
  coord_cartesian(ylim = c(30, 800)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.position = c(0.86, 0.88),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, units = "cm"),
        legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
        legend.spacing.y = unit(0.01, "cm"),
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

raw_algae_time

raw_algae_removal <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # raw data points
  geom_point(data = algae_continual_long %>% filter(treatment == "removal"), aes(x = time_since_end, y = algae_biomass), shape = 1, size = 1, alpha = 0.4, color = removal_col) +
  
  # prediction lines
  geom_line(data = predicted_raw_algae_during %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
  geom_line(data = predicted_raw_algae_recovery %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
  
  # confidence intervals
  geom_ribbon(data = predicted_raw_algae_during %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  geom_ribbon(data = predicted_raw_algae_recovery %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  scale_y_continuous(limits = c(-10, 800)) +
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
       y = "Biomass (dry g/m\U00B2)",
       title = "(a) Removal")

raw_algae_reference <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # raw data points
  geom_point(data = algae_continual_long %>% filter(treatment == "reference"), aes(x = time_since_end, y = algae_biomass), shape = 1, size = 1, alpha = 0.4, color = reference_col) +
  
  # prediction lines
  geom_line(data = predicted_raw_algae_during %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, color = reference_col) +
  geom_line(data = predicted_raw_algae_recovery %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, color = reference_col) +
  
  # confidence intervals
  geom_ribbon(data = predicted_raw_algae_during %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  geom_ribbon(data = predicted_raw_algae_recovery %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  scale_y_continuous(limits = c(-10, 800)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1) +
  labs(x = "Time since end of reference (years)", 
       y = "Biomass (dry g/m\U00B2)",
       title = "(b) Reference")

# data frame of predictions
delta_algae_predictions_during <- predicted_raw_algae_during %>% 
  as.data.frame() %>% 
  select(x, group, predicted) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta = removal - reference) %>% 
  mutate(exp_dates = "during")

delta_algae_predictions_after <- predicted_raw_algae_recovery %>% 
  as.data.frame() %>% 
  select(x, group, predicted) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta = removal - reference) %>% 
  mutate(exp_dates = "after")

overall_algae_predictions <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_algae_continual,
             aes(x = time_since_end, y = delta_continual_algae), shape = 2, size = 1, alpha = 0.4) +
  
  # overall
  geom_line(data = delta_algae_predictions_during, aes(x = x, y = delta), linewidth = 1) +
  geom_line(data = delta_algae_predictions_after, aes(x = x, y = delta), linewidth = 1) +
  
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-250, 800, by = 200), limits = c(-250, 775)) +
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
       y = "\U0394 biomass \n (removal \U2212 reference, dry g/m\U00B2)",
       fill = "Site",
       shape = "Site",
       title = "(c) Removal \U2212 reference")

overall_algae_predictions

algae_title/raw_algae_removal/raw_algae_reference/overall_algae_predictions +
  plot_layout(heights = c(3, 20, 20, 20))

##########################################################################-
# 4. epi. invert linear model ---------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
# lm_epi_during_lmer <- lmer(
#   delta_continual_epi ~ time_since_end + (1|site), 
#   data = delta_epi_continual %>% filter(exp_dates == "during"), 
#   na.action = na.pass)

# lm_raw_epi_during_lmer <- lmer(
#   epi_biomass ~ time_since_end*treatment + (1|site), 
#   data = epi_continual_long %>% filter(exp_dates == "during"), 
#   na.action = na.pass)

# lm_raw_epi_during_zigamma_01 <- glmmTMB(
#   epi_biomass ~ time_since_end*treatment + (1|site),
#   data = epi_continual_long %>% filter(exp_dates == "during"),
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1
# )

lm_raw_epi_during_zigamma_02 <- glmmTMB(
  epi_biomass ~ time_since_end*treatment + (1|site) + (1|year),
  data = epi_continual_long %>% filter(exp_dates == "during"),
  na.action = na.pass,
  family = ziGamma(link = "log"),
  ziformula = ~1
)

# diagnostics
# plot(simulateResiduals(lm_epi_during_lmer))
# check_model(lm_epi_during_lmer)

# plot(simulateResiduals(lm_raw_epi_during_lmer))
# check_model(lm_raw_epi_during_lmer)

# plot(simulateResiduals(lm_raw_epi_during_zigamma_01))

plot(simulateResiduals(lm_raw_epi_during_zigamma_02))

# R2
# r.squaredGLMM(lm_epi_during_lmer)
r.squaredGLMM(lm_raw_epi_during_zigamma_02)

# summary
# summary(lm_epi_during_lmer) 
# lm_epi_during_summary <- lm_epi_during_lmer %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**",
#     df = "**df**"
#   ) 
# lm_epi_during_summary

lm_raw_epi_during_zigamma_summary <- lm_raw_epi_during_zigamma_02 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Beta**"
  ) %>% 
  modify_column_indent(
    columns = label, 
    rows = variable %in% c("treatment", "time_since_end", "time_since_end:treatment"))

# filter out zero-inflated component
lm_raw_epi_during_zigamma_summary$table_body <- lm_raw_epi_during_zigamma_summary$table_body %>% 
  filter(component != "zi")
# change labels
lm_raw_epi_during_zigamma_summary$table_body$label <- c(
  time_since_end = "Time since end",
  treatmentremoval = "Treatment (removal)",
  `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
)

# final table 
lm_raw_epi_during_zigamma_summary
summary(lm_raw_epi_during_zigamma_02)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# predicted_epi_during <- ggpredict(lm_epi_during_lmer, terms = ~ time_since_end, type = "fixed")

predicted_raw_epi_during <- ggpredict(lm_raw_epi_during_zigamma_02, terms = c("time_since_end", "treatment"), type = "fixed")

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
# lm_epi_recovery_lmer <- lmer(
#   delta_continual_epi ~ time_since_end + (1|site), 
#   data = delta_epi_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass)

# lm_raw_epi_recovery_lmer <- lmer(
#   epi_biomass ~ time_since_end*treatment + (1|site), 
#   data = delta_epi_continual %>% filter(exp_dates == "after"), 
#   na.action = na.pass)

# lm_raw_epi_recovery_zigamma_01 <- glmmTMB(
#   epi_biomass ~ time_since_end*treatment + (1|site),
#   data = epi_continual_long %>% filter(exp_dates == "after"),
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1
# )

lm_raw_epi_recovery_zigamma_02 <- glmmTMB(
  epi_biomass ~ time_since_end*treatment + (1|site) + (1|year),
  data = epi_continual_long %>% filter(exp_dates == "after"),
  na.action = na.pass,
  family = ziGamma(link = "log"),
  ziformula = ~1
)


# diagnostics
# plot(simulateResiduals(lm_epi_recovery_lmer))
# check_model(lm_epi_recovery_lmer)

# plot(simulateResiduals(lm_raw_epi_recovery_zigamma_01))
plot(simulateResiduals(lm_raw_epi_recovery_zigamma_02))

# R2
MuMIn::r.squaredGLMM(lm_raw_epi_recovery_zigamma_02)

# summary table
# summary(lm_epi_recovery_lmer)
# lm_epi_recovery_summary <- lm_epi_recovery_lmer %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**",
#     df = "**df**"
#   ) 
# lm_epi_recovery_summary

lm_raw_epi_recovery_zigamma_summary <- lm_raw_epi_recovery_zigamma_02 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Beta**"
  ) %>% 
  modify_column_indent(
    columns = label, 
    rows = variable %in% c("treatment", "time_since_end", "time_since_end:treatment"))

# filter out zero-inflated component
lm_raw_epi_recovery_zigamma_summary$table_body <- lm_raw_epi_recovery_zigamma_summary$table_body %>% 
  filter(component != "zi")
# change labels
lm_raw_epi_recovery_zigamma_summary$table_body$label <- c(
  time_since_end = "Time since end",
  treatmentremoval = "Treatment (removal)",
  `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
)

# final table 
# lm_raw_epi_recovery_zigamma_summary
summary(lm_raw_epi_recovery_zigamma_02)

# ⟞ ⟞ ii. predictions -----------------------------------------------------

# predicted_epi_recovery <- ggpredict(lm_epi_recovery_lmer, terms = ~ time_since_end, type = "fixed")

predicted_raw_epi_recovery <- ggpredict(lm_raw_epi_recovery_zigamma_02, terms = c("time_since_end[0:6.75, by = 0.25]", "treatment"), type = "fixed")

# ⟞ c. figure ------------------------------------------------------------

# epi_time <- ggplot() +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_point(data = delta_epi_continual, 
#              aes(x = time_since_end, y = delta_continual_epi, fill = site, shape = site), 
#              size = 2, alpha = 0.9) +
#   scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
#   scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
#   # overall
#   # geom_line(data = predicted_epi_after, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
#   # geom_ribbon(data = predicted_epi_after, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   geom_line(data = predicted_epi_during, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_epi_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   geom_line(data = predicted_epi_recovery, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
#   geom_ribbon(data = predicted_epi_recovery, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
#   delta_timeseries_theme("epi") +
#   labs(x = "Time since end of removal (years)", 
#        y = "\U0394 biomass \n (treatment - control)",
#        subtitle = "(d)")
# epi_time

 # ⟞ ⟞ new model -----------------------------------------------------------

raw_epi_time <- ggplot() +
  # reference lines
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # raw data
  geom_point(data = epi_continual_long, 
             aes(x = time_since_end, y = epi_biomass, color = treatment, shape = treatment, size = treatment), alpha = 0.2) +
  
  # model predictions
  geom_line(data = predicted_raw_epi_recovery, aes(x = x, y = predicted, lty = group, color = group), linewidth = 1) +
  geom_ribbon(data = predicted_raw_epi_recovery, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
  geom_line(data = predicted_raw_epi_during, aes(x = x, y = predicted, lty = group, color = group), linewidth = 1) +
  geom_ribbon(data = predicted_raw_epi_during, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
  
  # colors and shapes
  scale_color_manual(values = c(reference = reference_col, removal = removal_col),
                     labels = c("Reference", "Removal")) +
  # scale_fill_manual(values = c(reference = "#6D5A1800", removal = "#CC754066"),
  #                    labels = c("Reference", "Removal")) +
  scale_linetype_manual(values = c(reference = 2, removal = 1),
                        labels = c("Reference", "Removal")) +
  scale_shape_manual(values = c(reference = 1, removal = 16),
                     labels = c("Reference", "Removal")) +
  scale_size_manual(values = c(reference = 1, removal = 1.3),
                    labels = c("Reference", "Removal")) +
  
  # theming
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
  coord_cartesian(ylim = c(5, 155)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1) +
  labs(x = "Time since end of removal (years)", 
       y = "Biomass (dry g/m\U00B2)", 
       title = "(c)")

raw_epi_time

raw_epi_removal <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # raw data points
  geom_point(data = epi_continual_long %>% filter(treatment == "removal"), aes(x = time_since_end, y = epi_biomass), shape = 1, size = 1, alpha = 0.4, color = removal_col) +
  
  # prediction lines
  geom_line(data = predicted_raw_epi_during %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
  geom_line(data = predicted_raw_epi_recovery %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
  
  # confidence intervals
  geom_ribbon(data = predicted_raw_epi_during %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  geom_ribbon(data = predicted_raw_epi_recovery %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  coord_cartesian(ylim = c(4, 155)) +
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
       y = "Biomass (dry g/m\U00B2)",
       title = "(d) Removal")
raw_epi_removal

raw_epi_reference <- ggplot() +
  # x at 0 and y at 0 lines
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # raw data points
  geom_point(data = epi_continual_long %>% filter(treatment == "reference"), aes(x = time_since_end, y = epi_biomass), shape = 1, size = 1, alpha = 0.4, color = reference_col) +
  
  # prediction lines
  geom_line(data = predicted_raw_epi_during %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
  geom_line(data = predicted_raw_epi_recovery %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
  
  # confidence intervals
  geom_ribbon(data = predicted_raw_epi_during %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  geom_ribbon(data = predicted_raw_epi_recovery %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
  
  theme_bw() + 
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(limits = c(-10, 155)) +
  coord_cartesian(ylim = c(4, 155)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10)) +
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1) +
  labs(x = "Time since end of reference (years)", 
       y = "Biomass (dry g/m\U00B2)",
       title = "(e) Reference")
raw_epi_reference

# data frame of predictions
delta_epi_predictions_during <- predicted_raw_epi_during %>% 
  as.data.frame() %>% 
  select(x, group, predicted) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta = removal - reference) %>% 
  mutate(exp_dates = "during")

delta_epi_predictions_after <- predicted_raw_epi_recovery %>% 
  as.data.frame() %>% 
  select(x, group, predicted) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta = removal - reference) %>% 
  mutate(exp_dates = "after")

overall_epi_predictions <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_epi_continual,
             aes(x = time_since_end, y = delta_continual_epi), shape = 2, size = 1, alpha = 0.4) +
  
  # overall
  geom_line(data = delta_epi_predictions_during, aes(x = x, y = delta), linewidth = 1) +
  geom_line(data = delta_epi_predictions_after, aes(x = x, y = delta), linewidth = 1) +
  
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  coord_cartesian(ylim = c(-40, 75)) +
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
       y = "\U0394 biomass \n (removal \U2212 reference, dry g/m\U00B2)",
       fill = "Site",
       shape = "Site",
       title = "(f) Removal \U2212 reference")

overall_epi_predictions

##########################################################################-
# 5. manuscript tables ----------------------------------------------------
##########################################################################-

# ⟞ a. model summary tables -----------------------------------------------

# individual group tables
# lm_algae_tables <- tbl_merge(tbls = list(lm_algae_during_summary, lm_algae_recovery_summary),
#                              tab_spanner = c("**Removal**", "**Recovery**")) 

# lm_epi_tables <- tbl_merge(tbls = list(lm_epi_during_summary, lm_epi_recovery_summary),
#                            tab_spanner = c("**Removal**", "**Recovery**")) 

# lm_endo_tables <- tbl_merge(tbls = list(lm_endo_during_summary, lm_endo_recovery_summary),
#                             tab_spanner = c("**Removal**", "**Recovery**")) 

# stack tables
# lm_summary_tables <- tbl_stack(
#   tbls = list(lm_kelp_tables, lm_algae_tables, lm_epi_tables, lm_endo_tables),
#   group_header = c("Kelp", "Understory macroalgae", "Epilithic invertebrates", "Endolithic invertebrates"),
#   quiet = TRUE) %>% 
#   as_flex_table() %>% 
#   font(fontname = "Times New Roman", part = "all")

# lm_summary_tables %>%
#   save_as_docx(path = here::here("tables", "ms-tables", paste("tbl-S1_", today(), ".docx", sep = "")))

# tbl_stack(
#   tbls = list(lm_kelp_tables, lm_algae_tables, lm_epi_tables, lm_endo_tables),
#   group_header = c("Kelp", "Understory macroalgae", "Epilithic invertebrates", "Endolithic invertebrates"),
#   quiet = TRUE) %>%
#   as_gt() %>%
#   tab_options(table.font.names = "Times New Roman") %>%
#   gtsave(here::here("tables", "ms-tables", paste("tbl-S1_", today(), ".png", sep = "")),
#          vwidth = 1500, vheight = 1000)

# individual group tables
lm_algae_zigamma_tables <- tbl_merge(tbls = list(lm_raw_algae_during_zigamma_summary, lm_raw_algae_recovery_zigamma_summary),
                             tab_spanner = c("**Kelp removal**", "**Recovery**")) 

lm_ep_zigamma_tables <- tbl_merge(tbls = list(lm_raw_epi_during_zigamma_summary, lm_raw_epi_recovery_zigamma_summary),
                           tab_spanner = c("**Kelp removal**", "**Recovery**")) 

# stack tables
lm_zigamma_summary_tables <- tbl_stack(
  tbls = list(lm_kelp_zigamma_tables, lm_algae_zigamma_tables, lm_ep_zigamma_tables),
  group_header = c("Giant kelp", "Understory macroalgae", "Sessile invertebrates"),
  quiet = TRUE) %>% 
  as_flex_table() %>% 
  font(fontname = "Times New Roman", part = "all")

# lm_zigamma_summary_tables %>%
#   save_as_docx(path = here::here("tables", "ms-tables", paste("tbl-S1_", today(), ".docx", sep = "")))

##########################################################################-
# 6. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. raw biomass through time -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_algae_raw,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S5_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_epi_raw,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S6_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_endo_raw,
#        height = 8, width = 16, units = "cm",
#        dpi = 300)

# ⟞ b. raw algae and epi model --------------------------------------------

fig2_v1 <-  (algae_title + epi_title) /
            (raw_algae_time + raw_epi_time) /
            (overall_algae_predictions + overall_epi_predictions) +
  plot_layout(heights = c(1, 10, 10))
fig2_v1

fig2_v2 <- (algae_title + epi_title) /
  (raw_algae_removal + raw_epi_removal) /
  (raw_algae_reference + raw_epi_reference) /
  (overall_algae_predictions + overall_epi_predictions) +
  plot_layout(heights = c(1, 10, 10, 10))
fig2_v2

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-2_new-model_", today(), ".jpg", sep = "")),
#        plot = raw_groups,
#        height = 18, width = 14, units = "cm",
#        dpi = 400)

# v1: legend position c(0.86, 0.88)
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-2_new-model_v1_", today(), ".jpg", sep = "")),
#        plot = fig2_v1,
#        height = 15, width = 18, units = "cm",
#        dpi = 400)

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-2_new-model_v2_", today(), ".jpg", sep = "")),
#        plot = fig2_v2,
#        height = 24, width = 18, units = "cm",
#        dpi = 400)



