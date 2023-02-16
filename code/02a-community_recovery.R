
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
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria"))

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
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria"))

# ⟞ c. endolithic invertebrates ------------------------------------------

# total biomass
endo_biomass <- biomass %>% 
  filter(taxon_family == "Pholadidae")

# delta biomass
delta_endo_biomass <- endo_biomass %>% 
  dplyr::select(site, year, month, treatment, date, dry_gm2) %>% 
  group_by(site, year, month, treatment, date) %>% 
  summarize(total_dry = sum(dry_gm2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = treatment, values_from = total_dry) %>% 
  mutate(delta_annual = annual - control,
         delta_continual = continual - control) 

# joining with kelp deltas and biomass
delta_endo_continual <- delta_endo_biomass %>% 
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
  rename("control_endo" = control,
         "continual_endo" = continual,
         "delta_continual_endo" = delta_continual) %>% 
  full_join(., delta_continual %>% dplyr::select(sample_ID, continual, control, delta_continual), by = "sample_ID") %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria"))

##########################################################################-
# 2. timeseries plots -----------------------------------------------------
##########################################################################-

# ⟞ a. raw biomass --------------------------------------------------------

# ⟞ ⟞ i. algae -----------------------------------------------------------

delta_continual_sites_algae_raw <- delta_algae_continual %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("A. ", site_full, sep = ""),
    site == "napl" ~ paste("B. ", site_full),
    site == "mohk" ~ paste("C. ", site_full),
    site == "carp" ~ paste("D. ", site_full)
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(x = time_since_end, y = control_algae, col = site), alpha = 0.5, linewidth = 2.5) +
  # control
  geom_point(aes(x = time_since_end, y = control_algae, shape = site), size = 1.5, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual_algae, col = site), linewidth = 2.5) +
  geom_point(aes(x = time_since_end, y = continual_algae, shape = site, col = site), size = 1.5, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  site_raw_biomass_theme() +
  labs(x = "Time since end of experiment (years)", 
       y = expression(Understory~algae~biomass~(dry~g/m^{"2"}))) +
  facet_wrap(~strip, scales = "free_y")
delta_continual_sites_algae_raw

# ⟞ ⟞ ii. epi inverts ----------------------------------------------------

delta_continual_sites_epi_raw <- delta_epi_continual %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("A. ", site_full, sep = ""),
    site == "napl" ~ paste("B. ", site_full),
    site == "mohk" ~ paste("C. ", site_full),
    site == "carp" ~ paste("D. ", site_full)
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  # control
  geom_line(aes(x = time_since_end, y = control_epi, col = site), alpha = 0.5, linewidth = 2.5) +
  geom_point(aes(x = time_since_end, y = control_epi, shape = site), size = 1.5, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual_epi, col = site), linewidth = 2.5) +
  geom_point(aes(x = time_since_end, y = continual_epi, shape = site, col = site), size = 1.5, fill = "#FFFFFF") +
  scale_shape_manual(values = shape_palette_site) +
  scale_color_manual(values = color_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  site_raw_biomass_theme() +
  labs(x = "Time since end of experiment (years)", 
       y = expression(Epilithic~invertebrate~biomass~(dry~g/m^{"2"}))) +
  facet_wrap(~strip, scales = "free_y")
delta_continual_sites_epi_raw

# ⟞ ⟞ iii. endo inverts --------------------------------------------------

delta_continual_sites_endo_raw <- delta_endo_continual %>% 
  mutate(strip = case_when(
    site == "aque" ~ paste("A. ", site_full, sep = ""),
    site == "napl" ~ paste("B. ", site_full),
    site == "mohk" ~ paste("C. ", site_full),
    site == "carp" ~ paste("D. ", site_full)
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(aes(x = time_since_end, y = control_endo, col = site), alpha = 0.5, linewidth = 2.5) +
  # control
  geom_point(aes(x = time_since_end, y = control_endo, shape = site), size = 1.5, alpha = 0.5, fill = "#FFFFFF") +
  # continual
  geom_line(aes(x = time_since_end, y = continual_endo, col = site), linewidth = 2.5) +
  geom_point(aes(x = time_since_end, y = continual_endo, shape = site, col = site), size = 1.5, fill = "#FFFFFF") +
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
       y = expression(Endolithic~invertebrate~biomass~(dry~g/m^{"2"}))) +
  facet_wrap(~strip, scales = "free_y")

delta_continual_sites_endo_raw

##########################################################################-
# 3. start-during-after comparisons ---------------------------------------
##########################################################################-

# ⟞ a. algae -------------------------------------------------------------

# anova with random effect
anova_algae_2yrs <- lmer(delta_continual_algae ~ comp_2yrs + (1|site), 
                    data = delta_algae_continual %>% drop_na(comp_2yrs))
anova_algae_3yrs <- lmer(delta_continual_algae ~ comp_3yrs + (1|site), 
                         data = delta_algae_continual %>% drop_na(comp_3yrs))

# diagnostics
plot(simulateResiduals(anova_algae_2yrs))
check_model(anova_algae_2yrs)

plot(simulateResiduals(anova_algae_3yrs))
check_model(anova_algae_3yrs)

# summary
summary(anova_algae_2yrs)
summary(anova_algae_3yrs)

# least squares comparison
difflsmeans(anova_algae_2yrs, test.effs = "Group", ddf = "Kenward-Roger")
difflsmeans(anova_algae_3yrs, test.effs = "Group", ddf = "Kenward-Roger")

# extract predicted values
anova_algae_2yrs_df <- ggpredict(anova_algae_2yrs, terms = "comp_2yrs", type = "fixed") %>% 
  mutate(x = case_when(
    x == "start" ~ "Start of removal",
    x == "during" ~ "End of removal",
    x == "after" ~ "Recovery period"
  )) %>% 
  mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))

# ⟞ b. epilithic invertebrates -------------------------------------------

# anova with random effect
anova_epi_2yrs <- lmer(delta_continual_epi ~ comp_2yrs + (1|site), 
                         data = delta_epi_continual %>% drop_na(comp_2yrs))
anova_epi_3yrs <- lmer(delta_continual_epi ~ comp_3yrs + (1|site), 
                         data = delta_epi_continual %>% drop_na(comp_3yrs))

# diagnostics
plot(simulateResiduals(anova_epi_2yrs))
check_model(anova_epi_2yrs)

plot(simulateResiduals(anova_epi_3yrs))
check_model(anova_epi_3yrs)

# summary
summary(anova_epi_2yrs)
summary(anova_epi_3yrs)

# least squares comparison
difflsmeans(anova_epi_2yrs, test.effs = "Group", ddf = "Kenward-Roger")
difflsmeans(anova_epi_3yrs, test.effs = "Group", ddf = "Kenward-Roger")

# anova with random effect
anova_epi_2yrs <- lmer(delta_continual_epi ~ comp_2yrs + (1|site), 
                       data = delta_epi_continual %>% drop_na(comp_2yrs))
anova_epi_3yrs <- lmer(delta_continual_epi ~ comp_3yrs + (1|site), 
                       data = delta_epi_continual %>% drop_na(comp_3yrs))

# diagnostics
plot(simulateResiduals(anova_epi_2yrs))
check_model(anova_epi_2yrs)

plot(simulateResiduals(anova_epi_3yrs))
check_model(anova_epi_3yrs)

# summary
summary(anova_epi_2yrs)
summary(anova_epi_3yrs)

# least squares comparison
difflsmeans(anova_epi_2yrs, test.effs = "Group", ddf = "Kenward-Roger")
difflsmeans(anova_epi_3yrs, test.effs = "Group", ddf = "Kenward-Roger")

# extract predicted values
anova_epi_2yrs_df <- ggpredict(anova_epi_2yrs, terms = "comp_2yrs", type = "fixed") %>% 
  mutate(x = case_when(
    x == "start" ~ "Start of removal",
    x == "during" ~ "End of removal",
    x == "after" ~ "Recovery period"
  )) %>% 
  mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))

# ⟞ c. endolithic invertebrates ------------------------------------------

# anova with random effect
anova_endo_2yrs <- lmer(delta_continual_endo ~ comp_2yrs + (1|site), 
                       data = delta_endo_continual %>% drop_na(comp_2yrs))
anova_endo_3yrs <- lmer(delta_continual_endo ~ comp_3yrs + (1|site), 
                       data = delta_endo_continual %>% drop_na(comp_3yrs))

# diagnostics
plot(simulateResiduals(anova_endo_2yrs))
check_model(anova_endo_2yrs)

plot(simulateResiduals(anova_endo_3yrs))
check_model(anova_endo_3yrs)

# summary
summary(anova_endo_2yrs)
summary(anova_endo_3yrs)

# extract predicted values
anova_endo_2yrs_df <- ggpredict(anova_endo_2yrs, terms = "comp_2yrs", type = "fixed") %>% 
  mutate(x = case_when(
    x == "start" ~ "Start of removal",
    x == "during" ~ "End of removal",
    x == "after" ~ "Recovery period"
  )) %>% 
  mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))

# ⟞ d. figures -----------------------------------------------------------

# ⟞ ⟞ i. algae -----------------------------------------------------------

sda_algae_biomass <- ggplot(anova_algae_2yrs_df) +
  # horizontal line at 0
  geom_hline(yintercept = 0, lty = 2) +
  # annotate("rect", xmin = 0, xmax = 4, ymin = 125.5, ymax = 140, fill = "#FFFFFF") +
  # points and error bars
  geom_point(aes(x = x, y = predicted), size = 2) +
  geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
  # path between time periods
  geom_line(aes(x = x, y = predicted, group = 2)) +
  # add in post-hoc comparisons
  annotate("text", x = 1, y = 122, label = "a", size = 3) +
  annotate("text", x = 2, y = 122, label = "b", size = 3) +
  annotate("text", x = 3, y = 122, label = "a", size = 3) +
  # annotate("text", x = 0.55, y = 135, label = "Understory algae", size = 10) +
  # aesthetics
  scale_x_discrete(labels = wrap_format(10)) +
  scale_y_continuous(limits = c(-30, 130), breaks = c(0, 50, 100, 150), expand = c(0, 0)) +
  sda_biomass_theme() +
  labs(x = "Time period", 
       y = expression(Biomass~(dry~g/m^{"2"})),
       title = "Understory algae", 
       subtitle = "A")
sda_algae_biomass

# ⟞ ⟞ ii. epi inverts ----------------------------------------------------

sda_epi_biomass <- ggplot(anova_epi_2yrs_df) +
  # horizontal line at 0
  geom_hline(yintercept = 0, lty = 2) +
  # points and error bars
  geom_point(aes(x = x, y = predicted), size = 2) +
  # annotate("rect", xmin = 0, xmax = 4, ymin = 20.1, ymax = 24, fill = "#FFFFFF") +
  geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
  # path between time periods
  geom_line(aes(x = x, y = predicted, group = 1)) +
  # add in post-hoc comparisons
  annotate("text", x = 1, y = 20, label = "a", size = 3) +
  annotate("text", x = 2, y = 20, label = "b", size = 3) +
  annotate("text", x = 3, y = 20, label = "b", size = 3) +
  # annotate("text", x = 0.68, y = 23, label = "Epilithic invertebrates", size = 10) +
  # aesthetics
  scale_x_discrete(labels = wrap_format(10)) +
  scale_y_continuous(limits = c(-13, 22), expand = c(0, 0)) +
  sda_biomass_theme() +
  labs(x = "Time period", 
       y = expression(Biomass~(dry~g/m^{"2"})),
       title = "Epilithic invertebrates",
       subtitle = "C")
sda_epi_biomass

# ⟞ ⟞ iii. endo inverts --------------------------------------------------

sda_endo_biomass <- ggplot(anova_endo_2yrs_df) +
  # horizontal line at 0
  geom_hline(yintercept = 0, lty = 2) +
  # annotate("rect", xmin = 0, xmax = 4, ymin = 551, ymax = 650, fill = "#FFFFFF") +
  # points and error bars
  geom_point(aes(x = x, y = predicted), size = 2) +
  geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
  # path between time periods
  geom_line(aes(x = x, y = predicted, group = 1)) +
  # annotate("text", x = 0.75, y = 625, label = "Endolithic invertebrates", size = 10) +
  # aesthetics
  scale_x_discrete(labels = wrap_format(10)) +
  scale_y_continuous(limits = c(-100, 550), expand = c(0, 0), breaks = c(0, 200, 400, 600)) +
  sda_biomass_theme() +
  labs(x = "Time period", 
       y = expression(Biomass~(dry~g/m^{"2"})),
       title = "Endolithic invertebrates",
       subtitle = "E")
sda_endo_biomass

##########################################################################-
# 2. algae linear model ---------------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_algae_during_lmer <- lmer(
  delta_continual_algae ~ time_since_end + (1|site), 
  data = delta_algae_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass)

# diagnostics
plot(simulateResiduals(lm_algae_during_lmer))
check_model(lm_algae_during_lmer)

# Rsquared
r.squaredGLMM(lm_algae_during_lmer)

# summary
summary(lm_algae_during_lmer)
lm_algae_during_summary <- lm_algae_during_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_algae_during_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_algae_during <- ggpredict(lm_algae_during_lmer, terms = ~ time_since_end, type = "fixed")

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_algae_recovery_lmer <- lmer(
  delta_continual_algae ~ time_since_end + (1|site), 
  data = delta_algae_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass)

# diagnostics
plot(simulateResiduals(lm_algae_recovery_lmer))
check_model(lm_algae_recovery_lmer)

# Rsquared
r.squaredGLMM(lm_algae_recovery_lmer)

# summary
summary(lm_algae_recovery_lmer)
lm_algae_recovery_summary <- lm_algae_recovery_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_algae_recovery_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_algae_recovery <- ggpredict(lm_algae_recovery_lmer, terms = ~time_since_end, type = "fixed")

# algae decreases to control in 4.7 years on average
ggpredict(lm_algae_recovery_lmer, terms = "time_since_end [4.6:4.7 by = 0.001]", type = "fixed")

# lower bound of 95% conf int: 2.6 years
ggpredict(lm_algae_recovery_lmer, terms = "time_since_end [2.6:2.7 by = 0.001]", type = "fixed")

# upper bound of 95% conf int: 7.7 years
ggpredict(lm_algae_recovery_lmer, terms = "time_since_end [7.6:7.7 by = 0.001]", type = "fixed")

# ⟞ c. figure ------------------------------------------------------------

algae_time <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_algae_continual, 
             aes(x = time_since_end, y = delta_continual_algae, fill = site, shape = site), size = 2, alpha = 0.9) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  # new_scale("color") + 
  # overall
  geom_line(data = predicted_algae_recovery, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_algae_recovery, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  geom_line(data = predicted_algae_during, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_algae_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
  delta_timeseries_theme("algae") +
  labs(x = "Time since end of removal (years)",
       y = "\U0394 biomass \n (treatment - control)",
       title = " ", subtitle = "B",
       fill = "Site", shape = "Site")

algae_time

##########################################################################-
# 3. epi. invert linear model ---------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_epi_during_lmer <- lmer(
  delta_continual_epi ~ time_since_end + (1|site), 
  data = delta_epi_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass)

# diagnostics
plot(simulateResiduals(lm_epi_during_lmer))
check_model(lm_epi_during_lmer)

# R2
r.squaredGLMM(lm_epi_during_lmer)

# summary
summary(lm_epi_during_lmer) 
lm_epi_during_summary <- lm_epi_during_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_epi_during_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_epi_during <- ggpredict(lm_epi_during_lmer, terms = ~ time_since_end, type = "fixed")

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_epi_recovery_lmer <- lmer(
  delta_continual_epi ~ time_since_end + (1|site), 
  data = delta_epi_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass)

# diagnostics
plot(simulateResiduals(lm_epi_recovery_lmer))
check_model(lm_epi_recovery_lmer)

# R2
MuMIn::r.squaredGLMM(lm_epi_recovery_lmer)

# summary table
summary(lm_epi_recovery_lmer)
lm_epi_recovery_summary <- lm_epi_recovery_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_epi_recovery_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_epi_recovery <- ggpredict(lm_epi_recovery_lmer, terms = ~ time_since_end, type = "fixed")

# ⟞ c. figure ------------------------------------------------------------

epi_time <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_epi_continual, 
             aes(x = time_since_end, y = delta_continual_epi, fill = site, shape = site), 
             size = 2, alpha = 0.9) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  # overall
  # geom_line(data = predicted_epi_after, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  # geom_ribbon(data = predicted_epi_after, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  geom_line(data = predicted_epi_during, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
  geom_ribbon(data = predicted_epi_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
  delta_timeseries_theme("epi") +
  labs(x = "Time since end of removal (years)", 
       y = "\U0394 biomass \n (treatment - control)",
       title = " ", subtitle = "D")
epi_time

##########################################################################-
# 4. endo. invert linear model --------------------------------------------
##########################################################################-

# ⟞ a. during removal -----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_endo_during_lmer <- lmer(
  delta_continual_endo ~ time_since_end + (1|site), 
  data = delta_endo_continual %>% filter(exp_dates == "during"), 
  na.action = na.pass
  )

# check
plot(simulateResiduals(lm_endo_during_lmer))
check_model(lm_endo_during_lmer)

# R2
r.squaredGLMM(lm_endo_during_lmer)

# summmary
summary(lm_endo_during_lmer)
lm_endo_during_summary <- lm_endo_during_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_endo_during_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_endo_during <- ggpredict(lm_endo_during_lmer, terms = ~ time_since_end, type = "fixed")

# ⟞ b. recovery period ----------------------------------------------------

# ⟞ ⟞ i. model and diagnostics  -------------------------------------------

# model
lm_endo_recovery_lmer <- lmer(
  delta_continual_endo ~ time_since_end + (1|site), 
  data = delta_endo_continual %>% filter(exp_dates == "after"), 
  na.action = na.pass
  )

# check
plot(simulateResiduals(lm_endo_recovery_lmer))
check_model(lm_endo_recovery_lmer)

# R2
r.squaredGLMM(lm_endo_recovery_lmer)

# summary
summary(lm_endo_recovery_lmer)
lm_endo_recovery_summary <- lm_endo_recovery_lmer %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)
lm_endo_recovery_summary

# ⟞ ⟞ ii. predictions -----------------------------------------------------

predicted_endo_recovery <- ggpredict(lm_endo_recovery_lmer, terms = ~ time_since_end, type = "fixed")

# ⟞ c. figure ------------------------------------------------------------

endo_time <- ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(data = delta_endo_continual, 
             aes(x = time_since_end, y = delta_continual_endo, fill = site, shape = site), size = 2, alpha = 0.9) +
  scale_shape_manual(values = shape_palette_site) +
  scale_fill_manual(values = color_palette_site) +
  # new_scale("color") + 
  # overall
  # geom_line(data = predicted_clam_after, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  # geom_ribbon(data = predicted_clam_after, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  # geom_line(data = predicted_clam_during, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
  # geom_ribbon(data = predicted_clam_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
  delta_timeseries_theme("endo") +
  labs(x = "Time since end of removal (years)", 
       y = "\U0394 biomass \n (treatment - control)",
       title = " ", subtitle = "F")
endo_time


##########################################################################-
# 5. manuscript tables ----------------------------------------------------
##########################################################################-

# individual group tables
lm_algae_tables <- tbl_merge(tbls = list(lm_algae_during_summary, lm_algae_recovery_summary),
                             tab_spanner = c("**Removal**", "**Recovery**")) 

lm_epi_tables <- tbl_merge(tbls = list(lm_epi_during_summary, lm_epi_recovery_summary),
                           tab_spanner = c("**Removal**", "**Recovery**")) 

lm_endo_tables <- tbl_merge(tbls = list(lm_endo_during_summary, lm_endo_recovery_summary),
                            tab_spanner = c("**Removal**", "**Recovery**")) 

# stack tables
lm_summary_tables <- tbl_stack(
  tbls = list(lm_kelp_tables, lm_algae_tables, lm_epi_tables, lm_endo_tables),
  group_header = c("Kelp", "Algae", "Epilithic invertebrates", "Endolithic invertebrates"),
  quiet = TRUE) %>% 
  as_gt()
lm_summary_tables

# gtsave(lm_summary_tables,
#        here::here("tables", "ms-tables", paste("lm_summary_tables_", today(), ".png", sep = "")),
#        vwidth = 1500, vheight = 1000)

##########################################################################-
# 6. manuscript figures ---------------------------------------------------
##########################################################################-

# ⟞ a. s-d-a and model predictions ----------------------------------------

sda_time_algae <- (sda_algae_biomass + algae_time) +
  plot_layout(widths = c(1.5, 2.5)) 

sda_time_epi <- (sda_epi_biomass + epi_time) +
  plot_layout(widths = c(1.5, 2.5))

sda_time_endo <- (sda_endo_biomass + endo_time) +
  plot_layout(widths = c(1.5, 2.5)) 

sda_time_together_v2 <- sda_time_algae/
                        sda_time_epi/
                        sda_time_endo

ggsave(here::here("figures", "ms-figures",
                  paste("fig-3_", today(), ".jpg", sep = "")),
       plot = sda_time_together_v2,
       height = 18, width = 16, units = "cm", 
       dpi = 300)

sda_time_together <- (sda_algae_biomass + algae_time) /
  (sda_epi_biomass + epi_time) /
  (sda_endo_biomass + endo_time) +
  plot_layout(widths = c(2, 3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30))

sda_time_together

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-3_", today(), ".jpg", sep = "")),
#        plot = sda_time_together,
#        height = 17, width = 16, dpi = 150)

# ⟞ b. raw biomass through time -------------------------------------------

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S4_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_algae_raw,
#        height = 8, width = 16, dpi = 150)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S5_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_epi_raw,
#        height = 8, width = 16, dpi = 150)
# 
# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-S6_", today(), ".jpg", sep = "")),
#        plot = delta_continual_sites_endo_raw,
#        height = 8, width = 16, dpi = 150)








