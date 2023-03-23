
##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "01a-kelp_recovery.R"))

##########################################################################-
# 1. predictions ----------------------------------------------------------
##########################################################################-

# ⟞ a. Arroyo Quemado ----------------------------------------------------

# predicted line crosses 0 at 4.0
plot(predicted_kelp_after_aque)
# upper bound crosses 0 at -3.5 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-3.5:-3.4 by = 0.001]", type = "random", condition = c(site = "aque"))
# lower bound crosses 0 at 12.0 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [12.0:12.1 by = 0.001]", type = "random", condition = c(site = "aque"))

# ⟞ b. Naples ------------------------------------------------------------

# predicted line crosses 0 at 3.4
plot(predicted_kelp_after_napl)
# upper bound crosses 0 at -4.0 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-4.0:-3.9 by = 0.001]", type = "random", condition = c(site = "napl"))
# lower bound crosses 0 at 11.5 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [11.5:11.6 by = 0.001]", type = "random", condition = c(site = "napl"))
# confidence interval at 3.4
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.3:3.5 by = 0.001]", type = "fixed")

# ⟞ c. Mohawk ------------------------------------------------------------

# predicted line crosses 0 at 5.4
plot(predicted_kelp_after_mohk)
# upper bound crosses 0 at -1.6 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-1.6:-1.5 by = 0.01]", type = "random", condition = c(site = "mohk"))
# lower bound crosses 0 at 14 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [13.9:14 by = 0.01]", type = "random", condition = c(site = "mohk"))

# ⟞ d. Carpinteria -------------------------------------------------------

# predicted line crosses 0 at 3.3
plot(predicted_kelp_after_carp)
# upper bound crosses 0 at -4.1 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [-4.1:-4.0 by = 0.01]", type = "random", condition = c(site = "carp"))
# lower bound crosses 0 at 11.4 years
ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [11.3:11.4 by = 0.01]", type = "random", condition = c(site = "carp"))

##########################################################################-
# 2. model without CARP ---------------------------------------------------
##########################################################################-

lm_kelp_recovery_nocarp_lmer <- lmerTest::lmer(
  delta_continual ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "after" & site != "carp"), 
  na.action = na.pass)

plot(DHARMa::simulateResiduals(lm_kelp_recovery_nocarp_lmer))
performance::check_model(lm_kelp_recovery_nocarp_lmer)

check_normality(lm_kelp_recovery_nocarp_lmer)
check_heteroscedasticity(lm_kelp_recovery_nocarp_lmer)

summary(lm_kelp_recovery_nocarp_lmer)

predicted_kelp_recovery_nocarp_overall <- ggpredict(lm_kelp_recovery_nocarp_lmer, terms = ~ time_since_end, type = "fixed") 
predicted_kelp_recovery_nocarp_aque <- ggpredict(lm_kelp_recovery_nocarp_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque")) 
predicted_kelp_recovery_nocarp_napl <- ggpredict(lm_kelp_recovery_nocarp_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl")) 
predicted_kelp_recovery_nocarp_mohk <- ggpredict(lm_kelp_recovery_nocarp_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk")) 

kelp_recovery_nocarp_overall_plot <- delta_continual %>% 
  filter(exp_dates == "after" & site != "carp") %>% 
  ggplot(aes(x = time_since_end, y = delta_continual)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(aes(fill = site, shape = site)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  geom_line(data = predicted_kelp_recovery_nocarp_overall, aes(x = x, y = predicted)) +
  geom_ribbon(data = predicted_kelp_recovery_nocarp_overall, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  # coord_cartesian(ylim = c(-25, 1600), xlim = c(-0.2, 5.7)) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
  # scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Time since end of removal", 
       y = expression(Giant~kelp~biomass~"(dry"~g/m^{"2"}~")"),
       fill = "Site", shape = "Site") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.position = c(0.8, 0.12),
        legend.key.size = unit(0.3, units = "cm")) 
kelp_recovery_nocarp_overall_plot

##########################################################################-
# 3. model of control biomass in post-removal period ----------------------
##########################################################################-

lm_kelp_control_during1 <- lmerTest::lmer(
  control ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "during"),
  na.action = na.pass
)

lm_kelp_control_during2 <- lmerTest::lmer(
  control ~ time_since_end*quality + (1|site),
  data = delta_continual %>% filter(exp_dates == "during"),
  na.action = na.pass
)

lm_kelp_control_during2 %>% 
  tbl_regression()

AICc(lm_kelp_control_during1, lm_kelp_control_during2, lm_kelp_control_during3) %>% 
  arrange(AICc)

plot(DHARMa::simulateResiduals(lm_kelp_control_during2))
performance::check_model(lm_kelp_control_during2)

summary(lm_kelp_control_during1)
# significant effect of time

summary(lm_kelp_control_during2)

r.squaredGLMM(lm_kelp_control_during)

predicted_control_kelp_during_overall <- ggpredict(lm_kelp_control_during, terms = ~ time_since_end, type = "fixed")

lm_kelp_control_after1 <- lmerTest::lmer(
  control ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"),
  na.action = na.pass
)

lm_kelp_control_after2 <- lmerTest::lmer(
  control ~ time_since_end*quality + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"),
  na.action = na.pass
)

plot(DHARMa::simulateResiduals(lm_kelp_control_after2))
performance::check_model(lm_kelp_control_after2)

summary(lm_kelp_control_after2)
# significant effect of time

lm_kelp_control_after2 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)

r.squaredGLMM(lm_kelp_control_after)

AICc(lm_kelp_control_after1, lm_kelp_control_after2) %>% 
  arrange(AICc) %>% 
  as.data.frame() %>% 
  rownames_to_column("model") %>% 
  gt()

lm_kelp_continual_after1 <- lmerTest::lmer(
  continual ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"),
  na.action = na.pass
)

lm_kelp_continual_after2 <- lmerTest::lmer(
  continual ~ time_since_end*quality + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"),
  na.action = na.pass
)




plot(DHARMa::simulateResiduals(lm_kelp_continual_after2))
performance::check_model(lm_kelp_continual_after2)

AICc(lm_kelp_continual_after1, lm_kelp_continual_after2) %>% 
  arrange(AICc)

lm_kelp_continual_after2 %>% 
  tbl_regression() %>% 
  bold_p(t = 0.05)

summary(lm_kelp_continual_after2)
# no significant effect of time

r.squaredGLMM(lm_kelp_continual_after)

predicted_control_kelp_after_overall <- ggpredict(lm_kelp_control_after, terms = ~ time_since_end, type = "fixed") 
predicted_continual_kelp_after_overall <- ggpredict(lm_kelp_continual_after, terms = ~ time_since_end, type = "fixed") 

ggpredict(lm_kelp_continual_after, terms = ~ time_since_end, type = "random", condition = c(site = "carp")) %>% 
  plot()

continual_kelp_after_overall_plot <- delta_continual %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = time_since_end, y = continual)) +
  geom_point(aes(fill = site, shape = site)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  # geom_line(data = predicted_continual_kelp_after_overall, aes(x = x, y = predicted)) +
  # geom_ribbon(data = predicted_continual_kelp_after_overall, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  coord_cartesian(ylim = c(-25, 1600), xlim = c(-0.2, 5.7)) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Time since end of removal", 
       y = expression(Giant~kelp~biomass~"(dry"~g/m^{"2"}~")"),
       fill = "Site", shape = "Site") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8),
        plot.title.position = "plot",
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.position = c(0.8, 0.82),
        legend.key.size = unit(0.3, units = "cm")) 

control_kelp_overall_plot <- delta_continual %>% 
  ggplot(aes(x = time_since_end, y = control)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(fill = site, shape = site)) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  geom_line(data = predicted_control_kelp_during_overall, aes(x = x, y = predicted)) +
  geom_ribbon(data = predicted_control_kelp_during_overall, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = predicted_control_kelp_after_overall, aes(x = x, y = predicted)) +
  geom_ribbon(data = predicted_control_kelp_after_overall, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  coord_cartesian(ylim = c(-25, 1600), xlim = c(-7.8, 5.7)) +
  scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Time since end of removal", 
       y = expression(Giant~kelp~biomass~"(dry"~g/m^{"2"}~")"),
       fill = "Site", shape = "Site") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8),
        plot.title.position = "plot",
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.position = c(0.91, 0.87),
        legend.margin = margin(0.3),
        legend.key.size = unit(0.3, units = "cm")) 

control_kelp_overall_plot

continual_kelp_after_overall_plot

##########################################################################-
# 4. recovery time vs biomass ---------------------------------------------
##########################################################################-

# ⟞ a. data frame ---------------------------------------------------------

# Not sure how to do this reproducibly. But pulled time to recovery from above predicted values with prediction intervals. Then calculated mean kelp biomass in the kelp year of recovery from `delta_continual`. For MOHK, calculated biomass using the most recent kelp year (2021-2022).

mean_kelp_all_sites <- delta_continual %>% 
  filter(exp_dates == "after") %>% 
  group_by(site) %>% 
  summarize(mean_control = mean(control),
            sd_control = sd(control),
            se_control = se(control)) %>% 
  ungroup() %>% 
  mutate(time_to_recovery = case_when(
    site == "aque" ~ 3.8,
    site == "napl" ~ 3.4,
    site == "mohk" ~ 5.4,
    site == "carp" ~ 3.3
  ), time_pred_low = case_when(
    site == "aque" ~ -3.5,
    site == "napl" ~ -4,  
    site == "mohk" ~ -1.6,
    site == "carp" ~ -4.1
  ), time_pred_high = case_when(
    site == "aque" ~ 12.0,  
    site == "napl" ~ 11.5,
    site == "mohk" ~ 14,
    site == "carp" ~ 11.4
  ), time_pred_full = paste(time_pred_low, ", ", time_pred_high, sep = ""),
  mean_se = paste(round(mean_control, 1), "\U00B1", round(se_control, 1), sep = " ")) %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  mutate(value = fct_relevel(value, "Carpinteria", "Naples", "Arroyo Quemado", "Mohawk"))

mean_kelp_all_sites %>% 
  select(value, mean_se, time_to_recovery, time_pred_full) %>% 
  arrange(value) %>% 
  gt() %>% 
  cols_label(
    value = "Site",
    mean_se = "Mean \U00B1 SE",
    time_to_recovery = "Time to recovery",
    time_pred_full = "95% PI"
  )

ggplot(mean_kelp_all_sites, aes(x = mean_control, y = time_to_recovery)) +
  geom_point() +
  geom_errorbar(aes(ymin = time_pred_low, ymax = time_pred_high))


# confidence intervals of biomass in kelp control plot from specific site recovery times
ci_3.8 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [3.7:3.9 by = 0.01]", type = "fixed") %>% 
  filter(x == "3.8") %>% 
  as_tibble() %>% 
  mutate(site = "aque")
ci_3.42 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [3.41:3.43 by = 0.01]", type = "fixed") %>% 
  filter(x == "3.42") %>% 
  as_tibble() %>% 
  mutate(site = "napl")
ci_5.4 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [5.3:5.4 by = 0.01]", type = "fixed") %>% 
  filter(x == "5.4") %>% 
  as_tibble() %>% 
  mutate(site = "mohk")
ci_3.33 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [3.3:3.4 by = 0.01]", type = "fixed") %>% 
  filter(x == "3.33") %>% 
  as_tibble() %>% 
  mutate(site = "carp")
ci_all <- rbind(ci_3.8, ci_3.42, ci_5.4, ci_3.33) %>% 
  select(site, conf.low, conf.high) 

# prediction intervals for biomass in kelp control plot
pi_3.8 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [3.7:3.9 by = 0.01]", type = "random", condition = c(site = "aque")) %>% 
  filter(x == "3.8") %>% 
  as_tibble() %>% 
  mutate(site = "aque") 
pi_3.42 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [3.41:3.43 by = 0.01]", type = "random", condition = c(site = "napl")) %>% 
  filter(x == "3.42") %>% 
  as_tibble() %>% 
  mutate(site = "napl") 
pi_5.4 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [5.3:5.4 by = 0.01]", type = "random", condition = c(site = "mohk")) %>% 
  filter(x == "5.4") %>% 
  as_tibble() %>% 
  mutate(site = "mohk") 
pi_3.33 <- ggpredict(lm_kelp_recovery_raw, terms = "time_since_end [3.3:3.4 by = 0.01]", type = "random", condition = c(site = "carp")) %>% 
  filter(x == "3.33") %>% 
  as_tibble() %>% 
  mutate(site = "carp") 

rec_time_2 <- rbind(pi_3.8, pi_3.42, pi_5.4, pi_3.33) %>% 
  rename(pred.low = conf.low,
         pred.high = conf.high) %>% 
  # select(site, pred.low, pred.high) %>% 
  left_join(., ci_all, by = "site") %>% 
  left_join(., rec_time, by = "site")


rec_time <- tribble(
  ~site, ~time_to_recovery, ~pi_time_low, ~pi_time_high, ~pi_kelp_low, ~pi_kelp_high, 
  "aque",       3.7,          -3.5,         12,           -542.53,       734.94,
  "napl",       3.4,          -4,           11.5,         -373.56,       900.80,
  "mohk",       5.4,           -1.6,         14.4,         -309.19,       991.41,
  "carp",       4,             -4.1,         11.4,         -611.64,       662.11
) %>% 
  left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
  rename("site_full" = value) %>% 
  mutate(site_full = fct_relevel(site_full, "Arroyo Quemado", "Naples", "Mohawk", "Carpinteria")) %>% 
  left_join(., mean_kelp_all_sites, by = "site") 

# ⟞ b. figure -------------------------------------------------------------

rec_time_plot <- ggplot(rec_time, aes(x = mean_control, y = time_to_recovery, shape = site, fill = site)) +
  geom_errorbar(aes(xmin = mean_control - se_control, xmax = mean_control + se_control), width = 1) +
  geom_errorbar(aes(ymin = pi_time_low, ymax = pi_time_high)) +
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

rec_time_plot_2 <- ggplot(rec_time_2, aes(x = predicted, y = time_to_recovery)) +
  # geom_point() +
  # prediction and confidence intervals for x axis: predicted mean kelp biomass
  geom_errorbar(aes(xmin = pred.low, xmax = pred.high)) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), color = "red") +
  geom_errorbar(aes(ymin = pi_time_low, ymax = pi_time_high))
rec_time_plot_2


##########################################################################-
# 5. control vs removal kelp biomass plot ---------------------------------
##########################################################################-

aque_arrows <- arrow_plot_fxn("aque") +
  labs(title = "(c) Arroyo Quemado") 
napl_arrows <- arrow_plot_fxn("napl") +
  labs(title = "(b) Naples") 
mohk_arrows <- arrow_plot_fxn("mohk") +
  labs(title = "(d) Mohawk") 
carp_arrows <- arrow_plot_fxn("carp") +
  labs(title = "(a) Carpinteria") 


delta_vs_biomass <- delta_continual %>% 
  filter(exp_dates == "after") %>% 
  ggplot(aes(x = control, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(fill = site, shape = site), size = 2) +
  scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
  # geom_smooth(method = "lm", color = "#000000", se = FALSE, linewidth = 1) +
  labs(x = expression(Control~biomass~"(dry"~g/m^{"2"}~")", sep = ""), 
       y = "\U0394 biomass (treatment - control)",
       fill = "Site", shape = "Site",
       title = "(a)") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.title = element_text(size = 8),
        plot.title.position = "plot",
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

delta_vs_biomass

##########################################################################-
# 6. manuscript tables ----------------------------------------------------
##########################################################################-

tor_table <- mean_kelp_all_sites %>% 
  select(value, mean_se, time_to_recovery, time_pred_full) %>% 
  arrange(value) %>% 
  gt() %>% 
  cols_label(
    value = "Site",
    mean_se = "Mean \U00B1 SE",
    time_to_recovery = "Time to recovery",
    time_pred_full = "95% PI"
  ) %>% 
  tab_options(table.font.names = "Times New Roman") 

# gtsave(tor_table,
#        here::here("tables", "ms-tables", paste("tbl-appendix2-S1_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)

##########################################################################-
# 6. manuscript figures ---------------------------------------------------
##########################################################################-

kelp_recovery_nocarp_overall_plot

ggsave(here::here("figures", "ms-figures",
                  paste("fig-appendix2-S2_", today(), ".jpg", sep = "")),
       plot = kelp_recovery_nocarp_overall_plot,
       height = 8, width = 8, units = "cm",
       dpi = 400)

panel_fig <- delta_vs_biomass + control_kelp_after_overall_plot

ggsave(here::here("figures", "ms-figures",
                  paste("fig-appendix2-S4_", today(), ".jpg", sep = "")),
       plot = panel_fig,
       height = 8, width = 16, units = "cm",
       dpi = 400)

control_label <- ggplot(data.frame(l = "Control biomass", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 3) + 
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(clip = "off")

continual_label <- ggplot(data.frame(l = "Treatment biomass", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = 3) + 
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(clip = "off")

top_arrows <- plot_grid(carp_arrows, napl_arrows, ncol = 2)
bottom_arrows <- plot_grid(aque_arrows, mohk_arrows, ncol = 2)
grid_arrows <- plot_grid(top_arrows, bottom_arrows, nrow = 2) %>% 
  plot_grid(., control_label, nrow = 2, rel_heights = c(12, 1)) %>% 
  plot_grid(continual_label, ., ncol = 2, rel_widths = c(1, 12))
grid_arrows

ggsave(here::here("figures", "ms-figures",
                  paste("fig-appendix2-S5_", today(), ".jpg", sep = "")),
       plot = grid_arrows,
       height = 8, width = 10, units = "cm",
       dpi = 400)





ggsave(here::here("figures", "ms-figures",
                  paste("fig-appendix2-S6_", today(), ".jpg", sep = "")),
       plot = control_kelp_overall_plot,
       height = 8, width = 14, units = "cm",
       dpi = 400)

ggsave(here::here("figures", "ms-figures",
                  paste("fig-appendix2-S7_", today(), ".jpg", sep = "")),
       plot = continual_kelp_after_overall_plot,
       height = 8, width = 8, units = "cm",
       dpi = 400)
