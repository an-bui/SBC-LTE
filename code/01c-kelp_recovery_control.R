##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "01a-kelp_recovery.R"))

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
    site == "aque" ~ 3.7,
    site == "napl" ~ 3.4,
    site == "mohk" ~ 5.4,
    site == "carp" ~ 4
  ), time_pred_low = case_when(
    site == "aque" ~ -3.5,
    site == "napl" ~ -3.8,  
    site == "mohk" ~ -1.8,
    site == "carp" ~ -3
  ), time_pred_high = case_when(
    site == "aque" ~ 11.9,  
    site == "napl" ~ 11.5,
    site == "mohk" ~ 14,
    site == "carp" ~ 12.4
  )) 

ggplot(mean_kelp_all_sites, aes(x = mean_control, y = time_to_recovery)) +
  geom_point() +
  geom_errorbar(aes(ymin = time_pred_low, ymax = time_pred_high))

lm_kelp_recovery_raw <- lmerTest::lmer(
  control ~ time_since_end + (1|site),
  data = delta_continual %>% filter(exp_dates == "after"),
  na.action = na.pass
)

plot(ggpredict(lm_kelp_recovery_raw, terms = ~ time_since_end, type = "fixed")) +
  geom_point(data = delta_continual %>% filter(exp_dates == "after"), aes(x = time_since_end, y = control, color = site)) +
  scale_color_manual(values = color_palette_site)
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
  "carp",       4,          -4.1,         11.4,         -611.64,       662.11
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

# arrows
delta_continual %>% 
  mutate(site = fct_relevel(site, "napl", "aque", "carp", "mohk")) %>% 
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
  ) +
  theme(legend.position = "none") +
  facet_wrap(~site, ncol = 2)

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
       fill = "Site", shape = "Site") +
  theme_bw() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

delta_vs_biomass