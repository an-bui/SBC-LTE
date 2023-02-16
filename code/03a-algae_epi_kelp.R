##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "02a-community_recovery.R"))

##########################################################################-
# 1. data frames and wranging functions -----------------------------------
##########################################################################-

cor.test(delta_algae_continual$delta_continual_algae, delta_algae_continual$delta_continual,
         method = "spearman")

lm_algae_kelp_during_m1 <- lm(
  delta_continual_algae ~ delta_continual, 
  data = delta_algae_continual %>% filter(exp_dates == "during"),
  na.action = na.omit
)
lm_algae_kelp_during_m2 <- lmer(
  delta_continual_algae ~ delta_continual + (1|year) + (1|site), 
  data = delta_algae_continual %>% filter(exp_dates == "during")
)
lm_algae_kelp_during_m3 <- lmer(
  delta_continual_algae ~ delta_continual + (1|site), 
  data = delta_algae_continual %>% filter(exp_dates == "during")
)

# diagnostics (NOT GOOD)
plot(simulateResiduals(lm_algae_kelp_during_m1))
plot(simulateResiduals(lm_algae_kelp_during_m2))
check_model(lm_algae_kelp_during_m1)
check_model(lm_algae_kelp_during_m2)
check_model(lm_algae_kelp_during_m3)

# model checks
check_autocorrelation(lm_algae_kelp_during_m1) # autocorrelation of residuals
check_autocorrelation(lm_algae_kelp_during_m2) # no autocorrelation
check_heteroscedasticity(lm_algae_kelp_during_m1) # herteroscedastic
check_normality(lm_algae_kelp_during_m3)

# Rsquared
r.squaredGLMM(lm_algae_kelp_during_m1)
r.squaredGLMM(lm_algae_kelp_during_m2)

# summary
summary(lm_algae_kelp_during_m1)
summary(lm_algae_kelp_during_m2)

# model selection
AICc(lm_algae_kelp_during_m1, lm_algae_kelp_during_m2)


predicted_delta_algae_vs_kelp <- ggpredict(delta_algae_vs_kelp_lmer, terms = ~ delta_continual, type = "fixed")


algae_vs_kelp <- delta_algae_continual %>% 
  mutate(exp_dates = case_when(
    exp_dates == "during" ~ "During removal",
    exp_dates == "after" ~ "Post-removal"
  )) %>% 
  # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
  ggplot(aes(x = delta_continual, y = delta_continual_algae)) +
  geom_hline(aes(yintercept = 0), lty = 3, color = "grey") +
  geom_vline(aes(xintercept = 0), lty = 3, color = "grey") +
  geom_point(aes(shape = exp_dates, fill = exp_dates), size = 5, shape = 21) + 
  geom_line(data = predicted_delta_algae_vs_kelp, aes(x = x, y = predicted), size = 2) +
  scale_linetype_manual(values = c("During removal" = 2, "Post-removal" = 1)) +
  scale_fill_manual(values = c("During removal" = "#FFFFFF", "Post-removal" = under_col)) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
  labs(x = "\U0394 kelp biomass (treatment - control)",
       y = "\U0394 understory algae biomass (treatment - control)") +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

algae_vs_kelp_spearman <- delta_algae_continual %>% 
  mutate(exp_dates = case_when(
    exp_dates == "during" ~ "During removal",
    exp_dates == "after" ~ "Recovery period"
  )) %>% 
  ggplot(aes(x = delta_continual, y = delta_continual_algae)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  geom_point(aes(shape = exp_dates, fill = exp_dates), size = 5, shape = 21) +
  stat_cor(method = "spearman", cor.coef.name = "rho") +
  geom_smooth(aes(linetype = exp_dates), method = "lm", color = "black", linewidth = 3, se = FALSE) +
  scale_linetype_manual(values = c("During removal" = 2, "Recovery period" = 1)) +
  scale_fill_manual(values = c("During removal" = "#FFFFFF", "Recovery period" = under_col)) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
  labs(x = "\U0394 kelp biomass (treatment - control)",
       y = "\U0394 understory algae biomass (treatment - control)") +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = c(0.83, 0.93),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

delta_algae_continual %>% 
  ggplot(aes(x = continual, y = continual_algae)) +
  geom_hline(aes(yintercept = 0), lty = 3, color = "grey") +
  geom_vline(aes(xintercept = 0), lty = 3, color = "grey") +
  geom_point(aes(shape = exp_dates, fill = exp_dates), size = 5, shape = 21) + 
  scale_fill_manual(values = c("during" = "#FFFFFF", after = under_col)) +
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
  labs(x = "kelp biomass",
       y = "understory algae biomass") +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

delta_algae_continual %>% 
  group_by(site, kelp_year) %>% 
  mutate(mean_delta_kelp = mean(delta_continual),
         mean_algae = mean(delta_continual_algae)) %>% 
  ggplot(aes(x = mean_delta_kelp, y = mean_algae)) +
  geom_hline(aes(yintercept = 0), lty = 3, color = "grey") +
  geom_vline(aes(xintercept = 0), lty = 3, color = "grey") +
  geom_point(size = 5, aes(fill = year, shape = site), alpha = 0.8) +
  scale_shape_manual(values = shape_palette_site) +
  geom_line(data = predicted_delta_algae_vs_kelp, aes(x = x, y = predicted), size = 2) +
  labs(x = "\U0394 kelp biomass (treatment - control)",
       y = "\U0394 understory algae biomass (treatment - control)",
       caption = "marginal R2 = 0.028, conditional R2 = 0.519") +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14),
        legend.position = c(0.9, 0.7),
        legend.text = element_text(size = 18),
        legend.title = element_blank())





