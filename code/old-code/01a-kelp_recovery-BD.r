# This just sources in your 01a-kelp_recovery.R file

df <- delta_continual %>% 
  pivot_longer(cols = c(control, continual)) %>% # Here I've pushed it to long form so that I can model control vs. continual removel as a categorical variable.
  filter(exp_dates == "during") %>% # Just the during time period
  mutate(value.adjusted = value + 0.001) # add a litte to be able to use a Gamma (I don't end up using this) 


lm1 <- lmer(value ~ time_since_end*name + (1|site),
  data = df) # So for this model I've included an interaction between your core predictor (time_since_end) and "name" which is just a dummy for control vs. continual 
summary(lm1) 
model_checks <- DHARMa::simulateResiduals(lm1) 
plot(model_checks) # pretty damn bad

lm2 <- glmmTMB::glmmTMB(value.adjusted ~ scale(time_since_end)*name + (1|site),
            data = df, family = Gamma(link = "log")) 
summary(lm2)
model_checks <- DHARMa::simulateResiduals(lm2) # pretty damn bad
plot(model_checks)

lm3 <- glmmTMB::glmmTMB(value ~ time_since_end*name + (1|site),
                        data = df, family = ziGamma(link = "log"), ziformula = ~.) # So this is a Gamma hurdle model. I've included all the covariates in the zero-component. But you could play aroudn with this if you wanted. 
summary(lm3)
model_checks <- DHARMa::simulateResiduals(lm3) # pretty good! Not perfect but I feel like this is acceptable
plot(model_checks)

lm3_plot <- plot(ggpredict(lm3, terms = ~time_since_end*name), add.data = T) # Fast way to get a plot

m.emm <- emmeans::emmeans (lm2,  ~ name | time_since_end, at = list(time_since_end = seq(min(df$time_since_end), max(df$time_since_end), length.out = 100))) # This allows you to actually generate the predicted marginal means. 
df_after <- delta_continual %>% 
  pivot_longer(cols = c(control, continual)) %>% # Here I've pushed it to long form so that I can model control vs. continual removel as a categorical variable.
  filter(exp_dates == "after")

lm4 <- glmmTMB::glmmTMB(value ~ scale(time_since_end)*name + (1|site),
                        data = df_after, family = ziGamma(link = "log"), ziformula = ~.) 
plot(simulateResiduals(lm4))
lm4_plot <- plot(ggpredict(lm4, terms = ~time_since_end*name), add.data = T)

df_long <- df_after <- delta_continual %>% 
  pivot_longer(cols = c(control, continual))

lm5 <- glmmTMB(value ~ exp_dates*name + (1|site),
               data = df_long, family = ziGamma(link = "log"), ziformula = ~.)
plot(simulateResiduals(lm5))
summary(lm5)
lm5_plot <- plot(ggpredict(lm5, terms = ~exp_dates*name), add.data = T)
lm5_plot

lm6 <- glmmTMB(value ~ scale(time_since_end)*exp_dates*name + (1|site),
               data = df_long, family = ziGamma(link = "log"), ziformula = ~.)
plot(simulateResiduals(lm6))
summary(lm6)
lm6_plot <- plot(ggpredict(lm6, terms = ~scale(time_since_end)*name*exp_dates), add.data = T) 
lm6_plot

emmeans::contrast(m.emm, "pairwise", type = "response", infer = T) %>% # This gives you the contrasts. Because the model has a log-link. The output is the log-ratio (log(X) - log(y) = log(x/y))
  broom::tidy(conf.int = T) %>%
  ggplot(aes(x = as.numeric(time_since_end), y = ratio))+
  geom_line()+
  # geom_point(data = df, aes(x = time_since_end, y = ))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25)


# Not sure this is correct. I believe you still have to force time_since_end into a factor to run the autoregressive model correctly. This seems real complicated and I wouldn't necessarily recommend this approach unless Max has more experience.
lm4 <- glmmTMB::glmmTMB(value ~ ar1(scale(time_since_end)*name + 0 | site),
                        data = df, family = ziGamma(link = "log"), ziformula = ~.) 
summary(lm4)
model_checks <- DHARMa::simulateResiduals(lm4) # pretty good
plot(model_checks)

AIC(lm3, lm4)



df2 <- delta_continual %>% 
  pivot_longer(cols = c(control, continual)) %>%
  filter(exp_dates == "after")

lm5 <- glmmTMB::glmmTMB(value ~ time_since_end*name + (1|site),
                        data = df2, family = ziGamma(link = "log"), ziformula = ~.) 
summary(lm5)
model_checks <- DHARMa::simulateResiduals(lm5) # pretty good?
plot(model_checks)

lm5_plot <- plot(ggpredict(lm5, terms = ~time_since_end*name), add.data = T)

lm3_plot + lm5_plot
delta_calc_during <- ggpredict(lm3, terms = ~time_since_end*name) %>% 
  select(group, predicted, x) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta_calc = continual - control)

delta_calc_after <- ggpredict(lm5, terms = ~time_since_end*name) %>% 
  select(group, predicted, x) %>% 
  pivot_wider(names_from = group, values_from = predicted) %>% 
  mutate(delta_calc = continual - control)

ggplot(data = delta_continual, aes(x = time_since_end, y = delta_continual)) +
  geom_point() +
  geom_line(data = delta_calc_during, aes(x = x, y = delta_calc)) +
  geom_line(data = delta_calc_after, aes(x = x, y = delta_calc))


