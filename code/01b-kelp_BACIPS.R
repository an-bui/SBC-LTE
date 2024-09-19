
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. set up -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only have to run this once per session
source(here::here("code", "01a-kelp_recovery.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. functions -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to conduct progressive change before-after-
# control-impact paired design analysis following Thiault et al. 2017 
# (doi: 10.1111/2041-210X.12655) The `ProgressiveChangeBACIPS()` function is 
# taken from the material presented in the paper.

# I made some edits to the function, namely: comments to clarify what each part
# was doing and making sure the results did not print out in the quartz viewer.
# Then, I created a new function called `biomass.pcbacips` that would use
# `ProgressiveChangeBACIPS` with a data frame in a standard format.
# That data frame would be derived using `time.model.fxn.new` from the 
# `delta_continual` object created in the `01a-kelp_recovery.R` script.

ProgressiveChangeBACIPS <- function(control, impact, time.true, time.model) {
  
  ### STEP 2 - Calculate the delta at each sampling date
  delta <- impact - control
  
  # Plot delta against time.true
  # dev.new(width=10, height=5)
  # par(mfrow=c(1,2))
  plot(delta~time.true, type="n")
  
  # counts the number of observations at "time = 0" (i.e. during disturbance)
  time.model.of.impact = max(which(time.model == 0))
  rect(time.model.of.impact, 
       min(delta) - 100, 
       max(time.model) + 10, 
       max(delta) + 100, 
       col = "grey")
  points(delta ~ time.true, 
         pch = 24,
         bg = "white", 
         cex = 2)
  
  ### STEP 3 - Fit and compete models
  ## Create a 'period' variable: if the time model == 0, assign "before"
  ## otherwise, assign "after"
  ## this is the only predictor in the aov() call below
  period <- ifelse(time.model == 0, "Before", "After")
  
  ## Fit a step model
  # analysis of variance: comparing means
  step.Model<-aov(delta ~ period)
  
  ## Fit a linear model
  # linear model: delta as a function of time.model
  # 0s in "before", sequence of numbers in "after"
  linear.Model<-lm(delta ~ time.model)
  
  ## Fit an asymptotic model
  # Create an asymptotic function
  myASYfun<-function(delta, time.model)
  {
    funAsy <- function(parS, time.model)	(parS$M * time.model) / (parS$L + time.model) + parS$B
    residFun <- function(p, observed, time.model) observed + funAsy(p,time.model)
    parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=1)
    nls_ASY_out <- nls.lm(par=parStart, 
                          fn= residFun, 
                          observed=delta, 
                          time.model=time.model, 
                          control = nls.lm.control(maxfev = integer(), 
                                                   maxiter = 1000))
    foAsy <- delta~(M * time.model) / (L + time.model) + B
    startPar <- c(-coef(nls_ASY_out)[1], coef(nls_ASY_out)[2], coef(nls_ASY_out)[3])
    asyFit <- nls2(foAsy, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
    asyFit
  }
  # Fit the asymptotic model
  asymptotic.Model<-myASYfun(delta=delta,time.model=time.model)
  
  
  ## Fit a sigmoid model
  ## Create a sigmoid function
  mySIGfun<-function(delta, time.model)
  {
    funSIG<-function(parS, time.model)	(parS$M * (time.model/parS$L)^parS$K) / (1 + (time.model/parS$L) ^ parS$K) + parS$B
    residFun<-function(p, observed, time.model) observed + funSIG(p,time.model)
    parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=mean(time.model), K=5)
    nls_SIG_out <- nls.lm(par=parStart, fn= residFun, observed=delta, time.model=time.model, control = nls.lm.control(maxfev = integer(), maxiter = 1000))
    foSIG<-delta~(M * (time.model/L) ^ K) / (1 + (time.model/L) ^ K) + B
    startPar<-c(-coef(nls_SIG_out)[1],-coef(nls_SIG_out)[2],coef(nls_SIG_out)[3],coef(nls_SIG_out)[4])
    sigFit<-nls2(foSIG, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
    sigFit
  }
  # Fit the sigmoid model
  sigmoid.Model<-mySIGfun(delta=delta,time.model=time.model)
  
  
  ## Compete models
  # Perform AIC tests
  AIC.test<-AIC(step.Model, linear.Model, asymptotic.Model, sigmoid.Model)
  AICc.test<-as.data.frame(cbind(AIC.test[,1], c(AICc(step.Model), AICc(linear.Model), AICc(asymptotic.Model), AICc(sigmoid.Model))))
  rownames(AICc.test)<-rownames(AIC.test)
  names(AICc.test)<-names(AIC.test)
  
  # Calculate AICc weight and selected the best model
  for(i in 1:dim(AICc.test)[1])
  {
    AICc.test$diff[i]<-AICc.test$AIC[i]-min(AICc.test$AIC)
  }
  AICc.test$RL<-exp(-0.5* AICc.test$diff)
  RL_sum<-sum(AICc.test$RL)
  AICc.test$aicWeights<-(AICc.test$RL/RL_sum)*100
  w<-AICc.test$aicWeights
  names(w)<-rownames(AICc.test)
  
  # Display raw AIC values
  print(AICc.test)
  
  # Display AICc weights
  # print(w)
  barplot(w, col="white", ylab="Relative likelihood (%)", cex.names = 0.9, names.arg =c("Step","Linear","Asymptotic","Sigmoid"))
  best.Model<-which(w==max(w))
  
  model.formula <- if(best.Model==1) {
    step.Model
  } else if(best.Model==2) {
    linear.Model
  } else if(best.Model==3) {
    asymptotic.Model
  } else if(best.Model==4) {
    sigmoid.Model
  }
  
  message <- if(best.Model==1) {
    "Step"
  } else if(best.Model==2) {
    "Linear"
  } else if(best.Model==3) {
    "Asymptotic"
  } else if(best.Model==4) {
    "Sigmoid"
  }
  
  likelihood <- if(best.Model==1) {
    paste(round(w[1],1), "%", sep="")
  } else if(best.Model==2) {
    paste(round(w[2],1), "%", sep="")
  } else if(best.Model==3) {
    paste(round(w[3],1), "%", sep="")
  } else if(best.Model==4) {
    paste(round(w[4],1), "%", sep="")
  }
  
  ### STEP 4 - Derive inference based on the best model (i.e., with the higher AICc weight)
  if(best.Model==1) {writeLines(paste("\n\nSTEP MODEL SELECTED - Likelihood = ", round(w[1],1), "%\n\n", sep=""))
    print(summary(step.Model))}
  if(best.Model==2) {writeLines(paste("\n\nLINEAR MODEL SELECTED - Likelihood = ", round(w[2],1), "%\n\n", sep=""))
    print(summary(linear.Model))}
  if(best.Model==3) {writeLines(paste("\n\nASYMPTOTIC MODEL SELECTED - Likelihood = ", round(w[3],1), "%\n\n", sep=""))
    print(asymptotic.Model)}
  if(best.Model==4) {writeLines(paste("\n\nSIGMOID MODEL SELECTED - Likelihood = ", round(w[4],1), "%\n\n", sep=""))
    print(sigmoid.Model)}
  
  return(c(message = message,
           likelihood = likelihood,
           aicc.test.results = AICc.test, 
           step.model.summary = summary(step.Model), 
           linear.model.summary = summary(linear.Model), 
           asy.model.summary = asymptotic.Model, 
           sigmoid.model.summary = sigmoid.Model))
}



# This is a function to wrangle the `delta_continual` data frame into a format
# that fits what `ProgressiveChangeBACIPS` expects from its input.

time.model.fxn.new <- function(site, treatment, season_choice) {
  # select a data frame based on site and treatment
  df <- if(treatment == "continual") {
    delta_continual
  } else {
    warning("Check your arguments! You might be missing site, treatment, or season.")
  }
  
  season_choice <- if(season_choice == "all") {
    c("winter", "spring", "summer", "fall")
  } else if(season_choice == "summer") {
    c("summer")
  } else if (season_choice == "fall") {
    c("fall")
  } else if (season_choice == "winter") {
    c("winter")
  } else if (season_choice == "spring") {
    c("spring")
  } else {
    warning("Check your arguments! You might be missing site, treatment, or season.")
    return(NA)
  }
  
  df %>% 
    filter(site == {{ site }}) %>% 
    dplyr::select(site, year, month, date, 
                  control, continual, 
                  exp_dates, contains("delta_"), time_since_end) %>% 
    mutate(exp_dates = fct_relevel(exp_dates, c("after", "during"))) %>% 
    arrange(exp_dates) %>% 
    # assign everything a "time step number" using the row numbers...
    rownames_to_column("time.model") %>% 
    # only keep the time step numbers for the "after dates"
    # make everything else 0 (i.e. everything else is before the "after")
    mutate(time.model = case_when(
      exp_dates == "during" ~ 0,
      TRUE ~ as.numeric(as.numeric(time.model))
    )) %>% 
    # rearrange the sampling dates to be in logical order
    mutate(exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
    arrange(exp_dates) 
}

# This is a function to run ProgressiveChangeBACIPS with the data frame created
# from `time.model.fxn.new`.

biomass.pcbacips <- function(df) {
  ProgressiveChangeBACIPS(
    control = pull(df, 6), # control column
    impact = pull(df, 7), # treatment column
    time.true = pull(df, 5), # date column
    time.model = pull(df, 1) # time.model column
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 2. model selection --------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. Arroyo Quemado -----------------------------------------------------

aque_biomass_continual <- time.model.fxn.new("aque", "continual", "all")

aque_continual_bacips_results <- biomass.pcbacips(aque_biomass_continual)

# Best model is linear model.

# fitting linear model to extract model predictions
aque_bacips_linear <- lm(delta_continual ~ time.model, 
                         data = aque_biomass_continual)

aque_bacips_predictions <- ggpredict(aque_bacips_linear, 
                                     terms = "time.model",
                                     type = "fixed")

# fitting model with time since end as a predictor for visualization
aque_time_since_end <- lm(delta_continual ~ time_since_end, 
                                 data = aque_biomass_continual %>% 
                                   filter(exp_dates == "after"))

aque_time_since_end_predictions <- ggpredict(aque_time_since_end, 
                                     terms = "time_since_end",
                                     type = "fixed")

# ⟞ b. Naples -------------------------------------------------------------

napl_biomass_continual <- time.model.fxn.new("napl", "continual", "all")

napl_continual_bacips_results <- biomass.pcbacips(napl_biomass_continual)

# Equal support for step, linear, and sigmoid models. Visualizations depict 
# step model (immediate increase in removal plot relative to reference), which 
# has the lowest AIC.

# ⟞ c. Mohawk -------------------------------------------------------------

mohk_biomass_continual <- time.model.fxn.new("mohk", "continual", "all")

mohk_continual_bacips_results <- biomass.pcbacips(mohk_biomass_continual)

# Equal support for linear and sigmoid model. Visualizations depict linear 
# model (gradual increase in removal plot relative to reference), which has the
# lowest AIC.

# fitting linear model to extract model predictions
mohk_bacips_linear <- lm(delta_continual ~ time.model, 
                         data = mohk_biomass_continual)

mohk_bacips_predictions <- ggpredict(mohk_bacips_linear, 
                                     terms = "time.model",
                                     type = "fixed")

# fitting model with time since end as a predictor for visualization
mohk_time_since_end <- lm(delta_continual ~ time_since_end, 
                          data = mohk_biomass_continual %>% 
                            filter(exp_dates == "after"))

mohk_time_since_end_predictions <- ggpredict(mohk_time_since_end, 
                                             terms = "time_since_end",
                                             type = "fixed")

# ⟞ d. Carpinteria --------------------------------------------------------

carp_biomass_continual <- time.model.fxn.new("carp", "continual", "all")

carp_continual_bacips_results <- biomass.pcbacips(carp_biomass_continual)

# equal support for step and linear models. Visualizations depict linear model,
# which has the lowest AIC.

# fitting linear model to extract model predictions
carp_bacips_linear <- lm(delta_continual ~ time.model, 
                         data = carp_biomass_continual)

carp_bacips_predictions <- ggpredict(carp_bacips_linear, 
                                     terms = "time.model",
                                     type = "fixed")

# fitting model with time since end as a predictor for visualization
carp_time_since_end <- lm(delta_continual ~ time_since_end, 
                          data = carp_biomass_continual %>% 
                            filter(exp_dates == "after"))

carp_time_since_end_predictions <- ggpredict(carp_time_since_end, 
                                     terms = "time_since_end",
                                     type = "fixed")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 3. model visualizations -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section contains code to create figures in Supplemental Material ____.
# For each site, there are two figures depicting delta biomass on the y-axis. 
# The first figure (sections i) has `time.model`, or the timesteps in the
# original `ProgressiveChangeBACIPS` function, as the x-axis. The second figure
# (sections ii) has `time_since_end`, or the time since the end of the removal
# in years as the x-axis. We created two figures to represent the results from
# the Progressive Change BACIPS procedure (using `time.model`), but also to 
# represent the observations as they appear in all other timeseries-style plots
# in the paper (with time since the end of removal).

# For all figures depicting timesteps (`time.model`), the "during" removal
# points are shown at the 0 point on the x-axis to represent the way the
# Progressive Change BACIPS procedure handles observations during a 
# disturbance.

# ⟞ a. figure aesthetics --------------------------------------------------

bacips_theme <- theme_bw() + 
    theme(axis.title = element_text(size = 15),
          plot.title = element_text(size = 18),
          plot.title.position = "plot",
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16), 
          legend.position = "none",
          panel.grid = element_blank()) 

bacips_hline <- geom_hline(yintercept = 0, lty = 2, alpha = 0.5)

bacips_labs <- list(
  labs(x = "Timesteps after removal", 
       y = "\U0394 giant kelp biomass\n(removal \U2212 reference, dry g/m\U00B2)")
)

time_since_end_labs <- list(
  labs(x = "Time since end of removal (years)", 
       y = "\U0394 giant kelp biomass\n(removal \U2212 reference, dry g/m\U00B2)")
)

bacips_aesthetics <- list(
  scale_alpha_discrete(range = c("during" = 0.3, "after" = 1))
)

# ⟞ b. Arroyo Quemado -----------------------------------------------------

# ⟞ ⟞ i. time.model -------------------------------------------------------

aque_bacips_plot <- aque_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2017-08-16" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time.model, y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, fill = aque_col, size = 4) +
  geom_line(data = aque_bacips_predictions,
            aes(x = x, y = predicted),
            linewidth = 1) + 
  geom_ribbon(data = aque_bacips_predictions,
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
              linewidth = 1, 
              alpha = 0.1) +
  bacips_aesthetics + 
  bacips_labs + 
  labs(title = "(a) Arroyo Quemado") +
  bacips_theme
aque_bacips_plot 

# ⟞ ⟞ ii. time_since_end --------------------------------------------------

aque_time_since_end <- aque_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2017-08-16" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time_since_end_model, 
             y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, 
             fill = aque_col, 
             size = 4) +
  geom_line(data = aque_time_since_end_predictions,
            aes(x = x, 
                y = predicted),
            linewidth = 1) + 
  geom_ribbon(data = aque_time_since_end_predictions,
              aes(x = x, 
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high),
              linewidth = 1, 
              alpha = 0.1) +
  bacips_aesthetics +
  time_since_end_labs + 
  labs(title = "(a) Arroyo Quemado") +
  bacips_theme

aque_time_since_end

# ⟞ c. Naples -------------------------------------------------------------

# ⟞ ⟞ i. time.model -------------------------------------------------------

napl_bacips_plot <- napl_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2016-05-17" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time_since_end_model, y = delta_continual)) +
  geom_point(aes(alpha = exp_dates), 
             fill = napl_col, size = 4, shape = 21) +
  geom_hline(yintercept = 0, linewidth = 1) +
  # scale_x_continuous(breaks = time_since_end) +
  bacips_aesthetics + 
  bacips_labs + 
  labs(title = "(b) Naples") +
  bacips_theme

napl_bacips_plot 

# ⟞ ⟞ ii. time_since_end --------------------------------------------------

napl_time_since_end <- napl_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2016-05-17" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time_since_end_model, 
             y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, 
             fill = napl_col, 
             size = 4) +
  geom_hline(yintercept = 0, linewidth = 1) +
  bacips_aesthetics +
  time_since_end_labs + 
  labs(title = "(b) Naples") +
  bacips_theme

napl_time_since_end


# ⟞ d. Mohawk -------------------------------------------------------------

# ⟞ ⟞ i. time.model -------------------------------------------------------

mohk_bacips_plot <- mohk_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2017-08-16" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time.model, y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, fill = mohk_col, size = 4) +
  geom_line(data = mohk_bacips_predictions,
            aes(x = x, y = predicted),
            linewidth = 1) + 
  geom_ribbon(data = mohk_bacips_predictions,
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
              linewidth = 1, 
              alpha = 0.1) +
  bacips_aesthetics + 
  bacips_labs + 
  labs(title = "(c) Mohawk") +
  bacips_theme
mohk_bacips_plot 

# ⟞ ⟞ ii. time_since_end --------------------------------------------------

mohk_time_since_end <- mohk_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2017-08-16" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time_since_end_model, 
             y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, 
             fill = mohk_col, 
             size = 4) +
  geom_line(data = mohk_time_since_end_predictions,
            aes(x = x, 
                y = predicted),
            linewidth = 1) + 
  geom_ribbon(data = mohk_time_since_end_predictions,
              aes(x = x, 
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high),
              linewidth = 1, 
              alpha = 0.1) +
  bacips_aesthetics +
  time_since_end_labs + 
  labs(title = "(c) Mohawk") +
  bacips_theme

mohk_time_since_end

# ⟞ e. Carpinteria --------------------------------------------------------

# ⟞ ⟞ i. time.model -------------------------------------------------------

carp_bacips_plot <- carp_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2017-08-16" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time.model, y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, fill = carp_col, size = 4) +
  geom_line(data = carp_bacips_predictions,
            aes(x = x, y = predicted),
            linewidth = 1) + 
  geom_ribbon(data = carp_bacips_predictions,
              aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
              linewidth = 1, 
              alpha = 0.1) +
  bacips_aesthetics + 
  bacips_labs + 
  labs(title = "(d) Carpinteria") +
  bacips_theme
carp_bacips_plot 

# ⟞ ⟞ ii. time_since_end --------------------------------------------------

carp_time_since_end <- carp_biomass_continual %>% 
  mutate(time_since_end_model = case_when(
    date < "2017-08-16" ~ 0,
    TRUE ~ time_since_end
  )) %>% 
  ggplot(aes(x = time_since_end_model, 
             y = delta_continual)) +
  bacips_hline +
  geom_point(aes(alpha = exp_dates), 
             shape = 21, 
             fill = carp_col, 
             size = 4) +
  geom_line(data = carp_time_since_end_predictions,
            aes(x = x, 
                y = predicted),
            linewidth = 1) + 
  geom_ribbon(data = carp_time_since_end_predictions,
              aes(x = x, 
                  y = predicted, 
                  ymin = conf.low,
                  ymax = conf.high),
              linewidth = 1, 
              alpha = 0.1) +
  bacips_aesthetics +
  time_since_end_labs + 
  labs(title = "(d) Carpinteria") +
  bacips_theme

carp_time_since_end


# ⟞ f. saving output ------------------------------------------------------

# ⟞ ⟞ i. time.model -------------------------------------------------------

bacips_plots <- plot_grid(
  aque_bacips_plot, napl_bacips_plot,
  mohk_bacips_plot, carp_bacips_plot,
  nrow = 2, ncol = 2
)

bacips_plots

# ggsave(here::here("figures", "ms-figures", paste("bacips_plots-", today(), ".jpg", sep = "")),
#        plot = bacips_plots,
#        height = 8, width = 12, dpi = 300)

# ⟞ ⟞ ii. time_since_end --------------------------------------------------

time_since_end_plots <- plot_grid(
  aque_time_since_end, napl_time_since_end,
  mohk_time_since_end, carp_time_since_end,
  nrow = 2, ncol = 2
)

time_since_end_plots

# ggsave(here::here("figures", "ms-figures", paste0("time_since_end_plots-", today(), ".jpg")),
#        plot = time_since_end_plots,
#        height = 8, width = 12, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 4. tables --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to create the summary table in Supplmentary 
# Material _____.

# ⟞ a. table formatting ---------------------------------------------------

# function to organize model selection statistics into a data frame
bacips_table <- function(site) {
  
  # select the right output for a given site
  results <- if(site == "aque") {
    aque_continual_bacips_results
  } else if(site == "napl") {
    napl_continual_bacips_results
  } else if(site == "mohk") {
    mohk_continual_bacips_results
  } else if(site == "carp") {
    carp_continual_bacips_results
  } else {
    warning("Check your arguments! You might be missing site.")
  }
  
  # select the right name for a given site
  site_full <- if(site == "aque") {
    "Arroyo Quemado"
  } else if(site == "napl") {
    "Naples"
  } else if(site == "mohk") {
    "Mohawk"
  } else if(site == "carp") {
    "Carpinteria"
  } else {
    warning("Check your arguments! You might be missing site.")
  }
  
  # create the data frame by extracting components from the model selection 
  # results objects created in section 2: model selection
  bind_rows(
    model = c("Step", "Linear", "Asymptotic", "Sigmoid"),
    deg.free = results$aicc.test.results.df,
    aic.val = round(results$aicc.test.results.AIC, 2),
    aic.diff = round(results$aicc.test.results.diff, 2),
    aic.rl = round(results$aicc.test.results.RL, 2),
    aic.weights = round(results$aicc.test.results.aicWeights, 2)
  ) %>% 
    mutate(site = site_full) %>% 
    relocate(site, .before = model)
}

# create the table
model_selection_summary_table <- bind_rows(
  bacips_table("aque"),
  bacips_table("napl"),
  bacips_table("mohk"),
  bacips_table("carp")
) %>% 
  flextable::flextable() %>% 
  merge_v(j = ~ site) %>% 
  valign(j = ~ site,
         i = NULL,
         valign = "top") %>% 
  style(i = ~ aic.diff < 2,
        j = c("model", "aic.diff"),
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  set_header_labels(site = "Site",
                    model = "Model",
                    deg.free = "Degrees of freedom",
                    aic.val = "AIC",
                    aic.diff = "\U0394 AIC",
                    aic.rl = "Relative likelihood",
                    aic.weights = "AIC weight") %>% 
  autofit %>% 
  fit_to_width(7) %>% 
  font(fontname = "Times New Roman", part = "all") %>% 
  fontsize(size = 11, part = "all")


# ⟞ b. saving output ------------------------------------------------------

# model_selection_summary_table %>%
#   save_as_docx(path = here::here(
#     "tables", 
#     "ms-tables", 
#     paste("tbl-S2_", today(), ".docx", sep = "")))

