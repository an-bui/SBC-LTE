
##########################################################################-
# 0. set up ---------------------------------------------------------------
##########################################################################-

# only have to run this once per session
source(here::here("code", "01a-kelp_recovery.R"))

##########################################################################-
# 1. functions ------------------------------------------------------------
##########################################################################-

# Using function from the Thiault script below with edits to not print out results in the quartz viewer.

ProgressiveChangeBACIPS<-function(control, impact, time.true, time.model) 
{
  ### STEP 2 - Calculate the delta at each sampling date
  delta <- impact - control
  
  # Plot delta against time.true
  # dev.new(width=10, height=5)
  # par(mfrow=c(1,2))
  plot(delta~time.true, type="n")
  time.model.of.impact=max(which(time.model==0))
  rect(time.model.of.impact, min(delta)-100, max(time.model)+10, max(delta)+100, col = "grey")
  points(delta~time.true, pch=24, bg="white", cex=2)
  
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
    funAsy<-function(parS, time.model)	(parS$M * time.model) / (parS$L + time.model) + parS$B
    residFun<-function(p, observed, time.model) observed + funAsy(p,time.model)
    parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=1)
    nls_ASY_out <- nls.lm(par=parStart, fn= residFun, observed=delta, time.model=time.model, control = nls.lm.control(maxfev = integer(), maxiter = 1000))
    foAsy<-delta~(M * time.model) / (L + time.model) + B
    startPar<-c(-coef(nls_ASY_out)[1], coef(nls_ASY_out)[2], coef(nls_ASY_out)[3])
    asyFit<-nls2(foAsy, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
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

# Function to add a column in the biomass data frame for each site called `time.model` that has the timesteps labelled the way the function wants them

time.model.fxn.new <- function(site, treatment, season_choice) {
  # select a data frame based on site and treatment
  df <- if(treatment == "annual") {
    delta_annual
  } else if(treatment == "continual") {
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
    dplyr::select(site, year, month, date, control, 7, exp_dates, contains("delta_"), time_since_end) %>% 
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

# function to run ProgressiveChangeBACIPS() using a data frame

biomass.pcbacips <- function(df) {
  ProgressiveChangeBACIPS(
    control = pull(df, 6), # control column
    impact = pull(df, 7), # treatment column
    time.true = pull(df, 5), # date column
    time.model = pull(df, 1) # time.model column
  )
}

##########################################################################-
# 2. models ---------------------------------------------------------------
##########################################################################-

# ⟞ a. Arroyo Quemado -----------------------------------------------------

aque_biomass_continual <- time.model.fxn.new("aque", "continual", "all")

aque_continual_bacips_results <- biomass.pcbacips(aque_biomass_continual)

aque_bacips_plot <- delta_continual %>% 
  filter(site == "aque" & exp_dates == "after") %>% 
  ggplot(aes(x = time_since_end, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_point(shape = aque_shape, fill = aque_col, size = 4) +
  geom_function(fun = function(x) -249.35+(15.8*x), size = 2) +
  labs(x = "Time since end of removal (years)", y = "\U0394 giant kelp biomass",
       title = aque_full) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16), 
        legend.position = "none") 
aque_bacips_plot 

# ⟞ b. Naples -------------------------------------------------------------

napl_biomass_continual <- time.model.fxn.new("napl", "continual", "all")

napl_continual_bacips_results <- biomass.pcbacips(napl_biomass_continual)

napl_bacips_plot <- delta_continual %>% 
  filter(site == "napl" & exp_dates == "after") %>% 
  ggplot(aes(x = time_since_end, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_point(shape = napl_shape, fill = napl_col, size = 4) +
  geom_function(fun = function(x) -196.417+(11.9*x), size = 2) +
  labs(x = "Time since end of removal (years)", y = "\U0394 giant kelp biomass",
       title = napl_full) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16), 
        legend.position = "none") 
napl_bacips_plot 

# ⟞ c. Mohawk -------------------------------------------------------------

mohk_biomass_continual <- time.model.fxn.new("mohk", "continual", "all")

mohk_continual_bacips_results <- biomass.pcbacips(mohk_biomass_continual)

mohk_bacips_plot <- delta_continual %>% 
  filter(site == "mohk" & exp_dates == "after") %>% 
  ggplot(aes(x = time_since_end, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_point(shape = mohk_shape, fill = mohk_col, size = 4) +
  geom_function(fun = function(x) -824.75+(54.85*x), size = 2) +
  labs(x = "Time since end of removal (years)", y = "\U0394 giant kelp biomass",
       title = mohk_full) +
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16), 
        legend.position = "none") 
mohk_bacips_plot 

# ⟞ d. Carpinteria --------------------------------------------------------

carp_biomass_continual <- time.model.fxn.new("carp", "continual", "all")

carp_continual_bacips_results <- biomass.pcbacips(carp_biomass_continual)

TukeyHSD(aov(delta_continual ~ exp_dates, data = carp_biomass_continual))

carp_bacips_plot <- delta_continual %>% 
  filter(site == "carp" & exp_dates == "after") %>% 
  ggplot(aes(x = time_since_end, y = delta_continual)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_point(shape = carp_shape, fill = carp_col, size = 4) +
  geom_function(fun = function(x) 0, size = 2) +
  labs(x = "Time since end of removal (years)", y = "\U0394 giant kelp biomass",
       title = carp_full) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16), 
        legend.position = "none") 
carp_bacips_plot 


##########################################################################-
# 3. manuscript tables ----------------------------------------------------
##########################################################################-

model_selection_summary_table <- cbind(
  model = list("Step", "Linear", "Asymptotic", "Sigmoid"),
  deg.free = aque_continual_bacips_results$aicc.test.results.df,
  aic.val = round(aque_continual_bacips_results$aicc.test.results.AIC, 2),
  aic.diff = round(aque_continual_bacips_results$aicc.test.results.diff, 2),
  aic.rl = round(aque_continual_bacips_results$aicc.test.results.RL, 2),
  aic.weights = round(aque_continual_bacips_results$aicc.test.results.aicWeights, 2)
) %>% 
  as.data.frame() %>% 
  mutate(site = aque_full) %>% 
  rbind(
    cbind(
      model = list("Step", "Linear", "Asymptotic", "Sigmoid"),
      deg.free = napl_continual_bacips_results$aicc.test.results.df,
      aic.val = round(napl_continual_bacips_results$aicc.test.results.AIC, 2),
      aic.diff = round(napl_continual_bacips_results$aicc.test.results.diff, 2),
      aic.rl = round(napl_continual_bacips_results$aicc.test.results.RL, 2),
      aic.weights = round(napl_continual_bacips_results$aicc.test.results.aicWeights, 2)
    ) %>% 
      as.data.frame() %>% 
      mutate(site = napl_full) 
  ) %>% 
  rbind(
    cbind(
      model = list("Step", "Linear", "Asymptotic", "Sigmoid"),
      deg.free = mohk_continual_bacips_results$aicc.test.results.df,
      aic.val = round(mohk_continual_bacips_results$aicc.test.results.AIC, 2),
      aic.diff = round(mohk_continual_bacips_results$aicc.test.results.diff, 2),
      aic.rl = round(mohk_continual_bacips_results$aicc.test.results.RL, 2),
      aic.weights = round(mohk_continual_bacips_results$aicc.test.results.aicWeights, 2)
    ) %>% 
      as.data.frame() %>% 
      mutate(site = mohk_full) 
  ) %>% 
  rbind(
    cbind(
      model = list("Step", "Linear", "Asymptotic", "Sigmoid"),
      deg.free = carp_continual_bacips_results$aicc.test.results.df,
      aic.val = round(carp_continual_bacips_results$aicc.test.results.AIC, 2),
      aic.diff = round(carp_continual_bacips_results$aicc.test.results.diff, 2),
      aic.rl = round(carp_continual_bacips_results$aicc.test.results.RL, 2),
      aic.weights = round(carp_continual_bacips_results$aicc.test.results.aicWeights, 2)
    ) %>% 
      as.data.frame() %>% 
      mutate(site = carp_full) 
  ) %>% 
  gt(groupname_col = "site") %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = aic.diff < 0.0001)
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>% 
  cols_label(model = "Model",
             deg.free = "Degrees of freedom",
             aic.val = "AIC",
             aic.diff = "\U0394 AIC",
             aic.rl = "Relative likelihood",
             aic.weights = "AIC weight") %>% 
  tab_options(table.font.names = "Times New Roman") 
model_selection_summary_table

# gtsave(model_selection_summary_table,
#        here::here("tables", "ms-tables", paste("tbl-S2_", today(), ".png", sep = "")),
#        vwidth = 500, vheight = 1000)

##########################################################################-
# 4. manuscript figures ---------------------------------------------------
##########################################################################-

bacips_plots <- (aque_bacips_plot + napl_bacips_plot) /
  (mohk_bacips_plot + carp_bacips_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20))
bacips_plots

# ggsave(here::here("figures", paste("bacips_plots-", today(), ".jpg", sep = "")), 
#        plot = bacips_plots, 
#        height = 8, width = 12, dpi = 150)







