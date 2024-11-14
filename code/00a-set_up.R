
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. libraries -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. general organization and cleaning ----------------------------------

library(tidyverse)
library(here)
library(janitor)
library(lubridate)

# ⟞ b. visualization ------------------------------------------------------

library(cowplot)
library(patchwork)

# ⟞ c. model tools --------------------------------------------------------

library(minpack.lm) # Fitting non-linear models
library(nls2) # Fitting non-linear models
library(AICcmodavg) # calculate second order AIC (AICc)
library(MuMIn)
library(boot)
library(lmerTest) # also loads `lme4`
library(glmmTMB) 
library(DHARMa)
library(performance)

# ⟞ d. model predictions, means, etc. -------------------------------------

library(ggeffects)
library(broom.mixed)

# ⟞ e. community analysis -------------------------------------------------

library(vegan)

# ⟞ f. table making -------------------------------------------------------

library(flextable)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 2. start and end dates ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. Arroyo Quemado (AQUE) ---------------------------------------------

# date of first removal: 2010-01-29
# first survey after first removal: 2010-04-26
aque_start_date_continual <- as_date("2010-04-26")
# date of last removal: 2017-03-02
# first survey after last removal: 2017-05-18
# start of recovery period: 2017-08-16
aque_after_date_continual <- as_date("2017-08-16")


# ⟞ b. Naples (NAPL) -----------------------------------------------------

# date of first removal: 2010-01-28
# first survey after first removal: 2010-04-27
napl_start_date_continual <- as_date("2010-04-27")

# date of last removal: 2016-02-09
# first survey after last removal: 2016-05-17
# start of recovery period: 2016-08-16
napl_after_date_continual <- as_date("2016-08-16")

# ⟞ c. Mohawk (MOHK) ------------------------------------------------------

# date of first removal: 2010-05-05
# first survey after first removal: 2010-06-14
mohk_start_date_continual <- as_date("2010-07-23")

# date of last removal: 2017-02-13
# first survey after last removal: 2017-05-17
# start of recovery period: 2017-08-11
mohk_after_date_continual <- as_date("2017-08-11")

# ⟞ d. Carpinteria (CARP)  ------------------------------------------------

# date of first removal: 2010-02-04 
# first survey after first removal: 2010-04-23
carp_start_date_continual <- as_date("2010-04-23")
# first survey after first removal: 2010-03-11

# date of last removal: 2017-02-15
# first survey after last removal: 2017-05-19
# start of recovery period: 2017-08-10
carp_after_date_continual <- as_date("2017-08-10")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------- 3. useful wrangling functions ---------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a new column for season and set factor levels
season_column <- function(df) {
  df %>% 
    mutate(season = case_when(
      month %in% c(12, 1, 2, 3) ~ "winter",
      month %in% c(4, 5) ~ "spring",
      month %in% c (6, 7, 8) ~ "summer",
      month %in% c(9, 10, 11) ~ "fall"
    ),
    season = fct_relevel(season, "spring", "summer", "fall", "winter")) 
}

# function to calculate standard error
se <- function(x,...){
  sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
}

time_since_columns_continual <- function(df) {
  df %>% 
    # create a column for quarter
    mutate(quarter = case_when(
      month <= 3 ~ "Q1",
      month >= 4 & month <= 6 ~ "Q2",
      month >= 7 & month <= 9 ~ "Q3",
      month >= 10 & month <= 12 ~ "Q4"
    )) %>% 
    # calculate time since start of experiment
    mutate(time_yrs = case_when(
      # AQUE, NAPL, MOHK, CARP: continual started in 2010
      quarter == "Q1" ~ year + 0.125 - 2010,
      quarter == "Q2" ~ year + 0.375 - 2010,
      quarter == "Q3" ~ year + 0.625 - 2010,
      quarter == "Q4" ~ year + 0.875 - 2010
    )) %>% 
    group_by(site) %>% 
    mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
    ungroup() %>% 
    # calculate time since end of experiment
    group_by(site, exp_dates) %>% 
    mutate(time_since_end = case_when(
      exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
      exp_dates == "after" ~ time_yrs - min(time_yrs)
    )) %>% 
    mutate(test_min_time_yrs = min(time_yrs)) %>% 
    ungroup()
}

# kelp year
# expects data frame that already has quarter from time_since_columns functions
kelp_year_column <- function(df) {
  df %>% 
  # create a new column for "kelp year"
  mutate(quarter = fct_relevel(quarter, "Q2", "Q3", "Q4", "Q1")) %>% 
    mutate(kelp_year = case_when(
      quarter %in% c("Q2", "Q3", "Q4") ~ year,
      quarter == "Q1" ~ year - 1
    )) %>% 
    mutate(kelp_year = paste("kelp_", kelp_year, "-", kelp_year + 1, sep = "")) 
}

comparison_column_continual <- function(df) {
  df %>% 
    mutate(comp_1yrs = case_when(
      site %in% c("aque", "mohk", "carp") & between(time_since_end, -7.25, -6.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -5.25) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -1.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 4.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 5.75, 6.75) ~ "after"
      
    )) %>% 
    mutate(comp_1yrs = fct_relevel(comp_1yrs, "start", "during", "after")) %>% 
    mutate(comp_2yrs = case_when(
      site %in% c("aque", "mohk", "carp") & between(time_since_end, -7.25, -5.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -4.25) ~ "start",
      site == "mohk" & between(time_since_end, -7.25, -5.25) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -2.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 3.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 4.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
    mutate(comp_3yrs = case_when(
      site %in% c("aque", "mohk", "carp") & between(time_since_end, -7.25, -4.25) ~ "start",
      site == "napl" & between(time_since_end, -6.25, -3.25) ~ "start",
      site == "mohk" & between(time_since_end, -7.25, -4.25) ~ "start",
      
      site %in% c("aque", "napl", "mohk", "carp") & between(time_since_end, -3.25, -0.25) ~ "during",
      
      site %in% c("aque", "mohk", "carp") & between(time_since_end, 2.75, 5.75) ~ "after",
      site == "napl" & between(time_since_end, 3.75, 6.75) ~ "after"
    )) %>% 
    mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>% 
    # unite("sample_ID", site, date, quarter, remove = FALSE) %>% 
    unite("sample_ID_short", site, date, remove = FALSE)
}

# function to extract model summaries
model_summary_fxn <- function(model) {
  model %>% 
    # use tidy to get model summary and calculate 95% CI
    tidy(conf.int = TRUE) %>% 
    # only include fixed conditional effects
    filter(effect == "fixed" & component == "cond") %>%
    select(term, estimate, p.value, conf.low, conf.high) %>% 
    # create a new column that indicates whether an effect is significant
    mutate(signif = case_when(
      p.value <= 0.05 ~ "yes",
      TRUE ~ "no"
    )) %>% 
    # create a p-value column that converts very small values to < 0.001
    # and rounds all other values to relevant digits
    mutate(p.value = case_when(
      between(p.value, 0, 0.001) ~ "<0.001",
      between(p.value, 0.001, 0.01) ~ as.character(round(p.value, digits = 3)),
      between(p.value, 0.01, 1) ~ as.character(round(p.value, digits = 2))
    )) %>%
    # round other numeric values to two digits
    mutate(across(where(is.numeric), ~ round(., digits = 2))) %>%
    # create a confidence interval column
    unite(ci_interval, conf.low, conf.high, sep = ", ") %>%
    # rename the terms to be neater
    mutate(term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "time_since_end" ~ "Time since end",
      term == "treatmentremoval" ~ "Treatment (removal)",
      term == "time_since_end:treatmentremoval" ~ "Time since end × treatment (removal)"
    ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------------- 4. data --------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. Max's guild data ---------------------------------------------------

guilds <- read_csv(here::here("code", "resources", "castorani", "LTE_guild_data.csv")) %>% 
  mutate(sp.code = replace_na(sp.code, "Nandersoniana")) %>% 
  rename("new_group" = biomass.guild) %>% 
  # for whatever reason Yellowtail Rockfish are not in the guild csv
  add_row(sp.code = "SFLA", new_group = "fish", diversity.guild = "fish")


# ⟞ b. LTE all species biomass --------------------------------------------

biomass <- read_csv(here::here("data",
                               "all-species-biomass",
                               "knb-lter-sbc.119.12",
                               "LTE_All_Species_Biomass_at_transect_20241014.csv")) %>%
  clean_names() %>%
  # replace NA sp_code with Nandersoniana
  mutate(sp_code = case_when(
    scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
    TRUE ~ sp_code
  )) %>%
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 == -99999, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 == -99999, NA),
         density = replace(density, density == -99999, NA)) %>%
  # change to lower case
  mutate_at(c("group", "mobility", "growth_morph", "treatment", "site"), str_to_lower) %>%
  # create a sample_ID for each sampling date at each treatment at each site
  unite("sample_ID", site, treatment, date, remove = FALSE) %>%
  # filter to only include continual removal plots and control plots
  filter(treatment %in% c("continual", "control")) %>%
  left_join(., guilds, by = c("sp_code" = "sp.code")) %>%
  mutate(exp_dates = case_when(
    site == "aque" & date >= aque_start_date_continual & date < aque_after_date_continual ~ "during",
    site == "aque" & date >= aque_after_date_continual ~ "after",
    site == "napl" & date >= napl_start_date_continual & date < napl_after_date_continual ~ "during",
    site == "napl" & date >= napl_after_date_continual ~ "after",
    site == "mohk" & date >= mohk_start_date_continual & date < mohk_after_date_continual ~ "during",
    site == "mohk" & date >= mohk_after_date_continual ~ "after",
    site == "carp" & date >= carp_start_date_continual & date < carp_after_date_continual ~ "during",
    site == "carp" & date >= carp_after_date_continual ~ "after"
  ),
  exp_dates = fct_relevel(exp_dates, "during", "after")) %>%
  # take out all surveys that were before the removal experiment started
  drop_na(exp_dates) %>%
  # add in time since the end of the experiment
  time_since_columns_continual() %>%
  # take average biomass for surveys that were done 2x seasonally
  group_by(site, year, treatment, quarter, sp_code) %>%
  mutate(dry_gm2 = mean(dry_gm2),
         wm_gm2 = mean(wm_gm2),
         density = mean(density)) %>%
  # take out the "duplicates": only one sampling date per quarter in the dataframe, with values averaged across the two sampling dates
  slice(1L) %>%
  ungroup() %>%
  # take out extraneous columns from time_since_columns_continual()
  select(!test_min_time_yrs) %>%
  # add the 1, 2, and 3 year comparisons
  comparison_column_continual() %>%
  # add the "kelp year" column
  kelp_year_column()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------- 5. useful vectors and data frames -------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. species -----------------------------------------------------------

# data frame of each species including scientific name, common name, and group
spp_names <- biomass %>% 
  dplyr::select(group, sp_code, scientific_name, common_name) %>% 
  unique() 


# ⟞ b. sites --------------------------------------------------------------

# vector of sites
LTE_sites <- c("aque", "napl", "ivee", "mohk", "carp")

site_quality <- tribble(
  ~site, ~quality,
  "mohk", "high",
  "ivee", "high",
  "aque", "medium",
  "napl", "medium",
  "carp", "low"
) %>% 
  mutate(quality = fct_relevel(quality, c("low", "medium", "high")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 7. plot aesthetics --------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. site colors and shapes --------------------------------------------
          
aque_col <- '#D46F10'
napl_col <- '#4CA49E'
mohk_col <- '#6B6D9F'
carp_col <- '#4B8FF7'

color_palette_site <- c("aque" = aque_col, 
                        "napl" = napl_col, 
                        "mohk" = mohk_col, 
                        "carp" = carp_col)

aque_shape <- 21
napl_shape <- 22
mohk_shape <- 23
carp_shape <- 24

shape_palette_site <- c("aque" = aque_shape,
                        "napl" = napl_shape,
                        "mohk" = mohk_shape,
                        "carp" = carp_shape)

removal_col <- "#CC7540"
reference_col <- "#6D5A18"
under_col <- "#6B6D9F"

model_predictions_theme <- theme_bw() +
  # coord_cartesian(ylim = c(30, 800)) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.position = "none",
        # legend.position = c(0.86, 0.88),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, units = "cm"),
        legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
        legend.spacing.y = unit(0.01, "cm"),
        panel.grid = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 10))

model_predictions_aesthetics <- list(
  scale_color_manual(values = c(reference = reference_col, 
                                removal = removal_col),
                     labels = c(reference = "Reference", 
                                removal = "Removal")),
  scale_linetype_manual(values = c(reference = "22", 
                                   removal = "solid"),
                        labels = c(reference = "Reference", 
                                   removal = "Removal")),
  scale_x_continuous(limits = c(-8, 7), 
                     breaks = seq(-8, 7, by = 1), 
                     minor_breaks = NULL),
  guides(color = guide_legend(keyheight = 0.6),
         shape = guide_legend(keyheight = 0.6),
         lty = guide_legend(keyheight = 0.6),
         keyheight = 1),
  labs(x = "Time since end of removal (years)", 
       y = "Biomass (dry g/m\U00B2)", 
       linetype = "Treatment",
       color = "Treatment",
       shape = "Treatment",
       size = "Treatment",
       group = "Treatment")
)

delta_aesthetics <- list(
  scale_x_continuous(limits = c(-8, 7), 
                     breaks = seq(-8, 7, by = 1), 
                     minor_breaks = NULL),
  labs(x = "Time since end of removal (years)", 
       y = "\U0394 biomass (removal \U2212 reference, dry g/m\U00B2)")
)

model_predictions_background <- list(
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey"),
    geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey"),
    annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.3)
)

# ⟞ b. site full names ----------------------------------------------------

aque_full <- "Arroyo Quemado"
napl_full <- "Naples"
mohk_full <- "Mohawk"
carp_full <- "Carpinteria"

sites_full <- setNames(c("Arroyo Quemado",
                         "Naples",
                         "Isla Vista",
                         "Mohawk",
                         "Carpinteria"),
                       c("aque",
                         "napl",
                         "ivee",
                         "mohk",
                         "carp"))

# ⟞ c. raw biomass plots --------------------------------------------------

# This theme is used to create the plots of biomass through time for giant
# kelp, understory macroalgae, and sessile invertebrates.

raw_biomass_plot_theme <-
    theme_bw() +
    theme(axis.title = element_text(size = 8),
          plot.title = element_text(size = 8),
          axis.text = element_text(size = 7),
          strip.text = element_text(size = 8, hjust = 0),
          strip.background = element_blank(),
          legend.text = element_text(size = 7),
          legend.position = "none",
          panel.grid = element_blank())

# ⟞ d. column titles ------------------------------------------------------

algae_title <- ggplot(data.frame(l = "Understory macroalgae", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

epi_title <- ggplot(data.frame(l = "Sessile invertebrates", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

kelp_title <- ggplot(data.frame(l = "Giant kelp", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

