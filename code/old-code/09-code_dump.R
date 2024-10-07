
# 1. kelp fronds ----------------------------------------------------------

# - `delta_fronds`: calculates difference in frond density between control and removal plots (removal - control)  
# - `control_fronds`: kelp fronds from control plots only  
# - `total_fronds`: kelp fronds from treatment + control plots

# total fronds

# total_fronds <- kelp_fronds %>% 
#   group_by(site, year, month, treatment, date) %>% 
#   summarize(total_fronds = sum(fronds)) %>% 
#   ungroup() %>% 
#   # # make a new column for after dates
#   # after_dates_column() %>% 
#   # make a new column for during and after and set factor levels
#   exp_dates_column() %>% 
#   # create a new column for season and set factor levels
#   season_column() %>% 
#   filter(total_fronds > -1) 
# 
# # delta fronds
# delta_fronds <- total_fronds %>% 
#   pivot_wider(names_from = treatment, values_from = total_fronds) %>% 
#   mutate(delta_annual = annual - control,
#          delta_continual = continual - control) 
# 
# # just control
# control_fronds <- total_fronds %>% 
#   filter(treatment == "control")


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


# 2. delta annual kelp biomass --------------------------------------------

# delta_annual <- delta_biomass %>% 
#   dplyr::select(site, year, month, date, control, annual, delta_annual) %>% 
#   drop_na(delta_annual) %>% 
#   mutate(exp_dates = case_when(
#     # after for annual removal:
#     site == "aque" & date >= aque_after_date_annual ~ "after",
#     site == "napl" & date >= napl_after_date_annual ~ "after",
#     site == "ivee" & date >= ivee_after_date_annual ~ "after",
#     site == "mohk" & date >= mohk_after_date_annual ~ "after",
#     site == "carp" & date >= carp_after_date_annual ~ "after",
#     # everything else is "during" the experiment
#     TRUE ~ "during"
#   ),
#   exp_dates = fct_relevel(exp_dates, c("during", "after"))) %>% 
#   # create a column for quarter
#   mutate(quarter = case_when(
#     month <= 3 ~ "Q1",
#     month <= 6 ~ "Q2",
#     month <= 9 ~ "Q3",
#     TRUE ~ "Q4"
#   )) %>% 
#   # calculate time since start of experiment
#   mutate(time_yrs = case_when(
#     # AQUE, NAPL, MOHK, CARP: control and annual started in 2008
#     site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q1" ~ year + 0.125 - 2008,
#     site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q2" ~ year + 0.375 - 2008,
#     site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q3" ~ year + 0.625 - 2008,
#     site %in% c("aque", "napl", "mohk", "carp") & quarter == "Q4" ~ year + 0.875 - 2008, 
#     # IVEE control and annual started in 2011
#     site == "ivee" & quarter == "Q1" ~ year + 0.125 - 2011,
#     site == "ivee" & quarter == "Q2" ~ year + 0.375 - 2011,
#     site == "ivee" & quarter == "Q3" ~ year + 0.625 - 2011,
#     site == "ivee" & quarter == "Q4" ~ year + 0.875 - 2011
#   )) %>% 
#   group_by(site) %>% 
#   mutate(time_since_start = time_yrs - min(time_yrs)) %>% 
#   ungroup() %>% 
#   # calculate time since end of experiment
#   group_by(site, exp_dates) %>% 
#   # if "after", then simple: the time in years - the minimum time in years
#   # if "during", then more complex: take the max time in years and add 0.25, then subtract the time in years
#   mutate(time_since_end = case_when(
#     exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
#     exp_dates == "after" ~ time_yrs - min(time_yrs)
#   )) %>% 
#   ungroup() %>% 
#   left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
#   rename(site_full = value) %>% 
#   mutate(site_full = fct_relevel(site_full, "Arroyo Quemado (AQUE)", "Naples (NAPL)", "Isla Vista (IVEE)", "Mohawk (MOHK)", "Carpinteria (CARP)")) %>% 
#   # create a new column for "kelp year"
#   mutate(quarter = fct_relevel(quarter, "Q2", "Q3", "Q4", "Q1")) %>% 
#   mutate(kelp_year = case_when(
#     quarter %in% c("Q2", "Q3", "Q4") ~ year,
#     quarter == "Q1" ~ year - 1
#   )) %>% 
#   mutate(kelp_year = paste("kelp_", kelp_year, "-", kelp_year + 1, sep = "")) %>% 
#   # create a column for the points to compare for "2 year interval"
#   mutate(comp_2yrs = case_when(
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010") ~ "start",
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during", 
#     site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
#     site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022") ~ "after"
#   )) %>% 
#   mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>% 
#   # create a column for the points to compare for "3 year interval"
#   mutate(comp_3yrs = case_when(
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2008-2009", "kelp_2009-2010", "kelp_2010-2011") ~ "start",
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during", 
#     site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
#     site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     kelp_year %in% c("kelp_2019-2020", "kelp_2020-2021", "kelp_2021-2022") ~ "after"
#   )) %>% 
#   mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after"))


# 3. delta continual kelp biomass -----------------------------------------

# delta_continual <- delta_biomass %>%
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
#   # create a column for quarter
#   mutate(quarter = case_when(
#     month <= 3 ~ "Q1",
#     month <= 6 ~ "Q2",
#     month <= 9 ~ "Q3",
#     TRUE ~ "Q4"
#   )) %>%
#   # calculate time since start of experiment
#   mutate(time_yrs = case_when(
#     # AQUE, NAPL, MOHK, CARP: continual started in 2010
#     quarter == "Q1" ~ year + 0.125 - 2010,
#     quarter == "Q2" ~ year + 0.375 - 2010,
#     quarter == "Q3" ~ year + 0.625 - 2010,
#     quarter == "Q4" ~ year + 0.875 - 2010
#   )) %>%
#   group_by(site) %>%
#   mutate(time_since_start = time_yrs - min(time_yrs)) %>%
#   ungroup() %>%
#   # calculate time since end of experiment
#   group_by(site, exp_dates) %>%
#   mutate(time_since_end = case_when(
#     exp_dates == "during" ~ -(max(time_yrs) + 0.25 - time_yrs),
#     exp_dates == "after" ~ time_yrs - min(time_yrs)
#   )) %>%
#   ungroup() %>%
#   left_join(., enframe(sites_full), by = c("site" = "name")) %>%
#   rename(site_full = value) %>%
#   mutate(site_full = fct_relevel(site_full, "Arroyo Quemado (AQUE)", "Naples (NAPL)", "Mohawk (MOHK)", "Carpinteria (CARP)")) %>%
#   # create a new column for "kelp year"
#   mutate(quarter = fct_relevel(quarter, "Q2", "Q3", "Q4", "Q1")) %>%
#   mutate(kelp_year = case_when(
#     quarter %in% c("Q2", "Q3", "Q4") ~ year,
#     quarter == "Q1" ~ year - 1
#   )) %>%
#   mutate(kelp_year = paste("kelp_", kelp_year, "-", kelp_year + 1, sep = "")) %>%
#   # create a column for the points to compare for "2 year interval"
#   mutate(comp_2yrs = case_when(
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012") ~ "start",
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013") ~ "start",
#     site == "ivee" & kelp_year %in% c("kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     kelp_year %in% c("kelp_2020-2021", "kelp_2021-2022") ~ "after"
#   )) %>%
#   mutate(comp_2yrs = fct_relevel(comp_2yrs, "start", "during", "after")) %>%
#   # create a column for the points to compare for "3 year interval"
#   mutate(comp_3yrs = case_when(
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2010-2011", "kelp_2011-2012", "kelp_2012-2013") ~ "start",
#     site %in% c("aque", "napl", "mohk", "carp") & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     site == "ivee" & kelp_year %in% c("kelp_2011-2012", "kelp_2012-2013", "kelp_2013-2014") ~ "start",
#     site == "ivee" & kelp_year %in% c("kelp_2014-2015", "kelp_2015-2016", "kelp_2016-2017") ~ "during",
#     kelp_year %in% c("kelp_2019-2020", "kelp_2020-2021", "kelp_2021-2022") ~ "after"
#   )) %>%
#   mutate(comp_3yrs = fct_relevel(comp_3yrs, "start", "during", "after")) %>%
#   # create a new sample ID that is site, year, quarter
#   unite("sample_ID", site, year, quarter, remove = FALSE) %>%
#   # add a new column for the delta continual at the last experimental period
#   mutate(delta_continual_end = case_when(
#     site == "aque" ~ -383.9406750,
#     site == "napl" ~ -1009.020417,
#     site == "mohk" ~ -1625.65538,
#     site == "carp" ~ -709.42300000
#   )) %>%
#   left_join(., site_quality, by = "site")
# # add a new column for the mean delta continual for the last 2 years of the experiment
# group_by(site, comp_3yrs) %>%
# mutate(delta_continual_end_3yrs = case_when(
#   comp_3yrs == "during" ~ mean(delta_continual)
# ))
# mutate(delta_continual_end_2yrs_test = case_when(
#   site == "aque" ~ -257.79005,
#   site == "napl" ~ -360.13564, 
#   site == "mohk" ~ -722.17182,
#   site == "carp" ~ -73.77815
# ))


# 4. comparisons ----------------------------------------------------------

### Before/after comparison

# This creates plots where the mean kelp frond density is plotted for the "before" and "after" time points ("before" = during the experiment, and "after" = after the experiment ended). 
# 
# #### Frond density
# ```{r}
# comparison_frond_fxn <- function(site_choice) {
#   plot_title <- if(site_choice == "aque") {
#     aque_full
#   } else if(site_choice == "napl") {
#     napl_full
#   } else if(site_choice == "carp") {
#     carp_full
#   } else if(site_choice == "ivee") {
#     ivee_full
#   } else if(site_choice == "mohk") {
#     mohk_full
#   }
#   
#   total_fronds %>% 
#     filter(site == site_choice) %>% 
#     ggplot(aes(x = exp_dates, y = total_fronds, col = treatment, shape = treatment)) +
#     # stat summary for point and SE for mean frond density during and after with line
#     stat_summary(aes(group = treatment), fun = mean, geom = "line", size = 1, alpha = 0.8) +
#     stat_summary(fun = mean, geom = "point", size = 4) +
#     stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.8, width = 0.2) +
#     scale_color_manual(values = color_palette) +
#     scale_shape_manual(values = shape_palette) +
#     theme_bw() +
#     theme(plot.title = element_text(size = 20)) +
#     labs(x = " ",
#          y = "Mean frond density",
#          title = plot_title) 
# }
# 
# aque_comparison_frond <- comparison_frond_fxn("aque")
# mohk_comparison_frond <- comparison_frond_fxn("mohk")
# carp_comparison_frond <- comparison_frond_fxn("carp")
# ivee_comparison_frond <- total_fronds %>% 
#   filter(site == "ivee") %>% 
#   ggplot(aes(x = exp_dates, y = total_fronds, col = treatment, shape = treatment)) +
#   stat_summary(aes(group = treatment), fun = mean, geom = "line", size = 1, alpha = 0.8) +
#   stat_summary(fun = mean, geom = "point", size = 4) +
#   stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.8, width = 0.2) +
#   scale_color_manual(values = color_palette) +
#   scale_shape_manual(values = shape_palette) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 20)) +
#   labs(x = " ",
#        y = "Mean frond density",
#        title = ivee_full)
# napl_comparison_frond <- comparison_frond_fxn("napl")
# 
# aque_comparison_frond
# mohk_comparison_frond
# carp_comparison_frond
# ivee_comparison_frond
# napl_comparison_frond
# ```
# 
# #### Kelp biomass
# ```{r}
# comparison_biomass_fxn <- function(site_choice) {
#   plot_title <- pluck(sites_full, site_choice)
#   
#   kelp_biomass %>% 
#     filter(site == site_choice) %>% 
#     ggplot(aes(x = exp_dates, y = wm_gm2, col = treatment, shape = treatment)) +
#     stat_summary(aes(group = treatment), fun = mean, geom = "line", size = 1, alpha = 0.8, na.rm = TRUE) +
#     stat_summary(fun = mean, geom = "point", size = 5, na.rm = TRUE) +
#     
#     stat_summary(fun.data = mean_se, geom = "errorbar", size = 1, width = 0.2, na.rm = TRUE) +
#     scale_color_manual(values = color_palette) +
#     scale_shape_manual(values = shape_palette) +
#     theme_bw() +
#     theme(plot.title = element_text(size = 20),
#           axis.text = element_text(size = 16),
#           axis.title = element_text(size = 18),
#           legend.text = element_text(size = 14),
#           legend.title = element_text(size = 14)) +
#     labs(x = " ",
#          y = "Mean kelp biomass (wm g/m2)",
#          title = plot_title,
#          col = "Treatment", shape = "Treatment")
# }
# 
# aque_comparison_biomass <- comparison_biomass_fxn("aque")
# napl_comparison_biomass <- comparison_biomass_fxn("napl")
# ivee_comparison_biomass <- comparison_biomass_fxn("ivee")
# mohk_comparison_biomass <- comparison_biomass_fxn("mohk")
# carp_comparison_biomass <- comparison_biomass_fxn("carp")
# 
# aque_comparison_biomass
# napl_comparison_biomass
# ivee_comparison_biomass
# mohk_comparison_biomass
# carp_comparison_biomass
# ```
# 
# ### best model for kelp frond density
# ```{r}
# summary(aov(data = total_fronds,
#             formula = total_fronds ~ treatment + site))
# 
# TukeyHSD(aov(data = total_fronds,
#              formula = total_fronds ~ site + treatment + site*treatment))
# 
# model.dredge <- dredge(
#   lm(total_fronds ~ site + treatment + site*treatment, data = total_fronds, na.action = na.pass)
# )
# 
# model.sel(model.dredge) %>% 
#   gt() %>% 
#   tab_header(title = "all model subsets",
#              subtitle = "total_fronds ~ site + treatment + site*treatment") 
# # %>% 
# #   gtsave("tab_1.png", expand = 10,
# #          path = here::here("tables", "kelp_frond_model.png"))
# ```
# 
# ### best model for kelp biomass
# 
# ```{r}
# summary(aov(data = kelp_biomass,
#             formula = wm_gm2 ~ site + treatment + site*treatment))
# 
# TukeyHSD(aov(data = kelp_biomass,
#              formula = wm_gm2 ~ site + treatment + site*treatment))
# 
# model.dredge <- dredge(
#   lm(wm_gm2 ~ site + treatment + site*treatment, data = kelp_biomass, na.action = na.pass)
# )
# 
# model.sel(model.dredge) %>% 
#   gt() %>% 
#   tab_header(title = "all model subsets",
#              subtitle = "kelp_biomass ~ site + treatment + site*treatment") 
# ```
# 
# ### Total kelp frond timeseries
# 
# This creates a plot where kelp frond density timeseries are plotted for each LTE site with control in open circles, annual in green circles, and continual in blue triangles. There is a `ggsave()` command at the end - **remember to change the date in the file name before saving**.
# 
# ```{r}
# total_frond_timeseries_fxn <- function(site_choice) {
#   after_date <- if(site_choice == "aque") {
#     aque_after_date_continual
#   } else if(site_choice == "napl") {
#     napl_after_date_continual
#   } else if(site_choice == "carp") {
#     carp_after_date_continual
#   } else if(site_choice == "ivee") {
#     ivee_after_date_continual
#   } else if(site_choice == "mohk") {
#     mohk_after_date_continual
#   }
#   
#   plot_title <- if(site_choice == "aque") {
#     aque_full
#   } else if(site_choice == "napl") {
#     napl_full
#   } else if(site_choice == "carp") {
#     carp_full
#   } else if(site_choice == "ivee") {
#     ivee_full
#   } else if(site_choice == "mohk") {
#     mohk_full
#   }
#   
#   rect_ymax <- total_fronds %>% filter(site == site_choice) %>% pull(total_fronds) %>% max()
#   
#   df <- total_fronds %>% filter(site == site_choice)
#   
#   ggplot(data = df,
#          aes(x = date, y = total_fronds, color = treatment, shape = treatment)) +
#     scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#     geom_hline(yintercept = 0, color = "grey") +
#     annotate("rect", xmin = after_date, xmax = as.Date("2022-02-02"),
#              ymin = -1, ymax = rect_ymax*1.05,
#              alpha = 0.2) +
#     geom_point(size = 3) +
#     geom_line(size = 1) +
#     scale_color_manual(values = color_palette) +
#     scale_shape_manual(values = shape_palette)  +
#     theme_bw() + 
#     labs(title = plot_title,
#          x = " ",
#          y = " ")
# }
# 
# aque_total_frond_timeseries <- total_frond_timeseries_fxn("aque")
# mohk_total_frond_timeseries <- total_frond_timeseries_fxn("mohk")
# carp_total_frond_timeseries <- total_frond_timeseries_fxn("carp")
# ivee_total_frond_timeseries <- ggplot(data = total_fronds %>% filter(site == "ivee"),
#                                       aes(x = date, y = total_fronds, color = treatment, shape = treatment)) +
#   scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#   geom_hline(yintercept = 0, color = "grey") +
#   annotate("rect", xmin = ivee_after_date, xmax = as.Date("2022-02-02"),
#            ymin = -1, ymax = (total_fronds %>% filter(site == "ivee") %>% pull(total_fronds) %>% max())*1.05,
#            alpha = 0.2) +
#   geom_point(size = 3) +
#   geom_line(size = 1) +
#   scale_color_manual(values = c("control" = control_col, "annual" = annual_col)) +
#   scale_shape_manual(values = c("control" = control_shape, "annual" = annual_shape))  +
#   theme_bw() + 
#   labs(title = ivee_full,
#        x = " ",
#        y = " ")
# napl_total_frond_timeseries <- total_frond_timeseries_fxn("napl")
# 
# aque_total_frond_timeseries
# mohk_total_frond_timeseries
# ivee_total_frond_timeseries
# napl_total_frond_timeseries
# carp_total_frond_timeseries
# 
# # y-axis label for patchwork side eye emoji
# yaxislabel <- ggplot(data.frame(l = "kelp frond density", x = 1, y = 1)) +
#   geom_text(aes(x, y, label = l), angle = 90) + 
#   theme_void() +
#   coord_cartesian(clip = "off")
# 
# plots_together <- yaxislabel + (aque_total_frond_timeseries / mohk_total_frond_timeseries / napl_total_frond_timeseries / ivee_total_frond_timeseries / carp_total_frond_timeseries) +
#   plot_layout(guides = "collect", widths = c(1, 25)) &
#   theme(legend.position = "bottom")
# 
# plots_together
# 
# # ggsave(here::here("figures", "kelp_total_frond_ts-20211012.jpg"), plots_together, height = 9, width = 5, dpi = 200)
# ```
# 
# #### Total frond timeseries for report
# 
# ```{r}
# aque_total_frond_timeseries_labs <- aque_total_frond_timeseries + 
#   labs(x = "Sampling date", y = "Total kelp fronds",
#        title = "Total kelp fronds") +
#   theme(legend.position = "none")
# mohk_total_frond_timeseries_labs <- mohk_total_frond_timeseries + 
#   labs(x = "Sampling date", y = "Total kelp fronds",
#        title = "Total kelp fronds") +
#   theme(legend.position = "none")
# ivee_total_frond_timeseries_labs <- ivee_total_frond_timeseries + 
#   labs(x = "Sampling date", y = "Total kelp fronds",
#        title = "Total kelp fronds") +
#   theme(legend.position = "none")
# napl_total_frond_timeseries_labs <- napl_total_frond_timeseries + 
#   labs(x = "Sampling date", y = "Total kelp fronds",
#        title = "Total kelp fronds") +
#   theme(legend.position = "none")
# carp_total_frond_timeseries_labs <- carp_total_frond_timeseries + 
#   labs(x = "Sampling date", y = "Total kelp fronds",
#        title = "Total kelp fronds") +
#   theme(legend.position = "none")
# ```


# 5. plotting, misc. timeseries -------------------------------------------

## Plotting

### Total kelp biomass timeseries

# ```{r}
# kelp_biomass_timeseries_fxn <- function(site_choice) {
#   after_date <- if(site_choice == "aque") {
#     aque_after_date_continual
#   } else if(site_choice == "napl") {
#     napl_after_date_continual
#   } else if(site_choice == "carp") {
#     carp_after_date_continual
#   } else if(site_choice == "ivee") {
#     ivee_after_date_continual
#   } else if(site_choice == "mohk") {
#     mohk_after_date_continual
#   }
#   
#   plot_title <- if(site_choice == "aque") {
#     aque_full
#   } else if(site_choice == "napl") {
#     napl_full
#   } else if(site_choice == "carp") {
#     carp_full
#   } else if(site_choice == "ivee") {
#     ivee_full
#   } else if(site_choice == "mohk") {
#     mohk_full
#   }
#   
#   rect_ymax <- kelp_biomass %>% filter(site == site_choice) %>% pull(wm_gm2) %>% na.omit %>%  max()
#   
#   df <- kelp_biomass %>% filter(site == site_choice) %>% drop_na(wm_gm2)
#   
#   ggplot(data = df,
#          aes(x = date, y = wm_gm2, color = treatment, shape = treatment)) +
#     scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#     geom_hline(yintercept = 0, color = "grey") +
#     annotate("rect", xmin = after_date, xmax = as.Date("2022-02-02"),
#              ymin = -1, ymax = rect_ymax*1.05,
#              alpha = 0.2) +
#     geom_point(size = 3) +
#     geom_line(size = 1) +
#     scale_color_manual(values = color_palette) +
#     scale_shape_manual(values = shape_palette)  +
#     theme_bw() + 
#     theme(legend.position = "none") +
#     labs(title = plot_title,
#          x = "Sampling date",
#          y = "Kelp biomass (wm g/m2)")
# }
# 
# aque_kelp_biomass_timeseries <- kelp_biomass_timeseries_fxn("aque")
# napl_kelp_biomass_timeseries <- kelp_biomass_timeseries_fxn("napl")
# ivee_kelp_biomass_timeseries <- ggplot(data = kelp_biomass %>% 
#                                          filter(site == "ivee"),
#                                        aes(x = date, y = wm_gm2, color = treatment, shape = treatment)) +
#   scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#   geom_hline(yintercept = 0, color = "grey") +
#   annotate("rect", xmin = ivee_after_date_annual, xmax = as.Date("2022-02-02"),
#            ymin = -1, ymax = (kelp_biomass %>% filter(site == "ivee") %>% pull(wm_gm2) %>% max())*1.05,
#            alpha = 0.2) +
#   geom_point(size = 3) +
#   geom_line(size = 1) +
#   scale_color_manual(values = c("control" = "grey", "annual" = annual_col, "continual" = continual_col)) +
#   scale_shape_manual(values = c("control" = control_shape, "annual" = annual_shape, "continual" = continual_shape))  +
#   theme_bw() + 
#   theme(legend.position = "none") +
#   labs(title = ivee_full,
#        x = "Sampling date",
#        y = "Total kelp biomass (wm g/m2)")
# mohk_kelp_biomass_timeseries <- kelp_biomass_timeseries_fxn("mohk")
# carp_kelp_biomass_timeseries <- kelp_biomass_timeseries_fxn("carp")
# 
# aque_kelp_biomass_timeseries
# napl_kelp_biomass_timeseries
# ivee_kelp_biomass_timeseries
# mohk_kelp_biomass_timeseries
# carp_kelp_biomass_timeseries
# ```
# 
# ### Annual kelp biomass timeseries
# 
# ```{r}
# annual_biomass_timeseries <- kelp_biomass %>% 
#   filter(treatment == "annual") %>% 
#   group_by(site) %>% 
#   mutate(site_label = case_when(
#     date == last(date) ~ sites_full[as.character(site)],
#     TRUE ~ ""
#   )) %>% 
#   ggplot(aes(x = date, y = wm_gm2)) +
#   annotate("rect", xmin = as.Date("2017-01-01"), xmax = as.Date("2022-02-02"),
#            ymin = -1, ymax = 10000,
#            alpha = 0.2) +
#   geom_point(aes(col = site), size = 3) +
#   geom_line(aes(col = site), size = 1) +
#   # geom_text(aes(label = site_label, col = site), hjust = 0.1, size = 4, vjust = 1, 
#   #           label.padding = unit(0.25, "lines")) +
#   scale_color_manual(values = site_palette) +
#   labs(x = "Sampling date", y = "Kelp biomass (wm g/m2)",
#        title = "Annual only") +
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         plot.title = element_text(size = 16)) 
# 
# # geom_text(data = . %>% filter(date == max(date)),
# #           aes(color = Province_State, x = as.Date(Inf),
# #               y = deaths_roll7_100k),
# #           hjust = 0, size = 4, vjust = 0.7,
# #           label = c("Arizona\n", "North Carolina"))
# 
# 
# annual_biomass_timeseries
# 
# # ggsave(here::here("figures", "annual_biomass_timeseries.jpg"), annual_biomass_timeseries, width = 14, height = 8, dpi = 300)
# ```
# 
# 
# ### Delta kelp frond timeseries
# 
# This creates a plot where delta kelp frond density (difference between control and removal density) timeseries are plotted for each LTE site with annual in green circles and continual in blue circles. This was supposed to be a recreation of the plots that Adrian sent a long time ago. There is a `ggsave()` command at the end - **remember to change the date in the file name before saving**.
# 
# ```{r}
# delta_frond_timeseries_fxn <- function(site_choice) {
#   after_date <- if(site_choice == "aque") {
#     aque_after_date
#   } else if(site_choice == "napl") {
#     napl_after_date
#   } else if(site_choice == "carp") {
#     carp_after_date
#   } else if(site_choice == "ivee") {
#     ivee_after_date
#   } else if(site_choice == "mohk") {
#     mohk_after_date
#   }
#   
#   delta_continual_max <- delta_fronds %>% filter(site == site_choice) %>% pull(delta_continual) %>% max(na.rm = TRUE)
#   delta_continual_min <- delta_fronds %>% filter(site == site_choice) %>% pull(delta_continual) %>% min(na.rm = TRUE)
#   delta_annual_max <- delta_fronds %>% filter(site == site_choice) %>% pull(delta_annual) %>% max(na.rm = TRUE)  
#   delta_annual_min <- delta_fronds %>% filter(site == site_choice) %>% pull(delta_annual) %>% min(na.rm = TRUE)
#   
#   rect_ymax <- ifelse(delta_continual_max > delta_annual_max, delta_continual_max, delta_annual_max)
#   rect_ymin <- ifelse(delta_continual_min < delta_annual_max, delta_continual_min, delta_annual_max)
#   
#   ggplot() +
#     scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#     annotate("rect", xmin = after_date, xmax = as.Date("2022-02-02"),
#              ymin = rect_ymin*1.1, ymax = rect_ymax*1.1,
#              alpha = 0.2) +
#     annotate("text", x = as.Date("2009-04-01"), y = rect_ymax*1.1, label = "removal > control", size = 3) +
#     annotate("text", x = as.Date("2009-04-01"), y = rect_ymin*1.1, label = "control > removal", size = 3) +
#     geom_point(data = delta_fronds %>% filter(site == site_choice & delta_annual != "NA"), 
#                aes(x = date, y = delta_annual, color = "annual"), 
#                size = 3, shape = annual_shape) +
#     geom_line(data = delta_fronds %>% filter(site == site_choice & delta_annual != "NA"),
#               aes(x = date, y = delta_annual, color = "annual"), size = 1) +
#     geom_point(data = delta_fronds %>% filter(site == site_choice & delta_continual != "NA"), 
#                aes(x = date, y = delta_continual, color = "continual"), 
#                size = 3, shape = continual_shape) +
#     geom_line(data = delta_fronds %>% filter(site == site_choice & delta_continual != "NA"),
#               aes(x = date, y = delta_continual, color = "continual"), size = 1) +
#     scale_color_manual(name = NULL,
#                        breaks = c("annual", "continual"),
#                        values = c("annual" = annual_col, "continual" = continual_col)) +
#     theme_bw() + 
#     labs(title = site_choice,
#          x = " ",
#          y = " ")
# }
# 
# aque_delta_kelp_timeseries <- delta_frond_timeseries_fxn("aque")
# mohk_delta_kelp_timeseries <- delta_frond_timeseries_fxn("mohk")
# carp_delta_kelp_timeseries <- delta_frond_timeseries_fxn("carp")
# ivee_delta_kelp_timeseries <- ggplot() +
#   scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#   annotate("rect", xmin = ivee_after_date, xmax = as.Date("2022-02-02"),
#            ymin = delta_fronds %>% filter(site == "ivee") %>% pull(delta_annual) %>% min()*1.1, 
#            ymax = delta_fronds %>% filter(site == "ivee") %>% pull(delta_annual) %>% max()*1.1,
#            alpha = 0.2) +
#   annotate("text", x = as.Date("2009-04-01"), y = delta_fronds %>% filter(site == "ivee") %>% pull(delta_annual) %>% max()*1.1, label = "removal > control", size = 3) +
#   annotate("text", x = as.Date("2009-04-01"), y = delta_fronds %>% filter(site == "ivee") %>% pull(delta_annual) %>% min()*1.1, label = "control > removal", size = 3) +
#   geom_point(data = delta_fronds %>% filter(site == "ivee" & delta_annual != "NA"), 
#              aes(x = date, y = delta_annual), color = annual_col, size = 3) +
#   geom_line(data = delta_fronds %>% filter(site == "ivee" & delta_annual != "NA"),
#             aes(x = date, y = delta_annual), color = annual_col, size = 1) +
#   theme_bw() + 
#   labs(title = "ivee",
#        x = "Sampling date",
#        y = " ")
# napl_delta_kelp_timeseries <- delta_frond_timeseries_fxn("napl")
# 
# aque_delta_kelp_timeseries
# mohk_delta_kelp_timeseries
# ivee_delta_kelp_timeseries
# napl_delta_kelp_timeseries
# carp_delta_kelp_timeseries
# 
# # y-axis label for patchwork side eye emoji
# yaxislabel <- ggplot(data.frame(l = "kelp fronds (removal - control)", x = 1, y = 1)) +
#   geom_text(aes(x, y, label = l), angle = 90) + 
#   theme_void() +
#   coord_cartesian(clip = "off")
# 
# plots_together <- yaxislabel + (aque_delta_kelp_timeseries / mohk_delta_kelp_timeseries / napl_delta_kelp_timeseries / ivee_delta_kelp_timeseries / carp_delta_kelp_timeseries) +
#   plot_layout(guides = "collect", widths = c(1, 25)) &
#   theme(legend.position = "bottom")
# 
# plots_together
# 
# # ggsave(here::here("figures", "kelp_frond_ts-20220201.jpg"), plots_together, height = 9, width = 5, dpi = 200)
# ```
# 
# #### Delta frond timeseries for report
# 
# ```{r}
# aque_delta_kelp_timeseries_labs <- aque_delta_kelp_timeseries +
#   labs(x = "Sampling date", y = "\u0394 fronds \n (removal - control)",
#        title = "\u0394 fronds (removal - control)")  +
#   theme(legend.position = "none")
# mohk_delta_kelp_timeseries_labs <- mohk_delta_kelp_timeseries +
#   labs(x = "Sampling date", y = "\u0394 fronds \n (removal - control)",
#        title = "\u0394 fronds (removal - control)")  +
#   theme(legend.position = "none")
# ivee_delta_kelp_timeseries_labs <- ivee_delta_kelp_timeseries +
#   labs(x = "Sampling date", y = "\u0394 fronds \n (removal - control)",
#        title = "\u0394 fronds (removal - control)")  +
#   theme(legend.position = "none")
# napl_delta_kelp_timeseries_labs <- napl_delta_kelp_timeseries +
#   labs(x = "Sampling date", y = "\u0394 fronds \n (removal - control)",
#        title = "\u0394 fronds (removal - control)")  +
#   theme(legend.position = "none")
# carp_delta_kelp_timeseries_labs <- carp_delta_kelp_timeseries +
#   labs(x = "Sampling date", y = "\u0394 fronds \n (removal - control)",
#        title = "\u0394 fronds (removal - control)") +
#   theme(legend.position = "none")
# ```
# 
# ### Delta biomass timeseries
# 
# ```{r}
# delta_biomass_timeseries_fxn <- function(site_choice) {
#   after_date <- if(site_choice == "aque") {
#     aque_after_date
#   } else if(site_choice == "napl") {
#     napl_after_date
#   } else if(site_choice == "carp") {
#     carp_after_date
#   } else if(site_choice == "ivee") {
#     ivee_after_date
#   } else if(site_choice == "mohk") {
#     mohk_after_date
#   }
#   
#   plot_title <- pluck(sites_full, site_choice)
#   
#   delta_continual_max <- delta_biomass %>% filter(site == site_choice) %>% pull(delta_continual) %>% max(na.rm = TRUE)
#   delta_continual_min <- delta_biomass %>% filter(site == site_choice) %>% pull(delta_continual) %>% min(na.rm = TRUE)
#   delta_annual_max <- delta_biomass %>% filter(site == site_choice) %>% pull(delta_annual) %>% max(na.rm = TRUE)  
#   delta_annual_min <- delta_biomass %>% filter(site == site_choice) %>% pull(delta_annual) %>% min(na.rm = TRUE)
#   
#   rect_ymax <- ifelse(delta_continual_max > delta_annual_max, delta_continual_max, delta_annual_max)
#   rect_ymin <- ifelse(delta_continual_min < delta_annual_max, delta_continual_min, delta_annual_max)
#   
#   ggplot() +
#     scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#     annotate("rect", xmin = after_date, xmax = as.Date("2022-02-02"),
#              ymin = rect_ymin*1.1, ymax = rect_ymax*1.1,
#              alpha = 0.2) +
#     annotate("text", x = as.Date("2009-04-01"), y = rect_ymax*1.1, label = "removal > control", size = 3) +
#     annotate("text", x = as.Date("2009-04-01"), y = rect_ymin*1.1, label = "control > removal", size = 3) +
#     geom_point(data = delta_biomass %>% filter(site == site_choice & delta_annual != "NA"), 
#                aes(x = date, y = delta_annual, color = "annual"), 
#                size = 3, shape = annual_shape) +
#     geom_line(data = delta_biomass %>% filter(site == site_choice & delta_annual != "NA"),
#               aes(x = date, y = delta_annual, color = "annual"), size = 1) +
#     geom_point(data = delta_biomass %>% filter(site == site_choice & delta_continual != "NA"), 
#                aes(x = date, y = delta_continual, color = "continual"), 
#                size = 3, shape = continual_shape) +
#     geom_line(data = delta_biomass %>% filter(site == site_choice & delta_continual != "NA"),
#               aes(x = date, y = delta_continual, color = "continual"), size = 1) +
#     scale_color_manual(name = NULL,
#                        breaks = c("annual", "continual"),
#                        values = c("annual" = annual_col, "continual" = continual_col)) +
#     theme_bw() + 
#     theme(legend.position = "none") +
#     labs(x = "Sampling date",
#          y = "\u0394 kelp biomass \n (wm g/m2, removal - control)",
#          title = plot_title)
# }
# 
# aque_delta_biomass_timeseries <- delta_biomass_timeseries_fxn("aque")
# napl_delta_biomass_timeseries <- delta_biomass_timeseries_fxn("napl")
# ivee_delta_biomass_timeseries <- ggplot() +
#   scale_x_date(limits = c(napl_start_date, as.Date("2022-02-02"))) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#   annotate("rect", xmin = ivee_after_date, xmax = as.Date("2022-02-02"),
#            ymin = delta_biomass %>% filter(site == "ivee") %>% pull(delta_annual) %>% min()*1.1, 
#            ymax = delta_biomass %>% filter(site == "ivee") %>% pull(delta_annual) %>% max()*1.1,
#            alpha = 0.2) +
#   annotate("text", x = as.Date("2009-04-01"), y = delta_biomass %>% filter(site == "ivee") %>% pull(delta_annual) %>% max()*1.1, label = "removal > control", size = 3) +
#   annotate("text", x = as.Date("2009-04-01"), y = delta_biomass %>% filter(site == "ivee") %>% pull(delta_annual) %>% min()*1.1, label = "control > removal", size = 3) +
#   geom_point(data = delta_biomass %>% filter(site == "ivee" & delta_annual != "NA"), 
#              aes(x = date, y = delta_annual, color = "annual"), 
#              size = 3, shape = annual_shape) +
#   geom_line(data = delta_biomass %>% filter(site == "ivee" & delta_annual != "NA"),
#             aes(x = date, y = delta_annual, color = "annual"), size = 1) +
#   scale_color_manual(name = NULL,
#                      breaks = c("annual"),
#                      values = c("annual" = annual_col)) +
#   theme_bw() + 
#   theme(legend.position = "none") +
#   labs(x = "Sampling date",
#        y = "\u0394 kelp biomass \n (wm g/m2, removal - control)",
#        title = "ivee")
# mohk_delta_biomass_timeseries <- delta_biomass_timeseries_fxn("mohk")
# carp_delta_biomass_timeseries <- delta_biomass_timeseries_fxn("carp")
# 
# 
# aque_delta_biomass_timeseries
# napl_delta_biomass_timeseries
# ivee_delta_biomass_timeseries
# mohk_delta_biomass_timeseries
# carp_delta_biomass_timeseries
# ```
# 
# ### 2022-02-04 report figures
# 
# ```{r}
# aque_panels <- (aque_comparison_frond / aque_total_frond_timeseries_labs / aque_delta_kelp_timeseries_labs) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 16))
# 
# napl_panels <- (napl_comparison_frond / napl_total_frond_timeseries_labs / napl_delta_kelp_timeseries_labs) +
#   plot_layout(guides = "collect") 
# 
# ivee_panels <- (ivee_comparison_frond / ivee_total_frond_timeseries_labs / ivee_delta_kelp_timeseries_labs) +
#   plot_layout(guides = "collect") 
# 
# carp_panels <- (carp_comparison_frond / carp_total_frond_timeseries_labs / carp_delta_kelp_timeseries_labs) +
#   plot_layout(guides = "collect") 
# 
# mohk_panels <- (mohk_comparison_frond / mohk_total_frond_timeseries_labs / mohk_delta_kelp_timeseries_labs) +
#   plot_layout(guides = "collect") 
# 
# # ggsave(here::here("figures", "panels_aque_2022-02-03.jpg"), aque_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_napl_2022-02-03.jpg"), napl_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_ivee_2022-02-03.jpg"), ivee_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_carp_2022-02-03.jpg"), carp_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_mohk_2022-02-03.jpg"), mohk_panels, height = 10, width = 9, dpi = 300)
# ```
# 
# ### 2022-02-11 report figures
# 
# ```{r}
# aque_panels <- (aque_comparison_biomass / aque_kelp_biomass_timeseries / aque_delta_biomass_timeseries) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 16))
# 
# napl_panels <- (napl_comparison_biomass / napl_kelp_biomass_timeseries / napl_delta_biomass_timeseries) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 16))
# 
# ivee_panels <- (ivee_comparison_biomass / ivee_kelp_biomass_timeseries / ivee_delta_biomass_timeseries) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 16))
# 
# carp_panels <- (carp_comparison_biomass / carp_kelp_biomass_timeseries / carp_delta_biomass_timeseries) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 16))
# 
# mohk_panels <- (mohk_comparison_biomass / mohk_kelp_biomass_timeseries / mohk_delta_biomass_timeseries) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 16))
# 
# # ggsave(here::here("figures", "panels_biomass_aque_2022-02-08.jpg"), aque_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_biomass_napl_2022-02-08.jpg"), napl_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_biomass_ivee_2022-02-08.jpg"), ivee_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_biomass_carp_2022-02-08.jpg"), carp_panels, height = 10, width = 9, dpi = 300)
# # 
# # ggsave(here::here("figures", "panels_biomass_mohk_2022-02-08.jpg"), mohk_panels, height = 10, width = 9, dpi = 300)
# ```
# 

# 6. biomass by season figures ---------------------------------------------

# 
# ## Biomass by season figure
# 
# x-axis: time from winter to fall  
# y-axis: kelp biomass  
# lines: colored by sampling date and year  
# 
# Going to try with aque first:
#   ```{r}
# aque_test <- kelp_biomass %>% 
#   filter(site == "aque" & year %in% c(2013, 2014, 2015, 2016, 2017, 2018, 2019)) %>% 
#   unite("color", treatment, year, sep = "_", remove = FALSE) %>% 
#   # make a new column for the year label lol i am a clown
#   mutate(year_label = case_when(
#     season == "winter" ~ as.character(year),
#     TRUE ~ ""
#   )) %>% 
#   ggplot(aes(x = season, y = total_wm_gm2, group = color)) +
#   geom_point(aes(color = color)) +
#   geom_line(aes(color = color)) +
#   geom_text_repel(aes(label = year_label, color = color)) +
#   theme(legend.position = "none") +
#   facet_wrap(~ treatment + exp_dates, nrow = 3, ncol = 2)
# 
# 
# aque_test
# ```
# 
# Ok great. Going to do the rest now:
#   
#   ```{r}
# biomass_by_season <- function(site_choice) {
#   plot_title <- pluck(sites_full, site_choice)
#   
#   kelp_biomass %>% 
#     # filter(site == site_choice) %>% 
#     unite("color", treatment, year, sep = "_", remove = FALSE) %>% 
#     # make a new column for the year label
#     mutate(year_label = case_when(
#       season == "winter" ~ as.character(year),
#       TRUE ~ ""
#     )) %>%
#     mutate(season = fct_relevel(season, c("winter", "spring", "summer", "fall"))) %>% 
#     ggplot(aes(x = season, y = wm_gm2, group = color)) +
#     # stat summary for point and SE for mean biomass during and after with line
#     stat_summary(aes(group = season, color = site), fun = mean, geom = "line", size = 1, alpha = 0.8) +
#     stat_summary(aes(group = season, color = site), fun = mean, geom = "point", size = 4) +
#     stat_summary(aes(group = season, color = site), fun.data = mean_se, geom = "errorbar", alpha = 0.8, width = 0.2) +
#     # geom_point(aes(color = year)) +
#     #   geom_line(aes(color = year)) +
#     # geom_text_repel(aes(label = year_label, color = year),
#     #                 nudge_x = 1.5, nudge_y = 0.1,
#     #                 segment.curvature = -1e-20, point.padding = 0.2,
#     #                 force_pull = 10) +
#     theme_bw() +
#     theme(legend.position = "none",
#           plot.title = element_text(size = 20),
#           strip.background = element_rect(fill = "#FFFFFF")) +
#     facet_wrap(~ treatment + exp_dates, nrow = 3, ncol = 2) +
#     labs(x = "Season", y = "Total biomass (wet g/m2)",
#          title = plot_title)
# }
# 
# aque_biomass_season <- biomass_by_season("aque")
# napl_biomass_season <- biomass_by_season("napl")
# ivee_biomass_season <- biomass_by_season("ivee")
# mohk_biomass_season <- biomass_by_season("mohk")
# carp_biomass_season <- biomass_by_season("carp")
# 
# delta_annual %>% 
#   filter(site == "carp") %>% 
#   mutate(season = case_when(
#     quarter == "Q2" ~ "spring",
#     quarter == "Q3" ~ "summer",
#     quarter == "Q4" ~ "fall",
#     quarter == "Q1" ~ "winter"
#   )) %>% 
#   mutate(season = fct_relevel(season, "spring", "summer", "fall", "winter")) %>% 
#   filter(exp_dates == "during") %>% 
#   mutate(quarter = fct_relevel(quarter, "Q2", "Q3", "Q4", "Q1")) %>% 
#   mutate(kelp_year = as.factor(kelp_year)) %>% 
#   ggplot(aes(x = season, y = annual)) +
#   stat_summary(aes(color = kelp_year, group = kelp_year), fun.data = mean_se, geom = "errorbar", alpha = 0.5, width = 0.2) +
#   # geom_point(aes(shape = site, color = kelp_year), size = 2) +
#   stat_summary(aes(color = kelp_year, group = kelp_year), fun = mean, geom = "line", size = 2) +
#   scale_color_manual(values = c("#BBD1FA", "#A2B9E4", "#89A1CE", "#7089B8",
#                                 "#5771A2", "#3E598C", "#254176", "#0D2A61", "darkblue", "navy", "black")) +
#   theme_bw() +
#   labs(x = "Season",
#        y = "Kelp biomass",
#        caption = "carp") +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16),
#         panel.grid = element_blank())
# # ggsave(here::here("figures", paste("biomass-by-season_carp_", todays_date, ".jpg", sep = "")),  
# #        height = 10, width = 8, dpi = 150)
# 
# # ggsave(here::here("figures", "biomass-by-season_aque_2022-02-08.jpg"), aque_biomass_season, height = 10, width = 8, dpi = 200)
# 
# # ggsave(here::here("figures", "biomass-by-season_napl_2022-02-08.jpg"), napl_biomass_season, height = 10, width = 8, dpi = 200)
# 
# # ggsave(here::here("figures", "biomass-by-season_ivee_2022-02-08.jpg"), ivee_biomass_season, height = 10, width = 8, dpi = 200)
# 
# # ggsave(here::here("figures", "biomass-by-season_mohk_2022-02-08.jpg"), mohk_biomass_season, height = 10, width = 8, dpi = 200)
# 
# # ggsave(here::here("figures", "biomass-by-season_carp_2022-02-08.jpg"), carp_biomass_season, height = 10, width = 8, dpi = 200)
# ```
# 
# ```{r}
# annual_removal_biomass_by_season <- kelp_biomass %>% 
#   filter(treatment == "annual") %>% 
#   unite("color", treatment, year, sep = "_", remove = FALSE) %>% 
#   # make a new column for the year label
#   mutate(year_label = case_when(
#     season == "winter" ~ as.character(year),
#     TRUE ~ ""
#   )) %>% 
#   ggplot(aes(x = season, y = wm_gm2)) +
#   # stat summary for point and SE for mean biomass during and after with line
#   stat_summary(aes(group = site, color = site), fun = mean, geom = "line", size = 1, alpha = 0.8) +
#   stat_summary(aes(group = site, color = site), fun = mean, geom = "point", size = 4) +
#   stat_summary(aes(group = site, color = site), fun.data = mean_se, geom = "errorbar", alpha = 0.8, width = 0.2) +
#   scale_color_manual(values = site_palette) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 40),
#         strip.background = element_rect(fill = "#FFFFFF"),
#         strip.text = element_text(size = 18)) +
#   facet_wrap(~ exp_dates) +
#   labs(x = "Season", y = "Mean biomass (wet g/m2)",
#        title = "Annual removal biomass") 
# 
# # ggsave(here::here("figures", "annual_biomass_timeseries.jpg"), annual_removal_biomass_by_season, width = 14, height = 8, dpi = 300)
# ```
# 
# ```{r}
# biomass_by_season_all <- kelp_biomass %>% 
#   unite("color", treatment, year, sep = "_", remove = FALSE) %>% 
#   # make a new column for the year label
#   mutate(year_label = case_when(
#     season == "winter" ~ as.character(year),
#     TRUE ~ ""
#   )) %>%
#   mutate(season = fct_relevel(season, c("winter", "spring", "summer", "fall"))) %>% 
#   ggplot(aes(x = season, y = wm_gm2, group = site)) +
#   # stat summary for point and SE for mean biomass during and after with line
#   stat_summary(aes(group = site, color = site), fun = mean, geom = "line", size = 1, alpha = 0.8) +
#   stat_summary(aes(group = site, color = site), fun = mean, geom = "point", size = 4) +
#   stat_summary(aes(group = site, color = site), fun.data = mean_se, geom = "errorbar", alpha = 0.8, width = 0.2) +
#   scale_color_manual(values = site_palette) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 20),
#         strip.background = element_rect(fill = "#FFFFFF"),
#         strip.text = element_text(size = 14),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 12)) +
#   facet_wrap(~ treatment + exp_dates, nrow = 3, ncol = 2) +
#   labs(x = "Season", y = "Mean biomass (wet g/m2)",
#        title = "Mean kelp biomass within season across sampling years")
# 
# # ggsave(here::here("figures", "annual_biomass_timeseries.jpg"), biomass_by_season_all, width = 11, height = 12, dpi = 300)
# ```
# 
# 
# ## Fall biomass figure
# 
# x-axis: time  
# y-axis: biomass  
# only including fall sampling dates
# 
# ```{r}
# kelp_biomass_fall <- kelp_biomass %>% 
#   filter(season == "fall") %>% 
#   ggplot(aes(x = date, y = wm_gm2)) +
#   geom_line(aes(color = treatment), size = 1) +
#   geom_point(aes(color = treatment, shape = treatment), size = 3) +
#   scale_color_manual(values = color_palette) +
#   scale_shape_manual(values = shape_palette) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         title = element_text(size = 20),
#         legend.text = element_text(size = 16)) +
#   facet_wrap(~site) +
#   labs(x = "Sampling date", y = "Total biomass (wet mass/gm2)",
#        title = "Fall kelp biomass")
# 
# kelp_biomass_summer <- kelp_biomass %>% 
#   filter(season == "summer") %>% 
#   ggplot(aes(x = date, y = wm_gm2)) +
#   geom_line(aes(color = treatment), size = 1) +
#   geom_point(aes(color = treatment, shape = treatment), size = 3) +
#   scale_color_manual(values = color_palette) +
#   scale_shape_manual(values = shape_palette) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         title = element_text(size = 20),
#         legend.text = element_text(size = 16)) +
#   facet_wrap(~site) +
#   labs(x = "Sampling date", y = "Total biomass (wet mass/gm2)",
#        title = "Summer kelp biomass")
# 
# kelp_biomass_fall
# kelp_biomass_summer
# 
# # ggsave(here::here("figures/timeseries", paste("kelp-biomass_summer_", todays_date, ".jpg", sep = "")), kelp_biomass_summer)
# 
# # ggsave(here::here("figures/timeseries", paste("kelp-biomass_fall_", todays_date, ".jpg", sep = "")), kelp_biomass_fall)
# ```
# 
# ## Kelp "recovery" by productivity
# 
# Plot with "productivity" on x-axis (CARP, AQUE, NAPL, MOHK, IVEE) and time to recovery (delta = 0) on y axis.
# 
# ```{r}
# df <- delta_biomass %>% 
#   filter(site == "mohk" & exp_dates == "after") %>% 
#   group_by(site) %>% 
#   mutate(time_to_zero = row_number())
# 
# ```
# 
# ## Removal kelp as a function of control
# 
# ### test
# 
# What I want is to have a plot where kelp biomass in the control is the x-axis, and kelp biomass in the removal is the y-axis, and there is a 1:1 line as the diagonal, and each point is connected by an arrow indicating time across sampling dates.
# 
# ```{r}
# delta_biomass %>% 
#   filter(site == "aque") %>% 
#   ggplot(aes(x = control, y = annual)) +
#   geom_text(aes(col = exp_dates, label = date)) +
#   geom_abline(slope = 1, intercept = 0) +
#   coord_cartesian()
# ```
# 
# ## Ratio of kelp in removal:control plots
# 
# Connell et al. 1997:  
#   - damage: 100 x (difference between abundance in census before with that just after disturbance)/(abundance in census just before disturbance)
# - short term rate of recovery: 100 x (change in abundance between first two censuses)/(damage during disturbance), which is then standardized to 2 year period by "linear interpolation of period between two post-storm censuses"  
# - extent of recovery: 100 x (highest abundance in period before next disturbance)/(damage during first disturbance)  
# 
# So, for kelp:
#   - impact of removal: 100 x (difference in biomass in control - removal)/(biomass in control)
# - rate of recovery: 100 x (change in abundance between end of experiment and sampling date)/(impact of removal)  
# - extent of recovery: 100 x (highest abundance in period after disturbance)/(impact of removal)
# 
# ```{r}
# df <- delta_biomass %>% 
#   filter(site == "aque") %>% 
#   mutate(delta_annual = control - annual,
#          delta_continual = control - continual, 
#          impact_annual = round((delta_annual/control), 4),
#          impact_continual = round((delta_continual/control), 4)) %>% 
#   filter(impact_annual > -Inf)
# 
# ggplot(df) +
#   geom_line(aes(x = date, y = impact_annual), col = "red") +
#   # geom_line(aes(x = date, y = impact_continual), col = "blue") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#   annotate("text", x = as.Date("2018-01-01"), y = 10, label = "more kelp in control", size = 3) +
#   annotate("text", x = as.Date("2018-01-01"), y = -10, label = "more kelp in removal", size = 3) +
#   theme_bw() +
#   labs(x = "Date", y = "Removal impact [(control-annual)/control]",
#        title = "Removal impact - AQUE")
# 
# ggplot(df, aes(x = control, y = delta_annual)) +
#   geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#   geom_smooth(method = "lm", aes(color = exp_dates)) +
#   theme_bw() +
#   labs(x = "Control kelp biomass", y = "Biomass difference (control - removal)",
#        title = "AQUE - delta biomass as a function of control kelp biomass") +
#   annotate("text", x = 6000, y = 200, label = "more kelp in control", size = 3) +
#   annotate("text", x = 6000, y = -200, label = "more kelp in removal", size = 3) 
# 
# ggplot(df %>% filter(impact_annual > -200), aes(x = control, y = impact_annual)) +
#   geom_point() +
#   # geom_point(data = df %>% filter(impact_annual < -500), aes(x = control, y = impact_annual), col = "red") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#   geom_smooth(method = "lm", aes(color = exp_dates)) +
#   theme_bw() +
#   labs(x = "Control kelp biomass", y = "Removal impact [(control-annual)/control]",
#        title = "AQUE - removal impact as a function of control kelp biomass") +
#   annotate("text", x = 6000, y = 50, label = "more kelp in control", size = 3) +
#   annotate("text", x = 6000, y = -50, label = "more kelp in removal", size = 3) 
# ```
# 
# Or, could try impact/rate/recovery with means...  
# 
# Just as a note, there probably needs to be two different data frames? Or maybe this can be done in vectors... regardless, maybe it could be a matrix that is 5 rows long, with columns for mean biomass annual during, mean biomass continual during... then deltas would be biomass
# 
# ```{r}
# 
# # - impact of removal: 100 x (difference in biomass in control - removal)/(biomass in control)
# # issues: when removal = 0 --> solution: change 0 to 0.000001? or some very small number?
# # - rate of recovery: 100 x (change in abundance between end of experiment and sampling date)/(impact of removal)  
# # - extent of recovery: 100 x (highest abundance in period after disturbance)/(impact of removal)
# 
# # Connell 1997: extent of recovery: ratio of increase to the amount of decline
# # e.g. if cover fell from 75% to 25%, then rose to 55%, then degree of recovery = 60% (rise 30%/fall 50%)
# 
# 
# recovery_df <- biomass %>% 
#   filter(sp_code == "MAPY") %>% 
#   select(-sp_code) %>% 
#   # # make a new column for after dates
#   after_dates_column() %>% 
#   # make a new column for during and after and set factor levels
#   exp_dates_column() %>% 
#   # select(site, treatment, date, exp_dates, wm_gm2) %>% 
#   group_by(exp_dates, site, treatment) %>% 
#   summarize(mean_biomass = mean(wm_gm2, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = c(treatment, exp_dates), values_from = mean_biomass) %>% 
#   # difference between control and removal during and after experiment ends
#   mutate(during_delta_annual = control_during - annual_during,
#          during_delta_continual = control_during - continual_during,
#          # "rise" is the increase after the experiment ends (difference between after and during)
#          both_rise_annual = annual_after - annual_during,
#          both_rise_continual = continual_after - continual_during,
#          # impact is the ratio of delta:control biomass
#          both_impact_annual = (during_delta_annual/control_during)*100,
#          both_impact_continual = (during_delta_continual/control_during)*100,
#          # rate of recovery is the ratio of rise:impact
#          both_rate_annual = both_rise_annual/both_impact_annual,
#          both_rate_continual = both_rise_continual/both_impact_continual,
#          # extent is the ratio of rise:delta
#          both_extent_annual = both_rise_annual/during_delta_annual,
#          both_extent_continual = both_rise_continual/during_delta_continual,
#          after_delta_annual = control_after - annual_after,
#          after_delta_continual = control_after - continual_after
#   ) %>% 
#   # new column with full site names
#   mutate(site_full = case_when(
#     site == "aque" ~ aque_full,
#     site == "napl" ~ napl_full,
#     site == "ivee" ~ ivee_full,
#     site == "mohk" ~ mohk_full,
#     site == "carp" ~ carp_full
#   )) %>% 
#   pivot_longer(cols = 8:19, names_to = "calc", values_to = "value") %>% 
#   separate(col = calc, into = c("exp_dates", "calc", "treatment")) %>% 
#   select(site_full, exp_dates, calc, treatment, value) %>% 
#   mutate(site_full = fct_relevel(site_full, aque_full, napl_full, ivee_full, mohk_full, carp_full)) 
# 
# aque_comparison_biomass
# napl_comparison_biomass
# ivee_comparison_biomass
# mohk_comparison_biomass
# carp_comparison_biomass
# 
# impact_plot <- ggplot(recovery_df %>% filter(calc == "impact"), aes(x = site_full, y = value)) +
#   geom_col(aes(fill = treatment), position = "dodge") +
#   scale_fill_manual(values = c(annual_col, continual_col)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
#   theme_bw() +
#   labs(x = "Site", y = "Mean kelp biomass impact [(\u0394 control - removal)/control]",
#        title = "Impact of removal during experiment") +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         title = element_text(size = 18),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 16))
# impact_plot
# # ggsave(here::here("figures/recovery-metrics", paste("impact_", todays_date, ".jpg", sep = "")), impact_plot)
# 
# extent_plot <- ggplot(recovery_df %>% filter(calc == "extent"), aes(x = site_full, y = value)) +
#   geom_col(aes(fill = treatment), position = "dodge") +
#   scale_fill_manual(values = c(annual_col, continual_col)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
#   theme_bw() +
#   labs(x = "Site", y = "Extent of recovery (mean kelp biomass) \n rise/\u0394",
#        title = "Extent of recovery following experiment end") +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         title = element_text(size = 18),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 16))
# extent_plot
# # ggsave(here::here("figures/recovery-metrics", paste("extent_", todays_date, ".jpg", sep = "")), extent_plot)
# 
# rate_plot <- ggplot(recovery_df %>% filter(calc == "rate"), aes(x = site_full, y = value)) +
#   geom_col(aes(fill = treatment), position = "dodge") +
#   scale_fill_manual(values = c(annual_col, continual_col)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
#   theme_bw() +
#   labs(x = "Site", y = "Rate of recovery (mean kelp biomass between 'during' and 'after') \n rise/impact",
#        title = "Rate of recovery following experiment end") +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         title = element_text(size = 18),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 16))
# rate_plot
# # ggsave(here::here("figures/recovery-metrics", paste("rate_", todays_date, ".jpg", sep = "")), rate_plot)
# ```


# 7. start-during-after comparisons ---------------------------------------

# summary(aov(delta_continual ~ comp_2yrs, data = delta_continual))
# TukeyHSD(aov(delta_continual ~ comp_2yrs, data = delta_continual))
# 
# summary(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "aque"))))
# TukeyHSD(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "aque"))))
# 
# summary(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "napl"))))
# TukeyHSD(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "napl"))))
# 
# summary(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "mohk"))))
# TukeyHSD(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "mohk"))))
# 
# summary(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "carp"))))
# TukeyHSD(aov(delta_continual ~ comp_2yrs, data = (delta_continual %>% filter(site == "carp"))))


# 
# 
# ```{r}
# ##### all sites #####
# # analysis of variance
# annual_comp_2yrs_aov <- aov(delta_annual ~ comp_2yrs, data = delta_annual)
# 
# # summary
# summary(annual_comp_2yrs_aov)
# 
# # Tukey and compact letter display
# annual_comp_2yrs_tukey <- TukeyHSD(aov(delta_annual ~ comp_2yrs, data = delta_annual))
# annual_comp_2yrs_cld <- multcompLetters4(annual_comp_2yrs_aov, annual_comp_2yrs_tukey)
# annual_comp_2yrs_df <- enframe(annual_comp_2yrs_cld[[1]]$Letters)
# 
# ##### AQUE #####
# # analysis of variance
# annual_aque_comp_2yrs_aov <- aov(delta_annual ~ comp_2yrs, data = (delta_annual %>% filter(site == "aque")))
# 
# # summary
# summary(annual_aque_comp_2yrs_aov)
# 
# # Tukey and compact letter display
# annual_aque_comp_2yrs_tukey <- TukeyHSD(annual_aque_comp_2yrs_aov)
# annual_aque_comp_2yrs_cld <- multcompLetters4(annual_aque_comp_2yrs_aov, annual_aque_comp_2yrs_tukey)
# annual_aque_comp_2yrs_df <- enframe(annual_aque_comp_2yrs_cld[[1]]$Letters)
# 
# ##### NAPL #####
# # analysis of variance
# annual_napl_comp_2yrs_aov <- aov(delta_annual ~ comp_2yrs, data = (delta_annual %>% filter(site == "napl")))
# 
# # summary
# summary(annual_napl_comp_2yrs_aov)
# 
# # Tukey and compact letter display
# annual_napl_comp_2yrs_tukey <- TukeyHSD(annual_napl_comp_2yrs_aov)
# annual_napl_comp_2yrs_cld <- multcompLetters4(annual_napl_comp_2yrs_aov, annual_napl_comp_2yrs_tukey)
# annual_napl_comp_2yrs_df <- enframe(annual_napl_comp_2yrs_cld[[1]]$Letters)
# 
# ##### IVEE #####
# # analysis of variance
# annual_ivee_comp_2yrs_aov <- aov(delta_annual ~ comp_2yrs, data = (delta_annual %>% filter(site == "ivee")))
# 
# # summary
# summary(annual_ivee_comp_2yrs_aov)
# 
# # Tukey and compact letter display
# annual_ivee_comp_2yrs_tukey <- TukeyHSD(annual_ivee_comp_2yrs_aov)
# annual_ivee_comp_2yrs_cld <- multcompLetters4(annual_ivee_comp_2yrs_aov, annual_ivee_comp_2yrs_tukey)
# annual_ivee_comp_2yrs_df <- enframe(annual_ivee_comp_2yrs_cld[[1]]$Letters)
# 
# ##### MOHK #####
# # analysis of variance
# annual_mohk_comp_2yrs_aov <- aov(delta_annual ~ comp_2yrs, data = (delta_annual %>% filter(site == "mohk")))
# 
# # summary
# summary(annual_mohk_comp_2yrs_aov)
# 
# # Tukey and compact letter display
# annual_mohk_comp_2yrs_tukey <- TukeyHSD(annual_mohk_comp_2yrs_aov)
# annual_mohk_comp_2yrs_cld <- multcompLetters4(annual_mohk_comp_2yrs_aov, annual_mohk_comp_2yrs_tukey)
# annual_mohk_comp_2yrs_df <- enframe(annual_mohk_comp_2yrs_cld[[1]]$Letters)
# 
# ##### CARP #####
# # analysis of variance
# annual_carp_comp_2yrs_aov <- aov(delta_annual ~ comp_2yrs, data = (delta_annual %>% filter(site == "carp")))
# 
# # summary
# summary(annual_carp_comp_2yrs_aov)
# 
# # Tukey and compact letter display
# annual_carp_comp_2yrs_tukey <- TukeyHSD(annual_carp_comp_2yrs_aov)
# annual_carp_comp_2yrs_cld <- multcompLetters4(annual_carp_comp_2yrs_aov, annual_carp_comp_2yrs_tukey)
# annual_carp_comp_2yrs_df <- enframe(annual_carp_comp_2yrs_cld[[1]]$Letters)
# 
# 
# ggplot(delta_annual %>% drop_na(comp_2yrs), aes(x = comp_2yrs, y = delta_annual)) +
#   geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
#   geom_jitter(width = 0.2, alpha = 0.4) +
#   stat_summary(aes(group = site), fun = mean, geom = "line") +
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
#   stat_summary(aes(col = site), fun = mean, geom = "point", size = 6) +
#   scale_color_manual(values = site_palette) +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 18),
#         axis.text = element_text(size = 10),
#         strip.text = element_text(size = 10)) +
#   labs(x = "Time period", y = "\U0394 biomass (treatment - control)",
#        title = "Annual removal, 2 year comparison") +
#   facet_wrap(~site_full)
# 
# 
# 
# 
# ggplot(delta_annual %>% drop_na(comp_3yrs), aes(x = comp_3yrs, y = delta_annual)) +
#   geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
#   geom_jitter(width = 0.2, alpha = 0.4) +
#   stat_summary(aes(group = site), fun = mean, geom = "line") +
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
#   stat_summary(aes(col = site), fun = mean, geom = "point", size = 6) +
#   scale_color_manual(values = site_palette) +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 18),
#         axis.text = element_text(size = 10),
#         strip.text = element_text(size = 10)) +
#   labs(x = "Time period", y = "\U0394 biomass (treatment - control)",
#        title = "Annual removal, 3 year comparison") +
#   facet_wrap(~site_full)
# 
# summary(aov(delta_annual ~ comp_3yrs, data = delta_annual))
# TukeyHSD(aov(delta_annual ~ comp_3yrs, data = delta_annual))
# 
# summary(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "aque"))))
# TukeyHSD(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "aque"))))
# 
# summary(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "napl"))))
# TukeyHSD(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "napl"))))
# 
# summary(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "ivee"))))
# TukeyHSD(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "ivee"))))
# 
# summary(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "mohk"))))
# TukeyHSD(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "mohk"))))
# 
# summary(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "carp"))))
# TukeyHSD(aov(delta_annual ~ comp_3yrs, data = (delta_annual %>% filter(site == "carp"))))
# 
# ```
# 
# 
# 
# 
# 
# 
# # aque_comm <- ggplot() +
#   coord_cartesian() +
#   # points in experiment period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "aque" & exp_dates == "during"), 
#              aes(x = NMDS1, y = NMDS2, shape = treatment), color = "grey", size = 2, alpha = 0.5) +
#   # points in after period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "aque" & exp_dates == "after" & treatment == "continual"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   geom_point(data = bc_continual_plotdf %>% filter(site == "aque" & exp_dates == "after" & treatment == "control"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2, shape = 17) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in control
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "aque" & treatment == "control"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                        yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                        color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in continual removal
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "aque" & treatment == "continual"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                       yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                       color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   
#   # geom_text(data = bc_continual_species, aes(x = NMDS1, y = NMDS2, label = scientific_name), col = "orange") +
#   theme_bw() +
#   theme(panel.grid = element_line(color = "white")) +
#   geom_hline(aes(yintercept = 0), lty = 2, color = "grey") +
#   geom_vline(aes(xintercept = 0), lty = 2, color = "grey") +
#   labs(title = "AQUE")
# 
# napl_comm <- ggplot() +
#   coord_cartesian() +
#   # points in experiment period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "napl" & exp_dates == "during"), 
#              aes(x = NMDS1, y = NMDS2, shape = treatment), color = "grey", size = 2, alpha = 0.5) +
#   # points in after period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "napl" & exp_dates == "after" & treatment == "continual"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   geom_point(data = bc_continual_plotdf %>% filter(site == "napl" & exp_dates == "after" & treatment == "control"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2, shape = 17) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in control
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "napl" & treatment == "control"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                        yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                        color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in continual removal
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "napl" & treatment == "continual"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                       yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                       color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   
#   # geom_text(data = bc_continual_species, aes(x = NMDS1, y = NMDS2, label = scientific_name), col = "orange") +
#   theme_bw() +
#   theme(panel.grid = element_line(color = "white"),
#         legend.position = "none") +
#   geom_hline(aes(yintercept = 0), lty = 2, color = "grey") +
#   geom_vline(aes(xintercept = 0), lty = 2, color = "grey") +
#   labs(title = "NAPL")
# 
# mohk_comm <- ggplot() +
#   coord_cartesian() +
#   # points in experiment period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "mohk" & exp_dates == "during"), 
#              aes(x = NMDS1, y = NMDS2, shape = treatment), color = "grey", size = 2, alpha = 0.5) +
#   # points in after period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "mohk" & exp_dates == "after" & treatment == "continual"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   geom_point(data = bc_continual_plotdf %>% filter(site == "mohk" & exp_dates == "after" & treatment == "control"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2, shape = 17) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in control
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "mohk" & treatment == "control"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                        yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                        color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in continual removal
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "mohk" & treatment == "continual"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                       yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                       color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   
#   # geom_text(data = bc_continual_species, aes(x = NMDS1, y = NMDS2, label = scientific_name), col = "orange") +
#   theme_bw() +
#   theme(panel.grid = element_line(color = "white"),
#         legend.position = "none") +
#   geom_hline(aes(yintercept = 0), lty = 2, color = "grey") +
#   geom_vline(aes(xintercept = 0), lty = 2, color = "grey") +
#   labs(title = "MOHK")
# 
# carp_comm <- ggplot() +
#   coord_cartesian() +
#   # points in experiment period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "carp" & exp_dates == "during"), 
#              aes(x = NMDS1, y = NMDS2, shape = treatment), color = "grey", size = 2, alpha = 0.5) +
#   # points in after period
#   geom_point(data = bc_continual_plotdf %>% filter(site == "carp" & exp_dates == "after" & treatment == "continual"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   geom_point(data = bc_continual_plotdf %>% filter(site == "carp" & exp_dates == "after" & treatment == "control"), 
#              aes(x = NMDS1, y = NMDS2, color = year), size = 2, shape = 17) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in control
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "carp" & treatment == "control"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                        yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                        color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#cdced1", high = "#99adcf", name = "Control years") +
#   new_scale_color() +
#   # arrows for after in continual removal
#   geom_segment(data = bc_continual_plotdf %>% filter(exp_dates == "after" & site == "carp" & treatment == "continual"), aes(x = NMDS1, y = NMDS2, xend = c(tail(NMDS1, n = -1), NA), 
#                                                                                                                       yend = c(tail(NMDS2, n = -1), NA),
#                                                                                                                       color = year),
#                arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
#                size = 0.5) +
#   scale_color_gradient(low = "#7590e0", high = "#032a9c", name = "Removal years") +
#   new_scale_color() +
#   
#   # geom_text(data = bc_continual_species, aes(x = NMDS1, y = NMDS2, label = scientific_name), col = "orange") +
#   theme_bw() +
#   theme(panel.grid = element_line(color = "white"),
#         legend.position = "none") +
#   geom_hline(aes(yintercept = 0), lty = 2, color = "grey") +
#   geom_vline(aes(xintercept = 0), lty = 2, color = "grey") +
#   labs(title = "CARP")
# 
# plots_together <- (aque_comm + napl_comm)/(mohk_comm + carp_comm + plot_spacer()) +
#   plot_layout(guides = "collect")
# plots_together
# 
##### old delta continual

# delta_continual <- delta_biomass %>% 
#   dplyr::select(site, year, month, date, control, continual, delta_continual) %>% 
#   # take out years where continual removal hadn't happened yet
#   drop_na(delta_continual) %>% 
#   mutate(exp_dates = case_when(
#     # after for continual removal:
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
#   left_join(., enframe(sites_full), by = c("site" = "name")) %>% 
#   rename("site_full" = value)

# total biomass
# kelp_biomass <- biomass %>% 
#   filter(sp_code == "MAPY") %>% 
#   dplyr::select(-sp_code) %>% 
#   # make a new column for after dates
#   after_dates_column() %>% 
#   # make a new column for during and after and set factor levels
#   exp_dates_column()
# 
# # delta biomass
# delta_biomass <- kelp_biomass %>% 
#   dplyr::select(site, year, month, treatment, date, exp_dates, dry_gm2) %>% 
#   pivot_wider(names_from = treatment, values_from = dry_gm2) %>% 
#   mutate(delta_annual = annual - control,
#          delta_continual = continual - control) 

#### code for centroid segment arrows

# geom_segment(data = aque_centroids_contin, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, color = aque_col) +
# geom_segment(data = napl_centroids_contin, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, color = napl_col) +
# # geom_point(data = mohk_centroids_contin, 
# #            aes(x = NMDS1, y = NMDS2, shape = site_full), size = 4) +
# geom_segment(data = mohk_centroids_contin, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, color = mohk_col) +
# geom_point(data = carp_centroids_contin, 
#            aes(x = NMDS1, y = NMDS2, shape = site_full), size = 4) +
# geom_segment(data = carp_centroids_contin, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, color = carp_col) +
#   geom_segment(data = aque_centroids_cont, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, lty = 2, color = aque_col) +
#     geom_segment(data = napl_centroids_cont, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, lty = 2, color = napl_col) +
#       geom_segment(data = mohk_centroids_cont, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, lty = 2, color = mohk_col) +
#         geom_segment(data = carp_centroids_cont, 
#              aes(x = NMDS1, y = NMDS2, 
#                  xend = c(tail(NMDS1, n = -1), NA), 
#                  yend = c(tail(NMDS2, n = -1), NA)), 
#              arrow = arrow(length = unit(0.5, "cm")),
#              size = 1, lty = 2, color = carp_col) +
# stat_ellipse(aes(group = comp_2yrs)) +

##########################################################################-
# 3. start-during-after comparisons ---------------------------------------
##########################################################################-

# ⟞ a. algae -------------------------------------------------------------

# anova with random effect
# anova_algae_1yr <- lmer(delta_continual_algae ~ comp_1yr + (1|site), 
#                          data = delta_algae_continual %>% drop_na(comp_1yr))
# anova_algae_2yrs <- lmer(delta_continual_algae ~ comp_2yrs + (1|site), 
#                     data = delta_algae_continual %>% drop_na(comp_2yrs))
# anova_algae_3yrs <- lmer(delta_continual_algae ~ comp_3yrs + (1|site), 
#                          data = delta_algae_continual %>% drop_na(comp_3yrs))
# 
# anova_raw_algae_1yr <- lmer(algae_biomass ~ comp_1yr*treatment + (1|site),
#                              data = algae_continual_long %>% drop_na(comp_1yr))
# anova_raw_algae_2yrs <- lmer(algae_biomass ~ comp_2yrs*treatment + (1|site),
#                              data = algae_continual_long %>% drop_na(comp_2yrs))
# anova_raw_algae_3yrs <- lmer(algae_biomass ~ comp_3yrs*treatment + (1|site),
#                              data = algae_continual_long %>% drop_na(comp_3yrs))

# diagnostics
# plot(simulateResiduals(anova_algae_1yr))
# check_model(anova_algae_1yr)
# 
# plot(simulateResiduals(anova_algae_2yrs))
# check_model(anova_algae_2yrs)
# 
# plot(simulateResiduals(anova_algae_3yrs))
# check_model(anova_algae_3yrs)
# 
# plot(simulateResiduals(anova_raw_algae_1yr))
# check_model(anova_raw_algae_1yr)
# 
# plot(simulateResiduals(anova_raw_algae_2yrs))
# check_model(anova_raw_algae_2yrs)
# 
# plot(simulateResiduals(anova_raw_algae_3yrs))
# check_model(anova_raw_algae_3yrs)

# summary
# summary(anova_algae_1yr) # difference when comparing start and after
# summary(anova_algae_2yrs)
# summary(anova_algae_3yrs) # same as 2 years
# 
# summary(anova_raw_algae_2yrs)
# summary(anova_raw_algae_3yrs)

# plot predictions
# plot(ggpredict(anova_raw_algae_1yr, terms = c("comp_1yr", "treatment"), type = "fixed")) 
# plot(ggpredict(anova_raw_algae_2yrs, terms = c("comp_2yrs", "treatment"), type = "fixed")) 
# plot(ggpredict(anova_raw_algae_3yrs, terms = c("comp_3yrs", "treatment"), type = "fixed")) 

# anova table
# anova_algae_1yr_aovt <- anova(anova_algae_1yr, ddf = c("Kenward-Roger"))
# 
# anova_algae_2yrs_aovt <- anova(anova_algae_2yrs, ddf = c("Kenward-Roger"))
# 
# anova_algae_3yrs_aovt <- anova(anova_algae_3yrs, ddf = c("Kenward-Roger")) 
# 
# anova_algae_1yr_diff <- delta_algae_continual %>% 
#   filter(comp_1yr %in% c("start", "during", "after")) %>% 
#   group_by(comp_1yr) %>% 
#   summarize(mean = mean(delta_continual_algae),
#             se = se(delta_continual_algae))
# 
# anova_algae_2yrs_diff <- delta_algae_continual %>% 
#   filter(comp_2yrs %in% c("start", "during", "after")) %>% 
#   group_by(comp_2yrs) %>% 
#   summarize(mean = mean(delta_continual_algae),
#             se = se(delta_continual_algae))
# 
# anova_algae_3yrs_diff <- delta_algae_continual %>% 
#   filter(comp_3yrs %in% c("start", "during", "after")) %>% 
#   group_by(comp_3yrs) %>% 
#   summarize(mean = mean(delta_continual_algae),
#             se = se(delta_continual_algae))

# least squares comparison

# anova_algae_2yrs_summary <- difflsmeans_summary_fxn(anova_algae_2yrs)
# 
# anova_algae_3yrs_summary <- difflsmeans_summary_fxn(anova_algae_3yrs)

# extract predicted values
# anova_algae_2yrs_df <- ggpredict(anova_algae_2yrs, terms = "comp_2yrs", type = "fixed") %>% 
#   mutate(x = case_when(
#     x == "start" ~ "Start of removal",
#     x == "during" ~ "End of removal",
#     x == "after" ~ "Recovery period"
#   )) %>% 
#   mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))
# anova_algae_3yrs_df <- ggpredict(anova_algae_3yrs, terms = "comp_3yrs", type = "fixed") %>% 
#   mutate(x = case_when(
#     x == "start" ~ "Start of removal",
#     x == "during" ~ "End of removal",
#     x == "after" ~ "Recovery period"
#   )) %>% 
#   mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))

# ⟞ b. epilithic invertebrates -------------------------------------------

# anova with random effect
# anova_epi_1yr <- lmer(delta_continual_epi ~ comp_1yr + (1|site), 
#                         data = delta_epi_continual %>% drop_na(comp_1yr))
# anova_epi_2yrs <- lmer(delta_continual_epi ~ comp_2yrs + (1|site), 
#                          data = delta_epi_continual %>% drop_na(comp_2yrs))
# anova_epi_3yrs <- lmer(delta_continual_epi ~ comp_3yrs + (1|site), 
#                          data = delta_epi_continual %>% drop_na(comp_3yrs))

# diagnostics
# plot(simulateResiduals(anova_epi_2yrs))
# check_model(anova_epi_2yrs)
# 
# plot(simulateResiduals(anova_epi_3yrs))
# check_model(anova_epi_3yrs)

# summary
# summary(anova_epi_2yrs)
# summary(anova_epi_3yrs)

# anova table
# anova_epi_2yrs_aovt <- anova(anova_epi_2yrs, ddf = c("Kenward-Roger")) 
# 
# anova_epi_3yrs_aovt <- anova(anova_epi_3yrs, ddf = c("Kenward-Roger"))
# 
# anova_epi_2yrs_diff <- emmeans(anova_epi_2yrs, "comp_2yrs", lmer.df = "kenward-roger") %>% 
#   test() %>% 
#   as.data.frame() 
# 
# anova_epi_3yrs_diff <- emmeans(anova_epi_3yrs, "comp_3yrs", lmer.df = "kenward-roger") %>% 
#   test() %>% 
#   as.data.frame() 

# least squares comparison
# anova_epi_2yrs_summary <- difflsmeans_summary_fxn(anova_epi_2yrs)
# 
# anova_epi_3yrs_summary <- difflsmeans_summary_fxn(anova_epi_3yrs)

# extract predicted values
# anova_epi_2yrs_df <- ggpredict(anova_epi_2yrs, terms = "comp_2yrs", type = "fixed") %>% 
#   mutate(x = case_when(
#     x == "start" ~ "Start of removal",
#     x == "during" ~ "End of removal",
#     x == "after" ~ "Recovery period"
#   )) %>% 
#   mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))
# 
# anova_epi_3yrs_df <- ggpredict(anova_epi_3yrs, terms = "comp_3yrs", type = "fixed") %>% 
#   mutate(x = case_when(
#     x == "start" ~ "Start of removal",
#     x == "during" ~ "End of removal",
#     x == "after" ~ "Recovery period"
#   )) %>% 
#   mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))

# ⟞ c. endolithic invertebrates ------------------------------------------

# anova with random effect
# anova_endo_2yrs <- lmer(delta_continual_endo ~ comp_2yrs + (1|site), 
#                        data = delta_endo_continual %>% drop_na(comp_2yrs))
# anova_endo_3yrs <- lmer(delta_continual_endo ~ comp_3yrs + (1|site), 
#                        data = delta_endo_continual %>% drop_na(comp_3yrs))

# diagnostics
# plot(simulateResiduals(anova_endo_2yrs))
# check_model(anova_endo_2yrs)
# 
# plot(simulateResiduals(anova_endo_3yrs))
# check_model(anova_endo_3yrs)

# summary
# summary(anova_endo_2yrs)
# summary(anova_endo_3yrs) # same as 2 years

# anova table
# anova_endo_2yrs_aovt <- anova(anova_endo_2yrs, ddf = c("Kenward-Roger"))
# 
# anova_endo_3yrs_aovt <- anova(anova_endo_3yrs, ddf = c("Kenward-Roger")) 
# 
# anova_endo_2yrs_diff <- emmeans(anova_endo_2yrs, "comp_2yrs", lmer.df = "kenward-roger") %>% 
#   test() %>% 
#   as.data.frame() 
# 
# anova_endo_3yrs_diff <- emmeans(anova_endo_3yrs, "comp_3yrs", lmer.df = "kenward-roger") %>% 
#   test() %>% 
#   as.data.frame() 

# least squares comparison
# anova_endo_2yrs_summary <- difflsmeans_summary_fxn(anova_endo_2yrs)
# 
# anova_endo_3yrs_summary <- difflsmeans_summary_fxn(anova_endo_3yrs)

# extract predicted values
# anova_endo_2yrs_df <- ggpredict(anova_endo_2yrs, terms = "comp_2yrs", type = "fixed") %>% 
#   mutate(x = case_when(
#     x == "start" ~ "Start of removal",
#     x == "during" ~ "End of removal",
#     x == "after" ~ "Recovery period"
#   )) %>% 
#   mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))
# anova_endo_3yrs_df <- ggpredict(anova_endo_3yrs, terms = "comp_3yrs", type = "fixed") %>% 
#   mutate(x = case_when(
#     x == "start" ~ "Start of removal",
#     x == "during" ~ "End of removal",
#     x == "after" ~ "Recovery period"
#   )) %>% 
#   mutate(x = fct_relevel(x, "Start of removal", "End of removal", "Recovery period"))

# ⟞ d. figures -----------------------------------------------------------

# ⟞ ⟞ i. algae -----------------------------------------------------------

# sda_algae_biomass <- ggplot(anova_algae_2yrs_df) +
#   # horizontal line at 0
#   geom_hline(yintercept = 0, lty = 2) +
#   # annotate("rect", xmin = 0, xmax = 4, ymin = 125.5, ymax = 140, fill = "#FFFFFF") +
#   # points and error bars
#   geom_point(aes(x = x, y = predicted), size = 2) +
#   geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
#   # path between time periods
#   geom_line(aes(x = x, y = predicted, group = 2)) +
#   # add in post-hoc comparisons
#   annotate("text", x = 1, y = 122, label = "a", size = 3) +
#   annotate("text", x = 2, y = 122, label = "b", size = 3) +
#   annotate("text", x = 3, y = 122, label = "a", size = 3) +
#   # annotate("text", x = 0.55, y = 135, label = "Understory algae", size = 10) +
#   # aesthetics
#   scale_x_discrete(labels = wrap_format(10)) +
#   scale_y_continuous(limits = c(-30, 130), breaks = c(0, 50, 100, 150), expand = c(0, 0)) +
#   sda_biomass_theme() +
#   labs(x = "Time period", 
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(a)")
# sda_algae_biomass

# sda_algae_biomass_3yrs <- ggplot(anova_algae_3yrs_df) +
#   # horizontal line at 0
#   geom_hline(yintercept = 0, lty = 2) +
#   # annotate("rect", xmin = 0, xmax = 4, ymin = 125.5, ymax = 140, fill = "#FFFFFF") +
#   # points and error bars
#   geom_point(aes(x = x, y = predicted), size = 2) +
#   geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
#   # path between time periods
#   geom_line(aes(x = x, y = predicted, group = 2)) +
#   # add in post-hoc comparisons
#   annotate("text", x = 1, y = 135, label = "a", size = 3) +
#   annotate("text", x = 2, y = 135, label = "b", size = 3) +
#   annotate("text", x = 3, y = 135, label = "a", size = 3) +
#   # annotate("text", x = 0.55, y = 135, label = "Understory algae", size = 10) +
#   # aesthetics
#   scale_x_discrete(labels = wrap_format(10)) +
#   scale_y_continuous(limits = c(-30, 140), breaks = c(0, 50, 100, 150), expand = c(0, 0)) +
#   sda_biomass_theme() +
#   labs(x = "Time period", 
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(a) Understory macroalgae")
# sda_algae_biomass_3yrs

# ⟞ ⟞ ii. epi inverts ----------------------------------------------------

# sda_epi_biomass <- ggplot(anova_epi_2yrs_df) +
#   # horizontal line at 0
#   geom_hline(yintercept = 0, lty = 2) +
#   # points and error bars
#   geom_point(aes(x = x, y = predicted), size = 2) +
#   # annotate("rect", xmin = 0, xmax = 4, ymin = 20.1, ymax = 24, fill = "#FFFFFF") +
#   geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
#   # path between time periods
#   geom_line(aes(x = x, y = predicted, group = 1)) +
#   # add in post-hoc comparisons
#   annotate("text", x = 1, y = 20, label = "a", size = 3) +
#   annotate("text", x = 2, y = 20, label = "b", size = 3) +
#   annotate("text", x = 3, y = 20, label = "b", size = 3) +
#   # annotate("text", x = 0.68, y = 23, label = "Epilithic invertebrates", size = 10) +
#   # aesthetics
#   scale_x_discrete(labels = wrap_format(10)) +
#   scale_y_continuous(limits = c(-13, 22), expand = c(0, 0)) +
#   sda_biomass_theme() +
#   labs(x = "Time period", 
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(c)")
# sda_epi_biomass

# sda_epi_biomass_3yrs <- ggplot(anova_epi_3yrs_df) +
#   # horizontal line at 0
#   geom_hline(yintercept = 0, lty = 2) +
#   # points and error bars
#   geom_point(aes(x = x, y = predicted), size = 2) +
#   # annotate("rect", xmin = 0, xmax = 4, ymin = 20.1, ymax = 24, fill = "#FFFFFF") +
#   geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
#   # path between time periods
#   geom_line(aes(x = x, y = predicted, group = 1)) +
#   # add in post-hoc comparisons
#   annotate("text", x = 1, y = 20, label = "a", size = 3) +
#   annotate("text", x = 2, y = 20, label = "b", size = 3) +
#   annotate("text", x = 3, y = 20, label = "c", size = 3) +
#   # annotate("text", x = 0.68, y = 23, label = "Epilithic invertebrates", size = 10) +
#   # aesthetics
#   scale_x_discrete(labels = wrap_format(10)) +
#   scale_y_continuous(limits = c(-13, 22), expand = c(0, 0)) +
#   sda_biomass_theme() +
#   labs(x = "Time period", 
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(b) Epilithic invertebrates")
# sda_epi_biomass_3yrs

# ⟞ ⟞ iii. endo inverts --------------------------------------------------

# sda_endo_biomass <- ggplot(anova_endo_2yrs_df) +
#   # horizontal line at 0
#   geom_hline(yintercept = 0, lty = 2) +
#   # annotate("rect", xmin = 0, xmax = 4, ymin = 551, ymax = 650, fill = "#FFFFFF") +
#   # points and error bars
#   geom_point(aes(x = x, y = predicted), size = 2) +
#   geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
#   # path between time periods
#   geom_line(aes(x = x, y = predicted, group = 1)) +
#   # aesthetics
#   scale_x_discrete(labels = wrap_format(10)) +
#   scale_y_continuous(limits = c(-100, 550), expand = c(0, 0), breaks = c(0, 200, 400, 600)) +
#   sda_biomass_theme() +
#   labs(x = "Time period", 
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(e)")
# sda_endo_biomass

# sda_endo_biomass_3yrs <- ggplot(anova_endo_3yrs_df) +
#   # horizontal line at 0
#   geom_hline(yintercept = 0, lty = 2) +
#   # annotate("rect", xmin = 0, xmax = 4, ymin = 551, ymax = 650, fill = "#FFFFFF") +
#   # points and error bars
#   geom_point(aes(x = x, y = predicted), size = 2) +
#   geom_errorbar(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2) +
#   # path between time periods
#   geom_line(aes(x = x, y = predicted, group = 1)) +
#   # aesthetics
#   scale_x_discrete(labels = wrap_format(10)) +
#   scale_y_continuous(limits = c(-100, 550), expand = c(0, 0), breaks = c(0, 200, 400, 600)) +
#   sda_biomass_theme() +
#   labs(x = "Time period", 
#        y = "\U0394 biomass \n (treatment - control)",
#        title = "(c) Endolithic invertebrates")
# sda_endo_biomass_3yrs

##########################################################################-
# 6. endo. invert linear model --------------------------------------------
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
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Slope**",
    df = "**df**"
  ) 
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
  bold_p(t = 0.05) %>% 
  modify_header(
    label = " ",
    estimate = "**Slope**",
    df = "**df**"
  ) 
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
       subtitle = "(f)")
endo_time

# ⟞ b. s-d-a summary tables -----------------------------------------------

sda_2yrs_tables <- rbind(anova_algae_2yrs_summary, anova_epi_2yrs_summary, anova_endo_2yrs_summary) %>% 
  rename_with(., .fn = ~paste(., "_2yrs", sep = "", .cols = everything(cols)))

sda_3yrs_tables <- rbind(anova_algae_3yrs_summary, anova_epi_3yrs_summary, anova_endo_3yrs_summary) %>% 
  rename_with(., .fn = ~paste(., "_3yrs", sep = "", .cols = everything(cols)))

sda_together_tables <- cbind(sda_2yrs_tables, sda_3yrs_tables) %>% 
  select(-levels_3yrs1) %>% 
  # turn the whole thing into a gt
  gt() %>% 
  # group labels
  tab_row_group(
    label = "Understory macroalgae", rows = 1:3
  ) %>% 
  tab_row_group(
    label = "Epilithic invertebrates", rows = 4:6
  ) %>% 
  tab_row_group(
    label = "Endolithic invertebrates", rows = 7:9
  ) %>% 
  row_group_order(groups = c("Understory macroalgae", "Epilithic invertebrates", "Endolithic invertebrates")) %>% 
  # spanner labels
  tab_spanner(
    label = "2 year comparison",
    id = "2 year comparison",
    columns = c(estimate_2yrs1, std_error_2yrs1, df_2yrs1, t_value_2yrs1, pr_t_2yrs1)
  ) %>%
  tab_spanner(
    label = "3 year comparison",
    id = "3 year comparison",
    columns = c(estimate_3yrs1, std_error_3yrs1, df_3yrs1, t_value_3yrs1, pr_t_3yrs1)
  ) %>% 
  # change column names
  cols_label(
    levels_2yrs1 = "",
    estimate_2yrs1 = "Estimated difference",
    std_error_2yrs1 = "SE",
    df_2yrs1 = "df",
    t_value_2yrs1 = "t-value", 
    pr_t_2yrs1 = "p-value",
    estimate_3yrs1 = "Estimated difference",
    std_error_3yrs1 = "SE",
    df_3yrs1 = "df",
    t_value_3yrs1 = "t-value", 
    pr_t_3yrs1 = "p-value",
  ) %>% 
  # align columns
  cols_align(columns = everything(),
             align = "center") %>% 
  # bold p < 0.05
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = pr_t_2yrs1,
      rows = pr_t_2yrs1 < 0.05
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = pr_t_3yrs1,
      rows = pr_t_3yrs1 < 0.05
    )
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>% 
  tab_options(table.font.names = "Times New Roman") 

# gtsave(sda_together_tables,
#        here::here("tables", "ms-tables", paste("tbl-S3_", today(), ".docx", sep = "")),
#        vwidth = 1500, vheight = 1000)

# ⟞ a. s-d-a and model predictions ----------------------------------------

# putting group plots together
top <- plot_grid(sda_algae_biomass, algae_time, ncol = 2)
middle <- plot_grid(sda_epi_biomass, epi_time, ncol = 2)
bottom <- plot_grid(sda_endo_biomass, endo_time, ncol = 2)

# putting group plots with labels
algae <- plot_grid(algae_label, top, ncol = 1, rel_heights = c(1, 12))
epi <- plot_grid(epi_label, middle, ncol = 1, rel_heights = c(1, 12))
endo <- plot_grid(endo_label, bottom, ncol = 1, rel_heights = c(1, 12))

# putting plots together
algae_epi <- plot_grid(algae, epi, ncol = 1)
sda_time_together <- plot_grid(algae_epi, endo, ncol = 1, rel_heights = c(2, 1))

# ggsave(here::here("figures", "ms-figures",
#                   paste("fig-2_", today(), ".jpg", sep = "")),
#        plot = sda_time_together,
#        height = 18, width = 16, units = "cm",
#        dpi = 400)

# ⟞ c. s-d-a with 3 year comparison ---------------------------------------

# sda_3yrs <- sda_algae_biomass_3yrs / sda_epi_biomass_3yrs / sda_endo_biomass_3yrs
# 
# ggsave(here::here("figures", "ms-figures",
#        paste("fig-S7_", today(), ".jpg", sep = ")),
#        plot = sda_3yrs,
#        height = 18, width = 9, units = "cm",
#        dpi = 400)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------OLD CODE BELOW HERE---------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ##########################################################################-
# # 3. linear models --------------------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. during removal -----------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # delta_hist <- delta_continual %>% 
# #   filter(exp_dates == "during") %>% 
# #   ggplot(aes(x = delta_continual)) +
# #   geom_histogram(bins = 10)
# # delta_hist
# 
# # model
# # normal model
# # lm_kelp_during_lmer <- lmerTest::lmer(
# #   delta_continual ~ time_since_end + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass)
# # normal model with zero inflated negative binomial structure
# # lm_kelp_during_zinb <- glmmTMB(
# #   delta_continual ~ time_since_end + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "during"),
# #   ziformula = ~time_since_end, 
# #   family = nbinom2,
# #   na.action = na.pass
# # )
# # normal model using glmmTMB with tweedie structure
# # lm_kelp_during_tweedie <- glmmTMB(
# #   delta_continual ~ time_since_end + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   family = tweedie(),
# #   na.action = na.pass, 
# #   REML = TRUE)
# # normal model with season
# # lm_kelp_during_season <- lmerTest::lmer(
# #   delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass)
# # normal model with season using glmmTM with tweedie structure
# # lm_kelp_during_season_tweedie <- glmmTMB(
# #   delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   family = tweedie(),
# #   na.action = na.pass)
# # normal model with season as random effect (failed to converge)
# # lm_kelp_during_season_re <- lmerTest::lmer(
# #   delta_continual ~ time_since_end + (1|site) + (1|quarter),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass)
# # continuous AR1
# # lm_kelp_during_lme_car1 <- nlme::lme(
# #   delta_continual ~ time_since_end, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corCAR1())
# # continuous AR1 with season
# # lm_kelp_during_lme_car1_season <- nlme::lme(
# #   delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corCAR1())
# # continuous AR1 with season as random effect
# # lm_kelp_during_lme_car1_season_re <- nlme::lme(
# #   delta_continual ~ time_since_end, random = list(site = ~1, quarter = ~1),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corCAR1())
# # AR1
# # lm_kelp_during_lme_ar1 <- nlme::lme(
# #   delta_continual ~ time_since_end, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corAR1())
# # AR1 with season
# # lm_kelp_during_lme_ar1_season <- nlme::lme(
# #   delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corAR1())
# # AR1 with season as random effect
# # lm_kelp_during_lme_ar1_season_re <- nlme::lme(
# #   delta_continual ~ time_since_end, random = list(site = ~1, quarter = ~1),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corAR1())
# # ARMA 4
# # lm_kelp_during_lme_ar4 <- nlme::lme(
# #   delta_continual ~ time_since_end, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corARMA(p = 4, q = 0)) 
# # with season no autoregressive
# # lm_kelp_during_glmmtmb_season <- glmmTMB(
# #   delta_continual ~ time_since_end + quarter + time_since_end*quarter + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   family = tweedie(),
# #   na.action = na.pass) 
# 
# # ARMA 4 with season
# # lm_kelp_during_lme_ar4_season <- nlme::lme(
# #   delta_continual ~ time_since_end + quarter + time_since_end*quarter, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corARMA(p = 4, q = 0)) 
# # with season
# # lm_kelp_during_glmmTMB_season <- glmmTMB(
# #   delta_continual ~ quarter + ar1(time_since_end + 0 |site),
# #   data = delta_continual %>% filter(exp_dates == "during") %>% mutate(time_since_end = as_factor(time_since_end)), 
# #   family = tweedie(),
# #   na.action = na.pass) 
# # ARMA 4 with season as random effect
# # lm_kelp_during_lme_ar4_season_re <- nlme::lme(
# #   delta_continual ~ time_since_end, random = list(site = ~1, quarter = ~1),
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   correlation = corARMA(p = 4, q = 0)) 
# # GLS CAR1
# # lm_kelp_during_gls_car1 <- nlme::gls(
# #   delta_continual ~ time_since_end,
# #   data = delta_continual %>% filter(exp_dates == "during"), 
# #   correlation = corCAR1(form = ~ 1|site)
# # )
# # linear transformation: multiply by -1, add 8
# # transform <- delta_continual %>% 
# #   mutate(delta_continual_tf = delta_continual*-1 + 8) %>% 
# #   mutate(logratio = log(continual)/log(control)) %>% 
# #   filter(exp_dates == "during")
# # ggplot(transform, aes(x = logratio)) +
# #   geom_histogram()
# # 
# # ggplot(transform, aes(x = control)) +
# #   geom_histogram(bins = 10)
# # transform_fxn <- function(x) x*-1 + 8
# # backtransform <- function(x) (x-8)*-1
# # taking out aque, mohk, or carp makes residuals look normal
# # lm_kelp_during_tf_gamma <- glmmTMB::glmmTMB(
# #   transform_fxn(delta_continual) ~ time_since_end + (1|site),
# #   data = transform,
# #   na.action = na.pass,
# #   family = Gamma(link = "log"))
# # simulateResiduals(lm_kelp_during_tf_gamma, plot = TRUE)
# # summary(lm_kelp_during_tf_gamma)
# # tf_pred <- ggpredict(lm_kelp_during_tf_gamma, terms = ~ time_since_end) %>% 
# #   mutate(transform = backtransform(predicted),
# #          tf_conf.low = backtransform(conf.low),
# #          tf_conf.high = backtransform(conf.high))
# # ggplot(delta_continual, aes(x = time_since_end, y = delta_continual)) +
# #   geom_point() +
# #   geom_line(data = tf_pred, aes(x = x, y = transform), color = "blue") +
# #   geom_ribbon(data = tf_pred, aes(x = x, y = transform, ymin = tf_conf.low, ymax = tf_conf.high), alpha = 0.2)
# # lm_kelp_during_tf_tweedie <- glmmTMB::glmmTMB(
# #   delta_continual_tf ~ time_since_end + (1|site),
# #   data = transform, 
# #   na.action = na.pass,
# #   family = tweedie(link = "log")) 
# # lm_kelp_during_tf_season <- glmmTMB::glmmTMB(
# #   delta_continual_tf ~ time_since_end + quarter + time_since_end*quarter + (1|site),
# #   data = transform, 
# #   na.action = na.pass,
# #   family = Gamma(link = "inverse"),
# #   start = list(psi = c(-1, 1))) 
# 
# # raw kelp biomass model
# ## site random effect
# # lm_kelp_during_zigamma_01 <- glmmTMB(
# #   kelp_biomass ~ time_since_end*treatment + (1|site),
# #   data = continual_long %>% filter(exp_dates == "during"), 
# #   family = ziGamma(link = "log"),
# #   ziformula = ~1)
# # 
# # ## site and year random effect
# # lm_kelp_during_zigamma_02 <- glmmTMB(
# #   kelp_biomass ~ time_since_end*treatment + (1|site) + (1|year),
# #   data = continual_long %>% filter(exp_dates == "during"), 
# #   family = ziGamma(link = "log"),
# #   ziformula = ~1)
# 
# # df <- delta_continual %>% 
# #   filter(exp_dates == "during") %>% 
# #   cbind(., residuals(lm_kelp_during_lmer)) %>% 
# #   dplyr::rename(resid_m1 = 'residuals(lm_kelp_during_lmer)') %>% 
# #   cbind(., residuals(lm_kelp_during_season)) %>% 
# #   dplyr::rename(resid_m2 = 'residuals(lm_kelp_during_season)') %>% 
# #   cbind(., residuals(lm_kelp_during_lme_ar4_season)) %>% 
# #   dplyr::rename(resid_m3 = 'residuals(lm_kelp_during_lme_ar4_season)')
# # 
# # ggplot(df, aes(x = time_since_end, y = resid_m3)) +
# #   geom_point() +
# #   geom_smooth(method = "lm", se = FALSE)
# 
# # diagnostics
# # normal model
# # DHARMa::simulateResiduals(lm_kelp_during_lmer, plot = T)
# # performance::check_model(lm_kelp_during_lmer)
# # performance::check_autocorrelation(lm_kelp_during_lmer) # Durbin-Watson-Test
# 
# # normal model with season
# # simulateResiduals(lm_kelp_during_season, plot = T)
# # check_model(lm_kelp_during_season)
# 
# # normal model with season with random effect
# # simulateResiduals(lm_kelp_during_season_re, plot = T)
# # check_model(lm_kelp_during_season_re)
# 
# # continuous AR1
# # resid_plot_fxn(lm_kelp_during_lme_car1)
# # plot(density(resid(lm_kelp_during_lme_car1)))
# # check_model(lm_kelp_during_lme_car1)
# 
# # continuous AR1 with season
# # resid_plot_fxn(lm_kelp_during_lme_car1_season)
# # plot(density(resid(lm_kelp_during_lme_car1_season)))
# # check_model(lm_kelp_during_lme_car1_season)
# 
# # continuous AR1 with season as random effect
# # resid_plot_fxn(lm_kelp_during_lme_car1_season_re)
# # plot(density(resid(lm_kelp_during_lme_car1_season_re)))
# # check_model(lm_kelp_during_lme_car1_season_re)
# 
# # AR1
# # resid_plot_fxn(lm_kelp_during_lme_ar1)
# # check_model(lm_kelp_during_lme_ar1)
# 
# # AR1 with season
# # resid_plot_fxn(lm_kelp_during_lme_ar1_season)
# 
# # AR1 with season as random effect
# # resid_plot_fxn(lm_kelp_during_lme_ar1_season_re)
# 
# # ARMA 4
# # resid_plot_fxn(lm_kelp_during_lme_ar4)
# # qqnorm(residuals(lm_kelp_during_lme_ar4))
# # qqline(residuals(lm_kelp_during_lme_ar4))
# # ks.test(residuals(lm_kelp_during_lme_ar4), "pnorm")
# 
# # ARMA 4 with season
# # resid_plot_fxn(lm_kelp_during_lme_ar4_season)
# # check_model(lm_kelp_during_lme_ar4_season)
# 
# # ARMA 4 with season as random effect
# # resid_plot_fxn(lm_kelp_during_lme_ar4_season_re)
# # check_model(lm_kelp_during_lme_ar4_season_re)
# 
# # GLS CAR1
# # check_model(lm_kelp_during_gls_car1)
# 
# # model checks
# # normal model
# # testUniformity(lm_kelp_during_lmer)
# # performance::check_convergence(lm_kelp_during_lmer)
# # performance::check_normality(lm_kelp_during_lmer) # Shapiro test
# # performance::check_homogeneity(lm_kelp_during_lmer) # Bartlett test
# # performance::check_heteroskedasticity(lm_kelp_during_lmer) # Breusch-Pagan test
# 
# # plot ACF/PACF
# # normal model
# # acf(resid(lm_kelp_during_lmer))
# # pacf(resid(lm_kelp_during_lmer))
# 
# # continuous AR1
# # acf(resid(lm_kelp_during_lme_car1))
# # pacf(resid(lm_kelp_during_lme_car1))
# 
# # ARMA 4
# # acf(resid(lm_kelp_during_lme_ar4))
# # pacf(resid(lm_kelp_during_lme_ar4))
# 
# # GLS CAR1
# # acf(resid(lm_kelp_during_gls_car1))
# # pacf(resid(lm_kelp_during_gls_car1))
# 
# # raw kelp biomass model
# # DHARMa::simulateResiduals(lm_kelp_during_zigamma_01, plot = T)
# # performance::check_model(lm_kelp_during_zigamma_01)
# 
# # DHARMa::simulateResiduals(lm_kelp_during_zigamma_02, plot = T)
# # performance::check_model(lm_kelp_during_zigamma_02)
# 
# # Rsquared
# # MuMIn::r.squaredGLMM(lm_kelp_during_lmer)
# # r.squaredGLMM(lm_kelp_during_lme_car1)
# # r.squaredGLMM(lm_kelp_during_lme_ar4)
# # r.squaredGLMM(lm_kelp_during_zigamma_01)
# # r.squaredGLMM(lm_kelp_during_zigamma_02)
# 
# # summaries
# # summary(lm_kelp_during_lmer)
# # summary(lm_kelp_during_zigamma_01)
# summary(lm_kelp_during_zigamma_02)
# 
# lm_kelp_during_zigamma_summary <- lm_kelp_during_zigamma_02 %>% 
#   tbl_regression(intercept = TRUE) %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Estimate**"
#   ) %>% 
#   modify_column_indent(
#     columns = label, 
#     rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))
# 
# # filter out zero-inflated component
# lm_kelp_during_zigamma_summary$table_body <- lm_kelp_during_zigamma_summary$table_body %>% 
#   filter(component != "zi")
# # change labels
# lm_kelp_during_zigamma_summary$table_body$label <- c(
#   `(Intercept)` = "(Intercept)",
#   time_since_end = "Time since end",
#   treatmentremoval = "Treatment (removal)",
#   `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
# )
# 
# # final table 
# lm_kelp_during_zigamma_summary
# 
# 
# # lm_kelp_during_summary <- lm_kelp_during_lmer %>% 
# #   tbl_regression() %>% 
# #   bold_p(t = 0.05) %>% 
# #   modify_header(
# #     label = " ",
# #     estimate = "**Slope**",
# #     df = "**df**"
# #   ) 
# # lm_kelp_during_summary
# 
# # AIC comparison
# # MuMIn::AICc(lm_kelp_during_lmer, 
# #             lm_kelp_during_season, lm_kelp_during_season_re,
# #             lm_kelp_during_lme_car1, lm_kelp_during_lme_car1_season,
# #             lm_kelp_during_lme_car1_season_re, lm_kelp_during_lme_ar1,
# #             lm_kelp_during_lme_ar1_season, lm_kelp_during_lme_ar1_season_re,
# #             lm_kelp_during_lme_ar4, lm_kelp_during_lme_ar4_season,
# #             lm_kelp_during_lme_ar4_season_re, lm_kelp_during_gls_car1, 
# #             lm_kelp_during_tweedie, lm_kelp_during_season_tweedie
# #             ) %>% 
# #   arrange(AICc)
# 
# # MuMIn::AICc(lm_kelp_during_zigamma_01, lm_kelp_during_zigamma_02)
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # predicted_kelp_during_overall <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "fixed")
# # predicted_kelp_during_aque <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
# # predicted_kelp_during_napl <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
# # predicted_kelp_during_mohk <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
# # predicted_kelp_during_carp <- ggpredict(lm_kelp_during_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))
# 
# # raw kelp biomass
# predicted_kelp_during_raw <- ggpredict(lm_kelp_during_zigamma_02,
#                                       terms = c("time_since_end [-7.25:0 by = 0.25]", "treatment"), type = "fixed")
# 
# 
# # ⟞ ⟞ iii. means ----------------------------------------------------------
# 
# # estimates predicted mean values of response for each level in predictor
# # using modelbased
# during_response_means <- estimate_means(lm_kelp_during_zigamma_02, 
#                                transform = "response") %>% 
#   mutate(exp_dates = "during")
# # using emmeans
# during_response_emmeans <- emmeans(lm_kelp_during_zigamma_02, 
#                           ~ treatment | time_since_end, type = "response") 
# 
# # calculates differences between levels of categorical predictor
# contrasts <- estimate_contrasts(lm_kelp_during_zigamma_02, 
#                                 transform = "response") 
# # when transforming to response scale, gives ratio of means (weird?)
# contrasts
# 
# # estimates slopes of numeric predictors for levels of categorical predictors
# # in summary object: estimate = difference in slope of time since end (removal - reference)
# during_slopes <- estimate_slopes(lm_kelp_during_zigamma_02, 
#                                  trend = "time_since_end",
#                                  at = "treatment") %>% 
#   as_tibble() %>% 
#   mutate(exp_dates = "Experimental removal")
# during_emtrends <- emmeans::emtrends(lm_kelp_during_zigamma_02, "treatment", var = "time_since_end")
# 
# estimate_relation(lm_kelp_during_zigamma_02) %>% plot()
# 
# # ⟞ b. recovery period ----------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # model
# # normal model
# # lm_kelp_recovery_lmer <- lmerTest::lmer(
# #   delta_continual ~ time_since_end + (1|site),
# #   data = delta_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass)
# # continuous AR1
# # lm_kelp_recovery_lme_ar1 <- nlme::lme(
# #   delta_continual ~ time_since_end, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass,
# #   correlation = corAR1())
# # ARMA 2
# # lm_kelp_recovery_lme_ar2 <- nlme::lme(
# #   delta_continual ~ time_since_end, random = ~1|site,
# #   data = delta_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass,
# #   correlation = corARMA(p = 2, q = 0))
# # GLS AR1
# # lm_kelp_recovery_gls_ar1 <- nlme::gls(
# #   delta_continual ~ time_since_end, 
# #   data = delta_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass,
# #   correlation = corAR1(form = ~1|site))
# 
# # raw kelp biomass model
# ## site random effect
# lm_kelp_recovery_zigamma_01 <- glmmTMB(
#   kelp_biomass ~ time_since_end*treatment + (1|site),
#   data = continual_long %>% filter(exp_dates == "after"), 
#   family = ziGamma(link = "log"), 
#   ziformula = ~1)
# 
# ## site and year random effect
# lm_kelp_recovery_zigamma_02 <- glmmTMB(
#   kelp_biomass ~ time_since_end*treatment + (1|site) + (1|year),
#   data = continual_long %>% filter(exp_dates == "after"), 
#   family = ziGamma(link = "log"), 
#   ziformula = ~1)
# 
# 
# # check for autocorrelation
# # performance::check_autocorrelation(lm_kelp_recovery_lmer)
# performance::check_autocorrelation(lm_kelp_recovery_zigamma_01)
# performance::check_autocorrelation(lm_kelp_recovery_zigamma_02)
# 
# # diagnostics
# # normal model
# # plot(DHARMa::simulateResiduals(lm_kelp_recovery_lmer))
# # performance::check_model(lm_kelp_recovery_lmer)
# 
# # continuous AR1
# # performance::check_model(lm_kelp_recovery_lme_ar1)
# 
# # GLS AR1
# # qqnorm(lm_kelp_recovery_gls_ar1)
# # plot(fitted(lm_kelp_recovery_gls_ar1), residuals(lm_kelp_recovery_gls_ar1))
# 
# # model checks
# # normal model
# # check_convergence(lm_kelp_recovery_lmer)
# # check_normality(lm_kelp_recovery_lmer)
# # check_heteroscedasticity(lm_kelp_recovery_lmer)
# 
# # plot ACF
# # normal model
# # acf(residuals(lm_kelp_recovery_lmer))
# # pacf(residuals(lm_kelp_recovery_lmer))
# 
# # ARMA 2
# # acf(residuals(lm_kelp_recovery_lme_ar2))
# # pacf(residuals(lm_kelp_recovery_lme_ar2))
# 
# # raw kelp biomass model
# plot(simulateResiduals(lm_kelp_recovery_zigamma_01))
# performance::check_model(lm_kelp_recovery_zigamma_01)
# 
# plot(simulateResiduals(lm_kelp_recovery_zigamma_02))
# performance::check_model(lm_kelp_recovery_zigamma_02)
# 
# # Rsquared
# # MuMIn::r.squaredGLMM(lm_kelp_recovery_lmer)
# # MuMIn::r.squaredGLMM(lm_kelp_recovery_lme_ar1)
# MuMIn::r.squaredGLMM(lm_kelp_recovery_zigamma_01)
# MuMIn::r.squaredGLMM(lm_kelp_recovery_zigamma_02)
# 
# # summary
# # summary(lm_kelp_recovery_lmer)
# # summary(lm_kelp_recovery_lme_ar1)
# # summary(lm_kelp_recovery_gls_ar1)
# # lm_kelp_recovery_summary <- lm_kelp_recovery_lmer %>% 
# #   tbl_regression() %>% 
# #   bold_p(t = 0.05) %>% 
# #   modify_header(
# #     label = " ",
# #     estimate = "**Slope**",
# #     df = "**df**"
# #   ) 
# # lm_kelp_recovery_summary
# 
# summary(lm_kelp_recovery_zigamma_02)
# parameters::model_parameters(lm_kelp_recovery_zigamma_02, exponentiate = TRUE)
# 
# lm_kelp_recovery_zigamma_summary <- lm_kelp_recovery_zigamma_02 %>% 
#   tbl_regression(intercept = TRUE) %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Estimate**"
#   ) %>% 
#   modify_column_indent(
#     columns = label, 
#     rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))
# 
# # filter out zero-inflated component
# lm_kelp_recovery_zigamma_summary$table_body <- lm_kelp_recovery_zigamma_summary$table_body %>% 
#   filter(component != "zi")
# # change labels
# lm_kelp_recovery_zigamma_summary$table_body$label <- c(
#   `(Intercept)` = "(Intercept)",
#   time_since_end = "Time since end",
#   treatmentremoval = "Treatment (removal)",
#   `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
# )
# 
# # final table 
# lm_kelp_recovery_zigamma_summary
# 
# 
# # AIC comparisons
# # MuMIn::AICc(lm_kelp_recovery_lmer, lm_kelp_recovery_lme_ar1, lm_kelp_recovery_lme_ar2, lm_kelp_recovery_gls_ar1) %>% 
# #   arrange(AICc)
# # gls best model when carp is not taken out
# 
# # MuMIn::AICc(lm_kelp_recovery_zigamma_01, lm_kelp_recovery_zigamma_02)
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # # all sites
# # predicted_kelp_after_overall <- ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end", type = "fixed")
# # # predicted line crosses 0 at 4.6
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [4:5 by = 0.01]", type = "fixed")
# # # upper bound crosses 0 at 3.1
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3:4 by = 0.01]", type = "fixed")
# # # lower bound crosses 0 at 6.8
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [6.7:7 by = 0.001]", type = "fixed")
# # 
# # # aque: 4.2 years
# # predicted_kelp_after_aque <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "aque"))
# # # predicted time to recovery: 4.2 years
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [4:5 by = 0.001]", type = "random", condition = c(site = "aque")) 
# # 
# # # napl: 4.1 years
# # predicted_kelp_after_napl <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "napl"))
# # # predicted time to recovery: 4.1 years
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [4:5 by = 0.001]", type = "random", condition = c(site = "napl"))
# # 
# # # mohk: 5.9 years
# # predicted_kelp_after_mohk <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "mohk"))
# # # predicted time to recovery: 5.9 years
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [5.5:6.5 by = 0.01]", type = "random", condition = c(site = "mohk"))
# # 
# # # carp: 3.9 years
# # predicted_kelp_after_carp <- ggpredict(lm_kelp_recovery_lmer, terms = ~ time_since_end, type = "random", condition = c(site = "carp"))
# # # predicted time to recovery: 3.9 years
# # ggpredict(lm_kelp_recovery_lmer, terms = "time_since_end [3.5:4 by = 0.01]", type = "random", condition = c(site = "carp"))
# 
# # raw kelp biomass
# predicted_kelp_after_raw <- ggpredict(lm_kelp_recovery_zigamma_02,
#                                       terms = c("time_since_end [0:6.75 by = 0.25]", "treatment"), type = "fixed")
# # kelp in treatment plot matches kelp in reference plot in about 4 years
# 
# # aque
# predicted_kelp_after_aque <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "aque"))
# # predicted time to recovery: 4 years
# ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.001]", "treatment"), type = "random", condition = c(site = "aque")) %>% plot()
# 
# # napl
# predicted_kelp_after_napl <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "napl"))
# # predicted time to recovery: 4 years
# ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.001]", "treatment"), type = "random", condition = c(site = "napl")) %>% plot()
# 
# # mohk
# predicted_kelp_after_mohk <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "mohk")) 
# # predicted time to recovery: 4 years
# ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.01]", "treatment"), type = "random", condition = c(site = "mohk")) %>% plot()
# 
# # carp
# predicted_kelp_after_carp <- ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end", "treatment"), type = "random", condition = c(site = "carp"))
# # predicted time to recovery: 4 years
# ggpredict(lm_kelp_recovery_zigamma_02, terms = c("time_since_end [3.5:4.5 by = 0.01]", "treatment"), type = "random", condition = c(site = "carp")) %>% plot()
# 
# # ⟞ ⟞ iii. means ----------------------------------------------------------
# 
# # estimates predicted mean values of response for each level in predictor
# recovery_response_means <- estimate_means(lm_kelp_recovery_zigamma_02) %>% 
#   mutate(exp_dates = "after")
# recovery_response_emmeans <- emmeans(lm_kelp_recovery_zigamma_02, ~ treatment | time_since_end, type = "response") %>% plot()
# 
# # calculates differences between levels of categorical predictor
# recovery_contrasts <- modelbased::estimate_contrasts(lm_kelp_recovery_zigamma_02, transform = "response") 
# # when transforming to response scale, gives ratio of means (weird?)
# recovery_contrasts
# 
# # estimates slopes of numeric predictors for levels of categorical predictors
# # in summary object: estimate = difference in slope of time since end removal - reference
# recovery_slopes <- modelbased::estimate_slopes(lm_kelp_recovery_zigamma_02, 
#                                                trend = "time_since_end",
#                                                at = "treatment") %>% 
#   as_tibble() %>% 
#   mutate(exp_dates = "Recovery") %>% 
#   rbind(during_slopes) %>% 
#   mutate(exp_dates = fct_relevel(exp_dates, "Experimental removal", "Recovery"))
# 
# recovery_emtrends <- emmeans::emtrends(lm_kelp_recovery_zigamma_02, "treatment", var = "time_since_end") 
# 
# # ggsave(filename = here::here("figures", "ms-figures", paste("fig-S11_", today(), ".jpg", sep = "")),
# #        effplot,
# #        width = 6, height = 4, dpi = 150)
# 
# # ⟞ c. figures -------------------------------------------------------------
# 
# # ⟞ ⟞ i. overall model predictions -----------------------------------------
# 
# # overall_ms <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_continual, 
# #              aes(x = time_since_end, y = delta_continual, fill = site, shape = site), size = 2, alpha = 0.9) +
# #   scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
# #   scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
# #   # new_scale("color") + 
# #   # overall
# #   geom_line(data = predicted_kelp_after_overall, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_after_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   geom_line(data = predicted_kelp_during_overall, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_during_overall, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
# #   scale_y_continuous(breaks = seq(-1500, 1000, by = 1000), limits = c(-1800, 900)) +
# #   theme_bw() + 
# #   theme(axis.title = element_text(size = 8),
# #         axis.text = element_text(size = 7),
# #         legend.text = element_text(size = 6), 
# #         legend.title = element_text(size = 6),
# #         # plot.margin = margin(0, 0, 0, 0),
# #         legend.position = c(0.1, 0.858),
# #         legend.key.size = unit(0.3, units = "cm")) +
# #   labs(x = "Time since end of removal (years)", 
# #        y = expression("\U0394"~giant~kelp~biomass~"(removal - reference, "~dry~g/m^{"2"}~")"), 
# #        fill = "Site",
# #        shape = "Site")
# 
# 
# # ⟞ ⟞ new model -----------------------------------------------------------
# 
# overall_kelp <- ggplot() +
#   # x at 0 and y at 0 lines
#   geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
#            fill = "grey", alpha = 0.3) +
#   
#   # raw data points
#   geom_point(data = continual_long, aes(x = time_since_end, y = kelp_biomass, color = treatment),
#              alpha = 0.15,
#              size = 0.75,
#              shape = 21) +
#   
#   # prediction lines
#   geom_line(data = predicted_kelp_during_raw, aes(x = x, y = predicted, color = group), linewidth = 1) +
#   geom_line(data = predicted_kelp_after_raw, aes(x = x, y = predicted, color = group), linewidth = 1) +
#   
#   # confidence intervals
#   geom_ribbon(data = predicted_kelp_during_raw, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
#   geom_ribbon(data = predicted_kelp_after_raw, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
#   
#   # colors and shapes
#   scale_color_manual(values = c(reference = reference_col, removal = removal_col),
#                      labels = c("Reference", "Removal")) +
#   scale_linetype_manual(values = c(reference = 2, removal = 1),
#                      labels = c("Reference", "Removal")) +
#   
#   # removal/recovery labels
#   annotate(geom = "text", x = -6.75, y = 1925, label = "Removal", size = 3) +
#   annotate(geom = "text", x = 5.5, y = 1925, label = "Recovery", size = 3) +
#   
#   theme_bw() + 
#   scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
#   coord_cartesian(ylim = c(-10, 2000)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         legend.background = element_blank(),
#         # plot.margin = margin(0, 0, 0, 0),
#         # legend.position = c(0.85, 0.9),
#         legend.position = "none",
#         legend.key.size = unit(0.5, units = "cm"),
#         legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
#         legend.spacing.y = unit(0.1, units = "cm"),
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   guides(color = guide_legend(keyheight = 0.6),
#          shape = guide_legend(keyheight = 0.6),
#          lty = guide_legend(keyheight = 0.6),
#          keyheight = 1) +
#   labs(x = "Time since end of removal (years)", 
#        y = "Biomass (dry g/m\U00B2)",
#        linetype = "Treatment",
#        color = "Treatment",
#        shape = "Treatment",
#        size = "Treatment",
#        title = "(a)")
# 
# overall_kelp
# 
# # napl_raw <- ggplot(data = continual_long %>% filter(site == "napl"),
# #        aes(x = time_since_end, y = kelp_biomass, color = treatment)) +
# #   # x at 0 and y at 0 lines
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   
# #   # raw data points
# #   geom_point(shape = 1, size = 1) +
# #   geom_smooth(method = "lm") +
# #   facet_wrap(~exp_dates)
# # 
# # napl_delta <- ggplot(data = delta_continual %>% filter(site == "napl"),
# #        aes(x = time_since_end, y = delta_continual)) +
# #   geom_point()
# 
# overall_kelp_removal <- ggplot() +
#   # x at 0 and y at 0 lines
#   geom_vline(xintercept = 0, lty = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
#   
#   # raw data points
#   geom_point(data = continual_long %>% filter(treatment == "removal"), aes(x = time_since_end, y = kelp_biomass), shape = 1, size = 1, alpha = 0.4, color = removal_col) +
#   
#   # prediction lines
#   geom_line(data = predicted_kelp_during_raw %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
#   geom_line(data = predicted_kelp_after_raw %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
#   
#   # confidence intervals
#   geom_ribbon(data = predicted_kelp_during_raw %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   geom_ribbon(data = predicted_kelp_after_raw %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   
#   theme_bw() + 
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(limits = c(-10, 2000)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         # plot.margin = margin(0, 0, 0, 0),
#         legend.position = c(0.88, 0.73),
#         legend.key.size = unit(0.5, units = "cm"),
#         legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
#         legend.spacing.y = unit(0.1, units = "cm"),
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   guides(color = guide_legend(keyheight = 0.6),
#          shape = guide_legend(keyheight = 0.6),
#          lty = guide_legend(keyheight = 0.6),
#          keyheight = 1) +
#   labs(x = "Time since end of removal (years)", 
#        y = "Giant kelp biomass (dry g/m\U00B2)",
#        title = "(a) Removal")
# 
# overall_kelp_removal
# 
# overall_kelp_reference <- ggplot() +
#   # x at 0 and y at 0 lines
#   geom_vline(xintercept = 0, lty = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
#   
#   # raw data points
#   geom_point(data = continual_long %>% filter(treatment == "reference"), aes(x = time_since_end, y = kelp_biomass), shape = 1, size = 1, alpha = 0.4, color = reference_col) +
#   
#   # prediction lines
#   geom_line(data = predicted_kelp_during_raw %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
#   geom_line(data = predicted_kelp_after_raw %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
#   
#   # confidence intervals
#   geom_ribbon(data = predicted_kelp_during_raw %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   geom_ribbon(data = predicted_kelp_after_raw %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   
#   theme_bw() + 
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(limits = c(-10, 2000)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         # plot.margin = margin(0, 0, 0, 0),
#         legend.position = c(0.88, 0.73),
#         legend.key.size = unit(0.5, units = "cm"),
#         legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
#         legend.spacing.y = unit(0.1, units = "cm"),
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   guides(color = guide_legend(keyheight = 0.6),
#          shape = guide_legend(keyheight = 0.6),
#          lty = guide_legend(keyheight = 0.6),
#          keyheight = 1) +
#   labs(x = "Time since end of removal (years)", 
#        y = "Giant kelp biomass (dry g/m\U00B2)",
#        title = "(b) Reference")
# 
# overall_kelp_reference
# 
# # ⟞ ⟞ ii. site level predictions -------------------------------------------
# 
# # aque <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_continual %>% filter(site == "aque"), aes(x = time_since_end, y = delta_continual), shape = aque_shape, fill = aque_col, size = 2, alpha = 0.9) +
# #   # during
# #   geom_line(data = predicted_kelp_during_aque, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_during_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   # after
# #   geom_line(data = predicted_kelp_after_aque, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_after_aque, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
# #   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1250, 1100), expand = c(0, 0)) +
# #   # geom_text(aes(x = -6.6, y = 600), label = "aque", size = 8) +
# #   theme_bw() + 
# #   theme(axis.title = element_text(size = 8),
# #         plot.title = element_text(size = 8),
# #         axis.text = element_text(size = 7),
# #         plot.title.position = "plot") +
# #   labs(x = "Time since end of removal", 
# #        y = "\U0394 giant kelp biomass",
# #        title = paste("(a) ", aque_full, sep = ""))
# # 
# # napl <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_continual %>% filter(site == "napl"), aes(x = time_since_end, y = delta_continual), shape = napl_shape, fill = napl_col, size = 2, alpha = 0.9) +
# #   # during
# #   geom_line(data = predicted_kelp_during_napl, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_during_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   # after
# #   geom_line(data = predicted_kelp_after_napl, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_after_napl, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
# #   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1250, 1100), expand = c(0, 0)) +
# #   # geom_text(aes(x = -6.8, y = 700), label = "napl", size = 8) +
# #   theme_bw() + 
# #   theme(axis.title = element_text(size = 8),
# #         plot.title = element_text(size = 8),
# #         axis.text = element_text(size = 7),
# #         plot.title.position = "plot") +
# #   labs(x = "Time since end of removal", 
# #        y = "\U0394 giant kelp biomass",
# #        title = paste("(b) ", napl_full, sep = ""))
# # 
# # mohk <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_continual %>% filter(site == "mohk"), aes(x = time_since_end, y = delta_continual), shape = mohk_shape, fill = mohk_col, size = 2, alpha = 0.9) +
# #   # during
# #   geom_line(data = predicted_kelp_during_mohk, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_during_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   # after
# #   geom_line(data = predicted_kelp_after_mohk, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_after_mohk, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
# #   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1800, 1100), expand = c(0, 0)) +
# #   # geom_text(aes(x = -6.6, y = 900), label = "mohk", size = 8) +
# #   theme_bw() + 
# #   theme(axis.title = element_text(size = 8),
# #         plot.title = element_text(size = 8),
# #         axis.text = element_text(size = 7),
# #         plot.title.position = "plot") +
# #   labs(x = "Time since end of removal", 
# #        y = "\U0394 giant kelp biomass",
# #        title = paste("(c) ", mohk_full, sep = ""))
# # 
# # carp <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_continual %>% filter(site == "carp"), aes(x = time_since_end, y = delta_continual), shape = carp_shape, fill = carp_col, size = 2, alpha = 0.9) +
# #   # during
# #   geom_line(data = predicted_kelp_during_carp, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_during_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   # after
# #   geom_line(data = predicted_kelp_after_carp, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_kelp_after_carp, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 5, by = 1), minor_breaks = NULL) +
# #   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1250, 1100), expand = c(0, 0)) +
# #   # geom_text(aes(x = -6.8, y = 1000), label = "carp", size = 8) +
# #   theme_bw() +
# #   theme(axis.title = element_text(size = 8),
# #         plot.title = element_text(size = 8),
# #         axis.text = element_text(size = 7),
# #         plot.title.position = "plot") +
# #   labs(x = "Time since end of removal", 
# #        y = "\U0394 giant kelp biomass",
# #        title = paste("(d) ", carp_full, sep = ""))
# # 
# # plots_together_sites <- (aque + napl) / (mohk + carp)
# # plots_together_sites
# 
# 
# # ⟞ ⟞ iii. delta from predictions -----------------------------------------
# 
# # data frame of predictions
# delta_predictions_during <- predicted_kelp_during_raw %>% 
#   as.data.frame() %>% 
#   select(x, group, predicted) %>% 
#   pivot_wider(names_from = group, values_from = predicted) %>% 
#   mutate(delta = removal - reference) %>% 
#   mutate(exp_dates = "during")
# 
# delta_predictions_after <- predicted_kelp_after_raw %>% 
#   as.data.frame() %>% 
#   select(x, group, predicted) %>% 
#   pivot_wider(names_from = group, values_from = predicted) %>% 
#   mutate(delta = removal - reference) %>% 
#   mutate(exp_dates = "after")
# 
# delta_kelp_predictions <- ggplot() +
#   geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
#            fill = "grey", alpha = 0.3) +
#   geom_point(data = delta_continual,
#              aes(x = time_since_end, y = delta_continual), 
#              shape = 2, 
#              alpha = 0.15,
#              size = 0.75) +
# 
#   # overall
#   geom_line(data = delta_predictions_during, aes(x = x, y = delta), linewidth = 1) +
#   geom_line(data = delta_predictions_after, aes(x = x, y = delta), linewidth = 1) +
# 
#   scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
#   scale_y_continuous(breaks = seq(-1500, 1000, by = 500), limits = c(-1800, 1000)) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         # plot.margin = margin(0, 0, 0, 0),
#         legend.position = "none",
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   labs(x = "Time since end of removal (years)", 
#        y = "\U0394 biomass\n(removal \U2212 reference, dry g/m\U00B2)",
#        fill = "Site",
#        shape = "Site",
#        title = "(b) Removal \U2212 reference")
# 
# delta_kelp_predictions
# 
# fig1_new <- overall_kelp / delta_kelp_predictions
# # fig1_new_v2 <- overall_kelp_removal / overall_kelp_reference / delta_kelp_predictions
# 
# 
# # ⟞ ⟞ iv. effect plots ----------------------------------------------------
# 
# kelp_means_df <- rbind(during_response_means, recovery_response_means) %>% 
#   as_tibble() %>% 
#   mutate(exp_dates = fct_relevel(exp_dates, "during", "after"))
# 
# # mean predicted responses: what are the differences in mean giant kelp biomass between treatments during each time period?
# kelp_means_plot <- continual_long %>% 
#   # filter(exp_dates == "during") %>% 
#   ggplot(aes(x = treatment, y = kelp_biomass, color = treatment)) +
#   facet_wrap(~exp_dates, 
#              labeller = labeller(exp_dates = c("during" = "Experimental removal",
#                                                "after" = "Recovery"))) +
#   geom_point(position = position_jitter(width = 0.1, seed = 666),
#              alpha = 0.3, shape = 21) +
#   geom_pointrange(data = kelp_means_df, 
#                   aes(x = treatment, y = Mean, ymin = CI_low, ymax = CI_high),
#                   size = 0.5, linewidth = 0.5) +
#   scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
#   scale_fill_manual(values = c(reference = reference_col, removal = removal_col)) +
#   scale_x_discrete(labels = c("reference" = "Reference", "removal" = "Removal")) +
#   scale_y_continuous(limits = c(-5, 2000), expand = c(0.01, 0.01)) +
#   labs(x = "Treatment", 
#        y = "Giant kelp biomass (dry g/m\U00B2)") +
#   theme_classic() +
#   theme(legend.position = "none",
#         strip.placement = "outer",
#         strip.text = element_text(hjust = 0.5),
#         strip.background = element_blank()) 
# 
# kelp_means_plot
# 
# kelp_means_nodata_plot <- ggplot(data = kelp_means_df,
#                             aes(x = treatment, y = kelp_biomass, color = treatment)) +
#   facet_wrap(~exp_dates, 
#              labeller = labeller(exp_dates = c("during" = "Experimental removal",
#                                                "after" = "Recovery"))) +
#   geom_pointrange(aes(x = treatment, y = Mean, ymin = CI_low, ymax = CI_high),
#                   size = 0.5, linewidth = 0.5) +
#   scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
#   scale_fill_manual(values = c(reference = reference_col, removal = removal_col)) +
#   scale_x_discrete(labels = c("reference" = "Reference", "removal" = "Removal")) +
#   scale_y_continuous(limits = c(-5, 550), 
#                      breaks = seq(from = 0, to = 550, by = 125),
#                      expand = c(0.01, 0.01)) +
#   labs(x = "Treatment", 
#        y = "Predicted giant kelp biomass \n(dry g/m\U00B2)") +
#   theme_classic() +
#   theme(legend.position = "none",
#         strip.placement = "outer",
#         strip.text = element_text(hjust = 0.5),
#         strip.background = element_blank()) 
# kelp_means_nodata_plot
# 
# # mean predicted slopes: what are the differences in rates of giant kelp biomass change?
# effplot <- ggplot(data = recovery_slopes, aes(x = treatment, y = Coefficient)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_pointrange(aes(ymin = CI_low, ymax = CI_high)) +
#   labs(x = "Treatment", 
#        y = "Effect of time since end (years)") +
#   theme_bw() +
#   raw_biomass_plot_theme +
#   facet_wrap(~exp_dates)
# effplot
# 
# # ⟞ ⟞ v. conditional means ------------------------------------------------
# 
# predicted_kelp_after_0_raw <- ggpredict(lm_kelp_recovery_zigamma_02,
#                                       terms = c("time_since_end [0]", "treatment"), type = "fixed")
# 
# predicted_kelp_after_both_raw <- ggpredict(lm_kelp_recovery_zigamma_02,
#                                         terms = c("time_since_end [4]", "treatment"), 
#                                         type = "fixed") %>% 
#   rbind(predicted_kelp_after_0_raw) %>% 
#   rename("treatment" = group, 
#          "kelp_biomass" = predicted,
#          "time_since_end" = x)
# 
# means <- continual_long %>% 
#   filter(time_since_end == 0 | time_since_end == 4) %>% 
#   ggplot(aes(x = treatment, y = kelp_biomass, color = treatment)) +
#   geom_point(position = position_jitter(width = 0.1, seed = 1),
#              shape = 21, alpha = 0.8, size = 1) +
#   geom_pointrange(data = predicted_kelp_after_both_raw,
#                   aes(ymin = conf.low, ymax = conf.high)) +
#   scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
#   scale_x_discrete(labels = c(reference = "Reference", removal = "Removal")) +
#   labs(x = "Treatment",
#        y = "Giant kelp biomass (dry g/m\U00B2)") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         strip.text = element_text(hjust = 0, size = 10),
#         strip.background = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   #geom_pointrange(data = predicted_kelp_after_0, ) +
#   facet_wrap(~time_since_end, labeller = labeller(time_since_end = c("0" = "(a) Time since end = 0", "4" = "(b) Time since end = 4")))
# means
# 
# 
# ##########################################################################-
# # 4. variation plots ------------------------------------------------------
# ##########################################################################-
# 
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
# # 
# var_df <- rbind(var_during_df, var_after_df) %>%
#   rename(treatment = x) %>%
#   mutate(exp_dates = fct_relevel(exp_dates, "during", "after"))
# # 
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
# # 
# ggplot(data = var_df, aes(x = treatment, y = predicted)) +
#   geom_pointrange(aes(x = treatment, y = predicted, ymin = conf.low, ymax = conf.high)) +
#   facet_wrap(~exp_dates) +
#   labs(x = "Treatment",
#        y = "Variation") +
#   theme_bw() +
#   theme(panel.grid = element_blank())
# 
# variation_site <- continual_long %>% 
#   select(exp_dates, treatment, site, kelp_biomass) %>% 
#   group_by(exp_dates, treatment, site) %>% 
#   summarize(mean = mean(kelp_biomass),
#             variation = sd(kelp_biomass)/mean(kelp_biomass)) %>% 
#   ungroup()
# 
# variation_plot <- ggplot(variation_site, aes(x = treatment, y = variation, color = treatment)) +
#   geom_point(position = position_jitter(width = 0.1, seed = 666),
#              alpha = 0.8, shape = 21) +
#   stat_summary(fun.data = mean_se, geom = "pointrange") +
#   scale_color_manual(values = c(reference = reference_col, removal = removal_col)) +
#   scale_x_discrete(labels = c(reference = "Reference", removal = "Removal")) +
#   facet_wrap(~exp_dates, labeller = labeller(exp_dates = c("during" = "(a) Experimental removal", "after" = "(b) Recovery period"))) +
#   labs(x = "Treatment",
#        y = "Coefficient of variation") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         strip.text = element_text(hjust = 0, size = 10),
#         strip.background = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "none") 
# 
# variation_plot
# 
# ##########################################################################-
# # 5. manuscript tables ----------------------------------------------------
# ##########################################################################-
# 
# # lm_kelp_tables <- tbl_merge(tbls = list(lm_kelp_during_summary, lm_kelp_recovery_summary), 
# #                             tab_spanner = c("**Removal**", "**Recovery**")) 
# # this table is compiled with others in the `02a-community_recovery.R` script
# 
# lm_kelp_zigamma_tables <- tbl_merge(tbls = list(lm_kelp_during_zigamma_summary, lm_kelp_recovery_zigamma_summary),
#                                     tab_spanner = c("**Experimental removal**", "**Recovery**"))
# 
# ##########################################################################-
# # 6. manuscript figures ---------------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. delta kelp through time -------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-1_", today(), ".jpg", sep = "")),
# #        plot = overall_ms,
# #        height = 8, width = 14, units = "cm",
# #        dpi = 400)
# 
# # ⟞ b. raw kelp biomass through time --------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S1_", today(), ".jpg", sep = "")),
# #        plot = continual_sites_raw,
# #        height = 12, width = 8, units = "cm",
# #        dpi = 300)
# 
# # ⟞ c. recovery time vs biomass -------------------------------------------
# 
# # s4_panels <- rec_time_plot + delta_vs_biomass +
# #   plot_layout(widths = c(2, 2.1)) +
# #   plot_annotation(tag_levels = "a", tag_suffix = ")") &
# #   theme(plot.tag = element_text(size = 12))
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S4_", today(), ".jpg", sep = "")),
# #        plot = rec_time_plot,
# #        height = 8, width = 12, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S4_panels_", today(), ".jpg", sep = "")),
# #        plot = s4_panels,
# #        height = 8, width = 16, units = "cm",
# #        dpi = 400)
# 
# # ⟞ d. predictions by site ------------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S2_", today(), ".jpg", sep = "")),
# #        plot = plots_together_sites,
# #        height = 10, width = 16, units = "cm",
# #        dpi = 300)
# 
# 
# # ⟞ e. new model ----------------------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-1_new-model_", today(), ".jpg", sep = "")),
# #        plot = fig1_new,
# #        height = 17, width = 13, units = "cm",
# #        dpi = 300)
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-1_new-model_v2_", today(), ".jpg", sep = "")),
# #        plot = fig1_new_v2,
# #        height = 17, width = 13, units = "cm",
# #        dpi = 300)
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-1_new-model_removal", today(), ".jpg", sep = "")),
# #        plot = overall_kelp_removal,
# #        height = 8, width = 14, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-1_new-model_reference", today(), ".jpg", sep = "")),
# #        plot = overall_kelp_reference,
# #        height = 8, width = 14, units = "cm",
# #        dpi = 300)
# 
# # Ecology max figure size: 18 x 24
# 
# # ⟞ f. means plots --------------------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("kelp_means_plot_", today(), ".jpg", sep = "")),
# #        plot = kelp_means_plot,
# #        height = 7, width = 14, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("kelp_means_nodata_plot_", today(), ".jpg", sep = "")),
# #        plot = kelp_means_nodata_plot,
# #        height = 7, width = 14, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("kelp-means-recovery-plot_", today(), ".jpg", sep = "")),
# #        plot = means,
# #        height = 7, width = 14, units = "cm",
# #        dpi = 300)
# 
# # ⟞ g. kelp density and fronds --------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("kelp-density-timeseries_", today(), ".jpg", sep = "")),
# #        plot = density_timeseries,
# #        height = 9, width = 16, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("kelp-frond-timeseries_", today(), ".jpg", sep = "")),
# #        plot = fronds_timeseries,
# #        height = 9, width = 16, units = "cm",
# #        dpi = 300)
# 
# # ⟞ h. variation plots --------------------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("cov_plot_", today(), ".jpg", sep = "")),
# #        plot = variation_plot,
# #        height = 7, width = 14, units = "cm",
# #        dpi = 300)
# 
# 


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # -------------------------- OLD CODE BELOW HERE --------------------------
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # ⟞ a. during removal -----------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # model
# # lm_algae_during_lmer <- lmer(
# #   delta_continual_algae ~ time_since_end + (1|site), 
# #   data = delta_algae_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass)
# 
# # lm_raw_algae_during_zigamma_01 <- glmmTMB(
# #   algae_biomass ~ time_since_end*treatment + (1|site), 
# #   data = algae_continual_long %>% filter(exp_dates == "during"), 
# #   na.action = na.pass,
# #   family = ziGamma(link = "log"),
# #   ziformula = ~1)
# 
# lm_raw_algae_during_zigamma_02 <- glmmTMB(
#   biomass ~ time_since_end*treatment + (1|site) + (1|year), 
#   data = algae_continual_long %>% filter(exp_dates == "during"), 
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1)
# 
# 
# # diagnostics
# # plot(simulateResiduals(lm_algae_during_lmer))
# # check_model(lm_algae_during_lmer)
# 
# # plot(simulateResiduals(lm_raw_algae_during_zigamma_01))
# # check_convergence(lm_raw_algae_during_zigamma_01)
# 
# plot(simulateResiduals(lm_raw_algae_during_zigamma_02))
# 
# # Rsquared
# # r.squaredGLMM(lm_algae_during_lmer)
# # r.squaredGLMM(lm_raw_algae_during_zigamma_01)
# r.squaredGLMM(lm_raw_algae_during_zigamma_02)
# 
# # summary
# # summary(lm_algae_during_lmer)
# # lm_algae_during_summary <- lm_algae_during_lmer %>% 
# #   tbl_regression() %>% 
# #   bold_p(t = 0.05) %>% 
# #   modify_header(
# #     label = " ",
# #     estimate = "**Slope**",
# #     df = "**df**"
# #   ) 
# # lm_algae_during_summary
# 
# lm_raw_algae_during_zigamma_summary <- lm_raw_algae_during_zigamma_02 %>% 
#   tbl_regression(intercept = TRUE) %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Estimate**"
#   ) %>% 
#   modify_column_indent(
#     columns = label, 
#     rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))
# 
# # filter out zero-inflated component
# lm_raw_algae_during_zigamma_summary$table_body <- lm_raw_algae_during_zigamma_summary$table_body %>% 
#   filter(component != "zi") %>% 
#   drop_na(term)
# # change labels
# lm_raw_algae_during_zigamma_summary$table_body$label <- c(
#   `(Intercept)` = "(Intercept)",
#   time_since_end = "Time since end",
#   treatmentremoval = "Treatment (removal)",
#   `time_since_end:treatmentremoval` = "Time since end × treatment (removal)" 
# )
# 
# # final table 
# lm_raw_algae_during_zigamma_summary
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # predicted_algae_during <- ggpredict(lm_algae_during_lmer, terms = ~ time_since_end, type = "fixed")
# 
# predicted_raw_algae_during <- ggpredict(lm_raw_algae_during_zigamma_02, terms = c("time_since_end", "treatment"), type = "fixed")
# 
# # ⟞ b. recovery period ----------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # model
# # lm_algae_recovery_lmer <- lmer(
# #   delta_continual_algae ~ time_since_end + (1|site), 
# #   data = delta_algae_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass)
# 
# # lm_raw_algae_recovery_zigamma_01 <- glmmTMB(
# #   algae_biomass ~ time_since_end*treatment + (1|site),
# #   data = algae_continual_long %>% filter(exp_dates == "after"),
# #   na.action = na.pass,
# #   family = ziGamma(link = "log"),
# #   ziformula = ~1)
# 
# lm_raw_algae_recovery_zigamma_02 <- glmmTMB(
#   algae_biomass ~ time_since_end*treatment + (1|site) + (1|year), 
#   data = algae_continual_long %>% filter(exp_dates == "after"), 
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1)
# 
# # diagnostics
# # plot(simulateResiduals(lm_algae_recovery_lmer))
# # check_model(lm_algae_recovery_lmer)
# 
# # plot(simulateResiduals(lm_raw_algae_recovery_zigamma_01))
# 
# plot(simulateResiduals(lm_raw_algae_recovery_zigamma_02))
# 
# # Rsquared
# # r.squaredGLMM(lm_algae_recovery_lmer)
# # r.squaredGLMM(lm_raw_algae_recovery_zigamma_01)
# r.squaredGLMM(lm_raw_algae_recovery_zigamma_02)
# 
# # summary
# # summary(lm_algae_recovery_lmer)
# # lm_algae_recovery_summary <- lm_algae_recovery_lmer %>% 
# #   tbl_regression() %>% 
# #   bold_p(t = 0.05) %>% 
# #   modify_header(
# #     label = " ",
# #     estimate = "**Slope**",
# #     df = "**df**"
# #   ) 
# # lm_algae_recovery_summary
# 
# lm_raw_algae_recovery_zigamma_summary <- lm_raw_algae_recovery_zigamma_02 %>% 
#   tbl_regression(intercept = TRUE) %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Estimate**"
#   ) %>% 
#   modify_column_indent(
#     columns = label, 
#     rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))
# 
# # filter out zero-inflated component
# lm_raw_algae_recovery_zigamma_summary$table_body <- lm_raw_algae_recovery_zigamma_summary$table_body %>% 
#   filter(component != "zi")
# # change labels
# lm_raw_algae_recovery_zigamma_summary$table_body$label <- c(
#   `(Intercept)` = "(Intercept)",
#   time_since_end = "Time since end",
#   treatmentremoval = "Treatment (removal)",
#   `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
# )
# 
# # final table 
# lm_raw_algae_recovery_zigamma_summary
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # predicted_algae_recovery <- ggpredict(lm_algae_recovery_lmer, terms = ~time_since_end, type = "fixed")
# 
# predicted_raw_algae_recovery <- ggpredict(lm_raw_algae_recovery_zigamma_02, terms = c("time_since_end[0:6.75, by = 0.25]", "treatment"), type = "fixed")
# 
# # ⟞ c. figure ------------------------------------------------------------
# 
# # algae_time <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_algae_continual, 
# #              aes(x = time_since_end, y = delta_continual_algae, fill = site, shape = site), size = 2, alpha = 0.9) +
# #   scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
# #   scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
# #   # new_scale("color") + 
# #   # overall
# #   geom_line(data = predicted_algae_recovery, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_algae_recovery, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   geom_line(data = predicted_algae_during, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_algae_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
# #   # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
# #   delta_timeseries_theme("algae") +
# #   labs(x = "Time since end of removal (years)",
# #        y = "\U0394 biomass \n (treatment - control)",
# #        title = "(b)",
# #        fill = "Site", shape = "Site")
# # 
# # algae_time
# 
# 
# # ⟞ ⟞ new model -----------------------------------------------------------
# 
# 
# 
# ##########################################################################-
# # 4. epi. invert linear model ---------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. during removal -----------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # model
# # lm_epi_during_lmer <- lmer(
# #   delta_continual_epi ~ time_since_end + (1|site), 
# #   data = delta_epi_continual %>% filter(exp_dates == "during"), 
# #   na.action = na.pass)
# 
# # lm_raw_epi_during_lmer <- lmer(
# #   epi_biomass ~ time_since_end*treatment + (1|site), 
# #   data = epi_continual_long %>% filter(exp_dates == "during"), 
# #   na.action = na.pass)
# 
# # lm_raw_epi_during_zigamma_01 <- glmmTMB(
# #   epi_biomass ~ time_since_end*treatment + (1|site),
# #   data = epi_continual_long %>% filter(exp_dates == "during"),
# #   na.action = na.pass,
# #   family = ziGamma(link = "log"),
# #   ziformula = ~1
# # )
# 
# lm_raw_epi_during_zigamma_02 <- glmmTMB(
#   epi_biomass ~ time_since_end*treatment + (1|site) + (1|year),
#   data = epi_continual_long %>% filter(exp_dates == "during"),
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1
# )
# 
# # diagnostics
# # plot(simulateResiduals(lm_epi_during_lmer))
# # check_model(lm_epi_during_lmer)
# 
# # plot(simulateResiduals(lm_raw_epi_during_lmer))
# # check_model(lm_raw_epi_during_lmer)
# 
# # plot(simulateResiduals(lm_raw_epi_during_zigamma_01))
# 
# plot(simulateResiduals(lm_raw_epi_during_zigamma_02))
# 
# # R2
# # r.squaredGLMM(lm_epi_during_lmer)
# r.squaredGLMM(lm_raw_epi_during_zigamma_02)
# 
# # summary
# # summary(lm_epi_during_lmer) 
# # lm_epi_during_summary <- lm_epi_during_lmer %>% 
# #   tbl_regression() %>% 
# #   bold_p(t = 0.05) %>% 
# #   modify_header(
# #     label = " ",
# #     estimate = "**Slope**",
# #     df = "**df**"
# #   ) 
# # lm_epi_during_summary
# 
# lm_raw_epi_during_zigamma_summary <- lm_raw_epi_during_zigamma_02 %>% 
#   tbl_regression(intercept = TRUE) %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Estimate**"
#   ) %>% 
#   modify_column_indent(
#     columns = label, 
#     rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))
# 
# # filter out zero-inflated component
# lm_raw_epi_during_zigamma_summary$table_body <- lm_raw_epi_during_zigamma_summary$table_body %>% 
#   filter(component != "zi")
# # change labels
# lm_raw_epi_during_zigamma_summary$table_body$label <- c(
#   `(Intercept)` = "(Intercept)",
#   time_since_end = "Time since end",
#   treatmentremoval = "Treatment (removal)",
#   `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
# )
# 
# # final table 
# lm_raw_epi_during_zigamma_summary
# summary(lm_raw_epi_during_zigamma_02)
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # predicted_epi_during <- ggpredict(lm_epi_during_lmer, terms = ~ time_since_end, type = "fixed")
# 
# predicted_raw_epi_during <- ggpredict(lm_raw_epi_during_zigamma_02, terms = c("time_since_end", "treatment"), type = "fixed")
# 
# # ⟞ b. recovery period ----------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # model
# # lm_epi_recovery_lmer <- lmer(
# #   delta_continual_epi ~ time_since_end + (1|site), 
# #   data = delta_epi_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass)
# 
# # lm_raw_epi_recovery_lmer <- lmer(
# #   epi_biomass ~ time_since_end*treatment + (1|site), 
# #   data = delta_epi_continual %>% filter(exp_dates == "after"), 
# #   na.action = na.pass)
# 
# # lm_raw_epi_recovery_zigamma_01 <- glmmTMB(
# #   epi_biomass ~ time_since_end*treatment + (1|site),
# #   data = epi_continual_long %>% filter(exp_dates == "after"),
# #   na.action = na.pass,
# #   family = ziGamma(link = "log"),
# #   ziformula = ~1
# # )
# 
# lm_raw_epi_recovery_zigamma_02 <- glmmTMB(
#   epi_biomass ~ time_since_end*treatment + (1|site) + (1|year),
#   data = epi_continual_long %>% filter(exp_dates == "after"),
#   na.action = na.pass,
#   family = ziGamma(link = "log"),
#   ziformula = ~1
# )
# 
# 
# # diagnostics
# # plot(simulateResiduals(lm_epi_recovery_lmer))
# # check_model(lm_epi_recovery_lmer)
# 
# # plot(simulateResiduals(lm_raw_epi_recovery_zigamma_01))
# plot(simulateResiduals(lm_raw_epi_recovery_zigamma_02))
# 
# # R2
# MuMIn::r.squaredGLMM(lm_raw_epi_recovery_zigamma_02)
# 
# # summary table
# # summary(lm_epi_recovery_lmer)
# # lm_epi_recovery_summary <- lm_epi_recovery_lmer %>% 
# #   tbl_regression() %>% 
# #   bold_p(t = 0.05) %>% 
# #   modify_header(
# #     label = " ",
# #     estimate = "**Slope**",
# #     df = "**df**"
# #   ) 
# # lm_epi_recovery_summary
# 
# lm_raw_epi_recovery_zigamma_summary <- lm_raw_epi_recovery_zigamma_02 %>% 
#   tbl_regression(intercept = TRUE) %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Estimate**"
#   ) %>% 
#   modify_column_indent(
#     columns = label, 
#     rows = variable %in% c("(Intercept)", "treatment", "time_since_end", "time_since_end:treatment"))
# 
# # filter out zero-inflated component
# lm_raw_epi_recovery_zigamma_summary$table_body <- lm_raw_epi_recovery_zigamma_summary$table_body %>% 
#   filter(component != "zi")
# # change labels
# lm_raw_epi_recovery_zigamma_summary$table_body$label <- c(
#   `(Intercept)` = "(Intercept)",
#   time_since_end = "Time since end",
#   treatmentremoval = "Treatment (removal)",
#   `time_since_end:treatmentremoval` = "Time since end * treatment (removal)" 
# )
# 
# # final table 
# lm_raw_epi_recovery_zigamma_summary
# summary(lm_raw_epi_recovery_zigamma_02)
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # predicted_epi_recovery <- ggpredict(lm_epi_recovery_lmer, terms = ~ time_since_end, type = "fixed")
# 
# predicted_raw_epi_recovery <- ggpredict(lm_raw_epi_recovery_zigamma_02, terms = c("time_since_end[0:6.75, by = 0.25]", "treatment"), type = "fixed")
# 
# # ⟞ c. figure ------------------------------------------------------------
# 
# # epi_time <- ggplot() +
# #   geom_vline(xintercept = 0, lty = 2) +
# #   geom_hline(yintercept = 0, lty = 2) +
# #   geom_point(data = delta_epi_continual, 
# #              aes(x = time_since_end, y = delta_continual_epi, fill = site, shape = site), 
# #              size = 2, alpha = 0.9) +
# #   scale_shape_manual(values = shape_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
# #   scale_fill_manual(values = color_palette_site, labels = c("aque" = aque_full, "napl" = napl_full, "mohk" = mohk_full, carp = carp_full)) +
# #   # overall
# #   # geom_line(data = predicted_epi_after, aes(x = x, y = predicted), size = 2, alpha = 0.7) +
# #   # geom_ribbon(data = predicted_epi_after, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   geom_line(data = predicted_epi_during, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_epi_during, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   geom_line(data = predicted_epi_recovery, aes(x = x, y = predicted), linewidth = 1, alpha = 0.7) +
# #   geom_ribbon(data = predicted_epi_recovery, aes(x = x, ymax = conf.high, ymin = conf.low), alpha = 0.2) +
# #   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
# #   # scale_y_continuous(breaks = seq(-250, 750, by = 250), limits = c(-250, 750)) +
# #   delta_timeseries_theme("epi") +
# #   labs(x = "Time since end of removal (years)", 
# #        y = "\U0394 biomass \n (treatment - control)",
# #        subtitle = "(d)")
# # epi_time
# 
#  # ⟞ ⟞ new model -----------------------------------------------------------
# 
# overall_epi_predictions <- ggplot() +
#   geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
#            fill = "grey", alpha = 0.3) +
#   
#   # raw data
#   geom_point(data = epi_continual_long, 
#              aes(x = time_since_end, y = epi_biomass, color = treatment), 
#              shape = 21,
#              alpha = 0.15,
#              size = 0.75) +
#   
#   # model predictions
#   geom_line(data = predicted_raw_epi_recovery, aes(x = x, y = predicted, color = group), linewidth = 1) +
#   geom_ribbon(data = predicted_raw_epi_recovery, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
#   geom_line(data = predicted_raw_epi_during, aes(x = x, y = predicted, color = group), linewidth = 1) +
#   geom_ribbon(data = predicted_raw_epi_during, aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.05) +
#   
#   # colors and shapes
#   scale_color_manual(values = c(reference = reference_col, removal = removal_col),
#                      labels = c("Reference", "Removal")) +
#   # scale_fill_manual(values = c(reference = "#6D5A1800", removal = "#CC754066"),
#   #                    labels = c("Reference", "Removal")) +
#   scale_linetype_manual(values = c(reference = 2, removal = 1),
#                         labels = c("Reference", "Removal")) +
#   scale_shape_manual(values = c(reference = 1, removal = 16),
#                      labels = c("Reference", "Removal")) +
#   scale_size_manual(values = c(reference = 1, removal = 1.3),
#                     labels = c("Reference", "Removal")) +
#   
#   # removal/recovery labels
#   annotate(geom = "text", x = -6.75, y = 150, label = "Removal", size = 3) +
#   annotate(geom = "text", x = 5.5, y = 150, label = "Recovery", size = 3) +
#   
#   # theming
#   theme_bw() + 
#   scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
#   coord_cartesian(ylim = c(5, 155)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.position = "none",
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   guides(color = guide_legend(keyheight = 0.6),
#          shape = guide_legend(keyheight = 0.6),
#          lty = guide_legend(keyheight = 0.6),
#          keyheight = 1) +
#   labs(x = "Time since end of removal (years)", 
#        y = "Biomass (dry g/m\U00B2)", 
#        title = "(e)")
# 
# overall_epi_predictions
# 
# raw_epi_removal <- ggplot() +
#   # x at 0 and y at 0 lines
#   geom_vline(xintercept = 0, lty = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
#   
#   # raw data points
#   geom_point(data = epi_continual_long %>% filter(treatment == "removal"), 
#              aes(x = time_since_end, y = epi_biomass), 
#              shape = 1, size = 1, alpha = 0.4, color = removal_col) +
#   
#   # prediction lines
#   geom_line(data = predicted_raw_epi_during %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
#   geom_line(data = predicted_raw_epi_recovery %>% filter(group == "removal"), aes(x = x, y = predicted), linewidth = 1, color = removal_col) +
#   
#   # confidence intervals
#   geom_ribbon(data = predicted_raw_epi_during %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   geom_ribbon(data = predicted_raw_epi_recovery %>% filter(group == "removal"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   
#   theme_bw() + 
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   coord_cartesian(ylim = c(4, 155)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         # plot.margin = margin(0, 0, 0, 0),
#         legend.position = c(0.88, 0.73),
#         legend.key.size = unit(0.5, units = "cm"),
#         legend.box.margin = margin(0.01, 0.01, 0.01, 0.01),
#         legend.spacing.y = unit(0.1, units = "cm"),
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   guides(color = guide_legend(keyheight = 0.6),
#          shape = guide_legend(keyheight = 0.6),
#          lty = guide_legend(keyheight = 0.6),
#          keyheight = 1) +
#   labs(x = "Time since end of removal (years)", 
#        y = "Biomass (dry g/m\U00B2)",
#        title = "(d) Removal")
# raw_epi_removal
# 
# raw_epi_reference <- ggplot() +
#   # x at 0 and y at 0 lines
#   geom_vline(xintercept = 0, lty = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
#   
#   # raw data points
#   geom_point(data = epi_continual_long %>% filter(treatment == "reference"), aes(x = time_since_end, y = epi_biomass), shape = 1, size = 1, alpha = 0.4, color = reference_col) +
#   
#   # prediction lines
#   geom_line(data = predicted_raw_epi_during %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
#   geom_line(data = predicted_raw_epi_recovery %>% filter(group == "reference"), aes(x = x, y = predicted), linewidth = 1, linetype = 2, color = reference_col) +
#   
#   # confidence intervals
#   geom_ribbon(data = predicted_raw_epi_during %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   geom_ribbon(data = predicted_raw_epi_recovery %>% filter(group == "reference"), aes(x = x, ymax = conf.high, ymin = conf.low, group = group), alpha = 0.2) +
#   
#   theme_bw() + 
#   scale_x_continuous(breaks = seq(-8, 6, by = 1), minor_breaks = NULL) +
#   # scale_y_continuous(limits = c(-10, 155)) +
#   coord_cartesian(ylim = c(4, 155)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   guides(color = guide_legend(keyheight = 0.6),
#          shape = guide_legend(keyheight = 0.6),
#          lty = guide_legend(keyheight = 0.6),
#          keyheight = 1) +
#   labs(x = "Time since end of reference (years)", 
#        y = "Biomass (dry g/m\U00B2)",
#        title = "(e) Reference")
# raw_epi_reference
# 
# # data frame of predictions
# delta_epi_predictions_during <- predicted_raw_epi_during %>% 
#   as.data.frame() %>% 
#   select(x, group, predicted) %>% 
#   pivot_wider(names_from = group, values_from = predicted) %>% 
#   mutate(delta = removal - reference) %>% 
#   mutate(exp_dates = "during")
# 
# delta_epi_predictions_after <- predicted_raw_epi_recovery %>% 
#   as.data.frame() %>% 
#   select(x, group, predicted) %>% 
#   pivot_wider(names_from = group, values_from = predicted) %>% 
#   mutate(delta = removal - reference) %>% 
#   mutate(exp_dates = "after")
# 
# delta_epi_predictions <- ggplot() +
#   geom_vline(xintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "grey") +
#   annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
#            fill = "grey", alpha = 0.3) +
#   geom_point(data = delta_epi_continual,
#              aes(x = time_since_end, y = delta_continual_epi), 
#              shape = 2, 
#              alpha = 0.15,
#              size = 0.75) +
#   
#   # overall
#   geom_line(data = delta_epi_predictions_during, aes(x = x, y = delta), linewidth = 1) +
#   geom_line(data = delta_epi_predictions_after, aes(x = x, y = delta), linewidth = 1) +
#   
#   scale_x_continuous(limits = c(-8, 7), breaks = seq(-8, 7, by = 1), minor_breaks = NULL) +
#   coord_cartesian(ylim = c(-40, 75)) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         legend.text = element_text(size = 6), 
#         legend.title = element_text(size = 6),
#         # plot.margin = margin(0, 0, 0, 0),
#         legend.position = "none",
#         panel.grid = element_blank(),
#         plot.title.position = "plot",
#         plot.title = element_text(size = 10)) +
#   labs(x = "Time since end of removal (years)", 
#        y = "\U0394 biomass \n (removal \U2212 reference, dry g/m\U00B2)",
#        fill = "Site",
#        shape = "Site",
#        title = "(f) Removal \U2212 reference")
# 
# delta_epi_predictions
# 
# ##########################################################################-
# # 5. manuscript tables ----------------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. model summary tables -----------------------------------------------
# 
# # individual group tables
# # lm_algae_tables <- tbl_merge(tbls = list(lm_algae_during_summary, lm_algae_recovery_summary),
# #                              tab_spanner = c("**Removal**", "**Recovery**")) 
# 
# # lm_epi_tables <- tbl_merge(tbls = list(lm_epi_during_summary, lm_epi_recovery_summary),
# #                            tab_spanner = c("**Removal**", "**Recovery**")) 
# 
# # lm_endo_tables <- tbl_merge(tbls = list(lm_endo_during_summary, lm_endo_recovery_summary),
# #                             tab_spanner = c("**Removal**", "**Recovery**")) 
# 
# # stack tables
# # lm_summary_tables <- tbl_stack(
# #   tbls = list(lm_kelp_tables, lm_algae_tables, lm_epi_tables, lm_endo_tables),
# #   group_header = c("Kelp", "Understory macroalgae", "Epilithic invertebrates", "Endolithic invertebrates"),
# #   quiet = TRUE) %>% 
# #   as_flex_table() %>% 
# #   font(fontname = "Times New Roman", part = "all")
# 
# # lm_summary_tables %>%
# #   save_as_docx(path = here::here("tables", "ms-tables", paste("tbl-S1_", today(), ".docx", sep = "")))
# 
# # tbl_stack(
# #   tbls = list(lm_kelp_tables, lm_algae_tables, lm_epi_tables, lm_endo_tables),
# #   group_header = c("Kelp", "Understory macroalgae", "Epilithic invertebrates", "Endolithic invertebrates"),
# #   quiet = TRUE) %>%
# #   as_gt() %>%
# #   tab_options(table.font.names = "Times New Roman") %>%
# #   gtsave(here::here("tables", "ms-tables", paste("tbl-S1_", today(), ".png", sep = "")),
# #          vwidth = 1500, vheight = 1000)
# 
# # individual group tables
# lm_algae_zigamma_tables <- tbl_merge(tbls = list(lm_raw_algae_during_zigamma_summary, lm_raw_algae_recovery_zigamma_summary),
#                              tab_spanner = c("**Kelp removal**", "**Recovery**")) 
# 
# lm_ep_zigamma_tables <- tbl_merge(tbls = list(lm_raw_epi_during_zigamma_summary, lm_raw_epi_recovery_zigamma_summary),
#                            tab_spanner = c("**Kelp removal**", "**Recovery**")) 
# 
# # stack tables
# lm_zigamma_summary_tables <- tbl_stack(
#   tbls = list(lm_kelp_zigamma_tables, lm_algae_zigamma_tables, lm_ep_zigamma_tables),
#   group_header = c("Giant kelp", "Understory macroalgae", "Sessile invertebrates"),
#   quiet = TRUE) %>% 
#   as_flex_table() %>% 
#   font(fontname = "Times New Roman", part = "all")
# 
# lm_zigamma_summary_tables %>%
#   save_as_docx(path = here::here("tables", "ms-tables", paste("tbl-S1_", today(), ".docx", sep = "")))
# 
# ##########################################################################-
# # 6. manuscript figures ---------------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. raw biomass through time -------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S4_", today(), ".jpg", sep = "")),
# #        plot = delta_continual_sites_algae_raw,
# #        height = 12, width = 8, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S5_", today(), ".jpg", sep = "")),
# #        plot = delta_continual_sites_epi_raw,
# #        height = 12, width = 8, units = "cm",
# #        dpi = 300)
# # 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S6_", today(), ".jpg", sep = "")),
# #        plot = delta_continual_sites_endo_raw,
# #        height = 8, width = 16, units = "cm",
# #        dpi = 300)
# 
# # ⟞ b. raw algae and epi model --------------------------------------------
# 
# # fig2_v1 <-  (kelp_title + algae_title + epi_title) /
# #             (overall_kelp + overall_algae_predictions + overall_epi_predictions) /
# #             (overall_predictions + overall_algae_predictions + delta_epi_predictions) +
# #   plot_layout(heights = c(1, 10, 10), widths = c(1, 1, 1))
# # fig2_v1
# 
# # fig2_v2 <- (algae_title + epi_title) /
# #   (raw_algae_removal + raw_epi_removal) /
# #   (raw_algae_reference + raw_epi_reference) /
# #   (overall_algae_predictions + delta_epi_predictions) +
# #   plot_layout(heights = c(1, 10, 10, 10))
# # fig2_v2
# 
# kelp_column <- plot_grid(overall_kelp, overall_predictions, nrow = 2) %>% 
#   plot_grid(kelp_title, ., 
#             nrow = 2,
#             rel_heights = c(1, 20))
# 
# algae_column <- plot_grid(overall_algae_predictions, overall_algae_predictions, nrow = 2) %>% 
#   plot_grid(algae_title, ., 
#             nrow = 2,
#             rel_heights = c(1, 20))
# 
# epi_column <- plot_grid(overall_epi_predictions, delta_epi_predictions, nrow = 2) %>% 
#   plot_grid(epi_title, ., 
#             nrow = 2,
#             rel_heights = c(1, 20))
# 
# fig2_v1 <- plot_grid(kelp_column, algae_column, ncol = 2) %>% 
#   plot_grid(., epi_column,
#             ncol = 2, 
#             rel_widths = c(2, 1))
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-2_new-model_", today(), ".jpg", sep = "")),
# #        plot = raw_groups,
# #        height = 18, width = 14, units = "cm",
# #        dpi = 400)
# 
# # v1: legend position c(0.86, 0.88)
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-2_new-model_v1_", today(), ".jpg", sep = "")),
# #        plot = fig2_v1,
# #        height = 15, width = 24, units = "cm",
# #        dpi = 400)
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-2_new-model_v2_", today(), ".jpg", sep = "")),
# #        plot = fig2_v2,
# #        height = 24, width = 18, units = "cm",
# #        dpi = 400)
# 
# 
# 


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # -------------------------- OLD CODE BELOW HERE --------------------------
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # lm_delta_algae_kelp_after_m3 <- lmer(
# #   delta_continual_algae ~ delta_continual + (1|site), 
# #   data = delta_algae_continual %>% filter(exp_dates == "after")
# # )
# # lm_algae_kelp_after_m1 <- lmer(
# #   continual_algae ~ continual + (1|year) + (1|site),
# #   data = delta_algae_continual %>% filter(exp_dates == "after") %>% filter(!(sample_ID %in% c("mohk_2019-11-19_Q4", "mohk_2021-05-13_Q2", "napl_2023-05-18_Q2")))
# # )
# # lm_algae_kelp_after_m2 <- lmer(
# #   continual_algae ~ continual*time_since_end + (1|site),
# #   data = delta_algae_continual %>% filter(exp_dates == "after")
# # )
# 
# # diagnostics
# # check_model(lm_delta_algae_kelp_after_m1)
# # simulateResiduals(lm_delta_algae_kelp_after_m1, plot = TRUE)
# 
# simulateResiduals(lm_delta_algae_kelp_after_m2, plot = TRUE)
# 
# outlier_check <- check_outliers(lm_delta_algae_kelp_after_m2_outlier) %>% plot()
# simulateResiduals(lm_delta_algae_kelp_after_m2_outlier, plot = TRUE)
# # check_model(lm_delta_algae_kelp_after_m3)
# 
# # check_model(lm_algae_kelp_after_m1)
# # plot(simulateResiduals(lm_algae_kelp_after_m1))
# 
# # Rsquared
# # r.squaredGLMM(lm_delta_algae_kelp_after_m1)
# r.squaredGLMM(lm_delta_algae_kelp_after_m2)
# # r.squaredGLMM(lm_delta_algae_kelp_after_m3)
# 
# # r.squaredGLMM(lm_algae_kelp_after_m1)
# # r.squaredGLMM(lm_algae_kelp_after_m2)
# 
# # summaries
# # summary(lm_delta_algae_kelp_after_m1)
# summary(lm_delta_algae_kelp_after_m2)
# # summary(lm_delta_algae_kelp_after_m3)
# 
# lm_delta_algae_kelp_after_summary <- lm_delta_algae_kelp_after_m2 %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**"
#   ) 
# 
# # AIC comparison
# # AICc(lm_delta_algae_kelp_after_m1, lm_delta_algae_kelp_after_m2, lm_delta_algae_kelp_after_m3) %>% 
# #   arrange(AICc)
# # AICc(lm_algae_kelp_after_m1, lm_algae_kelp_after_m2)
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # raw biomass
# # predicted_algae_kelp_after <- ggpredict(lm_algae_kelp_after_m1, terms = ~ continual, type = "fixed")
# 
# # deltas
# predicted_delta_algae_vs_kelp <- ggpredict(lm_delta_algae_kelp_after_m2, terms = ~ delta_continual, type = "fixed")
# 
# predicted_delta_algae_vs_kelp_outlier <- ggpredict(lm_delta_algae_kelp_after_m2_outlier, terms = ~ delta_continual, type = "fixed")
# 
# # ⟞ c. figures ------------------------------------------------------------
# 
# # ⟞ ⟞ i. raw biomass ------------------------------------------------------
# 
# # algae_kelp_after_plot <- ggplot(data = delta_algae_continual %>% 
# #                                   filter(exp_dates == "after"), 
# #                                 aes(x = continual, y = continual_algae)) +
# #   geom_point(shape = 21) +
# #   geom_line(data = predicted_algae_kelp_after, aes(x = x, y = predicted), lty = 1) +
# #   geom_ribbon(data = predicted_algae_kelp_after, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2)
# # 
# # algae_kelp_after_plot
# 
# # ⟞ ⟞ ii. deltas ----------------------------------------------------------
# 
# delta_algae_vs_kelp_lm <- delta_algae_continual %>% 
#   filter(exp_dates == "after" & sample_ID != "napl_2023-05-18_Q2") %>% 
#   # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
#   ggplot(aes(x = delta_continual, y = delta_continual_algae)) +
#   geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.5) +
#   geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.5) +
#   geom_point(size = 1, shape = 5, alpha = 0.4, color = under_col) + 
#   geom_ribbon(data = predicted_delta_algae_vs_kelp, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
#   geom_line(data = predicted_delta_algae_vs_kelp, aes(x = x, y = predicted), 
#             linewidth = 1,
#             color = under_col) +
#   scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
#   scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
#   labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
#        y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)", 
#        title = "(a) Understory macroalgae") +
#   annotate("text", x = -1100, y = -200,
#            label = "conditional R\U00B2 = 0.44\nmarginal R\U00B2 = 0.26\np < 0.001",
#            size = 1.5) +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 6),
#         axis.text = element_text(size = 5),
#         plot.title = element_text(size = 8),
#         plot.title.position = "plot",
#         panel.grid = element_blank()) 
# delta_algae_vs_kelp_lm
# 
# delta_algae_vs_kelp_lm_outlier <- delta_algae_continual %>% 
#   filter(exp_dates == "after") %>% 
#   # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
#   ggplot(aes(x = delta_continual, y = delta_continual_algae)) +
#   geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.5) +
#   geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.5) +
#   geom_point(size = 1, shape = 5, alpha = 0.4, color = under_col) + 
#   geom_ribbon(data = predicted_delta_algae_vs_kelp_outlier, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
#   geom_line(data = predicted_delta_algae_vs_kelp_outlier, aes(x = x, y = predicted), 
#             linewidth = 1,
#             color = under_col) +
#   scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
#   scale_y_continuous(breaks = seq(-200, 400, by = 200)) +
#   labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
#        y = "\U0394 understory macroalgae biomass\n(removal - reference, dry g/m\U00B2)", 
#        title = "(a) Understory macroalgae") +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 6),
#         axis.text = element_text(size = 5),
#         plot.title = element_text(size = 8),
#         plot.title.position = "plot",
#         panel.grid = element_blank()) 
# delta_algae_vs_kelp_lm_outlier
# 
# # delta_algae_vs_kelp_pearson <- delta_algae_continual %>% 
# #   mutate(exp_dates = case_when(
# #     exp_dates == "during" ~ "During removal",
# #     exp_dates == "after" ~ "Recovery period"
# #   )) %>% 
# #   ggplot(aes(x = delta_continual, y = delta_continual_algae, linetype = exp_dates, color = exp_dates, fill = exp_dates)) +
# #   geom_hline(aes(yintercept = 0), lty = 2) +
# #   geom_vline(aes(xintercept = 0), lty = 2) +
# #   geom_point(aes(shape = exp_dates, fill = exp_dates), size = 5, shape = 21, color = "#000000") +
# #   stat_cor(method = "pearson", cor.coef.name = "rho", inherit.aes = TRUE) +
# #   geom_smooth(aes(color = exp_dates), method = "lm", color = "black", linewidth = 3, se = FALSE) +
# #   scale_linetype_manual(values = c("During removal" = 2, "Recovery period" = 1)) +
# #   scale_color_manual(values = c("During removal" = "grey", "Recovery period" = under_col)) +
# #   scale_fill_manual(values = c("During removal" = "#FFFFFF", "Recovery period" = under_col)) +
# #   scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
# #   labs(x = "\U0394 kelp biomass (treatment - control)",
# #        y = "\U0394 understory macroalgae biomass (treatment - control)") +
# #   theme_bw() + 
# #   theme(axis.title = element_text(size = 18),
# #         plot.title = element_text(size = 18),
# #         axis.text = element_text(size = 16),
# #         legend.position = c(0.83, 0.93),
# #         legend.text = element_text(size = 18),
# #         legend.title = element_blank())
# # delta_algae_vs_kelp_pearson
# 
# ##########################################################################-
# # 2. epilithic inverts ----------------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. correlation --------------------------------------------------------
# 
# delta_epi_after <- delta_epi_continual %>% 
#   filter(exp_dates == "after")
# 
# cor.test(delta_epi_after$delta_continual_epi, delta_epi_after$delta_continual,
#          method = "pearson")
# 
# # ⟞ b. linear models ------------------------------------------------------
# 
# # ⟞ ⟞ i. model and diagnostics  -------------------------------------------
# 
# # lm_delta_epi_kelp_after_m1 <- lm(
# #   delta_continual_epi ~ delta_continual, 
# #   data = delta_epi_continual %>% filter(exp_dates == "after"),
# #   na.action = na.omit
# # )
# lm_delta_epi_kelp_after_m2 <- lmer(
#   delta_continual_epi ~ delta_continual + (1|year) + (1|site), 
#   data = delta_epi_continual %>% filter(exp_dates == "after")
# )
# # lm_delta_epi_kelp_after_m3 <- lmer(
# #   delta_continual_epi ~ delta_continual + (1|site), 
# #   data = delta_epi_continual %>% filter(exp_dates == "after")
# # )
# 
# # lm_epi_kelp_after_m1 <- lmer(
# #   continual_epi ~ continual + (1|year) + (1|site),
# #   data = delta_epi_continual %>% filter(exp_dates == "after")
# # )
# 
# # lm_epi_kelp_after_m2 <- lmer(
# #   continual_epi ~ continual*time_since_end + (1|site),
# #   data = delta_epi_continual %>% filter(exp_dates == "after")
# # )
# 
# # diagnostics
# # simulateResiduals(lm_delta_epi_kelp_after_m1, plot = T)
# # check_model(lm_delta_epi_kelp_after_m1)
# 
# simulateResiduals(lm_delta_epi_kelp_after_m2, plot = T)
# 
# # simulateResiduals(lm_delta_epi_kelp_after_m3, plot = T)
# # check_model(lm_delta_epi_kelp_after_m3)
# 
# # check_model(lm_epi_kelp_after_m1)
# # check_model(lm_epi_kelp_after_m2)
# 
# # Rsquared
# # r.squaredGLMM(lm_delta_epi_kelp_after_m1)
# r.squaredGLMM(lm_delta_epi_kelp_after_m2)
# # r.squaredGLMM(lm_delta_epi_kelp_after_m3)
# 
# # r.squaredGLMM(lm_epi_kelp_after_m1)
# # r.squaredGLMM(lm_epi_kelp_after_m2)
# 
# # summaries
# # summary(lm_delta_epi_kelp_after_m1)
# summary(lm_delta_epi_kelp_after_m2)
# # summary(lm_delta_epi_kelp_after_m3)
# 
# lm_delta_epi_kelp_after_summary <- lm_delta_epi_kelp_after_m2 %>% 
#   tbl_regression() %>% 
#   bold_p(t = 0.05) %>% 
#   modify_header(
#     label = " ",
#     estimate = "**Slope**"
#   ) 
# lm_delta_epi_kelp_after_summary
# 
# # AICc
# # AICc(lm_delta_epi_kelp_after_m1, lm_delta_epi_kelp_after_m2, lm_delta_epi_kelp_after_m3) %>% 
# #   arrange(AICc)
# 
# AICc(lm_epi_kelp_after_m1, lm_epi_kelp_after_m2) # same?
# 
# # ⟞ ⟞ ii. predictions -----------------------------------------------------
# 
# # raw biomass
# # predicted_epi_vs_kelp <- ggpredict(lm_epi_kelp_after_m1, terms = ~ continual, type = "fixed")
# 
# # deltas
# predicted_delta_epi_vs_kelp <- ggpredict(lm_delta_epi_kelp_after_m2, terms = ~ delta_continual, type = "fixed")
# 
# # ⟞ c. figures ------------------------------------------------------------
# 
# # ⟞ ⟞ i. raw biomass ------------------------------------------------------
# 
# # delta_epi_continual %>% 
# #   filter(exp_dates == "during") %>% 
# #   ggplot(aes(x = continual, y = continual_epi)) +
# #   geom_point() +
# #   geom_line(data = predicted_epi_vs_kelp, aes(x = x, y = predicted))
# 
# # ⟞ ⟞ ii. deltas ----------------------------------------------------------
# 
# delta_epi_vs_kelp_lm <- delta_epi_continual %>% 
#   filter(exp_dates == "after") %>% 
#   # two points missing from delta kelp: MOHK 2010-06-14, NAPL 2014-11-14
#   ggplot(aes(x = delta_continual, y = delta_continual_epi)) +
#   geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.5) +
#   geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.5) +
#   geom_point(size = 1, shape = 5, alpha = 0.4) + 
#   # geom_ribbon(data = predicted_delta_epi_vs_kelp, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
#   # geom_line(data = predicted_delta_epi_vs_kelp, aes(x = x, y = predicted), size = 2) +
#   # scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
#   labs(x = "\U0394 giant kelp biomass\n(removal - reference, dry g/m\U00B2)",
#        y = "\U0394 sessile invertebrate biomass\n(removal - reference, dry g/m\U00B2)", 
#        title = "(b) Sessile invertebrates") +
#   theme_bw() + 
#   theme(axis.title = element_text(size = 6),
#         axis.text = element_text(size = 5),
#         plot.title = element_text(size = 8),
#         plot.title.position = "plot",
#         panel.grid = element_blank()) 
# delta_epi_vs_kelp_lm
# 
# # epi_vs_kelp_pearson <- delta_epi_continual %>% 
# #   mutate(exp_dates = case_when(
# #     exp_dates == "during" ~ "During removal",
# #     exp_dates == "after" ~ "Post-removal"
# #   )) %>% 
# #   ggplot(aes(x = delta_continual, y = delta_continual_epi, linetype = exp_dates, color = exp_dates, fill = exp_dates)) +
# #   stat_cor(aes(color = exp_dates), method = "pearson", cor.coef.name = "rho") +
# #   geom_hline(aes(yintercept = 0), lty = 2) +
# #   geom_vline(aes(xintercept = 0), lty = 2) +
# #   geom_point(size = 5, shape = 25, color = "#000000") + 
# #   geom_smooth(aes(linetype = exp_dates), method = "lm", se = FALSE, color = "black", size = 3) +
# #   scale_linetype_manual(values = c("During removal" = 2, "Post-removal" = 1)) +
# #   scale_fill_manual(values = c("During removal" = "#FFFFFF", "Post-removal" = "#54662C")) +
# #   scale_color_manual(values = c("During removal" = "grey", "Post-removal" = "#54662C")) +
# #   scale_x_continuous(breaks = seq(-2000, 2000, by = 1000), minor_breaks = seq(-2000, 2000, by = 500)) +
# #   labs(x = "\U0394 kelp biomass (treatment - control)",
# #        y = "\U0394 epilithic invertebrate biomass (treatment - control)") +
# #   theme_bw() + 
# #   theme(axis.title = element_text(size = 18),
# #         plot.title = element_text(size = 18),
# #         axis.text = element_text(size = 16),
# #         legend.position = c(0.83, 0.93),
# #         legend.text = element_text(size = 18),
# #         legend.title = element_blank())
# # epi_vs_kelp_pearson
# 
# ##########################################################################-
# # 3. manuscript tables ----------------------------------------------------
# ##########################################################################-
# 
# lm_vs_kelp_summary_tables <- tbl_stack(
#   tbls = list(lm_delta_algae_kelp_after_summary, lm_delta_epi_kelp_after_summary),
#   group_header = c("Understory macroalgae", "Sessile invertebrates"),
#   quiet = TRUE) %>% 
#   as_flex_table() %>% 
#   font(fontname = "Times New Roman", part = "all")
# 
# # lm_vs_kelp_summary_tables %>%
# #   save_as_docx(path = here::here("tables", "ms-tables", paste("tbl-S5_", today(), ".docx", sep = "")))
# 
# 
# ##########################################################################-
# # 4. manuscript figures ---------------------------------------------------
# ##########################################################################-
# 
# # ⟞ a. correlation --------------------------------------------------------
# 
# algae_vs_kelp_spearman
# 
# # ⟞ b. linear models ------------------------------------------------------
# 
# group_vs_kelp <- plot_grid(delta_algae_vs_kelp_lm, delta_epi_vs_kelp_lm, ncol = 2)
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-4_", today(), ".jpg", sep = "")),
# #        plot = group_vs_kelp,
# #        height = 6, width = 12, units = "cm",
# #        dpi = 300)
# 
# 
# # ⟞ c. outlier check ------------------------------------------------------
# 
# # ggsave(here::here("figures", "ms-figures",
# #                   paste("fig-S12_", today(), ".jpg", sep = "")),
# #        plot = outlier_check,
# #        height = 8, width = 14, units = "cm",
# #        dpi = 300)

