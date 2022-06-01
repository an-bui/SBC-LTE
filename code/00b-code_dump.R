
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
# 
# 
# 
