
# 0. source ---------------------------------------------------------------

source(here::here("code", "00-set_up.R"))

# 1. biomass timeseries ---------------------------------------------------

############################################
# a. species specific function
############################################

# i. summary data frame to feed into function
biomass_summary <- biomass %>% 
  # calculate total dry mass by each survey
  group_by(site, year, month, treatment, date, sp_code) %>%
  mutate(total_biomass = sum(dry_gm2)) %>% 
  ungroup() %>% 
  # make a new column designating "start", "during" and "after" removal
  mutate(exp_dates = case_when(
    site == "NAPL" & sample_ID %in% napl_start_dates ~ "start",
    site == "NAPL" & date > napl_after_date ~ "after",
    site == "MOHK" & sample_ID %in% mohk_start_dates ~ "start",
    site == "MOHK" & date > mohk_after_date ~ "after",
    site == "AQUE" & sample_ID %in% aque_start_dates ~ "start",
    site == "AQUE" & date > aque_after_date ~ "after",
    site == "CARP" & sample_ID %in% carp_start_dates ~ "start",
    site == "CARP" & date > carp_after_date ~ "after",
    site == "IVEE" & sample_ID %in% ivee_start_dates ~ "start",
    site == "IVEE" & date > ivee_after_date ~ "after",
    TRUE ~ "during"
  ),
  exp_dates = factor(exp_dates, levels = c("start", "during", "after"))) 

spp_biomass_ts <- function(spp_code, this_site) {
  
  start_dates <- if(this_site == "NAPL") {
    napl_start_dates
  } else if(this_site == "MOHK") {
    mohk_start_dates
  } else if(this_site == "AQUE") {
    aque_start_dates
  } else if(this_site == "CARP") {
    carp_start_dates
  } else if(this_site == "IVEE") {
    ivee_start_dates
  }
  
  end_dates <- if(this_site == "NAPL") {
    napl_after_date
  } else if (this_site == "MOHK") {
    mohk_after_date
  } else if(this_site == "AQUE") {
    aque_after_date
  } else if(this_site == "CARP") {
    carp_after_date
  } else if(this_site == "IVEE") {
    ivee_after_date
  }
  
  full_sci_name <- spp_names %>% filter(sp_code == spp_code) %>% pull(scientific_name)
  
  full_common_name <- spp_names %>% filter(sp_code == spp_code) %>% pull(common_name)
  
  plot_label <- paste(full_sci_name, " (", full_common_name, ")", sep = "")
  
  biomass_summary %>% 
    filter(sp_code == spp_code & site == this_site) %>% 
    
    ggplot(., aes(x = date, y = total_biomass)) +
    geom_line(alpha = 0.5) +
    geom_point(aes(color = exp_dates)) +
    labs(title = plot_label,
         subtitle = this_site,
         x = "Year",
         y = "Total dry mass",
         color = "Experiment") +
    facet_wrap(~treatment) +
    theme_bw()
  
}

############################################
# b. tests
############################################

spp_biomass_ts("PTCA", "NAPL")
spp_biomass_ts("PTCA", "MOHK")
spp_biomass_ts("PTCA", "IVEE")
spp_biomass_ts("CYOS", "NAPL")
spp_biomass_ts("AB", "NAPL")
spp_biomass_ts("ANSP", "NAPL")

############################################
# c. loops
############################################

#### algae
for(j in 1:length(LTE_sites)) {
  # choose a site
  this_site <- LTE_sites[j]
  
  # list name
  list_name <- paste(this_site, "_algae_plot_list", sep = "")
  
  # create an empty list to hold plot names
  loop_plot_list <- rep(NaN, length(algae_spp))
  
  for(i in 1:length(algae_spp)) {
    # choose a species from the list
    this_sp <- algae_spp[i]
    
    # plot the species biomass timeseries
    plot <- spp_biomass_ts(this_sp, this_site)
    
    # plot name list in the loop
    loop_plot_list[[i]] <- paste(this_site, this_sp, sep = "_")
    
    # make a separate object for the plot
    assign(paste(this_site, this_sp, sep = "_"), plot)
  }
  
  # put that plot list into list_name
  assign(list_name, loop_plot_list)
  
}

#### fish

for(j in 1:length(LTE_sites)) {
  # choose a site
  this_site <- LTE_sites[j]
  
  # list name
  list_name <- paste(this_site, "_fish_plot_list", sep = "")
  
  # create an empty list to hold plot names
  loop_plot_list <- rep(NaN, length(fish_spp))
  
  for(i in 1:length(fish_spp)) {
    # choose a species from the list
    this_sp <- fish_spp[i]
    
    # plot the species biomass timeseries
    plot <- spp_biomass_ts(this_sp, this_site)
    
    # plot name list in the loop
    loop_plot_list[[i]] <- paste(this_site, this_sp, sep = "_")
    
    # make a separate object for the plot
    assign(paste(this_site, this_sp, sep = "_"), plot)
  }
  
  # put that plot list into list_name
  assign(list_name, loop_plot_list)
  
}

#### inverts

for(j in 1:length(LTE_sites)) {
  # choose a site
  this_site <- LTE_sites[j]
  
  # list name
  list_name <- paste(this_site, "_invert_plot_list", sep = "")
  
  # create an empty list to hold plot names
  loop_plot_list <- rep(NaN, length(invert_spp))
  
  for(i in 1:length(invert_spp)) {
    # choose a species from the list
    this_sp <- invert_spp[i]
    
    # plot the species biomass timeseries
    plot <- spp_biomass_ts(this_sp, this_site)
    
    # plot name list in the loop
    loop_plot_list[[i]] <- paste(this_site, this_sp, sep = "_")
    
    # make a separate object for the plot
    assign(paste(this_site, this_sp, sep = "_"), plot)
  }
  
  # put that plot list into list_name
  assign(list_name, loop_plot_list)
  
}

############################################
# d. group specific function
############################################

# data frame to feed into function
biomass_group_summary <- biomass %>% 
  # calculate total dry mass by each survey
  group_by(site, year, month, treatment, date, group) %>%
  mutate(total_biomass = sum(dry_gm2)) %>% 
  ungroup() %>% 
  # make a new column designating "start", "during" and "after" removal
  mutate(exp_dates = case_when(
    site == "NAPL" & sample_ID %in% napl_start_dates ~ "start",
    site == "NAPL" & date > napl_after_date ~ "after",
    site == "MOHK" & sample_ID %in% mohk_start_dates ~ "start",
    site == "MOHK" & date > mohk_after_date ~ "after",
    site == "AQUE" & sample_ID %in% aque_start_dates ~ "start",
    site == "AQUE" & date > aque_after_date ~ "after",
    site == "CARP" & sample_ID %in% carp_start_dates ~ "start",
    site == "CARP" & date > carp_after_date ~ "after",
    site == "IVEE" & sample_ID %in% ivee_start_dates ~ "start",
    site == "IVEE" & date > ivee_after_date ~ "after",
    TRUE ~ "during"
  ),
  exp_dates = factor(exp_dates, levels = c("start", "during", "after"))) 

group_biomass <- function(this_group, this_site){
  start_dates <- if(this_site == "NAPL") {
    napl_start_dates
  } else if(this_site == "MOHK") {
    mohk_start_dates
  } else if(this_site == "AQUE") {
    aque_start_dates
  } else if(this_site == "CARP") {
    carp_start_dates
  } else if(this_site == "IVEE") {
    ivee_start_dates
  }
  
  end_dates <- if(this_site == "NAPL") {
    napl_after_date
  } else if (this_site == "MOHK") {
    mohk_after_date
  } else if(this_site == "AQUE") {
    aque_after_date
  } else if(this_site == "CARP") {
    carp_after_date
  } else if(this_site == "IVEE") {
    ivee_after_date
  }
  
  group_title <- if(this_group == "invert") {
    "Invertebrates"
  } else if(this_group == "algae") {
    "Algae"
  } else if(this_group == "fish") {
    "Fish"
  }
  
  biomass_group_summary %>% 
    filter(group == this_group & site == this_site) %>%
    
    ggplot(., aes(x = date, y = total_biomass)) +
    geom_line(alpha = 0.5) +
    geom_point(aes(color = exp_dates)) +
    labs(title = group_title,
         subtitle = this_site,
         x = "Year",
         y = "Total dry mass",
         color = "Experiment") +
    facet_wrap(~treatment) +
    theme_bw()
  
}

############################################
# e. tests
############################################

group_biomass("invert", "NAPL")
group_biomass("algae", "NAPL")
group_biomass("fish", "NAPL")
