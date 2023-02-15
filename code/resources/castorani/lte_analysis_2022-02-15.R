# ====================================================================================
# SBC LTER LONG TERM EXPERIMENT - RECOVERY ANALYSIS
# ====================================================================================

# Max Castorani
# Updated Feb 15, 2022

# ====================================================================================
# LOAD LIBRARIES
# ====================================================================================

# ------------------------------------------------------------------------------------
# Clear environment
# rm(list = ls())

# ------------------------------------------------------------------------------------
# Set working environment
# setwd("/Users/castorani/Sync/Documents/Research/Active projects/SBC LTE recovery (Castorani, Reed, Stier)/Analyses/")

source(here::here("code", "00-set_up.R"))

# for (package in c('dplyr', 'tidyr', 'ggplot2', 'vegan', 'vegetarian', 'boot')) {
#   if (!require(package, character.only=T, quietly=T)) {
#     install.packages(package)
#     library(package, character.only=T)
#   }
# }

# ====================================================================================
# CUSTOM FUNCTIONS
# ====================================================================================

# ------------------------------------------------------------------------------------
sample.mean <- function(x, d) {
  return(mean(x[d]))
}

# ------------------------------------------------------------------------------------
boot.fun <- function(dat) {
  b <- boot::boot(data = dat,
            statistic = sample.mean, R = 1000)
  ci <- boot::boot.ci(b, type = "bca")$bca[4:5]
  return(ci)
}

# ====================================================================================
# IMPORT AND TIDY DATA
# ====================================================================================
# Package ID: knb-lter-sbc.119.7 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Long-term experiment: biomass of kelp forest species, ongoing since 2008.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Data set creator:  Robert J Miller -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/119/7/020f9b0e507561cb49fc8fd122da0a29" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "YEAR",     
                 "MONTH",     
                 "DATE",     
                 "SITE",     
                 "TRANSECT",     
                 "TREATMENT",     
                 "SP_CODE",     
                 "PERCENT_COVER",     
                 "DENSITY",     
                 "WM_GM2",     
                 "DRY_GM2",     
                 "SFDM",     
                 "AFDM",     
                 "SCIENTIFIC_NAME",     
                 "COMMON_NAME",     
                 "TAXON_KINGDOM",     
                 "TAXON_PHYLUM",     
                 "TAXON_CLASS",     
                 "TAXON_ORDER",     
                 "TAXON_FAMILY",     
                 "TAXON_GENUS",     
                 "GROUP",     
                 "MOBILITY",     
                 "GROWTH_MORPH"    ), check.names=TRUE)
unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings
if (class(dt1$MONTH)!="factor") dt1$MONTH<- as.factor(dt1$MONTH)                                   
# attempting to convert dt1$DATE dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1DATE<-as.Date(dt1$DATE,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1DATE) == length(tmp1DATE[!is.na(tmp1DATE)])){dt1$DATE <- tmp1DATE } else {print("Date conversion failed for dt1$DATE. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1DATE) 
if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$TREATMENT)!="factor") dt1$TREATMENT<- as.factor(dt1$TREATMENT)
if (class(dt1$SP_CODE)!="factor") dt1$SP_CODE<- as.factor(dt1$SP_CODE)
if (class(dt1$PERCENT_COVER)=="factor") dt1$PERCENT_COVER <-as.numeric(levels(dt1$PERCENT_COVER))[as.integer(dt1$PERCENT_COVER) ]               
if (class(dt1$PERCENT_COVER)=="character") dt1$PERCENT_COVER <-as.numeric(dt1$PERCENT_COVER)
if (class(dt1$DENSITY)=="factor") dt1$DENSITY <-as.numeric(levels(dt1$DENSITY))[as.integer(dt1$DENSITY) ]               
if (class(dt1$DENSITY)=="character") dt1$DENSITY <-as.numeric(dt1$DENSITY)
if (class(dt1$WM_GM2)=="factor") dt1$WM_GM2 <-as.numeric(levels(dt1$WM_GM2))[as.integer(dt1$WM_GM2) ]               
if (class(dt1$WM_GM2)=="character") dt1$WM_GM2 <-as.numeric(dt1$WM_GM2)
if (class(dt1$DRY_GM2)=="factor") dt1$DRY_GM2 <-as.numeric(levels(dt1$DRY_GM2))[as.integer(dt1$DRY_GM2) ]               
if (class(dt1$DRY_GM2)=="character") dt1$DRY_GM2 <-as.numeric(dt1$DRY_GM2)
if (class(dt1$SFDM)=="factor") dt1$SFDM <-as.numeric(levels(dt1$SFDM))[as.integer(dt1$SFDM) ]               
if (class(dt1$SFDM)=="character") dt1$SFDM <-as.numeric(dt1$SFDM)
if (class(dt1$AFDM)=="factor") dt1$AFDM <-as.numeric(levels(dt1$AFDM))[as.integer(dt1$AFDM) ]               
if (class(dt1$AFDM)=="character") dt1$AFDM <-as.numeric(dt1$AFDM)
if (class(dt1$SCIENTIFIC_NAME)!="factor") dt1$SCIENTIFIC_NAME<- as.factor(dt1$SCIENTIFIC_NAME)
if (class(dt1$COMMON_NAME)!="factor") dt1$COMMON_NAME<- as.factor(dt1$COMMON_NAME)
if (class(dt1$TAXON_KINGDOM)!="factor") dt1$TAXON_KINGDOM<- as.factor(dt1$TAXON_KINGDOM)
if (class(dt1$TAXON_PHYLUM)!="factor") dt1$TAXON_PHYLUM<- as.factor(dt1$TAXON_PHYLUM)
if (class(dt1$TAXON_CLASS)!="factor") dt1$TAXON_CLASS<- as.factor(dt1$TAXON_CLASS)
if (class(dt1$TAXON_ORDER)!="factor") dt1$TAXON_ORDER<- as.factor(dt1$TAXON_ORDER)
if (class(dt1$TAXON_FAMILY)!="factor") dt1$TAXON_FAMILY<- as.factor(dt1$TAXON_FAMILY)
if (class(dt1$TAXON_GENUS)!="factor") dt1$TAXON_GENUS<- as.factor(dt1$TAXON_GENUS)
if (class(dt1$GROUP)!="factor") dt1$GROUP<- as.factor(dt1$GROUP)
if (class(dt1$MOBILITY)!="factor") dt1$MOBILITY<- as.factor(dt1$MOBILITY)
if (class(dt1$GROWTH_MORPH)!="factor") dt1$GROWTH_MORPH<- as.factor(dt1$GROWTH_MORPH)

# Convert Missing Values to NA for non-dates
dt1$PERCENT_COVER <- ifelse((trimws(as.character(dt1$PERCENT_COVER))==trimws("-99999")),NA,dt1$PERCENT_COVER)               
suppressWarnings(dt1$PERCENT_COVER <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$PERCENT_COVER))==as.character(as.numeric("-99999"))),NA,dt1$PERCENT_COVER))
dt1$DENSITY <- ifelse((trimws(as.character(dt1$DENSITY))==trimws("-99999")),NA,dt1$DENSITY)               
suppressWarnings(dt1$DENSITY <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$DENSITY))==as.character(as.numeric("-99999"))),NA,dt1$DENSITY))
dt1$WM_GM2 <- ifelse((trimws(as.character(dt1$WM_GM2))==trimws("-99999")),NA,dt1$WM_GM2)               
suppressWarnings(dt1$WM_GM2 <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$WM_GM2))==as.character(as.numeric("-99999"))),NA,dt1$WM_GM2))
dt1$DRY_GM2 <- ifelse((trimws(as.character(dt1$DRY_GM2))==trimws("-99999")),NA,dt1$DRY_GM2)               
suppressWarnings(dt1$DRY_GM2 <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$DRY_GM2))==as.character(as.numeric("-99999"))),NA,dt1$DRY_GM2))
dt1$SFDM <- ifelse((trimws(as.character(dt1$SFDM))==trimws("-99999")),NA,dt1$SFDM)               
suppressWarnings(dt1$SFDM <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$SFDM))==as.character(as.numeric("-99999"))),NA,dt1$SFDM))
dt1$AFDM <- ifelse((trimws(as.character(dt1$AFDM))==trimws("-99999")),NA,dt1$AFDM)               
suppressWarnings(dt1$AFDM <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$AFDM))==as.character(as.numeric("-99999"))),NA,dt1$AFDM))
dt1$SCIENTIFIC_NAME <- as.factor(ifelse((trimws(as.character(dt1$SCIENTIFIC_NAME))==trimws("-99999")),NA,as.character(dt1$SCIENTIFIC_NAME)))

raw.dat <- dt1
rm(dt1, infile1, inUrl1)

# ------------------------------------------------------------------------------------
# Replace underscores with dots for convenience. Also convert to lowercase
colnames(raw.dat) <- tolower(gsub("_", ".", colnames(raw.dat)))  

raw.dat$group <- tolower(raw.dat$group) 
raw.dat$mobility <- tolower(raw.dat$mobility)
raw.dat$growth.morph <- tolower(raw.dat$growth.morph)

# ------------------------------------------------------------------------------------
# Replace "-99999" values with NA
raw.dat[raw.dat == -99999] <- NA

# ------------------------------------------------------------------------------------
# Import data describing guilds and add in
# guild.dat <- read.csv(file = "/Users/castorani/Sync/Data/SBC_LTER_data/Giant_kelp_LTE_data/LTE_guild_data.csv", header = TRUE) 

# read in LTE guild data
guild.dat <- read_csv(here::here("code/castorani", "LTE_guild_data.csv"))

raw.dat <- raw.dat %>%
  left_join(guild.dat, by = "sp.code")

# ------------------------------------------------------------------------------------
dat <- raw.dat %>%
  dplyr::mutate(month = as.numeric(as.character(month))) %>%
  dplyr::mutate(qt = ifelse(month <= 3, "Q1",
                     ifelse(month <= 6, "Q2",
                            ifelse(month <= 9, "Q3", "Q4")))) %>% # Create column for quarter
  mutate(yr.qt = year + ifelse(qt == "Q1", 0.125,
                               ifelse(qt == "Q2", 0.375,
                                      ifelse(qt == "Q3", 0.625, 0.875)))) %>%
  dplyr::mutate(time.temp = ifelse(month <= 3, 0.125,
                                        ifelse(month <= 6, 0.375,
                                               ifelse(month <= 9, 0.625, 0.875)))) %>%
  dplyr::mutate(time.yrs = year + time.temp - 2008) %>%
  dplyr::mutate(site = factor(plyr::revalue(site, c("MOHK" = "Mohawk",
                                             "AQUE" = "Arroyo Quemado",
                                             "NAPL" = "Naples",
                                             "IVEE" = "Isla Vista",
                                             "CARP" = "Carpinteria")),
                       levels = c("Mohawk",
                                  "Arroyo Quemado",
                                  "Naples",
                                  "Isla Vista",
                                  "Carpinteria")),
         treatment = factor(plyr::revalue(treatment, c("CONTROL" = "Control",
                                                       "ANNUAL" = "Annual disturbance",
                                                       "CONTINUAL" = "Quarterly disturbance")),
                            levels = c("Control", "Annual disturbance", "Quarterly disturbance"))) %>%
  dplyr::mutate(plot = paste(abbreviate(site, minlength = 2), abbreviate(treatment, minlength = 2), sep="_")) %>%
  dplyr::select(year, qt, yr.qt, time.yrs, month, date, site, treatment, plot, everything(), -transect, -time.temp, -month, -date) %>%
  dplyr::select(-group) %>%
  droplevels()

# ------------------------------------------------------------------------------------
# Create column for dry kg instead of dry g per m2
dat$biomass <- dat$dry.gm2 * 0.001  # CONVERT FROM g/m2 to kg/m2 in DRY MASS AS BIOMASS
  
# ------------------------------------------------------------------------------------
# Average early observations to quarterly
dat$yr.qt.plot <- paste(dat$yr.qt, dat$plot, sep = "_")

group.cols <- c("yr.qt.plot", "sp.code")
dots <-lapply(group.cols, as.symbol)

# Note: This next chunk of code takes a while to run (sorry! there is probably a more efficient way to do this, but it should work)
dat.temp <- dat %>%
  dplyr::group_by_(.dots = dots) %>%
  dplyr::summarise(year = unique(year),
                   qt = unique(qt),
                   time.yrs = unique(time.yrs),
                   yr.qt = unique(yr.qt),
                   site = unique(site),
                   treatment = unique(treatment),
                   plot = unique(plot),
                   taxon.genus = unique(taxon.genus),
                   scientific.name = unique(scientific.name),
                   sp.code = unique(sp.code),
                   biomass.guild = unique(biomass.guild),
                   diversity.guild = unique(diversity.guild),
                   biomass = mean(biomass, na.rm = TRUE)) %>%
  ungroup()

dat <- dat.temp %>%
  dplyr::select(year:plot, taxon.genus, scientific.name, sp.code, biomass.guild, diversity.guild, biomass)

remove(dat.temp, dots)

# ------------------------------------------------------------------------------------
# Add variable for time since experiment commenced
dat.temp <- dat %>%
  group_by(site) %>%
  mutate(time.since.start = time.yrs - min(time.yrs, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(year, qt, time.yrs, yr.qt, time.since.start, everything())

# Clean up
dat <- dat.temp
remove(dat.temp)

# ------------------------------------------------------------------------------------
# Add identifying code for each combination of site, treatment, and quarter
dat$id <- paste(dat$yr.qt, dat$plot, sep = ".")

# Reorder columns
dat <- dat %>%
  dplyr::select(id, everything())

# ------------------------------------------------------------------------------------
# Check number of taxa observations per year and per quarter
xtabs(~ site + year, data = dat)
xtabs(~ site + yr.qt, data = dat)

# Ensure that each plot (treatment x site) does not have multiple transects (so we are looking at plot mean density)
table(dat$plot, dat$site)


# ====================================================================================
# CONVERT FROM LONG TO WIDE FORMATS ... CREATE COMMUNITY MATRICES
# ====================================================================================

# ------------------------------------------------------------------------------------
# FULL COMMUNITY

dat.wide.full <- dat %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA) %>%
  dplyr::select(-`<NA>`)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.full.comm <- dat.wide.full %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.full.comm[is.na(dat.wide.full.comm)] <- 0

dat.wide.full <- dat.wide.full %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.full.comm)

# Check number of observations per year and per quarter
xtabs(~ site + year, data = dat.wide.full)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.full <- dat.wide.full.comm
rownames(comm.mat.full) <- dat.wide.full$treatment  
remove(dat.wide.full.comm)

# ------------------------------------------------------------------------------------
# CHECK TOTAL NUMBER OF OBSERVED TAXA IN EACH GUILD
no.taxa.obs <- dat %>%
  dplyr::select(id, sp.code, diversity.guild, biomass) %>%
  group_by(diversity.guild) %>%
  spread(., key = sp.code, value = biomass, fill = 0) %>%
  dplyr::select(-`<NA>`, -id) %>%
  dplyr::summarise_each(funs(sum)) %>%
  ungroup() %>%
  dplyr::mutate_each(funs(replace(., . > 0, 1)), -diversity.guild) %>%
  dplyr::mutate(no.taxa = rowSums(.[-1])) %>%
  dplyr::select(diversity.guild, no.taxa)
    
no.taxa.obs

unique(as.character(dat[is.na(dat$diversity.guild), ]$scientific.name))

# NOTE! THESE FOUR SPECIES, FOR SOME REASON, AREN'T GETTING PUT INTO THE RIGHT SPECIES GROUPS. PERHAPS THEY WERE ADDED TO THE MASTER SPECIES LIST AFTER 2018. EARLIER IN THE SPECIES LIST CSV FILE, THEY SHOULD BE CATEGORIZED INTO THE POSSIBLE GROUPINGS (E.G., UNDERSTORY ALGAE, EPILITHIC INVERTEBRATES, ETC.)

# ------------------------------------------------------------------------------------
# UNDERSTORY ALGAE
dat.wide.understory <- dat %>%
  dplyr::filter(diversity.guild == "algae") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA) %>%
  dplyr::select(-`<NA>`)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.understory.comm <- dat.wide.understory %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.understory.comm[is.na(dat.wide.understory.comm)] <- 0

dat.wide.understory <- dat.wide.understory %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.understory.comm)

xtabs(~ site + year, data = dat.wide.understory)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.understory <- dat.wide.understory.comm
rownames(comm.mat.understory) <- dat.wide.understory$treatment  
remove(dat.wide.understory.comm)

# ------------------------------------------------------------------------------------
# SESSILE INVERTEBRATES
dat.wide.sessile.inverts <- dat %>%
  dplyr::filter(diversity.guild == "sessile.invert") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA)
  
# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.sessile.inverts.comm <- dat.wide.sessile.inverts %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.sessile.inverts.comm[is.na(dat.wide.sessile.inverts.comm)] <- 0

dat.wide.sessile.inverts <- dat.wide.sessile.inverts %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.sessile.inverts.comm)

xtabs(~ site + year, data = dat.wide.sessile.inverts)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.sessile.inverts <- dat.wide.sessile.inverts.comm
rownames(comm.mat.sessile.inverts) <- dat.wide.sessile.inverts$treatment  
remove(dat.wide.sessile.inverts.comm)

# ------------------------------------------------------------------------------------
# EPILITHIC SESSILE INVERTEBRATES
dat.wide.epilithic.sessile.inverts <- dat %>%
  dplyr::filter(biomass.guild == "epilithic.sessile.invert") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.epilithic.sessile.inverts.comm <- dat.wide.epilithic.sessile.inverts %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.epilithic.sessile.inverts.comm[is.na(dat.wide.epilithic.sessile.inverts.comm)] <- 0

dat.wide.epilithic.sessile.inverts <- dat.wide.epilithic.sessile.inverts %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.epilithic.sessile.inverts.comm)

xtabs(~ site + year, data = dat.wide.epilithic.sessile.inverts)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.epilithic.sessile.inverts <- dat.wide.epilithic.sessile.inverts.comm
rownames(comm.mat.epilithic.sessile.inverts) <- dat.wide.epilithic.sessile.inverts$treatment  
remove(dat.wide.epilithic.sessile.inverts.comm)


# ------------------------------------------------------------------------------------
# ENDOLITHIC SESSILE INVERTS (CLAMS)
dat.wide.mobile.invert.carnivores <- dat %>%
  dplyr::filter(biomass.guild == "mobile.invert.carnivore") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.mobile.invert.carnivores.comm <- dat.wide.mobile.invert.carnivores %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.mobile.invert.carnivores.comm[is.na(dat.wide.mobile.invert.carnivores.comm)] <- 0

dat.wide.mobile.invert.carnivores <- dat.wide.mobile.invert.carnivores %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.mobile.invert.carnivores.comm)

xtabs(~ site + year, data = dat.wide.mobile.invert.carnivores)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.mobile.invert.carnivores <- dat.wide.mobile.invert.carnivores.comm
rownames(comm.mat.mobile.invert.carnivores) <- dat.wide.mobile.invert.carnivores$treatment  
remove(dat.wide.mobile.invert.carnivores.comm)

# ------------------------------------------------------------------------------------
# MOBILE INVERTEBRATE GRAZERS & DETRITIVORES

dat.wide.mobile.invert.grazers.detritivores <- dat %>%
  dplyr::filter(diversity.guild == "mobile.invert.grazer.detritivore") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.mobile.invert.grazers.detritivores.comm <- dat.wide.mobile.invert.grazers.detritivores %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.mobile.invert.grazers.detritivores.comm[is.na(dat.wide.mobile.invert.grazers.detritivores.comm)] <- 0

dat.wide.mobile.invert.grazers.detritivores <- dat.wide.mobile.invert.grazers.detritivores %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.mobile.invert.grazers.detritivores.comm)

xtabs(~ site + year, data = dat.wide.mobile.invert.grazers.detritivores)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.mobile.invert.grazers.detritivores <- dat.wide.mobile.invert.grazers.detritivores.comm
rownames(comm.mat.mobile.invert.grazers.detritivores) <- dat.wide.mobile.invert.grazers.detritivores$treatment  
remove(dat.wide.mobile.invert.grazers.detritivores.comm)

# ------------------------------------------------------------------------------------
# MOBILE INVERTEBRATE CARNIVORES

dat.wide.mobile.invert.carnivores <- dat %>%
  dplyr::filter(diversity.guild == "mobile.invert.carnivore") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.mobile.invert.carnivores.comm <- dat.wide.mobile.invert.carnivores %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.mobile.invert.carnivores.comm[is.na(dat.wide.mobile.invert.carnivores.comm)] <- 0

dat.wide.mobile.invert.carnivores <- dat.wide.mobile.invert.carnivores %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.mobile.invert.carnivores.comm)

xtabs(~ site + year, data = dat.wide.mobile.invert.carnivores)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.mobile.invert.carnivores <- dat.wide.mobile.invert.carnivores.comm
rownames(comm.mat.mobile.invert.carnivores) <- dat.wide.mobile.invert.carnivores$treatment  
remove(dat.wide.mobile.invert.carnivores.comm)

# ------------------------------------------------------------------------------------
# FISHES
dat.wide.fish <- dat %>%
  dplyr::filter(diversity.guild == "fish") %>%
  dplyr::select(-(taxon.genus:scientific.name), -biomass.guild, -diversity.guild) %>%
  spread(key = sp.code,  value = biomass, fill = NA)

# Replace NAs with zeros because data were not fully propagated with all species at all sites and times
dat.wide.fish.comm <- dat.wide.fish %>%
  dplyr::select(-(id:plot)) %>%
  as.matrix()

dat.wide.fish.comm[is.na(dat.wide.fish.comm)] <- 0

dat.wide.fish <- dat.wide.fish %>%
  dplyr::select(id:plot) %>%
  cbind(., dat.wide.fish.comm)

xtabs(~ site + year, data = dat.wide.fish)

# Rename community matrix and make row names disturbance frequency 'treatment'
comm.mat.fish <- dat.wide.fish.comm
rownames(comm.mat.fish) <- dat.wide.fish$treatment  
remove(dat.wide.fish.comm)


# ====================================================================================
# CALCULATE BIOMASS FOR EACH BIOMASS GUILD 
# ====================================================================================

biomass.dat.long <- dat %>%
  group_by(id, year, qt, time.yrs, yr.qt, time.since.start, site, treatment, plot, biomass.guild) %>%
  dplyr::summarise(biomass = sum(biomass, na.rm = TRUE)) %>%
  ungroup()

biomass.dat.wide <- spread(biomass.dat.long, key = biomass.guild, value = biomass)

# ====================================================================================
# PLOTTING
# ====================================================================================

ggplot(data = biomass.dat.wide, aes(x = treatment, y = giant.kelp, fill = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1, width = 0, fun.args = list(B = 10000)) +
  stat_summary(geom = "point", fun.y = "mean", size = 11.25) +
  stat_summary(geom = "point", fun.y = "mean", shape = 21, size = 10) +
  ylab(expression(paste("Giant kelp biomass (kg dry/", "m"^2, ")", sep = ""))) +
  guides(fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "cm"),
        aspect.ratio = 0.8) #0.4) 


# ====================================================================================
# CALCULATE DIVERSITY METRICS FOR EACH DIVERSITY GUILD
# ====================================================================================
# Species richness
dat.wide.full$richness <- specnumber(comm.mat.full)
dat.wide.understory$richness <- specnumber(comm.mat.understory)
dat.wide.sessile.inverts$richness <- specnumber(comm.mat.sessile.inverts)
dat.wide.mobile.invert.grazers.detritivores$richness <- specnumber(comm.mat.mobile.invert.grazers.detritivores)
dat.wide.mobile.invert.carnivores$richness <- specnumber(comm.mat.mobile.invert.carnivores)
dat.wide.fish$richness <- specnumber(comm.mat.fish)

dat.wide.full$initial.richness <- specnumber(comm.mat.full)
dat.wide.understory$richness <- specnumber(comm.mat.understory)
dat.wide.sessile.inverts$richness <- specnumber(comm.mat.sessile.inverts)
dat.wide.mobile.invert.grazers.detritivores$richness <- specnumber(comm.mat.mobile.invert.grazers.detritivores)
dat.wide.mobile.invert.carnivores$richness <- specnumber(comm.mat.mobile.invert.carnivores)
dat.wide.fish$richness <- specnumber(comm.mat.fish)

# Shannon diversity
dat.wide.full$shannon <- diversity(comm.mat.full, index = "shannon")
dat.wide.understory$shannon <- diversity(comm.mat.understory, index = "shannon")
dat.wide.sessile.inverts$shannon <- diversity(comm.mat.sessile.inverts, index = "shannon")
dat.wide.mobile.invert.grazers.detritivores$shannon <- diversity(comm.mat.mobile.invert.grazers.detritivores, index = "shannon")
dat.wide.mobile.invert.carnivores$shannon <- diversity(comm.mat.mobile.invert.carnivores, index = "shannon")
dat.wide.fish$shannon <- diversity(comm.mat.fish, index = "shannon")

# Species evenness (Pielou's J)
dat.wide.full$evenness <- dat.wide.full$shannon / log(dat.wide.full$richness)
dat.wide.understory$evenness <- dat.wide.understory$shannon / log(dat.wide.understory$richness)
dat.wide.sessile.inverts$evenness <- dat.wide.sessile.inverts$shannon / log(dat.wide.sessile.inverts$richness)
dat.wide.mobile.invert.grazers.detritivores$evenness <- dat.wide.mobile.invert.grazers.detritivores$shannon / log(dat.wide.mobile.invert.grazers.detritivores$richness)
dat.wide.mobile.invert.carnivores$evenness <-dat.wide.mobile.invert.carnivores$shannon / log(dat.wide.mobile.invert.carnivores$richness)
dat.wide.fish$evenness <- dat.wide.fish$shannon / log(dat.wide.fish$richness)

# Jost diversity
dat.wide.full$jost.d <- apply(X = comm.mat.full, MARGIN = 1, FUN = function(x){d(x, lev = "alpha")})
dat.wide.understory$jost.d <- apply(X = comm.mat.understory, MARGIN = 1, FUN = function(x){d(x, lev = "alpha")})
dat.wide.sessile.inverts$jost.d <- apply(X = comm.mat.sessile.inverts, MARGIN = 1, FUN = function(x){d(x, lev = "alpha")})
dat.wide.mobile.invert.grazers.detritivores$jost.d <- apply(X = comm.mat.mobile.invert.grazers.detritivores, MARGIN = 1, FUN = function(x){d(x, lev = "alpha")})
dat.wide.mobile.invert.carnivores$jost.d <- apply(X = comm.mat.mobile.invert.carnivores, MARGIN = 1, FUN = function(x){d(x, lev = "alpha")})
dat.wide.fish$jost.d <- apply(X = comm.mat.fish, MARGIN = 1, FUN = function(x){d(x, lev = "alpha")})

# ------------------------------------------------------------------------------------
# Bring richness and evenness data together into two dataframes
richness.dat.wide <- dat.wide.understory %>%
  dplyr::select(id:plot, understory = richness) %>%
  left_join(dat.wide.sessile.inverts[, c("id", "richness")], by = "id") %>%
  dplyr::rename(sessile.inverts = richness) %>%
  left_join(dat.wide.mobile.invert.grazers.detritivores[, c("id", "richness")], by = "id") %>%
  dplyr::rename(mobile.invert.grazers.detritivores = richness) %>%
  left_join(dat.wide.mobile.invert.carnivores[, c("id", "richness")], by = "id") %>%
  dplyr::rename(mobile.invert.carnivores = richness) %>%
  left_join(dat.wide.fish[, c("id", "richness")], by = "id") %>%
  dplyr::rename(fish = richness) %>%
  left_join(dat.wide.full[, c("id", "richness")], by = "id") %>%
  dplyr::rename(total = richness)

evenness.dat.wide <- dat.wide.understory %>%
  dplyr::select(id:plot, understory = evenness) %>%
  left_join(dat.wide.sessile.inverts[, c("id", "evenness")], by = "id") %>%
  dplyr::rename(sessile.inverts = evenness) %>%
  left_join(dat.wide.mobile.invert.grazers.detritivores[, c("id", "evenness")], by = "id") %>%
  dplyr::rename(mobile.invert.grazers.detritivores = evenness) %>%
  left_join(dat.wide.mobile.invert.carnivores[, c("id", "evenness")], by = "id") %>%
  dplyr::rename(mobile.invert.carnivores = evenness) %>%
  left_join(dat.wide.fish[, c("id", "evenness")], by = "id") %>%
  dplyr::rename(fish = evenness) %>%
  left_join(dat.wide.full[, c("id", "evenness")], by = "id") %>%
  dplyr::rename(total = evenness)

# ------------------------------------------------------------------------------------
# Summarize by guild / group
div.summary.fun <- function(df) {
  output <- df %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(mean.richness = mean(richness),
                     mean.shannon = mean(shannon),
                     mean.evenness = mean(evenness, na.rm = TRUE),
                     mean.jost.d = mean(jost.d),
                     se.richness = sd(richness) / sqrt(length(richness)),
                     se.shannon = sd(shannon) / sqrt(length(shannon)),
                     se.evenness = sd(evenness, na.rm = TRUE) / sqrt(length(na.omit(evenness))),
                     se.jost.d = sd(jost.d) / sqrt(length(jost.d)),
                     lower.richness = boot.fun(richness)[1],
                     upper.richness = boot.fun(richness)[2],
                     lower.shannon = boot.fun(shannon)[1],
                     upper.shannon = boot.fun(shannon)[2],
                     lower.evenness = boot.fun(na.omit(evenness))[1],
                     upper.evenness = boot.fun(na.omit(evenness))[2],
                     lower.jost.d = boot.fun(jost.d)[1],
                     upper.jost.d = boot.fun(jost.d)[2]) %>%
    dplyr::mutate(group = NA) %>%
    ungroup()
  return(output)
}

div.list <- list(dat.wide.full, dat.wide.understory, dat.wide.sessile.inverts, 
                 dat.wide.mobile.invert.grazers.detritivores, dat.wide.mobile.invert.carnivores, dat.wide.fish)
div.treatment.summary <- do.call(rbind, lapply(div.list, function(x) div.summary.fun(x)))

div.treatment.summary$group <- rep(c("Community\ntotal", "Understory\nalgae", "Sessile\ninvertebrates", "Grazers and\ndetritivores", "Carnivores", "Fishes"), each = length(unique(dat$treatment)))
div.treatment.summary$group <- factor(div.treatment.summary$group, levels = c("Understory\nalgae", "Sessile\ninvertebrates", "Grazers and\ndetritivores", "Carnivores", "Fishes", "Community\ntotal"))
div.treatment.summary$treatment <- factor(div.treatment.summary$treatment, levels = unique(div.treatment.summary$treatment))

div.treatment.summary <- div.treatment.summary %>%
  dplyr::select(group, treatment, everything())


# ====================================================================================
# PLOT DATA - RICHNESS AS A FUNCTION OF TREATMENT 
# ====================================================================================

p.base <- ggplot(data = richness.dat.wide %>% filter(year > 2016), aes(x = treatment, y = understory)) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), 
                                   color = "black", size = 26, angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), 
                                   color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = "none", fill = "none") 

p.base +
  ylab("Species richness") +
  stat_summary(aes(y = understory), fun.data = "mean_cl_boot", geom = "errorbar", size = 1, width = 0, fun.args = list(B = 10000)) +
  stat_summary(aes(y = understory), geom = "point", fun.y = "mean", size = 11.25) +
  stat_summary(aes(y = understory, fill = treatment), geom = "point", fun.y = "mean", shape = 21, size = 10) +
  ggtitle("Understory algae") +
  coord_cartesian(ylim = c(8, 11))
  
p.base +
  ylab("Species richness") +
  stat_summary(aes(y = sessile.inverts), fun.data = "mean_cl_boot", geom = "errorbar", size = 1, width = 0, fun.args = list(B = 10000)) +
  stat_summary(aes(y = sessile.inverts), geom = "point", fun.y = "mean", size = 11.25) +
  stat_summary(aes(y = sessile.inverts, fill = treatment), geom = "point", fun.y = "mean", shape = 21, size = 10) +
  ggtitle("Sessile invertebrates") +
  coord_cartesian(ylim = c(10, 16))

p.base +
  ylab("Species richness") +
  stat_summary(aes(y = mobile.invert.grazers.detritivores), fun.data = "mean_cl_boot", geom = "errorbar", size = 1, width = 0, fun.args = list(B = 10000)) +
  stat_summary(aes(y = mobile.invert.grazers.detritivores), geom = "point", fun.y = "mean", size = 11.25) +
  stat_summary(aes(y = mobile.invert.grazers.detritivores, fill = treatment), geom = "point", fun.y = "mean", shape = 21, size = 10) +
  ggtitle("Mobile invertebrate grazers and detritivores") +
  coord_cartesian(ylim = c(4, 6)) 

p.base +
  ylab("Species richness") +
  stat_summary(aes(y = mobile.invert.carnivores), fun.data = "mean_cl_boot", geom = "errorbar", size = 1, width = 0, fun.args = list(B = 10000)) +
  stat_summary(aes(y = mobile.invert.carnivores), geom = "point", fun.y = "mean", size = 11.25) +
  stat_summary(aes(y = mobile.invert.carnivores, fill = treatment), geom = "point", fun.y = "mean", shape = 21, size = 10) +
  ggtitle("Mobile invertebrate carnivores") +
  coord_cartesian(ylim = c(3, 5)) 

p.base +
  ylab("Species richness") +
  stat_summary(aes(y = fish), fun.data = "mean_cl_boot", geom = "errorbar", size = 1, width = 0, fun.args = list(B = 10000)) +
  stat_summary(aes(y = fish), geom = "point", fun.y = "mean", size = 11.25) +
  stat_summary(aes(y = fish, fill = treatment), geom = "point", fun.y = "mean", shape = 21, size = 10) +
  ggtitle("Fish") +
  coord_cartesian(ylim = c(4, 8))


# ============================================================================================================
# PLOT GUILD BIOMASS BY TREATMENT 
# ============================================================================================================

# -------------------------------------------------------------------------------------------------------------------
# Giant kelp

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = giant.kelp)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = giant.kelp)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = giant.kelp)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1)
p4

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = giant.kelp)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment)
p5

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide, aes(x = year, y = giant.kelp, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment)
p6

# -------------------------------------------------------------------------------------------------------------------
# Understory algae

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = algae)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = algae)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide %>% filter(year > 2016), aes(x = time.since.start, y = algae)) +
  # stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) +
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1)
p4

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide %>% filter(year > 2016), aes(x = time.since.start, y = algae)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment)
p5

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide %>% filter(year > 2016), aes(x = year, y = algae, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment)
p6


# -------------------------------------------------------------------------------------------------------------------
# Endolithic sessile inverts (clams)

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = mobile.invert.carnivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = mobile.invert.carnivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = mobile.invert.carnivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1)
p4

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = mobile.invert.carnivore)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p5

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide, aes(x = year, y = mobile.invert.carnivore, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p6



# -------------------------------------------------------------------------------------------------------------------
# Epilithic sessile inverts

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = epilithic.sessile.invert)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = epilithic.sessile.invert)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = epilithic.sessile.invert)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1)
p4

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = epilithic.sessile.invert)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p5

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide, aes(x = year, y = epilithic.sessile.invert, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p6


# -------------------------------------------------------------------------------------------------------------------
# Mobile invertebrate grazers and detritivores

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = mobile.invert.grazer.detritivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = mobile.invert.grazer.detritivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = mobile.invert.grazer.detritivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1)
p4

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = mobile.invert.grazer.detritivore)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p5

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide, aes(x = year, y = mobile.invert.grazer.detritivore, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p6

# -------------------------------------------------------------------------------------------------------------------
# Mobile invertebrate carnivores

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = mobile.invert.carnivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = mobile.invert.carnivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = mobile.invert.carnivore)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1)
p4

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = mobile.invert.carnivore)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p5

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide, aes(x = year, y = mobile.invert.carnivore, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p6


# -------------------------------------------------------------------------------------------------------------------
# Fish

# Plot mean by treatment
p1 <- ggplot(data = biomass.dat.wide, aes(x = treatment, y = fish)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Treatment") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  scale_x_discrete(labels = c("Control", "Annual", "Continual"), name = "") +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) 
p1

# Facet by season
p2 <- p1 + facet_wrap(~ qt, nrow = 1)
p2

# Facet by treatment
p3 <-  ggplot(data = biomass.dat.wide, aes(x = qt, y = fish)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 10) +
  xlab("Season") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1) +
  facet_wrap(~ treatment, ncol = 3)
p3

# Plot means over time
p4 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = fish)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_wrap(~treatment, ncol = 1) 
p4
# NOTE: LOOKS LIKE SOME OUTLIERS MAY NEED TO BE REMOVED

# Plot all data over time, faceted by site
p5 <- ggplot(data = biomass.dat.wide, aes(x = time.since.start, y = fish)) +
  geom_line() +
  geom_point(aes(fill = treatment), shape = 21, size = 2) +
  xlab("Years since start") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 26),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 26),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 26),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p5
# NOTE: LOOKS LIKE SOME OUTLIERS MAY NEED TO BE REMOVED

# Plot annual means over time, faceted by site
p6 <- ggplot(data = biomass.dat.wide, aes(x = year, y = fish, 
                                          group = paste(treatment, site, year))) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "black", width = 0, size = 0.5) + 
  stat_summary(aes(fill = treatment), fun = mean, geom = "point", shape = 21, size = 2) +
  scale_x_continuous(breaks = seq(2008, 2022, by = 4)) +
  xlab("Year") +
  ylab(expression(paste("Biomass (g dry/", "m"^2, ")", sep = ""))) +
  guides(color = FALSE, fill = FALSE) +
  theme(plot.title = element_text(size = 26, margin = margin(b = 20)),
        panel.border = element_rect(color="black", fill = NA, size=2),
        panel.background = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        panel.spacing = unit(0.5,"lines"),
        text = element_text(size = 15),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 15, 
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 15),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.ticks = element_line(size=1),
        axis.ticks.length = unit(0.4, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.5, "line"),
        aspect.ratio = 1/2) +
  facet_grid(site ~ treatment, scales = 'free_y')
p6
# NOTE: LOOKS LIKE SOME OUTLIERS MAY NEED TO BE REMOVED
