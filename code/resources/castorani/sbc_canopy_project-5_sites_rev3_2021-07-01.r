# ====================================================================================
# SBC LTER KELP CANOPY PROJECT
# ====================================================================================

# Max Castorani, University of Virginia, castorani@virginia.edu
# Shannon Harrer, University of California, Santa Barbara, harrer@ucsb.edu
# Dan Reed, University of California, Santa Barbara, danreed@ucsb.edu
# Robert Miller, University of California, Santa Barbara, rjmiller@ucsb.edu

# Revised July 1, 2021


# ====================================================================================
# LOAD LIBRARIES
# ====================================================================================

# ------------------------------------------------------------------------------------
# Clear environment
rm(list = ls())

# ------------------------------------------------------------------------------------
# Check for and install required packages

library('pacman')
library('plyr')
#devtools::install_github("easystats/easystats")
#devtools::install_github("easystats/report")
library('easystats')
library('report')

pacman::p_load('tidyverse', 'ggplot2', 'nlme', 'car', 'lme4',  'reshape', 'multcomp', 'segmented', 'merTools', 'GGally',
                  'ciTools', 'effects', 'here', 'RColorBrewer', 'ggeffects', 'emmeans', 'r2glmm', 'glmmTMB', 'bbmle',
                   'vegan', 'vegetarian', 'DHARMa', 'mgcv', 'lubridate', 'patchwork', 'lmerTest', 'tidymv', #'itsadug',
               'performance', 'see', 'svglite')

# ------------------------------------------------------------------------------------
# Use 'here' to tell R where your data are
here::here()
setwd(here::here())

# ------------------------------------------------------------------------------------
# Turn off scientific notation
options(scipen = 999)

# ------------------------------------------------------------------------------------
# Set custom colors
# col1 <- "#0868ac"
# col2 <- "#a1dab4"
# col3 <- "#ffffcc"

col1 <- "#6CC5AD" 
col2 <- "#4E96BF"
col3 <- "#EEB07F"

col4 <- "#065389"
col5 <- "#6ac488" 
col6 <- "#99994c"


# ====================================================================================
# IMPORT AND TIDY DATA
# ====================================================================================

# Light data
# Note: Benthic light data represent mean daily light (mol m-2 d-1) for each season within LTE transects

# ------------------------------------------------------------------------------------
# Package ID: knb-lter-sbc.36.18 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Kelp Removal Experiment: Hourly photon irradiance at the surface and seafloor.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Data set creator:  Shannon Harrer -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/36/18/d6e4640c963f78393c6f607601064f86" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "SITE",     
                 "TRANSECT",     
                 "SENSOR_LOCATION",     
                 "DATE_LOCAL",     
                 "TIME_LOCAL",     
                 "LIGHT_UMOL",     
                 "LIGHT_ADJ_CONV"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$SENSOR_LOCATION)!="factor") dt1$SENSOR_LOCATION<- as.factor(dt1$SENSOR_LOCATION)                                   
# attempting to convert dt1$DATE_LOCAL dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1DATE_LOCAL<-as.Date(dt1$DATE_LOCAL,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1DATE_LOCAL) == length(tmp1DATE_LOCAL[!is.na(tmp1DATE_LOCAL)])){dt1$DATE_LOCAL <- tmp1DATE_LOCAL } else {print("Date conversion failed for dt1$DATE_LOCAL. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1DATE_LOCAL) 
if (class(dt1$LIGHT_UMOL)=="factor") dt1$LIGHT_UMOL <-as.numeric(levels(dt1$LIGHT_UMOL))[as.integer(dt1$LIGHT_UMOL) ]               
if (class(dt1$LIGHT_UMOL)=="character") dt1$LIGHT_UMOL <-as.numeric(dt1$LIGHT_UMOL)
if (class(dt1$LIGHT_ADJ_CONV)=="factor") dt1$LIGHT_ADJ_CONV <-as.numeric(levels(dt1$LIGHT_ADJ_CONV))[as.integer(dt1$LIGHT_ADJ_CONV) ]               
if (class(dt1$LIGHT_ADJ_CONV)=="character") dt1$LIGHT_ADJ_CONV <-as.numeric(dt1$LIGHT_ADJ_CONV)

# Convert Missing Values to NA for non-dates

dt1$LIGHT_UMOL <- ifelse((trimws(as.character(dt1$LIGHT_UMOL))==trimws("-99999")),NA,dt1$LIGHT_UMOL)               
suppressWarnings(dt1$LIGHT_UMOL <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$LIGHT_UMOL))==as.character(as.numeric("-99999"))),NA,dt1$LIGHT_UMOL))
dt1$LIGHT_ADJ_CONV <- ifelse((trimws(as.character(dt1$LIGHT_ADJ_CONV))==trimws("-99999")),NA,dt1$LIGHT_ADJ_CONV)               
suppressWarnings(dt1$LIGHT_ADJ_CONV <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$LIGHT_ADJ_CONV))==as.character(as.numeric("-99999"))),NA,dt1$LIGHT_ADJ_CONV))

#Calculate average daily irradiance
light.dat.transect <-dt1 %>%
   dplyr::filter(!is.na(LIGHT_ADJ_CONV)) %>%
   dplyr::mutate(PAR_D=LIGHT_ADJ_CONV*3600/1000000) %>%
   dplyr::group_by(SITE,TRANSECT,SENSOR_LOCATION,DATE_LOCAL) %>%
   dplyr::summarise(LIGHT_MOL_DAY=sum(PAR_D),HR_DAILY=n()) %>%
   dplyr::ungroup() %>%
   #choose the dates with full-day data
   dplyr::filter(HR_DAILY==24) %>%
   dplyr::select(-HR_DAILY) %>%
   dplyr::arrange(SITE,TRANSECT,SENSOR_LOCATION,DATE_LOCAL) %>% 
   mutate(DOY =yday(DATE_LOCAL), year=substr(DATE_LOCAL,1,4),
          season = case_when (DOY>=355|DOY<79~"1-Winter",
                              DOY>=79&DOY<172~"2-Spring",
                              DOY>=172&DOY<265~"3-Summer",
                              DOY>=265&DOY<355~"4-Autumn") ) %>%
   dplyr::select(-DOY) %>%
   dplyr::filter(SENSOR_LOCATION == 'SEAFLOOR') %>%
   dplyr::filter(TRANSECT != "MKO" & TRANSECT != "MKI" & TRANSECT !="1") %>%
   dplyr::group_by(SITE,TRANSECT,year,season) %>%
   dplyr::summarise(Avg_seafloor_mol=mean(LIGHT_MOL_DAY),seafloor_n=n()) %>%
   dplyr::ungroup()

#Calculate hourly daily irradiance
light.dat.transect.daily <-dt1 %>%
  dplyr::filter(!is.na(LIGHT_ADJ_CONV)) %>%
  dplyr::mutate(PAR_D=LIGHT_ADJ_CONV*3600/1000000) %>%
  dplyr::group_by(SITE,TRANSECT,SENSOR_LOCATION,DATE_LOCAL) %>%
  dplyr::summarise(LIGHT_MOL_DAY=sum(PAR_D),HR_DAILY=n()) %>%
  dplyr::ungroup() %>%
  #choose the dates with full-day data
  dplyr::filter(HR_DAILY==24) %>%
  dplyr::select(-HR_DAILY) %>%
  dplyr::arrange(SITE,TRANSECT,SENSOR_LOCATION,DATE_LOCAL) %>% 
  mutate(DOY =yday(DATE_LOCAL), year=substr(DATE_LOCAL,1,4),
         season = case_when (DOY>=355|DOY<79~"1-Winter",
                             DOY>=79&DOY<172~"2-Spring",
                             DOY>=172&DOY<265~"3-Summer",
                             DOY>=265&DOY<355~"4-Autumn") ) %>%
  dplyr::filter(SENSOR_LOCATION == 'SEAFLOOR') %>%
  dplyr::filter(TRANSECT != "MKO" & TRANSECT != "MKI" & TRANSECT !="1") 

rm(dt1, infile1, inUrl1)

# ------------------------------------------------------------------------------------
# Light data QAQC
n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("n = ",length(x))))
}

ggplot(data = light.dat.transect, aes(x = season, y = Avg_seafloor_mol)) +
  geom_boxplot() +
  stat_summary(fun.data = n_fun, geom = "text", vjust = -1)

# ------------------------------------------------------------------------------------
# Depth data

# Package ID: knb-lter-sbc.44.8 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Long-term experiment: Kelp Removal: Transect depth data .
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/44/8/b33eb7e93b097b92adc48bd802451fb4" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "SITE",     
                 "TRANSECT",     
                 "DEPTH_MLLW_M",     
                 "SD_DEPTH",     
                 "CV_DEPTH",     
                 "LATITUDE",     
                 "LONGITUDE",     
                 "SITE_NAME"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$DEPTH_MLLW_M)=="factor") dt1$DEPTH_MLLW_M <-as.numeric(levels(dt1$DEPTH_MLLW_M))[as.integer(dt1$DEPTH_MLLW_M) ]               
if (class(dt1$DEPTH_MLLW_M)=="character") dt1$DEPTH_MLLW_M <-as.numeric(dt1$DEPTH_MLLW_M)
if (class(dt1$SD_DEPTH)=="factor") dt1$SD_DEPTH <-as.numeric(levels(dt1$SD_DEPTH))[as.integer(dt1$SD_DEPTH) ]               
if (class(dt1$SD_DEPTH)=="character") dt1$SD_DEPTH <-as.numeric(dt1$SD_DEPTH)
if (class(dt1$CV_DEPTH)=="factor") dt1$CV_DEPTH <-as.numeric(levels(dt1$CV_DEPTH))[as.integer(dt1$CV_DEPTH) ]               
if (class(dt1$CV_DEPTH)=="character") dt1$CV_DEPTH <-as.numeric(dt1$CV_DEPTH)
if (class(dt1$LATITUDE)=="factor") dt1$LATITUDE <-as.numeric(levels(dt1$LATITUDE))[as.integer(dt1$LATITUDE) ]               
if (class(dt1$LATITUDE)=="character") dt1$LATITUDE <-as.numeric(dt1$LATITUDE)
if (class(dt1$LONGITUDE)=="factor") dt1$LONGITUDE <-as.numeric(levels(dt1$LONGITUDE))[as.integer(dt1$LONGITUDE) ]               
if (class(dt1$LONGITUDE)=="character") dt1$LONGITUDE <-as.numeric(dt1$LONGITUDE)
if (class(dt1$SITE_NAME)!="factor") dt1$SITE_NAME<- as.factor(dt1$SITE_NAME)

# Convert Missing Values to NA for non-dates

depth.dat.transect <- dt1 %>%
  dplyr::select(SITE, TRANSECT, DEPTH_MLLW_M)

rm(dt1, infile1, inUrl1)

# Check data
ggplot(data = depth.dat.transect, aes(x = factor(TRANSECT), y = -DEPTH_MLLW_M)) +
  geom_point(size = 2) +
  facet_wrap(~ SITE, scales = "free_x") +
  ylim(-10, 0) +
  geom_hline(yintercept = 0) +
  guides(fill = F)

# ------------------------------------------------------------------------------------
# Sand percent cover data

# Package ID: knb-lter-sbc.28.28 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Long-term experiment: Kelp removal: Cover of sessile organisms, Uniform Point Contact.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/28/28/1a0798e8f4f37d2709419b460ff1e61d" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "YEAR",     
                 "MONTH",     
                 "DATE",     
                 "SITE",     
                 "TRANSECT",     
                 "TREATMENT",     
                 "QUAD",     
                 "SIDE",     
                 "SUBSTRATE_TYPE",     
                 "PERCENT_COVER"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$MONTH)!="factor") dt2$MONTH<- as.factor(dt2$MONTH)                                   
# attempting to convert dt2$DATE dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp2DATE<-as.Date(dt2$DATE,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2DATE) == length(tmp2DATE[!is.na(tmp2DATE)])){dt2$DATE <- tmp2DATE } else {print("Date conversion failed for dt2$DATE. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2DATE) 
if (class(dt2$SITE)!="factor") dt2$SITE<- as.factor(dt2$SITE)
if (class(dt2$TRANSECT)!="factor") dt2$TRANSECT<- as.factor(dt2$TRANSECT)
if (class(dt2$TREATMENT)!="factor") dt2$TREATMENT<- as.factor(dt2$TREATMENT)
if (class(dt2$QUAD)!="factor") dt2$QUAD<- as.factor(dt2$QUAD)
if (class(dt2$SIDE)!="factor") dt2$SIDE<- as.factor(dt2$SIDE)
if (class(dt2$SUBSTRATE_TYPE)!="factor") dt2$SUBSTRATE_TYPE<- as.factor(dt2$SUBSTRATE_TYPE)
if (class(dt2$PERCENT_COVER)=="factor") dt2$PERCENT_COVER <-as.numeric(levels(dt2$PERCENT_COVER))[as.integer(dt2$PERCENT_COVER) ]               
if (class(dt2$PERCENT_COVER)=="character") dt2$PERCENT_COVER <-as.numeric(dt2$PERCENT_COVER)

# Convert Missing Values to NA for non-dates

dt2$MONTH <- as.factor(ifelse((trimws(as.character(dt2$MONTH))==trimws("-99999")),NA,as.character(dt2$MONTH)))
dt2$SITE <- as.factor(ifelse((trimws(as.character(dt2$SITE))==trimws("-99999")),NA,as.character(dt2$SITE)))
dt2$TRANSECT <- as.factor(ifelse((trimws(as.character(dt2$TRANSECT))==trimws("-99999")),NA,as.character(dt2$TRANSECT)))
dt2$QUAD <- as.factor(ifelse((trimws(as.character(dt2$QUAD))==trimws("-99999")),NA,as.character(dt2$QUAD)))
dt2$SIDE <- as.factor(ifelse((trimws(as.character(dt2$SIDE))==trimws("-99999")),NA,as.character(dt2$SIDE)))
dt2$SUBSTRATE_TYPE <- as.factor(ifelse((trimws(as.character(dt2$SUBSTRATE_TYPE))==trimws("-99999")),NA,as.character(dt2$SUBSTRATE_TYPE)))
dt2$PERCENT_COVER <- ifelse((trimws(as.character(dt2$PERCENT_COVER))==trimws("-99999")),NA,dt2$PERCENT_COVER)               
suppressWarnings(dt2$PERCENT_COVER <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt2$PERCENT_COVER))==as.character(as.numeric("-99999"))),NA,dt2$PERCENT_COVER))

# ---------------------------------------------------------------------------------------------------
# Recode bottom substrate types

# B Bedrock
# BL Boulder Large: rock greater than 1 meter in diameter
# BM Boulder Medium: rock 50 to 100 cm in diameter
# BS Boulder Small: rock 25-50 cm in diameter
# C Cobble: rock less than 25 cm in diameter
# S Sand: sand greater than 2.5 cm deep
# SS Shallow Sand: sand less than 2.5 cm deep
# SH Shell Debris, broken

sand.dat.section <- dt2 %>%
  dplyr::mutate(SUBSTRATE_TYPE = dplyr::recode(.x = SUBSTRATE_TYPE,
                                               'B' = 'bedrock',
                                               'BL' = 'boulder large',
                                               'BM' = 'boulder medium',
                                               'BS' = 'boulder small',
                                               'C' = 'cobble',
                                               'S' = 'sand',
                                               'SH' = 'shell debris',
                                               'SS' = 'shallow sand')) %>%

  # Spread and gather data to propagate each class of substrate across all sections
  tidyr::pivot_wider(id_cols = c(YEAR, MONTH, DATE, SITE, TRANSECT, TREATMENT, QUAD, SIDE),
                     names_from = SUBSTRATE_TYPE,
                     values_from = PERCENT_COVER,
                     values_fill = list(PERCENT_COVER = 0)) %>%
  tidyr::pivot_longer(-c(YEAR:SIDE),
                      names_to = "SUBSTRATE_TYPE",
                      values_to = "PERCENT_COVER") %>%
  #dplyr::filter(SUBSTRATE_TYPE == 'sand') %>%
  dplyr::filter(SUBSTRATE_TYPE == 'sand' | SUBSTRATE_TYPE == 'shell debris') %>%
  #dplyr::filter(SUBSTRATE_TYPE == 'sand' | SUBSTRATE_TYPE == 'shell debris' | SUBSTRATE_TYPE == 'cobble') %>%
  #dplyr::filter(SUBSTRATE_TYPE == 'bedrock' | SUBSTRATE_TYPE == 'boulder large' | SUBSTRATE_TYPE == 'boulder medium' | SUBSTRATE_TYPE == 'boulder small') %>%
  dplyr::group_by(YEAR, MONTH, DATE, SITE, TRANSECT, TREATMENT, QUAD, SIDE) %>%
  dplyr::summarise(percent.sand = sum(PERCENT_COVER, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(DATE = as.Date(as.character(DATE)))

sand.dat.transect <- sand.dat.section %>%
  dplyr::group_by(YEAR, MONTH, DATE, SITE, TRANSECT, TREATMENT) %>%
  dplyr::summarise(percent.sand = mean(percent.sand, na.rm = T)) %>%
  dplyr::ungroup()

temp <- sand.dat.section %>%
  dplyr::group_by(YEAR, SITE) %>%
  dplyr::summarise(percent.sand = mean(percent.sand, na.rm = T)) %>%
  dplyr::ungroup()

ggplot(data = temp, aes(x = YEAR, y = percent.sand)) +
  geom_line() +
  facet_wrap(~ SITE) +
  ggtitle("Annual mean cover of sand and shallow sand") +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2))

rm(temp, dt2, infile2, inUrl2)

# ------------------------------------------------------------------------------------
# Biomass data

# Package ID: knb-lter-sbc.119.6 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Long-term experiment: biomass of kelp forest species, ongoing since 2008.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/119/6/020f9b0e507561cb49fc8fd122da0a29" 
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

biomass.dat.transect <- dt1
rm(dt1, infile1, inUrl1)

# ------------------------------------------------------------------------------------
# Isolate the sea urchin biomass data at transect scale
urchin.dat.transect <- biomass.dat.transect %>%
  dplyr::filter(COMMON_NAME == 'Red Urchin' | 
                COMMON_NAME == 'Purple Urchin' |
                COMMON_NAME == 'Crowned Sea Urchin' |
                COMMON_NAME == 'White Sea Urchin') %>%
  dplyr::select(YEAR, MONTH, DATE, SITE, TRANSECT, TREATMENT, SCIENTIFIC_NAME, DENSITY, DRY_GM2, SFDM) %>%
  
  # Sum the total number of urchins
  dplyr::group_by(YEAR, MONTH, DATE, SITE, TRANSECT, TREATMENT) %>%
  dplyr::summarise(urchin.dens.no.m2 = sum(DENSITY, na.rm = T),
                   urchin.dry.gm2    = sum(DRY_GM2, na.rm = T),
                   urchin.sfdm.gm2   = sum(SFDM, na.rm = T)) %>%
  dplyr::ungroup()

temp <- urchin.dat.transect %>%
  dplyr::group_by(YEAR, SITE, TREATMENT) %>%
  dplyr::summarise(mean.urchin.dry.gm2 = mean(urchin.dry.gm2, na.rm = T),
                   mean.urchin.dens.no.m2 = mean(urchin.dens.no.m2, na.rm = T)) %>%
  dplyr::ungroup()

ggplot(data = temp, aes(x = YEAR, y = mean.urchin.dens.no.m2)) +
  geom_line() +
  facet_grid(TREATMENT ~ SITE) +
  ggtitle("Annual mean density of sea urchins") +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2))

ggplot(data = temp, aes(x = YEAR, y = mean.urchin.dry.gm2)) +
  geom_line() +
  facet_grid(TREATMENT ~ SITE) +
  ggtitle("Annual mean dry mass of sea urchins") +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2))

rm(temp)

# ------------------------------------------------------------------------------------
# Import giant kelp frond density data

# Package ID: knb-lter-sbc.29.19 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Long-term experiment: Kelp removal: Abundance and size of Giant Kelp.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/29/19/83b290d3e780ebf95ecdbef239cdec48" 
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
                 "QUAD",     
                 "SIDE",     
                 "SP_CODE",     
                 "FRONDS",     
                 "AREA",     
                 "OBS_CODE",     
                 "NOTES",     
                 "SCIENTIFIC_NAME",     
                 "COMMON_NAME",     
                 "TAXON_KINGDOM",     
                 "TAXON_PHYLUM",     
                 "TAXON_CLASS",     
                 "TAXON_ORDER",     
                 "TAXON_FAMILY",     
                 "TAXON_GENUS",     
                 "GROUP",     
                 "SURVEY",     
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
if (class(dt1$QUAD)!="factor") dt1$QUAD<- as.factor(dt1$QUAD)
if (class(dt1$SIDE)!="factor") dt1$SIDE<- as.factor(dt1$SIDE)
if (class(dt1$SP_CODE)!="factor") dt1$SP_CODE<- as.factor(dt1$SP_CODE)
if (class(dt1$FRONDS)=="factor") dt1$FRONDS <-as.numeric(levels(dt1$FRONDS))[as.integer(dt1$FRONDS) ]               
if (class(dt1$FRONDS)=="character") dt1$FRONDS <-as.numeric(dt1$FRONDS)
if (class(dt1$AREA)=="factor") dt1$AREA <-as.numeric(levels(dt1$AREA))[as.integer(dt1$AREA) ]               
if (class(dt1$AREA)=="character") dt1$AREA <-as.numeric(dt1$AREA)
if (class(dt1$OBS_CODE)!="factor") dt1$OBS_CODE<- as.factor(dt1$OBS_CODE)
if (class(dt1$NOTES)!="factor") dt1$NOTES<- as.factor(dt1$NOTES)
if (class(dt1$SCIENTIFIC_NAME)!="factor") dt1$SCIENTIFIC_NAME<- as.factor(dt1$SCIENTIFIC_NAME)
if (class(dt1$COMMON_NAME)!="factor") dt1$COMMON_NAME<- as.factor(dt1$COMMON_NAME)
if (class(dt1$TAXON_KINGDOM)!="factor") dt1$TAXON_KINGDOM<- as.factor(dt1$TAXON_KINGDOM)
if (class(dt1$TAXON_PHYLUM)!="factor") dt1$TAXON_PHYLUM<- as.factor(dt1$TAXON_PHYLUM)
if (class(dt1$TAXON_CLASS)!="factor") dt1$TAXON_CLASS<- as.factor(dt1$TAXON_CLASS)
if (class(dt1$TAXON_ORDER)!="factor") dt1$TAXON_ORDER<- as.factor(dt1$TAXON_ORDER)
if (class(dt1$TAXON_FAMILY)!="factor") dt1$TAXON_FAMILY<- as.factor(dt1$TAXON_FAMILY)
if (class(dt1$TAXON_GENUS)!="factor") dt1$TAXON_GENUS<- as.factor(dt1$TAXON_GENUS)
if (class(dt1$GROUP)!="factor") dt1$GROUP<- as.factor(dt1$GROUP)
if (class(dt1$SURVEY)!="factor") dt1$SURVEY<- as.factor(dt1$SURVEY)
if (class(dt1$MOBILITY)!="factor") dt1$MOBILITY<- as.factor(dt1$MOBILITY)
if (class(dt1$GROWTH_MORPH)!="factor") dt1$GROWTH_MORPH<- as.factor(dt1$GROWTH_MORPH)

# Convert Missing Values to NA for non-dates

dt1$SIDE <- as.factor(ifelse((trimws(as.character(dt1$SIDE))==trimws("-99999")),NA,as.character(dt1$SIDE)))
dt1$SP_CODE <- as.factor(ifelse((trimws(as.character(dt1$SP_CODE))==trimws("-99999")),NA,as.character(dt1$SP_CODE)))
dt1$FRONDS <- ifelse((trimws(as.character(dt1$FRONDS))==trimws("-99999")),NA,dt1$FRONDS)               
suppressWarnings(dt1$FRONDS <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$FRONDS))==as.character(as.numeric("-99999"))),NA,dt1$FRONDS))
dt1$OBS_CODE <- as.factor(ifelse((trimws(as.character(dt1$OBS_CODE))==trimws("-99999")),NA,as.character(dt1$OBS_CODE)))
dt1$SCIENTIFIC_NAME <- as.factor(ifelse((trimws(as.character(dt1$SCIENTIFIC_NAME))==trimws("-99999")),NA,as.character(dt1$SCIENTIFIC_NAME)))


#Convert date column
dt1$DATE <- as.Date(dt1$DATE)

#Group by site, transect, date, and section
frond.dat.section <- dt1 %>% 
  dplyr::select(-SP_CODE) %>% 
  dplyr::group_by(DATE, SITE, TRANSECT, TREATMENT, QUAD, SIDE) %>%  
  dplyr::summarise(FRONDS_PER_SECTION = sum(FRONDS),
                   FRONDS_PER_M2 = sum(FRONDS)/20,                        # "FROND_DENSITY" = number of fronds per m2 in the 20 m2 swath of each section
                   FRONDS_PER_PLANT = mean(FRONDS),                       # "FRONDS_PER_PLANT" = the average number of fronds on each plant in the 20 m2 swath of each section
                   PLANTS_PER_SECTION = length(FRONDS[FRONDS != 0]),
                   PLANTS_PER_M2 = length(FRONDS[FRONDS != 0]) / 20) %>%  # "PLANTS_PER_M2" = number of plants per m2 in the 20 m2 swath of each section
  dplyr::ungroup()

# FIX ISSUE IN WHICH FRONDS WERE NOT MEASURED BUT CODE RECORDS PLANTS
frond.dat.section$PLANTS_PER_SECTION[is.na(frond.dat.section$FRONDS_PER_M2)] <- NA
frond.dat.section$PLANTS_PER_M2[is.na(frond.dat.section$FRONDS_PER_M2)] <- NA

rm(dt1, infile1, inUrl1)

frond.dat.transect <- frond.dat.section %>%
  dplyr::select(-QUAD, -SIDE) %>%
  dplyr::group_by(DATE, SITE, TRANSECT, TREATMENT) %>%
  dplyr::summarise(FRONDS_PER_TRANSECT = sum(FRONDS_PER_SECTION, na.rm = F),
                   PLANTS_PER_TRANSECT = sum(PLANTS_PER_SECTION, na.rm = F),
                   FRONDS_PER_M2 = mean(FRONDS_PER_M2, na.rm = F), 
                   FRONDS_PER_PLANT = mean(FRONDS_PER_PLANT, na.rm = F), 
                   PLANTS_PER_M2 = mean(PLANTS_PER_M2, na.rm = F)) %>% 
  dplyr::ungroup() %>%
  dplyr::rename(DATE_KELP = DATE)

# ------------------------------------------------------------------------------------
# Clean up light data
light.dat.transect.daily1 <- light.dat.transect.daily %>%
  dplyr::select(-SENSOR_LOCATION, -DOY) %>%
  dplyr::rename(YEAR = year, SEASON = season, DATE = DATE_LOCAL) %>%
  dplyr::mutate(YEAR = as.numeric(YEAR)) %>%
  droplevels()

light.dat.transect.daily2 <- light.dat.transect.daily1 %>%
  dplyr::rename(DATE_LIGHT = DATE) %>%
  dplyr::select(-SEASON, -YEAR)

# ------------------------------------------------------------------------------------
# Join kelp frond data with PAR data, and filter based on date match-ups within a range of days ("thresh")
thresh <- 7

daily.light.kelp.dat <- full_join(light.dat.transect.daily2, frond.dat.transect,
                                  by = c("SITE", "TRANSECT")) %>%
  dplyr::mutate(DELTA_DATE = as.numeric(DATE_LIGHT - DATE_KELP)) %>%
  #dplyr::filter(DELTA_DATE < 0) %>%
  #dplyr::filter(DELTA_DATE >= -thresh) %>%
  dplyr::filter(abs(DELTA_DATE) <= thresh) %>%
  dplyr::select(SITE, TREATMENT, TRANSECT, DATE_LIGHT, DATE_KELP, DELTA_DATE, LIGHT_MOL_DAY, FRONDS_PER_M2, FRONDS_PER_PLANT, PLANTS_PER_M2) 

# Filter only to control plots
daily.light.kelp.dat.all.treatments <- daily.light.kelp.dat

daily.light.kelp.dat <- daily.light.kelp.dat #%>%
 # dplyr::filter(SITE == "MOHK") %>%
 #  dplyr::filter(TREATMENT == "CONTROL" ) #| TREATMENT == "ANNUAL") 

# Summarized across window of time ("thresh")
daily.light.kelp.dat <- daily.light.kelp.dat %>%
  dplyr::group_by(SITE, TREATMENT, TRANSECT, DATE_KELP) %>%  
  dplyr::summarise_all(mean) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(YEAR = as.numeric(as.character(year(DATE_KELP)))) %>%
  dplyr::rename(DATE = DATE_KELP) %>%
  dplyr::mutate(MONTH = as.numeric(as.character(month(DATE)))) %>%
  dplyr::select(SITE, TREATMENT, TRANSECT, YEAR, MONTH, DATE, LIGHT_MOL_DAY, FRONDS_PER_M2, FRONDS_PER_PLANT, PLANTS_PER_M2)

add.season.column.to.light.dat <- function(x){
  x <- x %>%
    
    # Create column for day of YEAR (DOY)
    dplyr::mutate(doy = yday(DATE)) %>% 
    
    # Create column for quarter
    dplyr::mutate(SEASON = case_when(
      doy>=355 | doy<79  ~ "1-Winter",
      doy>=79  & doy<172 ~ "2-Spring",
      doy>=172 & doy<265 ~ "3-Summer",
      doy>=265 & doy<355 ~ "4-Autumn")) %>%
    
    mutate(TIME = YEAR + ifelse(SEASON == "1-Winter", 0.125,
                                ifelse(SEASON == "2-Spring", 0.375,
                                       ifelse(SEASON == "3-Summer", 0.625, 0.875)))) %>%
    dplyr::select(-doy) %>%
    dplyr::select(YEAR, SEASON, TIME, everything())
  
  return(x)
}

daily.light.kelp.dat <- add.season.column.to.light.dat(daily.light.kelp.dat)

ggplot(data = daily.light.kelp.dat, aes(x = sqrt(FRONDS_PER_M2), y = sqrt(LIGHT_MOL_DAY))) +
  geom_point() +
  geom_smooth(se = F, method = 'lm') +
  facet_wrap(~ SEASON, ncol = 4) +
  theme(aspect.ratio = 1) +
  xlim(0, 5) +
  ylim(0, 4) +
  ggtitle(thresh)
  #ggtitle("All sites - All plots")

# ------------------------------------------------------------------------------------
# Isolate the algal biomass data 
biomass.dat.transect <- biomass.dat.transect %>%
  dplyr::filter(GROUP == "ALGAE")

# ------------------------------------------------------------------------------------
# NPP data
# Package ID: knb-lter-sbc.58.13 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Long-term experiment: Taxon-specific seasonal net primary production (NPP) for macroalgae.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Shannon Harrer -  
# Data set creator:  Daniel C Reed -  
# Data set creator:  Robert J Miller -  
# Data set creator:  Sally J Holbrook -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/58/13/b737695757cb5f7a63d33611c2e3cc52" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "YEAR",     
                 "SEASON",     
                 "SITE",     
                 "TRANSECT",     
                 "TREATMENT",     
                 "QUAD",     
                 "SIDE",     
                 "SP_CODE",     
                 "NPP_season_gC_m2_day",     
                 "SCIENTIFIC_NAME",     
                 "COMMON_NAME",     
                 "TAXON_KINGDOM",     
                 "TAXON_PHYLUM",     
                 "TAXON_CLASS",     
                 "TAXON_ORDER",     
                 "TAXON_FAMILY"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$SEASON)!="factor") dt1$SEASON<- as.factor(dt1$SEASON)
if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$TREATMENT)!="factor") dt1$TREATMENT<- as.factor(dt1$TREATMENT)
if (class(dt1$QUAD)!="factor") dt1$QUAD<- as.factor(dt1$QUAD)
if (class(dt1$SP_CODE)!="factor") dt1$SP_CODE<- as.factor(dt1$SP_CODE)
if (class(dt1$SIDE)!="factor") dt1$SIDE<- as.factor(dt1$SIDE)
if (class(dt1$NPP_season_gC_m2_day)=="factor") dt1$NPP_season_gC_m2_day <-as.numeric(levels(dt1$NPP_season_gC_m2_day))[as.integer(dt1$NPP_season_gC_m2_day) ]               
if (class(dt1$NPP_season_gC_m2_day)=="character") dt1$NPP_season_gC_m2_day <-as.numeric(dt1$NPP_season_gC_m2_day)
if (class(dt1$SCIENTIFIC_NAME)!="factor") dt1$SCIENTIFIC_NAME<- as.factor(dt1$SCIENTIFIC_NAME)
if (class(dt1$COMMON_NAME)!="factor") dt1$COMMON_NAME<- as.factor(dt1$COMMON_NAME)
if (class(dt1$TAXON_KINGDOM)!="factor") dt1$TAXON_KINGDOM<- as.factor(dt1$TAXON_KINGDOM)
if (class(dt1$TAXON_PHYLUM)!="factor") dt1$TAXON_PHYLUM<- as.factor(dt1$TAXON_PHYLUM)
if (class(dt1$TAXON_CLASS)!="factor") dt1$TAXON_CLASS<- as.factor(dt1$TAXON_CLASS)
if (class(dt1$TAXON_ORDER)!="factor") dt1$TAXON_ORDER<- as.factor(dt1$TAXON_ORDER)
if (class(dt1$TAXON_FAMILY)!="factor") dt1$TAXON_FAMILY<- as.factor(dt1$TAXON_FAMILY)

# Convert Missing Values to NA for non-dates

dt1$NPP_season_gC_m2_day <- ifelse((trimws(as.character(dt1$NPP_season_gC_m2_day))==trimws("-99999")),NA,dt1$NPP_season_gC_m2_day)               
suppressWarnings(dt1$NPP_season_gC_m2_day <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$NPP_season_gC_m2_day))==as.character(as.numeric("-99999"))),NA,dt1$NPP_season_gC_m2_day))
dt1$COMMON_NAME <- as.factor(ifelse((trimws(as.character(dt1$COMMON_NAME))==trimws("-99999")),NA,as.character(dt1$COMMON_NAME)))
dt1$TAXON_KINGDOM <- as.factor(ifelse((trimws(as.character(dt1$TAXON_KINGDOM))==trimws("-99999")),NA,as.character(dt1$TAXON_KINGDOM)))
dt1$TAXON_PHYLUM <- as.factor(ifelse((trimws(as.character(dt1$TAXON_PHYLUM))==trimws("-99999")),NA,as.character(dt1$TAXON_PHYLUM)))
dt1$TAXON_CLASS <- as.factor(ifelse((trimws(as.character(dt1$TAXON_CLASS))==trimws("-99999")),NA,as.character(dt1$TAXON_CLASS)))
dt1$TAXON_ORDER <- as.factor(ifelse((trimws(as.character(dt1$TAXON_ORDER))==trimws("-99999")),NA,as.character(dt1$TAXON_ORDER)))
dt1$TAXON_FAMILY <- as.factor(ifelse((trimws(as.character(dt1$TAXON_FAMILY))==trimws("-99999")),NA,as.character(dt1$TAXON_FAMILY)))

npp.dat.section <- dt1
rm(dt1, infile1, inUrl1)

# ------------------------------------------------------------------------------------
# Replace underscores with dots for convenience. Also convert to lowercase
colnames(biomass.dat.transect) <- tolower(gsub("_", ".", colnames(biomass.dat.transect)))
colnames(depth.dat.transect)   <- tolower(gsub("_", ".", colnames(depth.dat.transect)))
colnames(frond.dat.section)   <- tolower(gsub("_", ".", colnames(frond.dat.section)))
colnames(frond.dat.transect)   <- tolower(gsub("_", ".", colnames(frond.dat.transect)))
colnames(light.dat.transect)   <- tolower(gsub("_", ".", colnames(light.dat.transect)))
colnames(npp.dat.section)      <- tolower(gsub("_", ".", colnames(npp.dat.section)))
colnames(sand.dat.section)     <- tolower(gsub("_", ".", colnames(sand.dat.section)))
colnames(sand.dat.transect)    <- tolower(gsub("_", ".", colnames(sand.dat.transect)))
colnames(urchin.dat.transect)  <- tolower(gsub("_", ".", colnames(urchin.dat.transect)))
colnames(daily.light.kelp.dat)  <- tolower(gsub("_", ".", colnames(daily.light.kelp.dat)))

# ------------------------------------------------------------------------------------
# Replace "-99999" values with NA
biomass.dat.transect[biomass.dat.transect== -99999] <- NA
depth.dat.transect[depth.dat.transect== -99999] <- NA
frond.dat.section[frond.dat.section== -99999] <- NA
frond.dat.transect[frond.dat.transect== -99999] <- NA
light.dat.transect[light.dat.transect== -99999] <- NA
npp.dat.section[npp.dat.section== -99999] <- NA
sand.dat.section[sand.dat.section== -99999] <- NA
sand.dat.transect[sand.dat.transect== -99999] <- NA
urchin.dat.transect[urchin.dat.transect== -99999] <- NA
daily.light.kelp.dat[daily.light.kelp.dat== -99999] <- NA

# ------------------------------------------------------------------------------------
# Add column for season
add.season.column <- function(x){
  x <- x %>%
    
    # Create column for day of year (DOY)
    dplyr::mutate(doy = yday(date)) %>% 
    
    # Create column for quarter
    dplyr::mutate(season = case_when(
      doy>=355 | doy<79  ~ "1-Winter",
      doy>=79  & doy<172 ~ "2-Spring",
      doy>=172 & doy<265 ~ "3-Summer",
      doy>=265 & doy<355 ~ "4-Autumn")) %>%
    
    mutate(time = year + ifelse(season == "1-Winter", 0.125,
                                ifelse(season == "2-Spring", 0.375,
                                       ifelse(season == "3-Summer", 0.625, 0.875)))) %>%
    dplyr::select(-doy) %>%
    dplyr::select(year, season, time, everything())
  
  return(x)
}

frond.dat.transect$date <- frond.dat.transect$date.kelp
frond.dat.transect <- frond.dat.transect %>%
  dplyr::select(-date.kelp)

frond.dat.section$year <- year(frond.dat.section$date)
frond.dat.transect$year <- year(frond.dat.transect$date)

biomass.dat.transect <- add.season.column(biomass.dat.transect)
daily.light.kelp.dat <- add.season.column(daily.light.kelp.dat)
frond.dat.section <- add.season.column(frond.dat.section)
frond.dat.transect <- add.season.column(frond.dat.transect)
sand.dat.section <- add.season.column(sand.dat.section)
sand.dat.transect <- add.season.column(sand.dat.transect)
urchin.dat.transect <- add.season.column(urchin.dat.transect)

# ------------------------------------------------------------------------------------
# Create ID data from biomass data to identify treatments for each transect and time point
id.dat.transect <- biomass.dat.transect %>%
  dplyr::select(site, transect, treatment, time, year, season) %>%
  dplyr::distinct()

# ------------------------------------------------------------------------------------
# Thin data by season

# For 2008-2012, data were gathered twice per season, potentially confounding comparisons to data collected after 2012. 
# For winter measurements, use only the earlier measurements that were taken before removing giant kelp in annual and continual removal transects. 

# Note that we don't have to worry about this for the NPP data, because it has already been averaged by season based on interpolations between dates of measurement.

# ---
# Biomass data
biomass.dat.transect.raw <- biomass.dat.transect

biomass.dat.transect <- biomass.dat.transect %>%
  dplyr::group_by(site, treatment, transect, year, season, time, sp.code, scientific.name) %>%
  dplyr::slice_min(date, n = 1) %>%
  dplyr::ungroup()

# ---
# Frond data
frond.dat.section <- frond.dat.section %>%
  dplyr::group_by(site, treatment, transect, quad, side, year, season, time) %>%
  dplyr::slice_min(date, n = 1) %>%
  dplyr::ungroup()

frond.dat.transect <- frond.dat.transect %>%
  dplyr::group_by(site, treatment, transect, year, season, time) %>%
  dplyr::slice_min(date, n = 1) %>%
  dplyr::ungroup()

# ---
# Sand data
sand.dat.section <- sand.dat.section %>%
  dplyr::group_by(site, treatment, transect, quad, side, year, season, time) %>%
  dplyr::slice_min(date, n = 1) %>%
  dplyr::ungroup()

sand.dat.transect <- sand.dat.transect %>%
  dplyr::group_by(site, treatment, transect, year, season, time) %>%
  dplyr::slice_min(date, n = 1) %>%
  dplyr::ungroup()

# ---
# Urchin data
urchin.dat.transect <- urchin.dat.transect %>%
  dplyr::group_by(site, treatment, transect, year, season, time) %>%
  dplyr::slice_min(date, n = 1) %>%
  dplyr::ungroup()

# ------------------------------------------------------------------------------------
# NPP DATA: Recode season data
npp.dat.section$season <- as.character(npp.dat.section$season)

npp.dat.section$season[npp.dat.section$season == "1-WINTER"] <- "1-Winter"
npp.dat.section$season[npp.dat.section$season == "2-SPRING"] <- "2-Spring"
npp.dat.section$season[npp.dat.section$season == "3-SUMMER"] <- "3-Summer"
npp.dat.section$season[npp.dat.section$season == "4-AUTUMN"] <- "4-Autumn"

# ------------------------------------------------------------------------------------
# Reformat columns in destination data frames
light.dat.transect$year <- as.numeric(as.character(light.dat.transect$year))
light.dat.transect$transect <- as.factor(as.character(light.dat.transect$transect))

npp.dat.section$year <- as.numeric(as.character(npp.dat.section$year))
npp.dat.section$transect <- as.factor(as.character(npp.dat.section$transect))

# Join light data
light.dat.transect2 <- dplyr::left_join(light.dat.transect, id.dat.transect) %>%
  mutate(time = year + ifelse(season == "1-Winter", 0.125,
                              ifelse(season == "2-Spring", 0.375,
                                     ifelse(season == "3-Summer", 0.625, 0.875))))

# Remove light measurements taken before some Continual transects were really started (winter prior to official start)
light.dat.transect2 <- light.dat.transect2[!is.na(light.dat.transect2$treatment), ]
light.dat.transect <- light.dat.transect2

# Join NPP data
npp.dat.section2 <- dplyr::left_join(npp.dat.section, id.dat.transect) %>%
  mutate(time = year + ifelse(season == "1-Winter", 0.125,
                              ifelse(season == "2-Spring", 0.375,
                                     ifelse(season == "3-Summer", 0.625, 0.875))))
npp.dat.section <- npp.dat.section2

# Join depth data (separate process because there is no time component)
id.dat.transect2 <- id.dat.transect %>%
  dplyr::select(site, transect, treatment) %>%
  dplyr::distinct()

depth.dat.transect$transect <- as.factor(as.character(depth.dat.transect$transect))

depth.dat.transect2 <- dplyr::left_join(depth.dat.transect, id.dat.transect2)
depth.dat.transect <- depth.dat.transect2

# Rearrange columns
light.dat.transect <- light.dat.transect %>%
  dplyr::select(year, season, time, site, transect, treatment, everything())

npp.dat.section <- npp.dat.section %>%
  dplyr::select(year, season, time, site, transect, treatment, quad, side, everything())

depth.dat.transect <- depth.dat.transect %>%
  dplyr::select(site, transect, treatment, depth.mllw.m)

# Clean up
rm(light.dat.transect2, depth.dat.transect2, npp.dat.section2, id.dat.transect2)

# ------------------------------------------------------------------------------------
# TREATMENT NAMES: Convert all treatment names to sentence case 
biomass.dat.transect.raw$treatment <- str_to_sentence(biomass.dat.transect.raw$treatment, locale = "en") 
biomass.dat.transect$treatment  <- str_to_sentence(biomass.dat.transect$treatment, locale = "en")
depth.dat.transect$treatment    <- str_to_sentence(depth.dat.transect$treatment, locale = "en")
frond.dat.section$treatment    <- str_to_sentence(frond.dat.section$treatment, locale = "en")
frond.dat.transect$treatment    <- str_to_sentence(frond.dat.transect$treatment, locale = "en")
id.dat.transect$treatment       <- str_to_sentence(id.dat.transect$treatment, locale = "en")
light.dat.transect$treatment    <- str_to_sentence(light.dat.transect$treatment, locale = "en")
npp.dat.section$treatment       <- str_to_sentence(npp.dat.section$treatment, locale = "en")
sand.dat.transect$treatment     <- str_to_sentence(sand.dat.transect$treatment, locale = "en")
sand.dat.section$treatment      <- str_to_sentence(sand.dat.section$treatment, locale = "en")
urchin.dat.transect$treatment   <- str_to_sentence(urchin.dat.transect$treatment, locale = "en")

# ------------------------------------------------------------------------------------
# Re-code sites
site.recode <- function(x){
  x$site = as.character(x$site)
  
  x$site[x$site == "AQUE"] <- "Arroyo Quemado"
  x$site[x$site == "CARP"] <- "Carpinteria"
  x$site[x$site == "IVEE"] <- "Isla Vista"
  x$site[x$site == "MOHK"] <- "Mohawk"
  x$site[x$site == "NAPL"] <- "Naples"
  
  return(x)
}

biomass.dat.transect.raw <- site.recode(biomass.dat.transect.raw)
biomass.dat.transect <- site.recode(biomass.dat.transect)
depth.dat.transect   <- site.recode(depth.dat.transect)
frond.dat.section   <- site.recode(frond.dat.section)
frond.dat.transect  <- site.recode(frond.dat.transect)
id.dat.transect      <- site.recode(id.dat.transect)
light.dat.transect   <- site.recode(light.dat.transect)
npp.dat.section      <- site.recode(npp.dat.section)
sand.dat.section     <- site.recode(sand.dat.section)
sand.dat.transect    <- site.recode(sand.dat.transect)
urchin.dat.transect  <- site.recode(urchin.dat.transect)

# ------------------------------------------------------------------------------------
# LIGHT DATA: Reorganize and rename columns
light.dat.transect <- light.dat.transect %>%
  dplyr::rename(mean.bottom.par = avg.seafloor.mol, 
                n.bottom.par = seafloor.n) %>% 
                #mean.surface.par = avg.surface.mol, 
                #n.surface.par = surface.n
  dplyr::select(site, transect, treatment, everything())

# ------------------------------------------------------------------------------------
# Recode month as a numeric, where appropriate
biomass.dat.transect.raw$month <- as.numeric(as.character(biomass.dat.transect.raw$month))
biomass.dat.transect$month <- as.numeric(as.character(biomass.dat.transect$month))


# ------------------------------------------------------------------------------------
# RAW BIOMASS TRANSECT DATA FOR GIANT KELP ONLY

# Sum juveniles and adults for each species
biomass.dat.transect.raw <- biomass.dat.transect.raw %>%
  dplyr::select(year, season, time, month, date, site, transect, treatment, scientific.name, dry.gm2) %>%
  dplyr::group_by(year, season, time, month, date, site, transect, treatment, scientific.name) %>%
  dplyr::summarise(dry.gm2 = sum(dry.gm2, na.rm = TRUE)) %>%
  ungroup() 

giant.kelp.biomass.transect.dat.raw <- biomass.dat.transect.raw %>%
  dplyr::filter(scientific.name == "Macrocystis pyrifera")
rm(biomass.dat.transect.raw)

# Add ID column
giant.kelp.biomass.transect.dat.raw <- giant.kelp.biomass.transect.dat.raw %>%
  dplyr::mutate(id = paste(year, season, site, treatment, sep = "_")) %>%
  dplyr::select(id, everything())

# ------------------------------------------------------------------------------------
# BIOMASS TRANSECT DATA: Clean up to get total biomass for groups

# Sum juveniles and adults for each species
biomass.dat.transect <- biomass.dat.transect %>%
  dplyr::select(year, season, time, month, date, site, transect, treatment, scientific.name, dry.gm2) %>%
  dplyr::group_by(year, season, time, month, date, site, transect, treatment, scientific.name) %>%
  dplyr::summarise(dry.gm2 = sum(dry.gm2, na.rm = TRUE)) %>%
  ungroup()

# Spread biomass data
biomass.dat.transect.wide <- biomass.dat.transect %>%
  tidyr::spread(key = scientific.name, value = dry.gm2, fill = 0)

# Check whether any species should be combined to genus level
colnames(biomass.dat.transect.wide)[-c(1:8)]
colnames(biomass.dat.transect.wide)[c(1:8)]

biomass.dat.transect <- biomass.dat.transect.wide %>%
  tidyr::gather(key = "taxon", value = "dry.gm2", -c(year:treatment))

biomass.dat.transect.wide <- biomass.dat.transect %>%
  tidyr::spread(key = taxon, value = dry.gm2, fill = 0)

# Add ID column
biomass.dat.transect <- biomass.dat.transect %>%
  dplyr::mutate(id = paste(year, season, site, treatment, sep = "_")) %>%
  dplyr::select(id, everything())

biomass.dat.transect.wide <- biomass.dat.transect.wide %>%
  dplyr::mutate(id = paste(year, season, site, treatment, sep = "_"))%>%
  dplyr::select(id, everything())

# Create understory biomass data frame and matrix
understory.transect.dat.wide <- biomass.dat.transect.wide %>%
  dplyr::select(-`Macrocystis pyrifera`)

understory.transect.mat <- understory.transect.dat.wide %>%
  dplyr::select(-c(id:treatment)) %>%
  as.matrix()

rownames(understory.transect.mat) <- understory.transect.dat.wide$id

# ------------------------------------------------------------------------------------
# BIOMASS TRANSECT DATA: Clean up column names, then aggregate understory biomass by species
biomass.sum.dat.transect <- biomass.dat.transect.wide %>%
  dplyr::select(id:treatment, `Macrocystis pyrifera`) %>%
  dplyr::rename(giant.kelp.dry.gm2 = `Macrocystis pyrifera`)

biomass.sum.dat.transect$understory.dry.gm2 <- rowSums(understory.transect.mat)
                            
# ------------------------------------------------------------------------------------
# Fix frond.dat.transect
frond.dat.transect <- frond.dat.transect %>%
  dplyr::select(year, season, time, date, everything())

# ------------------------------------------------------------------------------------
# JOIN TRANSECT-SCALE DATA: Biomass, light, depth, sand depth, and urchins
biomass.dat.transect$transect     <- as.factor(as.character(biomass.dat.transect$transect))
biomass.sum.dat.transect$transect <- as.factor(as.character(biomass.sum.dat.transect$transect))
depth.dat.transect$transect       <- as.factor(as.character(depth.dat.transect$transect))
frond.dat.transect$transect       <- as.factor(as.character(frond.dat.transect$transect))
light.dat.transect$transect       <- as.factor(as.character(light.dat.transect$transect))
npp.dat.section$transect          <- as.factor(as.character(npp.dat.section$transect))
sand.dat.section$transect         <- as.factor(as.character(sand.dat.section$transect))
sand.dat.transect$transect        <- as.factor(as.character(sand.dat.transect$transect))
urchin.dat.transect$transect      <- as.factor(as.character(urchin.dat.transect$transect))

urchin.dat.transect$month <- as.numeric(as.character(urchin.dat.transect$month))
sand.dat.transect$month <- as.numeric(as.character(sand.dat.transect$month))

depth.dat.transect2 <- depth.dat.transect %>%
  dplyr::select(site, transect, depth.mllw.m) %>%
  distinct()

transect.dat <- biomass.sum.dat.transect %>%
  #dplyr::select(-month, -date) %>%
  left_join(light.dat.transect) %>%
  left_join(depth.dat.transect2) %>%
  left_join(urchin.dat.transect) %>%
  left_join(sand.dat.transect) %>%
  left_join(frond.dat.transect) %>%
  dplyr::mutate(sqrt.giant.kelp.dry.gm2 = sqrt(giant.kelp.dry.gm2),
                sqrt.understory.dry.gm2 = sqrt(understory.dry.gm2),
                sqrt.mean.bottom.par = sqrt(mean.bottom.par))

transect.dat$season <- factor(transect.dat$season, levels = c("1-Winter", "2-Spring", "3-Summer", "4-Autumn"))
transect.dat$plot <- factor(paste(transect.dat$site, transect.dat$treatment, sep = "_"))

transect.dat <- transect.dat %>%
  dplyr::select(id:site, treatment,  plot, transect,
                depth.mllw.m, mean.bottom.par, sqrt.mean.bottom.par, n.bottom.par, 
                fronds.per.transect, plants.per.transect, fronds.per.m2, plants.per.m2, fronds.per.plant,
                giant.kelp.dry.gm2, sqrt.giant.kelp.dry.gm2, understory.dry.gm2, sqrt.understory.dry.gm2,
                percent.sand, urchin.dens.no.m2, urchin.dry.gm2, urchin.sfdm.gm2)

# ------------------------------------------------------------------------------------
# NPP DATA: Aggregate NPP by species

# Are there NAs in the NPP data?
table(is.na(npp.dat.section$npp.season.gc.m2.day))

# Sum NPP of juveniles and adults for each species
npp.dat.section <- npp.dat.section %>%
  dplyr::select(year, season, time, site, transect, treatment, quad, side, scientific.name, npp.season.gc.m2.day) %>%
  dplyr::group_by(year, season, time, site, transect, treatment, quad, side, scientific.name) %>%
  dplyr::summarise(npp.season.gc.m2.day = sum(npp.season.gc.m2.day, na.rm = TRUE)) %>%
  ungroup()

# Spread NPP data
npp.dat.section.wide <- npp.dat.section %>%
  tidyr::spread(key = scientific.name, value = npp.season.gc.m2.day, fill = 0)

# Check whether any species should be combined to genus level
colnames(npp.dat.section.wide)[-c(1:8)]
colnames(npp.dat.section.wide)[c(1:8)]

# Gather NPP data
npp.dat.section <- npp.dat.section.wide %>%
  tidyr::gather(key = "taxon", value = "npp.season.gc.m2.day", -c(year:side))

# Add ID column
npp.dat.section <- npp.dat.section %>%
  dplyr::mutate(id = paste(year, season, site, treatment, quad, side, sep = "_"),
                section = paste(quad, side, sep = "_")) %>%
  dplyr::select(id, year, season, time, site, treatment, transect, section, quad, side, taxon, npp.season.gc.m2.day)

npp.dat.section.wide <- npp.dat.section.wide %>%
  dplyr::mutate(id = paste(year, season, site, treatment, quad, side, sep = "_"),
                section = paste(quad, side, sep = "_")) %>%
  dplyr::select(id, year, season, time, site, treatment, transect, section, quad, side, everything()) 

# Check that there are equal numbers of sections
t(t(table(npp.dat.section$section)))
t(t(table(npp.dat.section.wide$section)))

# ------------------------------------------------------------------------------------
# For both NPP and biomass data at the section level, recode "section" as a unique section ID for each individual section
npp.dat.section$section          <- paste(npp.dat.section$site, npp.dat.section$treatment, npp.dat.section$section, sep = "_")
npp.dat.section.wide$section     <- paste(npp.dat.section.wide$site, npp.dat.section.wide$treatment, npp.dat.section.wide$section, sep = "_")

# ------------------------------------------------------------------------------------
# Calculate total NPP per section per quarter
npp.mat <- npp.dat.section.wide %>%
  dplyr::select(-c(id:side)) %>%
  as.matrix()
rownames(npp.mat) <- npp.dat.section.wide$id

npp.dat.section.sum <- npp.dat.section.wide
npp.dat.section.sum$total.npp.season.gc.m2.day <- rowSums(npp.mat)

# Clean up
npp.dat.section.sum <- npp.dat.section.sum %>%
  dplyr::select(id:side, total.npp.season.gc.m2.day)

# ------------------------------------------------------------------------------------
# Calculate total NPP per section per quarter - GIANT KELP ONLY
n <- which(colnames(npp.mat) %in% "Macrocystis pyrifera")
npp.dat.section.sum$giant.kelp.npp.season.gc.m2.day <- npp.mat[, n] 

rm(n)

# ------------------------------------------------------------------------------------
# Rename datasets
npp.dat.section.long <- npp.dat.section
npp.dat.section <- npp.dat.section.sum

biomass.dat.transect.long <- biomass.dat.transect
biomass.dat.transect <- biomass.sum.dat.transect

rm(npp.mat, npp.dat.section.sum, biomass.sum.dat.transect)

npp.dat.section$understory.npp.season.gc.m2.day = npp.dat.section$total.npp.season.gc.m2.day - 
                                                  npp.dat.section$giant.kelp.npp.season.gc.m2.day

# Make plotting version of NPP data 
npp.dat.section.plotting <- npp.dat.section %>%
  dplyr::rename(`Giant kelp` = giant.kelp.npp.season.gc.m2.day,
                `Understory` = understory.npp.season.gc.m2.day,
                `Total` = total.npp.season.gc.m2.day) %>%
  tidyr::gather(key = "guild", value = "guild.npp.season.gc.m2.day", -c(id:side))

# Reassign characters as factors
npp.dat.section.plotting$season <- factor(npp.dat.section.plotting$season, levels = c("1-Winter", "2-Spring", "3-Summer", "4-Autumn"))
npp.dat.section.plotting$site <- factor(npp.dat.section.plotting$site,
                                levels = c("Carpinteria", "Naples", "Arroyo Quemado", "Isla Vista", "Mohawk"))

npp.dat.section.plotting$treatment <- factor(npp.dat.section.plotting$treatment, 
                        levels = c("Control", "Annual", "Continual")) 
npp.dat.section.plotting$plot <- factor(paste(npp.dat.section.plotting$site, npp.dat.section.plotting$treatment, sep = "_"))
npp.dat.section.plotting$guild <- factor(npp.dat.section.plotting$guild, levels = c("Total", "Giant kelp", "Understory"))

# ------------------------------------------------------------------------------------
# Join section-scale NPP data with other section-scale data (and transect-scale data where there are no section-scale data)
# Note that the NPP data is longer because it is at the section level, not the plot level

light.dat.transect$year      <- as.numeric(light.dat.transect$year)
npp.dat.section$year         <- as.numeric(npp.dat.section$year)
sand.dat.section$year        <- as.numeric(sand.dat.section$year)
sand.dat.transect$year       <- as.numeric(sand.dat.transect$year)
frond.dat.section$year       <- as.numeric(frond.dat.section$year)
frond.dat.transect$year      <- as.numeric(frond.dat.transect$year)
transect.dat$year            <- as.numeric(transect.dat$year)
urchin.dat.transect$year     <- as.numeric(urchin.dat.transect$year)

depth.dat.transect$transect  <- as.character(depth.dat.transect$transect)
light.dat.transect$transect  <- as.character(light.dat.transect$transect)
npp.dat.section$transect     <- as.character(npp.dat.section$transect)
sand.dat.section$transect    <- as.character(sand.dat.section$transect)
sand.dat.transect$transect   <- as.character(sand.dat.transect$transect)
frond.dat.section$transect   <- as.character(frond.dat.section$transect)
frond.dat.transect$transect  <- as.character(frond.dat.transect$transect)
transect.dat$transect        <- as.character(transect.dat$transect)
urchin.dat.transect$transect <- as.character(urchin.dat.transect$transect)

npp.dat.section$quad     <- as.numeric(as.character(npp.dat.section$quad))
sand.dat.section$quad    <- as.numeric(as.character(sand.dat.section$quad))
frond.dat.section$quad   <- as.numeric(as.character(frond.dat.section$quad))

npp.dat.section$side     <- as.character(npp.dat.section$side)
sand.dat.section$side    <- as.character(sand.dat.section$side)
frond.dat.section$side   <- as.character(frond.dat.section$side)

# Create new "section" column for sand and frond data
sand.dat.section <- sand.dat.section %>%
  dplyr::mutate(id = paste(year, season, site, treatment, quad, side, sep = "_"),
                section = paste(site, treatment, quad, side, sep = "_")) %>%
  dplyr::select(id, year, season, time, site, treatment, transect, section, quad, side, everything()) 

frond.dat.section <- frond.dat.section %>%
  dplyr::mutate(id = paste(year, season, site, treatment, quad, side, sep = "_"),
                section = paste(site, treatment, quad, side, sep = "_")) %>%
  dplyr::select(id, year, season, time, site, treatment, transect, section, quad, side, everything()) 

# Join section-scale data
section.dat <- npp.dat.section %>%
  dplyr::left_join(light.dat.transect, by = c("year", "season", "time", "site", "treatment", "transect")) %>%
  dplyr::left_join(depth.dat.transect, by = c("site", "treatment", "transect")) %>%
  dplyr::left_join(sand.dat.section, by = c("id", "year", "season", "time", "site", "treatment", "transect", "section", "quad", "side")) %>%
  dplyr::left_join(frond.dat.section, by = c("id", "year", "season", "time", "site", "treatment", "transect", "section", "quad", "side")) %>%
  dplyr::select(id:side, depth.mllw.m, mean.bottom.par, n.bottom.par, 
                fronds.per.section, plants.per.section, fronds.per.m2, plants.per.m2, fronds.per.plant,
                total.npp.season.gc.m2.day, giant.kelp.npp.season.gc.m2.day, understory.npp.season.gc.m2.day,
                percent.sand) %>%
  dplyr::rename(depth.mllw.m_at.transect = depth.mllw.m,
                mean.bottom.par_at.transect = mean.bottom.par,
                n.bottom.par_at.transect = n.bottom.par,
                fronds.per.m2_at.section = fronds.per.m2, 
                plants.per.m2_at.section = plants.per.m2, 
                fronds.per.plant_at.section = fronds.per.plant,
                total.npp.season.gc.m2.day_at.section = total.npp.season.gc.m2.day,
                giant.kelp.npp.season.gc.m2.day_at.section = giant.kelp.npp.season.gc.m2.day,
                understory.npp.season.gc.m2.day_at.section = understory.npp.season.gc.m2.day,
                percent.sand_at.section = percent.sand)

# Reassign characters as factors
section.dat$season <- factor(section.dat$season, levels = c("1-Winter", "2-Spring", "3-Summer", "4-Autumn"))
section.dat$site <- factor(section.dat$site)
section.dat$treatment <- factor(section.dat$treatment, 
                        levels = c("Control", "Annual", "Continual")) 
section.dat$plot <- factor(paste(section.dat$site, section.dat$treatment, sep = "_"))
section.dat$section <- factor(section.dat$section)

section.dat <- section.dat %>%
  dplyr::select(id:treatment, plot, everything())

# ------------------------------------------------------------------------------------
# For transect scale data, average NPP data across quads and merge with other data

npp.dat.transect <- npp.dat.section %>%
  dplyr::group_by(year, season, time, site, treatment, transect) %>%
  dplyr::summarise(giant.kelp.npp.season.gc.m2.day_at.transect = mean(giant.kelp.npp.season.gc.m2.day, na.rm = T),
                  understory.npp.season.gc.m2.day_at.transect  = mean(understory.npp.season.gc.m2.day, na.rm = T),
                  total.npp.season.gc.m2.day_at.transect       = mean(total.npp.season.gc.m2.day, na.rm = T)) %>%
  dplyr::mutate(id = paste(year, season, site, treatment, sep = "_")) %>%
  dplyr::select(id, everything()) %>%
  dplyr::ungroup()

# Join data
transect.dat$season <- as.character(transect.dat$season)

transect.dat2 <- transect.dat %>%
  dplyr::left_join(npp.dat.transect, by = c("id", "year", "season", "time", "site", "treatment", "transect")) %>%
  dplyr::select(id:depth.mllw.m, 
                giant.kelp.npp.season.gc.m2.day_at.transect:total.npp.season.gc.m2.day_at.transect, 
                everything()) %>%
  dplyr::rename(giant.kelp.npp.season.gc.m2.day = giant.kelp.npp.season.gc.m2.day_at.transect,
                understory.npp.season.gc.m2.day = understory.npp.season.gc.m2.day_at.transect,
                total.npp.season.gc.m2.day      = total.npp.season.gc.m2.day_at.transect)

transect.dat <- transect.dat2
rm(transect.dat2)


# ====================================================================================
# CALCULATE DELTA-NPP DATA
# ====================================================================================

# ------------------------------------------------------------------------------------
# Using the transect-scale averages (because this would not be possible with section-level data), calculate "delta NPP" as the difference of NPP in the control plot and the annual and continual removal plots (X - Control)
delta.npp.dat <- npp.dat.transect %>%
  dplyr::select(-transect, -id) %>%
  tidyr::pivot_longer(-(year:treatment), 
                      names_to = "npp.component", 
                      values_to = "npp.season.gc.m2.day") %>%
  dplyr::mutate(npp.component = dplyr::recode(npp.component,
                                              `understory.npp.season.gc.m2.day_at.transect` = "understory",
                                              `giant.kelp.npp.season.gc.m2.day_at.transect` = "giant.kelp",
                                              `total.npp.season.gc.m2.day_at.transect` = "total")) %>%
  tidyr::pivot_wider(id_cols = c(year:site, npp.component),
                     names_from = treatment,
                     values_from = npp.season.gc.m2.day) %>%
  dplyr::mutate(delta.annual = Annual - Control,
                delta.continual = Continual - Control) %>%
  dplyr::select(-(Annual:Continual)) %>%
  tidyr::pivot_longer(-(year:npp.component),
                      names_to = "treatment",
                      values_to = "delta.npp") %>%
  dplyr::mutate(treatment = dplyr::recode(treatment,
                                          delta.annual = "Annual",
                                          delta.continual = "Continual"))

delta.npp.dat$treatment <- factor(delta.npp.dat$treatment, levels = c("Annual", "Continual"))

delta.npp.dat$site <- factor(delta.npp.dat$site,
                             levels = c("Carpinteria", "Naples", "Arroyo Quemado", "Isla Vista", "Mohawk"))

delta.npp.dat <- delta.npp.dat[!is.na(delta.npp.dat$delta.npp), ]
  

# ====================================================================================
# FILTER DATA BEFORE MOVING ON TO ADDITIONAL PLOTS AND ANALYSES
# ====================================================================================

ggplot(data = transect.dat[transect.dat$treatment == "Annual", ], 
       aes(x = time, y = giant.kelp.dry.gm2)) +
  geom_vline(xintercept = c(2012.375, 2017.125), color = 'red') +
  geom_line() +
  facet_wrap(~ site) +
  scale_x_continuous(breaks = 2007:2020)


# ------------------------------------------------------------------------------------
# For annual removal treatments, remove data > 1 year after last removal (in other words, if the last removal was in Winter 2016, we include data through Winter 2017 because winter measurements are made before cutting back kelp)

# For continual removal treatments, remove data > 1 quarter/season after last removal (in other words, if the last removal was in Winter 2017, we include data through Spring 2017 because measurements are made before cutting back kelp)

# For control treatments, match to the longest length of data at that site after making above changes

# Note, this approach retains the first 

treatment.filter <- function(x){
  x$flag.to.remove <- 0
  
  # Start dates
  x$flag.to.remove[x$site == "Mohawk"         & x$treatment == "Annual" & x$time < 2008.375] <- 1
  x$flag.to.remove[x$site == "Isla Vista"     & x$treatment == "Annual" & x$time < 2012.375] <- 1
  x$flag.to.remove[x$site == "Arroyo Quemado" & x$treatment == "Annual" & x$time < 2008.375] <- 1
  x$flag.to.remove[x$site == "Naples"         & x$treatment == "Annual" & x$time < 2008.375] <- 1
  x$flag.to.remove[x$site == "Carpinteria"    & x$treatment == "Annual" & x$time < 2008.375] <- 1
  
  x$flag.to.remove[x$site == "Mohawk"         & x$treatment == "Continual" & x$time < 2010.375] <- 1 #< 2010.625] <- 1
  # No continual removal at Isla Vista
  x$flag.to.remove[x$site == "Arroyo Quemado" & x$treatment == "Continual" & x$time < 2010.375] <- 1 #< 2010.625] <- 1
  x$flag.to.remove[x$site == "Naples"         & x$treatment == "Continual" & x$time < 2010.375] <- 1 #< 2010.625] <- 1
  x$flag.to.remove[x$site == "Carpinteria"    & x$treatment == "Continual" & x$time < 2010.375] <- 1 #< 2010.625] <- 1
  
  x$flag.to.remove[x$site == "Mohawk"         & x$treatment == "Control" & x$time < 2008.375] <- 1
  x$flag.to.remove[x$site == "Isla Vista"     & x$treatment == "Control" & x$time < 2012.375] <- 1
  x$flag.to.remove[x$site == "Arroyo Quemado" & x$treatment == "Control" & x$time < 2008.375] <- 1
  x$flag.to.remove[x$site == "Naples"         & x$treatment == "Control" & x$time < 2008.375] <- 1
  x$flag.to.remove[x$site == "Carpinteria"    & x$treatment == "Control" & x$time < 2008.375] <- 1
  
  # End dates
  x$flag.to.remove[x$site == "Mohawk"         & x$treatment == "Annual" & x$time > 2018.125] <- 1
  x$flag.to.remove[x$site == "Isla Vista"     & x$treatment == "Annual" & x$time > 2017.125] <- 1
  x$flag.to.remove[x$site == "Arroyo Quemado" & x$treatment == "Annual" & x$time > 2018.125] <- 1
  x$flag.to.remove[x$site == "Naples"         & x$treatment == "Annual" & x$time > 2017.125] <- 1
  x$flag.to.remove[x$site == "Carpinteria"    & x$treatment == "Annual" & x$time > 2018.125] <- 1
  
  x$flag.to.remove[x$site == "Mohawk"         & x$treatment == "Continual" & x$time > 2017.375] <- 1
  # No continual removal at Isla Vista
  x$flag.to.remove[x$site == "Arroyo Quemado" & x$treatment == "Continual" & x$time > 2017.375] <- 1
  x$flag.to.remove[x$site == "Naples"         & x$treatment == "Continual" & x$time > 2016.375] <- 1
  x$flag.to.remove[x$site == "Carpinteria"    & x$treatment == "Continual" & x$time > 2017.375] <- 1
  
  x$flag.to.remove[x$site == "Mohawk"         & x$treatment == "Control" & x$time > 2018.125] <- 1
  x$flag.to.remove[x$site == "Isla Vista"     & x$treatment == "Control" & x$time > 2017.125] <- 1
  x$flag.to.remove[x$site == "Arroyo Quemado" & x$treatment == "Control" & x$time > 2018.125] <- 1
  x$flag.to.remove[x$site == "Naples"         & x$treatment == "Control" & x$time > 2017.125] <- 1
  x$flag.to.remove[x$site == "Carpinteria"    & x$treatment == "Control" & x$time > 2018.125] <- 1
  
  x <- x %>%
    dplyr::filter(flag.to.remove != 1) %>%
    dplyr::select(-flag.to.remove)
  
  return(x)
}

giant.kelp.biomass.transect.dat.raw <- treatment.filter(giant.kelp.biomass.transect.dat.raw)
biomass.dat.transect          <- treatment.filter(biomass.dat.transect)
biomass.dat.transect.long     <- treatment.filter(biomass.dat.transect.long)
biomass.dat.transect.wide     <- treatment.filter(biomass.dat.transect.wide)
delta.npp.dat                 <- treatment.filter(delta.npp.dat)
frond.dat.section             <- treatment.filter(frond.dat.section)
frond.dat.transect            <- treatment.filter(frond.dat.transect)
light.dat.transect            <- treatment.filter(light.dat.transect)
npp.dat.section               <- treatment.filter(npp.dat.section)
npp.dat.section.long          <- treatment.filter(npp.dat.section.long)
npp.dat.section.plotting      <- treatment.filter(npp.dat.section.plotting)
npp.dat.section.wide          <- treatment.filter(npp.dat.section.wide)
npp.dat.transect              <- treatment.filter(npp.dat.transect)
sand.dat.section              <- treatment.filter(sand.dat.section)
sand.dat.transect             <- treatment.filter(sand.dat.transect)
section.dat                   <- treatment.filter(section.dat)
transect.dat                  <- treatment.filter(transect.dat)
understory.transect.dat.wide  <- treatment.filter(understory.transect.dat.wide)
urchin.dat.transect           <- treatment.filter(urchin.dat.transect)

# ------------------------------------------------------------------------------------
# Fix treatment levels 

giant.kelp.biomass.transect.dat.raw$treatment  <- factor(giant.kelp.biomass.transect.dat.raw$treatment, levels = c("Control", "Annual", "Continual"))
biomass.dat.transect$treatment          <- factor(biomass.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))
biomass.dat.transect.long$treatment     <- factor(biomass.dat.transect.long$treatment, levels = c("Control", "Annual", "Continual"))
biomass.dat.transect.wide$treatment     <- factor(biomass.dat.transect.wide$treatment, levels = c("Control", "Annual", "Continual"))
delta.npp.dat$treatment                 <- factor(delta.npp.dat$treatment, levels = c("Annual", "Continual"))
frond.dat.section$treatment             <- factor(frond.dat.section$treatment, levels = c("Control", "Annual", "Continual"))
frond.dat.transect$treatment            <- factor(frond.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))
id.dat.transect$treatment               <- factor(id.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))
light.dat.transect$treatment            <- factor(light.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))
npp.dat.section$treatment               <- factor(npp.dat.section$treatment, levels = c("Control", "Annual", "Continual"))
npp.dat.section.long$treatment          <- factor(npp.dat.section.long$treatment, levels = c("Control", "Annual", "Continual"))
npp.dat.section.plotting$treatment      <- factor(npp.dat.section.plotting$treatment, levels = c("Control", "Annual", "Continual"))
npp.dat.section.wide$treatment          <- factor(npp.dat.section.wide$treatment, levels = c("Control", "Annual", "Continual"))
npp.dat.transect$treatment              <- factor(npp.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))
sand.dat.section$treatment              <- factor(sand.dat.section$treatment, levels = c("Control", "Annual", "Continual"))
sand.dat.transect$treatment             <- factor(sand.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))
section.dat$treatment                   <- factor(section.dat$treatment, levels = c("Control", "Annual", "Continual"))
transect.dat$treatment                  <- factor(transect.dat$treatment, levels = c("Control", "Annual", "Continual"))
understory.transect.dat.wide$treatment  <- factor(understory.transect.dat.wide$treatment, levels = c("Control", "Annual", "Continual"))
urchin.dat.transect$treatment           <- factor(urchin.dat.transect$treatment, levels = c("Control", "Annual", "Continual"))

# ------------------------------------------------------------------------------------
# Remove particular sites here
# 
# npp.dat.section.all.sites <- npp.dat.section
# section.dat.all.sites <- section.dat
# 
# site.filter <- function(x){
#   x <- x %>%
#     dplyr::filter(site != "Carpinteria") %>%
#     droplevels()
#   return(x)
# }
# 
# biomass.dat.transect          <- site.filter(biomass.dat.transect)
# biomass.dat.transect.long     <- site.filter(biomass.dat.transect.long)
# biomass.dat.transect.wide     <- site.filter(biomass.dat.transect.wide)
# delta.dat                     <- site.filter(delta.dat)
# delta.npp.dat                 <- site.filter(delta.npp.dat)
# delta.npp.dat.annual          <- site.filter(delta.npp.dat.annual)
# frond.dat.section             <- site.filter(frond.dat.section)
# frond.dat.transect            <- site.filter(frond.dat.transect)
# id.dat.transect               <- site.filter(id.dat.transect)
# light.dat.transect            <- site.filter(light.dat.transect)
# npp.dat.section               <- site.filter(npp.dat.section)
# npp.dat.section.long          <- site.filter(npp.dat.section.long)
# npp.dat.section.plotting      <- site.filter(npp.dat.section.plotting)
# npp.dat.section.wide          <- site.filter(npp.dat.section.wide)
# npp.dat.transect              <- site.filter(npp.dat.transect)
# sand.dat.section              <- site.filter(sand.dat.section)
# sand.dat.transect             <- site.filter(sand.dat.transect)
# section.dat                   <- site.filter(section.dat)
# transect.dat                  <- site.filter(transect.dat)
# understory.section.dat.wide   <- site.filter(understory.section.dat.wide)
# understory.transect.dat.wide  <- site.filter(understory.transect.dat.wide)
# urchin.dat.transect           <- site.filter(urchin.dat.transect)

# # ------------------------------------------------------------------------------------
# # Remove data prior to 2010 at Naples because it was an urchin barren
# 
# naples.filter <- function(x){
#   x$flag <- 0
#   x$flag[x$site == "Naples" & x$year < 2009.875] <- 1
#   
#   # Note, we set the date above to 2009.875 instead of 2010.125 because farther down in the code we will
#   # be dropping the first quarter of all time series anyway, so that will cause the data to start in 2010
#   
#   x <- x %>%
#     dplyr::filter(flag == 0) %>%
#     droplevels() %>%
#     dplyr::select(-flag)
#   
#   return(x)
# }
# 
# biomass.dat.transect          <- naples.filter(biomass.dat.transect)
# biomass.dat.transect.long     <- naples.filter(biomass.dat.transect.long)
# biomass.dat.transect.wide     <- naples.filter(biomass.dat.transect.wide)
# delta.dat                     <- naples.filter(delta.dat)
# delta.npp.dat                 <- naples.filter(delta.npp.dat)
# delta.npp.dat.annual          <- naples.filter(delta.npp.dat.annual)
# frond.dat.section             <- naples.filter(frond.dat.section)
# frond.dat.transect            <- naples.filter(frond.dat.transect)
# id.dat.transect               <- naples.filter(id.dat.transect)
# light.dat.transect            <- naples.filter(light.dat.transect)
# npp.dat.section               <- naples.filter(npp.dat.section)
# npp.dat.section.long          <- naples.filter(npp.dat.section.long)
# npp.dat.section.plotting      <- naples.filter(npp.dat.section.plotting)
# npp.dat.section.wide          <- naples.filter(npp.dat.section.wide)
# npp.dat.transect              <- naples.filter(npp.dat.transect)
# sand.dat.section              <- naples.filter(sand.dat.section)
# sand.dat.transect             <- naples.filter(sand.dat.transect)
# section.dat                   <- naples.filter(section.dat)
# transect.dat                  <- naples.filter(transect.dat)
# understory.section.dat.wide   <- naples.filter(understory.section.dat.wide)
# understory.transect.dat.wide  <- naples.filter(understory.transect.dat.wide)
# urchin.dat.transect           <- naples.filter(urchin.dat.transect)

# # ------------------------------------------------------------------------------------
# # Check whether understory NPP is equivalent across treatments at the start of the continual removal
# start.dat <- section.dat %>%
#   dplyr::filter(time == 2008.125)
# 
# ggplot(data = start.dat, aes(x = treatment, y = understory.npp.season.gc.m2.day_at.section)) +
#   geom_point()
# 
# mod <- lm(understory.npp.season.gc.m2.day_at.section~treatment, data = start.dat)
# anova(mod)
# 
# # This tells us that we cannot discard the data prior to 2010.375
#
# # ------------------------------------------------------------------------------------
# # Remove first quarter from each treatment
# 
# drop.first.quarter <- function(x){
#   x1 <- x %>%
#     dplyr::filter(treatment == "Control" | treatment == "Annual") %>%
#     dplyr::filter(time > min(time))
# 
#   x2 <- x %>%
#     dplyr::filter(treatment == "Continual") %>%
#     dplyr::filter(time > min(time))
# 
#   x3 <- rbind(x1, x2)
# 
#   return(x3)
# }
# 
# # NOTE: Already did this for delta data!
# 
# biomass.dat.transect          <- drop.first.quarter(biomass.dat.transect)
# biomass.dat.transect.long     <- drop.first.quarter(biomass.dat.transect.long)
# biomass.dat.transect.wide     <- drop.first.quarter(biomass.dat.transect.wide)
# delta.npp.dat                 <- drop.first.quarter(delta.npp.dat)
# frond.dat.section             <- drop.first.quarter(frond.dat.section)
# frond.dat.transect            <- drop.first.quarter(frond.dat.transect)
# light.dat.transect            <- drop.first.quarter(light.dat.transect)
# npp.dat.section               <- drop.first.quarter(npp.dat.section)
# npp.dat.section.long          <- drop.first.quarter(npp.dat.section.long)
# npp.dat.section.plotting      <- drop.first.quarter(npp.dat.section.plotting)
# npp.dat.section.wide          <- drop.first.quarter(npp.dat.section.wide)
# npp.dat.transect              <- drop.first.quarter(npp.dat.transect)
# sand.dat.section              <- drop.first.quarter(sand.dat.section)
# sand.dat.transect             <- drop.first.quarter(sand.dat.transect)
# section.dat                   <- drop.first.quarter(section.dat)
# transect.dat                  <- drop.first.quarter(transect.dat)
# understory.transect.dat.wide  <- drop.first.quarter(understory.transect.dat.wide)
# urchin.dat.transect           <- drop.first.quarter(urchin.dat.transect)

# ------------------------------------------------------------------------------------
# Fix understory matrices to match
understory.transect.mat <- understory.transect.dat.wide %>%
  dplyr::select(-c(id:transect)) %>%
  as.matrix()
rownames(understory.transect.mat) <- understory.transect.dat.wide$id


# ====================================================================================
# CALCULATE DELTA-NPP ON AN ANNUAL SCALE
# ====================================================================================

# ------------------------------------------------------------------------------------
# Note that for annual removal, all data start in spring (20xx.375) and end in winter (20xx.125). Ready to calculate annual value. 
range(delta.npp.dat$time[delta.npp.dat$treatment == "Annual" & delta.npp.dat$site == "Mohawk"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Annual" & delta.npp.dat$site == "Isla Vista"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Annual" & delta.npp.dat$site == "Arroyo Quemado"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Annual" & delta.npp.dat$site == "Naples"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Annual" & delta.npp.dat$site == "Carpinteria"])

# Note that for continual removal, all data start in spring (20xx.375) and end in spring (20xx.375). To calculate an annual value, we must remove final spring or initial spring.
range(delta.npp.dat$time[delta.npp.dat$treatment == "Continual" & delta.npp.dat$site == "Mohawk"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Continual" & delta.npp.dat$site == "Arroyo Quemado"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Continual" & delta.npp.dat$site == "Naples"])
range(delta.npp.dat$time[delta.npp.dat$treatment == "Continual" & delta.npp.dat$site == "Carpinteria"])

# ------------------------------------------------------------------------------------
# Calculate mean delta-NPP for each year based on a 'kelp year', a  12 month period from spring of 20xx to winter of 20xx + 1

delta.npp.dat.annual <- delta.npp.dat %>%
  dplyr::mutate(time.temp = time - year,
                kelp.year = NA)

delta.npp.dat.annual$kelp.year[delta.npp.dat.annual$time.temp == 0.375] <- delta.npp.dat.annual$year[delta.npp.dat.annual$time.temp == 0.375] + 0.5
delta.npp.dat.annual$kelp.year[delta.npp.dat.annual$time.temp == 0.625] <- delta.npp.dat.annual$year[delta.npp.dat.annual$time.temp == 0.625] + 0.5
delta.npp.dat.annual$kelp.year[delta.npp.dat.annual$time.temp == 0.875] <- delta.npp.dat.annual$year[delta.npp.dat.annual$time.temp == 0.875] + 0.5
delta.npp.dat.annual$kelp.year[delta.npp.dat.annual$time.temp == 0.125] <- delta.npp.dat.annual$year[delta.npp.dat.annual$time.temp == 0.125] + 0.5 - 1

temp1 <- delta.npp.dat.annual %>%
  dplyr::filter(treatment == "Annual") %>%
  dplyr::select(year:time, kelp.year, everything(), -time.temp)

temp2 <- delta.npp.dat.annual %>%
  dplyr::filter(treatment == "Continual",
                site != "Naples",
                time >= 2010.375 & time < 2017.375) %>%
  dplyr::select(year:time, kelp.year, everything(), -time.temp)

temp3 <- delta.npp.dat.annual %>%
  dplyr::filter(treatment == "Continual",
                site == "Naples",
                time >= 2010.375 & time < 2016.375) %>%
  dplyr::select(year:time, kelp.year, everything(), -time.temp)

delta.npp.dat.annual <- rbind(temp1, temp2) %>%
  rbind(temp3)
rm(temp1, temp2, temp3)

# -----------------------------------------------------------------------------------------------------------------------------------
# Check

# Note that for annual removal, all data start in spring (20xx.375) and end in winter (20xx.125). Ready to calculate annual value. 
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Annual" & delta.npp.dat.annual$site == "Mohawk"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Annual" & delta.npp.dat.annual$site == "Isla Vista"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Annual" & delta.npp.dat.annual$site == "Arroyo Quemado"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Annual" & delta.npp.dat.annual$site == "Naples"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Annual" & delta.npp.dat.annual$site == "Carpinteria"]))/4

# Note that for continual removal, all data start in spring (20xx.375) and end in spring (20xx.375). To calculate an annual value, we must remove final spring or initial spring.
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Continual" & delta.npp.dat.annual$site == "Mohawk"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Continual" & delta.npp.dat.annual$site == "Arroyo Quemado"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Continual" & delta.npp.dat.annual$site == "Naples"]))/4
length(unique(delta.npp.dat.annual$time[delta.npp.dat.annual$treatment == "Continual" & delta.npp.dat.annual$site == "Carpinteria"]))/4

# -----------------------------------------------------------------------------------------------------------------------------------
# Calculate number of days per season
no.days.winter = (365 - 355) + (79 - 0)
no.days.spring = 172-79
no.days.summer = 265-172
no.days.autumn = 355-265

# Check
no.days.winter + no.days.spring + no.days.summer + no.days.autumn == 365

# -----------------------------------------------------------------------------------------------------------------------------------
# Calculate the total NPP per season
delta.npp.dat.annual$no.days.season <- NA
delta.npp.dat.annual$no.days.season[delta.npp.dat.annual$season == "1-Winter"] <- no.days.winter
delta.npp.dat.annual$no.days.season[delta.npp.dat.annual$season == "2-Spring"] <- no.days.spring
delta.npp.dat.annual$no.days.season[delta.npp.dat.annual$season == "3-Summer"] <- no.days.summer
delta.npp.dat.annual$no.days.season[delta.npp.dat.annual$season == "4-Autumn"] <- no.days.autumn

# Calculate annual summary metrics
delta.npp.dat.annual <- delta.npp.dat.annual %>%
  dplyr::mutate(delta.npp.season_gC.m2.season = delta.npp * no.days.season) %>%
  dplyr::rename(delta.npp.daily_gC.m2.d = delta.npp) %>%
  group_by(kelp.year, site, npp.component, treatment) %>%
  dplyr::summarize(no.quarters = length(delta.npp.daily_gC.m2.d),
                   annual.delta.npp_kgC.m2.y = sum(delta.npp.season_gC.m2.season, na.rm = F) / 1000,
                   daily.mean.delta.npp_gC.m2.d = mean(delta.npp.daily_gC.m2.d, na.rm = F),
                   daily.sd.delta.npp_gC.m2.d = sd(delta.npp.daily_gC.m2.d, na.rm = F)
                   ) %>%
  dplyr::mutate(daily.se.delta.npp = daily.sd.delta.npp_gC.m2.d / sqrt(no.quarters)) %>%
  dplyr::select(-daily.sd.delta.npp_gC.m2.d) %>%
  ungroup() %>%
  dplyr::filter(no.quarters == 4) %>% # Remove data for which there is a "yearly" observation for < 4 quarters. Should be none
  droplevels()

delta.npp.dat.annual$site <- factor(delta.npp.dat.annual$site,
                                    levels = c("Carpinteria", "Naples", "Arroyo Quemado", "Isla Vista", "Mohawk"))

delta.npp.dat.annual$treatment <- factor(delta.npp.dat.annual$treatment,
                                    levels = c("Continual", "Annual"))


# ====================================================================================
# CALCULATE ANNUAL SUM OR MEAN TRANSECT-SCALE DATA
# ====================================================================================

# ------------------------------------------------------------------------------------
# Calculate mean delta-NPP for each year based on a 'kelp year', a  12 month period from spring of 20xx to winter of 20xx + 1

transect.dat.annual <- transect.dat %>%
  dplyr::mutate(time.temp = time - year,
                kelp.year = NA)

transect.dat.annual$kelp.year[transect.dat.annual$time.temp == 0.375] <- transect.dat.annual$year[transect.dat.annual$time.temp == 0.375] + 0.5
transect.dat.annual$kelp.year[transect.dat.annual$time.temp == 0.625] <- transect.dat.annual$year[transect.dat.annual$time.temp == 0.625] + 0.5
transect.dat.annual$kelp.year[transect.dat.annual$time.temp == 0.875] <- transect.dat.annual$year[transect.dat.annual$time.temp == 0.875] + 0.5
transect.dat.annual$kelp.year[transect.dat.annual$time.temp == 0.125] <- transect.dat.annual$year[transect.dat.annual$time.temp == 0.125] + 0.5 - 1

temp1 <- transect.dat.annual %>%
  dplyr::filter(treatment == "Control") %>%
  dplyr::select(year:time, kelp.year, everything(), -month, -date, -id, -time.temp)

temp2 <- transect.dat.annual %>%
  dplyr::filter(treatment == "Annual") %>%
  dplyr::select(year:time, kelp.year, everything(), -month, -date, -id, -time.temp)

temp3 <- transect.dat.annual %>%
  dplyr::filter(treatment == "Continual",
                site != "Naples",
                time >= 2010.375 & time < 2017.375) %>%
  dplyr::select(year:time, kelp.year, everything(), -month, -date, -id, -time.temp)

temp4 <- transect.dat.annual %>%
  dplyr::filter(treatment == "Continual",
                site == "Naples",
                time >= 2010.375 & time < 2016.375) %>%
  dplyr::select(year:time, kelp.year, everything(), -month, -date, -id, -time.temp)

transect.dat.annual <- rbind(temp1, temp2) %>%
  rbind(temp3) %>%
  rbind(temp4)
rm(temp1, temp2, temp3, temp4)

# -----------------------------------------------------------------------------------------------------------------------------------
# Check

# Note that for control, all data start in spring (20xx.375) and end in winter (20xx.125). Ready to calculate annual value. 
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Control" & transect.dat.annual$site == "Mohawk"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Control" & transect.dat.annual$site == "Isla Vista"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Control" & transect.dat.annual$site == "Arroyo Quemado"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Control" & transect.dat.annual$site == "Naples"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Control" & transect.dat.annual$site == "Carpinteria"]))/4

# Note that for annual removal, all data start in spring (20xx.375) and end in winter (20xx.125). Ready to calculate annual value. 
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Annual" & transect.dat.annual$site == "Mohawk"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Annual" & transect.dat.annual$site == "Isla Vista"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Annual" & transect.dat.annual$site == "Arroyo Quemado"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Annual" & transect.dat.annual$site == "Naples"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Annual" & transect.dat.annual$site == "Carpinteria"]))/4

# Note that for continual removal, all data start in spring (20xx.375) and end in spring (20xx.375). To calculate an annual value, we must remove final spring or initial spring.
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Continual" & transect.dat.annual$site == "Mohawk"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Continual" & transect.dat.annual$site == "Arroyo Quemado"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Continual" & transect.dat.annual$site == "Naples"]))/4
length(unique(transect.dat.annual$time[transect.dat.annual$treatment == "Continual" & transect.dat.annual$site == "Carpinteria"]))/4

# -----------------------------------------------------------------------------------------------------------------------------------
# Calculate the total NPP per season
transect.dat.annual$no.days.season <- NA
transect.dat.annual$no.days.season[transect.dat.annual$season == "1-Winter"] <- no.days.winter
transect.dat.annual$no.days.season[transect.dat.annual$season == "2-Spring"] <- no.days.spring
transect.dat.annual$no.days.season[transect.dat.annual$season == "3-Summer"] <- no.days.summer
transect.dat.annual$no.days.season[transect.dat.annual$season == "4-Autumn"] <- no.days.autumn

# Calculate annual summary metrics
transect.dat.annual <- transect.dat.annual %>%
  dplyr::select(-depth.mllw.m, -mean.bottom.par, -sqrt.mean.bottom.par, -n.bottom.par, -fronds.per.transect, -plants.per.transect,
                -fronds.per.m2, -plants.per.m2, -fronds.per.plant, -sqrt.giant.kelp.dry.gm2, -understory.dry.gm2, -sqrt.understory.dry.gm2,
                -urchin.dens.no.m2, -urchin.sfdm.gm2) %>%
  dplyr::mutate(giant.kelp.npp_gC.m2.season = giant.kelp.npp.season.gc.m2.day * no.days.season,
                understory.npp_gC.m2.season = understory.npp.season.gc.m2.day * no.days.season,
                total.npp_gC.m2.season      = total.npp.season.gc.m2.day * no.days.season) %>%
  group_by(kelp.year, site, treatment) %>%
  dplyr::summarize(no.quarters = length(giant.kelp.npp_gC.m2.season),
                   giant.kelp.npp_kgC.m2.y = sum(giant.kelp.npp_gC.m2.season, na.rm = F) / 1000,
                   understory.npp_kgC.m2.y = sum(understory.npp_gC.m2.season, na.rm = F) / 1000,
                   total.npp_kgC.m2.y      = sum(total.npp_gC.m2.season, na.rm = F) / 1000,
                   mean.giant.kelp.dry.gm2 = mean(giant.kelp.dry.gm2, na.rm = F),
                   mean.urchin.dry.gm2     = mean(urchin.dry.gm2, na.rm = F),
                   mean.percent.sand       = mean(percent.sand, na.rm = F)
                   
  ) %>%
  ungroup() %>%
  dplyr::filter(no.quarters == 4) %>% # Remove data for which there is a "yearly" observation for < 4 quarters. Should be none
  droplevels()

transect.dat.annual$site <- factor(transect.dat.annual$site,
                                    levels = c("Carpinteria", "Naples", "Arroyo Quemado", "Isla Vista", "Mohawk"))

transect.dat.annual$treatment <- factor(transect.dat.annual$treatment,
                                         levels = c("Control", "Annual", "Continual"))


# ====================================================================================
# QAQC: PLOT DISTRIBUTION OF SAMPLING OVER SPACE AND TIME  
# ====================================================================================

ggplot(data = light.dat.transect, aes(x = time, y = transect, color = treatment)) +
  geom_point() +
  facet_wrap(~ site, scales = 'free_y') +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2))+
  ggtitle("Spatiotemporal sampling of benthic PAR data") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = npp.dat.transect, aes(x = time, y = transect, color = treatment)) +
  geom_point() +
  facet_wrap(~ site, scales = 'free_y') +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2)) +
  ggtitle("Spatiotemporal sampling of NPP data") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = delta.npp.dat, aes(x = time, y = treatment, color = treatment)) +
  geom_point() +
  facet_wrap(~ site, scales = 'free_y') +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2)) +
  ggtitle("Spatiotemporal sampling of delta NPP data") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = F)

ggplot(data = delta.npp.dat.annual, aes(x = kelp.year, y = treatment, color = treatment)) +
  geom_point() +
  facet_wrap(~ site, scales = 'free_y') +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2)) +
  ggtitle("Spatiotemporal sampling of annual delta NPP data") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = F)

ggplot(data = transect.dat.annual, aes(x = kelp.year, y = treatment, color = treatment)) +
  geom_point() +
  facet_wrap(~ site, scales = 'free_y') +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2)) +
  ggtitle("Spatiotemporal sampling of annual NPP data") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = F)


# =================================================================================================================================
# CHECK FOR COMPENSATION: ARE THERE TIMES WHEN DECREASES IN GIANT KELP NPP ARE MATCHED OR EXCEEDED BY INCREASES IN UNDERSTORY NPP?
# =================================================================================================================================
delta.npp.dat.annual.wide <- delta.npp.dat.annual %>%
  dplyr::select(-daily.se.delta.npp, -daily.mean.delta.npp_gC.m2.d, -no.quarters) %>%
  pivot_wider(id_cols = c(kelp.year, site, treatment),
              names_from = npp.component,
              values_from = annual.delta.npp_kgC.m2.y)

ribbon.dat <- data.frame(giant.kelp = c(-1.75, -0.5),
                         lwr = c(1.35, 0.1),
                         upr = c(2.15, 0.9),
                         understory = c(0, 0))


comp.dat <- delta.npp.dat.annual.wide %>%
  dplyr::filter(giant.kelp < -0.5 & giant.kelp > -1.75,
                understory > 0.4 & understory < 1.5)
comp.dat$ratio = comp.dat$understory / (-1 * comp.dat$giant.kelp)

delta.npp.dat.annual.wide$ratio = delta.npp.dat.annual.wide$understory / (-1 * delta.npp.dat.annual.wide$giant.kelp)


p <- ggplot(data = delta.npp.dat.annual.wide, aes(x = giant.kelp, y = understory)) +
  geom_ribbon(data = ribbon.dat, aes(ymin = lwr, ymax = upr),
           fill = "#FFE3EB",
           alpha = 0.6,
           color = "#FFE3EB") +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(23, 24), name = "Treatment") +
  geom_hline(yintercept = 0, size = 0.5) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_abline(intercept = 0, slope = -1, linetype = 'dashed') +
  geom_point(aes(shape = treatment, fill = treatment), size = 2) +
  geom_point(data = comp.dat, color = 'black', size = 0.25) +
  scale_x_continuous(limits = c(-3.25, 1), breaks = -3:3) +
  scale_y_continuous(limits = c(-1, 3), breaks = -3:3) +
  xlab(expression(atop("Giant kelp annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  ylab(expression(atop("Understory annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  ggtitle("NPP compensation by understory macroalgae") +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_line(color="black", size = 0.1),
        panel.grid.minor = element_line(color="black", size = 0.1),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
p
#ggsave(p, file = "Figures/Compensation plot.pdf", width = 6, height = 6)


# ====================================================================================
# CALCULATE RANGE OF TOTAL MACROALGAL NPP
# ====================================================================================

npp.summary <- transect.dat.annual$total.npp_kgC.m2.y[transect.dat.annual$treatment == "Control"]
hist(npp.summary)

k <- 1 + 3.322 * log(length(npp.summary))
h <- (max(npp.summary) - min(npp.summary)) / k

ggplot(data = transect.dat.annual[transect.dat.annual$treatment == "Control", ], 
       aes(x = total.npp_kgC.m2.y, y = ..density..)) +
  geom_histogram(breaks=seq(0, 4.5, by = 0.5), color = 'black', fill = 'grey') +
  ylab('Proportion of observations') +
  xlab(expression(paste("Total macroalgal NPP (kg C/", "m"^2, "/y)", sep = ""))) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5), expand = expansion(mult = c(0.003, 0.1))) +
  theme_classic() +
  theme(text = element_text(size = 16),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1)
#ggsave(filename = paste0(here::here("/Figures/"),  "total.npp.in.control.plots", ".pdf"), height = 4, width = 4)


# ====================================================================================
# PLOT 'HABITAT QUALITY' FOR TOTAL MACROALGAE IN CONTROL PLOTS
# ====================================================================================

# ----------------------------------------------------------------------------------------------------------------------
quality.dat.tot <- transect.dat %>%
  dplyr::filter(treatment == "Control") %>%
  dplyr::select(time, season,year, site, treatment, understory.dry.gm2, giant.kelp.dry.gm2, percent.sand, urchin.dens.no.m2) %>%
  dplyr::mutate(understory.dry.kgm2 = understory.dry.gm2 / 1000,
                giant.kelp.dry.kgm2 = giant.kelp.dry.gm2 / 1000) %>%
  dplyr::select(-understory.dry.gm2, -giant.kelp.dry.gm2) %>%
  dplyr::mutate(total.algae.dry.kgm2 = understory.dry.kgm2 + giant.kelp.dry.kgm2) 

algal.biomass.dat <- quality.dat.tot %>%
  dplyr::select(-percent.sand, -urchin.dens.no.m2) %>%
  pivot_longer(cols = c('understory.dry.kgm2', 'giant.kelp.dry.kgm2', 'total.algae.dry.kgm2'),
               names_to = "algal.group", 
               values_to = "biomass.dry.kgm2") %>%
  dplyr::group_by(site, treatment, algal.group) %>%
  dplyr::summarise(mean = mean(biomass.dry.kgm2),
                   sd = sd(biomass.dry.kgm2, na.rm = T),
                   n = length(!is.na(biomass.dry.kgm2))) %>%
  dplyr::mutate(se = sd / sqrt(n),
                ci = se * 1.96) %>%
  ungroup()

algal.biomass.dat$site <- factor(algal.biomass.dat$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

#algal.biomass.dat$vadj <- c(algal.biomass.dat$mean[2], 0, algal.biomass.dat$mean[4], 0, algal.biomass.dat$mean[6], 0, algal.biomass.dat$mean[8], 0, algal.biomass.dat$mean[10], 0, algal.biomass.dat$mean[12], 0, algal.biomass.dat$mean[14], 0)

# ----------------------------------------------------------------------------------------------------------------------
# Stacked bar chart of algal biomass
ggplot(data = algal.biomass.dat[algal.biomass.dat$algal.group != 'total.algae.dry.kgm2',], aes(x = site, y = mean)) +
  geom_bar(aes(fill = algal.group), position="stack", stat="identity", color = 'black') +
  #geom_errorbar(aes(ymax = mean + ci + vadj, ymin = mean - ci + vadj), width = 0) +
  xlab("\nSite") +
  ylab("Biomass (kg dry/m2)") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.2))) +
  theme(plot.title = element_text(size = 16),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/3,
        legend.position = c(0.88, 0.89),
        legend.background = element_blank()) +
  ggtitle("Macroalgal biomass in control plots (mean ± 95% CI)") +
  scale_fill_manual(values = c("darkgoldenrod3", "darkgreen"), name = "", labels = c("Giant kelp", "Understory"))
#ggsave(filename = paste0(here::here("/Figures/"),  "mean.algal.biomass_stack.barchart2", ".pdf"), height = 4, width = 7)

# ----------------------------------------------------------------------------------------------------------------------
quality.dat.tot.wide <- quality.dat.tot 

quality.dat.tot <- quality.dat.tot.wide %>%
  dplyr::select(-giant.kelp.dry.kgm2, -understory.dry.kgm2, -year, -time, -season) %>%
  pivot_longer(cols = c('total.algae.dry.kgm2', 'percent.sand', 'urchin.dens.no.m2'),
               names_to = "response.variable", 
               values_to = "value")
quality.dat.tot$site <- factor(quality.dat.tot$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

quality.dat.tot$response.variable <- factor(quality.dat.tot$response.variable, levels = c('urchin.dens.no.m2', 'percent.sand','total.algae.dry.kgm2'),
                                           labels = c(
                                             expression(paste("(A) Sea urchin density (no./m"^2, ")", sep = "")),
                                             expression(paste("(B) Sand cover (", '%', ")", sep = "")),
                                             expression(paste("(D) Macroalgal biomass (kg dry/m"^2, ")", sep = ""))
                                           ))

# ----------------------------------------------------------------------------------------------------------------------
ggplot(data = quality.dat.tot, aes(x = site, y = value)) +
  facet_wrap(~ response.variable, ncol = 1, scales = 'free_y', labeller = label_parsed) +
  stat_summary(fun.data = 'mean_cl_boot', aes(fill = response.variable), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  scale_fill_manual(values = c( 'mediumpurple1', 'moccasin','darkgoldenrod3'), name = "") +
  guides(fill = F) +
  xlab("\nSite") +
  ylab("Value") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.2))) +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/3) 
#ggsave(filename = paste0(here::here("/Figures/"),  "habitat.quality.plot1", ".pdf"), height = 9, width = 7)

urchin.col <- '#cc98a2'
a <- ggplot(data = quality.dat.tot %>% dplyr::filter(response.variable == levels(quality.dat.tot$response.variable)[1]), 
            aes(x = site, y = value)) +
  stat_summary(fun.data = 'mean_cl_boot', fill = urchin.col, geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  ylab( expression(paste("Sea urchin density (no./m"^2, ")", sep = ""))) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_x_discrete(labels = NULL) +
  theme_classic() +
  ggtitle("") +
  #ggtitle("(A) Herbivory by sea urchins") +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        aspect.ratio = 1/3) 
a

sand.col <- '#e0d7ac'
b <- ggplot(data = quality.dat.tot %>% dplyr::filter(response.variable == levels(quality.dat.tot$response.variable)[2]), 
            aes(x = site, y = value)) +
  stat_summary(fun.data = 'mean_cl_boot', fill = sand.col, geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  ylab("Sand cover (%)") +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1)), breaks = seq(0, 30, by = 10)) +
  coord_cartesian(ylim = c(0, 28)) +
  scale_x_discrete(labels = NULL) +
  theme_classic() +
  ggtitle("") +
  #ggtitle("(B) Sand burial and abrasion") +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        aspect.ratio = 1/3) 
b

kelp.col <- '#78bd8d'
c <- ggplot(data = quality.dat.tot %>% dplyr::filter(response.variable == levels(quality.dat.tot$response.variable)[3]), 
            aes(x = site, y = value)) +
  stat_summary(fun.data = 'mean_cl_boot', fill = kelp.col, geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  ylab(expression(paste("Macroalgal biomass (kg dry/m"^2, ")", sep = "")))+
  xlab("") +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.03)), breaks = seq(0, 1.2, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_classic() +
  ggtitle("") +
  #ggtitle("(C) Abundance of canopy-forming kelp species") +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        aspect.ratio = 1/3) 
c

(a / b / c)
#ggsave(filename = paste0(here::here("/Figures/"),  "habitat.quality.plot - experiment years", ".pdf"), height = 9.5, width = 7.75)

# ------------------------------------------------------------------------------------
# Test whether variables differ among sites -- urchins
quality.dat.tot.wide$time.factor <- factor(quality.dat.tot.wide$time, levels = seq(min(quality.dat.tot.wide$time),
                                                                                 max(quality.dat.tot.wide$time),
                                                                                 by = 0.25))

# site.urchin.mod <- glmmTMB(urchin.dens.no.m2 ~ site + poly(time, 2) +
#                            ar1(time.factor + 0 | site),
#                          family = ziGamma(link = "log"),
#                          ziformula = ~ 1,
#                          #dispformula = ~ factor(year),
#                          data = quality.dat.tot.wide)

site.urchin.mod <- glmmTMB(sqrt(urchin.dens.no.m2) ~ site + poly(time, 2) + 
                             ar1(time.factor + 0 | site),
                           family = gaussian(link = "identity"),
                           dispformula = ~ factor(year),
                           data = quality.dat.tot.wide)
AIC(site.urchin.mod)
summary(site.urchin.mod)
Anova(site.urchin.mod)
r2(site.urchin.mod)

em <- emmeans(site.urchin.mod, ~ site)
plot(em)

pairs <- pairs(em, simple = "site", reverse = TRUE, adjust = 'none')
pairs
plot(pairs)

# Validate site.urchin.mod
simulationOutput <- simulateResiduals(site.urchin.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

quality.dat.tot.wide$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(quality.dat.tot.wide$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(quality.dat.tot.wide$resid, xlab = "Pearson residuals", main = "")
plot(predict(site.urchin.mod, type = 'response'), quality.dat.tot.wide$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$time, quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$year), quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$season), quality.dat.tot.wide$resid, xlab = "Season", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$site), quality.dat.tot.wide$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(sqrt(quality.dat.tot.wide$urchin.dens.no.m2),  predict(site.urchin.mod, type = 'response'), ylab = "sqrt(Sea urchin density)", xlab = "Predicted values"); abline(c(0,1))

# Calculate ACF and PACF for each site
quality.dat.tot.wide$resid <- resid(site.urchin.mod)

acf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  acf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
      lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  pacf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
       lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
}
)

acf.dat <- data.frame((acf.dat))
pacf.dat <- data.frame((pacf.dat))

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "site", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "site", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "site"))

# Calculate critical value (based on the lowest length of kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(quality.dat.tot.wide[quality.dat.tot.wide$site == "Isla Vista", ]$year)))

# site ACF by site
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by site") +
  facet_wrap(~ site) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# site average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across sites") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Test whether variables differ among sites -- sand
site.sand.mod <- glmmTMB(percent.sand/100 ~ site, #+ 
                         #ar1(time.factor + 0 | site),
                         family = beta_family(link = "logit"),
                         ziformula = ~ 1,
                         #dispformula = ~ site,
                         data = quality.dat.tot.wide)
AIC(site.sand.mod)
summary(site.sand.mod)
Anova(site.sand.mod)
r2(site.sand.mod)

em <- emmeans(site.sand.mod, ~ site)
plot(em)

pairs <- pairs(em, simple = "site", reverse = TRUE, adjust = 'none')
pairs
plot(pairs)

# Validate site.sand.mod
simulationOutput <- simulateResiduals(site.sand.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

quality.dat.tot.wide$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(quality.dat.tot.wide$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(quality.dat.tot.wide$resid, xlab = "Pearson residuals", main = "")
plot(predict(site.sand.mod, type = 'response'), quality.dat.tot.wide$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$time, quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$year), quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$season), quality.dat.tot.wide$resid, xlab = "Season", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$site), quality.dat.tot.wide$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$percent.sand,  predict(site.sand.mod, type = 'response'), ylab = "Sea urchin density", xlab = "Predicted values"); abline(c(0,1))

# ----------------------------------------------------------------------------------------------------------
# Calculate ACF and PACF for each site
quality.dat.tot.wide$resid <- resid(site.sand.mod)

acf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  acf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
      lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  pacf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
       lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
}
)

acf.dat <- data.frame((acf.dat))
pacf.dat <- data.frame((pacf.dat))

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "site", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "site", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "site"))

# Calculate critical value (based on the lowest length of kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(quality.dat.tot.wide[quality.dat.tot.wide$site == "Isla Vista", ]$year)))

# site ACF by site
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by site") +
  facet_wrap(~ site) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# site average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across sites") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Test whether variables differ among sites -- macroalgal biomass
hist(quality.dat.tot.wide$total.algae.dry.kgm2)

site.tot.mod <- glmmTMB(sqrt(total.algae.dry.kgm2) ~ site + 
                         ar1(time.factor + 0 | site),
                       family = gaussian(link = "identity"),
                       #family = Gamma(link = "log"),
                       #ziformula = ~ 1,
                       #dispformula = ~ site,
                       data = quality.dat.tot.wide)
AIC(site.tot.mod)
summary(site.tot.mod)
Anova(site.tot.mod)
r2(site.tot.mod)

em <- emmeans(site.tot.mod, ~ site)
plot(em)

pairs <- pairs(em, simple = "site", reverse = TRUE, adjust = 'none')
pairs
plot(pairs)

# Validate site.tot.mod
simulationOutput <- simulateResiduals(site.tot.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

quality.dat.tot.wide$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(quality.dat.tot.wide$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(quality.dat.tot.wide$resid, xlab = "Pearson residuals", main = "")
plot(predict(site.tot.mod, type = 'response'), quality.dat.tot.wide$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$time, quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$year), quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$season), quality.dat.tot.wide$resid, xlab = "Season", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$site), quality.dat.tot.wide$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$total.algae.dry.kgm2,  predict(site.tot.mod, type = 'response'), ylab = "Total algal biomass", xlab = "Predicted values"); abline(c(0,1))

# ----------------------------------------------------------------------------------------------------------
# Calculate ACF and PACF for each site
quality.dat.tot.wide$resid <- resid(site.tot.mod)

acf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  acf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
      lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  pacf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
       lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
}
)

acf.dat <- data.frame((acf.dat))
pacf.dat <- data.frame((pacf.dat))

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "site", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "site", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "site"))

# Calculate critical value (based on the lowest length of kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(quality.dat.tot.wide[quality.dat.tot.wide$site == "Isla Vista", ]$year)))

# site ACF by site
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by site") +
  facet_wrap(~ site) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# site average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across sites") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Create linear model of total macroalgal biomass as a function of urchins and sand

#qual.mod <- lm(sqrt(total.algae.dry.kgm2) ~ urchin.dens.no.m2 + percent.sand, 
#               data = quality.dat.tot.wide)

# qual.mod <- glmmTMB(sqrt(total.algae.dry.kgm2) ~ urchin.dens.no.m2 + percent.sand, #+
#                     #ar1(time.factor + 0 | site),
#                     family = gaussian(link = "identity"),
#                     data = quality.dat.tot.wide)

qual.mod <- glmmTMB(total.algae.dry.kgm2 ~ urchin.dens.no.m2 + percent.sand, #+
                    #ar1(time.factor + 0 | site),
                    family = Gamma(link = "log"),
                    data = quality.dat.tot.wide)

AIC(qual.mod)
Anova(qual.mod)
summary(qual.mod)

r2(qual.mod)

qual.mod.glm <- glm(total.algae.dry.kgm2 ~ urchin.dens.no.m2 + percent.sand, 
                    family = Gamma(link = "log"),
                    data = quality.dat.tot.wide)
r2(qual.mod)

# Pseudo-R2 (squared correlation between the response and the predicted value; after Bolker)
cor(quality.dat.tot.wide$giant.kelp.dry.kgm2, predict(qual.mod, type="response"))^2

p <- ggpredict(qual.mod, "urchin.dens.no.m2", type = "fe")
plot(p)

p <- ggpredict(qual.mod, "percent.sand", type = "fe")
plot(p)

# Validate qual.mod
simulationOutput <- simulateResiduals(qual.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

quality.dat.tot.wide$resid <- simulationOutput$scaledResiduals
#quality.dat.tot.wide$resid <- resid(qual.mod)

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(quality.dat.tot.wide$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(quality.dat.tot.wide$resid, xlab = "Pearson residuals", main = "")
plot(predict(qual.mod, type = 'response'), quality.dat.tot.wide$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$time, quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$year), quality.dat.tot.wide$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$season), quality.dat.tot.wide$resid, xlab = "Season", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$site), quality.dat.tot.wide$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(quality.dat.tot.wide$treatment), quality.dat.tot.wide$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(quality.dat.tot.wide$giant.kelp.dry.kgm2,  predict(qual.mod, type = 'response'), ylab = "Giant kelp canopy biomass", xlab = "Predicted values"); abline(c(0,1))

# ----------------------------------------------------------------------------------------------------------
# Calculate ACF and PACF for each site
quality.dat.tot.wide$resid <- resid(qual.mod)

acf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  acf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
      lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(quality.dat.tot.wide$site), function(x){
  pacf(quality.dat.tot.wide$resid[quality.dat.tot.wide$site == x], 
       lag.max = length(unique(quality.dat.tot.wide$year)) / 2, plot = FALSE)$acf
}
)

acf.dat <- data.frame((acf.dat))
pacf.dat <- data.frame((pacf.dat))

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "site", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "site", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "site"))

# Calculate critical value (based on the lowest length of kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(quality.dat.tot.wide[quality.dat.tot.wide$site == "Isla Vista", ]$year)))

# site ACF by site
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by site") +
  facet_wrap(~ site) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# site average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across sites") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") + 
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals
preds.df <- expand.grid(urchin.dens.no.m2 = seq(min(quality.dat.tot.wide$urchin.dens.no.m2), 
                                                max(quality.dat.tot.wide$urchin.dens.no.m2),
                                                by = (max(quality.dat.tot.wide$urchin.dens.no.m2) - min(quality.dat.tot.wide$urchin.dens.no.m2)) / 50),
                        percent.sand = seq(min(quality.dat.tot.wide$percent.sand), 
                                           max(quality.dat.tot.wide$percent.sand),
                                           by = (max(quality.dat.tot.wide$percent.sand) - min(quality.dat.tot.wide$percent.sand)) / 50)
)

# Predict new data
preds.df$total.algae.dry.kgm2   <- predict(qual.mod, newdata = preds.df, type = "response", se.fit = FALSE)
preds.df$se              <- predict(qual.mod, newdata = preds.df, type = "response", se.fit = TRUE)$se.fit
preds.df$upr.ci          <- preds.df$total.algae.dry.kgm2 + (1.96 * preds.df$se)
preds.df$lwr.ci          <- preds.df$total.algae.dry.kgm2 - (1.96 * preds.df$se)
 
# Summarize predictions as a function of giant kelp frond density
preds.avg.urchins <- preds.df %>%
  dplyr::group_by(urchin.dens.no.m2) %>%
  dplyr::summarise_at(vars(total.algae.dry.kgm2:lwr.ci),
                      mean,
                      na.rm = T) %>%
  dplyr::ungroup()

preds.avg.sand <- preds.df %>%
  dplyr::group_by(percent.sand) %>%
  dplyr::summarise_at(vars(total.algae.dry.kgm2:lwr.ci),
                      mean,
                      na.rm = T) %>%
  dplyr::ungroup()


# ------------------------------------------------------------------------------------
# Plot predictions
ggplot(data = quality.dat.tot.wide, aes(x = urchin.dens.no.m2, y = total.algae.dry.kgm2)) +
  geom_ribbon(data = preds.avg.urchins, aes(y = total.algae.dry.kgm2, ymin = lwr.ci, ymax = upr.ci), fill = urchin.col, alpha = 0.33) +
  geom_point(size = 2, shape = 21, fill = urchin.col) +
  geom_line(data = preds.avg.urchins, aes(y = total.algae.dry.kgm2), color = 'black', size = 1.25) +
  xlab(expression("Sea urchin density (no./"*"m"^2*")")) +
  ylab(expression("Macroalgal biomass (kg dry/"*"m"^2*")")) +
  theme_classic() +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border= element_blank(), #element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = paste0(here::here("/Figures/"),  "urchin.effect.on.total.algal.biomass", ".pdf"), height = 4, width = 4)
#ggsave(filename = paste0(here::here("/Figures/"),  "urchin.effect.on.total.algal.biomass", ".svg"), height = 4, width = 4)

ggplot(data = quality.dat.tot.wide, aes(x = percent.sand, y = total.algae.dry.kgm2)) +
  geom_ribbon(data = preds.avg.sand, aes(y = percent.sand, ymin = lwr.ci, ymax = upr.ci), fill = sand.col, alpha = 0.33) +
  geom_point(size = 2, shape = 21, fill = sand.col) +
  geom_line(data = preds.avg.sand, aes(y = total.algae.dry.kgm2), color = 'black', size = 1.25) +
  xlab(expression("Sand cover (%)")) +
  ylab(expression("Macroalgal biomass (kg dry/"*"m"^2*")")) +
  theme_classic() +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border= element_blank(), #element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = paste0(here::here("/Figures/"),  "sand.effect.on.total.algal.biomass", ".pdf"), height = 4, width = 4)
#ggsave(filename = paste0(here::here("/Figures/"),  "sand.effect.on.total.algal.biomass", ".svg"), height = 4, width = 4)

# ----------------------------------------------------------------------------------------------------------
# Plot predicted habitat quality, rescaled to a max. of 1
quality.dat.tot.wide$preds <- predict(qual.mod, newdata = quality.dat.tot.wide, type = "response")

quality.dat.tot.wide$preds <- quality.dat.tot.wide$preds / max(quality.dat.tot.wide$preds )

quality.dat.tot.wide$site <- factor(quality.dat.tot.wide$site, 
                                   levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

blend.col <- '#D6B8A7'

d <- ggplot(data = quality.dat.tot.wide, aes(x = site, y = preds)) +
  stat_summary(fun.data = 'mean_cl_boot', fill = blend.col, geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  ylab("Habitat quality index") +
  xlab("\nSite") +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.03)), breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  ggtitle("") +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/3)   
d
#ggsave(filename = paste0(here::here("/Figures/"),  "habitat.qual.for.macroalgae-continuous.predictors", ".pdf"), height = 5, width = 7)

(c / a / b / d)
#ggsave(filename = paste0(here::here("/Figures/"),  "habitat.quality.plot - experiment years - 2", ".pdf"), height = 11.75, width = 7.75)
#ggsave(filename = paste0(here::here("/Figures/"),  "habitat.quality.plot - experiment years - 2", ".svg"), height = 11.75, width = 7.75)

# ------------------------------------------------------------------------------------
# Test whether variables differ among sites -- predicted habitat quality
site.preds.mod <- lm(preds ~ site, #+ 
                     #ar1(time.factor + 0 | site),
                     #family = ziGamma(link = "log"),
                     #ziformula = ~ 1,
                     #dispformula = ~ site,
                     data = quality.dat.tot.wide)

AIC(site.preds.mod)
summary(site.preds.mod)
Anova(site.preds.mod)
r2(site.preds.mod)

# Pseudo-R2 (squared correlation between the response and the predicted value; after Bolker)
cor(quality.dat.tot.wide$preds, predict(site.preds.mod, type="response"))^2

em <- emmeans(site.preds.mod, ~ site)
plot(em)

pairs <- pairs(em, simple = "site", reverse = TRUE, adjust = 'none')
pairs
plot(pairs)

# Validate site.preds.mod
simulationOutput <- simulateResiduals(site.preds.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)


# ====================================================================================
# EXPLORATORY PLOTTING
# ====================================================================================

# ------------------------------------------------------------------------------------
# Explore NPP data

# Boxplots of NPP by treatment
ggplot(data = npp.dat.section.plotting, aes(x = treatment, y = guild.npp.season.gc.m2.day, fill = guild)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 12)) +
  scale_fill_manual(values = c("grey", "darkgoldenrod", "darkolivegreen2")) +
  facet_wrap(~ site) +
  theme_classic()
#ggsave(filename = paste0(here::here("/Figures/"),  "npp.boxplot.by.treatment", ".pdf"), height = 5, width = 7)

# Mean +/- 95% CI of NPP by treatment
ggplot(data = npp.dat.section.plotting, aes(x = treatment, y = guild.npp.season.gc.m2.day, color = guild, group = guild)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1) +
  stat_summary(fun = "mean", geom = "line", size = 1) +
  ylab(expression(paste("NPP (g/", "m"^2, "/d)", sep = ""))) +
  xlab("Treatment") + 
  scale_color_manual(values = c("dimgrey", "darkgoldenrod", "darkolivegreen3")) +
  theme_classic()
#ggsave(filename = paste0(here::here("/Figures/"),  "npp.mean.95ci.by.treatment", ".pdf"), height = 5, width = 7)

# Mean +/- 95% CI of NPP by treatment, season, and site 
ggplot(data = npp.dat.section.plotting, 
       aes(x = season, y = guild.npp.season.gc.m2.day, color = guild, group = guild)) +
  stat_summary(data = npp.dat.section.plotting[npp.dat.section.plotting$guild == "Total", ], 
               fun = "mean", geom = "area", alpha = 0.25) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  stat_summary(fun = "mean", geom = "line", size = 0.6) +
  ylab(expression(paste("NPP (g/", "m"^2, "/d)", sep = ""))) +
  scale_y_continuous(breaks = seq(0,10, 2)) +
  scale_color_manual(values = c("darkgoldenrod3","dimgrey",  "darkolivegreen3"), name = "") +
  theme_bw() +
  facet_grid(site ~ treatment) +
  theme(axis.text.x=element_text(angle=-45, hjust=0)) +
  theme(aspect.ratio = 1) +
  xlab("Season")
#ggsave(filename = paste0(here::here("/Figures/"),  "Seasonal NPP by treatment and site", ".pdf"), height = 10, width = 7)

# ggplot(data = npp.dat.section.plotting[npp.dat.section.plotting$guild != "Total" &
#                                  npp.dat.section.plotting$season == "3-Summer", ], 
#        aes(x = site, y = guild.npp.season.gc.m2.day, color = guild, group = guild)) +
#   stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
#   ylab(expression(paste("NPP (g/", "m"^2, "/d)", sep = ""))) +
#   scale_color_manual(values = c("darkgoldenrod", "darkolivegreen3")) +
#   theme_classic() +
#   facet_wrap(~ treatment, ncol = 1) +
#   ggtitle("Summer data")
#ggsave(filename = paste0(here::here("/Figures/"),  "npp.mean.by.treatment.site.summer", ".pdf"), height = 6, width = 7)

# Relationship between understory and giant kelp NPP
ggplot(data = transect.dat, aes(x = fronds.per.m2, y = understory.npp.season.gc.m2.day)) +
  geom_point() +
  xlab(expression(paste("Giant kelp frond density (no./m"^2, ")", sep = ""))) +
  ylab(expression(paste("Understory NPP (g/", "m"^2, "/d)", sep = ""))) +
  theme_classic() +
  ggtitle("Transect-scale data")

ggplot(data = section.dat, aes(x = giant.kelp.npp.season.gc.m2.day_at.section, y = understory.npp.season.gc.m2.day_at.section)) +
  geom_point(alpha = 0.2) +
  xlab(expression(paste("Giant kelp frond density (no./m"^2, ")", sep = ""))) +
  ylab(expression(paste("Understory NPP (g/", "m"^2, "/d)", sep = ""))) +
  theme_classic() +
  ggtitle("Section-scale data")
#ggsave(filename = paste0(here::here("/Figures/"),  "understory.vs.kelp.npp", ".pdf"), height = 6, width = 6)

# Relationship between understory and giant kelp NPP by season
ggplot(data = npp.dat.section, aes(x = giant.kelp.npp.season.gc.m2.day, y = understory.npp.season.gc.m2.day)) +
  geom_point() +
  xlab(expression(paste("Giant kelp NPP (g/", "m"^2, "/d)", sep = ""))) +
  ylab(expression(paste("Understory NPP (g/", "m"^2, "/d)", sep = ""))) +
  theme_classic() +
  facet_wrap(~ season)
#ggsave(filename = paste0(here::here("/Figures/"),  "understory.vs.kelp.npp.by.season", ".pdf"), height = 7, width = 7)

# Changes in understory NPP over time as it relates to giant kelp frond density over time
ggplot(data = section.dat, 
       aes(x = time, y = understory.npp.season.gc.m2.day_at.section)) +
  geom_point(aes(shape = treatment, fill = treatment), size = 2) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(22, 23, 24), name = "Treatment") +
  geom_smooth(method = 'gam', se = F) +
  facet_grid(treatment ~ site) +
  xlab("Year") +
  ylab(expression(paste("Understory NPP (g/", "m"^2, "/d)", sep = ""))) +
  theme_classic() +
  ggtitle("Section-scale data")

ggplot(data = transect.dat, 
       aes(x = time, y = understory.npp.season.gc.m2.day)) +
  geom_point(aes(shape = treatment, fill = treatment), size = 2) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(22, 23, 24), name = "Treatment") +
  geom_smooth(method = 'gam', se = F) +
  facet_grid(treatment ~ site) +
  xlab("Year") +
  ylab(expression(paste("Understory NPP (g/", "m"^2, "/d)", sep = ""))) +
  theme_classic() +
  ggtitle("Transect-scale data")

ggplot(data = transect.dat, 
       aes(x = time, y = fronds.per.m2)) +
  geom_point(aes(shape = treatment, fill = treatment), size = 2) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(22, 23, 24), name = "Treatment") +
  geom_smooth(method = 'lm', se = F) +
  facet_grid(treatment ~ site) +
  xlab("Year") +
  ylab(expression(paste("Giant kelp frond density (no.", "m"^2, ")", sep = ""))) +
  theme_classic() +
  ggtitle("Transect-scale data")

transect.dat$site <- factor(transect.dat$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))


# ====================================================================================
# FURTHER EXPLORATORY PLOTS (AUGUST 2020)
# ====================================================================================

# Re-order sites and treatments
transect.dat$treatment <- factor(transect.dat$treatment, levels = c("Control", "Annual", "Continual"))
transect.dat$site <- factor(transect.dat$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

transect.dat$treatment2 <- dplyr::recode(transect.dat$treatment, 
                                         "Control" = "Control", 
                                         "Annual" = "Annual giant kelp removal", 
                                         "Continual" = "Quarterly giant kelp removal")
transect.dat$treatment2 <- factor(transect.dat$treatment2, 
                                  levels = c("Control", "Annual giant kelp removal", "Quarterly giant kelp removal"))


transect.dat$treatment4 <- dplyr::recode(transect.dat$treatment, 
                                        "Control" = "Control", 
                                        "Annual" = "Annual", 
                                        "Continual" = "Quarterly")
transect.dat$treatment4 <- factor(transect.dat$treatment4, 
                                 levels = c("Control", "Annual", "Quarterly"))

section.dat$treatment <- factor(section.dat$treatment, levels = c("Control", "Annual", "Continual"))
section.dat$site <- factor(section.dat$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

section.dat$treatment2 <- dplyr::recode(section.dat$treatment, 
                                         "Control" = "Control", 
                                         "Annual" = "Annual giant kelp removal", 
                                         "Continual" = "Quarterly giant kelp removal")
section.dat$treatment2 <- factor(section.dat$treatment2, 
                                  levels = c("Control", "Annual giant kelp removal", "Quarterly giant kelp removal"))

section.dat$treatment4 <- dplyr::recode(section.dat$treatment, 
                                         "Control" = "Control", 
                                         "Annual" = "Annual", 
                                         "Continual" = "Quarterly")
section.dat$treatment4 <- factor(section.dat$treatment4, 
                                  levels = c("Control", "Annual", "Quarterly"))

# -------------------------------------------------------------------------------------------------
# Plot time series of giant kelp frond density at all 5 sites
ggplot(data = transect.dat, aes(x = time, y = fronds.per.m2)) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Control", ], fill = "#1b9e77", alpha = 0.2) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Annual giant kelp removal", ], fill = "#d95f02", alpha = 0.2) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Quarterly giant kelp removal", ], fill = "#7570b3", alpha = 0.2) +
  geom_line(aes(linetype = treatment2, color = treatment2), size = 0.75) +
  scale_linetype_manual(values = c("solid", "solid", "solid"), name = "") +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), name = "") +
  facet_wrap(~site, ncol = 1, scale = 'fixed') +
  ylab(expression("Giant kelp frond density (no./"*"m"^2*")")) +
  #scale_y_continuous(breaks = seq(0,10, by = 5)) +
  xlab("Year") +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        #axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 10, angle = -45, hjust = 0),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/3) 
#ggsave(filename = "Figures/Giant kelp frond time series - 5 sites.pdf", height = 8, width = 7)

# -------------------------------------------------------------------------------------------------
# Plot time series of understory NPP at all 5 sites
ggplot(data = transect.dat, aes(x = time, y = understory.npp.season.gc.m2.day)) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Control", ], fill = "#1b9e77", alpha = 0.2) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Annual giant kelp removal", ], fill = "#d95f02", alpha = 0.2) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Quarterly giant kelp removal", ], fill = "#7570b3", alpha = 0.2) +
  geom_line(aes(linetype = treatment2, color = treatment2), size = 0.75) +
  scale_linetype_manual(values = c("solid", "solid", "solid"), name = "") +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), name = "") +
  facet_wrap(~site, ncol = 1, scale = 'fixed') +
  ylab(expression(paste("Understory NPP (g C/", "m"^2, "/d)", sep = ""))) +
  scale_y_continuous(breaks = seq(0,10, by = 5)) +
  xlab("Year") +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        #axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 10, angle = -45, hjust = 0),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/3) 
#ggsave(filename = "Figures/Understory NPP time series - 5 sites.pdf", height = 8, width = 7)

# -------------------------------------------------------------------------------------------------
# Plot time series of bottom light at all 5 sites
ggplot(data = transect.dat, aes(x = time, y = mean.bottom.par)) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Control", ], fill = "#1b9e77", alpha = 0.2) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Annual giant kelp removal", ], fill = "#d95f02", alpha = 0.2) +
  geom_area(data = transect.dat[transect.dat$treatment2 == "Quarterly giant kelp removal", ], fill = "#7570b3", alpha = 0.2) +
  geom_line(aes(linetype = treatment2, color = treatment2), size = 0.75) +
  scale_linetype_manual(values = c("solid", "solid", "solid"), name = "") +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), name = "") +
  facet_wrap(~site, ncol = 1, scale = 'fixed') +
  ylab(expression("Daily bottom irradiance (mol/"*"m"^2*"/d)")) +
  #scale_y_continuous(breaks = seq(0,10, by = 5)) +
  xlab("Year") +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        #axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 10, angle = -45, hjust = 0),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/3) 
#ggsave(filename = "Figures/Bottom light time series - 5 sites.pdf", height = 8, width = 7)

ggplot(data = transect.dat, aes(x = treatment4 , y = mean.bottom.par/(1/depth.mllw.m))) +
  stat_summary(fun.data = 'mean_cl_boot', aes(fill = treatment4), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  facet_wrap(~site, ncol = 5, scale = 'fixed') +
  ylab(expression("Daily bottom irradiance (mol/"*"m"^2*"/d) divided by (1 / depth)")) +
  #scale_y_continuous(breaks = seq(0,10, by = 5)) +
  xlab("Treatment") +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = "Figures/Bottom light div by depth - 5 sites.pdf", height = 8, width = 7)

# -------------------------------------------------------------------------------------------------
# Plot mean giant kelp frond density
frond.dens.plot <- ggplot(data = transect.dat, aes(x = treatment4, y = fronds.per.m2, group = paste(site, treatment4))) +
  stat_summary(fun.data = 'mean_cl_boot', aes(fill = treatment4), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  facet_wrap(~site, ncol = 5, scale = 'fixed') +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nTreatment") +
  ylab(expression("Giant kelp frond density (no./"*"m"^2*")")) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(0, 12)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
frond.dens.plot
#ggsave(filename = "Figures/Giant kelp frond means - 5 sites.pdf", height = 4, width = 12)

# -------------------------------------------------------------------------------------------------
# Plot mean giant kelp biomass
biomass.plot <- ggplot(data = transect.dat, aes(x = treatment4, y = giant.kelp.dry.gm2, 
                                                   group = paste(site, treatment4))) +
  stat_summary(fun.data = 'mean_cl_boot', aes(fill = treatment4), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0) +
  facet_wrap(~site, ncol = 5, scale = 'fixed') +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nTreatment") +
  ylab(expression("Giant kelp biomass (g dry/"*"m"^2*")")) +
  scale_y_continuous(breaks = seq(0, 1000, by = 250), expand = expansion(mult = c(0.004, 0.004))) +
  coord_cartesian(ylim = c(0, 1000)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
biomass.plot
#ggsave(filename = "Figures/Giant kelp biomass means - 5 sites.pdf", height = 4, width = 12)

# -------------------------------------------------------------------------------------------------
# Plot mean giant kelp NPP
ggplot(data = transect.dat, aes(x = treatment4, y = giant.kelp.npp.season.gc.m2.day, group = paste(site, treatment4))) +
  stat_summary(fun.data = 'mean_se', aes(fill = treatment4), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0) +
  facet_wrap(~site, ncol = 5, scale = 'fixed') +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nTreatment") +
  ylab(expression(atop(paste("Giant kelp NPP (g C/", "m"^2, "/d)", sep = "")), "")) +
  scale_y_continuous(breaks = seq(0, 8, by = 2), expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(0, 8)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = "Figures/Giant kelp NPP means - 5 sites.pdf", height = 4, width = 12)

# -------------------------------------------------------------------------------------------------
# Plot mean understory NPP
ggplot(data = transect.dat, aes(x = treatment4, y = understory.npp.season.gc.m2.day, group = paste(site, treatment4))) +
  stat_summary(fun.data = 'mean_se', aes(fill = treatment4), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0) +
  facet_wrap(~site, ncol = 5, scale = 'fixed') +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nTreatment") +
  ylab(expression(atop(paste("Understory NPP (g C/", "m"^2, "/d)", sep = "")), "")) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(0, 4)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = "Figures/Understory NPP means - 5 sites.pdf", height = 4, width = 12)

# -------------------------------------------------------------------------------------------------
# Plot mean total macroalgal NPP
ggplot(data = transect.dat, aes(x = treatment4, y = total.npp.season.gc.m2.day, group = paste(site, treatment4))) +
  stat_summary(fun.data = 'mean_se', aes(fill = treatment4), geom = 'bar', color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0) +
  facet_wrap(~site, ncol = 5, scale = 'fixed') +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nTreatment") +
  ylab(expression(atop(paste("Total macroalgal NPP (g C/", "m"^2, "/d)", sep = "")), "")) +
  scale_y_continuous(breaks = seq(0, 8, by = 2), expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(0, 8)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black"),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = "Figures/Total NPP means - 5 sites.pdf", height = 4, width = 12)

# ------------------------------------------------------------------------------------
# Recode season
transect.dat$season2 <- dplyr::recode(transect.dat$season, 
                                               "1-Winter" = "Winter",
                                               "2-Spring" = "Spring", 
                                               "3-Summer" = "Summer", 
                                               "4-Autumn" = "Autumn")
transect.dat$season2 <- factor(transect.dat$season2, levels = c("Winter", "Spring", "Summer", "Autumn"))

transect.dat$season3 <- dplyr::recode(transect.dat$season, 
                                      "2-Spring" = "Sp", 
                                      "3-Summer" = "Su", 
                                      "4-Autumn" = "Au", 
                                      "1-Winter" = "Wi")
transect.dat$season3 <- factor(transect.dat$season3, 
                               levels = c("Wi", "Sp", "Su", "Au"))

# -----------------------------------------------------------------------------------------------------------------
# Plot all NPP together in one plot (averaged across seasons)

transect.dat.plotting <- transect.dat %>%
  dplyr::select(season, season2, season3, site, treatment, treatment2, treatment4,
                giant.kelp.npp.season.gc.m2.day, understory.npp.season.gc.m2.day, total.npp.season.gc.m2.day)
transect.dat.plotting <- rbind(transect.dat.plotting, transect.dat.plotting[nrow(transect.dat.plotting), ])
transect.dat.plotting[nrow(transect.dat.plotting), ]$site <- "Isla Vista"
transect.dat.plotting[nrow(transect.dat.plotting), 8:10] <- 0

transect.dat.plotting <- transect.dat.plotting %>%
  pivot_longer(cols = c('giant.kelp.npp.season.gc.m2.day', 'understory.npp.season.gc.m2.day', 'total.npp.season.gc.m2.day'),
               names_to = "component", 
               values_to = "NPP")

transect.dat.plotting$component <- plyr::mapvalues(x = transect.dat.plotting$component,
                                                   from = unique(transect.dat.plotting$component),
                                                   to = c("(A) Giant kelp NPP", "(B) Understory NPP", "(C) Total macroalgal NPP"))
transect.dat.plotting$component <- factor(transect.dat.plotting$component,
                                          levels = c("(A) Giant kelp NPP", "(B) Understory NPP", "(C) Total macroalgal NPP"))

ggplot(data = transect.dat.plotting, aes(x = site, y = NPP, fill = treatment4)) +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  stat_summary(fun.data = 'mean_se', geom = 'bar', color = 'black', position = position_dodge(width=0.7), width = 0.7) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0,  position = position_dodge(width=0.7)) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nSite") +
  ylab(expression(paste("Net primary productivity (g C/", "m"^2, "/d)", sep = ""))) +
  #scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  scale_y_continuous(breaks = seq(0, 8, by = 2), expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(0, 8)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(1, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/4) 
#ggsave(filename = "Figures/NPP means - 5 sites - unified 2.pdf", height = 8, width = 12)

# -----------------------------------------------------------------------------------------------------------------
# Plot NPP patterns in the final year of the study
transect.dat.final.year <- transect.dat %>%
  dplyr::group_by(plot) %>%
  dplyr::arrange(desc(time)) %>%
  dplyr::slice(1:4) %>%
  dplyr::ungroup()


transect.dat.final.year.plotting <- transect.dat.final.year %>%
  dplyr::select(season, season2, season3, site, treatment, treatment2, treatment4,
                giant.kelp.npp.season.gc.m2.day, understory.npp.season.gc.m2.day, total.npp.season.gc.m2.day)
transect.dat.final.year.plotting <- rbind(transect.dat.final.year.plotting, transect.dat.final.year.plotting[nrow(transect.dat.final.year.plotting), ])
transect.dat.final.year.plotting[nrow(transect.dat.final.year.plotting), ]$site <- "Isla Vista"
transect.dat.final.year.plotting[nrow(transect.dat.final.year.plotting), ]$treatment <- "Continual"
transect.dat.final.year.plotting[nrow(transect.dat.final.year.plotting), ]$treatment2 <- factor("Quarterly giant kelp removal")
transect.dat.final.year.plotting[nrow(transect.dat.final.year.plotting), ]$treatment4 <- factor("Quarterly")
transect.dat.final.year.plotting[nrow(transect.dat.final.year.plotting), 8:10] <- 0

transect.dat.final.year.plotting <- transect.dat.final.year.plotting %>%
  pivot_longer(cols = c('giant.kelp.npp.season.gc.m2.day', 'understory.npp.season.gc.m2.day', 'total.npp.season.gc.m2.day'),
               names_to = "component", 
               values_to = "NPP")

transect.dat.final.year.plotting$component <- plyr::mapvalues(x = transect.dat.final.year.plotting$component,
                                                   from = unique(transect.dat.final.year.plotting$component),
                                                   to = c("(A) Giant kelp NPP", "(B) Understory NPP", "(C) Total macroalgal NPP"))

transect.dat.final.year.plotting$component <- factor(transect.dat.final.year.plotting$component,
                                          levels = c("(A) Giant kelp NPP", "(B) Understory NPP", "(C) Total macroalgal NPP"))

ggplot(data = transect.dat.final.year.plotting, aes(x = site, y = NPP, fill = treatment4)) +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  stat_summary(fun.data = 'mean_se', geom = 'bar', color = 'black', position = position_dodge(width=0.7), width = 0.7) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0,  position = position_dodge(width=0.7)) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("\nSite") +
  ylab(expression(paste("Net primary productivity (g C/", "m"^2, "/d)", sep = ""))) +
  #scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  scale_y_continuous(breaks = seq(0, 14, by = 2), expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(0, 14)) +
  guides(fill = F) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0),
        panel.spacing = unit(1, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        #axis.ticks.x = element_line(size = 0, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/4) 
#ggsave(filename = "Figures/NPP means - 5 sites - final year.pdf", height = 8, width = 12)

# -----------------------------------------------------------------------------------------------------------------
# Plot mean frond dens vs. mean understory NPP using section data
section.dat$treatment2 <- dplyr::recode(section.dat$treatment, 
                                         "Control" = "Control", 
                                         "Annual" = "Annual giant kelp removal", 
                                         "Continual" = "Quarterly giant kelp removal")
section.dat$treatment2 <- factor(section.dat$treatment2, 
                                  levels = c("Control", "Annual giant kelp removal", "Quarterly giant kelp removal"))

compare.dat <- section.dat %>%
  dplyr::group_by(site, treatment2) %>%
  dplyr::summarise(fronds.per.m2_mean = mean(fronds.per.m2_at.section, na.rm = T),
                      fronds.per.m2_sd   = sd(fronds.per.m2_at.section, na.rm = T),
                      fronds.per.m2_n    = length(fronds.per.m2_at.section),
                      understory.npp.season.gc.m2.day_mean = mean(understory.npp.season.gc.m2.day_at.section, na.rm = T),
                      understory.npp.season.gc.m2.day_sd   = sd(understory.npp.season.gc.m2.day_at.section, na.rm = T),
                      understory.npp.season.gc.m2.day_n    = length(understory.npp.season.gc.m2.day_at.section)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fronds.per.m2_se = fronds.per.m2_sd / sqrt(fronds.per.m2_n),
                understory.npp.season.gc.m2.day_se = understory.npp.season.gc.m2.day_sd / sqrt(understory.npp.season.gc.m2.day_n))

compare.dat$site <- factor(compare.dat$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

ggplot(data = compare.dat, aes(x = fronds.per.m2_mean, y = understory.npp.season.gc.m2.day_mean)) +
  geom_errorbar(aes(ymax = understory.npp.season.gc.m2.day_mean + 1.96 * understory.npp.season.gc.m2.day_se,
                ymin = understory.npp.season.gc.m2.day_mean - 1.96 * understory.npp.season.gc.m2.day_se)) +
  geom_errorbarh(aes(xmax = fronds.per.m2_mean + 1.96 * fronds.per.m2_se,
                    xmin = fronds.per.m2_mean - 1.96 * fronds.per.m2_se)) +
  geom_point(aes(shape = site, fill = treatment2), size = 3) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "Site") +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(expression("Giant kelp frond density (no./"*"m"^2*")")) +
  ylab(expression(paste("Understory NPP (g C/", "m"^2, "/d)", sep = ""))) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 3), expand = c(0.025, 0)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1), expand = c(0.025, 0)) +
  theme_classic() +
  theme(text = element_text(size = 16),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
#ggsave(filename = "Figures/All sites - Frond density vs. understory NPP - section data.pdf", height = 4.75, width = 6.75) 

# -----------------------------------------------------------------------------------------------------------------
# Plot mean frond dens vs. mean understory NPP using transect data
compare.dat <- transect.dat %>%
  dplyr::group_by(site, treatment2) %>%
  dplyr::summarise(fronds.per.m2_mean = mean(fronds.per.m2, na.rm = T),
                   fronds.per.m2_sd   = sd(fronds.per.m2, na.rm = T),
                   fronds.per.m2_n    = length(fronds.per.m2),
                   understory.npp.season.gc.m2.day_mean = mean(understory.npp.season.gc.m2.day, na.rm = T),
                   understory.npp.season.gc.m2.day_sd   = sd(understory.npp.season.gc.m2.day, na.rm = T),
                   understory.npp.season.gc.m2.day_n    = length(understory.npp.season.gc.m2.day)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fronds.per.m2_se = fronds.per.m2_sd / sqrt(fronds.per.m2_n),
                understory.npp.season.gc.m2.day_se = understory.npp.season.gc.m2.day_sd / sqrt(understory.npp.season.gc.m2.day_n))

compare.dat$site <- factor(compare.dat$site, levels = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"))

ggplot(data = compare.dat, aes(x = fronds.per.m2_mean, y = understory.npp.season.gc.m2.day_mean)) +
  geom_errorbar(aes(ymax = understory.npp.season.gc.m2.day_mean + 1.96 * understory.npp.season.gc.m2.day_se,
                    ymin = understory.npp.season.gc.m2.day_mean - 1.96 * understory.npp.season.gc.m2.day_se)) +
  geom_errorbarh(aes(xmax = fronds.per.m2_mean + 1.96 * fronds.per.m2_se,
                     xmin = fronds.per.m2_mean - 1.96 * fronds.per.m2_se)) +
  geom_point(aes(shape = site, fill = treatment2), size = 3) +
  #stat_smooth(method = 'lm', color = 'black', se = F) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "Site") +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(expression("Giant kelp frond density (no./"*"m"^2*")")) +
  ylab(expression(paste("Understory NPP (g C/", "m"^2, "/d)", sep = ""))) +
  #coord_cartesian(xlim = c(-0.1, 12), ylim = c(-0.1, 4)) +
  #scale_x_continuous(breaks = seq(0, 12, 3), expand = c(0.025, 0)) +
  #scale_y_continuous(breaks = seq(0, 4, 1), expand = c(0.025, 0)) +
  theme_classic() +
  theme(text = element_text(size = 16),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1)
#ggsave(filename = "Figures/All sites - Frond density vs. understory NPP - transect data.pdf", height = 4.75, width = 6.75) 


# ====================================================================================
# CALCULATE AND COMPARE RATIO OF CANOPY TO UNDERSTORY NPP
# ====================================================================================
transect.dat2 <- transect.dat %>%
  dplyr::mutate(npp.ratio = understory.npp.season.gc.m2.day / giant.kelp.npp.season.gc.m2.day)

ggplot(data = transect.dat2, aes(x = site, y = npp.ratio, fill = treatment4)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  xlab("") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 
rm(transect.dat2)


# ====================================================================================
# CALCULATE DELTA-LIGHT WHEN PAIRED OBSERVATIONS ARE AVAILABLE BETWEEN TREATMENTS
# ====================================================================================

# ------------------------------------------------------------------------------------
# Calculate from daily data
delta.light.dat <- daily.light.kelp.dat %>%
  dplyr::select(-transect, -month, -c(fronds.per.m2:plants.per.m2)) %>%
  dplyr::select(year, time, season, date, site, treatment, light.mol.day) %>%
  tidyr::pivot_wider(id_cols = c(year:site),
                     names_from = treatment,
                     values_from = light.mol.day) %>%
  dplyr::mutate(delta.annual = ANNUAL - CONTROL,
                delta.continual = CONTINUAL - CONTROL) %>%
  dplyr::select(-(ANNUAL:CONTROL)) %>%
  tidyr::pivot_longer(-(year:site),
                      names_to = "treatment",
                      values_to = "delta.light.mol.day") %>%
  dplyr::mutate(treatment = dplyr::recode(treatment,
                                          delta.annual = "Annual",
                                          delta.continual = "Continual")) 

delta.light.dat$treatment <- factor(delta.light.dat$treatment, levels = c("Annual", "Continual"))

delta.light.dat$site <- factor(mapvalues(delta.light.dat$site,
                                         from = c('CARP', 'NAPL', 'AQUE', 'IVEE', 'MOHK'),
                                         to = c("Carpinteria", "Naples", "Arroyo Quemado", "Isla Vista", "Mohawk")),
                             levels = c("Carpinteria", "Naples", "Arroyo Quemado", "Isla Vista", "Mohawk"))

delta.light.dat <- delta.light.dat[!is.na(delta.light.dat$delta.light.mol.day), ]

# ------------------------------------------------------------------------------------
ggplot(data = delta.light.dat, aes(x = time, y = delta.light.mol.day, color = treatment))  +
  facet_grid(treatment ~ site) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  geom_line() +
  ggtitle("Daily data")
#ggsave(filename = "Figures/delta.light.plot-all.dates.pdf", height = 6, width = 18)

# ------------------------------------------------------------------------------------
# Calculate mean delta-light for each year based on a 'kelp year', a  12 month period from spring of 20xx to winter of 20xx + 1

delta.light.dat.annual <- delta.light.dat %>%
  dplyr::mutate(time.temp = time - year,
                kelp.year = NA)

delta.light.dat.annual$kelp.year[delta.light.dat.annual$time.temp == 0.375] <- delta.light.dat.annual$year[delta.light.dat.annual$time.temp == 0.375] + 0.5
delta.light.dat.annual$kelp.year[delta.light.dat.annual$time.temp == 0.625] <- delta.light.dat.annual$year[delta.light.dat.annual$time.temp == 0.625] + 0.5
delta.light.dat.annual$kelp.year[delta.light.dat.annual$time.temp == 0.875] <- delta.light.dat.annual$year[delta.light.dat.annual$time.temp == 0.875] + 0.5
delta.light.dat.annual$kelp.year[delta.light.dat.annual$time.temp == 0.125] <- delta.light.dat.annual$year[delta.light.dat.annual$time.temp == 0.125] + 0.5 - 1

temp1 <- delta.light.dat.annual %>%
  dplyr::filter(treatment == "Annual") %>%
  dplyr::select(year:time, kelp.year, everything(), -time.temp)

temp2 <- delta.light.dat.annual %>%
  dplyr::filter(treatment == "Continual",
                site != "Naples",
                time >= 2010.375 & time < 2017.375) %>%
  dplyr::select(year:time, kelp.year, everything(), -time.temp)

temp3 <- delta.light.dat.annual %>%
  dplyr::filter(treatment == "Continual",
                site == "Naples",
                time >= 2010.375 & time < 2016.375) %>%
  dplyr::select(year:time, kelp.year, everything(), -time.temp)

delta.light.dat.annual <- rbind(temp1, temp2) %>%
  rbind(temp3) %>%
  dplyr::group_by(kelp.year, site, treatment) %>%
  dplyr::summarise(mean_delta.light.mol.day = mean(delta.light.mol.day)) %>%
  ungroup()

rm(temp1, temp2, temp3)

ggplot(data = delta.light.dat.annual, aes(x = kelp.year, y = mean_delta.light.mol.day, color = treatment))  +
  facet_grid(treatment ~ site) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  geom_line() +
  ggtitle("Annual data")
#ggsave(filename = "Figures/delta.light.plot-annual.pdf", height = 6, width = 18)


# ====================================================================================
# QUESTION 1: HOW DOES GIANT KELP BIOMASS AFFECT LIGHT? 
# ====================================================================================

# ------------------------------------------------------------------------------------
# Remove data beyond the time frame of the study
time.max <- max(transect.dat$time)
  
daily.light.kelp.dat2 <- daily.light.kelp.dat %>%
  dplyr::filter(time <= time.max) %>%
  dplyr::filter(treatment == "CONTROL")

# ------------------------------------------------------------------------------------
# Recode season
daily.light.kelp.dat2$season2 <- dplyr::recode(daily.light.kelp.dat2$season, 
                                               "1-Winter" = "Winter",
                                               "2-Spring" = "Spring", 
                                               "3-Summer" = "Summer", 
                                               "4-Autumn" = "Autumn")
daily.light.kelp.dat2$season2 <- factor(daily.light.kelp.dat2$season2, levels = c("Winter", "Spring", "Summer", "Autumn"))

# ------------------------------------------------------------------------------------
# Create statistical models of giant kelp effects on light using glmmTMB

daily.light.kelp.dat2$year2 <- daily.light.kelp.dat2$year - min(daily.light.kelp.dat2$year) + 1

daily.light.kelp.dat2$time.factor <- factor(daily.light.kelp.dat2$time, levels = seq(min(daily.light.kelp.dat2$time),
                                                                                     max(daily.light.kelp.dat2$time),
                                                                                     by = 0.25))

# Remove NA value
daily.light.kelp.dat2 <- daily.light.kelp.dat2 %>%
  dplyr::filter(!is.na(fronds.per.m2))

# ------------------------------------------------------------------------------------
# Convert from frond density to canopy biomass using monthly parameter estimates from Rassweiler et al.

daily.light.kelp.dat2$parm1 <- NA

# Note: These parameters convert frond density (no. / m2) to biomass density (kg dry / m2)
# convert to g dry / m2 by multiplying by 1000
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 1]  <- 0.07154
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 2]  <- 0.07578
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 3]  <- 0.0872
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 4]  <- 0.09642
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 5]  <- 0.09922
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 6]  <- 0.09355
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 7]  <- 0.08648
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 8]  <- 0.08359
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 9]  <- 0.07447
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 10] <- 0.07200
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 11] <- 0.06032
daily.light.kelp.dat2$parm1[daily.light.kelp.dat2$month == 12] <- 0.06594

daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2 <- daily.light.kelp.dat2$fronds.per.m2 * daily.light.kelp.dat2$parm1
daily.light.kelp.dat2$giant.kelp.biomass.dry.gm2  <- daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2 * 1000

ggplot(data = daily.light.kelp.dat2, aes(x = fronds.per.m2, y = giant.kelp.biomass.dry.kgm2)) +
  stat_smooth(aes(group = month), formula = y ~ 0 + x, method = 'lm', se = F) +
  geom_point() +
  facet_wrap(~ season2)

ggplot(data = daily.light.kelp.dat2, aes(x = giant.kelp.biomass.dry.kgm2, y = light.mol.day)) +
  geom_point() +
  facet_wrap(~ season2)

# ------------------------------------------------------------------------------------
# Create statistical models of giant kelp canopy biomass effects on light using glmmTMB

par.mod1   <- glmmTMB(light.mol.day ~ giant.kelp.biomass.dry.kgm2 * season + site + year2 +
                         ar1(time.factor + 0 | site),
                       family = Gamma("log"),
                       dispformula = ~ site + season,
                       data = daily.light.kelp.dat2)

par.mod2   <- glmmTMB(light.mol.day ~ giant.kelp.biomass.dry.kgm2 + season + site + year2 +
                       ar1(time.factor + 0 | site),
                     family = Gamma("log"),
                     dispformula = ~ site + season,
                     data = daily.light.kelp.dat2)

anova(par.mod1, par.mod2)
AICtab(par.mod1, par.mod2)

par.mod <- par.mod1

emtrends(par.mod, pairwise ~ season, var = "giant.kelp.biomass.dry.kgm2")$emtrends

summary(par.mod)
car::Anova(par.mod)

# Pseudo-R2 (squared correlation between the response and the predicted value; after Bolker)
cor(daily.light.kelp.dat2$light.mol.day, predict(par.mod, type="response"))^2

# ------------------------------------------------------------------------------------
# Validate par.mod - Normality, homogeneity of variance, dispersion, and zero inflation
simulationOutput <- simulateResiduals(par.mod, n=1000)
plot(simulationOutput)
testDispersion(simulationOutput) # Test for over/underdispersion based on simulated residuals
testZeroInflation(simulationOutput)

# Validate par.mod - Check for missing predictors
daily.light.kelp.dat2$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(daily.light.kelp.dat2$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(daily.light.kelp.dat2$resid, xlab = "Pearson residuals", main = "")
plot(predict(par.mod), daily.light.kelp.dat2$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(daily.light.kelp.dat2$time, daily.light.kelp.dat2$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(daily.light.kelp.dat2$year), daily.light.kelp.dat2$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(daily.light.kelp.dat2$season), daily.light.kelp.dat2$resid, xlab = "Season", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(daily.light.kelp.dat2$site), daily.light.kelp.dat2$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(daily.light.kelp.dat2$treatment), daily.light.kelp.dat2$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2, daily.light.kelp.dat2$resid, xlab = "Giant kelp canopy biomass", ylab = "Pearson residuals"); abline(c(0,0))

# ------------------------------------------------------------------------------------
# Validate par.mod - Check for temporal autocorrelation among par.model residuals
daily.light.kelp.dat2$resid <- resid(par.mod)

# Calculate ACF and PACF for each plot (transect)
resid.dat.transect <- daily.light.kelp.dat2 %>%
  dplyr::select(site, transect, time, resid) %>%
  dplyr::mutate(transect = paste(site, transect, sep = "_")) %>%
  dplyr::select(-site)
resid.dat.transect <- na.omit(resid.dat.transect)

acf.dat <- sapply(unique(resid.dat.transect$transect), function(x){
  acf(resid.dat.transect$resid[resid.dat.transect$transect == x], lag.max = length(unique(daily.light.kelp.dat2$time)) / 3, plot = FALSE)$acf
}
)

pacf.dat <- sapply(unique(resid.dat.transect$transect), function(x){
  pacf(resid.dat.transect$resid[resid.dat.transect$transect == x], lag.max = length(unique(daily.light.kelp.dat2$time)) / 3, plot = FALSE)$acf
}
)

acf.dat <- data.frame(acf.dat)
pacf.dat <- data.frame(pacf.dat)

# Use this code for multiple sites
acf.dat <- data.frame(t(plyr::ldply(acf.dat, rbind)))[-1, ]
pacf.dat <- data.frame(t(plyr::ldply(pacf.dat, rbind)))[-1, ]

colnames(acf.dat) <- unique(resid.dat.transect$transect)
colnames(pacf.dat) <- unique(resid.dat.transect$transect)

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "transect", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "transect", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "transect"))
rm(pacf.dat)

# Calculate critical value (based on the lowest length of time series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(daily.light.kelp.dat2[daily.light.kelp.dat2$transect == 4, ]$time)))

acf.dat$crit[acf.dat$plot == "Mohawk_Continual"] <- qnorm((1 + 0.95)/2) / sqrt(length(unique(daily.light.kelp.dat2[daily.light.kelp.dat2$transect == 4, ]$time)))

# Fix formatting issue
acf.dat$acf <- as.numeric(acf.dat$acf)
acf.dat$pacf <- as.numeric(acf.dat$pacf)

# Plot ACF by plot
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  facet_wrap(~ transect) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag (quarters)") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ylim(-1, 1)
p1

# Plot PACF by plot
p2 <- ggplot(data = acf.dat, aes(x = lag, y = pacf)) +
  ggtitle("Partial autocorrelation by plot") +
  facet_wrap(~ transect) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag (quarters)") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ylim(-1, 1)
p2

# For multiple sites: plot mean ACF
p3 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag (quarters)") +
  theme_classic() +
  theme(aspect.ratio = 1)

p4 <- ggplot(data = acf.dat, aes(x = lag, y = pacf)) +
  ggtitle("Partial autocorrelation by plot") +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag (quarters)") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p3 | p4)

# ------------------------------------------------------------------------------------
# Calculate effect sizes
em1 <- emtrends(par.mod, ~ season, var = "giant.kelp.biomass.dry.kgm2", transform = "response")
plot(em1)
summary(em1, infer=c(TRUE,TRUE),null=0)

em2 <- emtrends(par.mod, ~ 1, var = "giant.kelp.biomass.dry.kgm2", transform = "response")
plot(em2)
summary(em2, infer=c(TRUE,TRUE),null=0)

# ------------------------------------------------------------------------------------
# # Save results
# sink(here::here("Results/Results_Bottom-PAR-vs-kelp-biomass_All-sites.txt"))
# summary(par.mod)
# car::Anova(par.mod)
# 
# ref_grid(par.mod)
# emtrends(par.mod, ~ season, var = "giant.kelp.biomass.dry.kgm2", transform = "response")
# 
# sink()
# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals - based on mean giant kelp biomass observed in each treatment
preds.df <- daily.light.kelp.dat2 %>%
  dplyr::select(site, year2, season, time.factor) %>%
  dplyr::distinct()

n <- 3
kelp.breaks <- c(mean(transect.dat$giant.kelp.dry.kgm2[transect.dat$treatment == "Control"]), 
                 mean(transect.dat$giant.kelp.dry.kgm2[transect.dat$treatment == "Annual"]),
                 mean(transect.dat$giant.kelp.dry.kgm2[transect.dat$treatment == "Continual"]))

# Replicate kelp.breaks a many times as rows in pred.df ... then replicate preds.df 'n' times
kelp.vec <- sort(rep(kelp.breaks, nrow(preds.df)))
preds.df <- do.call("rbind", replicate(n, preds.df, simplify = FALSE))

# Add kelp data
preds.df$giant.kelp.biomass.dry.kgm2 <- kelp.vec
rm(kelp.vec)

# Predict new data
preds.df$light.mol.day   <- predict(par.mod, newdata = preds.df, type = "response", se.fit = FALSE)
preds.df$se              <- predict(par.mod, newdata = preds.df, type = "response", se.fit = TRUE)$se.fit
preds.df$upr.ci          <- preds.df$light.mol.day + (1.96 * preds.df$se)
preds.df$lwr.ci          <- preds.df$light.mol.day - (1.96 * preds.df$se)

# Summarize predictions as a function of giant kelp frond density
preds.avg.by.treatment <- preds.df %>%
  dplyr::group_by(giant.kelp.biomass.dry.kgm2) %>%
  dplyr::summarise_at(vars(light.mol.day:lwr.ci),
                      mean,
                      na.rm = T) %>%
  dplyr::ungroup()

preds.avg.by.treatment$treatment = factor(c("Quarterly\ndisturbance", "Annual\ndisturbance", "Control\n"),
                             levels = c("Control\n", "Annual\ndisturbance", "Quarterly\ndisturbance"))

preds.avg.by.treatment <- dplyr::select(preds.avg.by.treatment, treatment, everything())

ggplot(data = preds.avg.by.treatment, aes(x = treatment, y = light.mol.day)) +
  geom_bar(aes(fill = treatment), stat="identity", color = 'black') +
  #geom_errorbar(aes(ymax = upr.ci, ymin = lwr.ci), width = 0) +
  ylab(expression("Daily bottom irradiance (mol/"*"m"^2*"/d)")) +
  xlab("Treatment") +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        aspect.ratio = 1)

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals

preds.df <- daily.light.kelp.dat2 %>%
  dplyr::select(site, year2, season, time.factor) %>%
  dplyr::distinct()

n <- 50
min.val <- min(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2[!is.na(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2)])
max.val <- max(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2[!is.na(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2)])
by.val <- (max.val - min.val) / (n - 1)
kelp.breaks <-  seq(min.val, max.val, by = by.val)

# Replicate kelp.breaks a many times as rows in pred.df ... then replicate preds.df 'n' times
kelp.vec <- sort(rep(kelp.breaks, nrow(preds.df)))
preds.df <- do.call("rbind", replicate(n, preds.df, simplify = FALSE))

# Add kelp data
preds.df$giant.kelp.biomass.dry.kgm2 <- kelp.vec
rm(kelp.vec)

# Predict new data
preds.df$light.mol.day   <- predict(par.mod, newdata = preds.df, type = "response", se.fit = FALSE)
preds.df$se              <- predict(par.mod, newdata = preds.df, type = "response", se.fit = TRUE)$se.fit
preds.df$upr.ci          <- preds.df$light.mol.day + (1.96 * preds.df$se)
preds.df$lwr.ci          <- preds.df$light.mol.day - (1.96 * preds.df$se)

# Summarize predictions as a function of giant kelp frond density
preds.avg <- preds.df %>%
  dplyr::group_by(giant.kelp.biomass.dry.kgm2) %>%
  dplyr::summarise_at(vars(light.mol.day:lwr.ci),
                      mean,
                      na.rm = T) %>%
  dplyr::ungroup()

# Summarize predictions as a function of giant kelp biomass and season
preds.avg.season <- preds.df %>%
  dplyr::group_by(giant.kelp.biomass.dry.kgm2, season) %>%
  dplyr::summarise_at(vars(light.mol.day:lwr.ci),
                      mean,
                      na.rm = T) %>%
  dplyr::ungroup()

# For the non-averaged data, truncate predictions based on observed range of data for each unique combination of factors
min.max.dat <- daily.light.kelp.dat2[!is.na(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2), ] %>%
  dplyr::group_by(season, site) %>%
  dplyr::summarise(min.gk = min(giant.kelp.biomass.dry.kgm2),
                   max.gk = max(giant.kelp.biomass.dry.kgm2)) %>%
  ungroup()

preds.df <- dplyr::left_join(preds.df, min.max.dat) %>%
  dplyr::filter(giant.kelp.biomass.dry.kgm2 >= min.gk) %>%
  dplyr::filter(giant.kelp.biomass.dry.kgm2 <= max.gk) %>%
  dplyr::select( -min.gk, -max.gk)

# ------------------------------------------------------------------------------------
# For seasonal predictions, constrain predictions to maximum value of giant kelp biomass observed in each season
max1 <- max(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2[daily.light.kelp.dat2$season == "1-Winter"])
max2 <- max(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2[daily.light.kelp.dat2$season == "2-Spring"])
max3 <- max(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2[daily.light.kelp.dat2$season == "3-Summer"])
max4 <- max(daily.light.kelp.dat2$giant.kelp.biomass.dry.kgm2[daily.light.kelp.dat2$season == "4-Autumn"])

preds.avg.season$filter <- 1

preds.avg.season$filter[preds.avg.season$season == "1-Winter" &
                          preds.avg.season$giant.kelp.biomass.dry.kgm2 > max1] <- NA
preds.avg.season$filter[preds.avg.season$season == "2-Spring" &
                          preds.avg.season$giant.kelp.biomass.dry.kgm2 > max2] <- NA
preds.avg.season$filter[preds.avg.season$season == "3-Summer" &
                          preds.avg.season$giant.kelp.biomass.dry.kgm2 > max3] <- NA
preds.avg.season$filter[preds.avg.season$season == "4-Autumn" &
                          preds.avg.season$giant.kelp.biomass.dry.kgm2 > max4] <- NA

preds.avg.season <- preds.avg.season[!is.na(preds.avg.season$filter), ] %>%
  dplyr::select(-filter)

# ------------------------------------------------------------------------------------
# Plot predictions
light.plot <- ggplot(data = daily.light.kelp.dat2, aes(x = giant.kelp.biomass.dry.kgm2, y = light.mol.day)) +
  geom_ribbon(data = preds.avg, aes(y = light.mol.day, ymin = lwr.ci, ymax = upr.ci), fill = "dimgrey", alpha = 0.33) +
  geom_point(size = 1, shape = 21, fill = "black") +
  geom_line(data = preds.avg, aes(y = light.mol.day), color = 'black', size = 1.5) +
  ylab(expression("Daily bottom irradiance (mol/"*"m"^2*"/d)")) +
  xlab(expression("Giant kelp biomass (kg dry/"*"m"^2*")")) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        aspect.ratio = 1)
light.plot
#ggsave(file = here::here("Figures/giant.kelp.biomass.vs.light.mol.day.pdf"), width = 5.5, height = 4.5, units = "in")
#ggsave(file = here::here("Figures/giant.kelp.biomass.vs.light.mol.day.svg"), width = 5.5, height = 4.5, units = "in")

light.plot + 
  geom_point(data = preds.avg.by.treatment, aes(fill = treatment, shape = treatment)) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(22, 21, 24), name = "Treatment") 
#ggsave(file = here::here("Figures/giant.kelp.biomass.vs.light.mol.day_with_experiment_preds.pdf"), width = 5.5, height = 4.5, units = "in")

# ------------------------------------------------------------------------------------
# Recode season and re-order levels for plotting
daily.light.kelp.dat2$season2 <- dplyr::recode(daily.light.kelp.dat2$season,
                                               "2-Spring" = "Spring",
                                               "3-Summer" = "Summer",
                                               "4-Autumn" = "Autumn",
                                               "1-Winter" = "Winter")
daily.light.kelp.dat2$season2 <- factor(daily.light.kelp.dat2$season2, levels = c("Winter", "Spring", "Summer", "Autumn"))

preds.avg.season$season2 <- dplyr::recode(preds.avg.season$season,
                                          "2-Spring" = "Spring",
                                          "3-Summer" = "Summer",
                                          "4-Autumn" = "Autumn",
                                          "1-Winter" = "Winter")
preds.avg.season$season2 <- factor(preds.avg.season$season2, levels = c("Winter", "Spring", "Summer", "Autumn"))

daily.light.kelp.dat2$season4 <- dplyr::recode(daily.light.kelp.dat2$season,
                                               "2-Spring" = "(F) Spring",
                                               "3-Summer" = "(G) Summer",
                                               "4-Autumn" = "(H) Autumn",
                                               "1-Winter" = "(E) Winter")
daily.light.kelp.dat2$season4 <- factor(daily.light.kelp.dat2$season4, levels = c("(E) Winter", "(F) Spring", "(G) Summer", "(H) Autumn"))

preds.avg.season$season4 <- dplyr::recode(preds.avg.season$season,
                                          "2-Spring" = "(F) Spring",
                                          "3-Summer" = "(G) Summer",
                                          "4-Autumn" = "(H) Autumn",
                                          "1-Winter" = "(E) Winter")
preds.avg.season$season4 <- factor(preds.avg.season$season4, levels = c("(E) Winter", "(F) Spring", "(G) Summer", "(H) Autumn"))

ggplot(data = daily.light.kelp.dat2, aes(x = giant.kelp.biomass.dry.kgm2, y = light.mol.day)) +
  geom_ribbon(data = preds.avg.season, aes(y = light.mol.day, ymin = lwr.ci, ymax = upr.ci, fill = season2), alpha = 0.7) +
  geom_point(size = 1.5, shape = 21, aes(fill = season2)) +
  geom_line(data = preds.avg.season, aes(y = light.mol.day), size = 1.25) +
  scale_fill_manual(values = c("#80b1d3", "#fccde5", "#b3de69", "#e6ab02")) +
  scale_color_manual(values = c("#80b1d3", "#fccde5", "#b3de69", "#e6ab02")) +
  facet_wrap(~season2, ncol = 4) +
  ylab(expression("Daily bottom irradiance (mol/"*"m"^2*"/d)")) +
  xlab(expression("Giant kelp biomass (kg dry/"*"m"^2*")")) +
  #coord_cartesian(ylim = c(0, 4)) +
  #scale_x_continuous(breaks = seq(0, 1.2, by = 0.4)) +
  #scale_y_continuous(breaks = seq(0, 4, by = 1)) +
  theme_classic() +
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=15),
        panel.spacing = unit(0.2, "lines"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) +
  guides(fill = F, color = F)
#ggsave(file = here::here("Figures/giant.kelp.biomass.vs.light.mol.day.by.season.pdf"), width = 9, height = 3.5, units = "in")


# =====================================================================================================================
# CREATE TRANSFORMATIONS
# =====================================================================================================================

transect.dat$giant.kelp.dry.kgm2 <- transect.dat$giant.kelp.dry.gm2 / 1000
section.dat$giant.kelp.dry.kgm2_at.section <- section.dat$giant.kelp.dry.gm2_at.section / 1000

# Data wrangling
section.dat$season2 <- dplyr::recode(section.dat$season, 
                                     "2-Spring" = "Spring", 
                                     "3-Summer" = "Summer", 
                                     "4-Autumn" = "Autumn", 
                                     "1-Winter" = "Winter")
section.dat$season2 <- factor(section.dat$season2, 
                              levels = c("Winter", "Spring", "Summer", "Autumn"))

section.dat$season3 <- dplyr::recode(section.dat$season, 
                                     "2-Spring" = "Sp", 
                                     "3-Summer" = "Su", 
                                     "4-Autumn" = "Au", 
                                     "1-Winter" = "Wi")
section.dat$season3 <- factor(section.dat$season3, 
                              levels = c("Wi", "Sp", "Su", "Au"))

# Transform variables to satisfy normality and homoskedasticity requirements
section.dat$giant.kelp.npp.season.gc.m2.day_at.section_sq.root <- section.dat$giant.kelp.npp.season.gc.m2.day_at.section^(1/2)
section.dat$giant.kelp.npp.season.gc.m2.day_at.section_cube.root <- section.dat$giant.kelp.npp.season.gc.m2.day_at.section^(1/3)

section.dat$understory.npp.season.gc.m2.day_at.section_sq.root <- section.dat$understory.npp.season.gc.m2.day_at.section ^ (1/2)
section.dat$understory.npp.season.gc.m2.day_at.section_cube.root <- section.dat$understory.npp.season.gc.m2.day_at.section ^ (1/3)

section.dat$total.npp.season.gc.m2.day_at.section_sq.root <- section.dat$total.npp.season.gc.m2.day_at.section ^ (1/2)
section.dat$total.npp.season.gc.m2.day_at.section_cube.root <- section.dat$total.npp.season.gc.m2.day_at.section ^ (1/3)

section.dat$urchin.dry.gm2_at.section_cube.root <- section.dat$urchin.dry.gm2_at.section^(1/3)


# =================================================================================================================
# PLOT RAW DATA FOR ANALYSES TO COME FURTHER DOWN IN SCRIPT
# =================================================================================================================

# ---------------------------------------------------------------------------------------------------
ggplot(data = transect.dat, aes(x = urchin.thresh, y = giant.kelp.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  scale_x_discrete(name = "Sea urchin density") +
  ylab(expression("Giant kelp NPP (g C/"*"m"^2*"/d)")) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  #facet_wrap(~ site, ncol = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 

ggplot(data = transect.dat, aes(x = sand.thresh, y = giant.kelp.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  scale_x_discrete(name = "Sand cover") +
  ylab(expression("Giant kelp NPP (g C/"*"m"^2*"/d)")) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  #facet_wrap(~ site, ncol = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 

# ---------------------------------------------------------------------------------------------------
ggplot(data = transect.dat, aes(x = urchin.thresh, y = understory.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  scale_x_discrete(name = "Sea urchin density") +
  ylab(expression("Understory NPP (g C/"*"m"^2*"/d)")) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  #facet_wrap(~ site, ncol = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 

ggplot(data = transect.dat, aes(x = sand.thresh, y = understory.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  scale_x_discrete(name = "Sand cover") +
  ylab(expression("Understory NPP (g C/"*"m"^2*"/d)")) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  #facet_wrap(~ site, ncol = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 

# ---------------------------------------------------------------------------------------------------
ggplot(data = transect.dat, aes(x = urchin.thresh, y = total.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  scale_x_discrete(name = "Sea urchin density") +
  ylab(expression("Total NPP (g C/"*"m"^2*"/d)")) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  #facet_wrap(~ site, ncol = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 

ggplot(data = transect.dat, aes(x = sand.thresh, y = total.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) +
  scale_fill_manual(values = c(col1, col2, col3), name = "") +
  scale_x_discrete(name = "Sand cover") +
  ylab(expression("Total NPP (g C/"*"m"^2*"/d)")) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  #facet_wrap(~ site, ncol = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 

# ---------------------------------------------------------------------------------------------------
p1 <- ggplot(data = transect.dat.annual, 
                   aes(x = season2, y = giant.kelp.dry.kgm2, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar", position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col2), name = "") +
  xlab("") +
  ylab(expression("Giant kelp biomass (kg dry/"*"m"^2*")")) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.5), expand = expansion(mult = c(0.005, 0.006))) +
  coord_cartesian(ylim = c(0, 1.5)) +
  facet_wrap(~ site, ncol = 1) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 
p1
#ggsave(file = here::here("Figures/giant.kelp.biomass.by.season_all.sites_annual.v.control.pdf"), width = 4, height = 10, units = "in")

p1 +
  facet_wrap(~ site, ncol = 5)
#ggsave(file = here::here("Figures/giant.kelp.biomass.by.season_all.sites_annual.v.control-wide.pdf"), width = 15, height = 3.5, units = "in")

p2 <- ggplot(data = transect.dat.continual, 
                   aes(x = season2, y = giant.kelp.dry.kgm2, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col3), name = "") +
  xlab("") +
  ylab(expression("Giant kelp biomass (g dry/"*"m"^2*")")) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.5), expand = expansion(mult = c(0.005, 0.006))) +
  coord_cartesian(ylim = c(0, 1.5)) +
  facet_wrap(~ site, ncol = 1) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) 
p2
#ggsave(file = here::here("Figures/giant.kelp.biomass.by.season_all.sites_continual.v.control.pdf"), width = 4, height = 10, units = "in")

p2 +
  facet_wrap(~ site, ncol = 5)
#ggsave(file = here::here("Figures/giant.kelp.biomass.by.season_all.sites_continual.v.control-wide.pdf"), width = 15, height = 3.5, units = "in")

# --------------------------------------------------------------------------------------------------------------------
npp.p1a <- ggplot(data = transect.dat.annual, aes(x = season2, y = giant.kelp.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col2), name = "Treatment") +
  xlab("") +
  ylab(expression(paste("Giant kelp NPP (g C/", "m"^2, "/d)", sep = ""))) +
  coord_cartesian(ylim = c(-0.03, 12.03)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), expand = expansion(mult = c(0.004, 0.004))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) +
  facet_wrap(~ site, ncol = 5) 
npp.p1a
#ggsave(file = here::here("Figures/giant.kelp.NPP_control.vs.annual-all-sites.pdf"), width = 15, height = 3.5, units = "in")

npp.p1b <- ggplot(data = transect.dat.annual, aes(x = season2, y = understory.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col2), name = "Giant kelp treatment") +
  xlab("") +
  ylab(expression(paste("Understory NPP (g C/", "m"^2, "/d)", sep = ""))) +
  coord_cartesian(ylim = c(-0.03, 6.75)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = expansion(mult = c(0.004, 0.004))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) +
  facet_wrap(~ site, ncol = 5) +
  guides(fill = F)
npp.p1b
#ggsave(file = here::here("Figures/understory.NPP_control.vs.annual-all-sites.pdf"), width = 15, height = 3.5, units = "in")

npp.p1c <- ggplot(data = transect.dat.annual, aes(x = season2, y = total.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col2), name = "Giant kelp treatment") +
  xlab("") +
  ylab(expression(paste("Total macroalgal NPP (g C/", "m"^2, "/d)", sep = ""))) +
  coord_cartesian(ylim = c(-0.03, 12.03)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), expand = expansion(mult = c(0.004, 0.004))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) +
  facet_wrap(~ site, ncol = 5)  +
  guides(fill = F)
npp.p1c
#ggsave(file = here::here("Figures/total.NPP_control.vs.annual-all-sites.pdf"), width = 15, height = 3.5, units = "in")

# --------------------------------------------------------------------------------------------------------------------
npp.p2a <- ggplot(data = transect.dat.continual, aes(x = season2, y = giant.kelp.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col3), name = "Treatment") +
  xlab("") +
  ylab(expression(paste("Giant kelp NPP (g C/", "m"^2, "/d)", sep = ""))) +
  coord_cartesian(ylim = c(-0.03, 12.03)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), expand = expansion(mult = c(0.004, 0.004))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) +
  facet_wrap(~ site, ncol = 5) 
npp.p2a
#ggsave(file = here::here("Figures/giant.kelp.NPP_control.vs.continual-all-sites.pdf"), width = 15, height = 3.5, units = "in")

npp.p2b <- ggplot(data = transect.dat.continual, aes(x = season2, y = understory.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col3), name = "Giant kelp treatment") +
  xlab("") +
  ylab(expression(paste("Understory NPP (g C/", "m"^2, "/d)", sep = ""))) +
  coord_cartesian(ylim = c(-0.03, 6.75)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = expansion(mult = c(0.004, 0.004))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) +
  facet_wrap(~ site, ncol = 5) +
  guides(fill = F)
npp.p2b
#ggsave(file = here::here("Figures/understory.NPP_control.vs.continual-all-sites.pdf"), width = 15, height = 3.5, units = "in")

npp.p2c <- ggplot(data = transect.dat.continual, aes(x = season2, y = total.npp.season.gc.m2.day, fill = treatment2)) +
  stat_summary(fun = "mean", geom = "bar",position=position_dodge(0.65), color = "black", width = 0.65) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position=position_dodge(0.65), width = 0) + 
  scale_fill_manual(values = c(col1, col3), name = "Giant kelp treatment") +
  xlab("") +
  ylab(expression(paste("Total macroalgal NPP (g C/", "m"^2, "/d)", sep = ""))) +
  coord_cartesian(ylim = c(-0.03, 12.03)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), expand = expansion(mult = c(0.004, 0.004))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=12),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size =  unit(1, "lines"), 
        legend.spacing.x = unit(0.11, 'cm'),
        legend.spacing.y = unit(0.11, 'cm'),
        legend.background = element_blank(),
        legend.position = "top",
        aspect.ratio = 1/2) +
  facet_wrap(~ site, ncol = 5)  +
  guides(fill = F)
npp.p2c
#ggsave(file = here::here("Figures/total.NPP_control.vs.continual-all-sites.pdf"), width = 15, height = 3.5, units = "in")


# =============================================================================================================
# QUESTION 2: HOW DOES GIANT KELP BIOMASS DENSITY AFFECT UNDERSTORY NPP?
# =============================================================================================================

# ------------------------------------------------------------------------------------
# Create new treatment factor
transect.dat.annual$treatment2 <- factor(plyr::mapvalues(transect.dat.annual$treatment, 
                                                  from = c("Control", "Annual", "Continual"),
                                                  to = c("Control", "Annual disturbance", "Quarterly disturbance")),
                                         levels = c("Control", "Annual disturbance", "Quarterly disturbance"))

# ------------------------------------------------------------------------------------
# Create kelp.year as factor
transect.dat.annual$kelp.year.factor <- factor(transect.dat.annual$kelp.year, 
                                               levels = seq(min(transect.dat.annual$kelp.year),
                                                                   max(transect.dat.annual$kelp.year),
                                                                   by = 1))
transect.dat.annual$kelp.year2 <- transect.dat.annual$kelp.year - min(transect.dat.annual$kelp.year) + 1
section.dat$year2 <- section.dat$year - min(section.dat$year) + 1

# ------------------------------------------------------------------------------------
transect.dat.annual$mean.giant.kelp.dry.kgm2 <- transect.dat.annual$mean.giant.kelp.dry.gm2 / 1000 
transect.dat.annual$mean.urchin.dry.kgm2 <- transect.dat.annual$mean.urchin.dry.gm2 / 1000

# ------------------------------------------------------------------------------------
# Create statistical models

hist(transect.dat.annual$understory.npp_kgC.m2.y)
hist(sqrt(transect.dat.annual$understory.npp_kgC.m2.y))
hist(log(transect.dat.annual$understory.npp_kgC.m2.y))
range(transect.dat.annual$understory.npp_kgC.m2.y)

transect.dat.annual$sqrt.understory.npp_kgC.m2.y <- sqrt(transect.dat.annual$understory.npp_kgC.m2.y)

understory.npp.vs.giant.kelp.mod   <- glmmTMB(understory.npp_kgC.m2.y ~
                                                mean.giant.kelp.dry.kgm2 +
                                                mean.urchin.dry.kgm2 +
                                                mean.percent.sand +
                                                kelp.year2 +
                                                ar1(kelp.year.factor + 0 | site),
                                              family = Gamma(link = "log"),
                                              data = transect.dat.annual)

summary(understory.npp.vs.giant.kelp.mod)
car::Anova(understory.npp.vs.giant.kelp.mod)

check_collinearity(understory.npp.vs.giant.kelp.mod)

# ------------------------------------------------------------------------------------
# Validate understory.npp.vs.giant.kelp.mod - Normality, homogeneity of variance, dispersion, and zero inflation
simulationOutput <- simulateResiduals(understory.npp.vs.giant.kelp.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput) # Test for over/underdispersion based on simulated residuals
testZeroInflation(simulationOutput)

# Validate understory.npp.vs.giant.kelp.mod - Check for missing predictors
transect.dat.annual$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(transect.dat.annual$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(transect.dat.annual$resid, xlab = "Pearson residuals", main = "")
plot(predict(understory.npp.vs.giant.kelp.mod), transect.dat.annual$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$kelp.year, transect.dat.annual$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(transect.dat.annual$site), transect.dat.annual$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(transect.dat.annual$treatment), transect.dat.annual$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$mean.giant.kelp.dry.kgm2, transect.dat.annual$resid, xlab = "Giant kelp canopy biomass density", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$mean.urchin.dry.kgm2, transect.dat.annual$resid, xlab = "Sea urchin density", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$mean.percent.sand, transect.dat.annual$resid, xlab = "Sand cover (%)", ylab = "Pearson residuals"); abline(c(0,0))

# ------------------------------------------------------------------------------------
# Validate understory.npp.vs.giant.kelp.mod - Check for temporal autocorrelation among residuals
transect.dat.annual$resid <- resid(understory.npp.vs.giant.kelp.mod)

resid.dat.transect <- transect.dat.annual %>%
  dplyr::select(site, treatment, kelp.year, resid) %>%
  dplyr::mutate(transect = paste(site, treatment, sep = "_")) %>%
  dplyr::select(-site)
resid.dat.transect <- na.omit(resid.dat.transect)

acf.dat <- sapply(unique(resid.dat.transect$transect), function(x){
  acf(resid.dat.transect$resid[resid.dat.transect$transect == x], lag.max = length(unique(transect.dat.annual$kelp.year)) / 3, plot = FALSE)$acf
}
)

pacf.dat <- sapply(unique(resid.dat.transect$transect), function(x){
  pacf(resid.dat.transect$resid[resid.dat.transect$transect == x], lag.max = length(unique(transect.dat.annual$kelp.year)) / 3, plot = FALSE)$acf
}
)

acf.dat <- data.frame(acf.dat)
pacf.dat <- data.frame(pacf.dat)

# Use this code for multiple sites
acf.dat <- data.frame(t(plyr::ldply(acf.dat, rbind)))[-1, ]
pacf.dat <- data.frame(t(plyr::ldply(pacf.dat, rbind)))[-1, ]

colnames(acf.dat) <- unique(resid.dat.transect$transect)
colnames(pacf.dat) <- unique(resid.dat.transect$transect)

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "transect", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "transect", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "transect"))
rm(pacf.dat)

# Calculate critical value (based on the lowest length of kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(7)

# Fix formatting issue
acf.dat$acf <- as.numeric(acf.dat$acf)
acf.dat$pacf <- as.numeric(acf.dat$pacf)

# For multiple sites: plot mean ACF
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag (quarters)") +
  theme_classic() +
  theme(aspect.ratio = 1)

p2 <- ggplot(data = acf.dat, aes(x = lag, y = pacf)) +
  ggtitle("Partial autocorrelation by plot") +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 1:max(acf.dat$lag), name = "Lag (quarters)", limits = c(1, max(acf.dat$lag))) +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# # Save results
# sink(here::here("Results/Results_understory-NPP-vs-giant-kelp-biomass_All-sites.txt"))
# summary(understory.npp.vs.giant.kelp.mod)
# car::Anova(understory.npp.vs.giant.kelp.mod)
# 
# ref_grid(understory.npp.vs.giant.kelp.mod)
# em <- emtrends(understory.npp.vs.giant.kelp.mod, var = "mean.giant.kelp.dry.kgm2", transform = "response", infer = c(T,T), adjust = "none")
# summary(em)
# 
# sink()

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals
breaks.kelp <-  seq(min(transect.dat.annual$mean.giant.kelp.dry.kgm2), 
                    max(transect.dat.annual$mean.giant.kelp.dry.kgm2), 
                    by = (max(transect.dat.annual$mean.giant.kelp.dry.kgm2) - 
                            min(transect.dat.annual$mean.giant.kelp.dry.kgm2))/100)

preds <- ggemmeans(understory.npp.vs.giant.kelp.mod, 
                   terms = c("mean.giant.kelp.dry.kgm2[breaks.kelp]"),
                   type = "fixed", #)
                   back.transform = T)

preds <- data.frame(preds) %>%
  dplyr::select(-std.error, -group) %>%
  dplyr::rename(understory.npp_kgC.m2.y = predicted,
                mean.giant.kelp.dry.kgm2 = x)

ggplot(data = preds, aes(x = mean.giant.kelp.dry.kgm2, y = understory.npp_kgC.m2.y)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_line(size = 1)

# ------------------------------------------------------------------------------------
# Plot predictions - giant kelp
ggplot(data = transect.dat.annual, aes(x = mean.giant.kelp.dry.kgm2, y = understory.npp_kgC.m2.y)) +
  annotate("rect", 
           xmin = -0.02,
           xmax = 0.195,
           ymin = 1.5,
           ymax = 2.42,
           fill = "#FFE3EB",
           alpha = 0.6,
           color = "#FFE3EB") +
  geom_ribbon(data = preds, aes(y = understory.npp_kgC.m2.y, ymin = conf.low, ymax = conf.high), 
              fill = "dimgrey", alpha = 0.33) +
  geom_point(aes(shape = treatment2, fill = treatment2), size = 2, color = 'black') +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(22, 21, 24), name = "Treatment") +
  geom_line(data = preds, aes(y = understory.npp_kgC.m2.y), color = 'black', size = 1.5) +
  ylab(expression("Annual NPP (kg C/"*"m"^2*"/y)")) +
  xlab(expression("Giant kelp biomass (kg dry/"*"m"^2*")")) +
  ggtitle("(B) Understory NPP") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 2.5)) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  coord_cartesian(xlim = c(0, 1.23)) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        legend.text = element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 11, hjust = 0, face = 'bold'),
        legend.spacing.y = unit(0.15, 'lines'),
        legend.background = element_blank(),
        aspect.ratio = 1)
#ggsave(file = here::here("Figures/giant.kelp.vs.understory.NPP.pdf"), width = 7, height = 4.5, units = "in")
#ggsave(file = here::here("Figures/giant.kelp.vs.understory.NPP.svg"), width = 7, height = 4.5, units = "in")

# Calculate effect size of decline
round(100 * (0.07538344 - 0.55411252) / 0.55411252)

# Calculate percentile
percentile.dat <- sort(transect.dat.annual$understory.npp_kgC.m2.y)
k <- 0.95
n <- length(percentile.dat)
index.val <- round(k * n)
percentile.val <- round(percentile.dat[index.val + 1], 2)
percentile.val

# Find corresponding maximum giant kelp biomass values
round(max(transect.dat.annual$mean.giant.kelp.dry.kgm2[transect.dat.annual$understory.npp_kgC.m2.y >= percentile.val]), 2)

# Compare percentile to control treatment max
max(transect.dat.annual$understory.npp_kgC.m2.y[transect.dat.annual$treatment == "Control"])


# ===========================================================================================================================
# QUESTION 3: HOW DOES GIANT KELP BIOMASS DENSITY AFFECT TOTAL NPP?
# ===========================================================================================================================

# ------------------------------------------------------------------------------------
# Create statistical models

hist(transect.dat.annual$total.npp_kgC.m2.y)
hist(sqrt(transect.dat.annual$total.npp_kgC.m2.y))
hist(log(transect.dat.annual$total.npp_kgC.m2.y))
range(transect.dat.annual$total.npp_kgC.m2.y)

total.npp.vs.giant.kelp.mod   <- glmmTMB(total.npp_kgC.m2.y ~
                                           mean.giant.kelp.dry.kgm2 +
                                           mean.urchin.dry.kgm2 +
                                           mean.percent.sand +
                                           kelp.year2 +
                                           ar1(kelp.year.factor + 0 | site),
                                         dispformula = ~ site + treatment,
                                         family = gaussian(link = "identity"),
                                         data = transect.dat.annual)

summary(total.npp.vs.giant.kelp.mod)
car::Anova(total.npp.vs.giant.kelp.mod)

check_collinearity(total.npp.vs.giant.kelp.mod)

# ------------------------------------------------------------------------------------
# Validate total.npp.vs.giant.kelp.mod - Normality, homogeneity of variance, dispersion, and zero inflation
simulationOutput <- simulateResiduals(total.npp.vs.giant.kelp.mod, n=1000)
plot(simulationOutput, quantiles = NA)
testDispersion(simulationOutput) # Test for over/underdispersion based on simulated residuals
testZeroInflation(simulationOutput)

# Validate total.npp.vs.giant.kelp.mod - Check for missing predictors
transect.dat.annual$resid <- simulationOutput$scaledResiduals
#transect.dat.annual$resid <- resid(total.npp.vs.giant.kelp.mod, type = "pearson")

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(transect.dat.annual$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(transect.dat.annual$resid, xlab = "Pearson residuals", main = "")
plot(predict(total.npp.vs.giant.kelp.mod), transect.dat.annual$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$kelp.year, transect.dat.annual$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(transect.dat.annual$site), transect.dat.annual$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(transect.dat.annual$treatment), transect.dat.annual$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$mean.giant.kelp.dry.kgm2, transect.dat.annual$resid, xlab = "Giant kelp canopy biomass density", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$mean.urchin.dry.kgm2, transect.dat.annual$resid, xlab = "Sea urchin density", ylab = "Pearson residuals"); abline(c(0,0))
plot(transect.dat.annual$mean.percent.sand, transect.dat.annual$resid, xlab = "Sand cover (%)", ylab = "Pearson residuals"); abline(c(0,0))

# ------------------------------------------------------------------------------------
# Validate total.npp.vs.giant.kelp.mod - Check for temporal autocorrelation among residuals
transect.dat.annual$resid <- resid(total.npp.vs.giant.kelp.mod)

resid.dat.transect <- transect.dat.annual %>%
  dplyr::select(site, treatment, kelp.year, resid) %>%
  dplyr::mutate(transect = paste(site, treatment, sep = "_")) %>%
  dplyr::select(-site)
resid.dat.transect <- na.omit(resid.dat.transect)

acf.dat <- sapply(unique(resid.dat.transect$transect), function(x){
  acf(resid.dat.transect$resid[resid.dat.transect$transect == x], lag.max = length(unique(transect.dat.annual$kelp.year)) / 3, plot = FALSE)$acf
}
)

pacf.dat <- sapply(unique(resid.dat.transect$transect), function(x){
  pacf(resid.dat.transect$resid[resid.dat.transect$transect == x], lag.max = length(unique(transect.dat.annual$kelp.year)) / 3, plot = FALSE)$acf
}
)

acf.dat <- data.frame(acf.dat)
pacf.dat <- data.frame(pacf.dat)

# Use this code for multiple sites
acf.dat <- data.frame(t(plyr::ldply(acf.dat, rbind)))[-1, ]
pacf.dat <- data.frame(t(plyr::ldply(pacf.dat, rbind)))[-1, ]

colnames(acf.dat) <- unique(resid.dat.transect$transect)
colnames(pacf.dat) <- unique(resid.dat.transect$transect)

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "transect", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "transect", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "transect"))
rm(pacf.dat)

# Calculate critical value (based on the lowest length of kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(7)

# Fix formatting issue
acf.dat$acf <- as.numeric(acf.dat$acf)
acf.dat$pacf <- as.numeric(acf.dat$pacf)

# For multiple sites: plot mean ACF
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag (quarters)") +
  theme_classic() +
  theme(aspect.ratio = 1)

p2 <- ggplot(data = acf.dat, aes(x = lag, y = pacf)) +
  ggtitle("Partial autocorrelation by plot") +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 1:max(acf.dat$lag), name = "Lag (quarters)", limits = c(1, max(acf.dat$lag))) +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# # Save results
# sink(here::here("Results/Results_total-NPP-vs-giant-kelp-biomass_All-sites.txt"))
# summary(total.npp.vs.giant.kelp.mod)
# car::Anova(total.npp.vs.giant.kelp.mod)
# 
# ref_grid(total.npp.vs.giant.kelp.mod)
# em <- emtrends(total.npp.vs.giant.kelp.mod, var = "mean.giant.kelp.dry.kgm2", transform = "response", infer = c(T,T), adjust = "none")
# summary(em)
# 
# sink()

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals
preds <- ggemmeans(total.npp.vs.giant.kelp.mod, 
                   terms = c("mean.giant.kelp.dry.kgm2[breaks.kelp]"),
                   type = "fixed")

preds <- data.frame(preds) %>%
  dplyr::select(-std.error, -group) %>%
  dplyr::rename(total.npp_kgC.m2.y = predicted,
                mean.giant.kelp.dry.kgm2 = x) ## %>%
  # 
  # # Backtransform
  # dplyr::mutate(conf.low = conf.low ^ 2,
  #               conf.high = conf.high ^ 2,
  #               total.npp_kgC.m2.y = total.npp_kgC.m2.y ^ 2)

ggplot(data = preds, aes(x = mean.giant.kelp.dry.kgm2, y = total.npp_kgC.m2.y)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_line(size = 1)

# ------------------------------------------------------------------------------------
# Plot predictions - giant kelp
p <- ggplot(data = transect.dat.annual, aes(x = mean.giant.kelp.dry.kgm2, y = total.npp_kgC.m2.y)) +
  annotate("rect", 
           xmin = -0.02,
           xmax = 0.195,
           ymin = 1.87,
           ymax = 2.705,
           fill = "#FFE3EB",
           alpha = 0.6,
           color = "#FFE3EB") +
  # annotate("rect", 
  #          xmin = 0.484,
  #          xmax = 0.888,
  #          ymin = 1.87,
  #          ymax = 2.705,
  #          fill = "#C6EEE4",
  #          alpha = 0.6,
  #          color = "#C6EEE4") +
  geom_ribbon(data = preds, aes(y = total.npp_kgC.m2.y, ymin = conf.low, ymax = conf.high), 
              fill = "dimgrey", alpha = 0.5) +
  geom_point(aes(shape = treatment2, fill = treatment2), size = 2, color = 'black') +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_shape_manual(values = c(22, 21, 24), name = "Treatment") +
  geom_line(data = preds, aes(y = total.npp_kgC.m2.y), color = 'black', size = 1.5) +
  ylab(expression("Annual NPP (kg C/"*"m"^2*"/y)")) +
  xlab(expression("Giant kelp biomass (kg dry/"*"m"^2*")")) +
  ggtitle("(A) Total macroalgal NPP") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0, 1.23)) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        legend.text = element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 11, hjust = 0, face = 'bold'),
        legend.spacing.y = unit(0.15, 'lines'),
        legend.background = element_blank(),
        aspect.ratio = 1) 
p
#ggsave(file = here::here("Figures/giant.kelp.vs.total.NPP.pdf"), width = 7, height = 4.5, units = "in")
#ggsave(file = here::here("Figures/giant.kelp.vs.total.NPP.svg"), width = 7, height = 4.5, units = "in")

# Calculate effect size of increase
round(3.5325091/0.5285989, 1)
3.5/0.5

# Identify anomalously high values
high.vals <- transect.dat.annual %>%
  dplyr::filter(mean.giant.kelp.dry.kgm2 < 0.2) %>%
  dplyr::filter(treatment == "Annual" | treatment == "Continual") %>%
  dplyr::select(mean.giant.kelp.dry.kgm2, total.npp_kgC.m2.y) %>%
  dplyr::filter(total.npp_kgC.m2.y > 1.9) %>%
  dplyr::arrange(desc(total.npp_kgC.m2.y))

high.vals
p + geom_vline(xintercept = 0.2) + geom_hline(yintercept = c(1.94, 2.63)) + 
  geom_vline(xintercept = c(0.505, 0.865)) 

0.865/0.00849
0.865/0.167

# Find annual NPP at median giant kelp biomass value in control plots
mean(transect.dat.annual$mean.giant.kelp.dry.kgm2[transect.dat.annual$treatment == "Control"])
median(transect.dat.annual$mean.giant.kelp.dry.kgm2[transect.dat.annual$treatment == "Control"])


# =============================================================================================================================
# QUESTION 4: HOW DOES ANNUAL UNDERSTORY NPP CHANGE BASED ON TREATMENT * TIME-SINCE-START * HABITAT-QUALITY 
# =============================================================================================================================

# -----------------------------------------------------------------------------------------------------------------------------
# Add habitat quality category (high, medium, low)
delta.npp.dat.annual$habitat.quality.cat <- factor(plyr::mapvalues(delta.npp.dat.annual$site,
                                                   from = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"),
                                                   to = c("High", "High", "Medium", "Medium", "Low")),
                                   levels = c("Low", "Medium", "High"))

delta.npp.dat.annual <- delta.npp.dat.annual %>%
  dplyr::select(kelp.year, site, habitat.quality.cat, treatment, everything())

# -----------------------------------------------------------------------------------------------------------------------------
# Calculate time since start
delta.npp.dat.annual <- delta.npp.dat.annual %>%
  dplyr::mutate(plot = paste(site, treatment, sep = "_")) %>%
  dplyr::group_by(plot) %>%
  dplyr::mutate(min.time = min(kelp.year)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(time.since.start = kelp.year - min.time + 1,
                time.factor = factor(kelp.year)) %>%
  dplyr::select(-min.time, kelp.year, time.factor, time.since.start, site, treatment, plot, habitat.quality.cat, everything())

# -----------------------------------------------------------------------------------------------------------------------------
# Split annual delta NPP data into macroalgal components
delta.npp.dat.annual.understory <- delta.npp.dat.annual %>%
  dplyr::filter(npp.component == "understory") %>%
  dplyr::select(-npp.component)

delta.npp.dat.annual.total <- delta.npp.dat.annual %>%
  dplyr::filter(npp.component == "total") %>%
  dplyr::select(-npp.component)

delta.npp.dat.annual.giant.kelp <- delta.npp.dat.annual %>%
  dplyr::filter(npp.component == "giant.kelp") %>%
  dplyr::select(-npp.component)

# -----------------------------------------------------------------------------------------------------------------------------
# Fit models - Delta NPP data
hist(delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y)
hist((delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y + 1) ^ (1/2))
hist((delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y + 1) ^ (1/3))
hist(log(delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y + 1))

x <- 4.5
y <- log(x + 1)
exp(y) == x + 1
exp(y) - 1 == x

delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y_plus.one <- delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y + 1
delta.npp.dat.annual.understory$log.annual.delta.npp_kgC.m2.y <- log(delta.npp.dat.annual.understory$annual.delta.npp_kgC.m2.y_plus.one)
# 
# understory.delta.npp.mod <-  glmmTMB(log.annual.delta.npp_kgC.m2.y ~
#                   #treatment * time.since.start * habitat.quality.cat, #+
#                   treatment * poly(time.since.start, 2) * habitat.quality.cat, #+
#                   #ar1(time.factor + 0 | plot),
#                 family = gaussian(link = "identity"),
#                 data = delta.npp.dat.annual.understory)

understory.delta.npp.mod <-  glmmTMB(annual.delta.npp_kgC.m2.y_plus.one ~
                                     treatment * time.since.start * habitat.quality.cat, #+
                                     #treatment * poly(time.since.start, 2) * habitat.quality.cat, #+
                                     #ar1(time.factor + 0 | plot),
                                     family = Gamma(link = "log"),
                                     data = delta.npp.dat.annual.understory)

AIC(understory.delta.npp.mod)
car::Anova(understory.delta.npp.mod)
summary(understory.delta.npp.mod)

# ---
# Validate understory.delta.npp.mod - Check for missing predictors
simulationOutput <- simulateResiduals(understory.delta.npp.mod, n=1000)
plot(simulationOutput)
testDispersion(simulationOutput) # Test for over/underdispersion based on simulated residuals
testZeroInflation(simulationOutput)

# Validate understory.delta.npp.mod - Check for missing predictors
delta.npp.dat.annual.understory$resid <- resid(understory.delta.npp.mod)
delta.npp.dat.annual.understory$resid <- NA
delta.npp.dat.annual.understory$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(delta.npp.dat.annual.understory$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(delta.npp.dat.annual.understory$resid, xlab = "Pearson residuals", main = "")
plot(predict(understory.delta.npp.mod), delta.npp.dat.annual.understory$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.understory$habitat.quality.cat, delta.npp.dat.annual.understory$resid, xlab = "Habitat quality", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(delta.npp.dat.annual.understory$treatment), delta.npp.dat.annual.understory$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(delta.npp.dat.annual.understory$kelp.year), delta.npp.dat.annual.understory$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.understory$time.since.start, delta.npp.dat.annual.understory$resid, xlab = "Time since start", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.understory$site, delta.npp.dat.annual.understory$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))

# ---
# Validate understory.delta.npp.mod - Check for temporal autocorrelation among understory.delta.npp.model residuals
delta.npp.dat.annual.understory$resid <- resid(understory.delta.npp.mod)

# ---
# Calculate ACF and PACF for each plot
acf.dat <- sapply(unique(delta.npp.dat.annual.understory$plot), function(x){
  acf(delta.npp.dat.annual.understory$resid[delta.npp.dat.annual.understory$plot == x],
      lag.max = 4, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(delta.npp.dat.annual.understory$plot), function(x){
  pacf(delta.npp.dat.annual.understory$resid[delta.npp.dat.annual.understory$plot == x],
       lag.max = 4, plot = FALSE)$acf
}
)

acf.dat <- data.frame(acf.dat)
pacf.dat <- data.frame(pacf.dat)

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "plot", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "plot", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "plot"))

# Calculate critical value (based on the lowest length of kelp.kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(delta.npp.dat.annual.understory[delta.npp.dat.annual.understory$plot == "Isla Vista_Annual", ]$kelp.year)))

# Plot ACF by plot
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  facet_wrap(~ plot) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# Plot average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across plots") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals
breaks.hab <- c("Low", "Medium", "High")

breaks.time <-  seq(min(delta.npp.dat.annual.understory$time.since.start), 
                    max(delta.npp.dat.annual.understory$time.since.start), 
                    by = (max(delta.npp.dat.annual.understory$time.since.start)-min(delta.npp.dat.annual.understory$time.since.start))/100)

preds <- ggemmeans(understory.delta.npp.mod, 
                   terms = c("time.since.start[breaks.time]",
                             "treatment",
                             "habitat.quality.cat[breaks.hab]" #"habitat.quality[breaks.hab]", 
                             ),
                   type = "fixed", #)
                   back.transform = T)

preds <- data.frame(preds) %>%
  dplyr::select(-std.error) %>%
  dplyr::mutate(predicted = predicted - 1,
                conf.low = conf.low - 1,
                conf.high = conf.high - 1) %>%
  # dplyr::mutate(predicted = exp(predicted) - 1,
  #               conf.low = exp(conf.low) - 1,
  #               conf.high = exp(conf.high) - 1) %>%
  dplyr::rename(annual.delta.npp_kgC.m2.y = predicted,
                habitat.quality.cat = facet,
                treatment = group,
                time.since.start = x)

preds$time.since.start <- as.numeric(as.character(preds$time.since.start))

preds2 <- preds

preds2$annual.delta.npp_kgC.m2.y[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
preds2$conf.low[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
preds2$conf.high[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA

preds2$habitat.quality.cat <- factor(preds2$habitat.quality.cat, levels = c("Low", "Medium", "High"))
delta.npp.dat.annual.understory$habitat.quality.cat <- factor(delta.npp.dat.annual.understory$habitat.quality.cat, 
                                                              levels = c("Low", "Medium", "High"))

ggplot(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  #facet_grid(treatment ~ habitat.quality.cat) +  
  facet_wrap(~habitat.quality.cat) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = treatment, fill = treatment), alpha = 0.25) +
  geom_line(aes(colour = treatment, group = treatment), size = 1)

# # Compare trend lines to zero
em <- emtrends(understory.delta.npp.mod, ~ habitat.quality.cat | treatment, var = "time.since.start", max.degree = 1)
summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

# Remove trend lines that are not significant over time vs. zero 
preds2$annual.delta.npp_kgC.m2.y[preds$treatment == "Annual" & preds$habitat.quality.cat == "Medium"] <- NA
preds2$conf.low[preds$treatment == "Annual" & preds$habitat.quality.cat == "Medium"] <- NA
preds2$conf.high[preds$treatment == "Annual" & preds$habitat.quality.cat == "Medium"] <- NA

preds2$annual.delta.npp_kgC.m2.y[preds$habitat.quality.cat == "Low"] <- NA
preds2$conf.low[preds$habitat.quality.cat == "Low"] <- NA
preds2$conf.high[preds$habitat.quality.cat == "Low"] <- NA

# Plot
delta.npp.dat.annual.understory$treatment <- factor(delta.npp.dat.annual.understory$treatment, levels = c("Continual", "Annual"))
preds2$treatment <- factor(preds2$treatment, levels = c("Continual", "Annual"))

delta.npp.dat.annual.understory$treatment2 <- factor(mapvalues(delta.npp.dat.annual.understory$treatment, 
                                      from = c("Continual", "Annual"),
                                      to = c("Quarterly disturbance", "Annual disturbance")),
                                      levels = c("Annual disturbance", "Quarterly disturbance"))
preds2$treatment2 <- factor(mapvalues(preds2$treatment, 
                               from = c("Continual", "Annual"),
                               to = c("Quarterly disturbance", "Annual disturbance")),
                            levels = c("Annual disturbance", "Quarterly disturbance"))

col1.line <- "#2CA282"
col2.line <- "#2F85B7"
col3.line <- "#E88B42"

ggplot(data = delta.npp.dat.annual.understory, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_grid(treatment2 ~ habitat.quality.cat) +  
  geom_ribbon(data = preds2, 
              aes(ymin = conf.low, ymax = conf.high, group = treatment2, fill = treatment2), alpha = 0.33) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
  geom_line(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y,
                               color = treatment2, group = treatment), size = 1) +
  scale_shape_manual(values = c(21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col2.line, col3.line), name = "Treatment") +
  ylab(expression(atop("Understory annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  xlab("Years since start of experiment") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  scale_y_continuous(breaks = seq(-1, 3, by = 0.5)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-0.55, 2)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1,
        legend.position = "top", #c(0.81, 0.95),
        legend.text = element_text(size = 11, hjust = 0.5),
        legend.title = element_text(size = 11, hjust = 0.5, face = 'bold'),
        legend.spacing.y = unit(0.15, 'lines'),
        legend.background = element_blank()) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
#ggsave(filename = "Figures/new.understory.delta-npp.pdf", width = 8, height = 6)
ggsave(filename = "Figures/new.understory.delta-npp.svg", width = 8, height = 6)


# =============================================================================================================================
# QUESTION 5: HOW DOES ANNUAL TOTAL NPP CHANGE BASED ON TREATMENT * TIME-SINCE-START * HABITAT-QUALITY 
# =============================================================================================================================

# -----------------------------------------------------------------------------------------------------------------------------
# Fit models - Delta NPP data
hist(delta.npp.dat.annual.total$annual.delta.npp_kgC.m2.y)
hist((delta.npp.dat.annual.total$annual.delta.npp_kgC.m2.y + 3) ^ (1/2))
hist((delta.npp.dat.annual.total$annual.delta.npp_kgC.m2.y + 3) ^ (1/3))
hist(log(delta.npp.dat.annual.total$annual.delta.npp_kgC.m2.y + 1))

delta.npp.dat.annual.total$annual.delta.npp_kgC.m2.y_plus.three <- delta.npp.dat.annual.total$annual.delta.npp_kgC.m2.y + 3

total.delta.npp.mod <-  glmmTMB(annual.delta.npp_kgC.m2.y ~
                                #  treatment + time.since.start + habitat.quality.cat, #+
                                treatment * time.since.start * habitat.quality.cat, #+
                                #treatment + time.since.start + habitat.quality.cat, #+
                                #treatment * poly(time.since.start, 2) * habitat.quality.cat, #+
                                #ar1(time.factor + 0 | plot),
                                family = gaussian(link = "identity"),
                                dispformula = ~ habitat.quality.cat,
                                data = delta.npp.dat.annual.total)

car::Anova(total.delta.npp.mod)
summary(total.delta.npp.mod)

# ---
# Validate total.delta.npp.mod - Check for missing predictors
simulationOutput <- simulateResiduals(total.delta.npp.mod, n=1000)
plot(simulationOutput)
testDispersion(simulationOutput) # Test for over/underdispersion based on simulated residuals
testZeroInflation(simulationOutput)

# Validate total.delta.npp.mod - Check for missing predictors
delta.npp.dat.annual.total$resid <- resid(total.delta.npp.mod)
delta.npp.dat.annual.total$resid <- NA
delta.npp.dat.annual.total$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(delta.npp.dat.annual.total$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(delta.npp.dat.annual.total$resid, xlab = "Pearson residuals", main = "")
plot(predict(total.delta.npp.mod), delta.npp.dat.annual.total$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.total$habitat.quality.cat, delta.npp.dat.annual.total$resid, xlab = "Habitat quality", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(delta.npp.dat.annual.total$treatment), delta.npp.dat.annual.total$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(delta.npp.dat.annual.total$kelp.year), delta.npp.dat.annual.total$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.total$time.since.start, delta.npp.dat.annual.total$resid, xlab = "Time since start", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.total$site, delta.npp.dat.annual.total$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))

# ---
# Validate total.delta.npp.mod - Check for temporal autocorrelation among total.delta.npp.model residuals
delta.npp.dat.annual.total$resid <- resid(total.delta.npp.mod)

# ---
# Calculate ACF and PACF for each plot
acf.dat <- sapply(unique(delta.npp.dat.annual.total$plot), function(x){
  acf(delta.npp.dat.annual.total$resid[delta.npp.dat.annual.total$plot == x],
      lag.max = 4, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(delta.npp.dat.annual.total$plot), function(x){
  pacf(delta.npp.dat.annual.total$resid[delta.npp.dat.annual.total$plot == x],
       lag.max = 4, plot = FALSE)$acf
}
)

acf.dat <- data.frame(acf.dat)
pacf.dat <- data.frame(pacf.dat)

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "plot", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "plot", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "plot"))

# Calculate critical value (based on the lowest length of kelp.kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(delta.npp.dat.annual.total[delta.npp.dat.annual.total$plot == "Isla Vista_Annual", ]$kelp.year)))

# Plot ACF by plot
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  facet_wrap(~ plot) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# Plot average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across plots") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals
# preds <- ggemmeans(total.delta.npp.mod, 
#                    terms = c("time.since.start[breaks.time]",
#                              "treatment",
#                              "habitat.quality.cat[breaks.hab]" #"habitat.quality[breaks.hab]", 
#                    ),
#                    type = "fixed")
# 
# preds <- data.frame(preds) %>%
#   dplyr::select(-std.error) %>%
#   dplyr::mutate(predicted = predicted,
#                 conf.low = conf.low,
#                 conf.high = conf.high) %>%
#   dplyr::rename(annual.delta.npp_kgC.m2.y = predicted,
#                 habitat.quality.cat = facet,
#                 treatment = group,
#                 time.since.start = x)
# 
# preds$time.since.start <- as.numeric(as.character(preds$time.since.start))
# 
# preds2 <- preds
# 
# preds2$annual.delta.npp_kgC.m2.y[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
# preds2$conf.low[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
# preds2$conf.high[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
# 
# preds2$habitat.quality.cat <- factor(preds2$habitat.quality.cat, levels = c("Low", "Medium", "High"))
# delta.npp.dat.annual.total$habitat.quality.cat <- factor(delta.npp.dat.annual.total$habitat.quality.cat, 
#                                                               levels = c("Low", "Medium", "High"))
# 
# ggplot(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
#   #facet_grid(treatment ~ habitat.quality.cat) +  
#   facet_wrap(~habitat.quality.cat) +
#   geom_hline(yintercept = 0, linetype = 'dashed') +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = treatment, fill = treatment), alpha = 0.25) +
#   geom_line(aes(colour = treatment, group = treatment), size = 1)

# Compare means to each other and calculate effect sizes
em <- emmeans(total.delta.npp.mod, pairwise ~ treatment)
summary(em, infer=c(TRUE,TRUE), adjust = 'none')

# # Compare trend lines to zero
# em <- emtrends(total.delta.npp.mod, ~ habitat.quality.cat | treatment, var = "time.since.start", max.degree = 1)
# summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

# # Remove trend lines that are not significant over time vs. zero 
# preds2$annual.delta.npp_kgC.m2.y[preds$habitat.quality.cat == "Low"] <- NA
# preds2$conf.low[preds$habitat.quality.cat == "Low"] <- NA
# preds2$conf.high[preds$habitat.quality.cat == "Low"] <- NA
# 
# preds2$annual.delta.npp_kgC.m2.y[preds$habitat.quality.cat == "High"] <- NA
# preds2$conf.low[preds$habitat.quality.cat == "High"] <- NA
# preds2$conf.high[preds$habitat.quality.cat == "High"] <- NA
# 
# preds2$annual.delta.npp_kgC.m2.y[preds$treatment == "Continual"] <- NA
# preds2$conf.low[preds$treatment == "Continual"] <- NA
# preds2$conf.high[preds$treatment == "Continual"] <- NA

# Plot
delta.npp.dat.annual.total$treatment <- factor(delta.npp.dat.annual.total$treatment, levels = c("Continual", "Annual"))
delta.npp.dat.annual.total$treatment2 <- factor(mapvalues(delta.npp.dat.annual.total$treatment, 
                                                               from = c("Continual", "Annual"),
                                                               to = c("Quarterly disturbance", "Annual disturbance")),
                                                     levels = c("Annual disturbance", "Quarterly disturbance"))

#preds2$treatment <- factor(preds2$treatment, levels = c("Continual", "Annual"))
# preds2$treatment2 <- factor(mapvalues(preds2$treatment, 
#                                       from = c("Continual", "Annual"),
#                                       to = c("Quarterly disturbance", "Annual disturbance")),
#                             levels = c("Annual disturbance", "Quarterly disturbance"))

ggplot(data = delta.npp.dat.annual.total, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_grid(treatment2 ~ habitat.quality.cat) +  
#  geom_ribbon(data = preds2, 
#              aes(ymin = conf.low, ymax = conf.high, group = treatment2, fill = treatment2), alpha = 0.33) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
#  geom_line(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y,
#                               color = treatment2, group = treatment), size = 1) +
  scale_shape_manual(values = c(21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col2.line, col3.line), name = "Treatment") +
  ylab(expression(atop("Total macroalgal annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  xlab("Years since start of experiment") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  scale_y_continuous(breaks = seq(-3, 1, by = 1)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-3, 1)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1,
        legend.position = "top", #c(0.81, 0.95),
        legend.text = element_text(size = 11, hjust = 0.5),
        legend.title = element_text(size = 11, hjust = 0.5, face = 'bold'),
        legend.spacing.y = unit(0.15, 'lines'),
        legend.background = element_blank()) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
#ggsave(filename = "Figures/new.total.delta-npp--scatter.pdf", width = 8, height = 6)

# ------------------------------------------------------------------------------------
# Plot means
delta.npp.dat.annual.total$treatment2 <- factor(delta.npp.dat.annual.total$treatment2,
                                                levels = c("Annual disturbance", "Quarterly disturbance"))

delta.npp.dat.annual.total$treatment3 <- factor(mapvalues(delta.npp.dat.annual.total$treatment2,
                                                          from = c("Annual disturbance", "Quarterly disturbance"),
                                                          to = c("Annual", "Quarterly")),
                                                levels = c("Annual", "Quarterly"))

# Check significance of bars from zero
em <- emmeans(total.delta.npp.mod, ~ treatment | habitat.quality.cat, infer = T, adjust = 'none')
summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

ggplot(data = delta.npp.dat.annual.total, aes(x = habitat.quality.cat, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  #ylab(expression(paste("Change in total NPP from control (kg C/", "m"^2, "/y)", sep = ""))) +
  ylab(expression(atop("Total macroalgal annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  xlab("Site habitat quality") +
  scale_y_continuous(breaks = seq(-3, 1, by = 0.5)) + 
  coord_cartesian(ylim = c(-2, 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 11, hjust = 0, face = 'bold'),
        legend.spacing.y = unit(0.2, 'lines'),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        aspect.ratio = 1/2) 
#ggsave(filename = "Figures/new.total.delta-npp--bar.pdf", width = 6, height = 4)
#ggsave(filename = "Figures/new.total.delta-npp--bar.svg", width = 6, height = 4)


# =============================================================================================================================
# QUESTION 6: HOW DOES ANNUAL GIANT KELP NPP CHANGE BASED ON TREATMENT * TIME-SINCE-START * HABITAT-QUALITY 
# =============================================================================================================================

# -----------------------------------------------------------------------------------------------------------------------------
# Fit models - Delta NPP data
hist(delta.npp.dat.annual.giant.kelp$annual.delta.npp_kgC.m2.y)
hist((delta.npp.dat.annual.giant.kelp$annual.delta.npp_kgC.m2.y + 3) ^ (1/2))
hist((delta.npp.dat.annual.giant.kelp$annual.delta.npp_kgC.m2.y + 3) ^ (1/3))
hist(log(delta.npp.dat.annual.giant.kelp$annual.delta.npp_kgC.m2.y + 1))

delta.npp.dat.annual.giant.kelp$annual.delta.npp_kgC.m2.y_plus.three <- delta.npp.dat.annual.giant.kelp$annual.delta.npp_kgC.m2.y + 3.2

giant.kelp.delta.npp.mod <-  glmmTMB(annual.delta.npp_kgC.m2.y ~
                                  #  treatment + time.since.start + habitat.quality.cat, #+
                                  treatment * time.since.start * habitat.quality.cat, #+
                                #treatment + time.since.start + habitat.quality.cat, #+
                                #treatment * poly(time.since.start, 2) * habitat.quality.cat, #+
                                #ar1(time.factor + 0 | plot),
                                family = gaussian(link = "identity"),
                                dispformula = ~ habitat.quality.cat,
                                data = delta.npp.dat.annual.giant.kelp)

car::Anova(giant.kelp.delta.npp.mod)
summary(giant.kelp.delta.npp.mod)

emmeans(giant.kelp.delta.npp.mod, pairwise ~ treatment)
emmeans(giant.kelp.delta.npp.mod, pairwise ~ treatment, "habitat.quality.cat")

# ---
# Validate giant.kelp.delta.npp.mod - Check for missing predictors
simulationOutput <- simulateResiduals(giant.kelp.delta.npp.mod, n=1000)
plot(simulationOutput)
testDispersion(simulationOutput) # Test for over/underdispersion based on simulated residuals
testZeroInflation(simulationOutput)

# Validate giant.kelp.delta.npp.mod - Check for missing predictors
delta.npp.dat.annual.giant.kelp$resid <- resid(giant.kelp.delta.npp.mod)
delta.npp.dat.annual.giant.kelp$resid <- NA
delta.npp.dat.annual.giant.kelp$resid <- simulationOutput$scaledResiduals

par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))
qqPlot(delta.npp.dat.annual.giant.kelp$resid, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
hist(delta.npp.dat.annual.giant.kelp$resid, xlab = "Pearson residuals", main = "")
plot(predict(giant.kelp.delta.npp.mod), delta.npp.dat.annual.giant.kelp$resid, xlab = "Predicted values", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.giant.kelp$habitat.quality.cat, delta.npp.dat.annual.giant.kelp$resid, xlab = "Habitat quality", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(delta.npp.dat.annual.giant.kelp$treatment), delta.npp.dat.annual.giant.kelp$resid, xlab = "Treatment", ylab = "Pearson residuals"); abline(c(0,0))
plot(as.factor(delta.npp.dat.annual.giant.kelp$kelp.year), delta.npp.dat.annual.giant.kelp$resid, xlab = "Year", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.giant.kelp$time.since.start, delta.npp.dat.annual.giant.kelp$resid, xlab = "Time since start", ylab = "Pearson residuals"); abline(c(0,0))
plot(delta.npp.dat.annual.giant.kelp$site, delta.npp.dat.annual.giant.kelp$resid, xlab = "Site", ylab = "Pearson residuals"); abline(c(0,0))

# ---
# Validate giant.kelp.delta.npp.mod - Check for temporal autocorrelation among giant.kelp.delta.npp.model residuals
delta.npp.dat.annual.giant.kelp$resid <- resid(giant.kelp.delta.npp.mod)

# ---
# Calculate ACF and PACF for each plot
acf.dat <- sapply(unique(delta.npp.dat.annual.giant.kelp$plot), function(x){
  acf(delta.npp.dat.annual.giant.kelp$resid[delta.npp.dat.annual.giant.kelp$plot == x],
      lag.max = 4, plot = FALSE)$acf
})

pacf.dat <- sapply(unique(delta.npp.dat.annual.giant.kelp$plot), function(x){
  pacf(delta.npp.dat.annual.giant.kelp$resid[delta.npp.dat.annual.giant.kelp$plot == x],
       lag.max = 4, plot = FALSE)$acf
}
)

acf.dat <- data.frame(acf.dat)
pacf.dat <- data.frame(pacf.dat)

acf.dat <- acf.dat %>%
  dplyr::mutate(lag = 1:nrow(acf.dat) - 1) %>%
  tidyr::gather(key = "plot", value = "acf", -lag)

pacf.dat <- pacf.dat %>%
  dplyr::mutate(lag = 1:nrow(pacf.dat)) %>%
  tidyr::gather(key = "plot", value = "pacf", -lag)

acf.dat <- dplyr::left_join(acf.dat, pacf.dat, by = c("lag", "plot"))

# Calculate critical value (based on the lowest length of kelp.kelp.year series available)
acf.dat$crit <- qnorm((1 + 0.95)/2) / sqrt(length(unique(delta.npp.dat.annual.giant.kelp[delta.npp.dat.annual.giant.kelp$plot == "Isla Vista_Annual", ]$kelp.year)))

# Plot ACF by plot
p1 <- ggplot(data = acf.dat, aes(x = lag, y = acf)) +
  ggtitle("Autocorrelation by plot") +
  facet_wrap(~ plot) +
  geom_bar(stat = "identity", width = 0.1, color = "black", fill = "black") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)/10, name = "ACF") +
  scale_x_continuous(breaks = 0:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

# Plot average PACF
p2 <- ggplot(data = acf.dat[!is.na(acf.dat$pacf), ], aes(x = lag, y = pacf)) +
  ggtitle("Average partial autocorrelation across plots") +
  stat_summary(fun.data = mean_cl_boot) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = crit), linetype = "dashed") +
  geom_line(aes(y = -crit), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1), name = "PACF") +
  scale_x_continuous(limits = c(0.95, max(acf.dat$lag)), breaks = 1:max(acf.dat$lag), name = "Lag") +
  theme_classic() +
  theme(aspect.ratio = 1)

(p1 | p2)

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals
giant.kelp.delta.npp.mod.reduced <-  glmmTMB(annual.delta.npp_kgC.m2.y ~
                                       treatment + time.since.start + habitat.quality.cat,
                                     family = gaussian(link = "identity"),
                                     dispformula = ~ habitat.quality.cat,
                                     data = delta.npp.dat.annual.giant.kelp)

preds <- ggemmeans(giant.kelp.delta.npp.mod.reduced, 
                   terms = c("time.since.start[breaks.time]",
                             "treatment",
                             "habitat.quality.cat[breaks.hab]" #"habitat.quality[breaks.hab]", 
                   ),
                   type = "fixed")

preds <- data.frame(preds) %>%
  dplyr::select(-std.error) %>%
  dplyr::mutate(predicted = predicted,
                conf.low = conf.low,
                conf.high = conf.high) %>%
  dplyr::rename(annual.delta.npp_kgC.m2.y = predicted,
                habitat.quality.cat = facet,
                treatment = group,
                time.since.start = x)

preds$time.since.start <- as.numeric(as.character(preds$time.since.start))

preds2 <- preds

preds2$annual.delta.npp_kgC.m2.y[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
preds2$conf.low[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
preds2$conf.high[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA

preds2$habitat.quality.cat <- factor(preds2$habitat.quality.cat, levels = c("Low", "Medium", "High"))
delta.npp.dat.annual.giant.kelp$habitat.quality.cat <- factor(delta.npp.dat.annual.giant.kelp$habitat.quality.cat, 
                                                         levels = c("Low", "Medium", "High"))

ggplot(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  #facet_grid(treatment ~ habitat.quality.cat) +  
  facet_wrap(~habitat.quality.cat) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = treatment, fill = treatment), alpha = 0.25) +
  geom_line(aes(colour = treatment, group = treatment), size = 1)

# Compare trend lines to zero
em <- emtrends(giant.kelp.delta.npp.mod.reduced, ~ habitat.quality.cat | treatment, var = "time.since.start", max.degree = 1)
summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

em <- emtrends(giant.kelp.delta.npp.mod, pairwise ~ habitat.quality.cat | treatment, var = "time.since.start")
summary(em, adjust = 'none')


# Compare mean values to zero
em <- emmeans(giant.kelp.delta.npp.mod, ~ treatment | habitat.quality.cat)
summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

em <- emmeans(giant.kelp.delta.npp.mod.reduced, ~ treatment | habitat.quality.cat)
summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

# Plot
delta.npp.dat.annual.giant.kelp$treatment <- factor(delta.npp.dat.annual.giant.kelp$treatment, levels = c("Continual", "Annual"))
preds2$treatment <- factor(preds2$treatment, levels = c("Continual", "Annual"))

delta.npp.dat.annual.giant.kelp$treatment2 <- factor(mapvalues(delta.npp.dat.annual.giant.kelp$treatment, 
                                                          from = c("Continual", "Annual"),
                                                          to = c("Quarterly disturbance", "Annual disturbance")),
                                                levels = c("Annual disturbance", "Quarterly disturbance"))
preds2$treatment2 <- factor(mapvalues(preds2$treatment, 
                                      from = c("Continual", "Annual"),
                                      to = c("Quarterly disturbance", "Annual disturbance")),
                            levels = c("Annual disturbance", "Quarterly disturbance"))

ggplot(data = delta.npp.dat.annual.giant.kelp, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_grid(treatment2 ~ habitat.quality.cat) +  
  geom_ribbon(data = preds2, 
              aes(ymin = conf.low, ymax = conf.high, group = treatment2, fill = treatment2), alpha = 0.33) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
  geom_line(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y,
                               color = treatment2, group = treatment), size = 1) +
  scale_shape_manual(values = c(21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col2.line, col3.line), name = "Treatment") +
  ylab(expression(atop("Giant kelp annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  xlab("Years since start of experiment") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  scale_y_continuous(breaks = seq(-3, 1, by = 1)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-3.25, 1)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1,
        legend.position = "top", #c(0.81, 0.95),
        legend.text = element_text(size = 11, hjust = 0.5),
        legend.title = element_text(size = 11, hjust = 0.5, face = 'bold'),
        legend.spacing.y = unit(0.15, 'lines'),
        legend.background = element_blank()) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
#ggsave(filename = "Figures/new.giant.kelp.delta-npp--scatter.pdf", width = 8, height = 6)
#ggsave(filename = "Figures/new.giant.kelp.delta-npp--scatter.svg", width = 8, height = 6)

# ------------------------------------------------------------------------------------
# Calculate predictions and confidence intervals - all data pooled together (no facets)
preds <- ggemmeans(giant.kelp.delta.npp.mod.reduced, 
                   terms = "time.since.start[breaks.time]",
                   type = "fixed")

preds <- data.frame(preds) %>%
  dplyr::select(-std.error, -group) %>%
  dplyr::mutate(predicted = predicted,
                conf.low = conf.low,
                conf.high = conf.high) %>%
  dplyr::rename(annual.delta.npp_kgC.m2.y = predicted,
                time.since.start = x)

preds$time.since.start <- as.numeric(as.character(preds$time.since.start))

preds2 <- preds

preds2$annual.delta.npp_kgC.m2.y[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
preds2$conf.low[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA
preds2$conf.high[preds$time.since.start > 7.03 & preds$treatment == "Continual"] <- NA

ggplot(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_line(size = 1)

# Plot
delta.npp.dat.annual.giant.kelp$treatment <- factor(delta.npp.dat.annual.giant.kelp$treatment, levels = c("Continual", "Annual"))
delta.npp.dat.annual.giant.kelp$treatment2 <- factor(mapvalues(delta.npp.dat.annual.giant.kelp$treatment, 
                                                               from = c("Continual", "Annual"),
                                                               to = c("Quarterly disturbance", "Annual disturbance")),
                                                     levels = c("Annual disturbance", "Quarterly disturbance"))

ggplot(data = delta.npp.dat.annual.giant.kelp, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_ribbon(data = preds2, aes(ymin = conf.low, ymax = conf.high), fill = "dimgrey", alpha = 0.33) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
  geom_line(data = preds2, aes(x = time.since.start, y = annual.delta.npp_kgC.m2.y), size = 1) +
  scale_shape_manual(values = c(21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col2.line, col3.line), name = "Treatment") +
  ylab(expression(atop("Canopy annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  xlab("Years since start of experiment") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  scale_y_continuous(breaks = seq(-3, 1, by = 1)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-3.25, 1)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1,
        legend.position = "top", #c(0.81, 0.95),
        legend.text = element_text(size = 11, hjust = 0.5),
        legend.title = element_text(size = 11, hjust = 0.5, face = 'bold'),
        legend.spacing.y = unit(0.15, 'lines'),
        legend.background = element_blank()) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
#ggsave(filename = "Figures/new.giant.kelp.delta-npp--scatter-pooled.pdf", width = 5, height = 6)

# ------------------------------------------------------------------------------------
# Plot means
delta.npp.dat.annual.giant.kelp$treatment2 <- factor(delta.npp.dat.annual.giant.kelp$treatment2,
                                                levels = c("Annual disturbance", "Quarterly disturbance"))

delta.npp.dat.annual.giant.kelp$treatment3 <- factor(mapvalues(delta.npp.dat.annual.giant.kelp$treatment2,
                                                          from = c("Annual disturbance", "Quarterly disturbance"),
                                                          to = c("Annual", "Quarterly")),
                                                levels = c("Annual", "Quarterly"))

# Check significance of bars from zero
em <- emmeans(giant.kelp.delta.npp.mod, ~ treatment | habitat.quality.cat, infer = T, adjust = 'none')
summary(em, infer=c(TRUE,TRUE), null=0, adjust = 'none')

ggplot(data = delta.npp.dat.annual.giant.kelp, aes(x = habitat.quality.cat, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  #ylab(expression(paste("Change in giant.kelp NPP from control (kg C/", "m"^2, "/y)", sep = ""))) +
  ylab(expression(atop("Canopy annual NPP", paste("(difference from control; kg C/", "m"^2, "/y)", sep = "")))) +
  xlab("Site habitat quality") +
  scale_y_continuous(breaks = seq(-3, 1, by = 0.5)) + 
  coord_cartesian(ylim = c(-2.75, 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 11, hjust = 0, face = 'bold'),
        legend.spacing.y = unit(0.2, 'lines'),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        aspect.ratio = 1/2) 
#ggsave(filename = "Figures/new.giant.kelp.delta-npp--bar.pdf", width = 6, height = 4)

# ------------------------------------------------------------------------------------
# Plot means of giant kelp annual NPP and total macroalgal annual NPP together

a <- ggplot(data = delta.npp.dat.annual.giant.kelp, aes(x = habitat.quality.cat, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  ylab("") +
  scale_x_discrete(labels = NULL) +
  ggtitle("(A) Giant kelp NPP") +
  xlab("") +
  scale_y_continuous(breaks = seq(-3, 1, by = 0.5)) + 
  coord_cartesian(ylim = c(-2.75, 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 11, hjust = 0, face = 'bold'),
        legend.spacing.y = unit(0.2, 'lines'),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        aspect.ratio = 1/2) 
a

b <- ggplot(data = delta.npp.dat.annual.total, aes(x = habitat.quality.cat, y = annual.delta.npp_kgC.m2.y)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col2, col3), name = "Treatment") +
  ylab(expression(paste("Annual NPP (difference from control; kg C/", "m"^2, "/y)", sep = ""))) +
  ggtitle("(B) Total macroalgal NPP") +
  xlab("Site habitat quality") +
  scale_y_continuous(breaks = seq(-3, 1, by = 0.5)) + 
  coord_cartesian(ylim = c(-2.75, 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14, hjust = -0.46),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/2) +
  guides(fill = F)
b

(a/b)
#ggsave(filename = "Figures/new.combo.delta-npp--bar.pdf", width = 6, height = 8)

# ====================================================================================
# PLOT AVERAGE NPP BY LAYER AS A FUNCTION OF HABITAT QUALITY AND DISTURBANCE TREATMENT
# ====================================================================================

#----------------------------------------------------------------------------------------------
transect.dat.annual$habitat.quality.cat <- factor(plyr::mapvalues(transect.dat.annual$site,
                                                                  from = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"),
                                                                  to = c("High", "High", "Medium", "Medium", "Low")),
                                                  levels = c("Low", "Medium", "High"))

transect.dat.annual$treatment2 <- dplyr::recode(transect.dat.annual$treatment, 
                                                "Control" = "Control", 
                                                "Annual" = "Annual disturbance", 
                                                "Continual" = "Quarterly disturbance")
transect.dat.annual$treatment2 <- factor(transect.dat.annual$treatment2, 
                                         levels = c("Control", "Annual disturbance", "Quarterly disturbance"))

#----------------------------------------------------------------------------------------------
a <- ggplot(data = transect.dat.annual, aes(x = habitat.quality.cat, y = giant.kelp.npp_kgC.m2.y)) +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  ylab("") +
  scale_x_discrete(labels = NULL) +
  ggtitle("(A) Giant kelp NPP") +
  xlab("") +
  scale_y_continuous(breaks = seq(0, 3, by = 1), expand = expansion(mult = c(0.005, 0.005))) + 
  coord_cartesian(ylim = c(0, 3)) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(size = 14, hjust = 0),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        legend.position = c(0.3, 0.7),
        legend.text = element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 11, hjust = 0, face = 'bold'),
        legend.spacing.y = unit(0.2, 'lines'),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        aspect.ratio = 1/2) 
a

b <- ggplot(data = transect.dat.annual, aes(x = habitat.quality.cat, y = understory.npp_kgC.m2.y)) +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  ylab(expression(paste("Annual NPP (kg C/", "m"^2, "/y)", sep = ""))) +
  scale_x_discrete(labels = NULL) +
  ggtitle("(B) Understory NPP") +
  xlab("") +
  scale_y_continuous(breaks = seq(0, 4, by = 1), expand = expansion(mult = c(0.005, 0.005))) + 
  #scale_y_continuous(breaks = seq(0, 2, by = 0.5), expand = expansion(mult = c(0.005, 0.1))) + 
  coord_cartesian(ylim = c(0, 3)) +
  #coord_cartesian(ylim = c(0, 1.6)) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(size = 14, hjust = 0),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        aspect.ratio = 1/2) +
  guides(fill = F)
b

c <- ggplot(data = transect.dat.annual, aes(x = habitat.quality.cat, y = total.npp_kgC.m2.y)) +
  stat_summary(aes(fill = treatment2), geom = "bar", fun = "mean", color = "black", width = 0.65, 
               position=position_dodge(width = 0.65)) +
  stat_summary(aes(group = treatment2), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, fun.args = list(B = 1000), 
               position=position_dodge(width = 0.65)) +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  ylab("") +
  ggtitle("(C) Total macroalgal NPP") +
  xlab("Site habitat quality") +
  scale_y_continuous(breaks = seq(0, 4, by = 1), expand = expansion(mult = c(0.005, 0.005))) + 
  coord_cartesian(ylim = c(0, 3)) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(size = 14, hjust = 0),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(size = 14, hjust = -0.46),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1/2) +
  guides(fill = F)
c

(a/b/c)
#ggsave(filename = "Figures/new.combo.mean-npp--bar--fixed-scale.pdf", width = 4.25, height = 7.5)

#----------------------------------------------------------------------------------------------
# Plot over time
transect.dat.annual <- transect.dat.annual %>%
  dplyr::mutate(plot = paste(site, treatment, sep = "_")) %>%
  dplyr::group_by(plot) %>%
  dplyr::mutate(min.time = min(kelp.year)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(time.since.start = kelp.year - min.time + 1)
                
transect.dat.annual.mean <- transect.dat.annual %>%
  dplyr::group_by(time.since.start, habitat.quality.cat, treatment2) %>%
  dplyr::summarize_at(vars(giant.kelp.npp_kgC.m2.y, understory.npp_kgC.m2.y, total.npp_kgC.m2.y),
                     mean, na.rm = T) %>%
  dplyr::ungroup()

transect.dat.annual.error <- transect.dat.annual %>%
  dplyr::group_by(time.since.start, habitat.quality.cat, treatment2) %>%
  dplyr::summarize_at(vars(giant.kelp.npp_kgC.m2.y, understory.npp_kgC.m2.y, total.npp_kgC.m2.y),
                      function(x){sd(x) / sqrt(length(x))}) %>%
  dplyr::ungroup() %>%
  dplyr::rename(se.gk = giant.kelp.npp_kgC.m2.y,
                se.und = understory.npp_kgC.m2.y,
                se.tot = total.npp_kgC.m2.y)
transect.dat.annual.error <- left_join(transect.dat.annual.mean, transect.dat.annual.error)

transect.dat.annual.error <- transect.dat.annual.error %>%
  dplyr::mutate(upper.gk = giant.kelp.npp_kgC.m2.y + se.gk,
                lower.gk = giant.kelp.npp_kgC.m2.y - se.gk,
                upper.und = understory.npp_kgC.m2.y + se.und,
                lower.und = understory.npp_kgC.m2.y - se.und,
                upper.tot = total.npp_kgC.m2.y + se.tot,
                lower.tot = total.npp_kgC.m2.y - se.tot) 

#----------------------------------------------------------------------------------------------
d <- ggplot(data = transect.dat.annual.mean, aes(x = time.since.start, y = giant.kelp.npp_kgC.m2.y)) +
  facet_wrap(~ habitat.quality.cat, ncol = 3) +  
  geom_line(aes(color = treatment2, group = treatment2), size = 1) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
  scale_shape_manual(values = c(22, 21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col1.line, col2.line, col3.line), name = "Treatment") +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  #scale_y_continuous(breaks = seq(-3, 1, by = 1)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(0, 4.25)) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=14, color = "black", hjust = 0.5),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
d

e <- ggplot(data = transect.dat.annual.mean, aes(x = time.since.start, y = understory.npp_kgC.m2.y)) +
  facet_wrap(~ habitat.quality.cat, ncol = 3) +  
  geom_line(aes(color = treatment2, group = treatment2), size = 1) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
  scale_shape_manual(values = c(22, 21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col1.line, col2.line, col3.line), name = "Treatment") +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  #scale_y_continuous(breaks = seq(-3, 1, by = 1)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(0, 4.25)) +
  #coord_cartesian(xlim = c(0, 10.5), ylim = c(0, 2.5)) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_blank(),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
e

f <- ggplot(data = transect.dat.annual.mean, aes(x = time.since.start, y = total.npp_kgC.m2.y)) +
  facet_wrap(~ habitat.quality.cat, ncol = 3) +  
  geom_line(aes(color = treatment2, group = treatment2), size = 1) +
  geom_point(aes(fill = treatment2, shape = treatment2), size = 2) +
  scale_shape_manual(values = c(22, 21, 24), name = "Treatment") +
  scale_fill_manual(values = c(col1, col2, col3), name = "Treatment") +
  scale_color_manual(values = c(col1.line, col2.line, col3.line), name = "Treatment") +
  ylab("") +
  xlab("Years since start of experiment") +
  scale_x_continuous(breaks = seq(0, 10, by =2)) +
  #scale_y_continuous(breaks = seq(-3, 1, by = 1)) + 
  coord_cartesian(xlim = c(0, 10.5), ylim = c(0, 4.25)) +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_blank(),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) +
  guides(fill = F, color = F, shape = F)+
  ggtitle('Site habitat quality')
f

(d/e/f)
#ggsave(filename = "Figures/new.combo.mean-npp--over-time--fixed-scale.pdf", width = 4.25, height = 7.5)

((a/b/c) | (d/e/f)) +
  plot_layout(widths = c(1, 1.5))
#ggsave(filename = "Figures/new.combo.mean-npp--fixed-scale.pdf", width = 9.25, height = 7.5)


# ====================================================================================
# CALCULATE AMOUNT OF KELP BIOMASS REMOVED
# ====================================================================================

# Data wrangling
giant.kelp.biomass.transect.dat.raw$habitat.quality.cat <- factor(plyr::mapvalues(giant.kelp.biomass.transect.dat.raw$site,
                                                           from = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"),
                                                           to = c("High", "High", "Medium", "Medium", "Low")),
                                           levels = c("Low", "Medium", "High"))

transect.dat$habitat.quality.cat <- factor(plyr::mapvalues(transect.dat$site,
                                                                  from = c("Mohawk", "Isla Vista", "Arroyo Quemado", "Naples", "Carpinteria"),
                                                                  to = c("High", "High", "Medium", "Medium", "Low")),
                                                  levels = c("Low", "Medium", "High"))

giant.kelp.biomass.transect.dat.raw <- giant.kelp.biomass.transect.dat.raw %>%
  dplyr::select(-scientific.name) %>%
  dplyr::rename(giant.kelp.dry.gm2 = dry.gm2) %>%
  dplyr::select(treatment, year, season, time, site, habitat.quality.cat, giant.kelp.dry.gm2)

# Amount of kelp removed from quarterly removal plots = amount of kelp present at each quarterly sampling
# However, sampling was done twice per season in the early years. So, we need to SUM the total giant kelp biomass within each season, which will work for both early years and late years
kelp.removed.dat_quarterly <- giant.kelp.biomass.transect.dat.raw %>%
  dplyr::filter(treatment == "Continual") %>%
  droplevels() %>%
  dplyr::group_by(treatment, year, season, time, site, habitat.quality.cat) %>%
  dplyr::summarise(giant.kelp.dry.gm2 = sum(giant.kelp.dry.gm2, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Amount of kelp removed from quarterly removal plots = amount of kelp present at winter sampling (just before it was removed) 
kelp.removed.dat_annual <- transect.dat %>%
  dplyr::select(year:transect, habitat.quality.cat, giant.kelp.dry.gm2) %>%
  dplyr::filter(treatment == "Annual") #%>%
  dplyr::filter(season == "1-Winter") %>%
  droplevels()

# Plot
g1 <- ggplot(data = kelp.removed.dat_annual, aes(x = habitat.quality.cat, y = giant.kelp.dry.gm2/1000)) +
  geom_violin(fill = col2, alpha = 0.75) +
  geom_boxplot(fill = 'black', color = 'black', width = 0.05) +
  stat_summary(fun = median, na.rm = T, color = 'white', size = 0.25) +
  #geom_jitter(width = 0.1)
  theme_classic() +
  ylab(expression(atop("Giant kelp biomass removed ", paste("(kg dry/", "m"^2, "/year)", sep = "")))) +
  xlab("") +
  #xlab("\nSite habitat quality") +
  ggtitle("Annual disturbance plots") +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), expand = expansion(mult = c(0.003, 0.003))) +
  coord_cartesian(ylim = c(0, 1.2)) +
  theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
        plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
g1

g2 <- ggplot(data = kelp.removed.dat_quarterly, aes(x = habitat.quality.cat, y = giant.kelp.dry.gm2/1000)) +
  geom_violin(fill = col3, alpha = 0.75) +
  geom_boxplot(fill = 'black', color = 'black', width = 0.05) +
  stat_summary(fun = median, na.rm = T, color = 'white', size = 0.25) +
  #geom_jitter(width = 0.1)
  theme_classic() +
  ylab(expression(atop("Giant kelp biomass removed ", paste("(kg dry/", "m"^2, "/quarter)", sep = "")))) +
  xlab("\nSite habitat quality") +
  ggtitle("Quarterly disturbance plots") +
  scale_y_continuous(breaks = seq(0, 0.2, 0.05), expand = expansion(mult = c(0.003, 0.1))) +
  coord_cartesian(ylim = c(0, 0.21)) +
  theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
        plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", size=1, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.text.y = element_text(margin=margin(0.1,0.1,0.1,0.1, "cm"), 
                                   color = "black", size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        aspect.ratio = 1) 
g2

(g1 / g2) 
ggsave(filename = paste0(here::here("/Figures/"),  "giant.kelp.biomass.removed_variable.scale", ".pdf"), height = 8, width = 4)


g3 <- g2 +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), expand = expansion(mult = c(0.003, 0.003))) +
  coord_cartesian(ylim = c(0, 1.2)) 

(g1 / g3)
ggsave(filename = paste0(here::here("/Figures/"),  "giant.kelp.biomass.removed_fixed.scale", ".pdf"), height = 8, width = 4)


# ====================================================================================
# COMPARE ABSOLUTE NPP *WITHIN A TREATMENT* ACROSS HABITAT TYPES
# ====================================================================================

temp.dat_annual <- transect.dat.annual %>%
  dplyr::filter(treatment == "Annual") %>%
  droplevels()

temp.dat_continual <- transect.dat.annual %>%
  dplyr::filter(treatment == "Continual") %>%
  droplevels()
 
# ------------------------------------------------------------------------------------------------------------------------
# Compare total NPP among habitat quality types for ANNUAL DISTURBANCE
new.mod1 <- glmmTMB(total.npp_kgC.m2.y ~ habitat.quality.cat +
                      ar1(factor(kelp.year) + 0 | site),
                    data = temp.dat_annual)

new.mod1 <- glmmTMB(total.npp_kgC.m2.y ~ habitat.quality.cat * treatment +
                      ar1(factor(kelp.year) + 0 | site),
                    data = transect.dat.annual)

# Validate  
simulationOutput <- simulateResiduals(new.mod1, n=1000)
plot(simulationOutput)
testDispersion(simulationOutput)  
testZeroInflation(simulationOutput)

summary(new.mod1)
Anova(new.mod1)

em1 <- emmeans(new.mod1, pairwise ~ habitat.quality.cat * treatment)
summary(em1, infer=c(TRUE,TRUE), adjust = 'none')
plot(em1)

# ------------------------------------------------------------------------------------------------------------------------
# Compare total NPP among habitat quality types for ANNUAL DISTURBANCE
new.mod2 <- glmmTMB(total.npp_kgC.m2.y ~ habitat.quality.cat +
                      ar1(factor(kelp.year) + 0 | site),
                    data = temp.dat_continual)

# Validate  
simulationOutput <- simulateResiduals(new.mod1, n=1000)
plot(simulationOutput)
testDispersion(simulationOutput)  
testZeroInflation(simulationOutput)

summary(new.mod2)
Anova(new.mod2)

em2 <- emmeans(new.mod2, pairwise ~ habitat.quality.cat)
summary(em2, infer=c(TRUE,TRUE), adjust = 'none')
plot(em2)

# ====================================================================================
# CONSOLIDATE STATISTICAL RESULTS
# ====================================================================================

# -------------------------------------------------------------------------------------------
# Giant kelp habitat quality

#site.urchin.mod
p.vals1 <- data.frame('contrast' = NA,
                      'model.type' = 'GLS',
                      'family' = 'Gaussian',
                      'link' = 'identity',
                      'AR' = 1,
                      'trans' = 'sqrt',
                      'response.var' = 'urchin.dens.no.m2', 
                      'ZI.formula' = NA,
                      'disp.formula' = '~ factor(year)',
                      'predictor.var' = rownames(data.frame(Anova(site.urchin.mod))),
                      'num.df' = Anova(site.urchin.mod)$Df,
                      'den.df' = as.numeric(summary(site.urchin.mod)$AICtab["df.resid"]),
                      'test.statistic' = 'Chisq',
                      'test.val' = Anova(site.urchin.mod)$Chisq,
                      'p' = Anova(site.urchin.mod)$`Pr(>Chisq)`)

#site.sand.mod
p.vals2 <- data.frame('contrast' = NA,
                      'model.type' = 'GLM',
                         'family' = 'beta',
                         'link' = "logit",
                         'AR' = NA,
                         'trans' = NA,
                         'response.var' = 'percent.sand/100', 
                         'ZI.formula' = "~ 1",
                         'disp.formula' = '~ factor(year)',
                         'predictor.var' = rownames(data.frame(Anova(site.sand.mod))),
                         'num.df' = Anova(site.sand.mod)$Df,
                         'den.df' = as.numeric(summary(site.sand.mod)$AICtab["df.resid"]),
                         'test.statistic' = 'Chisq',
                         'test.val' = Anova(site.sand.mod)$Chisq,
                         'p' = Anova(site.sand.mod)$`Pr(>Chisq)`)

#site.tot.mod
p.vals3 <- data.frame('contrast' = NA,
                      'model.type' = 'GLS',
                         'family' = 'Gaussian',
                         'link' = "identity",
                         'AR' = 1,
                         'trans' = 'sqrt',
                         'response.var' = 'total.algae.dry.kgm2', 
                         'ZI.formula' = NA,
                         'disp.formula' = NA,
                         'predictor.var' = rownames(data.frame(Anova(site.tot.mod))),
                         'num.df' = Anova(site.tot.mod)$Df,
                         'den.df' = as.numeric(summary(site.tot.mod)$AICtab["df.resid"]),
                         'test.statistic' = 'Chisq',
                         'test.val' = Anova(site.tot.mod)$Chisq,
                         'p' = Anova(site.tot.mod)$`Pr(>Chisq)`)

#site.preds.mod
p.vals4 <- data.frame('contrast' = NA,
                      'model.type' = 'GLS',
                         'family' = 'Gaussian',
                         'link' = 'identity',
                         'AR' = NA,
                         'trans' = NA,
                         'response.var' = 'env.quality', 
                         'ZI.formula' = NA,
                         'disp.formula' = NA,
                         'predictor.var' = rownames(data.frame(Anova(site.preds.mod)))[1],
                         'num.df' = Anova(site.preds.mod)$Df[1],
                         'den.df' = Anova(site.preds.mod)$Df[2],
                         'test.statistic' = 'F',
                         'test.val' = Anova(site.preds.mod)$`F value`[1],
                         'p' = Anova(site.preds.mod)$`Pr(>F)`[1])

#qual.mod
p.vals5 <- data.frame('contrast' = NA,
                      'model.type' = 'GLM',
                         'family' = 'Gamma',
                         'link' = 'log',
                         'AR' = NA,
                         'trans' = NA,
                         'response.var' = 'total.algae.dry.kgm2', 
                         'ZI.formula' = '~ urchin.dens.no.m2',
                         'disp.formula' = NA,
                         'predictor.var' = rownames(data.frame(Anova(qual.mod))),
                         'num.df' = Anova(qual.mod)$Df,
                         'den.df' = as.numeric(summary(qual.mod)$AICtab["df.resid"]),
                         'test.statistic' = 'Chisq',
                         'test.val' = Anova(qual.mod)$Chisq,
                         'p' = Anova(qual.mod)$`Pr(>Chisq)`)

# -------------------------------------------------------------------------------------------
# Posthoc tests for site - habitat quality models
site.urchin.mod <- glmmTMB(sqrt(urchin.dens.no.m2) ~ site + poly(time, 2) + 
                             ar1(time.factor + 0 | site),
                           family = gaussian(link = "identity"),
                           dispformula = ~ factor(year),
                           data = quality.dat.tot.wide)

site.sand.mod <- glmmTMB(percent.sand/100 ~ site, 
                         family = beta_family(link = "logit"),
                         ziformula = ~ 1,
                         data = quality.dat.tot.wide)

site.tot.mod <- glmmTMB(sqrt(total.algae.dry.kgm2) ~ site + 
                          ar1(time.factor + 0 | site),
                        family = gaussian(link = "identity"),
                        data = quality.dat.tot.wide)

site.preds.mod <- lm(preds ~ site, 
                     data = quality.dat.tot.wide)

a <- data.frame(pairs(emmeans(site.urchin.mod, ~ site), simple = "site", reverse = TRUE, adjust = 'none'))
b <- data.frame(pairs(emmeans(site.sand.mod, ~ site), simple = "site", reverse = TRUE, adjust = 'none'))
c <- data.frame(pairs(emmeans(site.tot.mod, ~ site), simple = "site", reverse = TRUE, adjust = 'none'))
d <- data.frame(pairs(emmeans(site.preds.mod, ~ site), simple = "site", reverse = TRUE, adjust = 'none'))

p.vals6 <- data.frame('contrast' = c(a$contrast, b$contrast, c$contrast, d$contrast),
                      'model.type' = 'posthoc.test',
                      'family' = NA,
                      'link' = NA,
                      'AR' = NA,
                      'trans' = NA,
                      'response.var' = c(rep('urchin.dens.no.m2', 10),
                                         rep('percent.sand/100', 10),
                                         rep('total.algae.dry.kgm2', 10),
                                         rep('env.quality', 10)), 
                      'ZI.formula' = NA,
                      'disp.formula' = NA,
                      'predictor.var' = 'site',
                      'num.df' = 1,
                      'den.df' = c(a$df, b$df, c$df, d$df),
                      'test.statistic' = 't',
                      'test.val' = abs(c(a$t.ratio, b$t.ratio, c$t.ratio, d$t.ratio)),
                      'p' = c(a$p.value, b$p.value, c$p.value, d$p.value))

rm(a, b, c, d)

# -------------------------------------------------------------------------------------------
# Light model
p.vals7 <- data.frame('contrast' = NA,
                      'model.type' = 'GLM',
                       'family' = 'Gamma',
                       'link' = 'log',
                       'AR' = 1,
                       'trans' = NA,
                       'response.var' = 'light.mol.day', 
                       'ZI.formula' = NA,
                       'disp.formula' = '~ site + season',
                       'predictor.var' = rownames(data.frame(Anova(par.mod))),
                       'num.df' = Anova(par.mod)$Df,
                       'den.df' = as.numeric(summary(par.mod)$AICtab["df.resid"]),
                       'test.statistic' = 'Chisq',
                       'test.val' = Anova(par.mod)$Chisq,
                       'p' = Anova(par.mod)$`Pr(>Chisq)`)

# -------------------------------------------------------------------------------------------
# Delta-understory NPP as a function of treatment
p.vals8 <- data.frame('contrast' = NA,
                       'model.type' = 'GLM',
                       'family' = 'Gamma',
                       'link' = 'log',
                       'AR' = NA,
                       'trans' = NA,
                       'response.var' = 'understory.annual.delta.npp_kgC.m2.y', 
                       'ZI.formula' = NA,
                       'disp.formula' = NA,
                       'predictor.var' = rownames(data.frame(Anova(understory.delta.npp.mod))),
                       'num.df' = Anova(understory.delta.npp.mod)$Df,
                       'den.df' = as.numeric(summary(understory.delta.npp.mod)$AICtab["df.resid"]),
                       'test.statistic' = 'Chisq',
                       'test.val' = Anova(understory.delta.npp.mod)$Chisq,
                       'p' = Anova(understory.delta.npp.mod)$`Pr(>Chisq)`)

# Posthoc test to examine differences from zero
em <- data.frame(summary(emmeans(understory.delta.npp.mod, ~ treatment | habitat.quality.cat), 
                         infer=c(TRUE,TRUE), null=0, adjust = 'none'))

p.vals9 <- data.frame('contrast' = paste("Difference from zero: ", 
                                          em$treatment, " removal in ", em$habitat.quality.cat, 
                                          " quality habitat", sep = ""),
                       'model.type' = 'posthoc.test',
                       'family' = NA,
                       'link' = NA,
                       'AR' = NA,
                       'trans' = NA,
                       'response.var' = "understory.annual.delta.npp_kgC.m2.y", 
                       'ZI.formula' = NA,
                       'disp.formula' = NA,
                       'predictor.var' = 'habitat.quality and treatment',
                       'num.df' = 1,
                       'den.df' = em$df,
                       'test.statistic' = 't',
                       'test.val' = abs(em$t.ratio),
                       'p' = em$p.value)

# Posthoc test to examine changes over time-since-start
em <- data.frame(summary(emtrends(understory.delta.npp.mod, ~ habitat.quality.cat | treatment, 
                                  var = "time.since.start"), 
                         infer=c(TRUE,TRUE), null=0, adjust = 'none'))

p.vals10 <- data.frame('contrast' = paste("Trend over time (vs. zero): ",
                                         em$treatment, " removal in ", em$habitat.quality.cat, 
                                         " quality habitat", sep = ""),
                       'model.type' = 'posthoc.test',
                       'family' = NA,
                       'link' = NA,
                       'AR' = NA,
                       'trans' = NA,
                       'response.var' = "understory.annual.delta.npp_kgC.m2.y", 
                       'ZI.formula' = NA,
                       'disp.formula' = NA,
                       'predictor.var' = 'habitat.quality and treatment',
                       'num.df' = 1,
                       'den.df' = em$df,
                       'test.statistic' = 't',
                       'test.val' = abs(em$t.ratio),
                       'p' = em$p.value)

# -------------------------------------------------------------------------------------------
# Delta-giant kelp NPP as a function of treatment
p.vals11 <- data.frame('contrast' = NA,
                       'model.type' = 'GLS',
                       'family' = 'Gaussian',
                       'link' = 'identity',
                       'AR' = NA,
                       'trans' = NA,
                       'response.var' = 'giant.kelp.annual.delta.npp_kgC.m2.y', 
                       'ZI.formula' = NA,
                       'disp.formula' = '~ habitat.quality.cat',
                       'predictor.var' = rownames(data.frame(Anova(giant.kelp.delta.npp.mod))),
                       'num.df' = Anova(giant.kelp.delta.npp.mod)$Df,
                       'den.df' = as.numeric(summary(giant.kelp.delta.npp.mod)$AICtab["df.resid"]),
                       'test.statistic' = 'Chisq',
                       'test.val' = Anova(giant.kelp.delta.npp.mod)$Chisq,
                       'p' = Anova(giant.kelp.delta.npp.mod)$`Pr(>Chisq)`)

# Posthoc test to examine changes over time-since-start
em <- data.frame(summary(emtrends(giant.kelp.delta.npp.mod.reduced, ~ habitat.quality.cat | treatment, 
                                  var = "time.since.start"), 
                         infer=c(TRUE,TRUE), null=0, adjust = 'none'))

p.vals12 <- data.frame('contrast' = paste("Trend over time (vs. zero): ", 
                                         em$treatment, " removal in ", em$habitat.quality.cat, 
                                         " quality habitat", sep = ""),
                      'model.type' = 'posthoc.test',
                      'family' = NA,
                      'link' = NA,
                      'AR' = NA,
                      'trans' = NA,
                      'response.var' = "giant.kelp.annual.delta.npp_kgC.m2.y", 
                      'ZI.formula' = NA,
                      'disp.formula' = NA,
                      'predictor.var' = 'habitat.quality and treatment',
                      'num.df' = 1,
                      'den.df' = em$df,
                      'test.statistic' = 't',
                      'test.val' = abs(em$t.ratio),
                      'p' = em$p.value)

# Posthoc test to examine differences from zero
em <- data.frame(summary(emmeans(giant.kelp.delta.npp.mod, ~ treatment | habitat.quality.cat), 
                         infer=c(TRUE,TRUE), null=0, adjust = 'none'))

p.vals13 <- data.frame('contrast' = paste("Difference from zero: ", 
                                          em$treatment, " removal in ", em$habitat.quality.cat, 
                                          " quality habitat", sep = ""),
                       'model.type' = 'posthoc.test',
                       'family' = NA,
                       'link' = NA,
                       'AR' = NA,
                       'trans' = NA,
                       'response.var' = "giant.kelp.annual.delta.npp_kgC.m2.y", 
                       'ZI.formula' = NA,
                       'disp.formula' = NA,
                       'predictor.var' = 'habitat.quality and treatment',
                       'num.df' = 1,
                       'den.df' = em$df,
                       'test.statistic' = 't',
                       'test.val' = abs(em$t.ratio),
                       'p' = em$p.value)

# -------------------------------------------------------------------------------------------
# Delta-total NPP as a function of treatment

p.vals14 <- data.frame('contrast' = NA,
                      'model.type' = 'GLS',
                      'family' = 'Gaussian',
                      'link' = 'identity',
                      'AR' = NA,
                      'trans' = NA,
                      'response.var' = 'total.annual.delta.npp_kgC.m2.y', 
                      'ZI.formula' = NA,
                      'disp.formula' = '~ habitat.quality.cat',
                      'predictor.var' = rownames(data.frame(Anova(total.delta.npp.mod))),
                      'num.df' = Anova(total.delta.npp.mod)$Df,
                      'den.df' = as.numeric(summary(total.delta.npp.mod)$AICtab["df.resid"]),
                      'test.statistic' = 'Chisq',
                      'test.val' = Anova(total.delta.npp.mod)$Chisq,
                      'p' = Anova(total.delta.npp.mod)$`Pr(>Chisq)`)

# Posthoc test to examine differences from zero
em <- data.frame(summary(emmeans(total.delta.npp.mod, ~ treatment | habitat.quality.cat, 
                                 infer = T, adjust = 'none'), 
                         infer=c(TRUE,TRUE), null=0, adjust = 'none'))

p.vals15 <- data.frame('contrast' = paste("Difference from zero: ", em$treatment, 
                                          " removal in ", em$habitat.quality.cat, " quality habitat",
                                          sep = ""),
                      'model.type' = 'posthoc.test',
                      'family' = NA,
                      'link' = NA,
                      'AR' = NA,
                      'trans' = NA,
                      'response.var' = "total.annual.delta.npp_kgC.m2.y", 
                      'ZI.formula' = NA,
                      'disp.formula' = NA,
                      'predictor.var' = 'habitat.quality and treatment',
                      'num.df' = 1,
                      'den.df' = em$df,
                      'test.statistic' = 't',
                      'test.val' = abs(em$t.ratio),
                      'p' = em$p.value)

# -------------------------------------------------------------------------------------------
# Understory NPP as a function of giant kelp biomass
p.vals16 <- data.frame('contrast' = NA,
                       'model.type' = 'GLM',
                       'family' = 'Gamma',
                       'link' = 'log',
                       'AR' = 1,
                       'trans' = NA,
                       'response.var' = 'understory.npp_kgC.m2.y', 
                       'ZI.formula' = NA,
                       'disp.formula' = NA,
                       'predictor.var' = rownames(data.frame(Anova(understory.npp.vs.giant.kelp.mod))),
                       'num.df' = Anova(understory.npp.vs.giant.kelp.mod)$Df,
                       'den.df' = as.numeric(summary(understory.npp.vs.giant.kelp.mod)$AICtab["df.resid"]),
                       'test.statistic' = 'Chisq',
                       'test.val' = Anova(understory.npp.vs.giant.kelp.mod)$Chisq,
                       'p' = Anova(understory.npp.vs.giant.kelp.mod)$`Pr(>Chisq)`)

# -------------------------------------------------------------------------------------------
# Total NPP as a function of giant kelp biomass
p.vals17 <- data.frame('contrast' = NA,
                       'model.type' = 'GLS',
                       'family' = 'Gaussian',
                       'link' = 'identity',
                       'AR' = 1,
                       'trans' = NA,
                       'response.var' = 'total.npp.season.gc.m2.day', 
                       'ZI.formula' = NA,
                       'disp.formula' = "~ site + treatment",
                       'predictor.var' = rownames(data.frame(Anova(total.npp.vs.giant.kelp.mod))),
                       'num.df' = Anova(total.npp.vs.giant.kelp.mod)$Df,
                       'den.df' = as.numeric(summary(total.npp.vs.giant.kelp.mod)$AICtab["df.resid"]),
                       'test.statistic' = 'Chisq',
                       'test.val' = Anova(total.npp.vs.giant.kelp.mod)$Chisq,
                       'p' = Anova(total.npp.vs.giant.kelp.mod)$`Pr(>Chisq)`)

# -------------------------------------------------------------------------------------------
# Combine 
p.vals <- rbind(p.vals1, p.vals2, p.vals3, p.vals4, p.vals5, p.vals6, p.vals7,
                p.vals8, p.vals9, p.vals10, p.vals11, p.vals12, p.vals13, p.vals14,
                p.vals15, p.vals16, p.vals17)

p.vals$p.benj.hoch <- p.adjust(p.vals$p, method = "BH")

# Round and format
p.vals$test.val <- round(p.vals$test.val, 1)

p.vals$p[p.vals$p <= 0.1] <- round(p.vals$p[p.vals$p <= 0.1], 3) 
p.vals$p[p.vals$p > 0.1 & p.vals$p <= 0.2] <- round(p.vals$p[p.vals$p > 0.1 & p.vals$p <= 0.2], 2) 
p.vals$p[p.vals$p > 0.2] <- round(p.vals$p[p.vals$p > 0.2], 1) 

p.vals$p.benj.hoch[p.vals$p.benj.hoch <= 0.1] <- round(p.vals$p.benj.hoch[p.vals$p.benj.hoch <= 0.1], 3) 
p.vals$p.benj.hoch[p.vals$p.benj.hoch > 0.1 & p.vals$p.benj.hoch <= 0.2] <- round(p.vals$p.benj.hoch[p.vals$p.benj.hoch > 0.1 & p.vals$p.benj.hoch <= 0.2], 2) 
p.vals$p.benj.hoch[p.vals$p.benj.hoch > 0.2] <- round(p.vals$p.benj.hoch[p.vals$p.benj.hoch > 0.2], 1) 

# Which values change with p-value adjustment for multiple tests?
p.vals$flag <- 1
p.vals$flag[p.vals$p < 0.05 & p.vals$p.benj.hoch < 0.05] <- 0
p.vals$flag[p.vals$p > 0.05 & p.vals$p.benj.hoch > 0.05] <- 0

p.vals$p[p.vals$p == 0.000] <- "< 0.001"  
p.vals$p.benj.hoch[p.vals$p.benj.hoch == 0.000] <- "< 0.001"  

#write.csv(p.vals, file = here::here('Results/pvals.csv'))


# ====================================================================================
# CALCULATE THE PROPORTION OF TOTAL NPP COMPRISED BY UNDERSTORY
# ====================================================================================

transect.dat.annual$prop.understory <- transect.dat.annual$understory.npp_kgC.m2.y / transect.dat.annual$total.npp_kgC.m2.y

hist(transect.dat.annual$prop.understory)

range(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control' & 
                                            transect.dat.annual$habitat.quality.cat == 'Low'])
range(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control' & 
                                            transect.dat.annual$habitat.quality.cat == 'Medium'])
range(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control' & 
                                            transect.dat.annual$habitat.quality.cat == 'High'])

mean(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control' & 
                                            transect.dat.annual$habitat.quality.cat == 'Low'])
mean(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control' & 
                                            transect.dat.annual$habitat.quality.cat == 'Medium'])
mean(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control' & 
                                            transect.dat.annual$habitat.quality.cat == 'High'])
sd(transect.dat.annual$prop.understory[transect.dat.annual$treatment == 'Control'])


ggplot(data = transect.dat.annual, aes(x = prop.understory, y = ..density..)) +
  geom_histogram() +
  facet_grid(treatment2 ~ habitat.quality.cat)


