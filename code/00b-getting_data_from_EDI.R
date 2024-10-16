
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script includes code to get data from EDI. It only needs to be run
# once to get the data. For all subsequent uses of the code, you can start
# with `00a-set_up.R`.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. libraries -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# file path organization
library(here)

# pulling from the EDI
library(EDIutils)

# general use
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------ 2. SBC LTER data on EDI ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SBC LTE package ID
packageID <- "knb-lter-sbc.119.12"

# downloads zipped file
read_data_package_archive(packageID, path = here("data", "all-species-biomass"))

# unzips file
unzip(zipfile = here("data", "all-species-biomass", "knb-lter-sbc.119.12.zip"),
      exdir = here("data", "all-species-biomass", "knb-lter-sbc.119.12"))
