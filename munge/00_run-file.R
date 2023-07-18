#########################
# Suzanne M Dufault
# Run File
# Origin: 2023-06-16
# Revision: 2023-06-16
#########################

library(tidyverse)
library(here)
library(haven) # for calling in data that isn't in an RData format (e.g., .dta files)
library(labelled) # for working with the labelled structure from the .dta files

source(here("munge", "01_data-cleaning.R")) # cleans REMoxTB data

# Analysis
source(here("munge", "02_run-linear.R"))
source(here("munge", "03_run-complex-biphasic.R"))