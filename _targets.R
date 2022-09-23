# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse",
    "lubridate",
    "knitr",
    "targets",
    "tarchetypes",
    "nimble"
  ), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# Load the R scripts with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)
# source("other_functions.R") # Source other scripts as needed. # nolint

list(
  target_globalVariables,
  target_getData,
  target_prepareData,  
  target_getEH,
  target_getEH_main,
  target_getEH_trib,
  tt_trib,
  tt_main,
  tar_render(report, "main.Rmd"),
  tar_render(reportGrowth, "mainGrowth.Rmd")
)
