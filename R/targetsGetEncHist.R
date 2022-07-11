tar_option_set(packages = c("tidyverse"))

target_getEH <- 
  tar_plan(
    cols = list(), #list("Water"),
    ops = list(), #list("%in%"),
    vals = list(), #list("West Br Delaware River"),
    #eh <- getEHDataWide(tar_read(d), cols, ops, vals, date, valuesFill = NA)
    target_eh = getEH(tar_read(target_d), cols, ops, vals)
  )
