tar_option_set(packages = c("tidyverse"))

target_globalVariables <- 
  tar_plan(
    target_mainRiver = "West Br Delaware River",
    target_sizeCutoff1 = 230,
    target_sizeCutoff2 = 380
      
  )  
