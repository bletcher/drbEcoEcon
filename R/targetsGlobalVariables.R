tar_option_set(packages = c("tidyverse"))

target_globalVariables <- 
  tar_plan(
    target_mainRiver = "West Br Delaware River",
    target_sizeCutoff1 = 230,
    target_sizeCutoff2 = 380,
    target_mainOrTrib = "main" #"main" or "trib"
    
    # #not using until need sep cuttoffs for each location
    # target_sizeCutoffs = tribble(
    #   ~mainTrib, ~ cuttoff,
    #   "main", 230,
    #   "main", 380,
    #   "trib", 230,
    #   "trib", 380
    #   
    # )
      
  )  
