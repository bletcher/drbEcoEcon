tar_option_set(packages = c("tidyverse"))

getInits <- function(y) {
  zInits <- y + 1 # non-detection -> alive
  zInits[zInits == 2] = 1 # dead -> alive
  return(zInits)
}

target_getEH_InputData <- 
  tar_plan(
    target_y = tar_read(target_eh)$eh,
    target_state = tar_read(target_eh)$state,
    target_stateMatrix = tar_read(target_eh)$stateMatrix,
    target_first = tar_read(target_eh)$first, 
    target_last = tar_read(target_eh)$last, # check - this should be num Occ for all fish
    target_zInits = getInits(target_y) # non-detection -> alive
  )
