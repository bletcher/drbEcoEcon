# can move this section to ttt


tar_option_set(packages = c("tidyverse"))

getInits <- function(y) {
  zInits <- y + 1 # non-detection -> alive
  zInits[zInits == 2] = 1 # dead -> alive
  return(zInits)
}

target_getEH_InputData <- 
  tar_plan(
    #target_y = tar_read(target_eh)$eh,
    target_state = target_eh$state,
    target_stateMatrix = target_eh$stateMatrix,
    target_first = target_eh$first, 
    target_last = target_eh$last, # check - this should be num Occ for all fish
    target_zInits = getInits(target_eh$eh) # non-detection -> alive
  )
