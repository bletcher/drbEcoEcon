tar_option_set(packages = c("tidyverse"))

target_getEH <- 
  tar_plan(
    cols = list(), #list("Water"),
    ops = list(), #list("%in%"),
    vals = list(), #list("West Br Delaware River"),
    #eh <- getEHDataWide(tar_read(d), cols, ops, vals, date, valuesFill = NA)
    target_eh = getEH(target_d, cols, ops, vals)
  )

target_getEH_main <- 
  tar_plan(
    cols_main = list("mainTrib"), #list("Water"),
    ops_main = list("%in%"), #list("%in%"),
    vals_main = list("main"), #list("West Br Delaware River"),
    #eh <- getEHDataWide(tar_read(d), cols, ops, vals, date, valuesFill = NA)
    target_eh_main = getEH(target_d, cols_main, ops_main, vals_main)
  )

target_getEH_trib <- 
  tar_plan(
    cols_trib = list("mainTrib"), #list("Water"),
    ops_trib = list("%in%"), #list("%in%"),
    vals_trib = list("trib"), #list("West Br Delaware River"),
    #eh <- getEHDataWide(tar_read(d), cols, ops, vals, date, valuesFill = NA)
    target_eh_trib = getEH(target_d, cols_trib, ops_trib, vals_trib)
  )
#############################################
## Functions for creating encounter histories 
#############################################

getNeverCaptured <- function(d){
  d %>%
    #filter(ageInSamples > 0 & ageInSamples <= maxOccasionValue) %>%
    #filter(ageInSamples %in% 1:maxOccasionValue) %>%
    group_by(tag) %>%
    summarize(sumEnc = sum(enc, na.rm = TRUE)) %>%
    filter(sumEnc == 0) %>%
    dplyr::select(tag)
}

getFirstObservedOnLastOcc <- function(d) {
  d %>% 
    filter(firstObserved == max(d$occ)) %>%
    select(tag) %>%
    unique()
}

getGT1ObsPerOcc <- function(d) {
  d %>%
    dplyr::group_by(tag, dateYM) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n > 1L)
}

keepOnlyFirstObsPerOcc <- function(d) {
  # for fish that were observed gt 1 time per occasion (dateYM)
  rowsToRemove <- d %>%
    filter(tag %in% getGT1ObsPerOcc(d)$tag) %>%
    group_by(tag, dateYM) %>%
    summarise(dateTime = max(dateTime), n = n()) %>%
    filter(n > 1) %>%
    select(tag, dateTime)
  
  dOut <- d %>%
    anti_join(rowsToRemove)
}


#### getEH functions
op_call <- function(op, lhs, rhs) {
  call(op, sym(lhs), rhs)
}

ehFilter <- function(data, cols, ops, vals) {
  exprs <- purrr::pmap(list(ops, cols, vals), op_call)
  data %>% dplyr::filter(!!!exprs)
}

getEHDataWide <- function(d, cols, ops, vals, var, valuesFill = 0){
  #print(c("in getEHWide",d))
  
  dFiltered <- d %>%
    ehFilter(cols, ops, vals) %>% 
    filter(tag != "", tag != "ad") %>% # redundant
    
    filter(tag != "ad719472") %>%
    
    filter(species == "brown trout") %>%
    mutate(dateC = as.character(Date)) 
  
  dFiltered %>%
    arrange(occ) %>%
    pivot_wider(
      id_cols = tag,
      names_from = dateYM,
      names_prefix = "date_",
      values_from = eval(substitute(var)),
      #values_fill = as.character(valuesFill)
      values_fill = valuesFill
    )
}

getEH <- function(d, cols, ops, vals ){
  
  d <- d %>%
    filter(tag %notin% getNeverCaptured(d)$tag,             # Fish with no observed occasions
           tag %notin% getFirstObservedOnLastOcc(d)$tag     # exclude fish caught on the last capture occasion
    ) %>%
    keepOnlyFirstObsPerOcc()                               # just the first obs on a given day for fish caught gt 1 time per day
  
  encWide <- getEHDataWide(d, cols, ops, vals, "enc", valuesFill = 0)
  eh <- as.matrix(encWide %>% dplyr::select(-tag), nrow = nrow(encWide), ncol = ncol(encWide) - 1)
  
  riverWide <- getEHDataWide(d, cols, ops, vals, "Water", valuesFill = NA)
  riverMatrix <- as.matrix(riverWide %>% dplyr::select(-tag), nrow = nrow(riverWide), ncol = ncol(riverWide) - 1)
  
  riverNWide <- getEHDataWide(d, cols, ops, vals, "riverN", valuesFill = 0)
  riverNMatrix <- as.matrix(riverNWide %>% dplyr::select(-tag), nrow = nrow(riverNWide), ncol = ncol(riverNWide) - 1)
  
  stateWide <- getEHDataWide(d, cols, ops, vals, "state", valuesFill = 0)
  stateMatrix <- as.matrix(stateWide %>% dplyr::select(-tag), nrow = nrow(stateWide), ncol = ncol(stateWide) - 1)
  
  sizeStateWide <- getEHDataWide(d, cols, ops, vals, "sizeState", valuesFill = 0)
  sizeStateMatrix <- as.matrix(sizeStateWide %>% dplyr::select(-tag), nrow = nrow(sizeStateWide), ncol = ncol(sizeStateWide) - 1)
  
  tags <- encWide %>% dplyr::select(tag)
  
  data <- d %>%
    ehFilter(cols, ops, vals) %>% 
    #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue) %>%
    #filter(ageInSamples %in% 1:maxOccasionValue) %>%
    arrange(tag, date)
  
  first <- apply(eh, 1, function(x) min(which(x != 0)))
  #last <- apply(riverMatrix, 1, function(x) max(which(!is.na(x))))
  #last <- ifelse(last == ncol(riverMatrix), last, last - 1)
  
  # set all 'last' to the last occasion - we don't know age
  last <- rep(ncol(riverMatrix) - 0, nrow(riverMatrix))
  
  return(list(eh = eh,
              riverMatrix = riverMatrix,
              riverNMatrix = riverNMatrix,
              stateMatrix = stateMatrix,
              sizeStateMatrix = sizeStateMatrix,
              tags = tags, 
              first = first, 
              last = last, 
              data = data))
}

