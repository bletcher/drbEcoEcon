# R/functions.R
# Put code that want to source into "_targets.R" here.


  #####################################
  ## general functions 
  #####################################
  `%notin%` <- Negate(`%in%`)
  
  
  #####################################
  ## getData functions 
  #####################################
  
  addOccN <- function(d) {
    dates0 <- sort(unique(d$dateYM))
    dates <- data.frame(occ = 1:length(dates0), dateYM = dates0)
    left_join(d, dates)
  }
  
  getKnown <- function(x) {
    firstObs <- min(which(x == 1))
    lastObs <- max(which(x == 1))
    known <- rep(0, length(x))
    known[firstObs:lastObs] <- 1
    if (lastObs != length(known)) {
      known[(lastObs + 1):length(known)] <- NA
    }
    return(known)
  }
  
  addKnownZ2 <- function(d) {
    d %>% 
      group_by(tag) %>%
      arrange(occ) %>%
      mutate(knownZ = getKnown(enc)) %>%
      ungroup() %>%
      arrange(tag, occ)
  }
  
  addFirstLast <- function(d){
    firstLast <- d %>% 
      group_by(tag) %>%
      filter(knownZ == 1) %>%
      summarize(firstObserved = min(occ, na.rm = TRUE),
                lastObserved = max(occ, na.rm = TRUE)) %>%
      ungroup()
    
    left_join(d, firstLast) %>%
      mutate(isFirstObserved = occ == firstObserved,
             isLastObserved = occ == lastObserved)
  }    
  
  addRiverN <- function(d){
    rivers <- unique(d$Water) %>% sort()
    
    d %>% 
      mutate(riverN = as.numeric(factor(d$Water, levels = rivers)))
  }
  
  addMainTrib <- function(d) {
    d %>%
      mutate(
        mainTrib = ifelse(Water == tar_read(target_mainRiver), "main", "trib"),
        mainTribN = ifelse(Water == tar_read(target_mainRiver), 1, 2)
      )
  }
  
  addSizeState <- function(d) {
    d %>%
      mutate(
        sizeState = ifelse(Length < tar_read(target_sizeCutoff1), 1,
                    ifelse(Length < tar_read(target_sizeCutoff2), 2, 3))
      )
  }
  

  combineRiverSizeState <- function(d) {
    d %>%
      mutate(
        state = sizeState + (mainTribN - 1) * 3
      )
  }

  #####################################
  ## Create encounter histories 
  #####################################
  
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
      filter(tag != "", tag != "ad") %>%
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
    
    tags <- encWide %>% dplyr::select(tag)
    
    data <- d %>%
      ehFilter(cols, ops, vals) %>% 
      #filter(ageInSamples > 0, ageInSamples <= maxOccasionValue) %>%
      #filter(ageInSamples %in% 1:maxOccasionValue) %>%
      arrange(tag, date)
    
    first <- apply(eh, 1, function(x) min(which(x != 0)))
    #last <- apply(riverMatrix, 1, function(x) max(which(!is.na(x))))
    #last <- ifelse(last == ncol(riverMatrix), last, last - 1)
    last <- rep(ncol(riverMatrix) - 1, nrow(riverMatrix))
    
    return(list(eh = eh,
                riverMatrix = riverMatrix,
                riverNMatrix = riverNMatrix,
                stateMatrix = stateMatrix,
                tags = tags, 
                first = first, 
                last = last, 
                data = data))
  }
  
  
