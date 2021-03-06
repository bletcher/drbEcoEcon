tar_option_set(packages = c("tidyverse", "lubridate"))

#target_readData <- tar_file(d0, "./dataIn/JFIP_ElectrofishingData_Converted.csv", format = "file")

target_getData <-
  tar_plan(
    d0 = read.csv("./dataIn/JFIP_ElectrofishingData_Converted.csv"),
   
    dRaw0 = d0 %>%
      select(Longitude, Latitude, Date, Time, SiteLen, Discharge, WaterTemp, GearType, 
             TimeStart, Time.Stop, TimeFish, FishNum, SpeciesCd, Length, Weight, 
             Wild, Stage, MarkAppCd, Mark.PreCd, Water, tblSL_SL_ID, Rep, FishNum) %>%
      mutate(
        date = dmy(Date),
        #year = year(Date),
        #month = month(Date),
        dateC = as.character(Date),
        dateYM = substr(date, 0, 7),
        time = hms(Time),
        dateTime = date + time,
        species = recode(SpeciesCd, "326" = "rainbow trout", "328" = "brown trout", "329" = "brook trout", .default=NA_character_),
        tag = ifelse(MarkAppCd == "", Mark.PreCd, MarkAppCd),
        tag = tolower(tag),
        enc = 1
      ),
    
    dRaw = dRaw0 %>%
      removeUntagged() %>%
      removeUncaptured() %>%
      addOccN() %>%
      addKnownZ2() %>%
      addFirstLast() %>%
      addRiverN() %>%
      addMainTrib(target_mainRiver) %>%
      addSizeState(target_sizeCutoff1, target_sizeCutoff2) %>%
      combineRiverSizeState() %>%
      doNothing()
      
  )  


#####################################
## getData functions 
#####################################
doNothing <- function(d) {
  d
}
  
removeUntagged <- function(d) {
  d %>%
    filter(tag != "") %>%
    filter(tag != "ad")
}

# Fish with no captures - NONE
removeUncaptured <- function(d) {
  unCap <- d %>%
    group_by(date, tag) %>%
    summarise(n = n()) %>%
    filter(n == 0)
  d %>% anti_join(unCap)
}

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

addMainTrib <- function(d, r) {
  d %>%
    mutate(
      mainTrib = ifelse(Water == r, "main", "trib"),
      mainTribN = ifelse(Water == r, 1, 2)
    )
}

addSizeState <- function(d, s1, s2) {
  d %>%
    mutate(
      sizeState = ifelse(Length < s1, 1,
                         ifelse(Length < s2, 2, 3))
    )
}


combineRiverSizeState <- function(d) {
  d %>%
    mutate(
      state = sizeState + (mainTribN - 1) * 3
    )
}

