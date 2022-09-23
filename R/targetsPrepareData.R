tar_option_set(packages = c("tidyverse"))

target_prepareData <- 
  tar_plan(
    
    dRaw0 = d0 %>%
      add_row(d0_2021) %>%
      rename(MarkPreCd = `Mark PreCd`) %>%
      select(Longitude, Latitude, Date, Time, SiteLen, Discharge, WaterTemp, GearType,
             TimeFish, FishNum, SpeciesCd, Length, Weight,
             Wild, Stage, MarkAppCd, MarkPreCd, Water, tblSL_SL_ID, Rep, FishNum) %>%
      mutate(
        date = dmy(Date),
        #year = year(Date),
        #month = month(Date),
        dateC = as.character(Date),
        dateYM = substr(date, 0, 7),
        time = hms(Time),
        dateTime = date + time,
        species = recode(SpeciesCd, "326" = "rainbow trout", "328" = "brown trout", "329" = "brook trout", .default=NA_character_),
        tag = ifelse(is.na(MarkAppCd), MarkPreCd, MarkAppCd),
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
      doNothing(),
    
    ###########################################
    # find fish caught more than once on a day
    gtOneCapture = dRaw %>%
        filter(tag != "", tag != "ad") %>%
        group_by(date, tag) %>%
        summarise(n = n()) %>%
        filter(n > 1),
    
    # #all occasions of fish captured more than once on a day - sent to Fred Henson
    # gtOncePerDay = dRaw %>%
    #     filter(tag %in% gtOneCapture$tag) %>%
    #     arrange(tag, date) %>%
    #     select(Water, date, dateTime, tag, Length, Weight, tblSL_SL_ID, Rep, FishNum),
    # 
    # # based on email from Fred on xx, removing the following fish from gtOncePerDay
    # # added fish 12 and 17 because they are duplicates on the same day
    # toRemove = c(2,4,5,8,9,10,15,18,19,20,21, 12, 17),
    # gtOncePerDayToRemove = gtOncePerDay %>% filter(row_number() %in% toRemove),
    # 
    
    # all tags that were ever observed > once per day
    gtOncePerDay_tags = dRaw %>%
      filter(!is.na(date)) %>%
      filter(tag %in% gtOneCapture$tag, date %in% gtOneCapture$date) %>%
      distinct(tag) %>%
      select(tag),
    
    # remove all fish that were ever observed > once per day
    target_d = dRaw %>%
      # remove gtOncePerDayToRemove fish from target_dRaw
      # anti_join(gtOncePerDayToRemove)
      anti_join(gtOncePerDay_tags)
    ##########################################
    
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

