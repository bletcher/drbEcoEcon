tar_option_set(packages = c("tidyverse", "lubridate"))

target_getData <- 
  tar_plan(
    d0 = read.csv("./dataIn/JFIP_ElectrofishingData_Converted.csv"),
    # Simple drake-like syntax:
    #d0 = read_csv(d00),
    dRaw = d0 %>%
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
      ) %>%
      addOccN() %>%
      addKnownZ2() %>%
      addFirstLast() %>%
      addRiverN() %>%
      addMainTrib() %>%
      addSizeState() %>%
      combineRiverSizeState()
      
  )  


