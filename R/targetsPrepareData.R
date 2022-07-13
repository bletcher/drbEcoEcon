tar_option_set(packages = c("tidyverse"))

target_prepareData <- 
  tar_plan(
    
    ###########################################
    # find fish caught more than once on a day
    gtOneCapture = dRaw %>%
        filter(tag != "", tag != "ad") %>%
        group_by(date, tag) %>%
        summarise(n = n()) %>%
        filter(n > 1),
    
    #all occasions of fish captured more than once on a day - sent to Fred Henson
    gtOncePerDay = dRaw %>%
        filter(tag %in% gtOneCapture$tag) %>%
        arrange(tag, date) %>%
        select(Water, date, dateTime, tag, Length, Weight, tblSL_SL_ID, Rep, FishNum),
    
    # based on email from Fred on xx, removing the following fish from gtOncePerDay
    # added fish 12 and 17 because they are duplicates on the same day
    toRemove = c(2,4,5,8,9,10,15,18,19,20,21, 12, 17),
    gtOncePerDayToRemove = gtOncePerDay %>% filter(row_number() %in% toRemove),
    
    # remove gtOncePerDayToRemove fish from target_d
    target_d = dRaw %>%
      anti_join(gtOncePerDayToRemove)
    ##########################################
    
    
  )
