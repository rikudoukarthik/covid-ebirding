### data filters -----------------

# data file path
# path to list of group accounts to be filtered out
# path to classification of year-month as COVID categories


dataqual_filt <- function(datapath, groupaccspath, covidclasspath,
                          maxvel = 20, minsut = 2){

  ### importing from usual modified ebd RData
  load(datapath) 
  
  require(tidyverse)
  require(lubridate)

  ### list of group accounts to be filtered
  groupaccs <- read.csv(groupaccspath, 
                        na.strings = c(""," ",NA), quote = "", header = T, 
                        nrows = 401)  # excluding empty cells
  groupaccs <- groupaccs %>% 
    mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                                GA.2 == 1 ~ "GA.2", 
                                TRUE ~ "NG"))
  filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)
  
  ### COVID classification
  covidclass <- read_csv(covidclasspath)
  
  
  
  ### new observer data (to calculate no. of new observers metric) #######

  new_obsr_data <- data %>% 
    select(c("YEAR", "MONTH", "STATE", "SAMPLING.EVENT.IDENTIFIER",
             "LAST.EDITED.DATE", "OBSERVATION.DATE", "OBSERVER.ID")) %>% 
    mutate(LAST.EDITED.DATE = ymd_hms(LAST.EDITED.DATE)) %>% 
    group_by(OBSERVER.ID) %>% 
    arrange(LAST.EDITED.DATE) %>% 
    ungroup() %>% 
    distinct(OBSERVER.ID, .keep_all = TRUE) %>%
    mutate(YEAR = year(LAST.EDITED.DATE),
           MONTH = month(LAST.EDITED.DATE)) %>% 
    filter(YEAR >= 2019) %>% 
    left_join(covidclass) %>% 
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>% 
    rename(LE.YEAR = YEAR,
           LE.MONTH = MONTH) %>% 
    mutate(YEAR = year(OBSERVATION.DATE), 
           MONTH = month(OBSERVATION.DATE))
  
  # filtering
  new_obsr_data <- new_obsr_data %>% anti_join(filtGA) 
  save(new_obsr_data, file = "data/new_obsr_data.RData")
  
  
  
  ### main data filtering ######
  
  data0 <- data %>% 
    filter(YEAR >= 2019) %>% # retaining only data from 2019 onward
    anti_join(filtGA) %>% # removing data from group accounts
    # creating COVID factor
    left_join(covidclass) %>% 
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>% 
    group_by(GROUP.ID) %>% 
    mutate(NO.SP = n_distinct(COMMON.NAME)) %>%
    ungroup() %>% 
    mutate(TIME.D = hour(as_datetime(paste(OBSERVATION.DATE,
                                                  TIME.OBSERVATIONS.STARTED))),
           MIN = minute(as_datetime(paste(OBSERVATION.DATE,
                                          TIME.OBSERVATIONS.STARTED))),
           SPEED = EFFORT.DISTANCE.KM*60/DURATION.MINUTES, # kmph
           SUT = NO.SP*60/DURATION.MINUTES, # species per hour
           # calculate hour checklist ended
           END = floor((TIME.D*60 + MIN + DURATION.MINUTES)/60))
  
  # choose checklists without info on duration with 3 or fewer species
  temp <- data0 %>%
    filter(ALL.SPECIES.REPORTED == 1, PROTOCOL.TYPE != "Incidental") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter(NO.SP <= 3, is.na(DURATION.MINUTES)) %>%
    distinct(GROUP.ID)
  
  # exclude records based on various criteria 
  data0 <- data0 %>%
    mutate(ALL.SPECIES.REPORTED = 
             case_when(ALL.SPECIES.REPORTED == 1 & 
                         (GROUP.ID %in% temp | 
                            SPEED > maxvel |
                            (SUT < minsut & NO.SP <= 3) | 
                            PROTOCOL.TYPE == "Incidental" | 
                            (!is.na(TIME.D) & ((TIME.D <= 4 & END <= 4) | 
                                             (TIME.D >= 20 & END <= 28)
                                           )
                             )
                          ) ~ 0, 
                       ALL.SPECIES.REPORTED == 0 ~ 0,
                       TRUE ~ 1)) %>% 
    select(-SPEED, -SUT, -MIN, -END) %>% 
    filter(ALL.SPECIES.REPORTED == 1)
  
  assign("data0", data0, .GlobalEnv)
         
  save(data0, file = "data/data0.RData")
}


### bootstrapping standard error -----------------

# Useful resources: 
# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/variability-and-uncertainty-standard-deviations-standard-errors-confidence-intervals.html#bootstrap 
# https://websites.pmc.ucsc.edu/~mclapham/Rtips/resampling

# For some metrics like birding distance, it is possible to calculate SE from the data itself ("expected SE" = SD of sample / root sample size). But when there is a possibility to calculate the empirical SE instead (via bootstrapping) which is more accurate, there is no point in going for the former. 

boot_se = function(x, fn = mean, B = 1000) {
  1:B %>%
    # For each iteration, generate a sample of x with replacement
    map(~ x[sample(1:length(x), replace = TRUE)]) %>%
    # Obtain the fn estimate for each bootstrap sample
    map_dbl(fn) %>%
    # Obtain the standard error
    sd()
}

