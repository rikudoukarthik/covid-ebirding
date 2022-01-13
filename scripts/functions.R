### data filters -----------------

# data object
# list of group accounts to be filtered out


dataqual_filt <- function(data, groupaccspath){

  require(tidyverse)
  require(lubridate)

  groupaccs <- read.csv(groupaccspath, 
                        na.strings = c(""," ",NA), quote = "", header = T, 
                        nrows = 401)  # excluding empty cells
  groupaccs <- groupaccs %>% 
    mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                                GA.2 == 1 ~ "GA.2", 
                                TRUE ~ "NG"))
  filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)
  
  
  
  ### new observer data (to calculate no. of new observers metric)
  new_obsr_data <- data %>% 
    select(c("YEAR", "MONTH", "STATE", "SAMPLING.EVENT.IDENTIFIER",
             "LAST.EDITED.DATE", "OBSERVER.ID")) %>% 
    mutate(LAST.EDITED.DATE = ymd_hms(LAST.EDITED.DATE)) %>% 
    group_by(OBSERVER.ID) %>% 
    arrange(LAST.EDITED.DATE) %>% 
    ungroup() %>% 
    distinct(OBSERVER.ID, .keep_all = TRUE) %>%
    mutate(LE.YEAR = year(LAST.EDITED.DATE),
           LE.MONTH = month(LAST.EDITED.DATE)) %>% 
    filter(LE.YEAR >= 2019) %>% 
    mutate(COVID = factor(case_when(LE.YEAR == 2019 ~ "BEF_19", 
                                    LE.YEAR == 2020 ~ "DUR_20",
                                    LE.YEAR == 2021 & LE.MONTH < 9 ~ "DUR_21",
                                    LE.YEAR == 2021 & LE.MONTH >= 9 ~ "AFT_21"),
                          levels = c("BEF_19","DUR_20","DUR_21","AFT_21")))
  
  # filtering
  new_obsr_data <- new_obsr_data %>% anti_join(filtGA) 
  
  save(new_obsr_data, file = "data/new_obsr_data.RData")
  
  
  
  ### main data
  data <- data %>% filter(YEAR >= 2019) 
  save(data, file = "data/rawdata.RData")
  
  # filters
  data0 <- data %>% anti_join(filtGA) %>% 
    filter(PROTOCOL.TYPE %in% c("Traveling","Stationary"))
  
  rm(list = setdiff(ls(envir = .GlobalEnv), c("data0")), pos = ".GlobalEnv")
  save.image("data/data0.RData")
}