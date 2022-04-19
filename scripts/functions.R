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
  require(terra)
  require(sp)
  
  
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
    mutate(HOUR = hour(as_datetime(paste(OBSERVATION.DATE,
                                         TIME.OBSERVATIONS.STARTED))),
           MIN = minute(as_datetime(paste(OBSERVATION.DATE,
                                          TIME.OBSERVATIONS.STARTED))),
           SPEED = EFFORT.DISTANCE.KM*60/DURATION.MINUTES, # kmph
           SUT = NO.SP*60/DURATION.MINUTES, # species per hour
           # calculate hour checklist ended
           HOUR.END = floor((HOUR*60 + MIN + DURATION.MINUTES)/60))
  
  
  load("data/data0.RData")

  
  
  ### exclude records based on various criteria 
  
  # choose complete checklists without info on duration with 3 or fewer species
  temp1 <- data0 %>%
    filter(ALL.SPECIES.REPORTED == 1, PROTOCOL.TYPE != "Incidental") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter(NO.SP <= 3, is.na(DURATION.MINUTES)) %>%
    distinct(GROUP.ID)
  
  # getting list of GROUP.IDs inside IN boundary for pelagic filter
  temp2 <- data0 %>% 
    ungroup() %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    terra::vect(geom = c("LONGITUDE","LATITUDE"), crs = crs(india)) %>% 
    terra::intersect(india) %>% 
    terra::as.data.frame() %>% 
    select(GROUP.ID) # these are lists inside IN so need to be excluded (!) in filter step
  
  
  # true completeness + other filters
  data0 <- data0 %>%
    mutate(ALL.SPECIES.REPORTED = case_when(ALL.SPECIES.REPORTED == 1 & 
                                              (GROUP.ID %in% temp1 | 
                                                 SPEED > maxvel |
                                                 (SUT < minsut & NO.SP <= 3) | 
                                                 PROTOCOL.TYPE == "Incidental") ~ 0, 
                                            ALL.SPECIES.REPORTED == 0 ~ 0,
                                            TRUE ~ 1),
           OTHER.FILTERS = case_when(ALL.SPECIES.REPORTED == 1 & 
                                       (
                                         # nocturnal filter
                                         (!is.na(HOUR) & ((HOUR <= 4 & HOUR.END <= 4) | 
                                                             (HOUR >= 20 & HOUR.END <= 28)) ) |
                                           # pelagic filter
                                           !(GROUP.ID %in% temp2) |
                                           # distance filter
                                           (EFFORT.DISTANCE.KM > 50)
                                       ) ~ 0, 
                                     ALL.SPECIES.REPORTED == 0 ~ 0,
                                     TRUE ~ 1)) %>% 
    select(-SPEED, -SUT, -MIN, -HOUR.END) %>% 
    filter(ALL.SPECIES.REPORTED == 1 & OTHER.FILTERS == 1)
  
  assign("data0", data0, .GlobalEnv)
  
  save(data0, file = "data/data0.RData")
  
}


### bootstrapping standard error -----------------

# Useful resources: 
# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/variability-and-uncertainty-standard-deviations-standard-errors-confidence-intervals.html#bootstrap 
# https://websites.pmc.ucsc.edu/~mclapham/Rtips/resampling

# For some metrics like birding distance, it is possible to calculate SE from the data itself 
# ("expected SE" = SD of sample / root sample size). But when there is a possibility to 
# calculate the empirical SE instead (via bootstrapping) which is more accurate, there is 
# no point in going for the former. 

# here we are calculating SE but what we want is CIs, and we use +-1.96*SE to calculate them.
# this assumes Gaussian distribution about the mean/median (symmetric CIs). we can make this 
# assumption only because sample size is sufficiently large. 
# when we need to propagate from state to nation level, this assumption is necessary.
# thus, for all the initial metrics, this is okay.

boot_se = function(x, fn = mean, B = 1000) {
  
  require(tidyverse)
  require(rlang)
  
  a <- 1:B %>%
    # For each iteration, generate a sample of x with replacement
    map(~ x[sample(1:length(x), replace = TRUE)]) %>%
    # Obtain the fn estimate for each bootstrap sample
    map_dbl(fn)
  
  colnm <- paste0(unlist(str_split(deparse(substitute(x)), "\\$"))[2],
                  ".",
                  str_to_upper(deparse(substitute(fn))))
  
  return(data %>%
           dplyr::summarise(SE = sd(a),
                            !!colnm := fn(a)))
}

### bootstrapping confidence intervals -----------------

# when there is no issue of propagating SEs, always better to bootstrap CIs directly,
# without making the Gaussian assumption. this allows asymmetric CIs which are more likely
# with counts, proportions, etc.

boot_ci = function(data, x, fn = mean, B = 1000) {
  
  require(tidyverse)
  require(rlang)
  
  a <- 1:B %>%
    # For each iteration, generate a sample of x with replacement
    map(~ x[sample(1:length(x), replace = TRUE)]) %>%
    # Obtain the fn estimate for each bootstrap sample
    map_dbl(fn) 
  
  colnm <- paste0(unlist(str_split(deparse(substitute(x)), "\\$"))[2],
                  ".",
                  str_to_upper(deparse(substitute(fn))))

  return(data %>%
           dplyr::summarise(CI.L = stats::quantile(a, 0.025), # Obtain the CIs
                            CI.U = stats::quantile(a, 0.975),
                            !!colnm := fn(a)))
}


### raster aggregation -----------------

# Function to supply in raster::aggregate(), with 25% threshold for classification
# (works for any binary classification to be aggregated; here it is urban/non-)
# https://stackoverflow.com/questions/23369579/aggregate-raster-in-r-with-na-values 

rast_agg_fn <- function(x, ...){
  if_else(sum(x %in% 1) >= 0.25*length(x), 1, 0)
}

# x = c(rep(1, 2), rep(0, 8), rep(1, 10), rep(0, 30)) # for testing


### raster calc: log proportional change in number of lists -----------------

rast_logpropchange <- function(x, y, k = 1)  {
  
  # x is first, y is second, so change from x to y
  
  # # when no. of lists has been log transformed
  # # log(n2) - log(n1) = log(n2/n1)
  # return(y-x)
  
  # when value has not been transformed to prevent NA
  x <- x + k
  y <- y + k
  
  # proportion of urban birding in the cell changed by y/x [0,+n] times.
  return(log10(y/x))
  
}

### raster calc: proportional change in number of lists -----------------

rast_propchange <- function(x, y, k = 1)  {
  
  # x is first, y is second, so change from x to y
  
  # adding 1 to all rows to prevent NA (from zeroes)
  x <- x + k
  y <- y + k
  
  # proportion of urban birding in the cell changed by y/x [0,+n] times.
  return(y/x)
  
}
