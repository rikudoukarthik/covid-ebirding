# Step 1 of subsampling: generate 1000 versions of checklist-location-month-year 
# combinations (after data is filtered, cleaned, etc., which is already done by now)


# create the set of random locations (doesn't work inside a function)

require(tidyverse)
require(parallel)
require(foreach)
require(doParallel)


# doesn't return locations but rather their group IDs
createrandomlocs <- function(locs) {
  
  require(tidyverse)
  
  locs1 = locs %>% 
    group_by(LOCALITY.ID, MONTH, M.YEAR) %>% 
    # sampling 1 random list from the location for that month and year
    sample_n(1)
  
  return(locs1$GROUP.ID)
  
}


data0_MY_b <- data0_MY_b %>% 
  mutate(MONTHS.TYPE = case_when(MONTH %in% 4:5 ~ "LD",
                                 TRUE ~ "NL"))


# filtering separately for lockdown and non-lockdown months
for (mt in c("LD", "NL")) {
  
  data_locs <- data0_MY_b %>% 
    filter(ALL.SPECIES.REPORTED == 1,
           STATE == state_name,
           MONTHS.TYPE == mt) %>%
    distinct(LOCALITY.ID, GROUP.ID, MONTH, M.YEAR)
  
  
  # start parallel
  n.cores = parallel::detectCores()/2
  # create the cluster
  my.cluster = parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  # register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  
  randomgroupids = foreach(i = 1:1000, .combine = 'cbind') %dopar%
    createrandomlocs(data_locs)
  
  parallel::stopCluster(cl = my.cluster)
  
  save(randomgroupids, 
       file = glue("00_data/bird_models/{state_name}/ss01_groupids_{mt}.RData"))
  
  
  # cleaning up memory
  rm(data_locs, n.cores, my.cluster, randomgroupids)
  
  detach("package:doParallel", unload = TRUE)
  detach("package:foreach", unload = TRUE)
  detach("package:parallel", unload = TRUE)
  
  gc()
  
}