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


data0_MY_b <- bind_rows("LD" = data0_MY_b %>% filter(MONTH %in% 4:5), 
                        "ALL" = data0_MY_b, 
                        .id = "MONTHS.TYPE")


# creating new directory if it doesn't already exist
if (!dir.exists(glue("00_data/bird_models/{state_name}/"))) {
  dir.create(glue("00_data/bird_models/{state_name}/"), recursive = TRUE)
}


# filtering separately for lockdown and all months
for (mt in c("LD", "ALL")) {
  
  require(parallel)
  require(foreach)
  require(doParallel)

  data_locs <- data0_MY_b %>% 
    filter(ALL.SPECIES.REPORTED == 1,
           MONTHS.TYPE == mt) %>%
    {if (state_name == "India") {
      .
    } else {
      filter(., STATE == state_name)
    }} |> 
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
       file = glue("00_data/bird_models/{state_name}/b01_ss_groupids_{mt}.RData"))
  
  
  # cleaning up memory
  rm(data_locs, n.cores, my.cluster, randomgroupids)
  
  detach("package:doParallel", unload = TRUE)
  detach("package:foreach", unload = TRUE)
  detach("package:parallel", unload = TRUE)
  
}
