# Step 3 of subsampling: run 1000 models for each species and write 1000 CSVs of 
# prediction dfs.

require(tidyverse)
require(lme4)
require(VGAM)
require(parallel)
require(foreach)
require(doParallel)


cur_assignment <- 1:1000


path_folder <- "00_outputs/bird_models/{state_name}/ss03_models/"
datapath_folder <- "00_data/bird_models/{state_name}/ss02_datafiles/"


# creating new directory if it doesn't already exist
if (!dir.exists(path_folder)) {
  dir.create(path_folder, recursive = T)
}


for (k in cur_assignment) {
  
  # file names for individual files
  write_path <- glue("{path_folder}pred{k}.csv")
  data_path <- glue("{datapath_folder}data{k}.csv")

  
  tictoc::tic(glue("Species trends for {state_name}: {k}/{max(cur_assignment)}"))
  
  # read data files
  data_ss = read.csv(data_path) %>% 
    # filter for current month type
    filter(MONTHS.TYPE == cur_m)
  
  
  # start parallel
  n.cores = parallel::detectCores()/2
  # create the cluster
  my.cluster = parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  # register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  trends0 = foreach(i = cur_species_list$COMMON.NAME, 
                    .combine = 'cbind', .errorhandling = 'remove') %dopar%
    singlespeciesmodel(data = data_ss, 
                       species = i, 
                       specieslist = cur_species_list$COMMON.NAME)
  
  parallel::stopCluster(cl = my.cluster)
  
  trends = data.frame(trends0) %>% 
    # converting first row of species names (always true) to column names
    magrittr::set_colnames(.[1,]) %>% 
    slice(-1) %>% 
    # will always have 28 rows
    mutate(timegroupsf = rep(databins$timegroups, 2),
           timegroups = rep(databins$year, 2),
           type = rep(c("freq", "se"), each = 14),
           sl = k) %>%  # sim number
    # pivoting species names longer
    pivot_longer(-c(timegroups, timegroupsf, sl, type), 
                 values_to = "value", names_to = "COMMON.NAME") %>% 
    pivot_wider(names_from = type, values_from = value) %>% 
    # numerical ID for species names, for arranging
    mutate(sp = row_number(), .by = timegroupsf) %>%
    arrange(sl, sp) %>%
    dplyr::select(-sp) %>% 
    # reordering
    relocate(sl, COMMON.NAME, timegroupsf, timegroups, freq, se)
  
  write.csv(trends, file = write_path, row.names = F)
  
  tictoc::toc() 
  
  gc()
  
}