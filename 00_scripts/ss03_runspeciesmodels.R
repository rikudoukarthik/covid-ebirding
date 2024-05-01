# Step 3 of subsampling: run 1000 models for each species and write 1000 CSVs of 
# prediction dfs.

require(tidyverse)
require(lme4)
require(VGAM)
require(parallel)
require(foreach)
require(doParallel)

source("00_scripts/sstemp_functions.R")


for (mt in c("LD", "NL")) {
  
  cur_assignment <- 1:1000
  
  
  path_folder <- glue("00_outputs/bird_models/{state_name}/")
  datapath_folder <- glue("00_data/bird_models/{state_name}/ss02_datafiles_{mt}/")
  
  
  # creating new directory if it doesn't already exist
  if (!dir.exists(path_folder)) {
    dir.create(path_folder, recursive = T)
  }
  
  write_path <- glue("{path_folder}ss03_models_{mt}.csv")
  
  
  for (k in cur_assignment) {
    
    # file names for individual files
    data_path <- glue("{datapath_folder}data{k}.RData")
    
    
    tictoc::tic(glue("Species trends for {state_name}: {mt} {k}/{max(cur_assignment)}"))
    
    # read data files
    load(data_path)
    
    
    # start parallel
    n.cores = parallel::detectCores()/2
    # create the cluster
    my.cluster = parallel::makeCluster(
      n.cores, 
      type = "PSOCK"
    )
    # register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)
    
    birds_mod0 = foreach(i = cur_species_list$COMMON.NAME, 
                          .combine = "bind_rows", .errorhandling = "remove") %dopar%
      singlespeciesmodel(data = data_filt, 
                         species = i, 
                         specieslist = cur_species_list,
                         iter = k)
    
    parallel::stopCluster(cl = my.cluster)

    
    # storing as single object
    if (!exists("birds_mod", envir = .GlobalEnv)) {
      birds_mod <- birds_mod0
    } else {
      birds_mod <- birds_mod %>% 
        bind_rows(birds_mod0)
    }
    assign("birds_mod", birds_mod, envir = .GlobalEnv)
    
    tictoc::toc() 
    
    gc()
    
  }
  
  write.csv(birds_mod, file = write_path, row.names = F)
  
}