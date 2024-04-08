# Step 3 of subsampling: run 1000 models for each species and write 1000 CSVs of 
# prediction dfs.

require(tidyverse)
require(lme4)
require(VGAM)
require(parallel)
require(foreach)
require(doParallel)


for (mt in c("LD", "NL")) {
  
  cur_assignment <- 1:1000
  
  
  path_folder <- "00_outputs/bird_models/{state_name}/ss03_models_{mt}/"
  datapath_folder <- "00_data/bird_models/{state_name}/ss02_datafiles_{mt}/"
  
  
  # creating new directory if it doesn't already exist
  if (!dir.exists(path_folder)) {
    dir.create(path_folder, recursive = T)
  }
  
  
  for (k in cur_assignment) {
    
    # file names for individual files
    write_path <- glue("{path_folder}pred{k}.csv")
    data_path <- glue("{datapath_folder}data{k}.csv")
    
    
    tictoc::tic(glue("Species trends for {state_name}: {mt} {k}/{max(cur_assignment)}"))
    
    # read data files
    data_ss = read.csv(data_path)
    
    
    # start parallel
    n.cores = parallel::detectCores()/2
    # create the cluster
    my.cluster = parallel::makeCluster(
      n.cores, 
      type = "PSOCK"
    )
    # register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)
    
    birds_mod = foreach(i = cur_species_list$COMMON.NAME, 
                          .combine = 'cbind', .errorhandling = 'remove') %dopar%
      singlespeciesmodel(data = data_ss, 
                         species = i, 
                         specieslist = cur_species_list$COMMON.NAME)
    
    parallel::stopCluster(cl = my.cluster)
    

    write.csv(birds_mod, file = write_path, row.names = F)
    
    tictoc::toc() 
    
    gc()
    
  }
  
}