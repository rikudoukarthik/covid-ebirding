# Step 3 of subsampling: run 1000 models for each species and write 1000 CSVs of 
# prediction dfs.

require(tidyverse)
require(lme4)
require(VGAM)
require(parallel)
require(foreach)
require(doParallel)

source("00_scripts/b_model_functions.R")


for (mt in c("LD", "ALL")) {
  
  cur_assignment <- 2:500
  
  
  path_folder <- glue("00_outputs/bird_models/{state_name}/b03_models_{mt}/")
  datapath_folder <- glue("00_data/bird_models/{state_name}/b02_ss_datafiles_{mt}/")
  
  
  # creating new directory if it doesn't already exist
  if (!dir.exists(path_folder)) {
    dir.create(path_folder, recursive = T)
  }
  
  
  for (k in cur_assignment) {
    
    # file names for individual files  
    write_path <- glue("{path_folder}models{k}.csv")
    data_path <- glue("{datapath_folder}data{k}.RData")
    
    
    tictoc::tic(glue("Species trends for {state_name}: {mt} {k}/{max(cur_assignment)}"))
    
    # read data files
    load(data_path)
    
    
    # start parallel
    n.cores = if (state_name == "India" & mt == "ALL") {
      parallel::detectCores()/3
    } else {
      parallel::detectCores()/2
    }
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

    
    write.csv(birds_mod0, file = write_path, row.names = F)
    
    tictoc::toc() 

  }
  
  
  # need to combine all individual model results into a single CSV (as done for Assam)
  
  path_folder2 <- glue("00_outputs/bird_models/{state_name}/")
  write_path2 <- glue("{path_folder2}b03_models_{mt}.csv")
  
  if (dir.exists(path_folder)) {
    
    list_files <- list.files(path_folder)
    birds_mod <- map(list_files, ~ read.csv(glue("{path_folder}{.x}"))) %>% 
      list_rbind()
    
    write.csv(birds_mod, file = write_path2, row.names = FALSE)
    
  }
  
  
}
