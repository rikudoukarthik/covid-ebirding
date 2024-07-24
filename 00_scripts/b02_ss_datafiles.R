require(tidyverse)
require(furrr)
require(parallel)
require(data.table)
require(purrr)

# Step 2 of subsampling: generate 1000 versions of individual data files 
# subsampled using the subsampled GROUPIDs per location per month per year

# convert input data to data.table 
setDT(data0_MY_b)

for (mt in c("LD", "ALL")) {
  
  path_folder <- glue("00_data/bird_models/{state_name}/b02_ss_datafiles_{mt}/")
  
  # creating new directory if it doesn't already exist
  if (!dir.exists(path_folder)) {
    dir.create(path_folder, recursive = TRUE)
  }

  
  load(glue("00_data/bird_models/{state_name}/b01_ss_groupids_{mt}.RData"))
  cur_assignment <- 1:1000
  
  
  walk(cur_assignment, ~ {
    
    # file names for individual files
    write_path <- glue("{path_folder}data{.x}.RData")
    
    tictoc::tic(glue("{mt} ({.x}/{max(cur_assignment)}) Filtering data"))
    data_filt = data0_MY_b[GROUP.ID %in% randomgroupids[, .x]]
    tictoc::toc()
    
    tictoc::tic(glue("{mt} ({.x}/{max(cur_assignment)}) Writing data"))
    setDF(data_filt) # converting to data.frame for use further downstream
    save(data_filt, file = write_path)
    tictoc::toc()
    
  })

  # cleaning up memory
  rm(cur_assignment, write_path, data_filt, randomgroupids)
  gc()
  
}


# convert input data to data.frame for next state iterations
setDF(data0_MY_b)
