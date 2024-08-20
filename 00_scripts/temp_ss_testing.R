mt <- "ALL"
load(glue("00_data/bird_models/{state_name}/b01_ss_groupids_{mt}.RData"))

library(bench)

dt_filter <- function(data, ids) {
  
  require(data.table)
  
  setDT(data)
  # setDT(ids)
  
  # Perform the filtering
  result <- data[GROUP.ID %in% ids[, 1]]  
  
  return(result)
}

benchmarked <- mark(
  
  DPLYR = data0_MY_b %>% 
    filter(GROUP.ID %in% randomgroupids[, 1]),
  
  DT = dt_filter(data0_MY_b, randomgroupids),
  
  BASE = data0_MY_b[data0_MY_b$GROUP.ID %in% randomgroupids[, 1]]
  
  iterations = 1
  
)


tictoc::tic("dplyr")
data5 <- data0_MY_b |> 
  filter(GROUP.ID %in% randomgroupids[, 1])
tictoc::toc()

tictoc::tic("dt")
data2 <- setDT(data0_MY_b)
data3 <- data2[GROUP.ID %in% randomgroupids[, 1]]
tictoc::toc()






to_walk <- function(.x, folder, assignment, monthtype) {
  
  # file names for individual files
  write_path <- glue("{folder}data{.x}.RData")
  
  tictoc::tic(glue("{monthtype} ({.x}/{max(assignment)}) Filtering data"))
  data_filt = data0_MY_b[GROUP.ID %in% randomgroupids[, .x]]
  tictoc::toc()
  
  tictoc::tic(glue("{monthtype} ({.x}/{max(assignment)}) Writing data"))
  setDF(data_filt) # converting to data.frame for use further downstream
  save(data_filt, file = write_path)
  tictoc::toc()

}

cur_it <- 1:100

message("Activated future-walking using advanced Kenbunshoku Haki!")
tic(glue("Future-walked over {max(cur_it)} random data file generation"))
plan(multisession, workers = parallel::detectCores()/2)
options(future.globals.maxSize = 10000 * 1024^2) # 10GB
future_walk(cur_it, .progress = TRUE, 
            ~ to_walk(.x, path_folder, cur_it, mt))
plan(sequential)
toc()
gc()
