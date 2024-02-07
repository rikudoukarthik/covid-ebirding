# Step 2 of subsampling: generate 1000 versions of individual data files 
# subsampled using the subsampled GROUPIDs per location per month per year

randomgroupids<- load(glue("00_data/bird_models/{state_name}/ss01_groupids.RData"))

path_folder <- "00_data/bird_models/{state_name}/ss02_datafiles/"

# creating new directory if it doesn't already exist
if (!dir.exists(path_folder)) {
  dir.create(path_folder, recursive = T)
}

cur_assignment <- 1:1000


for (i in cur_assignment) {
  
  # file names for individual files
  write_path <- glue("{path_folder}data{i}.csv")
  
  
  tictoc::tic(glue("({i}/{max(cur_assignment)}) Filtering data"))
  data_filt = data0_MY_b %>% 
    filter(GROUP.ID %in% randomgroupids[,i]) %>% 
    mutate(MONTHS.TYPE = case_when(MONTH %in% 4:5 ~ "LD",
                                   TRUE ~ "NL"))
  tictoc::toc()
  
  tictoc::tic(glue("({i}/{max(cur_assignment)}) Writing data"))
  write.csv(data_filt, file = write_path, row.names = F)
  tictoc::toc()
  
  
  gc()
  
}


# cleaning up memory
rm(cur_assignment, write_path, data_filt, randomgroupids)

gc()