require(glue)
require(purrr)

# this analysis script need only be run once ###
# individual sections can be run independently, e.g., if subsampling has been
# run once already, in the next session modelling can be run directly


# subsampling data to model (parallel compute & main data object i --------

# importing main bird data
load("00_data/data0_MY_b.RData")
data0_MY_b <- bind_rows("LD" = data0_MY_b %>% filter(MONTH %in% 4:5), 
                        "ALL" = data0_MY_b, 
                        .id = "MONTHS.TYPE")

# iterate this step over all spatial units
walk(spat_iter, ~ {
  
  state_name <- .x
  cur_species_list <- bird_anal_spec_list() %>% filter(STATE == state_name)
  
  # creating subfolder for model outputs
  if (!dir.exists(glue("00_data/bird_models/{state_name}/"))) {
    dir.create(glue("00_data/bird_models/{state_name}/"), recursive = TRUE)
  }
  
  tictoc::tic(glue("Subsampling data (both steps) for {state_name} completed in:"))
  source("00_scripts/b03_overall_01_ss_groupids.R")
  source("00_scripts/b03_overall_02_ss_datafiles.R")
  tictoc::toc()
  
})

# removing main bird data
rm(data0_MY_b)


# modelling bird reporting ------------------------------------------------

# modelling directly with presence-absence data instead of relative abundance
# parallelisation is somehow much less efficient when within a function or loop

walk(spat_iter, ~ {
  
  state_name <- .x
  cur_species_list <- bird_anal_spec_list() %>% filter(STATE == state_name)

  tictoc::tic(glue("Modelling for {state_name} completed in:"))
  source("00_scripts/b03_overall_03_runspeciesmodels.R")
  tictoc::toc()
  
})


# bootstrapping for results -----------------------------------------------

walk(spat_iter, ~ {
  
  state_name <- .x
  cur_species_list <- bird_anal_spec_list() %>% filter(STATE == state_name)

  tictoc::tic(glue("Resolving & bootstrapping for {state_name} completed in:"))
  source("00_scripts/b03_overall_04_resolvebootstrap.R")
  tictoc::toc()
  
})

