# Modelling directly with presence-absence data instead of relative abundance.
#
# parallelisation is somehow much less efficient when within a function or loop.

tictoc::tic("Total time elapsed for b_overall_model():")

cur_species_list <- species_list %>% filter(STATE == state_name)

# creating subfolder for model outputs
if (!dir.exists(glue("00_data/bird_models/{state_name}/"))) {
  dir.create(glue("00_data/bird_models/{state_name}/"), recursive = TRUE)
}


###

# subsampling data to model
source("00_scripts/b01_ss_groupids.R")
source("00_scripts/b02_ss_datafiles.R")

# modelling and bootstrapping for results
source("00_scripts/b03_runspeciesmodels.R")
source("00_scripts/b04_resolvebootstrap.R")


###

tictoc::toc()