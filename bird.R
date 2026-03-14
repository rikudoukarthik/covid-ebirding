# pipeline for bird reporting analysis

# global setup
source("00_scripts/b01_setup.R")


# data import & prep (one-time) -------------------------------------------

# only needs to be run once, later analyses use saved output data/files
source("00_scripts/b02_import_one_time.R") 


# analyses ----------------------------------------------------------------

# full country analysis (4 steps)
source("00_scripts/b03_overall_00_pipeline.R")

# prolific birder analysis 
source("00_scripts/b03_prolific.R")


# plotting ----------------------------------------------------------------


