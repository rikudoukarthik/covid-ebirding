require(tidyverse)
require(lubridate)
require(glue)
require(spdep)
require(sf)

load("00_data/maps_sf.RData")
load("00_data/grids_st_sf.RData")
# removing unnecessary objects from maps_sf.RData
rm(g2_in_sf, g3_in_sf, g4_in_sf,
   g2_st_sf, g3_st_sf, g4_st_sf)


# one time imports --------------------------------------------------------

# Only needed once in the beginning to obtain the data.

# Getting MODIS data ###

if (!file.exists("00_data/in_LULC_MODIS/in_LULC_MODIS.tif") & 
    !file.exists("00_data/rast_UNU.RData")) {
  
  getmodisdata()
  
} else {
  print("MODIS data is ready to use!")
}


# #  Getting EBD data + filtering and preparing for current use (62 mins on server) ###
# tictoc::tic("Preparing data for analyses")
# rawdatapath <- "00_data/ebd_IN_prv_relMay-2022.txt"
# senspath <- "00_data/ebd_sensitive_relMay-2022_IN.txt"
# datapath <- "00_data/ebd_IN_relMay-2022.RData"
# groupaccspath <- "00_data/ebd_users_GA_relMay-2022.csv"
# covidclasspath <- "00_data/covid_classification.csv"
# rast_UNU_path <- "00_data/rast_UNU.RData"
# 
# data_qualfilt_prep(rawdatapath, senspath, datapath, groupaccspath, covidclasspath,
#                    rast_UNU_path,
#                    maxvel = 20, minsut = 2)
# tictoc::toc()


# processing for prolific analysis ----------------------------------------

# load required data
load("00_data/data0_MY_b.RData")
load("00_data/data0_MY_b_slice.RData")

# timeline metadata
timeline <- data0_MY_b_slice_G %>% distinct(COVID, YEAR, MONTH, M.YEAR)

# species list (all species included in list for other analyses, across all states)
species_prolific <- bird_anal_spec_list() %>% distinct(COMMON.NAME, SP.CATEGORY)


## highly birded personal locations -----------------------------------------------

# eBird user info
eBird_users <- read.delim("00_data/ebd_users_relMay-2022.txt", 
                          sep = "\t", header = T, quote = "", 
                          stringsAsFactors = F, na.strings = c(""," ",NA)) %>% 
  transmute(OBSERVER.ID = observer_id,
            FULL.NAME = paste(first_name, last_name, sep = " "))


# getting overall checklist leaders from data
prolific_a <- data0_MY_b_slice_S %>% 
  group_by(OBSERVER.ID, LOCALITY.ID, LOCALITY) %>% 
  dplyr::summarise(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  ungroup() %>% 
  filter(TOT.LISTS >= 800) %>% # 200 per year threshold (arbitrary)
  arrange(desc(TOT.LISTS)) %>% 
  left_join(eBird_users) %>% 
  dplyr::select(-TOT.LISTS)


# getting checklist leaders with equal effort across years from data
prolific_b <- data0_MY_b_slice_S %>% 
  group_by(OBSERVER.ID, LOCALITY.ID, LOCALITY, M.YEAR, MONTH) %>% 
  dplyr::summarise(MONTH.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  filter(MONTH.LISTS >= 15) %>% # 200 per year roughly translates to 16 per month
  dplyr::summarise(N = n_distinct(MONTH),
                   TOT.LISTS = sum(MONTH.LISTS)) %>% 
  filter(N >= 10) %>% # at least 10 months in a year
  dplyr::summarise(N = n_distinct(M.YEAR),
                   TOT.LISTS = sum(TOT.LISTS)) %>% 
  filter(N == 4) %>% # all four years
  ungroup() %>% 
  arrange(desc(TOT.LISTS)) %>% 
  left_join(eBird_users) %>% 
  dplyr::select(-N, -TOT.LISTS)

prolific <- inner_join(prolific_a, prolific_b)

rm(eBird_users, prolific_a, prolific_b, data0_MY_b_slice_G, data0_MY_b_slice_S)


## preparing data for models -----------------------------------------------

# downsizing main data to prolific
prol_data <- data0_MY_b |> 
  right_join(prolific, by = c("OBSERVER.ID", "LOCALITY.ID", "LOCALITY"))
rm(data0_MY_b)
gc()


# to join for presence-absence of various species
prol_list_spec_presabs <- prol_data %>% 
  group_by(SAMPLING.EVENT.IDENTIFIER, COMMON.NAME) %>% 
  reframe(OBSERVATION.COUNT = max(OBSERVATION.COUNT))

# to later join checklist metadata
prol_list_metadata <- prol_data |> 
  distinct(OBSERVER.ID, LOCALITY.ID, SAMPLING.EVENT.IDENTIFIER, NO.SP, M.YEAR, MONTH)

# join patches with patch-specific species
prol_patch_spec <- prol_data %>% 
  group_by(LOCALITY.ID) %>% 
  filter(COMMON.NAME %in% species_prolific$COMMON.NAME) %>% 
  distinct(LOCALITY.ID, COMMON.NAME)


prol_data <- prol_list_metadata %>% 
  # for every list from a person's patch, adding all possible species for the patch
  # many-to-many because we are expanding by species (y) AND checklists (x)
  left_join(prol_patch_spec,, by = "LOCALITY.ID",
            relationship = "many-to-many") %>% 
  # joining observation count to species, and will replace NAs with 0
  left_join(prol_list_spec_presabs, 
            by = c("COMMON.NAME", "SAMPLING.EVENT.IDENTIFIER")) %>% 
  mutate(REPORT = replace_na(OBSERVATION.COUNT, "0"),
         OBSERVATION.COUNT = NULL) %>% 
  # converting counts to presences
  mutate(REPORT = if_else(REPORT == "0", 0, 1))


## filtering data ----------------------------------------------------------

# - calculating overall (12 months, 4 years) repfreq for each patch-species combo to 
#   filter out very uncommon species based on threshold.
# - also removing species only "present" (same threshold) in one locality because 
#   cannot use random effect with only one level in model.
#
# - finally, removing locations where after above filters, both sp. categ. do not exist

prol_list_tot <- prol_list_metadata %>% 
  group_by(OBSERVER.ID, LOCALITY.ID) %>% 
  reframe(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) 

filt <- prol_data %>% 
  # prol_data is list-level, so some patch-species combos will be 0 in certain lists
  # interested in species repfreq for each patch, so here removing 0 to count presences
  filter(REPORT != 0) %>% 
  group_by(OBSERVER.ID, LOCALITY.ID, COMMON.NAME) %>% 
  reframe(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  right_join(prol_list_tot, by = c("OBSERVER.ID", "LOCALITY.ID")) %>% 
  mutate(PROP.LISTS = 100*NO.LISTS/TOT.LISTS) %>% 
  # filtering for species present in threshold lists 
  # (RoPi with 11% not converging in max iterations so setting threshold at 15, not 10)
  filter(PROP.LISTS >= 15) %>%
  # removing species present only in one locality
  group_by(COMMON.NAME) %>% 
  mutate(N.LOC = n_distinct(LOCALITY.ID)) %>%
  ungroup() %>% 
  filter(N.LOC > 1) %>% 
  distinct(COMMON.NAME, LOCALITY.ID) %>%
  arrange(COMMON.NAME, LOCALITY.ID) 

prol_data <- prol_data %>% 
  # filtering for patch-species combos present in threshold lists at >1 locality
  right_join(filt, by = c("LOCALITY.ID", "COMMON.NAME")) %>% 
  arrange(OBSERVER.ID, LOCALITY.ID, SAMPLING.EVENT.IDENTIFIER, COMMON.NAME) %>% 
  # joining species info
  left_join(species_prolific) %>% 
  # filtering out locations where both species categories do not exist (after species filters)
  group_by(LOCALITY.ID) %>% 
  mutate(N.SPC = n_distinct(SP.CATEGORY)) %>% 
  ungroup() %>% 
  filter(N.SPC == 2) %>% 
  mutate(N.SPC = NULL)

# Lockdown and all months modelled separately, since the month variable in
# model otherwise pulled to a very different average 
prol_data <- bind_rows("LD" = prol_data %>% filter(MONTH %in% 4:5), 
                       "ALL" = prol_data, 
                       .id = "MONTHS.TYPE")


## saving data -------------------------------------------------------------

save(timeline, prol_data, 
     file = "00_data/bird_prolific.RData")

