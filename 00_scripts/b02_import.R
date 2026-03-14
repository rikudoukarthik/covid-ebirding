library(lubridate)
library(glue)
library(spdep)
library(sf)

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

load("00_data/data0_MY_b_slice.RData")



# timeline metadata
timeline <- data0_MY_b_slice_G %>% distinct(COVID, YEAR, MONTH, M.YEAR)


# highly birded personal locations -----------------------------------------------

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

rm(eBird_users, prolific_a, prolific_b)
