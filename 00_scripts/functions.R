
# get labelled months ---------------------------------------------------------------

# and COVID colour scheme

get_labelled_months <- function() {

  month_lab <- data.frame(M.YEAR = 2018:2021) %>% 
    group_by(M.YEAR) %>% 
    summarise(MONTH = 1:12 %>% month()) %>% 
    mutate(MONTH.LAB = month(MONTH, label = T)) %>% 
    mutate(MONTH = factor(MONTH, levels = month(c(6:12, 1:5))),
           M.YEAR = factor(M.YEAR),
           MONTH.LAB = factor(MONTH.LAB, levels = month(c(6:12, 1:5), label = T))) %>% 
    mutate(COL1 = ifelse(MONTH.LAB %in% c("Apr", "May"), "red", "black"), # for normal plots
           COL2 = ifelse((M.YEAR %in% 2019:2020) & 
                           (MONTH.LAB %in% c("Apr", "May")), "red", "black"), # for timeline 
           COL3 = ifelse(M.YEAR %in% 2019:2020, "red", "black")) %>% # for bird model plots
    arrange(M.YEAR, MONTH.LAB)
  
  return(month_lab)
  
}

# error functions -------------------------------------------------------------------

# getting CI limits from SE
get_CI_lims <- function(data) {
  
  data %>% 
    mutate(CI.L = PRED - 1.96*SE,
           CI.U = PRED + 1.96*SE)
  
}


# propagating SE when averaging
propagate_se_formean <- function(x, N = n()) {
  
  sqrt(sum(x^2)/N)
  
}


# averaging an estimate with its uncertainty
summarise_mean_and_se <- function(data_grouped, est, se) {
  
  data_grouped %>% 
    reframe({{ se }} := (sd({{ est }})/sqrt(n())) + propagate_se_formean({{ se }}),
            {{ est }} := mean({{ est }}))
  
}


# dividing an estimate with its uncertainty
summarise_div_and_se <- function(data_grouped, 
                                 est_num, se_num, est_denom, se_denom,
                                 est_out = {{ est_num }}, se_out = {{ se_num }}) {
  
  data_grouped %>% 
    reframe(TEMP = sqrt(
      ({{ se_num }} / {{ est_num }}) ^ 2 + ({{ se_denom }} / {{ est_denom }}) ^ 2
    ),
    {{ est_out }} := {{ est_num }} / {{ est_denom }},
    {{ se_out }} := {{ est_out }} * TEMP) %>% 
    dplyr::select(-TEMP)
  
}



# data metrics composite ------------------------------------------------------------

# min-max or 0-1 scaling (normalisation)
scale01 <- function(x, x_min, x_max, se = FALSE) {
  
  # for SE, we scale only by slope and not intercept
  
  if (se == TRUE) {
    (x)/(x_max - x_min)
  } else {
    (x - x_min)/(x_max - x_min)
  }
  
}


# composite of metrics summarise across months (COVID month yes/no) within the metric
# - transforms estimates into proportional changes from t0
join_metric_composite <- function(base, tojoin, metric_colname, 
                                  scale = "country", state = NULL) {
  
  # tojoin is a list (actually character vector of the list name) 
  # containing country and state results
  if (scale == "country") {
    pluck_which <- 1
  } else if (scale == "state") {
    pluck_which <- 2
    if (is.null(state)) return("Please choose a state to filter for!")
  }
  
  tojoin <- get(tojoin) %>% 
    pluck(pluck_which) %>% 
    dplyr::select(-c(PRED.LINK, SE.LINK, SE.L, YEAR, COVID, MONTH.LAB, COL1, COL2, COL3)) %>% 
    {if (scale == "state") {
      filter(., STATE == state)
    } else {
      .
    }} %>% 
    {if (metric_colname == "URBAN") {
      # we want to show declines due to COVID in the vis, so we need 1 - urb. prop.
      mutate(., PRED = 1 - PRED) %>% 
        get_CI_lims()
    } else {
      .
    }}
  
  reference <- tojoin %>% 
    filter(M.YEAR == 2018) %>% 
    distinct(MONTH, PRED, SE) %>% 
    magrittr::set_colnames(c("MONTH", "PRED.REF", "SE.REF"))
  
  tojoin <- tojoin %>% 
    left_join(reference) %>% 
    group_by(MONTH, M.YEAR, COVID.MONTH) %>% 
    summarise_div_and_se(PRED, SE, PRED.REF, SE.REF) %>% 
    mutate(METRIC = metric_colname)
  
  if (is.null(base)) {
    return(tojoin)
  } else {
    return(bind_rows(base, tojoin))
  }
  
}


### joinmapvars ########################################

# adapted from Ashwin's "addmapvars" function to allow use within data preparation steps
# i.e., inputs and outputs are data objects in the environment, no writing of files involved

# data: main data object to which map vars will be added (needs to be sliced already!)
# admin = T if DISTRICT and ST_NM from shapefile required, else F
# grids = T if grid cells (four resolutions) from shapefile required, else F


#### maps.RData must already be loaded
#### column names here are all uppercase, unlike in Ashwin's function


joinmapvars = function(data, admin = T, grids = T){
  
  require(tidyverse)
  require(sp)
  require(rgeos)
  
  # single object at group ID level (same group ID, same grid/district/state)
  temp0 <- data 
  
  
  ### add columns with DISTRICT and ST_NM to main data 
  
  if (admin == T) {
    
    temp = temp0 # separate object to prevent repeated slicing (intensive step)
    
    rownames(temp) = temp$GROUP.ID # only to setup adding the group.id column for the future left_join
    coordinates(temp) = ~LONGITUDE + LATITUDE # convert to SPDF
    proj4string(temp) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    temp = sp::over(temp, districtmap) %>% # returns only ATTRIBUTES of districtmap (DISTRICT and ST_NM)
      dplyr::select(1, 2) %>% 
      rename(DISTRICT = dtname,
             ST_NM = stname) %>% 
      rownames_to_column("GROUP.ID") 
    
    data = left_join(temp, data)
    
  }
  
  ### add grid cell info (at four resolutions) to main data
  
  if (grids == T) {
    
    temp = temp0
    
    rownames(temp) = temp$GROUP.ID
    coordinates(temp) = ~LONGITUDE + LATITUDE
    
    temp = sp::over(temp, gridmapg1) %>% 
      rownames_to_column("GROUP.ID") %>% 
      rename(GRIDG1 = id)
    
    data = left_join(temp, data)
    
    
    temp = temp0
    
    rownames(temp) = temp$GROUP.ID
    coordinates(temp) = ~LONGITUDE + LATITUDE
    
    temp = sp::over(temp, gridmapg2) %>% 
      rownames_to_column("GROUP.ID") %>% 
      rename(GRIDG2 = id)
    
    data = left_join(temp, data)
    
    
    temp = temp0
    
    rownames(temp) = temp$GROUP.ID
    coordinates(temp) = ~LONGITUDE + LATITUDE
    
    temp = sp::over(temp, gridmapg3) %>% 
      rownames_to_column("GROUP.ID") %>% 
      rename(GRIDG3 = id)
    
    data = left_join(temp, data)
    
    
    temp = temp0
    
    rownames(temp) = temp$GROUP.ID
    coordinates(temp) = ~LONGITUDE + LATITUDE
    
    temp = sp::over(temp, gridmapg4) %>% 
      rownames_to_column("GROUP.ID") %>% 
      rename(GRIDG4 = id)
    
    data = left_join(temp, data)
    
  }
  
  return(data)
  
}

### MODIS land-use map to get urban-non-urban information -----------------

getmodisdata <- function(){
  
# MCD12Q1 LULC data from December 2021, being aggregated to 2kmx2km ("UNU")
# this will be retained as a raster and also joined to the EBD data

require(raster) # masks dplyr functions, so have used package::function()
require(terra)
require(sf)
require(geodata)
# install.packages("luna", repos = "https://rspatial.r-universe.dev") # requires Rtools 4.0
require(luna)


### Obtaining MODIS LULC maps for urban/non-urban classification #######

## Obtaining MODIS data from within R for reproducibility, then saving the cropped and 
## masked raster of area of interest as .tif file.
## Walkthrough at https://rspatial.org/terra/modis/2-download.html

# The MODIS Land Cover Type Product (MCD12Q1) provides a suite of science data sets (SDSs) 
# that map global land cover at 500 meter spatial resolution at annual time step for six 
# different land cover legends. The maps were created from classifications of spectro-temporal 
# features derived of data from the Moderate Resolution Imaging Spectroradiometer (MODIS).


# ## Below is to see list of MODIS products. We already know which one we want, so ignore.
# MODISprod <- getProducts("^MOD|^MYD|^MCD")
# MODISprod


prod <- "MCD12Q1" # short name of product of interest
start <- "2021-12-31" # time period of mapping
end <- "2021-12-31" # time period of mapping


# creating folder to store MODIS data that will be downloaded
MODISpath <- "00_data/in_LULC_MODIS/"
dir.create(MODISpath, showWarnings=FALSE)


# ## create RDS file containing login credentials for EOSDIS account
# cred <- data.frame(UN = "xxx",
#                    PW = "yyy")
# saveRDS(cred, "userpass.rds")
userpass <- readRDS("userpass.rds")


## downloading MODIS tiled data (16 tiles for our defined area of interest, India)
# resulting object is a list of file names
modis <- luna::getModis(prod, start, end,
                        aoi = india, download = T, path = MODISpath,
                        username = userpass$UN, password = userpass$PW)

print("Downloaded MODIS data")

# choose the band (SDS) that you want using sds[] option and write GTiff files.
# University of Maryland SDS (LC_Type2) is SDS1
for (i in (modis)) {
  sds <- terra::sds(i)
  modis2 <- terra::writeRaster(sds[1], filename = paste0(i, ".tif"))
}


## listing all the different tiles present as .tif files, to later merge into one
rast_list <- list.files(path = MODISpath,
                        pattern = "MCD12Q1*.*tif",
                        all.files = TRUE, full.names = TRUE)

# merge into one raster
rast_full <- rast_list %>%
  lapply(raster::raster) %>% # rasterise the different files
  do.call(what = raster::merge) %>% # merge the files
  terra::rast()

print("Merged multiple MODIS rasters into one")


## projecting from MODIS sinusoidal to WGS84 (flat 2D)
projto <- raster::crs(india)
# next step takes a few minutes
# "near" because LC values are categories despite being numerical
rast_proj <- terra::project(rast_full, projto, method = "near")


## cropping and masking the raster to our area of interest
rast_aoi <- terra::crop(rast_proj, india)
rast_aoi <- terra::mask(rast_aoi, india)

print("Projected to WGS84, and cropped and masked to AOI")


raster::writeRaster(rast_aoi, paste0(MODISpath, "/in_LULC_MODIS.tif"), overwrite = TRUE)

rm(prod, start, end, MODISpath, userpass, modis, sds, modis2, rast_list, rast_full, projto, rast_proj, rast_aoi)

print("Created merged TIFF file")

### Modifying the MODIS data #######

## read in cropped and masked .tif
rast <- raster::raster("00_data/in_LULC_MODIS/in_LULC_MODIS.tif")


## reclassifying the raster as urban vs. non-urban 

# (help from: http://www.wvview.org/spatial_analytics/Raster_Analysis/_site/index.html ; http://zevross.com/blog/2015/03/30/map-and-analyze-raster-data-in-r/ )
# category 13 corresponds to "urban and built-up lands"
categs <- unique(values(rast))
reclass <- data.frame(old = categs) %>% 
  mutate(new = case_when(old == 13 ~ 1, TRUE ~ 0)) %>% 
  as.matrix()
# ocean also becomes "non-urban" but that's fine since we won't be looking at those locations

rast_UNU <- raster::reclassify(rast, rcl = reclass)


## changing resolution of data (MODIS data has 500m*500m resolution)
# https://stackoverflow.com/questions/37956422/resample-raster

# aggregating to 2km*2km
rast_UNU <- rast_UNU %>% 
  raster::aggregate(fact = 2, fun = rast_agg_fn) %>% 
  raster::aggregate(fact = 2, fun = rast_agg_fn)

save(rast_UNU, file = "00_data/rast_UNU.RData")

print("Reclassified and aggregated MODIS data")

# 

}

### data quality filters and preparation -----------------

# data file path (first time .txt which then creates .RData for future uses; functionised)
# latest Ashwin's maps.RData file (having areas)
# path to list of group accounts to be filtered out
# path to classification of year-month as COVID categories
# path to raster data at 2x2 scale


data_qualfilt_prep <- function(rawdatapath, senspath,
                               datapath, groupaccspath, covidclasspath,
                               rast_UNU_path,
                               maxvel = 20, minsut = 2){
  
  require(tidyverse)
  require(lubridate)
  require(tictoc)

  # if first time, importing data from EBD .txt and saving as .RData for future
  # else directly loading .RData
  
  if (!file.exists(datapath) & file.exists(rawdatapath)) {
    
    # variables required in data object
    preimp <- c("CATEGORY","EXOTIC.CODE","COMMON.NAME","OBSERVATION.COUNT",
                "LOCALITY.ID","LOCALITY.TYPE","REVIEWED","APPROVED","STATE","COUNTY","LAST.EDITED.DATE",
                "LATITUDE","LONGITUDE","OBSERVATION.DATE","TIME.OBSERVATIONS.STARTED","OBSERVER.ID",
                "PROTOCOL.TYPE","DURATION.MINUTES","EFFORT.DISTANCE.KM","LOCALITY","BREEDING.CODE",
                "NUMBER.OBSERVERS","ALL.SPECIES.REPORTED","GROUP.IDENTIFIER","SAMPLING.EVENT.IDENTIFIER",
                "TRIP.COMMENTS","SPECIES.COMMENTS", "HAS.MEDIA")
    

    ### main EBD ###
    
    # this method using base R import takes only 373 sec with May 2022 release
    nms <- names(read.delim(rawdatapath, nrows = 1, sep = "\t", header = T, quote = "", 
                            stringsAsFactors = F, na.strings = c(""," ", NA)))
    nms[!(nms %in% preimp)] <- "NULL"
    nms[nms %in% preimp] <- NA
    data <- read.delim(rawdatapath, colClasses = nms, sep = "\t", header = T, quote = "",
                       stringsAsFactors = F, na.strings = c(""," ",NA)) 
    
    # # tidy import takes way longer, a total of 877 sec, but could be useful for smaller data
    # data <- read_delim(rawpath, col_select = preimp,
    #                    name_repair = make.names, # base R nomencl. with periods for spaces
    #                    quote = "", na = c(""," ", NA), show_col_types = F)
    
    
    ### sensitive species ###
    nms1 <- names(read.delim(senspath, nrows = 1, sep = "\t", header = T, quote = "", 
                             stringsAsFactors = F, na.strings = c(""," ", NA)))
    nms1[!(nms1 %in% preimp)] <- "NULL"
    nms1[nms1 %in% preimp] <- NA
    senssp <- read.delim(senspath, colClasses = nms1, sep = "\t", header = T, quote = "",
                         stringsAsFactors = F, na.strings = c(""," ",NA))
    
    ### combing the two ###
    data <- bind_rows(data, senssp) %>% 
      # filtering unvetted as well as exotic species (provisional and escapee)
      filter(APPROVED == 1 & !(EXOTIC.CODE %in% c("P", "X")))
    
    
    data <- data %>% 
      # trimming whitespace in breeding codes
      mutate(BREEDING.CODE = str_trim(BREEDING.CODE)) %>% 
      # group ID and dates
      mutate(GROUP.ID = ifelse(is.na(GROUP.IDENTIFIER), SAMPLING.EVENT.IDENTIFIER, GROUP.IDENTIFIER), 
             OBSERVATION.DATE = as.Date(OBSERVATION.DATE), 
             YEAR = year(OBSERVATION.DATE), 
             MONTH = month(OBSERVATION.DATE),
             DAY.M = day(OBSERVATION.DATE)) %>% 
      # migratory year and month information
      mutate(M.YEAR = if_else(MONTH > 5, YEAR, YEAR-1), # from June to May
             M.MONTH = if_else(MONTH > 5, MONTH-5, 12-(5-MONTH))) 
    
    
    tic("Completed data import, quality filtering and preparation in")
    
  } else if (!file.exists(datapath) & !file.exists(rawdatappath)) {
    
    print("Raw data file of EBD does not exist!")
    
  } else if (file.exists(datapath)) {
    
    tic("Completed data quality filtering and preparation in")
    
    ### importing from usual modified ebd RData
    load(datapath) 
    
  }
  

  ### list of group accounts to be filtered
  groupaccs <- read_csv(groupaccspath) %>% 
    mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                                GA.2 == 1 ~ "GA.2", 
                                TRUE ~ "NG"))
  # for bird analyses
  filtGA_b <- groupaccs %>% 
    filter(CATEGORY == "GA.1") %>% 
    dplyr::select(OBSERVER.ID)
  # for birder/data analyses
  filtGA_d <- groupaccs %>% 
    filter(CATEGORY == "GA.1" | CATEGORY == "GA.2") %>% 
    dplyr::select(OBSERVER.ID)
  
  
  ### COVID classification
  covidclass <- read_csv(covidclasspath)
  
  
  require(tidyverse)
  require(lubridate)
  require(terra)
  require(sp)
  
  
  ### adding UNU information ######
  
  load(rast_UNU_path)

  # adding extracted land cover values at eBird data points with grid cell information at both scales
  # Pakistan-administered Kashmir lists produce NA for URBAN
  lists_UNU <- data %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    mutate(URBAN = raster::extract(rast_UNU, cbind(LONGITUDE, LATITUDE)),
           NONURBAN = if_else(URBAN == 1, 0, 1),
           # 2km*2km
           SUBCELL.ID = raster::cellFromXY(rast_UNU, cbind(LONGITUDE, LATITUDE))) %>% 
    filter(!is.na(URBAN)) %>% 
    dplyr::select(-LONGITUDE, -LATITUDE)
  
  save(lists_UNU, file = "00_data/lists_UNU.RData")
  
  print("Added UNU information")
  
  ### new observer data (to calculate no. of new observers metric) #######

  new_obsr_data <- data %>%
    dplyr::select(c("YEAR", "MONTH", "STATE", "SAMPLING.EVENT.IDENTIFIER",
             "LAST.EDITED.DATE", "OBSERVATION.DATE", "OBSERVER.ID")) %>%
    mutate(LAST.EDITED.DATE = ymd_hms(LAST.EDITED.DATE)) %>%
    group_by(OBSERVER.ID) %>%
    arrange(LAST.EDITED.DATE) %>%
    ungroup() %>%
    distinct(OBSERVER.ID, .keep_all = TRUE) %>%
    mutate(YEAR = year(LAST.EDITED.DATE),
           MONTH = month(LAST.EDITED.DATE),
           M.YEAR = if_else(MONTH > 5, YEAR, YEAR-1), # from June to May
           M.MONTH = if_else(MONTH > 5, MONTH-5, 12-(5-MONTH))) %>%
    filter(M.YEAR >= 2018) %>%
    left_join(covidclass) %>%
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>%
    rename(LE.YEAR = M.YEAR,
           LE.MONTH = M.MONTH) %>%
    mutate(LE.YEAR = factor(LE.YEAR, levels = seq(2018, 2022, by = 1)),
           LE.MONTH = factor(LE.MONTH, levels = seq(1, 12, by = 1)),
           YEAR = year(OBSERVATION.DATE),
           MONTH = month(OBSERVATION.DATE),
           M.YEAR = if_else(MONTH > 5, YEAR, YEAR-1), # from June to May
           M.MONTH = if_else(MONTH > 5, MONTH-5, 12-(5-MONTH))) %>% 
    mutate(YEAR = factor(YEAR, levels = seq(2018, 2022, by = 1)),
           MONTH = factor(MONTH, levels = seq(1, 12, by = 1)),
           M.YEAR = factor(M.YEAR, levels = seq(2018, 2021, by = 1)),
           M.MONTH = factor(M.MONTH, levels = seq(1, 12, by = 1)))

  # filtering
  new_obsr_data <- new_obsr_data %>% anti_join(filtGA_d)
  save(new_obsr_data, file = "00_data/new_obsr_data.RData")
  
  print("Obtained new observer data")
  
  ### main data filtering ######
  
  # to calculate list length aka number of species
  temp0 <- data %>% 
    mutate(CATEGORY = case_when(CATEGORY == "domestic" & 
                                  COMMON.NAME == "Rock Pigeon" ~ "species",
                                TRUE ~ CATEGORY)) %>% 
    filter(CATEGORY %in% c("issf","species")) %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(NO.SP = n_distinct(COMMON.NAME))
  
  ## for bird analyses (group acc filter different) ####
  
  data0_MY_b <- data %>% 
    # this data (including latter months of 2018 only needed for bird behaviour section)
    # so will later save separate data filtering out 2018 months
    filter(M.YEAR >= 2018) %>%
    anti_join(filtGA_b) %>% # removing data from group accounts
    # creating COVID factor
    left_join(covidclass) %>% 
    # adding UNU information for every GROUP.ID
    left_join(lists_UNU) %>% 
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>% 
    # NO.SP column
    left_join(temp0) %>% 
    mutate(NO.SP = replace_na(NO.SP, 0),
           DATETIME = as_datetime(paste(OBSERVATION.DATE,
                                        TIME.OBSERVATIONS.STARTED)),
           HOUR = hour(DATETIME),
           MIN = minute(DATETIME),
           SPEED = EFFORT.DISTANCE.KM*60/DURATION.MINUTES, # kmph
           SUT = NO.SP*60/DURATION.MINUTES, # species per hour
           # calculate hour checklist ended
           HOUR.END = floor((HOUR*60 + MIN + DURATION.MINUTES)/60))

  
  ### exclude records based on various criteria 
  
  # false complete lists (without duration info & 3 or fewer species, <3min, low SUT)
  temp1 <- data0_MY_b %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE != "Incidental") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter((NO.SP <= 3 & is.na(DURATION.MINUTES)) | 
             (DURATION.MINUTES < 3) | 
             (SUT < minsut & NO.SP <= 3)) %>%
    distinct(GROUP.ID)
  
  # getting list of GROUP.IDs inside IN boundary for pelagic filter
  temp2 <- data0_MY_b %>% 
    ungroup() %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    terra::vect(geom = c("LONGITUDE","LATITUDE"), crs = crs(india)) %>% 
    terra::intersect(india) %>% 
    terra::as.data.frame() %>% 
    dplyr::select(GROUP.ID) 

  # speed and distance filter for travelling lists
  temp3 <- data0_MY_b %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE == "Traveling") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter((SPEED > maxvel) | (EFFORT.DISTANCE.KM > 50)) %>%
    distinct(GROUP.ID)
  
  # true completeness + other filters
  data0_MY_b <- data0_MY_b %>%
    # nocturnal filter
    mutate(NOCT.FILTER = case_when((!is.na(HOUR) & ((HOUR <= 4 & HOUR.END <= 4) | 
                                                      (HOUR >= 20 & HOUR.END <= 28))) ~ 0, 
                                   TRUE ~ 1)) %>% 
    filter((ALL.SPECIES.REPORTED == 1) & 
             (NOCT.FILTER == 1) &
             # true completeness
             !(GROUP.ID %in% temp1$GROUP.ID |
                GROUP.ID %in% temp3$GROUP.ID |
                (ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE == "Incidental")) &
             # pelagic filter
             (GROUP.ID %in% temp2$GROUP.ID)) %>% 
    dplyr::select(-BREEDING.CODE, -SPEED, -SUT, -MIN, -DATETIME, -HOUR.END, -NOCT.FILTER)
  
  
  # making month and year ordered factors
  data0_MY_b <- data0_MY_b %>% 
    mutate(M.YEAR = factor(M.YEAR, levels = seq(2018, 2021, by = 1)),
           MONTH = factor(MONTH, levels = seq(1, 12, by = 1)),
           M.MONTH = factor(M.MONTH, levels = seq(1, 12, by = 1)))
  
  
  # sliced data which is what is required for analyses
  data0_MY_b_slice_S <- data0_MY_b %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
    slice(1) %>% 
    ungroup()
  data0_MY_b_slice_G <- data0_MY_b %>% 
    group_by(GROUP.ID) %>% 
    slice(1) %>%
    ungroup()
  
  # adding map variables (CELL.ID) to main data
  load("00_data/maps.RData", envir = .GlobalEnv) # Ashwin's maps data
  
  lists_grids <- data0_MY_b_slice_G %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    joinmapvars() %>% 
    # 25x25 grid cells
    rename(CELL.ID = GRIDG1)

  data0_MY_b <- data0_MY_b %>% left_join(lists_grids)
  data0_MY_b_slice_S <- data0_MY_b_slice_S %>% left_join(lists_grids)
  data0_MY_b_slice_G <- data0_MY_b_slice_G %>% left_join(lists_grids)
  
  save(data0_MY_b, file = "00_data/data0_MY_b.RData")
  save(data0_MY_b_slice_S, data0_MY_b_slice_G, file = "00_data/data0_MY_b_slice.RData")
  
  print("Complete main data filtering for bird analyses!")
  
  
  ## for birder analyses (group acc filter different) ####
  
  data0_MY_d <- data %>% 
    # this data (including latter months of 2018 only needed for bird behaviour section)
    # so will later save separate data filtering out 2018 months
    filter(M.YEAR >= 2018) %>%
    anti_join(filtGA_d) %>% # removing data from group accounts
    # creating COVID factor
    left_join(covidclass) %>% 
    # adding UNU information for every GROUP.ID
    left_join(lists_UNU) %>% 
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>% 
    # NO.SP column
    left_join(temp0) %>% 
    mutate(NO.SP = replace_na(NO.SP, 0),
           DATETIME = as_datetime(paste(OBSERVATION.DATE,
                                        TIME.OBSERVATIONS.STARTED)),
           HOUR = hour(DATETIME),
           MIN = minute(DATETIME),
           SPEED = EFFORT.DISTANCE.KM*60/DURATION.MINUTES, # kmph
           SUT = NO.SP*60/DURATION.MINUTES, # species per hour
           # calculate hour checklist ended
           HOUR.END = floor((HOUR*60 + MIN + DURATION.MINUTES)/60))
  
  
  ### exclude records based on various criteria 
  
  # false complete lists (without duration info & 3 or fewer species, <3min, low SUT)
  temp1 <- data0_MY_d %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE != "Incidental") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter((NO.SP <= 3 & is.na(DURATION.MINUTES)) | 
             (DURATION.MINUTES < 3) | 
             (SUT < minsut & NO.SP <= 3)) %>%
    distinct(GROUP.ID)
  
  # getting list of GROUP.IDs inside IN boundary for pelagic filter
  temp2 <- data0_MY_d %>% 
    ungroup() %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    terra::vect(geom = c("LONGITUDE","LATITUDE"), crs = crs(india)) %>% 
    terra::intersect(india) %>% 
    terra::as.data.frame() %>% 
    dplyr::select(GROUP.ID) 
  
  # speed and distance filter for travelling lists
  temp3 <- data0_MY_d %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE == "Traveling") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter((SPEED > maxvel) | (EFFORT.DISTANCE.KM > 50)) %>%
    distinct(GROUP.ID)
  
  # true completeness + other filters
  data0_MY_d <- data0_MY_d %>%
    # nocturnal filter
    mutate(NOCT.FILTER = case_when((!is.na(HOUR) & ((HOUR <= 4 & HOUR.END <= 4) | 
                                                      (HOUR >= 20 & HOUR.END <= 28))) ~ 0, 
                                   TRUE ~ 1)) %>% 
    filter((ALL.SPECIES.REPORTED == 1) & 
             (NOCT.FILTER == 1) &
             # true completeness
             !(GROUP.ID %in% temp1$GROUP.ID |
                 GROUP.ID %in% temp3$GROUP.ID |
                 (ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE == "Incidental")) &
             # pelagic filter
             (GROUP.ID %in% temp2$GROUP.ID)) %>% 
    dplyr::select(-BREEDING.CODE, -SPEED, -SUT, -MIN, -DATETIME, -HOUR.END, -NOCT.FILTER)
  
  
  # making month and year ordered factors
  data0_MY_d <- data0_MY_d %>% 
    mutate(M.YEAR = factor(M.YEAR, levels = seq(2018, 2021, by = 1)),
           MONTH = factor(MONTH, levels = seq(1, 12, by = 1)),
           M.MONTH = factor(M.MONTH, levels = seq(1, 12, by = 1)))
  
  
  # sliced data which is what is required for analyses
  data0_MY_d_slice_S <- data0_MY_d %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
    slice(1) %>% 
    ungroup()
  data0_MY_d_slice_G <- data0_MY_d %>% 
    group_by(GROUP.ID) %>% 
    slice(1) %>%
    ungroup()
  
  # adding map variables (CELL.ID) to main data
  load("00_data/maps.RData", envir = .GlobalEnv) # Ashwin's maps data
  
  lists_grids <- data0_MY_d_slice_G %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    joinmapvars() %>% 
    # 25x25 grid cells
    rename(CELL.ID = GRIDG1)
  
  data0_MY_d <- data0_MY_d %>% left_join(lists_grids)
  data0_MY_d_slice_S <- data0_MY_d_slice_S %>% left_join(lists_grids)
  data0_MY_d_slice_G <- data0_MY_d_slice_G %>% left_join(lists_grids)
  
  save(data0_MY_d, file = "00_data/data0_MY_d.RData")
  save(data0_MY_d_slice_S, data0_MY_d_slice_G, file = "00_data/data0_MY_d_slice.RData")
  
  print("Completed main data filtering for birder analyses!")
  
  
  
  toc()
  
}

### bootstrapping confidence -----------------

# Useful resources: 
# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/variability-and-uncertainty-standard-deviations-standard-errors-confidence-intervals.html#bootstrap 
# https://websites.pmc.ucsc.edu/~mclapham/Rtips/resampling

# For some metrics like birding distance, it is possible to calculate SE from the data itself 
# ("expected SE" = SD of sample / root sample size). But when there is a possibility to 
# calculate the empirical SE instead (via bootstrapping) which is more accurate, there is 
# no point in going for the former. 

# here we are calculating SE but what we want is CIs, and we use +-1.96*SE to calculate them.
# this assumes Gaussian distribution about the mean/median (symmetric CIs). we can make this 
# assumption only because sample size is sufficiently large. 
# when we need to propagate from state to nation level, this assumption is necessary.
# thus, for all the initial metrics, this is okay.

###
# when there is no issue of propagating SEs, always better to bootstrap CIs directly,
# without making the Gaussian assumption. this allows asymmetric CIs which are more likely
# with counts, proportions, etc.

boot_conf = function(x, fn = mean, B = 1000) {
  
  1:B %>%
    # For each iteration, generate a sample of x with replacement
    map(~ x[sample(1:length(x), replace = TRUE)]) %>%
    # Obtain the fn estimate for each bootstrap sample
    map_dbl(fn)
  
}


# metrics that use this function
# fidelity, time, dist, duration, length, spatial (net change), new




### bootstrapping confidence from GLMMs (bootMer) -------------------------------------

# adapted from Ashwin's function for SoIB

boot_conf_GLMM = function(model, 
                          new_data, # separately specify dataframe with vars for model
                          new_data_string, # string for clusterExport()
                          re_form = NA,
                          nsim = 1000,
                          pred_type = "link")
{

  require(tidyverse)
  require(lme4)
  require(VGAM)
  require(parallel) # to parallelise bootstrap step

  pred_fun <- function(model) {
    predict(model, newdata = new_data, type = pred_type, re.form = re_form, 
            allow.new.levels = TRUE)
    # not specifying type = "response" because will later transform prediction along with SE
  }
  
  par_cores <- max(1, floor(detectCores()/2))
  par_cluster <- makeCluster(rep("localhost", par_cores), outfile = "log.txt")
  clusterEvalQ(par_cluster, library(lme4))
  clusterExport(par_cluster, varlist = new_data_string)

  print(glue("Using {par_cores} cores."))
  
  pred_bootMer <- bootMer(model, nsim = nsim, FUN = pred_fun,
                          parallel = "snow", 
                          use.u = FALSE, type = "parametric", 
                          ncpus = par_cores, cl = par_cluster)
  
  stopCluster(par_cluster)

  return(pred_bootMer$t)
  
}

### splitting simulations for bootstrapping (bootMer) -------------------------------------

# because doing all 1000 sims at once results in error

split_par_boot <- function(model, 
                           new_data, # separately specify dataframe with vars for model
                           new_data_string, # string for clusterExport()
                           re_form = NA, 
                           mode = "normal",
                           pred_type = "link") {
  
  if (mode == "extra") {
    
    count <- 0
    prediction1 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction2 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction3 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction4 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction5 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction6 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction7 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction8 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction9 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 50,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction10 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction11 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction12 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction13 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction14 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction15 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction16 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction17 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction18 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction19 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    prediction20 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 50,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/20 sets of 50 simulations completed"))
    
    
    prediction <- rbind(prediction1, prediction2, prediction3, prediction4, prediction5,
                        prediction6, prediction7, prediction8, prediction9, prediction10,
                        prediction11, prediction12, prediction13, prediction14, prediction15,
                        prediction16, prediction17, prediction18, prediction19, prediction20)
    
    return(prediction)
    
  }
  
  if (mode == "normal") {
    
    count <- 0
    prediction1 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction2 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction3 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction4 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction5 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction6 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction7 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction8 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction9 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    prediction10 <- boot_conf_GLMM(model, 
                                   new_data, 
                                   new_data_string, 
                                   nsim = 100,
                                   re_form = re_form,
                                   pred_type = pred_type)
    count <- count + 1
    print(glue("{count}/10 sets of 100 simulations completed"))
    
    
    prediction <- rbind(prediction1, prediction2, prediction3, prediction4, prediction5,
                        prediction6, prediction7, prediction8, prediction9, prediction10)
    
    return(prediction)
    
  }
  
  if (mode == "normal_low") {
    
    # only 500 sims in total
    
    count <- 0
    
    count <- count + 1
    tictoc::tic((glue("{count}/5 sets of 100 simulations completed")))
    prediction1 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    tictoc::toc()
    
    count <- count + 1
    tictoc::tic((glue("{count}/5 sets of 100 simulations completed")))
    prediction2 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    tictoc::toc()
    
    
    count <- count + 1
    tictoc::tic((glue("{count}/5 sets of 100 simulations completed")))
    prediction3 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    tictoc::toc()
    
    
    count <- count + 1
    tictoc::tic((glue("{count}/5 sets of 100 simulations completed")))
    prediction4 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    tictoc::toc()
    
    
    count <- count + 1
    tictoc::tic((glue("{count}/5 sets of 100 simulations completed")))
    prediction5 <- boot_conf_GLMM(model, 
                                  new_data, 
                                  new_data_string, 
                                  nsim = 100,
                                  re_form = re_form,
                                  pred_type = pred_type)
    tictoc::toc()

    
    prediction <- rbind(prediction1, prediction2, prediction3, prediction4, prediction5)
    
    return(prediction)
    
  }
  
}

### raster aggregation -----------------

# Function to supply in raster::aggregate(), with 25% threshold for classification
# (works for any binary classification to be aggregated; here it is urban/non-)
# https://stackoverflow.com/questions/23369579/aggregate-raster-in-r-with-na-values 

rast_agg_fn <- function(x, ...){
  if_else(sum(x %in% 1) >= 0.25*length(x), 1, 0)
}

# x = c(rep(1, 2), rep(0, 8), rep(1, 10), rep(0, 30)) # for testing


### raster calc: log proportional change in number of lists -----------------

rast_logpropchange <- function(x, y, k = 1, emptycheck = F)  {
  
  # To return NA if both values being compared are zero 
  # (mainly useful for maps of proportional change to visualise empty grid cells)
  # x is first, y is second, so change from x to y
  
  # # when no. of lists has been log transformed (this not vectorised)
  # # log(n2) - log(n1) = log(n2/n1)
  # return(y-x)
  
  # when value has not been transformed to prevent NA

  case_when(emptycheck == T & (x == 0 & y == 0) ~ NA_real_, 
            TRUE ~ log10((y + k)/(x + k)))
  
  # proportion of urban birding in the cell changed by y/x [0,+n] times.
  
}


### raster calc: proportional change in number of lists -----------------

rast_propchange <- function(x, y, k = 1, emptycheck = F)  {
  
  case_when(emptycheck == T & (x == 0 & y == 0) ~ NA_real_, 
            TRUE ~ (y + k)/(x + k)) # percentage
  
  # proportion of urban birding in the cell changed by y/x [0,+n] times.

}

### change in spatial spread/clustering metric -----------------

s_spread_propchange <- function(metric, x, y, k = 1)  {

  # when both values being compared are zero, converting to NA which can later be removed
  # when LHS is zero, adding k = 1 to both values (to avoid infinity due to 0 in denominator)
  # when RHS is zero, adding nothing

  # for the proportional change of number of lists, setting a threshold of 10 lists per district
  # to avoid large values of proportional change caused by them
  
  case_when(x == 0 & y == 0 ~ NA_real_,
            x < 10 & metric == "NO.LISTS" ~ NA_real_, 
            x == 0 & metric != "NO.LISTS" & y != 0 ~ ((y + k) - (x + k))/(x + k), 
            TRUE ~ (y - x)/x)

}

s_spread_change <- function(x, y, k = 1)  {

  # when both values being compared are zero, converting to NA which can later be removed
  case_when(x == 0 & y == 0 ~ NA_real_, 
            TRUE ~ (y - x))

}


### create nb object but omit grid cells having no neighbours -----------------

# these need to be removed when calculating weights for the Moran analyses
# also creates data object "nb_zero" in environment: list of cells with no neighbours
# BUT needs there to be CELL.ID column in input data object

# ref: https://stackoverflow.com/a/57378930/13000254


# only realised after creating this function that there was no need to do this as
# localmoran() has zero.policy argument as well like poly2nb()


require(spdep)

poly2omit0nb <- function(data) {
  
  nbobj <- data %>% poly2nb()
  
  # getting number of neighbours for each cell
  nbcount <- card(nbobj)
  
  # list of cells with no neighbours
  nb_zero <- which(nbcount == 0)
  # list of cells with neighbours
  nb_true <- which(nbcount != 0)
  
  
  # actual CELL.IDs of cells with zero neighbours
  cell_zero <- data.frame(CELL.ID = data$CELL.ID[nb_zero])
  assign("cell_zero", cell_zero, envir = .GlobalEnv)


  # only return the cell if it has neighbours
  return(subset.nb(nbobj, (1:length(nbobj) %in% nb_true)))
  
}


### ggplot for non-model bird reporting patterns -----------------

gg_b_nonmodel <- function(data, region, time) {
  
  # cur_city_list should already be in environment
  
  require(tidyverse)
  require(patchwork)
  
  if (region == "city"){
    plot_title <- glue("{unique(cur_city_list$CITY)} city")
  } else if (region == "state"){
    plot_title <- glue("{unique(cur_city_list$STATE)} state")
  }
  
  if (time == "monthly") {
    
    ((ggplot(filter(data, SP.CATEGORY == "U"), 
             aes(MONTH, REP.FREQ, colour = M.YEAR)) +
        geom_point(size = 2, position = position_dodge(0.8)) +
        geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                      size = 1.25, width = 0.6, position = position_dodge(0.8)) +
        scale_colour_manual(values = covid_palette, name = "Migratory\nyear") +
        facet_wrap(~ COMMON.NAME, dir = "h", ncol = 3, 
                   strip.position = "left", scales = "free_y") +
        labs(title = "Urban species",
             x = "Month", y = "Reporting frequency")) /
       (ggplot(filter(data, SP.CATEGORY == "R"), 
               aes(MONTH, REP.FREQ, colour = M.YEAR)) +
          geom_point(size = 2, position = position_dodge(0.8)) +
          geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                        size = 1.25, width = 0.6, position = position_dodge(0.8)) +
          scale_colour_manual(values = covid_palette, name = "Migratory\nyear") +
          facet_wrap(~ COMMON.NAME, dir = "h", ncol = 3, 
                     strip.position = "left", scales = "free_y") +
          labs(title = "Non-urban species",
               x = "Month", y = "Reporting frequency"))) +
      plot_layout(guides = "collect", heights = c(4, 4)) +
      plot_annotation(title = plot_title) &
      theme(strip.text = element_text(size = 7))
    
  } else if (time == "yearly") {
    
    ((ggplot(filter(data, SP.CATEGORY == "U"), aes(M.YEAR, REP.FREQ)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                      size = 1.25, width = 0.3) +
        facet_wrap(~ COMMON.NAME, dir = "h", ncol = 3, 
                   strip.position = "left", scales = "free_y") +
        labs(title = "Urban species",
             x = "Migratory year", y = "Reporting frequency")) /
       (ggplot(filter(data, SP.CATEGORY == "R"), aes(M.YEAR, REP.FREQ)) +
          geom_point(size = 3) +
          geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                        size = 1.25, width = 0.3) +
          facet_wrap(~ COMMON.NAME, dir = "h", ncol = 3, 
                     strip.position = "left", scales = "free_y") +
          labs(title = "Non-urban species",
               x = "Migratory year", y = "Reporting frequency"))) +
      plot_layout(guides = "collect", heights = c(4, 4)) +
      plot_annotation(title = plot_title) &
      theme(strip.text = element_text(size = 7))
    
  }
  
}

### ggplot for model bird reporting patterns -----------------

gg_b_model <- function(data, type, data_points) {
  

  require(tidyverse)
  require(patchwork)
  
  if (type == "overall") {
    
    plot_title <- glue("{state_name} state")
    plot_subtitle <- paste0(
      "Predicted reporting frequencies of ",
      n_distinct(data_occ$COMMON.NAME), " species (",
      n_distinct(filter(data_occ, SP.CATEGORY == "U")$COMMON.NAME), " urban, ",
      n_distinct(filter(data_occ, SP.CATEGORY == "R")$COMMON.NAME), " rural)",
      " in three separate monthwise models")
    
    # different y lims for different states
    if (state_name == "Karnataka") {
      plot_ylims <- c(0.13, 0.35)
    } else if (state_name == "Kerala") {
      plot_ylims <- c(0.05, 0.22)
    } else if (state_name == "Maharashtra") {
      plot_ylims <- c(0.175, 0.35)
    } else if (state_name == "Assam") {
      plot_ylims <- c(0.1, 0.5)
    }
    
    
    model_plot <- (ggplot(filter(data, MONTHS.TYPE == "LD"), 
                          aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = plot_ylims) +
                     labs(title = "For the months of April and May",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     geom_point(size = 1.75, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1, width = 0.15, position = position_dodge(0.5)) |
                     ggplot(filter(data, MONTHS.TYPE == "NL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = plot_ylims) +
                     labs(title = "For other ten months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     geom_point(size = 1.75, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1, width = 0.15, position = position_dodge(0.5)) |
                     ggplot(filter(data, MONTHS.TYPE == "ALL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = plot_ylims) +
                     labs(title = "For all twelve months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     geom_point(size = 1.75, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1, width = 0.15, position = position_dodge(0.5))) +
      plot_layout(guides = "collect") +
      plot_annotation(title = plot_title,
                      subtitle = plot_subtitle) 
    
    return(model_plot)
    
  } else if (type == "prolific") {
    
    plot_title <- "Change in bird species reporting from locations with consistent effort"
    plot_subtitle <- paste0(
      "Predicted reporting frequencies of ",
      n_distinct(data_prolific$COMMON.NAME), " species (",
      n_distinct(filter(data_prolific, SP.CATEGORY == "U")$COMMON.NAME), " urban, ",
      n_distinct(filter(data_prolific, SP.CATEGORY == "R")$COMMON.NAME), " rural)",
      " in three separate monthwise models")
    
    
    model_plot <- (ggplot(filter(data, MONTHS.TYPE == "LD"), 
                          aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = c(0.1, 1.0)) +
                     labs(title = "For the months of April and May",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     # data points
                     geom_point(data = filter(data_points, MONTHS.TYPE == "LD"), 
                                aes(group = PATH.GROUP),
                                size = 3, alpha = 0.25,
                                position = position_dodge(0.25)) + 
                     geom_path(data = filter(data_points, MONTHS.TYPE == "LD"), 
                               aes(x = as.numeric(M.YEAR), y = PRED, group = PATH.GROUP),
                               size = 1, alpha = 0.15, position = position_dodge(0.25)) + 
                     # main points
                     geom_point(size = 3, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1.5, width = 0.2, position = position_dodge(0.5)) |
                    ggplot(filter(data, MONTHS.TYPE == "NL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = c(0.1, 1.0)) +
                     labs(title = "For other ten months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     # data points
                     geom_point(data = filter(data_points, MONTHS.TYPE == "NL"), 
                                aes(group = PATH.GROUP),
                                size = 3, alpha = 0.25,
                                position = position_dodge(0.25)) + 
                     geom_path(data = filter(data_points, MONTHS.TYPE == "NL"), 
                               aes(x = as.numeric(M.YEAR), y = PRED, group = PATH.GROUP),
                               size = 1, alpha = 0.15, position = position_dodge(0.25)) + 
                     # main points
                     geom_point(size = 3, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1.5, width = 0.2, position = position_dodge(0.5)) |
                    ggplot(filter(data, MONTHS.TYPE == "ALL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = c(0.1, 1.0)) +
                     labs(title = "For all twelve months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     # data points
                     geom_point(data = filter(data_points, MONTHS.TYPE == "ALL"), 
                                aes(group = PATH.GROUP),
                                size = 3, alpha = 0.25,
                                position = position_dodge(0.25)) + 
                     geom_path(data = filter(data_points, MONTHS.TYPE == "ALL"), 
                               aes(x = as.numeric(M.YEAR), y = PRED, group = PATH.GROUP),
                               size = 1, alpha = 0.15, position = position_dodge(0.25)) + 
                     # main points
                     geom_point(size = 3, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1.5, width = 0.2, position = position_dodge(0.5))) +
      plot_layout(guides = "collect") +
      plot_annotation(title = plot_title,
                      subtitle = plot_subtitle) 
    
    return(model_plot)
    
  } else {print("Please select valid analysis type!")}
  
}

### iterative code for overall bird reporting patterns -----------------

b01_overall_monthly <- function(state_name) {
  
  # anal_name, city_list and species_list need to be in environment
  # sliced and non-sliced data as well
  
  cur_city_list <- city_list %>% filter(STATE == state_name)
  cur_species_list <- species_list %>% filter(STATE == state_name)
  
  assign("cur_city_list", cur_city_list, envir = .GlobalEnv)
  assign("cur_species_list", cur_species_list, envir = .GlobalEnv)
  
  
  data_a <- data0_MY_b_slice_G %>% 
    filter(COUNTY == unique(cur_city_list$COUNTY), 
           CELL.ID %in% cur_city_list$CELLS) %>% 
    group_by(M.YEAR, MONTH, CELL.ID) %>% 
    dplyr::summarise(GROUP.ID = GROUP.ID, # not mutating so as to remove unwanted join columns
              TOT.LISTS = n_distinct(GROUP.ID)) %>% 
    ungroup() %>% 
    left_join(data0_MY_b) %>% 
    filter(COMMON.NAME %in% cur_species_list$COMMON.NAME) %>% 
    group_by(COMMON.NAME, M.YEAR, MONTH, CELL.ID) %>% 
    dplyr::summarise(TOT.LISTS = min(TOT.LISTS),
              NO.LISTS = n_distinct(GROUP.ID),
              REP.FREQ = round(NO.LISTS/TOT.LISTS, 4)) %>% 
    dplyr::summarise(REP.FREQ = boot_conf(REP.FREQ)) %>% 
    group_by(COMMON.NAME, M.YEAR, MONTH) %>% 
    dplyr::summarise(CI.L = stats::quantile(REP.FREQ, 0.025), # Obtain the CIs
              CI.U = stats::quantile(REP.FREQ, 0.975),
              REP.FREQ = median(REP.FREQ)) %>% 
    left_join(cur_species_list) %>% 
    dplyr::select(-STATE)
  
  print("data_a completed.")
  
  data_b <- data0_MY_b_slice_G %>% 
    filter(STATE == unique(cur_city_list$STATE)) %>% 
    group_by(M.YEAR, MONTH, CELL.ID) %>% 
    dplyr::summarise(GROUP.ID = GROUP.ID, # not mutating so as to remove unwanted join columns
              TOT.LISTS = n_distinct(GROUP.ID)) %>% 
    ungroup() %>% 
    left_join(data0_MY_b) %>% 
    filter(COMMON.NAME %in% cur_species_list$COMMON.NAME) %>% 
    group_by(COMMON.NAME, M.YEAR, MONTH, CELL.ID) %>% 
    dplyr::summarise(TOT.LISTS = min(TOT.LISTS),
              NO.LISTS = n_distinct(GROUP.ID),
              REP.FREQ = round(NO.LISTS/TOT.LISTS, 4)) %>% 
    dplyr::summarise(REP.FREQ = boot_conf(REP.FREQ)) %>% 
    group_by(COMMON.NAME, M.YEAR, MONTH) %>% 
    dplyr::summarise(CI.L = stats::quantile(REP.FREQ, 0.025), # Obtain the CIs
              CI.U = stats::quantile(REP.FREQ, 0.975),
              REP.FREQ = median(REP.FREQ)) %>% 
    left_join(cur_species_list) %>% 
    dplyr::select(-STATE)
  
  print("data_b completed.")
  
  assign("data_a", data_a, envir = .GlobalEnv)
  assign("data_b", data_b, envir = .GlobalEnv)
  
  
  plot_a <- gg_b_nonmodel(data_a, region = "city", time = "monthly")
  plot_b <- gg_b_nonmodel(data_b, region = "state", time = "monthly")
  
  assign("plot_a", plot_a, envir = .GlobalEnv)
  assign("plot_b", plot_b, envir = .GlobalEnv)
  
  
  ggsave(filename = glue("03_wrap_figs/{anal_name}_a.png"), plot = plot_a,
         dpi = 300, width = 16, height = 15, units = "in")
  
  ggsave(filename = glue("03_wrap_figs/{anal_name}_b.png"), plot = plot_b,
         dpi = 300, width = 16, height = 15, units = "in")
  
  print("plot_a and plot_b created and written to disk.")
  
}


b01_overall_annual <- function(state_name) {
  
  # previous needs to be run before this for data objects!
  
  data_a <- data_a %>% 
    mutate(SE = (REP.FREQ - CI.L)/1.96) %>% 
    group_by(COMMON.NAME, M.YEAR) %>% 
    dplyr::summarise(REP.FREQ = mean(REP.FREQ),
              SE = sqrt(sum((SE)^2))/n(),
              CI.L = REP.FREQ - 1.96*SE,
              CI.U = REP.FREQ + 1.96*SE) %>% 
    left_join(cur_species_list) %>% 
    dplyr::select(-STATE)
  
  print("data_a completed.")
  
  data_b <- data_b %>% 
    mutate(SE = (REP.FREQ - CI.L)/1.96) %>% 
    group_by(COMMON.NAME, M.YEAR) %>% 
    dplyr::summarise(REP.FREQ = mean(REP.FREQ),
              SE = sqrt(sum((SE)^2))/n(),
              CI.L = REP.FREQ - 1.96*SE,
              CI.U = REP.FREQ + 1.96*SE) %>% 
    left_join(cur_species_list) %>% 
    dplyr::select(-STATE)
  
  print("data_b completed.")
  
  assign("data_a", data_a, envir = .GlobalEnv)
  assign("data_b", data_b, envir = .GlobalEnv)
  
  
  plot_a <- gg_b_nonmodel(data_a, region = "city", time = "yearly")
  plot_b <- gg_b_nonmodel(data_b, region = "state", time = "yearly")
  
  assign("plot_a", plot_a, envir = .GlobalEnv)
  assign("plot_b", plot_b, envir = .GlobalEnv)
  
  
  ggsave(filename = glue("03_wrap_figs/{anal_name}_a.png"), plot = plot_a,
         dpi = 300, width = 10, height = 12, units = "in")
  
  ggsave(filename = glue("03_wrap_figs/{anal_name}_b.png"), plot = plot_b,
         dpi = 300, width = 10, height = 12, units = "in")
  
  print("plot_a and plot_b created and written to disk.")
  
}

