# map EBD admin units to sf admin units ---------------------------------------------

# main use is to resolve issue of same CELL.ID or COUNTY.CODE getting mapped to multiple
# sf states/districts
# we use similarity to limit this mapping to 1:1


# function to calc similarity of two strings (https://stackoverflow.com/a/11535768/13000254)
str_similarity <- function(x, y) {
  
  # Levenshtein edit distance
  str_dist <- adist(x, y) %>% as.vector()
  str_longer_length <- max(str_length(c(x, y)))
  
  str_simil <- 1 - (str_dist/str_longer_length)
  return(str_simil)
  
}

map_admin_ebd_sf <- function(data) {

  mapping <- data %>% 
    group_by(STATE, COUNTY.CODE, COUNTY, STATE.NAME, DISTRICT.NAME) %>% 
    filter(!is.na(COUNTY.CODE)) %>% 
    reframe(NO.LISTS = n_distinct(GROUP.ID)) %>% 
    rowwise() %>% 
    mutate(NAME.SIMILARITY = str_similarity(COUNTY, DISTRICT.NAME)) %>% 
    arrange(desc(NAME.SIMILARITY), desc(NO.LISTS)) %>% 
    group_by(STATE, COUNTY.CODE, COUNTY) %>% 
    slice_head(n = 1) %>%
    ungroup() %>% 
    dplyr::select(-c(NAME.SIMILARITY, NO.LISTS)) %>% 
    arrange(STATE, COUNTY.CODE)
  
  return(mapping)
  
}


# calculate cell dim of extremes ----------------------------------------------------

calc_celldim_extremes <- function() {
  
  load("../../00_data/maps_sf.RData")
  load("../../00_data/grids_sf_full.RData")
  
  g1_dims <- map(g1_in_sf$GEOM.G1, ~ st_bbox(.x) %>% 
                   as.matrix() %>% 
                   t() %>% 
                   as.data.frame()) %>% 
    list_rbind() %>% 
    bind_cols(g1_in_sf %>% dplyr::select(GRID.G1) %>% st_drop_geometry())
  
  g1_extremes <- g1_dims %>% 
    mutate(LIMITS = case_when(
      xmin == min(xmin) ~ "West",
      xmax == max(xmax) ~ "East",
      ymin == min(ymin) ~ "South",
      ymax == max(ymax) ~ "North",
      TRUE ~ NA
    )) %>% 
    filter(!is.na(LIMITS))
  
  # we need the full cells for dimensions (g1_in_sf is clipped at borders)
  # so use g1_sf
  g1_sel <- map(g1_sf$GEOM.G1, ~ st_bbox(.x) %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame()) %>% 
    list_rbind() %>% 
    bind_cols(g1_sf %>% dplyr::select(GRID.G1) %>% st_drop_geometry()) %>% 
    filter(GRID.G1 %in% g1_extremes$GRID.G1) %>% 
    mutate(X1 = map2(xmin, ymin, ~ st_point(c(.x, .y)) %>% 
                       st_sfc(crs = st_crs(g1_in_sf)) %>% 
                       st_sf() %>% 
                       rename(X1 = geometry)) %>% 
             list_rbind(),
           X2 = map2(xmax, ymin, ~ st_point(c(.x, .y)) %>% 
                       st_sfc(crs = st_crs(g1_in_sf)) %>% 
                       st_sf() %>% 
                       rename(X2 = geometry)) %>% 
             list_rbind(),
           Y1 = map2(xmin, ymin, ~ st_point(c(.x, .y)) %>% 
                       st_sfc(crs = st_crs(g1_in_sf)) %>% 
                       st_sf() %>% 
                       rename(Y1 = geometry)) %>% 
             list_rbind(),
           Y2 = map2(xmin, ymax, ~ st_point(c(.x, .y)) %>% 
                       st_sfc(crs = st_crs(g1_in_sf)) %>% 
                       st_sf() %>% 
                       rename(Y2 = geometry)) %>% 
             list_rbind())
  
  
  temp_fn <- function(a, b, c, ...) {
    x <- tibble(DIST = c(a, b),
                GRID.G1 = rep(c, 2))
    return(x)
  }
  
  x <- pmap(list(a = g1_sel$X1$X1, b = g1_sel$X2$X2, c = g1_sel$GRID.G1), 
            temp_fn) %>% 
    list_rbind() %>% 
    rowwise() %>% 
    reframe(st_point(DIST) %>% 
              st_sfc(crs = st_crs(g1_in_sf)) %>% 
              st_sf(),
            GRID.G1 = GRID.G1) %>% 
    st_as_sf(crs = st_crs(g1_in_sf)) %>% 
    mutate(GRID.G1 = as.numeric(GRID.G1)) %>% 
    group_by(GRID.G1) %>% 
    reframe(X = st_distance(.)[2*cur_group_id() - 1, 2*cur_group_id()]) %>% 
    mutate(GRID.G1 = as.character(GRID.G1)) 
  
  y <- pmap(list(a = g1_sel$Y1$Y1, b = g1_sel$Y2$Y2, c = g1_sel$GRID.G1), 
            temp_fn) %>% 
    list_rbind() %>% 
    rowwise() %>% 
    reframe(st_point(DIST) %>% 
              st_sfc(crs = st_crs(g1_in_sf)) %>% 
              st_sf(),
            GRID.G1 = GRID.G1) %>% 
    st_as_sf(crs = st_crs(g1_in_sf)) %>% 
    mutate(GRID.G1 = as.numeric(GRID.G1)) %>% 
    group_by(GRID.G1) %>% 
    reframe(Y = st_distance(.)[2*cur_group_id() - 1, 2*cur_group_id()]) %>% 
    mutate(GRID.G1 = as.character(GRID.G1)) 
  
  g1_celldim_extremes <- g1_extremes %>% 
    dplyr::select(GRID.G1, LIMITS) %>% 
    left_join(x, by = "GRID.G1") %>% 
    left_join(y, by = "GRID.G1") %>% 
    mutate(across(c(X, Y),
                  ~ round(units::set_units(., "km"), 1)))
  
  return(g1_celldim_extremes)
  
}


# get labelled months ---------------------------------------------------------------

# and COVID colour scheme

get_labelled_months <- function() {

  month_lab <- data.frame(M.YEAR = 2018:2021) %>% 
    group_by(M.YEAR) %>% 
    summarise(MONTH = 1:12 %>% month()) %>% 
    mutate(MONTH.LAB = month(MONTH, label = T),
           M.YEAR.LAB = glue("{M.YEAR}\u2013{as.character(M.YEAR+1) %>% str_sub(3)}")) %>% 
    mutate(MONTH = factor(MONTH, levels = month(c(6:12, 1:5))),
           M.YEAR = factor(M.YEAR),
           MONTH.LAB = factor(MONTH.LAB, levels = month(c(6:12, 1:5), label = T)),
           M.YEAR.LAB = factor(M.YEAR.LAB, levels = c("2018\u201319", "2019\u201320",
                                                      "2020\u201321", "2021\u201322"))) %>% 
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


### map functions ########################################

# getting mapping functions
source(url("https://raw.githubusercontent.com/birdcountindia/bci-functions/main/01_functions/mapping.R"))


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
# maps_sf.RData file (having areas)
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
  sf_use_s2(FALSE) # not spherical geometry
  
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
  
  ### new observer data (to calculate no. of new observers metric---NO LONGER USED) #######

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
           HOUR.END = floor((HOUR*60 + MIN + DURATION.MINUTES)/60)) %>%
    # data quality
    filter(REVIEWED == 0 | APPROVED == 1)
  
  
  ### remove probable mistakes
  source(url("https://github.com/stateofindiasbirds/soib_2023/raw/master/00_scripts/rm_prob_mistakes.R"))
  data0_MY_b <- rm_prob_mistakes(data0_MY_b)

  
  ### exclude records based on various criteria 
  
  # false complete lists (without duration info & 3 or fewer species, <3min, low SUT)
  temp1 <- data0_MY_b %>%
    # only complete
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE != "Incidental") %>%
    # false complete
    distinct(GROUP.ID, .keep_all = TRUE) %>% 
    filter((NO.SP <= 3 & is.na(DURATION.MINUTES)) | 
             (DURATION.MINUTES < 3) | 
             (SUT < minsut & NO.SP <= 3))
  
  # getting list of GROUP.IDs inside IN boundary for pelagic filter
  temp2 <- data0_MY_b %>% 
    ungroup() %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    terra::vect(geom = c("LONGITUDE","LATITUDE"), crs = crs(india_sf)) %>% 
    terra::intersect(terra::vect(india_sf)) %>% 
    terra::as.data.frame() %>% 
    dplyr::select(GROUP.ID) 

  # speed and distance filter for travelling lists
  temp3 <- data0_MY_b %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE == "Traveling") %>%
    distinct(GROUP.ID, .keep_all = TRUE) %>%
    filter((SPEED > maxvel) | (EFFORT.DISTANCE.KM > 10)) 
  
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
    distinct(SAMPLING.EVENT.IDENTIFIER, .keep_all = TRUE)
  
  data0_MY_b_slice_G <- data0_MY_b %>% 
    distinct(GROUP.ID, .keep_all = TRUE)
  

  lists_grids <- data0_MY_b_slice_G %>% 
    join_map_sf() %>% 
    # 25x25 grid cells
    rename(CELL.ID = GRID.G1) %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE, CELL.ID, GRID.G3)
  
  data0_MY_b <- data0_MY_b %>% left_join(lists_grids)
  data0_MY_b_slice_S <- data0_MY_b_slice_S %>% left_join(lists_grids)
  data0_MY_b_slice_G <- data0_MY_b_slice_G %>% left_join(lists_grids)
  
  save(data0_MY_b, file = "00_data/data0_MY_b.RData")
  save(data0_MY_b_slice_S, data0_MY_b_slice_G, file = "00_data/data0_MY_b_slice.RData")
  
  print("Complete main data filtering for bird analyses!")
  
  
  ## for birder analyses (group acc filter different) ####
  
  data0_MY_d <- data %>% 
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
    # only complete
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE != "Incidental") %>%
    # false complete
    group_by(GROUP.ID) %>% 
    slice(1) %>%
    filter((NO.SP <= 3 & is.na(DURATION.MINUTES)) | 
             (DURATION.MINUTES < 3) | 
             (SUT < minsut & NO.SP <= 3)) %>%
    distinct(GROUP.ID)
  
  # getting list of GROUP.IDs inside IN boundary for pelagic filter
  temp2 <- data0_MY_d %>% 
    ungroup() %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    terra::vect(geom = c("LONGITUDE","LATITUDE"), crs = crs(india_sf)) %>% 
    terra::intersect(terra::vect(india_sf)) %>% 
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
  
  set.seed(10)
  data0_MY_d_slice_G <- data0_MY_d %>% 
    group_by(GROUP.ID) %>% 
    slice_sample(n = 1) %>%
    ungroup()

  lists_grids <- data0_MY_d_slice_G %>% 
    join_map_sf() %>% 
    # 25x25 grid cells
    rename(CELL.ID = GRID.G1) %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE, CELL.ID)
  
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
  
  pred_bootMer <- bootMer(model, FUN = pred_fun, 
                          nsim = nsim, seed = 20231010, 
                          parallel = "snow", use.u = FALSE, type = "parametric", 
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


# map of India and selected states --------------------------------------------------

map_sites <- function() {
  
  require(tidyverse)
  require(sf)
  require(mapview)
  # require(ggmap)
  
  require(rnaturalearth)
  require(rnaturalearthdata)
  require(rmapshaper) # to simplify sf
  require(ggspatial) # for north arrow and legend
  
  
  world_sf <- ne_countries(type = "countries", 
                           scale = "large", returnclass = "sf") %>% 
    # selecting columns of interest
    reframe(NAME = sovereignt, 
            geometry = geometry) %>% 
    st_as_sf()
  
  # replacing default India with Indian India boundaries
  world_sf <- world_sf %>% 
    filter(NAME != "India") %>% 
    # since the row is at the end, will lay over other disputed territories
    bind_rows(india_sf %>% 
                dplyr::select(geometry) %>% 
                mutate(NAME = "India"))
  
  world_sf <- world_sf %>% 
    # # converting to Robinson projection
    # st_transform(crs = "ESRI:54030") %>% 
    # # simplifying polygons
    # ms_simplify(keep = 0.1, keep_shapes = TRUE) %>% 
    # fill colours and line widths
    mutate(FILL = case_when(NAME == "India" ~ "#ffffff",
                                 TRUE ~ "#e5d8ca"),
           LINEWIDTH = case_when(NAME == "India" ~ 0.5,
                                 TRUE ~ 0.35)) %>% 
    # select country labels of interest
    mutate(NAME = case_when(NAME %in% c("India", "Pakistan", "Afghanistan", "Tajikistan",
                                        "China", "Nepal", "Bhutan", "Bangladesh", "Myanmar",
                                        "Sri Lanka") ~ NAME,
                            TRUE ~ NA_character_))
  
  ourstates_sf <- states_sf %>%
    filter(STATE.NAME %in%
             c("Karnataka", "Kerala", "Maharashtra", "Assam")) %>% 
    dplyr::select(-AREA) %>% 
    # # converting to Robinson projection
    # st_transform(crs = "ESRI:54030") %>% 
    # # simplifying polygons
    # ms_simplify(keep = 0.1, keep_shapes = TRUE) %>% 
    # capitalising labels
    mutate(STATE.NAME = str_to_upper(STATE.NAME)) %>% 
    arrange(STATE.NAME)
  
  
  # plot limits 
  plot_lims <- tibble(X = c(60, 60, 100, 100), 
                      Y = c(5, 40, 5, 40)) %>% 
    st_as_sf(coords = c("X", "Y")) %>% 
    dplyr::summarise() %>% 
    st_cast("POLYGON") %>% 
    # important that limits are also projected to Robinson
    st_set_crs(st_crs(india_sf)) %>% 
    # st_transform(crs = "ESRI:54030") %>% 
    st_bbox()
  
  map <- ggplot(world_sf) +
    geom_sf(aes(fill = FILL, linewidth = LINEWIDTH)) +
    geom_sf_text(aes(label = NAME), size = 3) +
    scale_fill_identity() +
    scale_linewidth_identity() + 
    geom_sf(data = ourstates_sf, fill = "grey90", linewidth = 0.55) +
    geom_sf_text(data = ourstates_sf, 
                 aes(label = STATE.NAME), 
                 size = 2,
                 nudge_x = c(2, 0, 0, 0.5),
                 nudge_y = c(0.5, 0, 0, 1)) +
    # need to set coord limits (plot zoom limits)
    coord_sf(xlim = c(plot_lims$xmin, plot_lims$xmax), 
             ylim = c(plot_lims$ymin, plot_lims$ymax)) +
    theme_void() +
    theme(panel.background = element_rect(colour = NA, fill = "#a9d5e0"),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(margin = margin(4, 0, 0, 0)),
          axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
          panel.border = element_rect(fill = NA),
          plot.margin = margin(0, 30, 0, 10), 
          plot.background = element_rect(fill = "white", colour = NA))
  
  
  map_scale <- annotation_scale(location = "bl", 
                                bar_cols = c("black", "white")) 
  
  map_arrow <- annotation_north_arrow(location = "bl", 
                                      which_north = "true",
                                      pad_x = unit(0.4, "in"), 
                                      pad_y = unit(0.4, "in"),
                                      style = north_arrow_fancy_orienteering(fill = c("grey20", "white"),
                                                                   line_col = "black"))
  
  map <- map + map_scale + map_arrow
  
  
  return(map)
  
  
  # map_data("world") %>% 
  #   st_as_sf(coords = c("long", "lat")) %>% 
  #   group_by(group, region, order) %>% 
  #   summarise()
  #   filter(region != "India") %>% 
  #   dplyr::select(region) %>% 
  #   
  #   ggplot() +
  #   geom_sf()
  # 
  # 
  # mapview::mapView(india_sf,
  #                  map.types = "Esri.WorldImagery",
  #                  legend = FALSE,
  #                  viewer.suppress = TRUE) %>% 
  #   mapview::mapshot(file = "test.png")
  
}
