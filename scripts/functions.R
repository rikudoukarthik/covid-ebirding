### MODIS land-use map to get urban-non-urban information -----------------

getmodisdata <- function(){
  
# MCD12Q1 LULC data from December 2019, being aggregated to 2kmx2km ("UNU") and 
# 24kmx24km (SoIB) scales.
# this will be retained as a raster and also joined to the EBD data

require(raster) # masks dplyr functions, so have used package::function()
require(terra)
require(gdalUtils)
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
start <- "2019-12-31" # time period of mapping
end <- "2019-12-31" # time period of mapping


# creating folder to store MODIS data that will be downloaded
MODISpath <- "data/in_LULC_MODIS/"
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

# choose the band (SDS) that you want using sds[] option and write GTiff files.
# University of Maryland SDS (LC_Type2) is SDS1
for (i in (modis)) {
  sds <- get_subdatasets(i)
  modis2 <- gdal_translate(sds[1], dst_dataset = paste0(i, ".tif"))
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


## projecting from MODIS sinusoidal to WGS84 (flat 2D)
projto <- raster::crs(india)
# next step takes a few minutes
# "near" because LC values are categories despite being numerical
rast_proj <- terra::project(rast_full, projto, method = "near")


## cropping and masking the raster to our area of interest
rast_aoi <- terra::crop(rast_proj, india)
rast_aoi <- terra::mask(rast_aoi, india)


raster::writeRaster(rast_aoi, paste0(MODISpath, "/in_LULC_MODIS.tif"), overwrite = TRUE)

rm(prod, start, end, MODISpath, userpass, modis, sds, modis2, rast_list, rast_full, projto, rast_proj, rast_aoi)


### Modifying the MODIS data #######

## read in cropped and masked .tif
rast <- raster::raster("data/in_LULC_MODIS/in_LULC_MODIS.tif")


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

save(rast_UNU, file = "data/rast_UNU.RData")

# aggregating to 24km*24km for SoIB inference
# (here, priority is getting grid cell information, not so much the UNU classification, as that is not directly going to be used; instead will use birding values from them which will be averaged to the SoIB scale)
rast_SoIB <- rast_UNU %>% 
  raster::aggregate(fact = 12, fun = rast_agg_fn)

save(rast_SoIB, file = "data/rast_SoIB.RData")


# 
# below needed only when finalising the correct code for this whole process. no need to repeat. 
# ### Leaflet overlay to verify urban/non- classification #######
# 
# # https://rstudio.github.io/leaflet/raster.html
# # https://rpubs.com/mgei/drivingtimes
# 
# library(leaflet)
# 
# leafcols <- colorFactor(c(NA, "#0C2C84"), values(rast_UNU),
#                         na.color = "#ffba00")
# 
# # leaflet uses different crs so must be projected; method should be made categorical
# 
# leaflet() %>% 
#   addTiles() %>% 
#   addRasterImage(rast_UNU, 
#                  method = "ngb", # nearest neighbour
#                  colors = leafcols,
#                  opacity = 0.7,
#                  maxBytes = 150 * 1024 * 1024) %>% 
#   addLegend(pal = leafcols, values = values(rast_UNU))



}

### data quality filters and preparation -----------------

# data file path
# path to list of group accounts to be filtered out
# path to classification of year-month as COVID categories
# path to raster data at 2x2 and 24x24 scales


data_qualfilt_prep <- function(datapath, groupaccspath, covidclasspath,
                               rast_UNU_path, rast_SoIB_path,
                               maxvel = 20, minsut = 2){
  
  ### importing from usual modified ebd RData
  load(datapath) 
  # adding migratory year column
  data <- data %>% 
    mutate(DAY.Y = yday(OBSERVATION.DATE),
           M.YEAR = if_else(DAY.Y <= 151, YEAR-1, YEAR)) # from 1st June to 31st May
  
  
  ### list of group accounts to be filtered
  groupaccs <- read_csv(groupaccspath) %>% 
    mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                                GA.2 == 1 ~ "GA.2", 
                                TRUE ~ "NG"))
  filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)
  
  
  ### COVID classification
  covidclass <- read_csv(covidclasspath)
  

  ### adding UNU information ######
  
  require(tidyverse)
  require(lubridate)
  require(terra)
  require(sp)
  
  load(rast_UNU_path)
  load(rast_SoIB_path)
  
  # adding extracted land cover values at eBird data points with grid cell information at both scales
  # Pakistan-administered Kashmir lists produce NA for URBAN
  lists_UNU <- data %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    mutate(URBAN = raster::extract(rast_UNU, cbind(LONGITUDE, LATITUDE)),
           NONURBAN = if_else(URBAN == 1, 0, 1),
           # 2km*2km
           SUBCELL.ID = raster::cellFromXY(rast_UNU, cbind(LONGITUDE, LATITUDE)), 
           # 24km*24km
           CELL.ID = raster::cellFromXY(rast_SoIB, cbind(LONGITUDE, LATITUDE))) %>% 
    filter(!is.na(URBAN)) %>% 
    select(-LONGITUDE, -LATITUDE)
  
  save(lists_UNU, file = "data/lists_UNU.RData")
  
  
  ### new observer data (to calculate no. of new observers metric) #######
  
  new_obsr_data <- data %>% 
    select(c("YEAR", "MONTH", "STATE", "SAMPLING.EVENT.IDENTIFIER",
             "LAST.EDITED.DATE", "OBSERVATION.DATE", "OBSERVER.ID")) %>% 
    mutate(LAST.EDITED.DATE = ymd_hms(LAST.EDITED.DATE)) %>% 
    group_by(OBSERVER.ID) %>% 
    arrange(LAST.EDITED.DATE) %>% 
    ungroup() %>% 
    distinct(OBSERVER.ID, .keep_all = TRUE) %>%
    mutate(YEAR = year(LAST.EDITED.DATE),
           MONTH = month(LAST.EDITED.DATE)) %>% 
    filter(YEAR >= 2019) %>% 
    left_join(covidclass) %>% 
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>% 
    rename(LE.YEAR = YEAR,
           LE.MONTH = MONTH) %>% 
    mutate(YEAR = year(OBSERVATION.DATE), 
           MONTH = month(OBSERVATION.DATE))
  
  # filtering
  new_obsr_data <- new_obsr_data %>% anti_join(filtGA) 
  save(new_obsr_data, file = "data/new_obsr_data.RData")
  
  
  
  ### main data filtering ######
  
  # to calculate list length aka number of species
  temp0 <- data %>% 
    mutate(CATEGORY = case_when(CATEGORY == "domestic" & 
                                  COMMON.NAME == "Rock Pigeon" ~ "species",
                                TRUE ~ CATEGORY)) %>% 
    filter(CATEGORY %in% c("issf","species")) %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
    summarise(NO.SP = n_distinct(COMMON.NAME))
  
  data0_MY <- data %>% 
    # this data (including latter months of 2018 only needed for bird behaviour section)
    # so will later save separate data filtering out 2018 months
    filter(M.YEAR >= 2018) %>%
    anti_join(filtGA) %>% # removing data from group accounts
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
  temp1 <- data0_MY %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE != "Incidental") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter((NO.SP <= 3 & is.na(DURATION.MINUTES)) | 
             (DURATION.MINUTES < 3) | 
             (SUT < minsut & NO.SP <= 3)) %>%
    distinct(GROUP.ID)
  
  # getting list of GROUP.IDs inside IN boundary for pelagic filter
  temp2 <- data0_MY %>% 
    ungroup() %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    terra::vect(geom = c("LONGITUDE","LATITUDE"), crs = crs(india)) %>% 
    terra::intersect(india) %>% 
    terra::as.data.frame() %>% 
    select(GROUP.ID) 

  # speed and distance filter for travelling lists
  temp3 <- data0_MY %>%
    filter(ALL.SPECIES.REPORTED == 1 & PROTOCOL.TYPE == "Traveling") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter((SPEED > maxvel) | (EFFORT.DISTANCE.KM > 50)) %>%
    distinct(GROUP.ID)
  
  # true completeness + other filters
  data0_MY <- data0_MY %>%
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
    select(-BREEDING.CODE, -SPEED, -SUT, -MIN, -DATETIME, -HOUR.END, -NOCT.FILTER)
  
  # sliced data which is what is required for analyses
  data0_MY_slice_S <- data0_MY %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
    slice(1) %>% 
    ungroup()
  data0_MY_slice_G <- data0_MY %>% 
    group_by(GROUP.ID) %>% 
    slice(1) %>%
    ungroup()
  
  save(data0_MY, file = "data/data0_MY.RData")
  save(data0_MY_slice_S, data0_MY_slice_G, file = "data/data0_MY_slice.RData")
  
  

  # retaining only data from 2019 onward (for birder behaviour section)
  
  data0 <- data0_MY %>% filter(YEAR >= 2019)
  data0_slice_S <- data0_MY_slice_S %>% filter(YEAR >= 2019)
  data0_slice_G <- data0_MY_slice_G %>% filter(YEAR >= 2019)
  
  # Getting list of months having data from all 3 COVID periods 
  month_compar <- data0_slice_S %>% 
    group_by(MONTH) %>% 
    mutate(COVID = as.character(COVID)) %>% 
    mutate(COVID = factor(case_when(COVID %in% c("DUR_20","DUR_21") ~ "DUR",
                                    TRUE ~ COVID),
                          levels = c("BEF","DUR","AFT"))) %>% 
    dplyr::summarise(N = n_distinct(COVID)) %>% 
    filter(N == 3) %>%  # all 3 COVID categories
    dplyr::select(MONTH)
  
  save(data0, file = "data/data0.RData")
  save(data0_slice_S, data0_slice_G, month_compar, file = "data/data0_slice.RData")

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



### raster aggregation -----------------

# Function to supply in raster::aggregate(), with 25% threshold for classification
# (works for any binary classification to be aggregated; here it is urban/non-)
# https://stackoverflow.com/questions/23369579/aggregate-raster-in-r-with-na-values 

rast_agg_fn <- function(x, ...){
  if_else(sum(x %in% 1) >= 0.25*length(x), 1, 0)
}

# x = c(rep(1, 2), rep(0, 8), rep(1, 10), rep(0, 30)) # for testing


### raster calc: log proportional change in number of lists -----------------

rast_logpropchange <- function(x, y, k = 1)  {
  
  # x is first, y is second, so change from x to y
  
  # # when no. of lists has been log transformed
  # # log(n2) - log(n1) = log(n2/n1)
  # return(y-x)
  
  # when value has not been transformed to prevent NA
  x <- x + k
  y <- y + k
  
  # proportion of urban birding in the cell changed by y/x [0,+n] times.
  return(log10(y/x))
  
}

### raster calc: proportional change in number of lists -----------------

rast_propchange <- function(x, y, k = 1)  {
  
  # x is first, y is second, so change from x to y
  
  # adding 1 to all rows to prevent NA (from zeroes)
  x <- x + k
  y <- y + k
  
  # proportion of urban birding in the cell changed by y/x [0,+n] times.
  return(y/x)
  
}
