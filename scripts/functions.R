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
  
# MCD12Q1 LULC data from December 2019, being aggregated to 2kmx2km ("UNU")
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
# latest Ashwin's maps.RData file (having areas)
# path to list of group accounts to be filtered out
# path to classification of year-month as COVID categories
# path to raster data at 2x2 scale


data_qualfilt_prep <- function(datapath, groupaccspath, covidclasspath,
                               rast_UNU_path,
                               maxvel = 20, minsut = 2){
  
  require(tictoc)
  tic("Completed data quality filtering and preparation in")
  
  ### importing from usual modified ebd RData
  load(datapath) 
  # adding migratory year column
  data <- data %>% 
    mutate(M.YEAR = if_else(MONTH > 5, YEAR, YEAR-1), # from June to May
           M.MONTH = if_else(MONTH > 5, MONTH-5, 12-(5-MONTH))) 
    
  
  ### list of group accounts to be filtered
  groupaccs <- read_csv(groupaccspath) %>% 
    mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                                GA.2 == 1 ~ "GA.2", 
                                TRUE ~ "NG"))
  filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)
  
  
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
    select(-LONGITUDE, -LATITUDE)
  
  save(lists_UNU, file = "data/lists_UNU.RData")
  
  print("added UNU information")
  
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
           MONTH = month(LAST.EDITED.DATE),
           M.YEAR = if_else(MONTH > 5, YEAR, YEAR-1), # from June to May
           M.MONTH = if_else(MONTH > 5, MONTH-5, 12-(5-MONTH))) %>%
    filter(M.YEAR >= 2018) %>%
    left_join(covidclass) %>%
    mutate(COVID = factor(COVID,
                          levels = c("BEF","DUR_20","DUR_21","AFT"))) %>%
    rename(LE.YEAR = M.YEAR,
           LE.MONTH = M.MONTH) %>%
    mutate(LE.YEAR = factor(LE.YEAR, levels = seq(2018, 2021, by = 1)),
           LE.MONTH = factor(LE.MONTH, levels = seq(1, 12, by = 1)),
           YEAR = factor(year(OBSERVATION.DATE), levels = seq(2018, 2022, by = 1)),
           MONTH = factor(month(OBSERVATION.DATE), levels = seq(1, 12, by = 1))) %>% 
    mutate(M.YEAR = if_else(MONTH > 5, YEAR, YEAR-1), # from June to May
           M.MONTH = if_else(MONTH > 5, MONTH-5, 12-(5-MONTH)))

  # filtering
  new_obsr_data <- new_obsr_data %>% anti_join(filtGA)
  save(new_obsr_data, file = "data/new_obsr_data.RData")
  
  print("obtained new observer data")
  
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
  
  
  # making month and year ordered factors
  data0_MY <- data0_MY %>% 
    mutate(M.YEAR = factor(M.YEAR, levels = seq(2018, 2021, by = 1)),
           MONTH = factor(MONTH, levels = seq(1, 12, by = 1)),
           M.MONTH = factor(M.MONTH, levels = seq(1, 12, by = 1)))
  
  
  # sliced data which is what is required for analyses
  data0_MY_slice_S <- data0_MY %>% 
    group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
    slice(1) %>% 
    ungroup()
  data0_MY_slice_G <- data0_MY %>% 
    group_by(GROUP.ID) %>% 
    slice(1) %>%
    ungroup()
  
  # adding map variables (CELL.ID) to main data
  load("data/maps.RData") # Ashwin's maps data
  
  lists_grids <- data0_MY_slice_G %>% 
    distinct(GROUP.ID, LONGITUDE, LATITUDE) %>% 
    joinmapvars() %>% 
    # 25x25 grid cells
    rename(CELL.ID = GRIDG1)

  data0_MY <- data0_MY %>% left_join(lists_grids)
  data0_MY_slice_S <- data0_MY_slice_S %>% left_join(lists_grids)
  data0_MY_slice_G <- data0_MY_slice_G %>% left_join(lists_grids)
  
  save(data0_MY, file = "data/data0_MY.RData")
  save(data0_MY_slice_S, data0_MY_slice_G, file = "data/data0_MY_slice.RData")
  
  print("completed main data filtering")
  
  
  
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


### local moran cluster and outlier analysis (COA) -----------------

# Deriving the cluster/outlier types (COType in ArcGIS term) for each spatial feature 

# manual implementation of COA using local Moran's I (Anselin), because no function 
# exists for this

# ARGUMENTS:
# data is the input spatial object (sf dataframe)
# data_weights is the list of weights
# sig.lvl <- 0.05 # 95% confidence
# correction <- "fdr" # False Discovery Rate, performs better in COA; other options in p.adjust()

localmoran_coa <- function(data, data_weights, 
                           sig.lvl = 0.05, correction = "fdr", 
                           highlight = "outlier", trans.resp = F){
  
  if (highlight == "outlier") {
    CO_levels <- c("HL", "HH", "NS", "LL", "LH")
    } else if (highlight == "cluster") {
      CO_levels <- c("HH", "HL", "NS", "LH", "LL")
    }
  
  if (trans.resp == T) {
    # log-transforming the number of lists (+1) in order to allow analysis to be more sensitive,
    # and to not let the high values overwhelm
    data <- data %>% mutate(NO.LISTS = log(NO.LISTS + 1))
    
    print("Response variable transformed:")
    print(summary(data$NO.LISTS))
  } else {
    print("Response variable not transformed.")
  }
  # is response variable normally distributed?
  hist(data$NO.LISTS)
  
  
  # calculating mean value because in cluster and outlier analysis, the reference to 
  # high and low is relative to the mean of the variable, and should not be interpreted 
  # in an absolute sense.
  mean.ref <- mean(data$NO.LISTS)
  print(paste("Mean reference value:", mean.ref))
  
  data_COA <- localmoran(data$NO.LISTS, data_weights) %>% 
    as_tibble() %>% 
    magrittr::set_colnames(c("Ii","E.Ii","Var.Ii","Z.Ii","Pr(z > 0)")) # for easy reference
  
  # adjusting the p-value based on correction method
  data_COA$P.ADJ <- p.adjust(data_COA$`Pr(z > 0)`, method = correction)
  
  data_COA <- data_COA %>% 
    # joining to main data
    bind_cols(data) %>% 
    mutate(CO.TYPE = factor(
      case_when(P.ADJ > sig.lvl ~ "NS",
                P.ADJ <= sig.lvl & Ii >= 0 & NO.LISTS >= mean.ref ~ "HH",
                P.ADJ <= sig.lvl & Ii >= 0 & NO.LISTS < mean.ref ~ "LL",
                P.ADJ <= sig.lvl & Ii < 0 & NO.LISTS >= mean.ref ~ "HL",
                P.ADJ <= sig.lvl & Ii < 0 & NO.LISTS < mean.ref ~ "LH"),
      levels = CO_levels
    )) %>% 
    # renaming Ii to MORAN for easy reference
    rename(MORAN = Ii) %>% 
    # reclassifying from "localmoran" type to double
    mutate(across(c("MORAN","E.Ii","Var.Ii","Z.Ii","Pr(z > 0)"), ~ as.double(.x)))
  
  return(data_COA)
  
}


### period-wise spatial spread analysis (local Moran cluster and outlier analysis) -------

# need to input three different data objects corresponding to three periods

pw_spread_LMCOA <- function(data_p1, data_p2, data_p3, trans.resp = F) {
  
  require(spdep)

  
  {
    clust_data1 <- data_p1 %>% 
      mutate(COVID = factor(COVID, levels = c("BEF", "DUR", "AFT"))) 
    
    clust_w1 <- clust_data1 %>% poly2omit0nb() %>% nb2listw()
    
    # removing 0nb cells from main data obj also
    clust_data1 <- clust_data1 %>% anti_join(cell_zero)
    
    clust_data1 <- localmoran_coa(clust_data1, clust_w1, 
                                  sig.lvl = 0.05, correction = "fdr", trans.resp = trans.resp)

  }
  
  {
    clust_data2 <- data_p2 %>% 
      mutate(COVID = factor(COVID, levels = c("BEF", "DUR", "AFT"))) 
    
    clust_w2 <- clust_data2 %>% poly2omit0nb() %>% nb2listw()
    
    # removing 0nb cells from main data obj also
    clust_data2 <- clust_data2 %>% anti_join(cell_zero)
    
    clust_data2 <- localmoran_coa(clust_data2, clust_w2, 
                                  sig.lvl = 0.05, correction = "fdr", trans.resp = trans.resp)
  }
  
  {
    clust_data3 <- data_p3 %>% 
      mutate(COVID = factor(COVID, levels = c("BEF", "DUR", "AFT"))) 
    
    clust_w3 <- clust_data3 %>% poly2omit0nb() %>% nb2listw()
    
    # removing 0nb cells from main data obj also
    clust_data3 <- clust_data3 %>% anti_join(cell_zero)
    
    clust_data3 <- localmoran_coa(clust_data3, clust_w3, 
                                  sig.lvl = 0.05, correction = "fdr", trans.resp = trans.resp)
  }
  
  
  clust_data <- clust_data1 %>% 
    bind_rows(clust_data2) %>% 
    bind_rows(clust_data3) %>% 
    # making it sf object
    st_as_sf()
  
  
  # should the weights objects be saved in environment? don't think so
  
  
  # returning overall data with moran values
  return(clust_data)
  
}

### state-level modelling of bird repfreq -----------------

bird_model_state <- function(data_full = data0_MY, 
                             data_sliceG = data0_MY_slice_G, 
                             state) {
  
  require(lme4)
  
  species <- species_list %>% filter(STATE == state)

  ### modelling directly with presence-absence data instead of relative abundance
  
  # to join for presence-absence of various species
  temp1 <- data_full %>% 
    filter(STATE == state) %>% 
    group_by(GROUP.ID, COMMON.NAME) %>% 
    summarise(OBSERVATION.COUNT = max(OBSERVATION.COUNT)) %>% 
    ungroup()
  
  # to later join checklist metadata
  temp2 <- data_full %>% 
    filter(STATE == state) %>% 
    arrange(SAMPLING.EVENT.IDENTIFIER) %>% 
    group_by(GROUP.ID) %>% 
    slice(1) %>% ungroup() %>% 
    select(GROUP.ID, STATE, COUNTY, LOCALITY, LATITUDE, LONGITUDE, OBSERVATION.DATE, 
           M.YEAR, MONTH, DAY.M, M.YEAR, URBAN, CELL.ID, SUBCELL.ID, NO.SP)
  
  data_occ <- data_sliceG %>% 
    filter(STATE == state) %>% 
    group_by(GROUP.ID) %>% 
    summarise(COMMON.NAME = species$COMMON.NAME) %>% 
    left_join(temp1) %>% 
    # for species not reported in lists, filling in NAs in COMMON.NAME and REPORT
    mutate(REPORT = replace_na(OBSERVATION.COUNT, "0")) %>% 
    select(-OBSERVATION.COUNT) %>% 
    mutate(REPORT = as.numeric(case_when(REPORT != "0" ~ "1", TRUE ~ REPORT))) %>% 
    # checklist metadata
    left_join(temp2, by = "GROUP.ID") %>% 
    arrange(GROUP.ID) %>% 
    # species categories
    left_join(species) %>% 
    ungroup()
  
  data_occ0 <- bind_rows("LD" = data_occ %>% filter(MONTH %in% 4:5), 
                         "NL" = data_occ %>% filter(!(MONTH %in% 4:5)), 
                         "ALL" = data_occ, 
                         .id = "MONTHS.TYPE") 
  
  # getting median list length for prediction later
  median_length <- data_occ0 %>% 
    distinct(MONTHS.TYPE, M.YEAR, MONTH, GROUP.ID, NO.SP) %>% 
    group_by(MONTHS.TYPE, MONTH) %>% 
    summarise(NO.SP.MED = floor(median(NO.SP)))

  # dataframe with empty column to populate with looped values
  # total rows: product of distinct values of predictors
  birds_pred <- data_occ0 %>% 
    group_by(MONTHS.TYPE) %>% 
    tidyr::expand(COMMON.NAME, nesting(MONTH), M.YEAR) %>% 
    left_join(species) %>% 
    mutate(REP.FREQ.PRED = NA)
  
  print("Completed preparations for modelling. Now starting modelling.")
  
  
  count <- 0
  for (m in 1:n_distinct(birds_pred$MONTHS.TYPE)) {
    
    data_mtype <- data_occ0 %>% 
      filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[m])
    
    birds_pred0 <- birds_pred %>% 
      filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[m]) %>% 
      rename(REP.FREQ.PRED2 = REP.FREQ.PRED)

    for (i in 1:n_distinct(birds_pred0$COMMON.NAME)) {
      
      data_spec <- data_mtype %>% 
        filter(COMMON.NAME == unique(birds_pred0$COMMON.NAME)[i]) %>% 
        # filtering for CELL.ID-MONTH (space-time) combos in which species occurs
        filter(REPORT == 1) %>% 
        distinct(COMMON.NAME, CELL.ID, MONTH) %>% 
        left_join(data_occ)
      
      
      model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:NO.SP + MONTH:M.YEAR + 
                            (1|CELL.ID),
                          data = data_spec, family = binomial(link = "cloglog"),
                          nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
      
      for (j in 1:n_distinct(birds_pred0$MONTH)) {
        
        for (k in 1:n_distinct(birds_pred0$M.YEAR)) {
          count <- count + 1
          
          birds_pred0$REP.FREQ.PRED2[count] = predict(
            model_spec,
            data.frame(COMMON.NAME = birds_pred0$COMMON.NAME[count],
                       MONTH = birds_pred0$MONTH[count],
                       M.YEAR = birds_pred0$M.YEAR[count],
                       NO.SP = median_length$NO.SP.MED[birds_pred0$MONTH[count]]),
            re.form = NA, 
            type = "response")
          
          print(glue::glue("Months type {m}: Completed predictions for {birds_pred0$COMMON.NAME[count]} in
                     month {birds_pred0$MONTH[count]}, migratory year {birds_pred0$M.YEAR[count]}"))
        }
      }
    } 
    
    birds_pred <- birds_pred %>% 
      left_join(birds_pred0) %>% 
      mutate(REP.FREQ.PRED = coalesce(REP.FREQ.PRED, REP.FREQ.PRED2)) %>% 
      select(-REP.FREQ.PRED2)
    
    # resetting counter for next MONTHS.TYPE
    count <- 0
  }
  
  birds_pred <- birds_pred %>% 
    group_by(MONTHS.TYPE, M.YEAR, SP.CATEGORY) %>% 
    summarise(SE = sd(REP.FREQ.PRED)/sqrt(n()),
              REP.FREQ.PRED = mean(REP.FREQ.PRED),
              CI.L = REP.FREQ.PRED - 1.96*SE,
              CI.U = REP.FREQ.PRED + 1.96*SE)


  (ggplot(filter(birds_pred, MONTHS.TYPE == "LD"), 
          aes(M.YEAR, REP.FREQ.PRED, col = SP.CATEGORY)) +
      scale_color_manual(values = c("#8F85C1", "#A3383C"),
                         name = "Species category",
                         labels = c("Rural", "Urban")) +
      scale_y_continuous(limits = c(0.1, 1.0)) +
      labs(title = "For the months of April and May",
           x = "Migratory year", y = "Predicted reporting frequency") +
      geom_point(size = 3, position = position_dodge(0.5)) +
      geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                    size = 1.5, width = 0.2, position = position_dodge(0.5)) |
    ggplot(filter(birds_pred, MONTHS.TYPE == "NL"), 
             aes(M.YEAR, REP.FREQ.PRED, col = SP.CATEGORY)) +
      scale_color_manual(values = c("#8F85C1", "#A3383C"),
                         name = "Species category",
                         labels = c("Rural", "Urban")) +
      scale_y_continuous(limits = c(0.1, 1.0)) +
      labs(title = "For other ten months",
           x = "Migratory year", y = "Predicted reporting frequency") +
      geom_point(size = 3, position = position_dodge(0.5)) +
      geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                    size = 1.5, width = 0.2, position = position_dodge(0.5)) |
    ggplot(filter(birds_pred, MONTHS.TYPE == "ALL"), 
             aes(M.YEAR, REP.FREQ.PRED, col = SP.CATEGORY)) +
      scale_color_manual(values = c("#8F85C1", "#A3383C"),
                         name = "Species category",
                         labels = c("Rural", "Urban")) +
      scale_y_continuous(limits = c(0.1, 1.0)) +
      labs(title = "For all twelve months",
           x = "Migratory year", y = "Predicted reporting frequency") +
      geom_point(size = 3, position = position_dodge(0.5)) +
      geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                    size = 1.5, width = 0.2, position = position_dodge(0.5))) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste(state, "state"),
                    subtitle = paste0(
                      "Predicted reporting frequencies of ",
                      n_distinct(data_occ$COMMON.NAME), " species (",
                      n_distinct(filter(data_occ, SP.CATEGORY == "U")$COMMON.NAME), " urban, ",
                      n_distinct(filter(data_occ, SP.CATEGORY == "R")$COMMON.NAME), " rural)",
                      " in three separate models")) -> birds_graph
  
  assign("birds_pred", birds_pred, envir = .GlobalEnv)
  assign("birds_graph", birds_graph, envir = .GlobalEnv)
  
  }
