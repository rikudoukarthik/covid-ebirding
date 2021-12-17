# Code for 'Reduced human activity during COVID-19 alters avian land use across North America'

# The first several sections describe how data were prepared from original sources.
# Links to those sources are provided, but several steps in the data acquisition process may differ from our approach
# We recommend contacting the lead author, Michael Schrimpf (michael.schrimpf@umanitoba.ca or michael.schrimpf@gmail.com) if you would like to learn more about how raw data were obtained.

# The data file used in the model itself is produced in the section 'Model data prep' and is also available to download from the University of Manitoba library (see paper for link)
# A description of each column is also providied in the 'Model data prep' section.

# The model section 'Running Model' provides the actual rstanarm regression model



# Mobility Data -----------------------------------------------------------


library(tidyverse)
library(here)


# The mobility data came from here:
# https://www.google.com/covid19/mobility


# raw mobility data
mobility <- here("data", "mobility.csv") %>% 
  read_csv() %>% 
  rename_all(tolower) %>% 
  select(-duplicate, -pop_prop)

# data must have no missing dates or values
n_missing <- mobility %>% 
  arrange(county_code, date) %>% 
  group_by(county_code) %>% 
  mutate(date_diff = c(diff(date), 1)) %>%  
  ungroup() %>% 
  filter(date_diff > 1) %>% 
  nrow()
stopifnot(n_missing == 0)
stopifnot(complete.cases(mobility))

# function to calculate rolling cumulative sum
rolling_cumsum <- function(x, days) {
  cumsum(x) - cumsum(c(rep(0, days), head(x, -days)))
}

# calculate for a set of days back
mobility_metrics <- mobility %>% 
  arrange(county_code, date) %>% 
  select(county_code, date, retail, grocery, parks, transit, work, home)
mobility_cumulative <- mobility %>% 
  distinct(city_num, city, county, state, county, county_code, date)
for (d in c(1, 3, 7, 14, 30)) {
  mobility_cumulative <- mobility_metrics %>% 
    group_by(county_code) %>% 
    mutate_if(is.numeric, rolling_cumsum, days = d) %>% 
    ungroup() %>% 
    rename_if(is.numeric, ~ paste(., d, sep = "_")) %>% 
    inner_join(mobility_cumulative, ., by = c("county_code", "date"))
}

# calculate for the whole period (mean)
mobility_cumulative_full_mean <- mobility_metrics %>%
  filter(date >= "2020-03-01" & date <= "2020-05-31") %>%
  group_by(county_code) %>%
  mutate_at(.vars = c("retail", "grocery", "parks", "transit", "work", "home"), .funs = mean) %>%
  select(-date) %>%
  distinct()


# export
here("data", "mobility_cumulative.csv") %>% 
  write_csv(mobility_cumulative, .)
here("data", "mobility_cumulative_full_mean.csv") %>% 
  write_csv(mobility_cumulative_full_mean, .)




# Extracting eBird data ---------------------------------------------------

# eBird data can be found at www.ebird.org

# Code written by M. Strimas-Mackey


### 01 ebd-extract

library(auk)

# county list
counties <- readRDS("data/county_code_list.rds")

# extract from ebd
f_ebd <- "data/ebd_2020-05.txt"
f_sed <- "data/sed_2020-05.txt"
# columns to retain
keep_cols <- c("sampling_event_identifier", "observer_id", "group_identifier", 
               "category", "scientific_name", "subspecies_scientific_name", 
               "observation_count", 
               "latitude", "longitude", 
               "country_code", "state_code", "county_code", 
               "locality_id", "locality_type", 
               "protocol_type", "protocol_code", 
               "observation_date", "time_observations_started", 
               "duration_minutes", "effort_distance_km", "effort_area_ha", 
               "number_observers", 
               "all_species_reported")
saveRDS(keep_cols, "data/ebd-columns.rds")
if (!file.exists(f_ebd)) {
  auk_ebd("ebd_relMay-2020.txt", "ebd_sampling_relMay-2020.txt") %>% 
    auk_county(county = counties) %>% 
    # march-may, 2017-2020
    auk_year(year = 2017:2020) %>% 
    auk_date(date = c("*-03-01", "*-05-31")) %>% 
    # < 5 h, < 5 km, traveling and stationary counts
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_distance(distance = c(0, 5)) %>% 
    auk_duration(duration = c(0, 5 * 60)) %>% 
    auk_complete() %>% 
    auk_filter(f_ebd, file_sampling = f_sed, keep = keep_cols)
}


### 02 process-counties

library(auk)
library(lubridate)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(detectCores() - 2)

# county list
counties <- readRDS("data/county_code_list.rds")

# pre-filtered ebd
keep_cols <- readRDS("data/ebd-columns.rds")
f_ebd <- "data/ebd_2020-05.txt"
f_sed <- "data/sed_2020-05.txt"

# process by county
td <- tempdir()
check_status <- function(x) {
  yr <- lubridate::year(x)
  if (any(yr < 2020)) {
    if (any(yr == 2020)) {
      "continuing"
    } else {
      "old"
    }
  } else {
    return("new")
  }
}
tax <- select(ebird_taxonomy, scientific_name, species_code)
county_stats <- foreach (cnty = counties, .combine = bind_rows) %dopar% {
  # extract county data
  f_ebd_tmp <- file.path(td, str_glue("ebd_{cnty}.txt"))
  f_sed_tmp <- file.path(td, str_glue("sed_{cnty}.txt"))
  auk_ebd(f_ebd, f_sed) %>% 
    auk_county(cnty) %>% 
    auk_filter(f_ebd_tmp, f_sed_tmp)
  
  # import and deal with taxonomy rollup
  ebd <- read_ebd(f_ebd_tmp, unique = FALSE)
  sed <- read_sampling(f_sed_tmp, unique = FALSE)
  unlink(f_ebd_tmp)
  unlink(f_sed_tmp)
  
  # identify observer status
  obs_status <- sed %>% 
    group_by(observer_id) %>% 
    summarise(observer_status = check_status(observation_date)) %>% 
    ungroup() %>% 
    mutate(experienced_obs = (observer_status != "new"))
  
  # identify checklist status
  checklist_status <- sed %>% 
    mutate(checklist_id = coalesce(group_identifier, 
                                   sampling_event_identifier)) %>% 
    inner_join(obs_status, by = "observer_id") %>% 
    group_by(checklist_id) %>% 
    summarise(experienced_obs = any(experienced_obs)) %>% 
    ungroup()
  status_counts <- table(obs_status$observer_status)
  
  # remove duplicates
  ebd <- auk_unique(ebd)
  sed <- auk_unique(sed, checklists_only = TRUE) %>% 
    inner_join(checklist_status, by = "checklist_id")
  
  # zero fill
  zf <- auk_zerofill(ebd, sed)
  
  # convert scientific names to codes
  zf$observations <- zf$observations %>% 
    inner_join(tax, by = "scientific_name") %>% 
    select(checklist_id, species_code, species_observed, observation_count)
  
  # prevalence
  prevalence <- checklist_status %>%
    filter(experienced_obs) %>% 
    semi_join(zf$observations, ., by = "checklist_id") %>% 
    group_by(species_code) %>%
    summarise(n_checklists = n(),
              n_pos_obs = sum(species_observed),
              prevalence = mean(species_observed)) %>%
    ungroup() %>% 
    mutate(county_code = cnty) %>% 
    select(county_code, species_code, 
           n_checklists, n_pos_obs, prevalence)
  
  # convert x to -1, drop zeros
  zf$observations <- zf$observations %>% 
    mutate(observation_count = if_else(observation_count == "X", -1L, 
                                       as.integer(observation_count))) %>% 
    filter(species_observed) %>% 
    select(checklist_id, species_code, observation_count)
  
  # nicer order
  cols <- intersect(keep_cols, names(zf$sampling_events))
  zf$sampling_events <- select(zf$sampling_events, 
                               checklist_id, 
                               any_of(cols),
                               experienced_obs)
  
  # output
  str_glue("checklists_{cnty}.rds") %>% 
    file.path("data", "county", .) %>% 
    write_rds(zf$sampling_events, .)
  str_glue("observations_{cnty}.rds") %>% 
    file.path("data", "county", .) %>% 
    write_rds(zf$observations, .)
  str_glue("prevalence_{cnty}.rds") %>% 
    file.path("data", "county", .) %>% 
    write_rds(prevalence, .)
  
  tibble(county = cnty,
         n_checklists = nrow(zf$sampling_events),
         n_included = sum(zf$sampling_events$experienced_obs),
         n_species = n_distinct(zf$observations$species_code),
         n_observers = sum(status_counts),
         n_continuing = coalesce(status_counts["continuing"], 0),
         n_old = coalesce(status_counts["old"], 0),
         n_new = coalesce(status_counts["new"], 0))
}
write_csv(county_stats, "data/county-stats.csv")

# combine counties
# checklists
checklists <- str_glue("checklists_{counties}.rds") %>% 
  file.path("data", "county", .) %>% 
  map_dfr(read_rds)
write_rds(checklists, "data/checklists.rds", compress = "gz")
# observations
observations <- str_glue("observations_{counties}.rds") %>% 
  file.path("data", "county", .) %>% 
  map_dfr(read_rds)
write_rds(observations, "data/observations.rds", compress = "gz")

# prevalence
prevalence <- str_glue("prevalence_{counties}.rds") %>% 
  file.path("data", "county", .) %>% 
  map_dfr(read_rds)
# total
prevalence_total <- prevalence %>% 
  group_by(species_code) %>% 
  summarize(n_checklists = sum(n_checklists),
            n_pos_obs = sum(n_pos_obs)) %>% 
  ungroup() %>% 
  mutate(county_code = "total",
         prevalence = n_pos_obs / n_checklists) %>% 
  select(any_of(names(prevalence)))
bind_rows(prevalence_total, prevalence) %>% 
  write_csv("data/prevalence.csv")

# taxonomy
ebird_taxonomy %>% 
  filter(category == "species") %>% 
  select(species_code, common_name, scientific_name, family) %>% 
  inner_join(prevalence_total %>% select(-county_code), by = "species_code") %>% 
  write_csv("data/ebird-taxonomy.csv")








# Subsetting checklists ---------------------------------------------------

# Spatial functions
library(dggridR) # Creating Hexagonal grid
library(sf) # for gis functions



# Directory for GIS data:
gis_dir <- "C:/Docs/GIS_raw"

# Directory for processed species data:
out_dir <- "C:/Docs/ebird_processed/by_species"



# Loading checklists, with a few additional filters for data quality:
checklists_f <- here("data", "checklists.rds") %>% 
  read_rds() %>%
  # Removing a couple of grouped checklists that evaded the first round
  distinct(checklist_id, .keep_all = T) %>%
  # give zero distance to stationary checklists
  mutate(effort_distance_km = if_else(protocol_type == "Stationary", 0, 
                                      effort_distance_km)) %>%
  # Additional filters
  filter(number_observers <= 10,
         effort_distance_km <= 3,
         experienced_obs == T) %>%
  filter(locality_type == "P" | locality_type == "H") %>%
  # Add year column
  mutate(observation_year = lubridate::year(observation_date)) %>%
  # Pre/post column
  mutate(pandemic = (observation_year == 2020))

# Create hexagonal grid system
# generate hexagonal grid with ~ 3 km betweeen cells
dggs <- dgconstruct(spacing = 3)

# For each checklist, get the gridcell ID
checklists_f <- checklists_f %>%
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)


# Function to subsample. This randomly selects an equal number of checklists from pre- and post-pandemic, up to a maximum threshold in each county (thresh). It selects them in pairs from each gridcell in the county, looping through gridcells until the threshold is reached (or until every gridcell in the county no longer has paired pre- and post- checklists to choose from.

SpatialSubsample <- function(dat, thresh) {
  # Object for full selection:
  select <- NULL
  
  # Run selection algorithm
  for (i in 1:length(unique(dat$county_code))) {
    
    # Selection from county
    cty_select <- NULL
    
    # Current county
    cty <- unique(dat$county_code)[i]
    
    # Grid cells covered
    cells <- dat %>% filter(county_code == cty) %>% distinct(cell)
    
    # Split data into pre and post pandemic
    pre <- dat %>% filter(county_code == cty,
                          pandemic == FALSE)
    post <- dat %>% filter(county_code == cty,
                           pandemic == TRUE)
    
    # Step through cells to sample
    while(length(cty_select) < thresh &
          intersect(
            (pre %>% group_by(cell, .drop = FALSE) %>% count())$cell,
            (post %>% group_by(cell, .drop = FALSE) %>% count())$cell
          ) %>% length() > 0) {
      
      # Loop through each cell
      for (j in 1:nrow(cells)) {
        
        # Current cell
        cl <- cells$cell[j]
        
        # Sample
        if(any(pre$cell == cl) &
           any(post$cell == cl) &
           length(cty_select) < thresh){
          # Sample one from each
          s1 <- sample((pre %>% filter(cell == cl))$checklist_id, size = 1)
          s2 <- sample((post %>% filter(cell == cl))$checklist_id, size = 1)
          
          # Add to county selection
          cty_select <- c(cty_select, s1, s2)
          
          # Remove from "supply"
          pre <- pre[-which(pre$checklist_id == s1),]
          post <- post[-which(post$checklist_id == s2),]
        }
      }
    }
    
    # Add to full selection
    select <- c(select, cty_select)
  }
  
  return(select)
  
}

# Run selection
set.seed(22)

tt <- proc.time() # For run time, later
selection <- SpatialSubsample(dat = checklists_f, thresh = 1000)
run_time <- (proc.time()[3] - tt[3]) / 60 # Run time in minutes
# takes ~ 7 minutes on my laptop

checklists_sample <- checklists_f %>% filter(checklist_id %in% selection)

# Save
saveRDS(checklists_sample, file = "checklists_1000_v2.rds")
saveRDS(selection, file = "selection_1000_v2.rds")








# I use the list of checklists to get the locations associated with those checklists

# Jessica will extract the landcover and distance info, and those columns can then be added to the subsetted checklists object to use in the zero-filling step for each species.


# checklist data, subset columns, get distinct locations
locations_sub <- here("data", "checklists_1000_v2.rds") %>% 
  read_rds() %>%
  select(latitude, longitude, country_code, state_code, county_code, locality_id, locality_type) %>%
  distinct(locality_id, .keep_all = T)

# Number of checklists per location
location_count_sub <- here("data", "checklists_1000_v2.rds") %>% 
  read_rds() %>%
  group_by(locality_id) %>%
  summarize(count = n())

# Merge objects
locations_sub <- merge(locations_sub, location_count_sub, by = "locality_id")

# Save list
write_csv(locations_sub, path = "locations_subset_1000_v2.csv")







# Extract data for focal species ------------------------------------------


# setting the threshold:
# drop observations from counties with lower than this prevalence
prevalence_thresh <- 0.05

# drop counties with less than this number of sampled checklists
check_thresh <- 200




# Load list of species
species_list <- read_csv("species_list_4let.csv")




# Function to Obtain and zero-fill a single species data product:
SpeciesPrep_1000 <- function(species) {
  # define species
  taxonomy <- here("data", "ebird-taxonomy.csv") %>% 
    read_csv() %>% 
    filter(common_name == species)
  
  # input data
  # prevalence
  prevalence <- here("data", "prevalence.csv") %>% 
    read_csv() %>% 
    filter(species_code == taxonomy$species_code,
           county_code != "total",
           prevalence > prevalence_thresh)
  
  # Get list of counties to keep 
  keep_cty <- here("data", "checklists_1000_v2.rds") %>% 
    read_rds() %>%
    count(county_code) %>%
    filter(n >= check_thresh)
  
  # checklist data
  checklists <- here("data", "checklists_1000_v2.rds") %>% 
    read_rds() %>% 
    # drop counties below threshold
    semi_join(prevalence, by = "county_code") %>%
    # drop counties with poor coverage
    filter(county_code %in% keep_cty$county_code) ### debug
  
  # counts
  observations <- here("data", "observations.rds") %>% 
    read_rds() %>% 
    # subset to species
    filter(species_code == taxonomy$species_code) %>% 
    mutate(species_observed = TRUE)
  
  # zero fill
  zf <- left_join(checklists, observations, by = "checklist_id") %>% 
    # deal with X/-1 and fill zeros
    mutate(species_observed = coalesce(species_observed, FALSE),
           observation_count = coalesce(observation_count, 0L),
           observation_count = if_else(observation_count == -1, NA_integer_, 
                                       observation_count),
           species_code = coalesce(species_code, taxonomy$species_code))
  
  return(zf)
}




# Run for all species on the list
for (i in 1:nrow(species_list)) {
  SpeciesPrep_1000(species = species_list$common_name[i]) %>%
    saveRDS(file = file.path(out_dir, paste(species_list$band_code[i],
                                            "_zf.rds", sep = "")))
}





# Model data prep ---------------------------------------------------------


library(tidyverse)
library(lubridate)



# Directories 

# Directory for processed species data:
out_dir <- "C:/Docs/ebird_processed/by_species"

# Directory for mobility data
cov_dir <- "C:/Docs/covariate_data"

# Directory for model output:
mod_dir <- "C:/Docs/model_output"



# Species List 

sp_list <- read_csv(file = file.path(out_dir, "species_list_4let.csv"))



# Load checklist metadata 

checklists <- readRDS(file = file.path(cov_dir, "checklists_1000_v2.rds"))


# Load county metadata 

cty_meta <- read_csv(file = file.path(cov_dir, "county_metadata.csv"))



# Time variables 

# Day of Year
checklists <- checklists %>% mutate(day_of_year = lubridate::yday(observation_date))

# Month
checklists <- checklists %>% mutate(month = lubridate::month(observation_date, label = TRUE))






# Load covariates


# Day-specific Mobility
mobility_day <- read_csv(file = file.path(cov_dir, "mobility_cumulative.csv"))


# Mean mobility for full period
mobility_full_mean <- read_csv(file = file.path(cov_dir, "mobility_cumulative_full_mean.csv"))


# Location-specific data (extracted from GIS)

# Road data came from:
# https://www.arcgis.com/home/item.html?id=06e71cbbefab401fb99b6c2bb5139487

# Airport data came from:
# https://open.canada.ca/data/en/dataset/3a1eb6ef-6054-4f9d-b1f6-c30322cd7abf
# and
# https://ais-faa.opendata.arcgis.com/datasets/e747ab91a11045e8b3f8a3efd093d3b5_0

# And landcover data came from:
# http://www.cec.org/north-american-environmental-atlas/land-cover-30m-2015-landsat-and-rapideye/

loc_cov_data <- read_csv(file = file.path(cov_dir, "location_covariate_data.csv"))






# Mobility/Delta-Traffic 

# Add centered delta-traffic
# Inverse of relative time-at-home
# Negative = county with stronger reduction in traffic,
# Positive = county with weaker reduction in traffic
checklists <- mobility_full_mean %>%
  select(county_code, home) %>%
  mutate(delta_traffic = home *-1) %>%
  mutate(c_delta_traffic = scale(delta_traffic, scale = FALSE)) %>%
  select(county_code, c_delta_traffic) %>%
  left_join(checklists, ., by = "county_code")





# Determine period of 'peak' lockdown in each county:
mobility_peak <- mobility_day %>%
  mutate(day_of_year = lubridate::yday(date)) %>%
  filter(date >= "2020-03-01" & date <= "2020-05-31") %>%
  group_by(county_code) %>%  # do everything below for each county separately
  arrange(day_of_year) %>%  # order by day
  mutate(cumul_lock = cumsum(home_1)) %>%  # calculate accumulation curve
  # Find count cutoff for lower quartile, and then last day counts were below
  mutate(area_25 = max(cumul_lock) * 0.25) %>%
  mutate(lock_peak_start = max(day_of_year[cumul_lock <= area_25])) %>%
  # Find count cutoff for upper quartile, and then last day counts were below
  mutate(area_75 = max(cumul_lock) * 0.75) %>%
  mutate(lock_peak_end = max(day_of_year[cumul_lock <= area_75])) %>%
  # Wrap up
  select(county_code, lock_peak_start, lock_peak_end) %>%
  distinct()




# Distance to Road/Airport 

# Add centered, logged distance to major road and distance airport
# Negative = closer than average to road/airport,
# Positive = further than average to road/airport
checklists <- loc_cov_data %>%
  mutate(air = scale(log(air_dist), scale = FALSE)) %>%
  mutate(road = scale(log(road_maj_dist), scale = FALSE)) %>%
  select(locality_id, air, road) %>%
  left_join(checklists, ., by = "locality_id")







# Landcover 

# Developed/undeveloped, based on 50m buffer of landcover
# Rounded and subtracted from 1 so that:
# 0 = mostly developed (i.e. 'urban')
# 1 = mostly undeveloped (i.e. 'rural')
checklists <- loc_cov_data %>%
  mutate(land = 1-round(urban_buff)) %>%
  select(locality_id, land) %>%
  left_join(checklists, ., by = "locality_id")






# Overlap by species 

# Create the functions that will calculate overlap by species

# Overlap value itself for each row
Overlap <- function(x) {
  lock <- seq(x$lock_peak_start, x$lock_peak_end)
  mig <- seq(x$mig_peak_start, x$mig_peak_end)
  return((sum(mig %in% lock) / length(mig)) - 1)
}

# Function to run on each species
CalcOverlap <- function(species, bird_dat) {
  # Determine period of 'peak' bird presence (i.e. interquartile range of cumulative abundance)
  mig_peak <- bird_dat %>%
    select(county_code, day_of_year, observation_count) %>%
    group_by(county_code) %>%  # do everything below for each county separately
    arrange(day_of_year) %>%  # order by day
    mutate(cumul_count = cumsum(observation_count)) %>%  # calculate accumulation curve
    # Find count cutoff for lower ~25% area, and then last day counts were below
    # (if cutoff is on the first day, it can throw an error, hence need for ifelse statement)
    mutate(area_25 = max(cumul_count) * 0.25) %>%
    mutate(mig_peak_start = ifelse(
      is.finite(max(day_of_year[cumul_count <= area_25])),
      max(day_of_year[cumul_count <= area_25]),
      min(day_of_year))) %>%
    # Find count cutoff for upper ~25%, and then last day counts were below
    mutate(area_75 = max(cumul_count) * 0.75) %>%
    mutate(mig_peak_end = max(day_of_year[cumul_count <= area_75])) %>%
    # Wrap up
    select(county_code, mig_peak_start, mig_peak_end) %>%
    distinct()
  
  # Combine with county metadata
  cty_meta_sp <- left_join(cty_meta, mig_peak, by = "county_code") %>%
    left_join(mobility_peak, by = "county_code") %>%
    mutate(used = !is.na(mig_peak_start))
  
  # Calculate Overlap:
  cty_meta_sp <- cty_meta_sp %>%
    filter(used == TRUE) %>%
    rowwise() %>%
    do(row = as_tibble(.)) %>%
    mutate(overlap = Overlap(row)) %>%
    unnest(cols = c(row)) %>%
    select(county_code, overlap) %>%
    inner_join(cty_meta_sp, by = "county_code")
  
  # Return the overlap by county
  return(cty_meta_sp)
}




# Prep Function 

# Function to prep a single species:
DatPrep <- function(species) {
  # Load species and remove presence-only rows
  zf <- readRDS(file = file.path(out_dir, paste(species, "_zf.rds", sep = ""))) %>%
    filter(!is.na(observation_count))
  # Convert Pandemic variable
  # -1 = Pre(2017-19); 0 = Post(2020)
  zf$pandemic <- (zf$pandemic %>% as.numeric()) - 1
  # Add other data
  zf <- checklists %>%
    select(checklist_id, day_of_year, c_delta_traffic, air, road, land) %>%
    left_join(zf, ., by = "checklist_id")
  # Remove any NA landcover checklists (should only be from one bad location)
  zf <- zf %>% filter(!is.na(land))
  # Add overlap by county
  zf <- CalcOverlap(species = species, bird_dat = zf) %>%
    select(county_code, overlap) %>%
    left_join(zf, ., by = "county_code")
  # Return only the data needed for the model
  dat <- zf %>%
    mutate(sp = species) %>%
    mutate(y = observation_count) %>%
    mutate(cty = county_code) %>%
    mutate(delta_traf = c_delta_traffic) %>%
    mutate(dist = effort_distance_km) %>%
    mutate(dur = duration_minutes) %>%
    select(sp, checklist_id, day_of_year,
           y, cty, overlap, delta_traf,
           pandemic, dist, dur, road, air, land)
  return(dat)
}





# Single Species 
# Running the data prep for one species:

species <- "WODU"

assign(paste("dat_", species, sep = ""),
       DatPrep(species = species))





# All Species 

# List of species:
sp_list$band_code

# Running the data prep for all species, and combining them into one file
# Note: will throw an error, but should still work
dat_full <- tibble()
for (i in 1:nrow(sp_list)) {
  dat_full <- DatPrep(species = sp_list$band_code[i]) %>%
    bind_rows(dat_full)
}

# Just some checks:
is.nan(dat_full$overlap) %>% sum()
is.na(dat_full$overlap) %>% sum()
is.infinite(dat_full$overlap) %>% sum()


# Export as rds
saveRDS(dat_full, file = file.path(out_dir, "dat_full.rds"))

sink(file = file.path(out_dir, "dat_full_readme.txt"))
cat("
This data file combines all of the filtered and zero-filled eBird data for the species we are using.
Each row corresponds to a checklist that was selected according to the filters we used and the spatial subsample (described elsewhere).
There should be between 200-1000 checklists per species, per county.
This data file is designed to use for individual species by county models.
Not all counties are represented in each species: only the counties for which the species occurs in at least 5% of checklists.

Column descriptions:

sp:
Species, identified by 4-letter banding code.

checklist_id:
ID number of checklist

day_of_year:
Day of observation

y:
The observed count of that species.

cty:
The county, identified by eBird county codes.

overlap:
A re-classfified proportion of days during the middle (aka 'peak') of migration that overlap the middle of the lockdown period for that county.
It was re-classified so that zero corresponded to 100% overlap, to make it easier to interpret other main effects in counties where the overlap was the highest.
Note: this measure is specific to the county, and should therefore NOT be included in models specific to each county.
-1 = 0% overlap
0 = 100% overlap

delta_traf:
A centered measure of county-specific change in traffic during the pandemic.
Negative = a county with stronger than average reduction in traffic
Positive = a county with weaker than average reduction in traffic

pandemic:
A binary variable for 'year' category.
-1 = Pre-pandemic (2017-19)
0 = 'Post'-pandemic (2020)

dist:
Distance traveled by the observers, in km.

dur:
Duration of the checklist, in minutes.

road:
Distance to closest major road, logged and then centered (among all checklists, prior to selecting which counties would be included for each species).
Negative = closer than average to road
Positive = further than average to road

air:
Distance to the city's major airport, logged and then centered (among all checklists, prior to selecting which counties would be included for each species).
Negative = closer than average to airport
Positive = further than average to airport

land:
Binary variable representing immediate landcover.
Based on percent of pixels within 50m of location that are developed/un-developed, then rounded.
0 = mostly developed (i.e. 'urban')
1 = mostly undeveloped (i.e. 'rural')
", fill = TRUE)
sink()









# Running Model -----------------------------------------------------------

# Designed to be run on the U of Manitoba compute cluster:



# Set species

species <- "PAWA"


# Start error sink 

errors <- file(paste("errors_", species, ".txt", sep = ""), open="wt")
sink(errors, type="message")

# Libraries 

# Include the library path to my workspace
myPaths <- .libPaths()
myPaths <- c(myPaths, '/local/workspace01/schrimpm/R_libs')
.libPaths(myPaths)

# Load libraries
library(dplyr)


# Stan settings
options(mc.cores = 5)
rstan::rstan_options(auto_write = TRUE)


# Set directories 

# Directory for data
data_dir <- "/home/u1/schrimpm/Documents/Data/"

# Directory for processed species data:
out_dir <- "/home/u1/schrimpm/Documents/Output/"



# Data 

dat <- readRDS(file = file.path(data_dir, "dat_full.rds")) %>%
  filter(sp == species)


# Run Model 

tt <- proc.time() # For run time, later
# Run:
fit <- rstanarm::stan_glmer(y ~ (1|cty) + dist + dur +
                              pandemic + overlap + overlap:pandemic +
                              delta_traf + delta_traf:pandemic + 
                              air + air:pandemic +
                              road + road:pandemic +
                              land + land:pandemic,
                            data = dat,
                            family = rstanarm::neg_binomial_2,
                            chains = 5,
                            warmup = 300,
                            iter = 1300,
                            thin = 1,
                            cores = 5,
                            seed = 121212,
                            refresh = 1)
run_time <- (proc.time()[3] - tt[3]) / 60 # Run time in minutes

# Save Output:
saveRDS(fit, file = file.path(out_dir, paste("fit_", species, ".rds", sep = "")))
saveRDS(run_time, file = file.path(out_dir, paste("time_", species, ".rds", sep = "")))


# Finish error sink 

sink(type="message")
close(errors)









# Interpreting results ----------------------------------------------------


# Libraries 

library(tidyverse)
library(MCMCvis)
library(bayesplot)


# Set directories 

# Directory for bird data
out_dir <- "C:/Docs/ebird_processed/by_species"

# Directory for GIS and mobility data
cov_dir <- "C:/Docs/covariate_data"

# Directory for model output:
mod_dir <- "C:/Docs/model_output"




# Species List 

sp_list <- read_csv(file = file.path(out_dir, "species_list_4let.csv"))

# Remove dropped species
'%notin%' <- Negate('%in%')
sp_list <- sp_list %>% filter(band_code %notin% c("NSHO","RNDU"))


sp_cats <- read_csv(file = file.path(out_dir, "species_list_categories.csv"))
sp_cats <- sp_cats %>% filter(band_code %notin% c("NSHO","RNDU"))



# Parameter posteriors to summarize
pars <- c("dist",
          "dur",
          "pandemic",
          "overlap",
          "delta_traf",
          "air",
          "road",
          "land",
          "pandemic:overlap",
          "pandemic:delta_traf",
          "pandemic:air",
          "pandemic:road",
          "pandemic:land",
          "reciprocal_dispersion")

# Summary stats desired for each posterior
coln <- c("mean",
          "sd",
          "lower",
          "median",
          "upper",
          "Rhat",
          "ESS")


# Function to extract data for each species, rounded to a certain number of significant figures
PostSum <- function(species_list, figs) {
  # Create object to store results:
  res_name_list <- vector(length = length(coln)*length(pars) + 3)
  res_name_list[1:3] <- c("Species_common", "Species_sci", "Band_code")
  for(i in 1:length(pars)) {
    for (j in 1:length(coln)) {
      res_name_list[3 + length(coln) * (i - 1) + j] <- paste0(pars[i],"-",coln[j])
    }
  }
  res <- matrix(nrow = nrow(sp_list), ncol = length(res_name_list))
  dimnames(res) <- list(
    NULL,
    res_name_list
  )
  res[,1:3] <- select(sp_list, common_name, scientific_name, band_code) %>% as.matrix()
  
  # Fill each row with appropriate species results
  for (i in 1:nrow(res)) {
    species <- res[i,"Band_code"]
    
    # Extract species results
    fit <- readRDS(file = file.path(mod_dir, paste("fit_", species, ".rds", sep = "")))
    
    # Fill row with summary stats, rounded
    res[i,-(1:3)] <- MCMCsummary(fit, params = pars, HPD = FALSE) %>%
      as.matrix() %>% t() %>% as.vector() %>% signif(digits = figs)
  }
  
  return(res)
  
}

post_summary <- PostSum(species_list = sp_list, figs = 4) %>% as.tibble
post_summary <- post_summary %>% as_tibble

# Save as csv
write_csv(post_summary, path = file.path(mod_dir, "post_summary.csv"))

# reload
post_summary <- read_csv(file = file.path(mod_dir, "post_summary.csv"))


sink(file = file.path(mod_dir, "post_summary_readme.txt"))
cat("
This file provides summary statistics for parameter posterior distrbutions from each of the single species models.
Each column represents a different summary statistic from a different parameter, and the columns are labeled with the parameter name followed by the statistic, separated by a hyphen (-).

The first three columns provide the species common name, scientific name, and 4-letter banding code.

All except the last parameter included (each of which is described in more detail in the methods) are the coefficients for each of the following variables:
dist = distance traveled
dur = duration
pandemic = main effect of pandemic (pre- vs post-pandemic)
overlap = main effect of overlap between peak traffic reduction and peak bird counts
delta_traf = main effect of traffic reduction
air = main effect of logged distance from airport
road = main effect of logged distance from road
land = main effect of local landcover (urban vs rural)
pandemic:overlap = interaction between pandemic and overlap
pandemic:delta_traf = interaction between pandemic and delta_traf
pandemic:air = interaction between pandemic and air
pandemic:road = interaction between pandemic and road
pandemic:land = interaction between pandemic and land

The last parameter is:
reciprocal_dispersion = the measure of overdispersion used by rstanarm to model the negative binomial distribution



Each parameter has seven summary statistics, which are:
mean = mean
sd = standard deviation
lower = 2.5% quantile, i.e., lower limit of 95% credible interval
median = median, i.e. 50% quantile
upper = 97.5% quantile, i.e., upper limit of 95% credible interval
Rhat = R-hat value, used to check convergence
ESS = crude measure of effective sample size of the posterior
", fill = TRUE)
sink()





# Function to plot all interaction CIs and medians in balck and white:
PostSumBW <- function(trm) {
  # Extract data
  dat <- post_summary %>%
    select(Band_code,
           paste("pandemic:", trm, "-median", sep = ""),
           paste("pandemic:", trm, "-upper", sep = ""),
           paste("pandemic:", trm, "-lower", sep = "")) %>%
    rename(band_code = Band_code,
           median = paste("pandemic:", trm, "-median", sep = ""),
           upper = paste("pandemic:", trm, "-upper", sep = ""),
           lower = paste("pandemic:", trm, "-lower", sep = "")) %>%
    left_join((sp_cats %>% select(band_code,
                                  common_name,
                                  cat_1)), ., by = "band_code") %>%
    arrange(median)
  
  # Conversion to the correct sign (due to reversing axis direction for all variables other than 'overlap')
  if (trm == "overlap") {
    dat <- dat
  } else {
    dat <- dat %>%
      mutate(median_old = median,
             lower_old = lower,
             upper_old = upper,
             median = median_old * -1,
             lower = upper_old * -1,
             upper = lower_old * -1) %>%
      arrange(median)
  }
  
  # Plotting info:
  y <- -(1:nrow(dat))
  xlims <- c(min(dat$lower), max(dat$upper))
  
  # Margin parameters
  par(mar = c(0.1, 6.1, 3.1, 0.1),
      mgp = c(2, 0.5, 0))
  
  # Blank plot
  plot(0,
       xaxt = 'n',
       yaxt = 'n',
       bty = 'n',
       pch = '',
       ylab = '', xlab = '',
       xlim = xlims, ylim = c(-1,-length(y))
  )
  
  # X axis
  axis(side = 3, tck = -0.01, cex.axis = 0.75)
  mtext(text = paste("pandemic:", trm, " interaction", sep = ""), side = 3, line = 2)
  
  # Species Guidelines
  for (i in 1:nrow(dat)) {
    lines(x = xlims,
          y = c(y[i], y[i]),
          lty = 3,
          col = rgb(0.8,0.8,0.8,1))
  }
  
  # Zero Line
  abline(v = 0, lty = 2)
  
  # Data
  for (i in 1:nrow(dat)) {
    points(x = dat$median[i], y = y[i],
           pch = 16, cex = 0.75,
           col = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                    (dat$upper[i] < 0 & dat$lower[i] < 0)) {
             rgb(0,0,0,1)} else {
               rgb(0.6,0.6,0.6,1)
             })
  }
  
  
  for (i in 1:nrow(dat)) {
    lines(x = c(dat$lower[i], dat$upper[i]),
          y = c(y[i], y[i]),
          lwd = 1.5,
          col = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                   (dat$upper[i] < 0 & dat$lower[i] < 0)) {
            rgb(0,0,0,1)} else {
              rgb(0.6,0.6,0.6,1)
            })
  }
  
  # Species names
  par(mgp = c(0, 0, 0))
  for (i in 1:nrow(dat)) {
    axis(labels = dat$common_name[i],
         at = y[i],
         side = 2, las = 2, tick = FALSE,
         cex.axis = 0.3,
         font = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                   (dat$upper[i] < 0 & dat$lower[i] < 0)) {
           2} else {1}
    )
  }
  
}



PostSumBW(trm = "air")






# Marginal effects/Interaction Plots

### Function to make the species plots (all panels)
SpeciesIntPlot <- function(species) {
  # Load results
  fit <- readRDS(file = file.path(mod_dir, paste("fit_", species, ".rds", sep = "")))
  
  # Panel arrangement
  par(mfrow = c(3,2),
      mar = c(2, 2, 0.3, 0.3),
      mgp = c(1.2, 0.3, 0),
      las = 1)
  
  # Plot name and legend
  plot(0,
       xaxt = 'n', yaxt = 'n', bty = 'n',
       pch = '',
       ylab = '', xlab = '',
       xlim = c(-1,1), ylim = c(-1.5,1.5))
  
  
  text(x = 0, y = 0.5,
       labels = sp_list %>% filter(band_code == species) %>%
         select(common_name) %>% paste(),
       cex = 1.5
  )
  text(x = 0, y = 0, font = 3,
       labels = sp_list %>% filter(band_code == species) %>%
         select(scientific_name) %>% paste(),
       cex = 1)
  legend(x = 0, y = -1, ncol = 2, cex = 0.75,
         legend = c("2017-19", "2020"),
         col = c(rgb(0.9,0.5,0,1),
                 rgb(0.5,0,0.8,1)),
         lwd = 2,
         xjust = 0.5, yjust = 0.5)
  
  
  # Custom functions for floor and ceiling that come from:
  # https://stackoverflow.com/questions/35807523/r-decimal-ceiling)
  # They will calc the floor and ceiling to 0.1:
  floor_dec <- function(x, level = 1) round(x - 5 * 10^(-level - 1), level)
  ceiling_dec <- function(x, level = 1) round(x + 5  *10^(-level - 1), level)
  
  # Function to predict and plot each variable
  IntPlot <- function(trm) {
    
    # Decide whether to show the plot:
    if ( # if the posterior does not occur all above or all below zero (i.e. overlaps zero)
      !(((post_summary %>% filter(Band_code == species) %>%
          select(paste("pandemic:", trm, "-lower", sep = "")) %>%
          pull > 0) &
         (post_summary %>% filter(Band_code == species) %>%
          select(paste("pandemic:", trm, "-upper", sep = "")) %>%
          pull > 0)) |
        ((post_summary %>% filter(Band_code == species) %>%
          select(paste("pandemic:", trm, "-lower", sep = "")) %>%
          pull < 0) &
         (post_summary %>% filter(Band_code == species) %>%
          select(paste("pandemic:", trm, "-upper", sep = "")) %>%
          pull < 0)))
    ) {
      # Plot for 'non-significant' variable
      plot(0,
           xaxt = 'n', yaxt = 'n', bty = 'n',
           pch = '',
           ylab = '', xlab = '',
           xlim = c(-1,1), ylim = c(-1,1))
      
      text(x = 0, y = 0,
           labels = paste("posterior distribution 95% CI \n for interaction with\n*",
                          trm, "* includes zero", sep = ""),
           cex = 1)
      
    } else { # if the posterior did not overlap zero
      
      # Set limits and resolution for ggpredict (auto limits are not great)
      if (trm == "land") {
        mn <- 0
        mx <- 1
      } else {
        mn <- fit$x[,trm] %>% min() %>% floor_dec(level = 1)
        mx <- fit$x[,trm] %>% max() %>% ceiling_dec(level = 1)
      }
      tms <- c(paste(trm, " [", mn, ":", mx, ", by = 0.1]", sep = ""), "pandemic")
      
      
      # Run predict function:
      temp <- ggeffects::ggpredict(fit, tms)
      
      
      
      ### Values to plot
      # Median bird values
      y_plot <- temp %>% filter(group == 0) %>% select(predicted) %>%
        as.data.frame() %>% rename(post = predicted)
      y_plot <- temp %>% filter(group == -1) %>% select(predicted) %>%
        as.data.frame() %>% rename(pre = predicted) %>% add_column(y_plot)
      
      # CI (upper)
      y_up_plot <- temp %>% filter(group == 0) %>% select(conf.high) %>%
        as.data.frame() %>% rename(post = conf.high)
      y_up_plot <- temp %>% filter(group == -1) %>% select(conf.high) %>%
        as.data.frame() %>% rename(pre = conf.high) %>% add_column(y_up_plot)
      
      # CI (lower)
      y_low_plot <- temp %>% filter(group == 0) %>% select(conf.low) %>%
        as.data.frame() %>% rename(post = conf.low)
      y_low_plot <- temp %>% filter(group == -1) %>% select(conf.low) %>%
        as.data.frame() %>% rename(pre = conf.low) %>% add_column(y_low_plot)
      
      # Axes limits and values
      ylims <- c(min(y_low_plot) - min(y_low_plot) * 0.1,
                 max(y_up_plot) + max(y_up_plot) * 0.1)
      
      # Predictor Values
      x <- temp %>% filter(group == 0) %>% select(x) %>%
        as.data.frame() %>% rename(post = x)
      x <- temp %>% filter(group == -1) %>% select(x) %>%
        as.data.frame() %>% rename(pre = x) %>% add_column(x)
      
      # X values for plotting
      if (trm == "overlap") {
        x_plot <- x
      } else {
        x_plot <- -x
      }
      xlims <- c(
        min(x_plot) %>% floor,
        max(x_plot) %>% ceiling
      )
      
      # Other X axis inputs:
      if (trm == "overlap") {
        x_lab_pos <- c(-1, -0.5, 0)
        x_lab_val <- c(0, 50, 100)
        x_lab_txt <- "Percent Overlap"
      } else {
        if (trm == "delta_traf") {
          
          # Center of delta_traf is: -14.05706
          
          x_lab_val <- seq(from = -25, to = -10, by = 5)
          x_lab_pos <- (x_lab_val - -14.05706) * -1
          
          # x_lab_pos <- seq(from = xlims[1], to = xlims[2],
          #                  by = 3)
          # x_lab_val <- x_lab_pos * -1
          x_lab_txt <- "Percent Change in Traffic"
        } else {
          if (trm == "land") {
            x_lab_pos <- c(-1, 0)
            x_lab_val <- c("Rural", "Urban")
            x_lab_txt <- "Landcover"
          } else {
            if (trm == "air") {
              
              # center of air is: exp(10.26713)
              
              # # To go from the x_plot limits to real distances (km):
              # (exp(-xlims[1] + 10.26713) / 1000) # max distance to airport (km) (i.e. min x_plot)
              # (exp(-xlims[2] + 10.26713) / 1000) # min distance to airport (max xplot)
              
              x_lab_val <- c(1, 10, 100, 500)
              x_lab_pos <- (c(1, 10, 100, 500) * 1000) %>% log() %>% -10.26713 * -1
              x_lab_txt <- "Distance to Airport (km)"
            } else {
              # trm == road
              # center of road is : exp(5.491703)
              
              # # To go from the x_plot limits to real distances (km):
              # (exp(-xlims[1] + 5.491703) / 1000) # max distance to road (km) (i.e. min x_plot)
              # (exp(-xlims[2] + 5.491703) / 1000) # min distance to road (max xplot)
              
              x_lab_val <- c(0.001, 0.01, 0.1, 1, 10, 100)
              x_lab_pos <- (c(0.001, 0.01, 0.1, 1, 10, 100) * 1000) %>% log() %>% -5.491703 * -1
              x_lab_txt <- "Distance to Road (km)"
              
            }
          }
        }
      }
      
      
      ### Draw panel
      # Draw blank plot with specific axes
      plot(0,
           xaxt = 'n',
           pch = '',
           ylab = '', xlab = '',
           xlim = xlims, ylim = ylims,
           cex.axis = 0.65, tck = -0.009)
      title(ylab = "Expected Count", line = 1.3, cex.lab = 0.65)
      title(xlab = x_lab_txt, line = 0.8, cex.lab = 0.65)
      axis(side = 1,
           at = x_lab_pos,
           labels = x_lab_val,
           cex.axis = 0.65, tck = -0.01, mgp = c(0.9, 0, 0)
      )
      # Add pre-CI
      polygon(x = c(x_plot$pre, rev(x_plot$pre)),
              y = c(y_low_plot$pre, rev(y_up_plot$pre)),
              border = NA, col = rgb(0.9,0.5,0,0.4))
      # Add pre-line
      points(x = x_plot$pre,
             y = y_plot$pre,
             type = "l",
             col = rgb(0.9,0.5,0,1),
             lwd = 2)
      
      # Add post-CI
      polygon(x = c(x_plot$post, rev(x_plot$post)),
              y = c(y_low_plot$post, rev(y_up_plot$post)),
              border = NA, col = rgb(0.5,0,0.8,0.25))
      # Add post-line
      points(x = x_plot$post,
             y = y_plot$post,
             type = "l",
             col = rgb(0.5,0,0.8,1),
             lwd = 2)
      
    }
  }
  
  
  # Loop through all five variables
  
  vars <- c("overlap",
            "delta_traf",
            "land",
            "air",
            "road")
  
  for (i in 1:length(vars)) {
    IntPlot(trm = vars[i])
  }
  
  
}



SpeciesIntPlot(species = "BAEA")










# Function that calculates marginal effects for one of the variables (and the pandemic variable) for all species with "significant" posteriors for that interaction.
PropPlotMultiple <- function(trm, highlight = NULL) {
  
  ### Getting data for plot
  
  # Select species lists for the increase panel and the decrease panel for that variable
  if (trm == "overlap") {
    spl_inc <- post_summary$Band_code[
      which(post_summary[,paste("pandemic:", trm, "-lower", sep = "")] > 0 &
              post_summary[,paste("pandemic:", trm, "-upper", sep = "")] > 0)
    ]
    spl_dec <- post_summary$Band_code[
      which(post_summary[,paste("pandemic:", trm, "-lower", sep = "")] < 0 &
              post_summary[,paste("pandemic:", trm, "-upper", sep = "")] < 0)
    ]
  } else {
    spl_inc <- post_summary$Band_code[
      which(post_summary[,paste("pandemic:", trm, "-lower", sep = "")] < 0 &
              post_summary[,paste("pandemic:", trm, "-upper", sep = "")] < 0)
    ]
    spl_dec <- post_summary$Band_code[
      which(post_summary[,paste("pandemic:", trm, "-lower", sep = "")] > 0 &
              post_summary[,paste("pandemic:", trm, "-upper", sep = "")] > 0)
    ]
  }
  
  
  # Custom functions for floor and ceiling that come from:
  # https://stackoverflow.com/questions/35807523/r-decimal-ceiling)
  # They will calc the floor and ceiling to 0.1:
  floor_dec <- function(x, level = 1) round(x - 5 * 10^(-level - 1), level)
  ceiling_dec <- function(x, level = 1) round(x + 5  *10^(-level - 1), level)
  
  
  # Function to extract proportional change in abundance by species for the current variable (trm)
  ExtrProp <- function(species) {
    fit <- readRDS(file = file.path(mod_dir, paste("fit_", species, ".rds", sep = "")))
    
    # Set limits and resolution for ggpredict (auto limits are not great)
    if (trm == "land") {
      mn <- 0
      mx <- 1
    } else {
      mn <- fit$x[,trm] %>% min() %>% floor_dec(level = 1)
      mx <- fit$x[,trm] %>% max() %>% ceiling_dec(level = 1)
    }
    
    tms <- c(paste(trm, " [", mn, ":", mx, ", by = 0.1]", sep = ""), "pandemic")
    
    # Run predict function:
    temp <- ggeffects::ggpredict(fit, tms)
    
    
    # Percent difference in predicted abundance
    y <- ( temp[which(temp$group == 0),"predicted"] -
             temp[which(temp$group == -1),"predicted"] ) /
      temp[which(temp$group == -1),"predicted"] * 100
    
    # actual X values
    x <- temp$x[which(temp$group == 0)]
    
    # X values for plotting
    if (trm == "overlap") {
      x_plot <- x
    } else {
      x_plot <- -x
    }
    
    temp_plot <- cbind(x, x_plot, y)
    
    return(temp_plot)
  }
  
  # Create empty list (for species with increases)
  temp_list_inc <- list()
  # Extract values for each species with increases
  for (i in 1:length(spl_inc)) {
    temp_list_inc[[spl_inc[i]]] <- ExtrProp(spl_inc[i])
  }
  
  
  # Create empty list (for species with decreases)
  temp_list_dec <- list()
  # Extract values for each species with increases
  for (i in 1:length(spl_dec)) {
    temp_list_dec[[spl_dec[i]]] <- ExtrProp(spl_dec[i])
  }
  
  # Create a total list (just to get axes limits)
  temp_list_tot <- c(temp_list_inc, temp_list_dec)
  
  
  ### Getting plots set up
  
  # Axes limits
  ylims <- c(
    lapply(temp_list_tot, function(x) { min(x[,"y"]) } ) %>% as.data.frame %>% min %>% floor_dec,
    lapply(temp_list_tot, function(x) { max(x[,"y"]) } ) %>% as.data.frame %>% max %>% ceiling_dec
  )
  xlims <- c(
    lapply(temp_list_tot, function(x) { min(x[,"x_plot"]) } ) %>% as.data.frame %>% min %>% floor,
    lapply(temp_list_tot, function(x) { max(x[,"x_plot"]) } ) %>% as.data.frame %>% max %>% ceiling
  )
  
  # Setting up axes
  if (trm == "overlap") {
    x_lab_pos <- c(-1, -0.5, 0)
    x_lab_val <- c(0, 50, 100)
    x_lab_txt <- "Percent Overlap"
  } else {
    if (trm == "delta_traf") {
      
      # Center of delta_traf is: -14.05706
      
      x_lab_val <- seq(from = -25, to = -10, by = 5)
      x_lab_pos <- (x_lab_val - -14.05706) * -1
      x_lab_txt <- "Percent Change in Traffic"
    } else {
      if (trm == "land") {
        x_lab_pos <- c(-1, 0)
        x_lab_val <- c("Rural", "Urban")
        x_lab_txt <- "Landcover"
      } else {
        if (trm == "air") {
          
          # center of air is: exp(10.26713)
          
          # # To go from the x_plot limits to real distances (km):
          # (exp(-xlims[1] + 10.26713) / 1000) # max distance to airport (km) (i.e. min x_plot)
          # (exp(-xlims[2] + 10.26713) / 1000) # min distance to airport (max xplot)
          
          x_lab_val <- c(1, 10, 100, 500)
          x_lab_pos <- (c(1, 10, 100, 500) * 1000) %>% log() %>% -10.26713 * -1
          x_lab_txt <- "Distance to Airport (km)"
        } else {
          # trm == road
          # center of road is : exp(5.491703)
          
          # # To go from the x_plot limits to real distances (km):
          # (exp(-xlims[1] + 5.491703) / 1000) # max distance to road (km) (i.e. min x_plot)
          # (exp(-xlims[2] + 5.491703) / 1000) # min distance to road (max xplot)
          
          x_lab_val <- c(0.001, 0.01, 0.1, 1, 10, 100)
          x_lab_pos <- (c(0.001, 0.01, 0.1, 1, 10, 100) * 1000) %>% log() %>% -5.491703 * -1
          x_lab_txt <- "Distance to Road (km)"
          
        }
      }
    }
  }
  
  
  ### Plot
  # Left panel (increases)
  par(mar = c(1.6, 1.8, 0.1, 0.1),
      mgp = c(1, 0.3, 0),
      las = 1,
      mfrow = c(1,2))
  # Draw blank plot with specific axes
  plot(0,
       xaxt = 'n',
       pch = '',
       ylab = '', xlab = '',
       xlim = xlims, ylim = ylims,
       cex.axis = 0.5, tck = -0.009)
  title(ylab = "Percent Change", line = 1, cex.lab = 0.75)
  title(xlab = x_lab_txt, line = 0.6, cex.lab = 0.75)
  axis(side = 1,
       at = x_lab_pos,
       labels = x_lab_val,
       cex.axis = 0.5, tck = -0.01, mgp = c(0.9, 0, 0)
  )
  # Add lines
  for (i in 1:length(spl_inc)) {
    points(x = temp_list_inc[[spl_inc[i]]][,"x_plot"],
           y = temp_list_inc[[spl_inc[i]]][,"y"],
           type = "l",
           col = rgb(0.6,0.6,0.6,1),
           lwd = 1.25
    )
    
    
    
  }
  abline(h = 0, lty = 2) # Zero line (Percent change)
  
  # Add highlighted species
  if(!is.null(highlight)) {
    if(highlight %in% spl_inc) {
      points(x = temp_list_inc[[highlight]][,"x_plot"],
             y = temp_list_inc[[highlight]][,"y"],
             type = "l",
             col = rgb(0,0,0,1),
             lwd = 2
      )
    }
  }
  
  # Right panel (decreases)
  par(mar = c(1.6, 0.1, 0.1, 1.8))
  # Draw blank plot with specific axes
  plot(0,
       xaxt = 'n', yaxt = 'n',
       pch = '',
       ylab = '', xlab = '',
       xlim = xlims, ylim = ylims,
       cex.axis = 0.5, tck = -0.009)
  title(xlab = x_lab_txt, line = 0.6, cex.lab = 0.75)
  axis(side = 1,
       at = x_lab_pos,
       labels = x_lab_val,
       cex.axis = 0.5, tck = -0.01, mgp = c(0.9, 0, 0)
  )
  # Add lines
  for (i in 1:length(spl_dec)) {
    points(x = temp_list_dec[[spl_dec[i]]][,"x_plot"],
           y = temp_list_dec[[spl_dec[i]]][,"y"],
           type = "l",
           col = rgb(0.6,0.6,0.6,1),
           lwd = 1.25
    )
  }
  abline(h = 0, lty = 2) # Zero line (Percent change)
  
  if(!is.null(highlight)) {
    if(highlight %in% spl_dec) {
      points(x = temp_list_dec[[highlight]][,"x_plot"],
             y = temp_list_dec[[highlight]][,"y"],
             type = "l",
             col = rgb(0,0,0,1),
             lwd = 2
      )
    }
  }
  
}




