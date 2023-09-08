# # for showing dimensions of data in slide
# load("00_data/ebd_IN_relMay-2022.RData")
# our_small_dataset <- as_tibble(data)
# 
# rm(our_small_dataset, data)

#

# group birding ####


# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0 <- data0_MY_slice_S %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

count <- 0
for (i in 1:4) { # for each state
  
  data_a <- data0_MY_slice_S %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(GROUP.BIRDING = if_else(NUMBER.OBSERVERS > 2, 1, 0)) %>%
    ungroup() 
  
  tictoc::tic(glue::glue("GLMM for {anal_states[, i]}"))
  model_a <- glmer(GROUP.BIRDING ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_a, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  for (j in 1:n_distinct(data0$MONTH)) {
    for (k in 1:n_distinct(data0$M.YEAR)) {
      count <- count + 1
      
      data0$PRED[count] = predict(model_a,
                                  data.frame(MONTH = data0$MONTH[count],
                                             M.YEAR = data0$M.YEAR[count]),
                                  re.form = NA, type = "response")
    }
  }
} 

data0 <- data0 %>% 
  left_join(timeline)

group_birding <- data0



# site fidelity ####


# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0 <- data0_MY_slice_S %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

count <- 0
for (i in 1:4) { # for each state
  
  data_a <- data0_MY_slice_S %>% 
    left_join(obs_mainstates, by = "OBSERVER.ID") %>% 
    filter(MAIN.STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID) %>% 
    dplyr::summarise(NO.SITES = n_distinct(CELL.ID)) %>% 
    mutate(NO.SITES = replace_na(NO.SITES, 0)) %>% 
    ungroup()
  
  tictoc::tic(glue::glue("GLMM for {anal_states[, i]}"))
  model_a <- glmer(NO.SITES ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_a, family = poisson,
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  for (j in 1:n_distinct(data0$MONTH)) {
    for (k in 1:n_distinct(data0$M.YEAR)) {
      count <- count + 1
      
      data0$PRED[count] = predict(model_a,
                                  data.frame(MONTH = data0$MONTH[count],
                                             M.YEAR = data0$M.YEAR[count]),
                                  re.form = NA, type = "response")
    }
  }
} 

data0 <- data0 %>% 
  left_join(timeline)

site_fidelity <- data0


# urban bias ####


data_urban <- data0_MY_slice_G %>% 
  group_by(M.YEAR, MONTH, COUNTY, CELL.ID) %>% 
  filter(URBAN == 1) %>% 
  dplyr::summarise(U.LISTS = n_distinct(GROUP.ID)) %>% 
  right_join(data0_MY_slice_G) %>% 
  group_by(M.YEAR, MONTH, COUNTY, CELL.ID) %>% 
  dplyr::summarise(TOT.LISTS = n_distinct(GROUP.ID),
                   U.LISTS = replace_na(min(U.LISTS), 0),
                   PROP.U = round(U.LISTS/TOT.LISTS, 4)) %>% 
  ungroup()

counties <- data_urban %>% distinct(COUNTY, CELL.ID)

# no MONTH main effect unlike in birds model
tictoc::tic("UNU GLMM using weighted proportions")
model_urban <- glmer(PROP.U ~ M.YEAR + M.YEAR:MONTH + (1|COUNTY/CELL.ID), 
                     weights = TOT.LISTS,
                     data = data_urban, family = binomial(link = "cloglog"),
                     nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) # 119 sec
tictoc::toc()

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
urban1 <- data_urban %>% 
  tidyr::expand(nesting(MONTH, M.YEAR, COUNTY, CELL.ID)) %>% 
  filter(!is.na(COUNTY)) %>% 
  mutate(PROP.U.PRED = predict(model_urban, data.frame(MONTH, M.YEAR, COUNTY, CELL.ID),
                               type = "response")) %>% 
  group_by(M.YEAR, MONTH) %>% 
  summarise(SE = sd(PROP.U.PRED)/sqrt(n()),
            PROP.U.PRED = mean(PROP.U.PRED),
            CI.L = PROP.U.PRED - 1.96*SE,
            CI.U = PROP.U.PRED + 1.96*SE)

# overall coverage and intensity ####


# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0 <- data0_MY_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

count <- 0
for (i in 1:4) { # for each state
  
  data_a <- data0_MY_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    select(M.YEAR, MONTH, CELL.ID, GROUP.ID) %>% 
    tidyr::complete(M.YEAR, MONTH, CELL.ID) %>% 
    group_by(M.YEAR, MONTH, CELL.ID) %>% 
    dplyr::summarise(COVERED = if_else(n_distinct(GROUP.ID) >= 10, 1, 0))
  
  tictoc::tic(glue::glue("GLMM for {anal_states[, i]}"))
  model_a <- glmer(COVERED ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
                   data = data_a, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))  
  tictoc::toc()
  
  for (j in 1:n_distinct(data0$MONTH)) {
    for (k in 1:n_distinct(data0$M.YEAR)) {
      count <- count + 1
      
      data0$PRED[count] = predict(model_a,
                                  data.frame(MONTH = data0$MONTH[count],
                                             M.YEAR = data0$M.YEAR[count]),
                                  re.form = NA, type = "response")
    }
  }
} 

data0 <- data0 %>% 
  left_join(timeline)

coverage <- data0


# bird repfreq ####


require(lme4)


### modelling directly with presence-absence data instead of relative abundance

# to join for presence-absence of various species
temp1 <- data0_MY %>% 
  filter(STATE == "Karnataka") %>% 
  group_by(GROUP.ID, COMMON.NAME) %>% 
  summarise(OBSERVATION.COUNT = max(OBSERVATION.COUNT)) %>% ungroup()

# to later join checklist metadata
temp2 <- data0_MY %>% 
  filter(STATE == "Karnataka") %>% 
  arrange(SAMPLING.EVENT.IDENTIFIER) %>% 
  group_by(GROUP.ID) %>% 
  slice(1) %>% ungroup() %>% 
  select(GROUP.ID, STATE, COUNTY, LOCALITY, LATITUDE, LONGITUDE, OBSERVATION.DATE, 
         M.YEAR, MONTH, DAY.M, M.YEAR, URBAN, CELL.ID, SUBCELL.ID, NO.SP)

data_b2 <- data0_MY_slice_G %>% 
  filter(STATE == "Karnataka") %>% 
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
  filter(SP.CATEGORY != "W") %>% # removing wetland spp. cos we want to compare UNU
  ungroup()


# getting median list length for prediction later
median_length <- data_b2 %>% 
  distinct(M.YEAR, MONTH, GROUP.ID, NO.SP) %>% 
  group_by(MONTH) %>% 
  summarise(NO.SP.MED = median(NO.SP))

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
bang0 <- data_b2 %>% 
  tidyr::expand(COMMON.NAME, MONTH, M.YEAR) %>% 
  left_join(species) %>% 
  mutate(REP.FREQ.PRED = NA)

count <- 0
for (i in 1:n_distinct(bang0$COMMON.NAME)) {
  
  data_spec <- data_b2 %>% 
    filter(COMMON.NAME == unique(bang0$COMMON.NAME)[i]) %>% 
    # filtering for CELL.ID-MONTH (space-time) combos in which species occurs
    filter(REPORT == 1) %>% 
    distinct(COMMON.NAME, CELL.ID, MONTH) %>% 
    left_join(data_b2)
  
  
  model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:NO.SP + MONTH:M.YEAR + 
                        (1|CELL.ID),
                      data = data_spec, family = binomial(link = "cloglog"),
                      nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
  
  for (j in 1:n_distinct(bang0$MONTH)) {
    for (k in 1:n_distinct(bang0$M.YEAR)) {
      count <- count + 1
      
      bang0$REP.FREQ.PRED[count] = predict(
        model_spec,
        data.frame(COMMON.NAME = bang0$COMMON.NAME[count],
                   MONTH = bang0$MONTH[count],
                   M.YEAR = bang0$M.YEAR[count],
                   NO.SP = median_length$NO.SP.MED[bang0$MONTH[count]]),
        re.form = NA, 
        type = "response")
    }
  }
} # 7 mins


bang1 <- bang0 %>% 
  filter(MONTH %in% 4:5) %>% 
  group_by(M.YEAR, SP.CATEGORY) %>% 
  summarise(SE = sd(REP.FREQ.PRED)/sqrt(n()),
            REP.FREQ.PRED = mean(REP.FREQ.PRED),
            CI.L = REP.FREQ.PRED - 1.96*SE,
            CI.U = REP.FREQ.PRED + 1.96*SE)

bang2 <- bang0 %>% 
  group_by(M.YEAR, SP.CATEGORY) %>% 
  summarise(SE = sd(REP.FREQ.PRED)/sqrt(n()),
            REP.FREQ.PRED = mean(REP.FREQ.PRED),
            CI.L = REP.FREQ.PRED - 1.96*SE,
            CI.U = REP.FREQ.PRED + 1.96*SE)


# SUPPL. ####
# birding time ####


# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0 <- data0_MY_slice_S %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

count <- 0
for (i in 1:4) { # for each state
  
  data_a <- data0_MY_slice_S %>% 
    left_join(obs_mainstates, by = "OBSERVER.ID") %>% 
    filter(MAIN.STATE == anal_states[, i]) %>% 
    group_by(MONTH, M.YEAR, OBSERVER.ID) %>% 
    dplyr::summarise(B.TIME = sum(DURATION.MINUTES)) %>% 
    mutate(B.TIME = replace_na(B.TIME, 0)) %>% 
    ungroup()
  
  tictoc::tic(glue::glue("GLMM for {anal_states[, i]}"))
  model_a <- glmer(B.TIME ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_a, family = poisson,
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  for (j in 1:12) {
    for (k in 1:4) {
      count <- count + 1
      
      data0$PRED[count] = predict(model_a,
                                  data.frame(MONTH = data0$MONTH[count],
                                             M.YEAR = data0$M.YEAR[count]),
                                  re.form = NA, type = "response")
    }
  }
} 

data0 <- data0 %>% 
  left_join(timeline)

birding_time <- data0


# hotspot birding ####


# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0 <- data0_MY_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

count <- 0
for (i in 1:4) { # for each state
  
  data_a <- data0_MY_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, CELL.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(HOTSPOT = if_else(LOCALITY.TYPE == "H", 1, 0)) %>%
    ungroup() 
  
  tictoc::tic(glue::glue("GLMM for {anal_states[, i]}"))
  model_a <- glmer(HOTSPOT ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
                   data = data_a, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  for (j in 1:12) {
    for (k in 1:4) {
      count <- count + 1
      
      data0$PRED[count] = predict(model_a,
                                  data.frame(MONTH = data0$MONTH[count],
                                             M.YEAR = data0$M.YEAR[count]),
                                  re.form = NA, type = "response")
    }
  }
} 

data0 <- data0 %>% 
  left_join(timeline)

hotspot_birding <- data0

# birding protocol ####


# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0 <- data0_MY_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

count <- 0
for (i in 1:4) { # for each state
  
  # using 300m as threshold to define travelling birding
  data_a <- data0_MY_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, CELL.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(TRAVELLING = if_else((PROTOCOL.TYPE == "Traveling" & 
                                             EFFORT.DISTANCE.KM > 0.3), 1, 0)) %>%
    ungroup() 
  
  tictoc::tic(glue::glue("GLMM for {anal_states[, i]}"))
  model_a <- glmer(TRAVELLING ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
                   data = data_a, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  for (j in 1:n_distinct(data0$MONTH)) {
    for (k in 1:n_distinct(data0$M.YEAR)) {
      count <- count + 1
      
      data0$PRED[count] = predict(model_a,
                                  data.frame(MONTH = data0$MONTH[count],
                                             M.YEAR = data0$M.YEAR[count]),
                                  re.form = NA, type = "response")
    }
  }
} 

data0 <- data0 %>% 
  left_join(timeline)

birding_protocol <- data0


# spatial spread ####


# tidy version of 25*25 grid to which values of interest can be joined (for spatial object)
gridmapg1_sf <- gridmapg1 %>% 
  st_as_sf() %>% 
  rename(CELL.ID = id)

load("00_data/maps_gridmapg1_IN.RData") # all cells within India


data_SoIBmap0 <- data0_MY_slice_G %>% 
  # selecting only nine months (DUR has only nine)
  filter(MONTH %in% month_compar$MONTH) %>% 
  group_by(COVID, CELL.ID, MONTH) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(GROUP.ID)) %>% 
  # summarise removes one grouping variable so no need to repeat group_by() 
  dplyr::summarise(NO.LISTS = mean(NO.LISTS)) %>% 
  ungroup() %>% 
  inner_join(gridmapg1_IN) %>% 
  st_as_sf()
# not filling zero lists for cell.id-month combinations because we are interested in the 
# average monthly effort (when effort present)

# all grid cells covered in whole time period (four years) 
# (not all grids in India)
# this is enough to show difference/change; other grids will be blank
foc_grids <- data_SoIBmap0 %>% distinct(CELL.ID)

# only grid cells covered in all four periods
foc_grids1 <- data_SoIBmap0 %>% 
  group_by(CELL.ID) %>% 
  summarise(N.PERIODS = n_distinct(COVID)) %>% 
  filter(N.PERIODS == 4) %>% 
  distinct(CELL.ID)

# only non-zero grid cells but period-wise
foc_grids2 <- data_SoIBmap0 %>% distinct(COVID, CELL.ID)


data_SoIBmap0_filt1 <- data_SoIBmap0 %>% right_join(foc_grids1)
data_SoIBmap0_filt2 <- data_SoIBmap0 %>% right_join(foc_grids2)


data_SoIBmap1 <- data_SoIBmap0_filt2 %>% filter(COVID == "BEF")

# # to plot and check
# data_SoIBmap1 %>%
#   ggplot() +
#   geom_polygon(data = indiamap,
#                aes(long, lat, group = group),
#                colour = "black", fill = NA, size = 0.2) +
#   geom_sf(aes(fill = NO.LISTS), col = NA)

data_SoIBmap2 <- data_SoIBmap0_filt2 %>% 
  filter(COVID %in% c("DUR_20","DUR_21")) %>% 
  mutate(COVID = "DUR") %>% 
  # averages for cells with lists both years, but returns single year value for others (n/1)
  group_by(COVID, CELL.ID) %>% 
  filter(NO.LISTS != 0) %>% 
  dplyr::summarise(NO.LISTS = mean(NO.LISTS)) %>% 
  ungroup()


data_SoIBmap3 <- data_SoIBmap0_filt2 %>% filter(COVID == "AFT")




require(spdep)

#### analysing three periods separately

clust_data <- pw_spread_LMCOA(data_SoIBmap1, data_SoIBmap2, data_SoIBmap3,
                              trans.resp = T)


# saving RData ####

save(anal_states, bang1, bang2, birding_protocol, birding_time, clust_data, coverage, 
     group_birding, hotspot_birding, indiamap, obs_mainstates, site_fidelity, 
     species, species_guwa, species_koch, urban1,
     file = "00_data/NCF_AM2022_covid.RData")
