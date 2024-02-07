expandbyspecies = function(data, species) {
  
  require(tidyverse)
  
  # considers only complete lists
  
  occinfo <- data %>% 
    group_by(GROUP.ID, COMMON.NAME) %>% 
    reframe(REPORT = as.numeric(max(OBSERVATION.COUNT))) %>% 
    ungroup()
  
  checklistinfo = data %>% 
    filter(ALL.SPECIES.REPORTED == 1) %>%
    distinct(GROUP.ID, M.YEAR, MONTH, DAY.M, URBAN, CELL.ID, SUBCELL.ID, NO.SP) %>%
    group_by(GROUP.ID) %>% 
    slice(1) %>% 
    ungroup()
  
  # expand data frame to include the bird species in every list
  expanded = checklistinfo %>% 
    # getting species of interest
    group_by(GROUP.ID) %>% 
    reframe(COMMON.NAME = species) %>% 
    left_join(occinfo, by = c("GROUP.ID", "COMMON.NAME")) %>%
    # deal with NAs (column is character)
    # for species not reported in lists, filling in NAs in COMMON.NAME and REPORT
    mutate(REPORT = replace_na(REPORT, "0"))

  return(expanded)
  
}



singlespeciesmodel = function(data, species, specieslist) {
  
  require(tidyverse)
  require(lme4)
  require(merTools)


 
  
  
  # getting median list length for prediction later
  median_length <- data_occ0 %>% 
    distinct(MONTHS.TYPE, M.YEAR, MONTH, GROUP.ID, NO.SP) %>% 
    # interested in seasonality so not concerned with M.YEAR
    group_by(MONTHS.TYPE, MONTH) %>% 
    dplyr::summarise(NO.SP.MED = floor(median(NO.SP)))
  
  # dataframe with empty column to populate with looped values
  # total rows: product of distinct values of predictors
  birds_pred <- data_occ0 %>% 
    group_by(MONTHS.TYPE) %>% 
    # nest month within species cos all species not present all-year
    tidyr::expand(COMMON.NAME, nesting(MONTH), M.YEAR) %>% 
    left_join(cur_species_list) %>% 
    # joining median list length
    left_join(median_length) %>% 
    rename(NO.SP = NO.SP.MED) %>% 
    mutate(PRED.LINK = NA,
           SE.LINK = NA)
  
  message("Completed preparations for modelling. Now starting modelling.")
    


  data1 = data1 %>%
    filter(COMMON.NAME == species) %>%
    distinct(gridg3, month) %>% 
    left_join(data1)
  
  tm = data1 %>% distinct(timegroups)
  #rm(data, pos = ".GlobalEnv")
  
  datay = data1 %>%
    group_by(gridg3, gridg1, group.id) %>% 
    slice(1) %>% 
    group_by(gridg3, gridg1) %>% 
    reframe(medianlla = median(no.sp)) %>%
    group_by(gridg3) %>% 
    reframe(medianlla = mean(medianlla)) %>%
    reframe(medianlla = round(mean(medianlla)))
  # getting median list length for prediction later
  median_length <- data_occ0 %>% 
    distinct(M.YEAR, MONTH, GROUP.ID, NO.SP) %>% 
    # interested in seasonality so not concerned with M.YEAR
    group_by(MONTH) %>% 
    dplyr::summarise(NO.SP.MED = floor(median(NO.SP)))
  
  medianlla = datay$medianlla
  
  
  # expand dataframe to include absences as well
  ed = expandbyspecies(data1, species) %>% 
    # join species categories
    left_join(specieslist)
  
  
  # the model ---------------------------------------------------------------
  
    m1 = glmer(OBSERVATION.COUNT ~ month + month:log(no.sp) + timegroups + (1|gridg1), 
               data = ed, family = binomial(link = 'cloglog'), 
               nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))

  
  # predicting from model ---------------------------------------------------
  
  # prepare a new data file to predict
  ltemp <- ed %>% 
    group_by(month) %>% 
    reframe(timegroups = unique(tm$timegroups)) %>% 
    mutate(no.sp = medianlla,
           # <annotation_pending_AV> why taking 1st value?
           gridg1 = data1$gridg1[1], 
           gridg3 = data1$gridg3[1])
  
  f2 <- ltemp %>% 
    dplyr::select(timegroups) %>% 
    # this is not actually needed
    mutate(freq = 0, se = 0)
  
  
  if (flag != 2)
  {
    #pred = predict(m1, newdata = ltemp, type = "response", re.form = NA, allow.new.levels=TRUE)
    pred = predictInterval(m1, newdata = ltemp, which = "fixed",
                           level = 0.48, type = "linear.prediction")
    f2$freqt = pred$fit
    f2$set = pred$fit-pred$lwr
  }
  
  if (flag == 2)
  {
    pred = predict(m1, newdata = ltemp, type = "link", se.fit = T)
    f2$freqt = pred$fit
    f2$set = pred$se.fit
  }
  
  f1 = f2 %>%
    filter(!is.na(freqt) & !is.na(se)) %>%
    group_by(timegroups) %>% 
    reframe(freq = mean(freqt), se = mean(set)) %>% 
    right_join(tm) %>% 
    left_join(databins %>% distinct(timegroups, year)) %>% 
    rename(timegroupsf = timegroups,
           timegroups = year) %>% 
    mutate(timegroupsf = factor(timegroupsf, 
                                levels = c("before 2000","2000-2006","2007-2010",
                                           "2011-2012","2013","2014","2015","2016",
                                           "2017","2018","2019","2020","2021","2022"))) %>% 
    complete(timegroupsf) %>% 
    arrange(timegroupsf)
  
  
  tocomb = c(species, f1$freq, f1$se)
  return(tocomb)
  # each species's tocomb becomes one column in final trends0 output object
  
}