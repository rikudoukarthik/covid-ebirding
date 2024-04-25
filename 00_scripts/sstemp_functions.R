expandbyspecies = function(data, species) {
  
  require(tidyverse)
  
  # considers only complete lists
  
  occinfo <- data %>% 
    group_by(GROUP.ID, COMMON.NAME) %>% 
    reframe(REPORT = 1) %>% 
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
    left_join(checklistinfo, by = "GROUP.ID") %>% 
    # for species not reported in lists, filling in NAs in COMMON.NAME and REPORT
    mutate(REPORT = replace_na(REPORT, 0))

  return(expanded)
  
}



singlespeciesmodel = function(data, species, specieslist) {
  
  require(tidyverse)
  require(lme4)
  require(merTools)


  # getting median list length for prediction later
  median_length <- data %>% 
    distinct(M.YEAR, MONTH, GROUP.ID, NO.SP) %>% 
    group_by(MONTH) %>% 
    reframe(NO.SP.MED = floor(median(NO.SP)))
  
  message("Completed preparations for modelling. Now starting modelling.")

  
  data_filt = data %>%
    filter(COMMON.NAME == species) %>%
    # using only CELL.ID-MONTH (space-time) combos in which species occurs
    distinct(CELL.ID, MONTH) %>% 
    # joining median list length
    left_join(median_length, by = "MONTH") %>% 
    # join rest of dataset back
    left_join(data)

  # expand dataframe to include absences as well
  data_exp = expandbyspecies(data_filt, species) %>% 
    # join species categories
    left_join(specieslist %>% distinct(COMMON.NAME, SP.CATEGORY),
              by = "COMMON.NAME")
  
  
  # the model ---------------------------------------------------------------
  
  tictoc::tic(glue("GLMM for {mt}, {species}"))
  
  # for some species in KL and MH, using cloglog link is resulting in "PIRLS step-halvings
  # failed to reduce deviance in pwrssUpdate"
  
  if ((state_name == "Kerala" & species %in% fail_spec_KL) |
      (state_name == "Maharashtra" & species %in% fail_spec_MH)){
    
    model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:log(NO.SP) + MONTH:M.YEAR +
                          (1|CELL.ID),
                        data = data_exp, family = binomial,
                        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
    
  } else {
    
    model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:log(NO.SP) + MONTH:M.YEAR +
                          (1|CELL.ID),
                        data = data_exp, family = binomial(link = "cloglog"),
                        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
    
  }
  
  tictoc::toc()
  
  # predicting from model ---------------------------------------------------
  
  # prepare a new data file to predict
  birds_pred <- data_exp %>% 
    # selecting random CELL.ID because predictInterval needs input even if which == "fixed"
    mutate(CELL.ID = sample(unique(CELL.ID), 1)) %>% 
    distinct(MONTH, M.YEAR, CELL.ID) %>% 
    # joining median list length
    left_join(median_length, by = "MONTH") %>% 
    rename(NO.SP = NO.SP.MED)
  
  
  pred = predictInterval(model_spec, newdata = birds_pred, which = "fixed",
                         level = 0.48, type = "linear.prediction")
  birds_pred$PRED.LINK = pred$fit
  birds_pred$SE.LINK = pred$fit - pred$lwr
  
  
  birds_pred = birds_pred %>%
    filter(!is.na(PRED.LINK) & !is.na(SE.LINK)) %>%
    group_by(M.YEAR) %>% 
    reframe(PRED.LINK = mean(PRED.LINK), 
            SE.LINK = mean(SE.LINK)) 
  
  return(birds_pred)
  
}

