# trends files
trends <- list.files(path = trends_pathonly, 
                     # Generate the full file paths
                     full.names = T) %>% 
  # Read each CSV file and combine them into a single data frame
  map_df(read.csv) 



totsims = n_distinct(trends$sl)


# data filtering: problem species -----------------------------------------

# long-term (problematic species) ###

trendsa = trends %>%
  group_by(sl, COMMON.NAME) %>%  
  # is any simulation of any species problematic: unable to calc. SE or SE > |mean|
  filter(any(is.na(se) | se > abs(freq))) %>%
  distinct(sl, COMMON.NAME)





###

birds_pred <- birds_pred %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = TRUE),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = TRUE)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  left_join(timeline) %>% 
  # summarising for species categories
  group_by(STATE, MONTHS.TYPE, M.YEAR, SP.CATEGORY) %>% 
  dplyr::summarise(PRED = mean(PRED),
                   # propagating SE across species of a category
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PRED - 1.96*SE,
                   CI.U = PRED + 1.96*SE)
