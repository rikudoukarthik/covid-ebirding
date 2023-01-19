


## month type 2

data_mtype <- data_occ0 %>% 
  filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[2])

birds_pred0 <- birds_pred %>% 
  filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[2]) 

for (i in 1:n_distinct(birds_pred0$COMMON.NAME)) {
  
  model_test(birds_pred0 = birds_pred0,
             data_mtype = data_mtype,
             data_occ = data_occ,
             cur_i = i, cur_m = 2)    
  
  birds_pred <- birds_pred %>% 
    left_join(birds_pred0, by = c("MONTHS.TYPE", "COMMON.NAME", "MONTH", "M.YEAR", "STATE",
                                  "SP.CATEGORY", "NO.SP")) %>% 
    mutate(PRED.LINK = coalesce(PRED.LINK.x, PRED.LINK.y),
           SE.LINK = coalesce(SE.LINK.x, SE.LINK.y)) %>% 
    dplyr::select(-PRED.LINK.x, -PRED.LINK.y, -SE.LINK.x, -SE.LINK.y)
  
}

## month type 3

data_mtype <- data_occ0 %>% 
  filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[3])

birds_pred0 <- birds_pred %>% 
  filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[3]) 

for (i in 1:n_distinct(birds_pred0$COMMON.NAME)) {
  
  model_test(birds_pred0 = birds_pred0,
             data_mtype = data_mtype,
             data_occ = data_occ,
             cur_i = i, cur_m = 3)    
  
  birds_pred <- birds_pred %>% 
    left_join(birds_pred0, by = c("MONTHS.TYPE", "COMMON.NAME", "MONTH", "M.YEAR", "STATE",
                                  "SP.CATEGORY", "NO.SP")) %>% 
    mutate(PRED.LINK = coalesce(PRED.LINK.x, PRED.LINK.y),
           SE.LINK = coalesce(SE.LINK.x, SE.LINK.y)) %>% 
    dplyr::select(-PRED.LINK.x, -PRED.LINK.y, -SE.LINK.x, -SE.LINK.y)
  
}

## combining all three

  birds_pred <- birds_pred %>% 
    mutate(PRED = clogloglink(PRED.LINK, inverse = T),
           # to transform lower bound of SE (not CI.L! think "mean +- SE")
           SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
    mutate(SE = PRED - SE.L) %>% 
    left_join(timeline) %>% 
    # summarising for species categories
    group_by(STATE, MONTHS.TYPE, M.YEAR, SP.CATEGORY) %>% 
    summarise(PRED = mean(PRED),
              # propagating SE across species of a category
              SE = sqrt(sum((SE)^2))/n(),
              CI.L = PRED - 1.96*SE,
              CI.U = PRED + 1.96*SE)
  
  tictoc::toc()
  


