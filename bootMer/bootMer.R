data0 <- data0_MY_slice_S %>% 
  filter(STATE %in% anal_states[1, 4]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED = NA,
         SE = NA,
         STATE = factor(STATE, levels = anal_states[1,4])) %>% 
  arrange(STATE)

count <- 0
for (i in 4) { # for Assam
  
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
      
      tictoc::tic(glue::glue("Bootstrapped predictions for {anal_states[, i]} in month {j}, year {k}"))
      prediction <- boot_conf_GLMM(model_a,
                                   new_data = data.frame(MONTH = data0$MONTH[count],
                                                         M.YEAR = data0$M.YEAR[count]),
                                   nsim = 100)
      tictoc::toc()
      
      data0$PRED[count] <- median(na.omit(prediction$t[count]))
      data0$SE[count] <- sd(na.omit(prediction$t[count]))
      
    }
  }
} 

save(data0, file = "bootMer/data0.RData")
