# Step 4 of subsampling: resolve the 1000 predictions to get bootstrapped
# estimates

require(tidyverse)
require(glue)
require(VGAM) #clogloglink()

source("00_scripts/b_model_functions.R")


for (mt in c("LD", "NL")) {

  mod_path <- glue("00_outputs/bird_models/{state_name}/b03_models_{mt}.csv")
  
  # read data of model predictions
  data_mod <- read_csv(mod_path) %>% 
    mutate(M.YEAR = factor(M.YEAR))
  
  
  # resolve bootstrapped estimates
  
  data_res = data_mod %>%
    # across iterations
    group_by(COMMON.NAME, SP.CATEGORY, M.YEAR) %>% 
    summarise_mean_and_se(PRED.LINK, SE.LINK) |> 
    group_by(COMMON.NAME, SP.CATEGORY, M.YEAR) %>% 
    mutate(CI.L = inverse_link(PRED.LINK - 1.96*SE.LINK, state_name, COMMON.NAME),
           PRED = inverse_link(PRED.LINK, state_name, COMMON.NAME),
           CI.U = inverse_link(PRED.LINK + 1.96*SE.LINK, state_name, COMMON.NAME)) %>%
    ungroup()
  
  data_res <- data_res %>% 
    arrange(COMMON.NAME, M.YEAR) %>%
    # getting trends values of first year (2018--19)
    group_by(COMMON.NAME, SP.CATEGORY) %>% 
    mutate(PRED.LINK.Y1 = first(PRED.LINK),
           PRED.Y1 = first(PRED),
           SE.LINK.Y1 = first(SE.LINK)) %>% 
    ungroup() %>% 
    # for calculating as % change
    #   back-transformed so value is % of year1 value
    mutate(PRED.PERC = 100*PRED/PRED.Y1) 
  
  # calculting simulated CIs for the % change
  set.seed(10) 
  data_res <- data_res %>% 
    group_by(COMMON.NAME, SP.CATEGORY, M.YEAR) %>%
    # 1000 simulations of transformed ratio of present:original values
    # quantiles*100 from these gives us our CI limits for PRED.PERC
    reframe(SIM.RATIOS = simerrordiv(PRED.LINK, PRED.LINK.Y1, SE.LINK, SE.LINK.Y1,
                                     state_name, COMMON.NAME) %>% 
              pull(rat)) %>% 
    group_by(COMMON.NAME, SP.CATEGORY, M.YEAR) %>%
    reframe(CI.L.PERC = 100*as.numeric(quantile(SIM.RATIOS, 0.025)),
            CI.U.PERC = 100*as.numeric(quantile(SIM.RATIOS, 0.975))) %>% 
    right_join(data_res, by = c("COMMON.NAME", "SP.CATEGORY", "M.YEAR")) %>%
    group_by(COMMON.NAME, SP.CATEGORY) %>%
    # making CI band zero for first year
    mutate(CI.L.PERC = case_when(M.YEAR == first(M.YEAR) ~ PRED.PERC,
                                 TRUE ~ CI.L.PERC),
           CI.U.PERC = case_when(M.YEAR == first(M.YEAR) ~ PRED.PERC,
                                 TRUE ~ CI.U.PERC)) %>%
    ungroup() %>%
    dplyr::select(M.YEAR, COMMON.NAME, SP.CATEGORY, CI.L, PRED, CI.U,
                  CI.L.PERC, PRED.PERC, CI.U.PERC)
  
  
  # averaging per species category
  data_res_avg <- data_res %>% 
    # calculating SE to propagate
    mutate(SE.PERC = (CI.U.PERC - PRED.PERC)/1.96,
           MONTHS.TYPE = mt) %>% 
    group_by(MONTHS.TYPE, M.YEAR, SP.CATEGORY) %>% 
    reframe(PRED.PERC = mean(PRED.PERC),
            # propagating SE across species of a category
            # ADD sd of PRED.PERC to sqrt()
            SE.PERC = sqrt(sum((SE.PERC)^2))/n(),
            CI.L.PERC = PRED.PERC - 1.96*SE.PERC,
            CI.U.PERC = PRED.PERC + 1.96*SE.PERC)
  
  
  # writing
  save(data_res_avg, data_res, 
       file = glue("00_outputs/bird_models/{state_name}/b04_results_{mt}.RData"))
  
}
