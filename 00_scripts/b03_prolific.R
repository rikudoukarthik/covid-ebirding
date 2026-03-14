# Model run separately for each location, and species variables included in model 


require(tidyverse)

source("00_scripts/b03_model_functions.R")
source("00_scripts/functions.R")
load("00_data/bird_prolific.RData")



# getting median list length (patch-wise) for prediction later
median_length <- prol_data %>% 
  distinct(MONTHS.TYPE, LOCALITY.ID, M.YEAR, MONTH, SAMPLING.EVENT.IDENTIFIER, NO.SP) %>% 
  group_by(MONTHS.TYPE, LOCALITY.ID, MONTH) %>% 
  reframe(NO.SP.MED = floor(median(NO.SP)))

# dataframe with empty column to populate with looped values
prolific_pred <- prol_data %>% 
  group_by(MONTHS.TYPE, LOCALITY.ID) %>% 
  tidyr::expand(nesting(MONTH), M.YEAR, SP.CATEGORY) %>% 
  ungroup() %>% 
  # joining median list length
  left_join(median_length, by = c("MONTHS.TYPE", "LOCALITY.ID", "MONTH")) %>% 
  rename(NO.SP = NO.SP.MED) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         ID = row_number()) # for correct join of iters

to_iter <- prolific_pred %>% distinct(MONTHS.TYPE, LOCALITY.ID)
# # to test with subset of iters
# %>% 
# group_by(MONTHS.TYPE) %>% 
# slice(1:2)


prolific_pred <- map2(to_iter$MONTHS.TYPE, to_iter$LOCALITY.ID, 
                      ~ fit_pred_prolific(.x, .y, 
                                          prol_data, prolific_pred)) %>%
  list_rbind() %>%
  # merge back into main prediction table
  right_join(prolific_pred, 
             by = c("MONTHS.TYPE","LOCALITY.ID","MONTH","M.YEAR","SP.CATEGORY", "NO.SP", "ID")) %>%
  mutate(PRED.LINK = coalesce(PRED.LINK, PRED.LINK2),
         SE.LINK   = coalesce(SE.LINK, SE.LINK2)) %>%
  select(-PRED.LINK2, -SE.LINK2)


# for data points

prolific_points <- prolific_pred %>% 
  left_join(timeline) %>% 
  # averaging across months of year
  group_by(MONTHS.TYPE, M.YEAR, SP.CATEGORY, LOCALITY.ID) %>% 
  summarise_mean_and_se(PRED.LINK, SE.LINK, n_is_sim = FALSE) %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         CI.L = clogloglink((PRED.LINK - 1.96*SE.LINK), inverse = T),
         CI.U = clogloglink((PRED.LINK + 1.96*SE.LINK), inverse = T)) %>% 
  mutate(SE = (PRED - CI.L)/1.96)

prolific_points <- prolific_points %>% 
  arrange(MONTHS.TYPE, SP.CATEGORY, LOCALITY.ID, M.YEAR) %>%
  # getting trends values of first year (2018--19)
  group_by(MONTHS.TYPE, SP.CATEGORY, LOCALITY.ID) %>% 
  mutate(PRED.LINK.Y1 = first(PRED.LINK),
         PRED.Y1 = first(PRED),
         SE.LINK.Y1 = first(SE.LINK)) %>% 
  ungroup() %>% 
  # for calculating as % change
  #   back-transformed so value is % of year1 value
  mutate(PRED.PERC = 100*PRED/PRED.Y1) 

# calculting simulated CIs for the % change
set.seed(10) 
prolific_points <- prolific_points %>% 
  group_by(MONTHS.TYPE, SP.CATEGORY, LOCALITY.ID, M.YEAR) %>% 
  # 1000 simulations of transformed ratio of present:original values
  # quantiles*100 from these gives us our CI limits for PRED.PERC
  reframe(SIM.RATIOS = simerrordiv(PRED.LINK, PRED.LINK.Y1, SE.LINK, SE.LINK.Y1,
                                   state = NULL, species = NULL) %>% 
            pull(rat)) %>% 
  group_by(MONTHS.TYPE, SP.CATEGORY, LOCALITY.ID, M.YEAR) %>% 
  reframe(CI.L.PERC = 100*as.numeric(quantile(SIM.RATIOS, 0.025)),
          CI.U.PERC = 100*as.numeric(quantile(SIM.RATIOS, 0.975))) %>% 
  right_join(prolific_points, by = c("MONTHS.TYPE", "SP.CATEGORY", "LOCALITY.ID", "M.YEAR")) %>%
  group_by(MONTHS.TYPE, SP.CATEGORY, LOCALITY.ID) %>% 
  # making CI band zero for first year
  mutate(CI.L.PERC = case_when(M.YEAR == first(M.YEAR) ~ PRED.PERC,
                               TRUE ~ CI.L.PERC),
         CI.U.PERC = case_when(M.YEAR == first(M.YEAR) ~ PRED.PERC,
                               TRUE ~ CI.U.PERC)) %>%
  ungroup() %>%
  dplyr::select(MONTHS.TYPE, SP.CATEGORY, LOCALITY.ID, M.YEAR, CI.L, PRED, CI.U,
                CI.L.PERC, PRED.PERC, CI.U.PERC) %>% 
  # calculating SE to propagate
  mutate(SE.PERC = (CI.U.PERC - PRED.PERC)/1.96)

# to make lines connect per location                      
prolific_points <- prolific_points %>% 
  arrange(MONTHS.TYPE, LOCALITY.ID, SP.CATEGORY, M.YEAR) %>% 
  mutate(PATH.GROUP = glue("{LOCALITY.ID}{SP.CATEGORY}"))


# summarising and propagating SE for species categories (across locations)
prolific_pred <- prolific_points %>% 
  group_by(MONTHS.TYPE, M.YEAR, SP.CATEGORY) %>% 
  summarise_mean_and_se(PRED.PERC, SE.PERC, n_is_sim = FALSE) |> 
  mutate(CI.L.PERC = PRED.PERC - 1.96*SE.PERC,
         CI.U.PERC = PRED.PERC + 1.96*SE.PERC)


# save results
save(prolific_pred, prolific_points, 
     file = glue("00_outputs/bird_models/prolific.RData"))
