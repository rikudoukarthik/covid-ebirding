# this script will be sourced instead of functionising some of these steps and directly 
# calling them from the main script, because parallelisation is somehow much less
# efficient when within a function or loop.



tictoc::tic("Total time elapsed for bird_model_state():")

require(lme4)

cur_species_list <- species_list %>% filter(STATE == state_name)

### modelling directly with presence-absence data instead of relative abundance

# to join for presence-absence of various species
temp1 <- data_full %>% 
  filter(STATE == state_name) %>% 
  group_by(GROUP.ID, COMMON.NAME) %>% 
  summarise(OBSERVATION.COUNT = max(OBSERVATION.COUNT)) %>% 
  ungroup()

# to later join checklist metadata
temp2 <- data_sliceG %>% 
  filter(STATE == state_name) %>% 
  arrange(SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::select(GROUP.ID, STATE, COUNTY, LOCALITY, LATITUDE, LONGITUDE, OBSERVATION.DATE, 
                M.YEAR, MONTH, DAY.M, M.YEAR, URBAN, CELL.ID, SUBCELL.ID, NO.SP)

data_occ <- data_sliceG %>% 
  filter(STATE == state_name) %>% 
  group_by(GROUP.ID) %>% 
  summarise(COMMON.NAME = cur_species_list$COMMON.NAME) %>% 
  left_join(temp1) %>% 
  # for species not reported in lists, filling in NAs in COMMON.NAME and REPORT
  mutate(REPORT = replace_na(OBSERVATION.COUNT, "0")) %>% 
  dplyr::select(-OBSERVATION.COUNT) %>% 
  mutate(REPORT = as.numeric(case_when(REPORT != "0" ~ "1", TRUE ~ REPORT))) %>% 
  # checklist metadata
  left_join(temp2, by = "GROUP.ID") %>% 
  arrange(GROUP.ID) %>% 
  # species categories
  left_join(cur_species_list) %>% 
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
  left_join(cur_species_list) %>% 
  # joining median list length
  left_join(median_length) %>% 
  rename(NO.SP = NO.SP.MED) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

print("Completed preparations for modelling. Now starting modelling.")




# month type 1 ------------------------------------------------------------



data_mtype <- data_occ0 %>% 
  filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[1])

birds_pred0 <- birds_pred %>% 
  filter(MONTHS.TYPE == unique(birds_pred$MONTHS.TYPE)[1]) 


birds_pred0_b <- birds_pred0 %>% 
  filter(COMMON.NAME == unique(birds_pred0$COMMON.NAME)[1]) %>% 
  rename(PRED.LINK2 = PRED.LINK,
         SE.LINK2 = SE.LINK)

assign("birds_pred0_b", birds_pred0_b, envir = .GlobalEnv)

data_spec <- data_mtype %>% 
  filter(COMMON.NAME == unique(birds_pred0$COMMON.NAME)[1]) %>% 
  # filtering for CELL.ID-MONTH (space-time) combos in which species occurs
  filter(REPORT == 1) %>% 
  distinct(COMMON.NAME, CELL.ID, MONTH) %>% 
  left_join(data_occ)


tictoc::tic(glue("GLMM for months type {1}, {unique(birds_pred0$COMMON.NAME)[1]}"))
model_spec <- glmer(REPORT ~ M.YEAR + MONTH:NO.SP + MONTH:M.YEAR + 
                      (1|CELL.ID),
                    data = data_spec, family = binomial(link = "cloglog"),
                    nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
tictoc::toc() 

source("scripts/b01_overall_model.R")


tictoc::tic(glue("Bootstrapped predictions for months type {1}, {unique(birds_pred0$COMMON.NAME)[1]}"))
prediction <- split_par_boot(model = model_spec, 
                             new_data = birds_pred0_b, 
                             new_data_string = "birds_pred0_b", 
                             mode = "normal")
tictoc::toc() 


count <- 0
for (j in 1:n_distinct(birds_pred0_b$MONTH)) {
  
  for (k in 1:n_distinct(birds_pred0_b$M.YEAR)) {
    
    count <- count + 1
    
    birds_pred0_b$PRED.LINK2[count] = median(na.omit(prediction[,count]))
    birds_pred0_b$SE.LINK2[count] = sd(na.omit(prediction[,count]))
    
  }
}

birds_pred0 <- birds_pred0 %>% 
  left_join(birds_pred0_b) %>% 
  mutate(PRED.LINK = coalesce(PRED.LINK, PRED.LINK2),
         SE.LINK = coalesce(SE.LINK, SE.LINK2)) %>% 
  dplyr::select(-PRED.LINK2, -SE.LINK2)

assign("birds_pred0", birds_pred0, envir = .GlobalEnv)  

birds_pred <- birds_pred %>% 
  left_join(birds_pred0, by = c("MONTHS.TYPE", "COMMON.NAME", "MONTH", "M.YEAR", "STATE",
                                "SP.CATEGORY", "NO.SP")) %>% 
  mutate(PRED.LINK = coalesce(PRED.LINK.x, PRED.LINK.y),
         SE.LINK = coalesce(SE.LINK.x, SE.LINK.y)) %>% 
  dplyr::select(-PRED.LINK.x, -PRED.LINK.y, -SE.LINK.x, -SE.LINK.y)
