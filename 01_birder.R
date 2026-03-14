# Setup -------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(patchwork)
library(boot)
library(lme4)
library(glue)
library(spdep)
library(sf)
library(parallel) # for bootstrapping error from GLMMs
library(colorspace) # for ability to set midpoint in ggplot scale colour
library(ggallin) # for inverse arcsine transformation of ggplot scale

theme_set(theme_classic())

covid_palette <- c("#1B9E77", "#EF4050", "#E89005", "#9678B6")
# covid_palette2 <- c("#1B9E77", "#D95F02", "#7570B3", "#555555")
covid_palette3 <- c("#1B9E77", "#E89005", "#9678B6")


anal_states <- data.frame(a = "Karnataka",
                          b = "Kerala",
                          c = "Maharashtra",
                          d = "Assam")


load("00_data/maps_sf.RData")



source("00_scripts/functions.R")

# Only needed once in the beginning to obtain the data.

# Getting MODIS data ###
if (!file.exists("00_data/in_LULC_MODIS/in_LULC_MODIS.tif") & 
    !file.exists("00_data/rast_UNU.RData")) {
  
  getmodisdata()
  
} else {
  print("MODIS data is ready to use!")
}

# #  Getting EBD data + filtering and preparing for current use (62 mins on server) ###
# tictoc::tic("Preparing data for analyses")
# rawdatapath <- "00_data/ebd_IN_prv_relMay-2022.txt"
# senspath <- "00_data/ebd_sensitive_relMay-2022_IN.txt"
# datapath <- "00_data/ebd_IN_relMay-2022.RData"
# groupaccspath <- "00_data/ebd_users_GA_relMay-2022.csv"
# covidclasspath <- "00_data/covid_classification.csv"
# rast_UNU_path <- "00_data/rast_UNU.RData"
# 
# data_qualfilt_prep(rawdatapath, senspath, datapath, groupaccspath, covidclasspath,
#                    rast_UNU_path,
#                    maxvel = 20, minsut = 2)
# tictoc::toc()


load("00_data/data0_MY_d_slice.RData")
load("00_data/grids_st_sf.RData")



# timeline metadata
timeline <- data0_MY_d_slice_G %>% distinct(COVID, YEAR, MONTH, M.YEAR)

# finding main state of birders (mainly for observer-level metrics)
obs_mainstates <- data0_MY_d_slice_S %>% 
  group_by(OBSERVER.ID, STATE) %>% 
  dplyr::summarise(LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  arrange(desc(LISTS)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  rename(MAIN.STATE = STATE) %>% 
  dplyr::select(-LISTS)

# choosing only months with data from all 3 COVID categories
month_compar <- data0_MY_d_slice_G %>% 
  group_by(MONTH) %>% 
  mutate(COVID = factor(case_when(as.character(COVID) %in% c("DUR_20","DUR_21") ~ "DUR",
                                  TRUE ~ as.character(COVID)),
                        levels = c("BEF","DUR","AFT"))) %>% 
  dplyr::summarise(N = n_distinct(COVID)) %>% 
  filter(N == 3) %>%  # all 3 COVID categories
  dplyr::select(MONTH)

# removing unnecessary objects from maps_sf.RData
rm(g2_in_sf, g3_in_sf, g4_in_sf,
   g2_st_sf, g3_st_sf, g4_st_sf)

# Group birding (\>2 observer) per observer -------------------------------

anal_name <- "d01_group_po"

# Final analysis: national GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(GROUP.BIRDING = if_else(NUMBER.OBSERVERS > 2, 1, 0)) %>%
  ungroup()

tictoc::tic(glue("GLMM for India"))
model_a <- glmer(GROUP.BIRDING ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                 data = data_a, family = binomial(link = "cloglog"),
                 nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
tictoc::toc()
# 160 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 137 min


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Final analysis: statewise GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(GROUP.BIRDING = if_else(NUMBER.OBSERVERS > 2, 1, 0)) %>%
    ungroup() 
  
  tictoc::tic(glue("GLMM for {anal_states[, i]}"))
  model_b <- glmer(GROUP.BIRDING ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_b, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc()
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted group birding proportion") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted group birding proportion") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")


# Site fidelity per observer ----------------------------------------------

anal_name <- "d02_fidelity_po"

# Final analysis: national GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID) %>% 
  dplyr::summarise(NO.SITES = log(n_distinct(CELL.ID))) %>% 
  # log-transforming because Poisson GLMM won't be appropriate (zero is not a possible 
  # value in our response variable); so will run LMM instead
  ungroup() 
# not filling zeroes for M.YEAR-MONTH combos where observer had no list, because
# our aim is to see, from the people who were birding, how many sites did a birder 
# visit on average?
# number of sites is the focus here, not the individual observer (anyway a random effect)


tictoc::tic(glue("LMM for India"))
model_a <- lmer(NO.SITES ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                data = data_a, 
                control = lmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 6 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "normal")
tictoc::toc() # 476 sec (8 min)

for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  # this step doesn't change even though we log-transform prior to fitting the model
  # since the output of prediction is on log scale either way. 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Final analysis: statewise GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID) %>% 
    dplyr::summarise(NO.SITES = log(n_distinct(CELL.ID))) %>% 
    ungroup() 
  # not filling zeroes for M.YEAR-MONTH combos where observer had no list, because
  # our aim is to see, from the people who were birding, how many sites did a birder 
  # visit on average?
  # number of sites is the focus here, not the individual observer (anyway a random effect)
  
  
  tictoc::tic(glue("LMM for {anal_states[, i]}"))
  model_b <- lmer(NO.SITES ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                  data = data_b, 
                  control = lmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc()
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted no. of sites visited") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted no. of sites visited") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")



# Birding time per observer -----------------------------------------------

anal_name <- "d03_time_po"

# Final analysis: national LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID) %>% 
  dplyr::summarise(BIRDING.TIME = log(sum(DURATION.MINUTES))) %>% 
  # log-transforming because Poisson GLMM won't be appropriate (zero is not a possible 
  # value in our response variable); so will run LMM instead
  ungroup() 
# not filling zeroes for M.YEAR-MONTH combos where observer had no list, because
# our aim is to see, from the people who were birding, how much time was spent?



tictoc::tic(glue("LMM for India"))
model_a <- lmer(BIRDING.TIME ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                data = data_a, 
                control = lmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 6 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "normal")
tictoc::toc() # 471 sec (8 min)

for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Final analysis: statewise LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE %in% anal_states[1,]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID) %>% 
    dplyr::summarise(BIRDING.TIME = log(sum(DURATION.MINUTES))) %>% 
    ungroup()
  # not filling zeroes for M.YEAR-MONTH combos where observer had no list, because
  # our aim is to see, from the people who were birding, how much time was spent?
  
  
  tictoc::tic(glue("LMM for {anal_states[, i]}"))
  model_b <- lmer(BIRDING.TIME ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                  data = data_b, 
                  control = lmerControl(optimizer = "bobyqa")) 
  tictoc::toc()
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc()
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted birding time per observer") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted birding time per observer") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")



# Hotspot birding ---------------------------------------------------------

anal_name <- "d04_hotspot"

# Final analysis: national GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(HOTSPOT = if_else(LOCALITY.TYPE == "H", 1, 0)) %>%
  ungroup()

tictoc::tic(glue("GLMM for India"))
model_a <- glmer(HOTSPOT ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                 data = data_a, family = binomial(link = "cloglog"),
                 nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
tictoc::toc()
# 150 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 7520 sec (2h5m)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Final analysis: statewise GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(HOTSPOT = if_else(LOCALITY.TYPE == "H", 1, 0)) %>%
    ungroup() 
  
  tictoc::tic(glue("GLMM for {anal_states[, i]}"))
  model_b <- glmer(HOTSPOT ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_b, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc() #
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() #
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted hotspot birding proportion") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted hotspot birding proportion") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")


# Birding protocol --------------------------------------------------------


anal_name <- "d05_protocol"

# Final analysis: national GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(TRAVELLING = if_else((PROTOCOL.TYPE == "Traveling" & 
                                           EFFORT.DISTANCE.KM > 0.3), 1, 0)) %>%
  ungroup()

tictoc::tic(glue("GLMM for India"))
model_a <- glmer(TRAVELLING ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                 data = data_a, family = binomial(link = "cloglog"),
                 nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 110 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 6245 sec (104 min)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Final analysis: statewise GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(TRAVELLING = if_else((PROTOCOL.TYPE == "Traveling" & 
                                             EFFORT.DISTANCE.KM > 0.3), 1, 0)) %>%
    ungroup()
  
  tictoc::tic(glue("GLMM for {anal_states[, i]}"))
  model_b <- glmer(TRAVELLING ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_b, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc() #
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() #
  
  
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted travelling birding proportion") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted travelling birding proportion") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")



# Birding distance --------------------------------------------------------

anal_name <- "d06_distance"

# Final analysis: national LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  filter(PROTOCOL.TYPE == "Traveling", !is.na(EFFORT.DISTANCE.KM), EFFORT.DISTANCE.KM > 0.3) %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(DISTANCE = log(EFFORT.DISTANCE.KM)) %>% 
  # because no other distribution family fits 
  ungroup() 


tictoc::tic(glue("LMM for India"))
model_a <- lmer(DISTANCE ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                data = data_a, 
                control = lmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 18 sec

tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 1603 sec (26 min)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  # this step doesn't change even though we log-transform prior to fitting the model
  # since the output of prediction is on log scale either way. 
  mutate(PRED = exp(PRED.LINK),
         SE.L = exp(PRED.LINK - SE.LINK)) %>%
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Final analysis: statewise LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    filter(PROTOCOL.TYPE == "Traveling", !is.na(EFFORT.DISTANCE.KM), EFFORT.DISTANCE.KM > 0.3) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(DISTANCE = log(EFFORT.DISTANCE.KM)) %>% 
    # because no other distribution family fits 
    ungroup() 
  
  tictoc::tic(glue("LMM for {anal_states[, i]}"))
  model_b <- lmer(DISTANCE ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                  data = data_b, 
                  control = lmerControl(optimizer = "bobyqa")) 
  tictoc::toc() #
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() #
  
  
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = exp(PRED.LINK),
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)



# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted birding distance") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted birding distance") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")



# List duration -----------------------------------------------------------

anal_name <- "d07_duration"

# Final analysis: national LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(DURATION = log(DURATION.MINUTES)) %>% 
  # because no other distribution family fits 
  ungroup() 


tictoc::tic(glue("LMM for India"))
model_a <- lmer(DURATION ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                data = data_a, 
                control = lmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 35 sec

tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 2937 sec (49 min)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  # this step doesn't change even though we log-transform prior to fitting the model
  # since the output of prediction is on log scale either way. 
  mutate(PRED = exp(PRED.LINK),
         SE.L = exp(PRED.LINK - SE.LINK)) %>%
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Final analysis: statewise LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(DURATION = log(DURATION.MINUTES)) %>% 
    # because no other distribution family fits 
    ungroup() 
  
  tictoc::tic(glue("LMM for {anal_states[, i]}"))
  model_b <- lmer(DURATION ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                  data = data_b, 
                  control = lmerControl(optimizer = "bobyqa")) 
  tictoc::toc() #
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() #
  
  
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = exp(PRED.LINK),
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted list duration") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted list duration") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")


# List length -------------------------------------------------------------

anal_name <- "d08_length"

# Final analysis: national LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
# taking data from all states instead of just chosen 4
data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  # one list has 1 spuh and no species (due to slice_sample)
  filter(NO.SP != 0) %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(LENGTH = log(NO.SP)) %>%
  # because no other distribution family fits 
  ungroup() 

tictoc::tic(glue("LMM for India"))
model_a <- lmer(LENGTH ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                data = data_a, 
                control = lmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 38 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 3021 sec (50 min)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  # this step doesn't change even though we log-transform prior to fitting the model
  # since the output of prediction is on log scale either way. 
  mutate(PRED = exp(PRED.LINK),
         SE.L = exp(PRED.LINK - SE.LINK)) %>%
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Final analysis: statewise LMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
    dplyr::summarise(LENGTH = log(NO.SP)) %>%
    # because no other distribution family fits 
    ungroup() 
  
  tictoc::tic(glue("LMM for {anal_states[, i]}"))
  model_b <- lmer(LENGTH ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                  data = data_b, 
                  control = lmerControl(optimizer = "bobyqa")) 
  tictoc::toc() 
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() 
  
  
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = exp(PRED.LINK),
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted list length") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted list length") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")


# Spatial spread and coverage ---------------------------------------------

anal_name <- "d09a_s_UNU"

# https://github.com/lme4/lme4/issues/281
# Issue with high number of zeroes (or hence low proportion values) and back-transformation of linear variable without accounting for random variation (Jensen's inequality effects)
# So, need to include random effects in prediction, and need to predict at response scale.

# Final analysis: national GLMMs #

data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  distinct(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER, URBAN) %>%
  ungroup()

tictoc::tic(glue("GLMM for India"))
model_a <- glmer(URBAN ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                 data = data_a, family = binomial(link = "cloglog"),
                 nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 300 sec (5 min)


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a", 
                             mode = "extra")
tictoc::toc() # 8324 sec (138 min)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Final analysis: state GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)

for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    distinct(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER, URBAN) %>%
    ungroup()
  
  tictoc::tic(glue("GLMM for {anal_states[, i]}"))
  model_b <- glmer(URBAN ~ M.YEAR + M.YEAR:MONTH + (1|OBSERVER.ID),
                   data = data_b, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc() #
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() #
  
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted urban birding proportion") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b,
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted urban birding proportion") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")



anal_name <- "d09b_s_cover"

# Final analysis: national GLMMs #

data0_a <- data0_MY_d_slice_G %>% 
  tidyr::expand(nesting(MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA)

data_a <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, CELL.ID) %>% 
  dplyr::summarise(COVERED = if_else(n_distinct(SAMPLING.EVENT.IDENTIFIER) >= 5, 1, 0)) %>% 
  # threshold to consider cell "covered"
  group_by(M.YEAR, MONTH) %>% 
  complete(CELL.ID = g1_in_sf$GRID.G1, fill = list(COVERED = 0)) %>% 
  ungroup()

tictoc::tic(glue("GLMM for India"))
model_a <- glmer(COVERED ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
                 data = data_a, family = binomial(link = "cloglog"),
                 nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
tictoc::toc() # 114 sec


tictoc::tic(glue("Bootstrapped predictions for India"))
prediction <- split_par_boot(model = model_a, 
                             new_data = data0_a, 
                             new_data_string = "data0_a")
tictoc::toc() # 2494 sec (41 min)


for (i in 1:length(data0_a$MONTH)) {
  
  data0_a$PRED.LINK[i] <- median(na.omit(prediction[,i]))
  data0_a$SE.LINK[i] <- sd(na.omit(prediction[,i]))
  
}

data0_a <- data0_a %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)

# Final analysis: statewise GLMMs #

# dataframe with empty column to populate with looped values
# total rows: product of distinct values of predictors
data0_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  tidyr::expand(nesting(STATE, MONTH, M.YEAR)) %>% 
  mutate(PRED.LINK = NA,
         SE.LINK = NA,
         STATE = factor(STATE, levels = anal_states[1,])) %>% 
  arrange(STATE)


for (i in 1:4) { # for each state
  
  data0_b2 <- data0_b %>% filter(STATE == anal_states[, i])
  
  data_b <- data0_MY_d_slice_G %>% 
    filter(STATE == anal_states[, i]) %>% 
    group_by(M.YEAR, MONTH, CELL.ID) %>% 
    dplyr::summarise(COVERED = if_else(n_distinct(SAMPLING.EVENT.IDENTIFIER) >= 5, 1, 0)) %>% 
    # threshold to consider cell "covered"
    group_by(M.YEAR, MONTH) %>% 
    complete(CELL.ID = g1_st_sf %>% 
               filter(STATE.NAME == anal_states[, i]) %>% 
               pull(GRID.G1), 
             fill = list(COVERED = 0)) %>% 
    ungroup()
  
  tictoc::tic(glue("GLMM for {anal_states[, i]}"))
  model_b <- glmer(COVERED ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
                   data = data_b, family = binomial(link = "cloglog"),
                   nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 
  tictoc::toc() #
  
  
  tictoc::tic(glue("Bootstrapped predictions for {anal_states[, i]}"))
  prediction <- boot_conf_GLMM(model_b,
                               new_data = data0_b2,
                               new_data_string = "data0_b2",
                               nsim = 1000)
  tictoc::toc() #
  
  
  for (j in 1:length(data0_b2$MONTH)) {
    
    data0_b2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
    data0_b2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
    
  }
  
  data0_b <- data0_b %>% 
    left_join(data0_b2, by = c("STATE", "MONTH", "M.YEAR")) %>% 
    transmute(MONTH = MONTH,
              M.YEAR = M.YEAR,
              STATE = STATE,
              PRED.LINK = ifelse(!is.na(PRED.LINK.y), PRED.LINK.y, PRED.LINK.x),
              SE.LINK = ifelse(!is.na(SE.LINK.y), SE.LINK.y, SE.LINK.x))
  
} 

data0_b <- data0_b %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  get_CI_lims() %>% 
  left_join(timeline)


# Saving analysis objects #

save(data0_a, data_a, model_a, data0_b, data_b, model_b, 
     file = glue("00_outputs/{anal_name}.RData"))

# Graphs #

plot_a <- ggplot(data0_a, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  labs(title = "National-level models",
       x = "Month", y = "Predicted proportion of grid cells covered") +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1, width = 0.3, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")

plot_b <- ggplot(data0_b, 
                 aes(MONTH, PRED, colour = M.YEAR)) +
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  labs(title = "State-level models",
       x = "Month", y = "Predicted proportion of grid cells covered") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U),
                size = 1.5, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")



anal_name <- "d09c_s_spread"

load(url("https://github.com/birdcountindia/ebird-datasets/raw/main/region_codes.RData"))
region_codes <- region_codes %>% 
  filter(COUNTRY == "India", 
         !is.na(COUNTY.CODE)) %>% 
  distinct(COUNTY.CODE, STATE, COUNTY)

# Preparing data #

# linking ebd state, district, latlong to grid cell (already done) and sf state, district
sf_use_s2(FALSE)
data_spat <- data0_MY_d_slice_G %>% 
  st_as_sf(coords = c("LONGITUDE", "LATITUDE")) %>%
  st_set_crs(st_crs(dists_sf)) %>% 
  st_join(dists_sf) %>% 
  st_drop_geometry() %>% 
  dplyr::select(STATE, COUNTY, 
                SAMPLING.EVENT.IDENTIFIER, GROUP.IDENTIFIER, GROUP.ID,
                DISTRICT.NAME, STATE.NAME, CELL.ID) %>% 
  left_join(region_codes, by = c("STATE", "COUNTY"))


map_ebd_sf <- map_admin_ebd_sf(data_spat)


# linking sf districts to EBD districts and grid cells
dists_ebd_sf <- data_spat %>% 
  distinct(CELL.ID, COUNTY.CODE) %>% 
  filter(!is.na(CELL.ID), !is.na(COUNTY.CODE)) %>% 
  left_join(map_ebd_sf, by = "COUNTY.CODE") %>% 
  left_join(dists_sf %>% distinct(STATE.NAME, DISTRICT.NAME, DISTRICT.GEOM),
            by = c("STATE.NAME", "DISTRICT.NAME"))


set.seed(31)
# for number of lists 
grid_dist_NL <- dists_ebd_sf %>% 
  st_drop_geometry() %>% 
  # Ensuring that a grid cell is only linked with one district
  group_by(CELL.ID) %>% 
  slice_sample(n = 1) %>% 
  group_by(COUNTY.CODE) %>% 
  mutate(NO.CELLS.TOT = n_distinct(CELL.ID)) %>% 
  ungroup() %>% 
  st_as_sf()

# for grid coverage
grid_dist_GC <- grid_dist_NL %>% 
  # considering only districts with at least 5 grid cells 
  filter(NO.CELLS.TOT >= 5)


data_slice_G <- data0_MY_d_slice_G %>% 
  # combining the two DUR years
  mutate(COVID = case_when(M.YEAR == 2018 ~ "BEF",
                           M.YEAR %in% 2019:2020 ~ "DUR",
                           M.YEAR == 2021 ~ "AFT")) %>% 
  mutate(COVID = factor(COVID, levels = c("BEF", "DUR", "AFT")))


# starting data 
temp0_NL <- data_slice_G %>% 
  inner_join(grid_dist_NL, by = c("CELL.ID")) 

temp0_GC <- data_slice_G %>% 
  inner_join(grid_dist_GC, by = c("CELL.ID"))


# districts being considered for each
dist_cons_NL <- temp0_NL %>% 
  distinct(COUNTY.CODE, NO.CELLS.TOT) %>% 
  arrange(COUNTY.CODE)

dist_cons_GC <- temp0_GC %>% 
  distinct(COUNTY.CODE, NO.CELLS.TOT) %>% 
  arrange(COUNTY.CODE)


# for number of lists per district
temp_NL <- temp0_NL %>% 
  # summarising metrics of interest (more districts for no.lists than for coverage)
  group_by(COVID, COUNTY.CODE) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(GROUP.ID)) %>% 
  ungroup() %>% 
  # only considering districts covered at least once in the four years
  inner_join(dist_cons_NL %>% distinct(COUNTY.CODE),
             by = "COUNTY.CODE") %>% 
  group_by(COUNTY.CODE) %>% 
  complete(COVID, 
           fill = list(NO.LISTS = 0)) %>% 
  ungroup() %>% 
  # since DUR is two years
  mutate(NO.LISTS = ifelse(COVID == "DUR", floor(NO.LISTS/2), NO.LISTS)) %>% 
  left_join(dist_cons_NL, by = "COUNTY.CODE")

# for grid coverage per district without threshold
temp_GC1 <- temp0_GC %>% 
  # average gridcov across two DUR years (since DUR is two years)
  group_by(M.YEAR, COUNTY.CODE, NO.CELLS.TOT) %>% 
  dplyr::summarise(NO.CELLS = n_distinct(CELL.ID),
                   GRID.COV = 100*NO.CELLS/min(NO.CELLS.TOT),
                   COVID = first(COVID)) %>% 
  group_by(COVID, COUNTY.CODE) %>% 
  dplyr::summarise(GRID.COV = mean(GRID.COV)) %>% 
  ungroup() %>% 
  # only considering districts covered at least once in the four years
  inner_join(dist_cons_GC %>% distinct(COUNTY.CODE),
             by = "COUNTY.CODE") %>% 
  group_by(COUNTY.CODE) %>% 
  complete(COVID, 
           fill = list(GRID.COV = 0)) %>% 
  ungroup() %>% 
  left_join(dist_cons_GC, by = "COUNTY.CODE")

# for grid coverage per district with threshold
temp_GC2 <- temp0_GC %>% 
  # average gridcov across two DUR years (since DUR is two years)
  group_by(M.YEAR, COUNTY.CODE, NO.CELLS.TOT, CELL.ID) %>% 
  # threshold of 5 lists in a grid cell (per year) for the cell to be considered "covered"
  dplyr::summarise(CELL.LISTS = n_distinct(GROUP.ID),
                   COVID = first(COVID)) %>% 
  filter(CELL.LISTS >= 5) %>% 
  # grid coverage with threshold
  dplyr::summarise(NO.CELLS = n_distinct(CELL.ID),
                   GRID.COV.T = 100*NO.CELLS/min(NO.CELLS.TOT),
                   COVID = first(COVID)) %>% 
  group_by(COVID, COUNTY.CODE) %>% 
  dplyr::summarise(GRID.COV.T = mean(GRID.COV.T)) %>% 
  ungroup() %>% 
  # only considering districts covered at least once in the four years
  inner_join(dist_cons_GC %>% distinct(COUNTY.CODE),
             by = "COUNTY.CODE") %>% 
  group_by(COUNTY.CODE) %>% 
  complete(COVID, 
           fill = list(GRID.COV = 0)) %>% 
  ungroup() %>% 
  left_join(dist_cons_GC, by = "COUNTY.CODE")


data_spread0 <- temp_NL %>% # NL has most districts
  left_join(temp_GC1, by = c("COUNTY.CODE", "COVID", "NO.CELLS.TOT")) %>% 
  left_join(temp_GC2, by = c("COUNTY.CODE", "COVID", "NO.CELLS.TOT")) %>% 
  pivot_longer(cols = c("NO.LISTS", "GRID.COV", "GRID.COV.T"),
               names_to = "METRIC", values_to = "VALUE") %>% 
  # some NAs produced for grid coverage because those districts don't qualify for coverage
  # analysis but present in number of lists analysis, so need to remove them
  filter(!(is.na(VALUE))) %>% 
  pivot_wider(names_from = "COVID", values_from = "VALUE") %>% 
  # calculating change (produces NA for x == 0 & y == 0 which need to be removed)
  # number of lists proportional change also has threshold of 10 lists (in function)
  mutate(RAW.CHANGE1 = s_spread_change(BEF, DUR),
         RAW.CHANGE2 = s_spread_change(DUR, AFT),
         PROP.CHANGE1a = s_spread_propchange(METRIC, BEF, DUR),
         PROP.CHANGE1b = s_spread_propchange(METRIC, DUR, AFT),
         # reverse direction
         PROP.CHANGE2a = s_spread_propchange(METRIC, DUR, BEF),
         PROP.CHANGE2b = s_spread_propchange(METRIC, AFT, DUR)) %>% 
  # for district geometry
  left_join(dists_ebd_sf %>% distinct(COUNTY.CODE, STATE, COUNTY, DISTRICT.GEOM),
            by = "COUNTY.CODE")


transition_names <- data.frame(T.CODE = c("T1", "T2", "T3", "T4"),
                               T.LABEL = c("Before to During", "During to After",
                                           "During to Before", "After to During"))


# Change in number of lists per district #

# 1. map of raw change
# 2. with SEs of raw change
# 3. map of proportional change

# 1. map (to give an idea of overall clustering)

data_spread1 <- data_spread0 %>% 
  filter(METRIC == "NO.LISTS") %>% 
  distinct(COUNTY.CODE, DISTRICT.GEOM, RAW.CHANGE1, RAW.CHANGE2) %>% 
  # pivoting longer for facetting
  pivot_longer(cols = c("RAW.CHANGE1", "RAW.CHANGE2"),
               names_to = "METRIC", values_to = "VALUE") %>% 
  mutate(T.CODE = case_when(str_detect(METRIC, "1") ~ "T1",
                            str_detect(METRIC, "2") ~ "T2")) %>% 
  left_join(transition_names, by = "T.CODE")

s_spread_map <- ggplot(data_spread1) +
  geom_sf(data = india_sf,
          fill = "#ACACAC", colour = "black", size = 0.4) +
  geom_sf(aes(fill = VALUE, geometry = DISTRICT.GEOM), col ="#ACACAC", size = 0.1) +
  facet_wrap(~ T.LABEL) +
  # scale_fill_continuous_divergingx(palette = "RdBu", mid = 0,
  #                                  n_interp = 7, na.value = "grey") +
  scale_fill_stepsn(colours = RColorBrewer::brewer.pal(n = 9, name = "RdBu"),
                    breaks = c(-1000, -500, -100, -50, 50, 100, 500, 1000),
                    values = scales::rescale(c(-11000, -1000, -500, -100, -50, 
                                               50, 100, 500, 1000, 11000)), 
                    limits = c(-11000, 11000),
                    na.value = "#ACACAC") +
  labs(x = "Longitude", y = "Latitude",
       fill = "Difference in\nno. of lists",
       title = "Change in no. of lists per district in transition periods",
       subtitle = "Grey: districts not considered in analysis")
# not considered can be either because both time periods are zero, or if generally no 
# district-grid cell link present

ggsave(filename = glue("01_birder_figs/{anal_name}_1effort_a.png"), 
       plot = s_spread_map,
       dpi = 300, width = 14, height = 8, units = "in")



# 2. with SEs

data_spread2 <- data_spread1 %>% 
  # removing NA cells (when both periods of the transition are zero)
  filter(!is.na(VALUE)) %>% 
  group_by(T.CODE, T.LABEL) %>% 
  # this gives B number of means of bootstrapped samples
  reframe(VALUE = boot_conf(x = VALUE)) %>% 
  group_by(T.CODE, T.LABEL) %>% 
  reframe(CI.L = stats::quantile(VALUE, 0.025), # Obtain the CIs
          CI.U = stats::quantile(VALUE, 0.975),
          VALUE = median(VALUE)) 

s_spread_rawchange <- ggplot(data_spread2, aes(x = T.LABEL, y = VALUE)) + 
  geom_point(size = 4, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                linewidth = 1.75, width = 0.08, position = position_dodge(0.5)) +
  # scale_y_continuous(limits = c(1.1, 2.3), breaks = seq(1, 3, 0.2)) +
  labs(title = "Net change in birding effort (no. of lists) across districts\nin transition periods",
       x = "Transition", y = "Difference in no. of lists") 

ggsave(filename = glue("01_birder_figs/{anal_name}_1effort_b.png"), 
       plot = s_spread_rawchange,
       dpi = 300, width = 6, height = 9, units = "in")

# is DUR-AFT really lower than BEF-DUR?
ggplot(data_spread1) +
  geom_histogram(aes(VALUE)) +
  facet_wrap(~ T.LABEL)


# 3. map

# forward direction
data_spread3a <- data_spread0 %>% 
  filter(METRIC == "NO.LISTS") %>% 
  distinct(COUNTY.CODE, DISTRICT.GEOM, PROP.CHANGE1a, PROP.CHANGE1b) %>% 
  # pivoting longer for facetting
  pivot_longer(cols = c("PROP.CHANGE1a", "PROP.CHANGE1b"),
               names_to = "METRIC", values_to = "VALUE") %>% 
  mutate(T.CODE = case_when(str_detect(METRIC, "a") ~ "T1",
                            str_detect(METRIC, "b") ~ "T2")) %>% 
  left_join(transition_names) %>% 
  # keeping only negative change
  mutate(VALUE = ifelse(VALUE < 0, VALUE, NA_real_))

# backward direction
data_spread3b <- data_spread0 %>% 
  filter(METRIC == "NO.LISTS") %>% 
  distinct(COUNTY.CODE, DISTRICT.GEOM, PROP.CHANGE2a, PROP.CHANGE2b) %>% 
  # pivoting longer for facetting
  pivot_longer(cols = c("PROP.CHANGE2a", "PROP.CHANGE2b"),
               names_to = "METRIC", values_to = "VALUE") %>% 
  mutate(T.CODE = case_when(str_detect(METRIC, "a") ~ "T3",
                            str_detect(METRIC, "b") ~ "T4")) %>% 
  left_join(transition_names) %>% 
  # keeping only negative change
  mutate(VALUE = ifelse(VALUE < 0, VALUE, NA_real_)) %>% 
  # order
  mutate(T.LABEL = factor(T.LABEL, levels = c("After to During", "During to Before")))


s_spread_mapprop <- (ggplot(data_spread3a) +
                       geom_sf(data = india_sf,
                               fill = "white", colour = "black", size = 0.4) +
                       geom_sf(aes(fill = VALUE, geometry = DISTRICT.GEOM), col ="#ACACAC", size = 0.1) +
                       facet_wrap(~ T.LABEL) +
                       # scale_fill_continuous_divergingx(palette = "RdBu", mid = 0,
                       #                                  n_interp = 7, na.value = "grey") +
                       scale_fill_steps2(low = viridisLite::cividis(5, direction = -1)[1],
                                         mid = viridisLite::cividis(5, direction = -1)[3],
                                         high = viridisLite::cividis(5, direction = -1)[5],
                                         midpoint = -0.45,
                                         breaks = c(-1, -0.9, -0.65, -0.45, -0.25, -0.1, -0),
                                         limits = c(-1, 0),
                                         na.value = "#ACACAC") +
                       labs(x = "Longitude", y = "Latitude",
                            fill = "Proportional change")) /
  (ggplot(data_spread3b) +
     geom_sf(data = india_sf,
             fill = "white", colour = "black", size = 0.4) +
     geom_sf(aes(fill = VALUE, geometry = DISTRICT.GEOM), col ="#ACACAC", size = 0.1) +
     facet_wrap(~ T.LABEL) +
     # scale_fill_viridis_b(option = "cividis", direction = -1, na.value = "grey",
     #                      values = c(0.1, 0.25, 0.45, 0.65, 0.9)) +
     scale_fill_steps2(low = viridisLite::cividis(5, direction = -1)[1],
                       mid = viridisLite::cividis(5, direction = -1)[3],
                       high = viridisLite::cividis(5, direction = -1)[5],
                       midpoint = -0.45,
                       breaks = c(-1, -0.9, -0.65, -0.45, -0.25, -0.1, -0),
                       limits = c(-1, 0),
                       na.value = "#ACACAC") +
     labs(x = "Longitude", y = "Latitude",
          fill = "Proportional change")) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Proportional declines in no. of lists per district in transition periods (only when N1 >= 10)",
                  subtitle = "Grey: districts with positive proportional change\nWhite: districts not considered in analysis")

ggsave(filename = glue("01_birder_figs/{anal_name}_1effort_c.png"), 
       plot = s_spread_mapprop,
       dpi = 300, width = 14, height = 15, units = "in")


# Change in grid coverage per district (without & with threshold) #

# 1. map of raw change
# 2. with SEs of raw change
# 3. map of proportional change

coverage_types <- data.frame(METRIC = c("GRID.COV", "GRID.COV.T"),
                             METRIC.LABEL = c("Coverage (%)", "Coverage with threshold (%)"))

# 1. map (to give an idea of overall clustering)

data_spread4 <- data_spread0 %>% 
  filter(METRIC == "GRID.COV" | METRIC == "GRID.COV.T") %>% 
  distinct(METRIC, COUNTY.CODE, DISTRICT.GEOM, RAW.CHANGE1, RAW.CHANGE2) %>% 
  # pivoting longer for facetting
  pivot_longer(cols = c("RAW.CHANGE1", "RAW.CHANGE2"),
               names_to = "CHANGE.TYPE", values_to = "VALUE") %>% 
  mutate(T.CODE = case_when(str_detect(CHANGE.TYPE, "1") ~ "T1",
                            str_detect(CHANGE.TYPE, "2") ~ "T2")) %>% 
  left_join(transition_names) %>% 
  left_join(coverage_types)

s_spread_gridcov_map <- ggplot(data_spread4) +
  geom_sf(data = india_sf,
          fill = "#ACACAC", colour = "black", size = 0.4) +
  geom_sf(aes(fill = VALUE, geometry = DISTRICT.GEOM), col ="#ACACAC", size = 0.1) +
  facet_grid(METRIC.LABEL ~ T.LABEL, switch = "y") +
  scale_fill_stepsn(colours = RColorBrewer::brewer.pal(n = 9, name = "RdBu"),
                    breaks = c(-40, -25, -10, -1, 1, 10, 25, 40),
                    values = scales::rescale(c(-100, -40, -25, -10, -1, 1, 10, 25, 40, 100)), 
                    limits = c(-100, 100),
                    na.value = "#ACACAC") +
  labs(x = "Longitude", y = "Latitude",
       fill = "Difference in\ngrid coverage",
       title = "Change in grid coverage per district in transition periods",
       subtitle = "Grey: districts not considered in analysis")

ggsave(filename = glue("01_birder_figs/{anal_name}_2gridcov_a.png"), 
       plot = s_spread_gridcov_map,
       dpi = 300, width = 14, height = 14, units = "in")


# 2. with SEs

data_spread5 <- data_spread4 %>% 
  # removing NA cells (when both periods of the transition are zero)
  filter(!is.na(VALUE)) %>% 
  group_by(METRIC, METRIC.LABEL, T.CODE, T.LABEL) %>% 
  # this gives B number of means of bootstrapped samples
  reframe(VALUE = boot_conf(x = VALUE)) %>% 
  group_by(METRIC, METRIC.LABEL, T.CODE, T.LABEL) %>% 
  reframe(CI.L = stats::quantile(VALUE, 0.025), # Obtain the CIs
          CI.U = stats::quantile(VALUE, 0.975),
          VALUE = median(VALUE)) 

s_spread_gridcov_rawchange <- ggplot(data_spread5, aes(x = T.LABEL, y = VALUE)) + 
  facet_wrap(~ METRIC.LABEL, nrow = 2) +
  geom_point(size = 4, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                size = 1.75, width = 0.08, position = position_dodge(0.5)) +
  # scale_y_continuous(limits = c(1.1, 2.3), breaks = seq(1, 3, 0.2)) +
  labs(title = "Net change in grid coverage across districts\nin transition periods",
       x = "Transition", y = "Difference in grid coverage") 

ggsave(filename = glue("01_birder_figs/{anal_name}_2gridcov_b.png"), 
       plot = s_spread_gridcov_rawchange,
       dpi = 300, width = 6, height = 9, units = "in")


# 3. map

data_spread6 <- data_spread0 %>% 
  filter(METRIC == "GRID.COV" | METRIC == "GRID.COV.T") %>% 
  distinct(METRIC, COUNTY.CODE, DISTRICT.GEOM, 
           PROP.CHANGE1a, PROP.CHANGE1b, PROP.CHANGE2a, PROP.CHANGE2b) %>% 
  # pivoting longer for facetting
  pivot_longer(cols = c("PROP.CHANGE1a", "PROP.CHANGE1b",
                        "PROP.CHANGE2a", "PROP.CHANGE2b"),
               names_to = "CHANGE.TYPE", values_to = "VALUE") %>% 
  mutate(T.CODE = case_when(str_detect(CHANGE.TYPE, "1a") ~ "T1",
                            str_detect(CHANGE.TYPE, "1b") ~ "T2",
                            str_detect(CHANGE.TYPE, "2a") ~ "T3",
                            str_detect(CHANGE.TYPE, "2b") ~ "T4")) %>% 
  left_join(transition_names) %>% 
  left_join(coverage_types) %>% 
  # keeping only negative change
  mutate(VALUE = ifelse(VALUE < 0, VALUE, NA_real_)) %>% 
  mutate(T.LABEL = factor(T.LABEL, levels = c("Before to During", "During to After",
                                              "After to During", "During to Before")))

# with the threshold, the IQR of proportional change has increased even though range is smaller
# i.e., the distribution is much more rounded
# this results in more districts with higher declines, but what's important here is to
# compare bef-dur with aft-dur and see which districts overlap---not many
# then compare these two with dur-aft and dur-bef respectively to check whether changes are 
# only in one direction---evidently not

s_spread_gridcov_mapprop_a <- ggplot(data_spread6 %>% 
                                       filter(METRIC == "GRID.COV")) +
  geom_sf(data = india_sf,
          fill = "white", colour = "black", size = 0.4) +
  geom_sf(aes(fill = VALUE, geometry = DISTRICT.GEOM), col ="#ACACAC", size = 0.1) +
  facet_wrap(~ T.LABEL, ncol = 2) +
  scale_fill_steps2(low = viridisLite::cividis(5, direction = -1)[1],
                    mid = viridisLite::cividis(5, direction = -1)[3],
                    high = viridisLite::cividis(5, direction = -1)[5],
                    midpoint = -0.45,
                    breaks = c(-1, -0.9, -0.65, -0.45, -0.25, -0.1, -0),
                    limits = c(-1, 0),
                    na.value = "#ACACAC") +
  labs(x = "Longitude", y = "Latitude",
       fill = "Proportional change",
       title = "Proportional declines in grid coverage per district in transition periods",
       subtitle = "Grey: districts with positive proportional change\nWhite: districts not considered in analysis")

ggsave(filename = glue("01_birder_figs/{anal_name}_2gridcov_c1.png"), 
       plot = s_spread_gridcov_mapprop_a,
       dpi = 300, width = 14, height = 14, units = "in")

s_spread_gridcov_mapprop_b <- ggplot(data_spread6 %>% 
                                       filter(METRIC == "GRID.COV.T")) +
  geom_sf(data = india_sf,
          fill = "white", colour = "black", size = 0.4) +
  geom_sf(aes(fill = VALUE, geometry = DISTRICT.GEOM), col ="#ACACAC", size = 0.1) +
  facet_wrap(~ T.LABEL, ncol = 2) +
  scale_fill_steps2(low = viridisLite::cividis(5, direction = -1)[1],
                    mid = viridisLite::cividis(5, direction = -1)[3],
                    high = viridisLite::cividis(5, direction = -1)[5],
                    midpoint = -0.45,
                    breaks = c(-1, -0.9, -0.65, -0.45, -0.25, -0.1, -0),
                    limits = c(-1, 0),
                    na.value = "#ACACAC") +
  labs(x = "Longitude", y = "Latitude",
       fill = "Proportional change",
       title = "Proportional declines in grid coverage per district (with threshold) in transition periods",
       subtitle = "Grey: districts with positive proportional change\nWhite: districts not considered in analysis")

ggsave(filename = glue("01_birder_figs/{anal_name}_2gridcov_c2.png"), 
       plot = s_spread_gridcov_mapprop_b,
       dpi = 300, width = 14, height = 14, units = "in")

# all these values are zero in one of the two transition periods.

# Saving analysis objects #

map_dists_ebd_sf <- dists_ebd_sf %>% 
  distinct(COUNTY.CODE, STATE.NAME, DISTRICT.NAME)

temp1 <- data_spread1 %>% st_drop_geometry() %>% dplyr::select(-DISTRICT.GEOM)
temp2 <- data_spread2
temp3 <- data_spread3a %>% st_drop_geometry() %>% dplyr::select(-DISTRICT.GEOM)
temp4 <- data_spread3b %>% st_drop_geometry() %>% dplyr::select(-DISTRICT.GEOM)
temp5 <- data_spread4 %>% st_drop_geometry() %>% dplyr::select(-DISTRICT.GEOM)
temp6 <- data_spread5 
temp7 <- data_spread6 %>% st_drop_geometry() %>% dplyr::select(-DISTRICT.GEOM)


save(map_dists_ebd_sf, temp1, temp2, temp3, temp4, temp5, temp6, temp7,
     file = glue("00_outputs/{anal_name}.RData"))


# Temporal spread ---------------------------------------------------------


anal_name <- "d10a_t_dow"

# Change in DoW spread #

met_week <- function(dates) {
  require(lubridate)
  
  normalyear <- c((0:363 %/% 7 + 1), 52)
  leapyear   <- c(normalyear[1:59], 9, normalyear[60:365])
  yearday    <- yday(dates)
  
  return(ifelse(leap_year(dates), leapyear[yearday], normalyear[yearday])) 
}


data_a <- data0_MY_d_slice_G %>% 
  mutate(WEEK.Y = met_week(OBSERVATION.DATE)) %>% 
  mutate(DAY.W = wday(OBSERVATION.DATE, 
                      # start week from Monday
                      week_start = getOption("lubridate.week.start", 1)),
         DAY.W.LABEL = wday(OBSERVATION.DATE, 
                            week_start = getOption("lubridate.week.start", 1),
                            label = T)) %>% 
  group_by(M.YEAR, WEEK.Y) %>% 
  mutate(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  group_by(M.YEAR, WEEK.Y, DAY.W, DAY.W.LABEL) %>% 
  dplyr::summarise(TOT.LISTS = min(TOT.LISTS),
                   DAY.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER)),
                   PROP.LISTS = DAY.LISTS/TOT.LISTS) %>% 
  group_by(M.YEAR, DAY.W, DAY.W.LABEL) %>% 
  dplyr::summarise(PROP.LISTS = boot_conf(PROP.LISTS)) %>% 
  group_by(M.YEAR, DAY.W, DAY.W.LABEL) %>% 
  dplyr::summarise(CI.L = stats::quantile(PROP.LISTS, 0.025), # Obtain the CIs
                   CI.U = stats::quantile(PROP.LISTS, 0.975),
                   PROP.LISTS = median(PROP.LISTS)) 

data_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  mutate(STATE = factor(STATE, levels = anal_states[1,]),
         WEEK.Y = met_week(OBSERVATION.DATE)) %>% 
  mutate(DAY.W = wday(OBSERVATION.DATE, 
                      # start week from Monday
                      week_start = getOption("lubridate.week.start", 1)),
         DAY.W.LABEL = wday(OBSERVATION.DATE, 
                            week_start = getOption("lubridate.week.start", 1),
                            label = T)) %>% 
  group_by(M.YEAR, STATE, WEEK.Y) %>% 
  mutate(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  group_by(M.YEAR, STATE, WEEK.Y, DAY.W, DAY.W.LABEL) %>% 
  dplyr::summarise(TOT.LISTS = min(TOT.LISTS),
                   DAY.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER)),
                   PROP.LISTS = DAY.LISTS/TOT.LISTS) %>% 
  group_by(M.YEAR, STATE, DAY.W, DAY.W.LABEL) %>% 
  dplyr::summarise(PROP.LISTS = boot_conf(PROP.LISTS)) %>% 
  group_by(M.YEAR, STATE, DAY.W, DAY.W.LABEL) %>% 
  dplyr::summarise(CI.L = stats::quantile(PROP.LISTS, 0.025), # Obtain the CIs
                   CI.U = stats::quantile(PROP.LISTS, 0.975),
                   PROP.LISTS = median(PROP.LISTS)) 


plot_a <- ggplot(data_a, aes(DAY.W.LABEL, PROP.LISTS, colour = M.YEAR)) + 
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                size = 1.25, width = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear") +
  labs(title = "National level",
       x = "Day of week", y = "Proportion of lists")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 11, height = 6, units = "in")


plot_b <- ggplot(data_b, aes(DAY.W.LABEL, PROP.LISTS, colour = M.YEAR)) + 
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  geom_point(size = 3.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                size = 1.5, width = 0.5, position = position_dodge(0.5)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear") +
  labs(title = "State-level",
       x = "Day of week", y = "Proportion of lists")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 22, height = 13, units = "in")

# Saving analysis objects #

save(data_a, data_b, 
     file = glue("00_outputs/{anal_name}.RData"))





anal_name <- "d10b_t_tod"

# Change in ToD spread #

data_a <- data0_MY_d_slice_G %>% 
  mutate(DAY.Y = yday(OBSERVATION.DATE)) %>% 
  group_by(M.YEAR, DAY.Y) %>% 
  mutate(TOT.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER))) %>%
  group_by(M.YEAR, DAY.Y, HOUR) %>% 
  dplyr::summarise(TOT.LISTS = min(TOT.LISTS),
                   HOUR.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER)),
                   PROP.LISTS = HOUR.LISTS/TOT.LISTS) %>% 
  group_by(M.YEAR) %>% 
  complete(HOUR = 0:23, fill = list(PROP.LISTS = 0)) %>% 
  group_by(M.YEAR, HOUR) %>% 
  dplyr::summarise(PROP.LISTS = boot_conf(PROP.LISTS)) %>% 
  group_by(M.YEAR, HOUR) %>% 
  dplyr::summarise(CI.L = stats::quantile(PROP.LISTS, 0.025), # Obtain the CIs
                   CI.U = stats::quantile(PROP.LISTS, 0.975),
                   PROP.LISTS = median(PROP.LISTS)) 

data_b <- data0_MY_d_slice_G %>% 
  filter(STATE %in% anal_states[1,]) %>% 
  mutate(STATE = factor(STATE, levels = anal_states[1,]),
         DAY.Y = yday(OBSERVATION.DATE)) %>% 
  group_by(M.YEAR, STATE, DAY.Y) %>% 
  mutate(TOT.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER))) %>%
  group_by(M.YEAR, STATE, DAY.Y, HOUR) %>% 
  dplyr::summarise(TOT.LISTS = min(TOT.LISTS),
                   HOUR.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER)),
                   PROP.LISTS = HOUR.LISTS/TOT.LISTS) %>% 
  group_by(M.YEAR, STATE) %>% 
  complete(HOUR = 0:23, fill = list(PROP.LISTS = 0)) %>% 
  group_by(M.YEAR, STATE, HOUR) %>% 
  dplyr::summarise(PROP.LISTS = boot_conf(PROP.LISTS)) %>% 
  group_by(M.YEAR, STATE, HOUR) %>% 
  dplyr::summarise(CI.L = stats::quantile(PROP.LISTS, 0.025), # Obtain the CIs
                   CI.U = stats::quantile(PROP.LISTS, 0.975),
                   PROP.LISTS = median(PROP.LISTS)) 


plot_a <- ggplot(data_a, aes(HOUR, PROP.LISTS, colour = M.YEAR)) + 
  geom_point(size = 2, position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                size = 1.2, width = 0.6, position = position_dodge(0.8)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear") +
  scale_x_continuous(limits = c(-1, 24), breaks = 0:23) +
  labs(title = "National level",
       x = "Time of day (hours)", y = "Proportion of lists")

ggsave(filename = glue("01_birder_figs/{anal_name}_a.png"), plot = plot_a,
       dpi = 300, width = 16, height = 8, units = "in")


plot_b <- ggplot(data_b, aes(HOUR, PROP.LISTS, colour = M.YEAR)) + 
  facet_wrap(~ STATE, ncol = 2, scales = "free_y") +
  geom_point(size = 2.5, position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                size = 1.4, width = 0.6, position = position_dodge(0.8)) +
  scale_colour_manual(values = covid_palette, name = "Migratory\nyear") +
  scale_x_continuous(limits = c(-1, 24), breaks = 0:23) +
  labs(title = "State-level",
       x = "Time of day (hours)", y = "Proportion of lists")

ggsave(filename = glue("01_birder_figs/{anal_name}_b.png"), plot = plot_b,
       dpi = 300, width = 32, height = 17, units = "in")


# Saving analysis objects #

save(data_a, data_b, 
     file = glue("00_outputs/{anal_name}.RData"))



