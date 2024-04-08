# Modelling directly with presence-absence data instead of relative abundance.
#
# parallelisation is somehow much less efficient when within a function or loop.

tictoc::tic("Total time elapsed for b_overall_model():")

cur_species_list <- species_list %>% filter(STATE == state_name)

# creating subfolder for model outputs
if (!dir.exists(glue("00_data/bird_models/{state_name}/"))) {
  dir.create(glue("00_data/bird_models/{state_name}/"), recursive = TRUE)
}


###

source("00_scripts/ss01_groupids.R")
source("00_scripts/ss02_datafiles.R")


# month type 1 ------------------------------------------------------------

cur_m <- "LD"

source("00_scripts/ss03_runspeciesmodels.R")



for (cur_sp in 1:n_distinct(birds_pred0$COMMON.NAME)) {

  birds_pred0_b <- birds_pred0 %>%
    filter(COMMON.NAME == unique(birds_pred0$COMMON.NAME)[cur_sp]) %>%
    rename(PRED.LINK2 = PRED.LINK,
           SE.LINK2 = SE.LINK)

  assign("birds_pred0_b", birds_pred0_b, envir = .GlobalEnv)


  data_spec <- data_mtype %>%
    filter(COMMON.NAME == unique(birds_pred0$COMMON.NAME)[cur_sp]) %>%
    # using only CELL.ID-MONTH (space-time) combos in which species occurs
    filter(REPORT == 1) %>%
    distinct(COMMON.NAME, CELL.ID, MONTH) %>%
    # joining both presences and absences for space-time combos of interest
    left_join(data_occ)


  tictoc::tic(glue("GLMM for months type {cur_m}, {unique(birds_pred0$COMMON.NAME)[cur_sp]}"))

  # for some species in KL and MH, using cloglog link is resulting in "PIRLS step-halvings
  # failed to reduce deviance in pwrssUpdate"

  if ((state_name == "Kerala" &
       unique(birds_pred0$COMMON.NAME)[cur_sp] %in% fail_spec_KL) |
      (state_name == "Maharashtra" &
       unique(birds_pred0$COMMON.NAME)[cur_sp] %in% fail_spec_MH)){
    model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:NO.SP + MONTH:M.YEAR +
                          (1|CELL.ID),
                        data = data_spec, family = binomial,
                        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
  } else {
    model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:NO.SP + MONTH:M.YEAR +
                          (1|CELL.ID),
                        data = data_spec, family = binomial(link = "cloglog"),
                        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
  }

  tictoc::toc()


  tictoc::tic(glue("Bootstrapped predictions for months type {cur_m}, {unique(birds_pred0$COMMON.NAME)[cur_sp]}"))
  prediction <- split_par_boot(model = model_spec,
                               new_data = birds_pred0_b,
                               new_data_string = "birds_pred0_b",
                               mode = "normal_low")
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
    # coalesce takes first non-NA so retains NA for non-current species
    mutate(PRED.LINK = coalesce(PRED.LINK, PRED.LINK2),
           SE.LINK = coalesce(SE.LINK, SE.LINK2)) %>%
    dplyr::select(-PRED.LINK2, -SE.LINK2)

  assign("birds_pred0", birds_pred0, envir = .GlobalEnv)

}

birds_pred <- birds_pred %>%
  left_join(birds_pred0, by = c("MONTHS.TYPE", "COMMON.NAME", "MONTH", "M.YEAR", "STATE",
                                "SP.CATEGORY", "NO.SP")) %>%
  mutate(PRED.LINK = coalesce(PRED.LINK.x, PRED.LINK.y),
         SE.LINK = coalesce(SE.LINK.x, SE.LINK.y)) %>%
  dplyr::select(-PRED.LINK.x, -PRED.LINK.y, -SE.LINK.x, -SE.LINK.y)


# month type 2 ------------------------------------------------------------

cur_m <- "NL"

source("00_scripts/ss03_runspeciesmodels.R")



# backtransform, summarise, CIs -----------------------------------------------------

birds_pred <- birds_pred %>% 
  mutate(PRED = clogloglink(PRED.LINK, inverse = T),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = clogloglink((PRED.LINK - SE.LINK), inverse = T)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  left_join(timeline) %>% 
  # summarising for species categories
  group_by(STATE, MONTHS.TYPE, M.YEAR, SP.CATEGORY) %>% 
  dplyr::summarise(PRED = mean(PRED),
                   # propagating SE across species of a category
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PRED - 1.96*SE,
                   CI.U = PRED + 1.96*SE)

tictoc::toc()

