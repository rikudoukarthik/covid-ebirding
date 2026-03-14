# data objects required for summaries, tables and figures in manuscript

# first run all setup and data import steps in 03_wrap.Rmd, 04_birder.Rmd and 05_bird.Rmd



# data quantity (supp.) ---------------------------------------------------

dataquant_a <- data0_MY_d_slice_S %>% 
  group_by(M.YEAR, MONTH) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER),
                   NO.EBIRDERS = n_distinct(OBSERVER.ID)) %>%
  ungroup()

dataquant_b <- data0_MY_d_slice_S %>% 
  filter(STATE %in% anal_states) %>% 
  group_by(STATE, M.YEAR, MONTH) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER),
                   NO.EBIRDERS = n_distinct(OBSERVER.ID)) %>%
  ungroup()


# months timeline ---------------------------------------------------------

months_timeline <- data0_MY_b_slice_S %>% 
  arrange(OBSERVATION.DATE) %>% 
  distinct(M.YEAR, MONTH) 

# total no. of unique lists -----------------------------------------------

tot_uniq_lists <- n_distinct(data0_MY_b_slice_S$GROUP.ID)

# b_slice_S because b data has more unique than d data (different groupacc filters)


# prolific qualifications -------------------------------------------------

# locations, species, lists that qualified for prolific models

tot_qual_prolific <- data.frame(
  LOCATIONS = n_distinct(data_prolific$LOCALITY.ID),
  SPECIES = n_distinct(data_prolific$COMMON.NAME),
  NU.LISTS = n_distinct(data_prolific$SAMPLING.EVENT.IDENTIFIER)
)


# saving ------------------------------------------------------------------

save(covid_palette, # timeline graph 
     totalcells, # how many 25 x 25 cells
     dataquant_a, dataquant_b, # data quantity
     months_timeline,
     tot_uniq_lists, # how many total unique lists
     species_list, # how many species (and urban, rural)
     fail_spec_KL, fail_spec_MH, # for how many spp in overall bird models used logit
     tot_qual_prolific,
     file = "00_outputs/data_for_ms.RData")
