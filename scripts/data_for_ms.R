# data objects required for summaries, tables and figures in manuscript

# first run all setup and data import steps in 03_wrap.Rmd, 04_birder.Rmd and 05_bird.Rmd


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
     tot_uniq_lists, # how many total unique lists
     species_list, # how many species (and urban, rural)
     fail_spec_KL, fail_spec_MH, # for how many spp in overall bird models used logit
     tot_qual_prolific,
     file = "outputs/data_for_ms.RData")