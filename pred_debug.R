# model_a <- glmer(URBAN ~ M.YEAR + M.YEAR:MONTH + (1|COUNTY/CELL.ID),
#                  data = data_a, family = binomial(link = "cloglog"),
#                  nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 


x <- data0_MY_d_slice_G %>% 
  distinct(M.YEAR, MONTH, COUNTY, CELL.ID, SAMPLING.EVENT.IDENTIFIER, URBAN) %>%
  group_by(M.YEAR, MONTH, COUNTY, CELL.ID) %>% 
  mutate(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  filter(URBAN == 1) %>% 
  summarise(U.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER),
            PROP = U.LISTS/min(TOT.LISTS)) %>% 
  summarise(PROP = mean(PROP)) %>% 
  summarise(PROP = mean(PROP))

hist(x$PROP)

n_distinct((data0_MY_d_slice_G %>% filter(URBAN == 1))$SAMPLING.EVENT.IDENTIFIER) /
  n_distinct(data0_MY_d_slice_G$SAMPLING.EVENT.IDENTIFIER)




# model_a <- glmer(HOTSPOT ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
#                  data = data_a, family = binomial(link = "cloglog"),
#                  nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 

y <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, CELL.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(HOTSPOT = if_else(LOCALITY.TYPE == "H", 1, 0)) %>%
  mutate(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  filter(HOTSPOT == 1) %>% 
  summarise(H.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER),
            PROP = H.LISTS/min(TOT.LISTS)) %>% 
  summarise(PROP = mean(PROP))

hist(y$PROP)



# model_a <- glmer(TRAVELLING ~ M.YEAR + M.YEAR:MONTH + (1|CELL.ID),
#                  data = data_a, family = binomial(link = "cloglog"),
#                  nAGQ = 0, control = glmerControl(optimizer = "bobyqa")) 

z <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, CELL.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(TRAVELLING = if_else((PROTOCOL.TYPE == "Traveling" & 
                                           EFFORT.DISTANCE.KM > 0.3), 1, 0)) %>%
  mutate(TOT.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  filter(TRAVELLING == 1) %>% 
  summarise(T.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER),
            PROP = T.LISTS/min(TOT.LISTS)) %>% 
  summarise(PROP = mean(PROP))

hist(z$PROP)



anal_name <- "d04_hotspot"
anal_name <- "d05_protocol"
anal_name <- "d09b_s_cover"
anal_name <- "d09a_s_UNU"
anal_name <- "d01_group_po"
anal_name <- "d03_time_po"


load(glue("outputs/{anal_name}.RData"))

summary(model_a)






a <- data0_MY_d_slice_S %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID, SAMPLING.EVENT.IDENTIFIER) %>% 
  dplyr::summarise(GROUP.BIRDING = if_else(NUMBER.OBSERVERS > 2, 1, 0)) %>%
  mutate(TOT = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  filter(GROUP.BIRDING == 1) %>% 
  summarise(GROUP.BIRDING = n_distinct(SAMPLING.EVENT.IDENTIFIER),
            PROP = GROUP.BIRDING/min(TOT)) %>% 
  summarise(PROP = mean(PROP))

b <- data0_MY_d_slice_S %>% 
  group_by(M.YEAR, MONTH, OBSERVER.ID) %>% 
  dplyr::summarise(NO.SITES = log(n_distinct(CELL.ID))) %>% 
  dplyr::summarise(NO.SITES = exp(mean(NO.SITES))) 

c <- data0_MY_d_slice_G %>% 
  group_by(M.YEAR, MONTH, CELL.ID) %>% 
  dplyr::summarise(COVERED = if_else(n_distinct(SAMPLING.EVENT.IDENTIFIER) >= 5, 1, 0)) %>% 
  # threshold to consider cell "covered"
  group_by(M.YEAR, MONTH) %>% 
  complete(CELL.ID = gridmapg1_IN$CELL.ID, fill = list(COVERED = 0)) %>% 
  mutate(TOT = n_distinct(CELL.ID)) %>% 
  filter(COVERED == 1) %>% 
  summarise(COVERED = n_distinct(CELL.ID),
            PROP = COVERED/min(TOT)) 
