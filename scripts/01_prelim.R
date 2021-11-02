library(tidyverse)
library(lubridate)


# ### importing from raw ebd
# load("data/ebd_IN_relSep-2021.RData") # 2 mins!
# data <- data %>% filter(YEAR >= 2018) # 20 secs!
# save.image("data/rawdata.RData") # 3 mins
# # load("data/rawdata.RData")
# 
# ### filters
# groupaccs <- read.csv("data/ebd_users_GA_relSep-2021.csv", na.strings = c(""," ",NA),
#                         quote = "", header = T, nrows = 401)  # excluding empty cells
# 
# groupaccs <- groupaccs %>% 
#   mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", GA.2 == 1 ~ "GA.2", TRUE ~ "NG"))
# 
# filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)
# 
# data0 <- data %>% anti_join(filtGA) %>% filter(PROTOCOL.TYPE %in% c("Traveling","Stationary"))
# 
# rm(list = setdiff(ls(envir = .GlobalEnv), c("data0")), pos = ".GlobalEnv")
# save.image("data/data0.RData")
# 
# # random subset for analysis + adding COVID category
# set.seed(10)
# data_sub <- data0 %>% filter(GROUP.ID %in% sample(unique(data0$GROUP.ID),100000)) %>% 
#   mutate(COVID = case_when(YEAR %in% c(2017,2018,2019) ~ "BEF", TRUE ~ "DUR"))
# save(data_sub, file = "data/data_sub.RData")



data_obs <- data0 %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% 
  slice(1) %>% ungroup() %>% group_by(YEAR,OBSERVER.ID,OBSERVATION.DATE) %>%
  summarise(BIRDING.TIME = sum(DURATION.MINUTES), NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  ungroup %>% mutate(COVID = case_when(YEAR == 2017:2019 ~ "BEF", TRUE ~ "DUR"))


ggplot(data_obs) + theme_bw() +
  scale_y_log10() +
  geom_violin(aes(x = as.factor(YEAR), y = BIRDING.TIME))

ggplot(data_obs) + theme_bw() +
  scale_y_log10() +
  geom_violin(aes(x = as.factor(YEAR), y = NO.LISTS))


data_sf <- data0 %>% 
  group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1) %>% 
  ungroup() %>% mutate(LATITUDE = round(LATITUDE, 2), LONGITUDE = round(LONGITUDE, 2)) %>% #1kmx1km
  group_by(YEAR,OBSERVER.ID,MONTH) %>%
  summarise(NO.SITES = n_distinct(c(LATITUDE, LONGITUDE))) %>% # monthly no. of sites
  ungroup() %>% mutate(COVID = case_when(YEAR == 2017:2019 ~ "BEF", TRUE ~ "DUR"))

ggplot(data_sf) + theme_bw() +
  scale_y_log10() +
  geom_boxplot(aes(x = COVID, y = NO.SITES))



rm(list = setdiff(ls(envir = .GlobalEnv), c("data0","data_obs","data_sf")), pos = ".GlobalEnv")


library(lme4)

m1 <- glmer(BIRDING.TIME ~ COVID + (1|OBSERVER.ID) + (1|OBSERVATION.DATE),
            control = glmerControl(optimizer = "bobyqa"),
            data = data_obs, family=poisson()) # very large eigenvalue
m2 <- glmer(NO.LISTS ~ COVID + (1|OBSERVER.ID) + (1|OBSERVATION.DATE),
            data = data_obs, family=poisson())
m3 <- glmer(NO.SITES ~ COVID + (1|OBSERVER.ID) + (1|MONTH), 
            data = data_sf, family=poisson())



data_cl <- data0 %>% 
  group_by(GROUP.ID) %>% slice(1) %>% ungroup() %>% 
  mutate(DURATION.MINUTES = DURATION.MINUTES,
         PROTOCOL.TYPE = PROTOCOL.TYPE,
         # MEDIA = MEDIA
         COVID = case_when(YEAR %in% c(2017,2018,2019) ~ "BEF", TRUE ~ "DUR"),
         SHARED = case_when(NUMBER.OBSERVERS == 1 ~ 0, TRUE ~ 1),
         RARITY = REVIEWED)



save.image("temp.RData")
load("temp.RData")

