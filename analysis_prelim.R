############ covid-ebirding


# importing from raw ebd ###
#
# require(tidyverse)
# require(dtplyr)
# require(data.table)
# library(lubridate)
# 
# memory.limit(size = 20000)
# load("ebd_IN_relAug-2021.RData") # 8-11 minutes
# ### ad hoc filters to produce tibble from dtplyr
# data <- data %>% filter(YEAR >= 2017)
# data <- collect(data)
# rm("data_DT")
# gc()
# save.image("data.RData") 
# 
# detach("package:dtplyr", unload = T)
# detach("package:data.table", unload = T)
# 
# ###

library(tidyverse)
library(lubridate)
memory.limit(size = 20000)
load("data.RData")

data <- data %>% filter(PROTOCOL.TYPE %in% c("Traveling","Stationary")) %>%
  filter(DURATION.MINUTES < 1500) %>% 
  filter(!(OBSERVER.ID %in% c("obsr852457","obsr949737"))) # Salem group acc.s

data_obs <- data %>% filter(PROTOCOL.TYPE %in% c("Traveling","Stationary")) %>% 
  group_by(YEAR,OBSERVER.ID,OBSERVATION.DATE,SAMPLING.EVENT.IDENTIFIER) %>% slice(1) %>% 
  ungroup() %>% group_by(YEAR,OBSERVER.ID,OBSERVATION.DATE) %>%
  summarise(BIRDING.TIME = sum(DURATION.MINUTES), NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>% 
  ungroup()

eBird_users <- read.delim("ebd_users_relJun-2021.txt", sep = "\t", header = T, quote = "", 
                          stringsAsFactors = F, na.strings = c(""," ",NA))
names(eBird_users) <- c("OBSERVER.ID","FIRST.NAME","LAST.NAME")
eBird_users <- eBird_users %>% transmute(OBSERVER.ID = OBSERVER.ID,
                                         FULL.NAME = paste(FIRST.NAME, LAST.NAME))


data_x <- data_obs %>% filter(NO.LISTS > 50) %>% select(OBSERVER.ID) %>% distinct()
x <- left_join(data_x, eBird.users)


ggplot(data_obs) + theme_bw() +
  scale_y_log10() +
  geom_violin(aes(x = as.factor(YEAR), y = BIRDING.TIME))

ggplot(data_obs) + theme_bw() +
  geom_violin(aes(x = as.factor(YEAR), y = NO.LISTS))




data_sf <- data %>% filter(PROTOCOL.TYPE %in% c("Traveling","Stationary")) %>% 
  group_by(YEAR,OBSERVER.ID,OBSERVATION.DATE,SAMPLING.EVENT.IDENTIFIER) %>% slice(1) %>% 
  ungroup() %>% group_by(YEAR,OBSERVER.ID,MONTH) %>%
  summarise(NO.SITES = n_distinct(LOCALITY.ID)) %>% # monthly no. of sites
  ungroup()

ggplot(data_sf) + theme_bw() +
  scale_y_log10() +
  geom_boxplot(aes(x = as.factor(YEAR), y = NO.SITES))



ggplot(data_obs) + theme_bw() +
  scale_y_log10() +
  geom_violin(aes(x = as.factor(YEAR), y = DURATION.MINUTES))




