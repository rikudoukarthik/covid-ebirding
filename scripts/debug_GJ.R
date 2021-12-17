library(tidyverse)
library(lubridate)
library(patchwork)
theme_set(theme_bw())

#### load EBD Sep release RData

load("data/ebd_IN_relSep-2021.RData") 

### list of group accounts to be filtered

groupaccs <- read.csv("data/ebd_users_GA_relSep-2021.csv", 
                      na.strings = c(""," ",NA), quote = "", header = T, 
                      nrows = 401)  # excluding empty cells
groupaccs <- groupaccs %>% 
  mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                              GA.2 == 1 ~ "GA.2", 
                              TRUE ~ "NG"))
filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)



### new eBirder data

new_obsr_data <- data %>% # name of data object from RData file
  select(c("YEAR", "MONTH", "STATE", "SAMPLING.EVENT.IDENTIFIER",
           "LAST.EDITED.DATE", "OBSERVER.ID")) %>% 
  mutate(LAST.EDITED.DATE = ymd_hms(LAST.EDITED.DATE)) %>% 
  group_by(OBSERVER.ID) %>% 
  arrange(LAST.EDITED.DATE) %>% 
  ungroup() %>% 
  distinct(OBSERVER.ID, .keep_all = TRUE) %>% # keeps first row (rows are already arranged by date)
  mutate(LE.YEAR = year(LAST.EDITED.DATE),
         LE.MONTH = month(LAST.EDITED.DATE)) %>% 
  mutate(COVID = case_when(LE.YEAR %in% c(2018,2019) ~ "BEF", 
                           LE.YEAR %in% c(2020,2021) ~ "DUR",
                           TRUE ~ "NA")) 

### filtering

new_obsr_data <- new_obsr_data %>% anti_join(filtGA) 

# save(new_obsr_data, file = "data/new_obsr_data.RData")
# 
# load("data/new_obsr_data.RData")


new_obsr_nw <- new_obsr_data %>% 
  filter(COVID != "NA") %>% 
  group_by(COVID, LE.YEAR, LE.MONTH) %>% 
  summarise(NEW = n_distinct(OBSERVER.ID))

# not new to the state, but new overall, just grouped by state of first eBirding
new_obsr_sw <- new_obsr_data %>% 
  filter(COVID != "NA") %>% 
  group_by(COVID, LE.YEAR, LE.MONTH, STATE) %>% 
  summarise(NEW = n_distinct(OBSERVER.ID))


(ggplot(new_obsr_nw, aes(LE.MONTH, NEW, colour = COVID)) + 
    scale_x_continuous(breaks = 1:12) + scale_y_log10() +
    geom_point(size = 2) +
    ggtitle("Monthly number of new eBirders at national level")) /
((ggplot(filter(new_obsr_sw, STATE == "Kerala"), 
           aes(LE.MONTH, NEW, colour = COVID)) + 
      scale_x_continuous(breaks = 1:12) + scale_y_log10() +
      geom_point(size = 2) +
      ggtitle("Kerala"))  |
(ggplot(filter(new_obsr_sw, STATE == "Karnataka"), 
             aes(LE.MONTH, NEW, colour = COVID)) + 
        scale_x_continuous(breaks = 1:12) + scale_y_log10() +
        geom_point(size = 2) +
        ggtitle("Karnataka"))) /
((ggplot(filter(new_obsr_sw, STATE == "Gujarat"), 
           aes(LE.MONTH, NEW, colour = COVID)) + 
      scale_x_continuous(breaks = 1:12) + scale_y_log10() +
      geom_point(size = 2) +
      ggtitle("Gujarat")) |
(ggplot(filter(new_obsr_sw, STATE == "Arunachal Pradesh"), 
           aes(LE.MONTH, NEW, colour = COVID)) + 
      scale_x_continuous(breaks = 1:12) + scale_y_log10() +
      geom_point(size = 2) +
      ggtitle("Arunachal Pradesh")))

