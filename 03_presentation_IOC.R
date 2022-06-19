library(tidyverse)
library(lubridate)
library(patchwork)
library(kableExtra)
library(rasterVis)
library(RColorBrewer)

# boundary of India 
india <- rgdal::readOGR(dsn = "data/in_2011", layer = "India_2011")
sp::proj4string(india) <- "+proj=longlat +datum=WGS84"
india <- terra::vect(india)
# state boundaries
states <- rgdal::readOGR(dsn = "data/in_states_2019", layer = "in_states_2019") %>% 
  terra::vect()

load("data/data0_slice.RData")
load("data/rast_UNU.RData")
load("data/rast_SoIB.RData")

source("scripts/functions.R")

# choosing only months with data from all 3 COVID categories
month_compar <- data0_slice_S %>% 
  group_by(MONTH) %>% 
  mutate(COVID = as.character(COVID)) %>% 
  mutate(COVID = factor(case_when(COVID %in% c("DUR_20","DUR_21") ~ "DUR",
                                  TRUE ~ COVID),
                        levels = c("BEF","DUR","AFT"))) %>% 
  dplyr::summarise(N = n_distinct(COVID)) %>% 
  filter(N == 3) %>%  # all 3 COVID categories
  dplyr::select(MONTH)



### group birding per observer ####

temp1 <- data0_slice_S %>% 
  group_by(COVID, YEAR, MONTH, STATE, OBSERVER.ID) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  ungroup() 

temp2 <- data0_slice_S %>% 
  filter(NUMBER.OBSERVERS > 1) %>% 
  group_by(COVID, YEAR, MONTH, STATE, OBSERVER.ID) %>% 
  dplyr::summarise(NO.SLISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  ungroup()

nl_po_sw <- temp1 %>% 
  left_join(temp2) %>% 
  mutate(NO.SLISTS = replace_na(NO.SLISTS, 0)) %>% 
  mutate(PROP.SLISTS = NO.SLISTS/NO.LISTS,
         SE = sqrt((PROP.SLISTS)*(1 - PROP.SLISTS)/NO.LISTS)) %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(PROP.SLISTS = mean(PROP.SLISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.SLISTS - 1.96*SE,
                   CI.U = PROP.SLISTS + 1.96*SE)

nl_po_nw <- nl_po_sw %>% 
  group_by(COVID, YEAR, MONTH) %>% 
  dplyr::summarise(PROP.SLISTS = mean(PROP.SLISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.SLISTS - 1.96*SE,
                   CI.U = PROP.SLISTS + 1.96*SE)


rm(list = c("temp1","temp2","nl_po_sw"))

### hotspot birding ####

temp1 <- data0_slice_G %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  ungroup() 

temp2 <- data0_slice_G %>% 
  filter(LOCALITY.TYPE == "H") %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(NO.HLISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  ungroup()

hot_sw <- temp1 %>% 
  left_join(temp2) %>% 
  mutate(NO.HLISTS = replace_na(NO.HLISTS, 0)) %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  mutate(PROP.HLISTS = NO.HLISTS/NO.LISTS,
         SE = sqrt((PROP.HLISTS)*(1 - PROP.HLISTS)/NO.LISTS),
         CI.L = PROP.HLISTS - 1.96*SE,
         CI.U = PROP.HLISTS + 1.96*SE) 

hot_nw <- hot_sw %>%
  filter(NO.LISTS != 0) %>% 
  group_by(COVID, YEAR, MONTH) %>% 
  dplyr::summarise(PROP.HLISTS = mean(PROP.HLISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.HLISTS - 1.96*SE,
                   CI.U = PROP.HLISTS + 1.96*SE) 

rm(list = c("temp1","temp2","hot_sw"))

### protocol ####

temp1 <- data0_slice_G %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  ungroup() 

temp2 <- data0_slice_G %>% 
  filter(PROTOCOL.TYPE == "Traveling" & EFFORT.DISTANCE.KM >= 0.1) %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(NO.TLISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  ungroup()

prot_sw <- temp1 %>% 
  left_join(temp2) %>% 
  mutate(NO.TLISTS = replace_na(NO.TLISTS, 0)) %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  mutate(PROP.TLISTS = NO.TLISTS/NO.LISTS,
         SE = sqrt((PROP.TLISTS)*(1 - PROP.TLISTS)/NO.LISTS),
         CI.L = PROP.TLISTS - 1.96*SE,
         CI.U = PROP.TLISTS + 1.96*SE) 

prot_nw <- prot_sw %>%
  filter(NO.LISTS != 0) %>% 
  group_by(COVID, YEAR, MONTH) %>% 
  dplyr::summarise(PROP.TLISTS = mean(PROP.TLISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.TLISTS - 1.96*SE,
                   CI.U = PROP.TLISTS + 1.96*SE) 

rm(list = c("temp1","temp2","prot_sw"))

### site fidelity per observer ####

set.seed(213) # for bootstrap
fidel_sw <- data0_slice_S %>% 
  # grouping in terms of 25kmx25km cells (not 1degx1deg)
  mutate(LAT.25K = LATITUDE*111/25, 
         LON.25K = LONGITUDE*111/25) %>% 
  # pasting lat-lon as string to find truly unique sites
  mutate(COORD.25K = paste0(ceiling(LON.25K), ",", ceiling(LAT.25K))) %>% 
  group_by(COVID, YEAR, MONTH, STATE, OBSERVER.ID) %>% 
  dplyr::summarise(NO.SITES = n_distinct(COORD.25K)) %>% 
  mutate(NO.SITES = replace_na(NO.SITES, 0)) %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(NO.SITES = boot_conf(x = NO.SITES)) %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  dplyr::summarise(SE = sd(NO.SITES),
                   NO.SITES = mean(NO.SITES),
                   CI.L = NO.SITES - 1.96*SE,
                   CI.U = NO.SITES + 1.96*SE)

fidel_nw <- fidel_sw %>% 
  ungroup() %>% 
  group_by(COVID, YEAR, MONTH) %>% 
  dplyr::summarise(NO.SITES = mean(NO.SITES),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = NO.SITES - 1.96*SE,
                   CI.U = NO.SITES + 1.96*SE)

rm(fidel_sw)


### Spatial cover ####

cover <- data0_slice_G %>% 
  ungroup() %>% 
  mutate(TOT.CELLS = n_distinct(CELL.ID)) %>% 
  group_by(COVID, YEAR, MONTH, CELL.ID) %>% 
  dplyr::summarise(TOT.CELLS = min(TOT.CELLS),
                   N.LISTS = n_distinct(GROUP.ID)) %>% 
  filter(N.LISTS >= 10) %>% 
  group_by(COVID, YEAR, MONTH) %>% 
  dplyr::summarise(TOT.CELLS = min(TOT.CELLS),
                   NO.CELLS = n_distinct(CELL.ID),
                   PROP.CELLS = NO.CELLS/TOT.CELLS,
                   SE = sqrt((PROP.CELLS)*(1 - PROP.CELLS)/TOT.CELLS),
                   CI.L = PROP.CELLS - 1.96*SE,
                   CI.U = PROP.CELLS + 1.96*SE)


### Spatial net change in effort ####

temp <- data0_slice_G %>% 
  ungroup() %>% 
  filter(MONTH %in% month_compar$MONTH) %>% 
  group_by(COVID, CELL.ID, MONTH) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(GROUP.ID)) %>% 
  # summarise removes one grouping variable so no need to repeat group_by() 
  dplyr::summarise(NO.LISTS = mean(NO.LISTS)) %>% 
  ungroup()

temp0 <- temp %>% distinct(CELL.ID)

temp1 <- temp %>% 
  filter(COVID == "BEF") %>% 
  # to include all grid cells and not just those with birding in BEF
  right_join(temp0) %>% 
  mutate(NO.LISTS = replace_na(NO.LISTS, replace = 0)) %>% 
  rename(BEF = NO.LISTS) %>% 
  dplyr::select(CELL.ID, BEF)

temp2 <- temp %>% 
  filter(COVID %in% c("DUR_20","DUR_21")) %>% 
  group_by(CELL.ID) %>% 
  # averages cells with lists from both years, and returns single year value for others (n/1)
  dplyr::summarise(NO.LISTS = mean(NO.LISTS)) %>% 
  ungroup() %>% 
  right_join(temp0) %>% 
  mutate(NO.LISTS = replace_na(NO.LISTS, replace = 0)) %>% 
  rename(DUR = NO.LISTS) %>% 
  dplyr::select(CELL.ID, DUR)

temp3 <- temp %>% 
  filter(COVID == "AFT") %>% 
  # to include all grid cells and not just those with birding in BEF
  right_join(temp0) %>% 
  mutate(NO.LISTS = replace_na(NO.LISTS, replace = 0)) %>% 
  rename(AFT = NO.LISTS) %>% 
  dplyr::select(CELL.ID, AFT)


set.seed(23) # for bootstrap
net_effort <- temp0 %>% 
  left_join(temp1) %>% 
  left_join(temp2) %>% 
  left_join(temp3) %>% 
  group_by(CELL.ID) %>% 
  dplyr::summarise(BeforeDuring = rast_propchange(BEF, DUR, k = 1),
                   DuringAfter = rast_propchange(DUR, AFT, k = 1)) %>% 
  pivot_longer(cols = c(BeforeDuring, DuringAfter),
               names_to = "PERIOD", 
               values_to = "PROP.CHANGE") %>%  # log. prop. change in monthly no.lists per 24x24
  group_by(PERIOD) %>% 
  # this gives B number of means/medians of bootstrapped samples
  # median resulting in just 1s so mean is fine
  dplyr::summarise(PROP.CHANGE = boot_conf(x = PROP.CHANGE)) %>% 
  group_by(PERIOD) %>% 
  dplyr::summarise(CI.L = stats::quantile(PROP.CHANGE, 0.025), # Obtain the CIs
                   CI.U = stats::quantile(PROP.CHANGE, 0.975),
                   PROP.CHANGE = mean(PROP.CHANGE))


rm(temp1, temp2, temp3, temp0, temp)

### Spatial proportion of urban birding ####

temp <- data0_slice_G %>% 
  group_by(COVID, YEAR, MONTH, CELL.ID) %>% # state is meaninglesss cos CELL.ID
  filter(URBAN == 1) %>% 
  dplyr::summarise(U.LISTS = n_distinct(GROUP.ID)) %>% 
  ungroup()

urban_lists <- data0_slice_G %>% 
  group_by(COVID, YEAR, MONTH, CELL.ID) %>%
  dplyr::summarise(TOT.LISTS = n_distinct(GROUP.ID)) %>% 
  left_join(temp) %>% 
  group_by(COVID, YEAR, MONTH, CELL.ID) %>%
  dplyr::summarise(U.LISTS = replace_na(U.LISTS, 0),
                   PROP.U = U.LISTS/TOT.LISTS,
                   SE = sqrt((PROP.U)*(1 - PROP.U)/TOT.LISTS)) %>% 
  group_by(COVID, YEAR, MONTH) %>% 
  dplyr::summarise(PROP.U = mean(PROP.U),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.U - 1.96*SE,
                   CI.U = PROP.U + 1.96*SE)


rm(temp)

### Spatial spread by period (map) ####

data_SoIBmap0 <- data0_slice_G %>% 
  filter(MONTH %in% month_compar$MONTH) %>% 
  mutate(CELL.LONG = raster::xyFromCell(rast_SoIB, CELL.ID)[,"x"],
         CELL.LAT = raster::xyFromCell(rast_SoIB, CELL.ID)[,"y"]) %>%
  group_by(COVID, CELL.ID, CELL.LONG, CELL.LAT, MONTH) %>% 
  dplyr::summarise(NO.LISTS = n_distinct(GROUP.ID)) %>% 
  # summarise removes one grouping variable so no need to repeat group_by() 
  dplyr::summarise(NO.LISTS = mean(NO.LISTS)) %>% 
  ungroup()


temp0 <- data_SoIBmap0 %>% distinct(CELL.ID, CELL.LONG, CELL.LAT)


data_SoIBmap1 <- data_SoIBmap0 %>% 
  filter(COVID == "BEF") %>% 
  # to include all grid cells and not just those with birding in BEF
  right_join(temp0) %>% 
  mutate(COVID = "BEF", 
         NO.LISTS = replace_na(NO.LISTS, replace = 0)) %>% 
  dplyr::select(CELL.LONG, CELL.LAT, NO.LISTS)

rast_SoIBmap1 <- raster::rasterFromXYZ(xyz = data_SoIBmap1,
                                       crs = raster::crs(rast_SoIB))


data_SoIBmap2 <- data_SoIBmap0 %>% 
  filter(COVID %in% c("DUR_20","DUR_21")) %>% 
  group_by(CELL.ID, CELL.LONG, CELL.LAT) %>% 
  # averages cells with lists from both years, and returns single year value for others (n/1)
  dplyr::summarise(NO.LISTS = mean(NO.LISTS)) %>% 
  ungroup() %>% 
  right_join(temp0) %>% 
  mutate(COVID = "DUR", 
         NO.LISTS = replace_na(NO.LISTS, replace = 0)) %>% 
  dplyr::select(CELL.LONG, CELL.LAT, NO.LISTS)

rast_SoIBmap2 <- raster::rasterFromXYZ(xyz = data_SoIBmap2,
                                       crs = raster::crs(rast_SoIB))


data_SoIBmap3 <- data_SoIBmap0 %>% 
  filter(COVID == "AFT") %>% 
  right_join(temp0) %>% 
  mutate(COVID = "AFT", 
         NO.LISTS = replace_na(NO.LISTS, replace = 0)) %>% 
  dplyr::select(CELL.LONG, CELL.LAT, NO.LISTS)

rast_SoIBmap3 <- raster::rasterFromXYZ(xyz = data_SoIBmap3,
                                       crs = raster::crs(rast_SoIB))



rast_SoIBmap <- raster::stack(
  raster::overlay(x = rast_SoIBmap1, y = rast_SoIBmap2, 
                  k = 1, # constant to prevent NAs
                  fun = rast_logpropchange),
  raster::overlay(x = rast_SoIBmap2, y = rast_SoIBmap3, 
                  k = 1, 
                  fun = rast_logpropchange))
names(rast_SoIBmap) <- c("Before-During", "During-After")



library(RColorBrewer)


in_bound <- rbind(india, states) %>% as("Spatial") # converting for sp::sp.polygons()

# for graph legends
plotfrom <- raster::cellStats(rast_SoIBmap, "min")[1]
plotto <- raster::cellStats(rast_SoIBmap, "max")[1]
plotlegend <- seq(round(plotfrom, 2), round(plotto, 2), 0.1)

labellegend <- data.frame(NORMAL = round(seq(exp(plotfrom), 
                                             exp(plotto),  
                                             length.out = 10), 3)) %>% 
  mutate(LOG = round(log(NORMAL), 2),
         LABEL = NORMAL)

mytheme <- rasterVis::rasterTheme(region = viridis::viridis(100),
                                  box.rectangle = list(col = "white", fill = "black"),
                                  box.umbrella = list(col = 'black', lty = 1),
                                  panel.background = list(col = "#F0F0F0"),
                                  strip.background = list(col = "transparent"),
                                  strip.border = list(col = "transparent"))

levelplot(rast_SoIBmap, 
          par.settings = mytheme, 
          margin = F, 
          at = plotlegend,
          main = list(label = "Proportional change in monthly (seven months) no. of lists per 24km \u00d7 24km cell",
                      cex = 1.5),
          xlab = list(label = "Longitude", cex = 1.3),
          ylab = list(label = "Latitude", cex = 1.3),
          scales = list(cex = 1),
          names.attr = c("Before to During", "During to After"),
          par.strip.text = list(cex = 1.2),
          colorkey = list(labels = list(at = labellegend$LOG,
                                        labels = labellegend$LABEL))) +
  latticeExtra::layer(sp::sp.polygons(in_bound, alpha = 0.4)) -> s_spread_SoIBmap


png(file = "pics/03_presentation_IOC_map.png", 
    width = 16, height = 9, units = "in", res = 300)
print(s_spread_SoIBmap)
dev.off()



### Temporal patterns ####

t_dow_sw <- data0_slice_G %>% 
  filter(MONTH %in% month_compar$MONTH) %>% 
  mutate(DAY.W = wday(OBSERVATION.DATE, 
                      week_start = getOption("lubridate.week.start", 1))) %>% 
  group_by(COVID, STATE, YEAR, MONTH) %>% 
  mutate(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  group_by(COVID, DAY.W, STATE, YEAR, MONTH) %>% 
  dplyr::summarise(NO.LISTS = min(NO.LISTS),
                   DAY.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER)),
                   PROP.LISTS = DAY.LISTS/NO.LISTS,
                   SE = sqrt((PROP.LISTS)*(1 - PROP.LISTS)/NO.LISTS)) %>% 
  group_by(COVID, DAY.W, STATE) %>% 
  dplyr::summarise(PROP.LISTS = mean(PROP.LISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.LISTS - 1.96*SE,
                   CI.U = PROP.LISTS + 1.96*SE)

t_dow_nw <- t_dow_sw %>% 
  group_by(COVID, DAY.W) %>% 
  dplyr::summarise(PROP.LISTS = mean(PROP.LISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.LISTS - 1.96*SE,
                   CI.U = PROP.LISTS + 1.96*SE)



t_tod_sw <- data0_slice_G %>% 
  filter(MONTH %in% month_compar$MONTH) %>% 
  group_by(COVID, STATE, YEAR, MONTH) %>% 
  mutate(NO.LISTS = n_distinct(SAMPLING.EVENT.IDENTIFIER)) %>%
  group_by(COVID, STATE, HOUR, YEAR, MONTH) %>% 
  dplyr::summarise(NO.LISTS = min(NO.LISTS),
                   TIME.LISTS = replace_na(n_distinct(SAMPLING.EVENT.IDENTIFIER)),
                   PROP.LISTS = TIME.LISTS/NO.LISTS,
                   SE = sqrt((PROP.LISTS)*(1 - PROP.LISTS)/NO.LISTS)) %>% 
  group_by(COVID, STATE, HOUR) %>% 
  dplyr::summarise(PROP.LISTS = mean(PROP.LISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.LISTS - 1.96*SE,
                   CI.U = PROP.LISTS + 1.96*SE)

t_tod_nw <- t_tod_sw %>%
  group_by(COVID, HOUR) %>% 
  dplyr::summarise(PROP.LISTS = mean(PROP.LISTS),
                   SE = sqrt(sum((SE)^2))/n(),
                   CI.L = PROP.LISTS - 1.96*SE,
                   CI.U = PROP.LISTS + 1.96*SE)

rm(t_dow_nw, t_tod_nw)

### saving RData ####

save(month_compar, nl_po_nw, hot_nw, prot_nw, fidel_nw, urban_lists, cover,
     net_effort, t_dow_sw, t_tod_sw, file = "data/03_presentation_IOC.RData")

