## from github.com/ashwinv2005/trend-analyses/blob/master/functions.R

#################### create maps and grids ###################################

## requires shapefiles and several packages
## provide path to folder and (common) name of files within
## saves a workspace image called "maps.RData"

## this can be edited for more flexibility with grid sizes; current default is 25,50,100,200


create_maps <- function(g1 = 25, g2 = 50, g3 = 100, g4 = 200,
                        path1 = "data/in_2011", name1 = "India_2011",
                        path2 = "data/in_states_2019", name2 = "in_states_2019",
                        path3 = "data/in_dists_2019", name3 = "in_dist_2019")
{
  require(tidyverse)
  require(rgdal)
  require(sp)
  require(sf)
  
  # reading maps
  
  assign("indiamap",readOGR(path1,name1),.GlobalEnv)
  
  assign("statemap",readOGR(path2,name2),.GlobalEnv)
  names(statemap@data)[2] <- "ST_NM"
  assign("statemap",statemap,.GlobalEnv)
  
  assign("districtmap",readOGR(path3,name3),.GlobalEnv)
  names(districtmap@data)[1:2] <- c("DISTRICT","ST_NM")
  assign("districtmap",districtmap,.GlobalEnv)

  # creating SPDF grids below that can be intersected with various maps and overlaid on to data
  
  bb <- bbox(indiamap) # creates a box with extents from map
  cs <- c(g1*1000/111111,g1*1000/111111)  # cell size g1 km x g1 km
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd) # create required grids
  sp_grd <- SpatialGridDataFrame(grd, data = data.frame(id = 1:prod(cd))) # create spatial grid data frame
  sp_grd_poly <- as(sp_grd, "SpatialPolygonsDataFrame") # SGDF to SPDF
  assign("gridmapg1",sp_grd_poly,.GlobalEnv)
  
  bb <- bbox(indiamap)
  cs <- c(g2*1000/111111,g2*1000/111111) 
  cc <- bb[, 1] + (cs/2)  
  cd <- ceiling(diff(t(bb))/cs)  
  grd <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
  sp_grd <- SpatialGridDataFrame(grd, data = data.frame(id = 1:prod(cd)))
  sp_grd_poly <- as(sp_grd, "SpatialPolygonsDataFrame")
  assign("gridmapg2",sp_grd_poly,.GlobalEnv)
  
  bb <- bbox(indiamap)
  cs <- c(g3*1000/111111,g3*1000/111111) 
  cc <- bb[, 1] + (cs/2)  
  cd <- ceiling(diff(t(bb))/cs)  
  grd <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
  sp_grd <- SpatialGridDataFrame(grd, data = data.frame(id = 1:prod(cd)))
  sp_grd_poly <- as(sp_grd, "SpatialPolygonsDataFrame")
  assign("gridmapg3",sp_grd_poly,.GlobalEnv)
  
  bb <- bbox(indiamap)
  cs <- c(g4*1000/111111,g4*1000/111111) 
  cc <- bb[, 1] + (cs/2)  
  cd <- ceiling(diff(t(bb))/cs)  
  grd <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
  sp_grd <- SpatialGridDataFrame(grd, data = data.frame(id = 1:prod(cd)))
  sp_grd_poly <- as(sp_grd, "SpatialPolygonsDataFrame")
  assign("gridmapg4",sp_grd_poly,.GlobalEnv)
  
  # indiamap = spTransform(indiamap,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # not required here, CRS is NA
  
  assign("gridlevels",c(g1,g2,g3,g4),.GlobalEnv)
  
  rm(list = setdiff(ls(envir = .GlobalEnv), 
                    c("districtmap", "statemap", "indiamap", "gridmapg1", "gridmapg2", 
                      "gridmapg3", "gridmapg4", "gridlevels")), 
     pos = ".GlobalEnv")
  
  save.image("data/maps.RData")
  
  rm(districtmap, statemap, indiamap, gridmapg1, gridmapg2, gridmapg3, gridmapg4, gridlevels, 
     pos = ".GlobalEnv")
}



#################### combine map and grid data with raw data ###########################

## prepare data for analyses, add map variables, grids
## place the 'maps' workspace in working directory

add_mapvars <- function(datapath = "data/data_sub.RData", mappath = "data/maps.RData")
{
  require(tidyverse)
  require(data.table)
  require(sp)
  require(rgeos)
  
  load(datapath)
  load(mappath)
  load("data/clips.RData")
  data_map <- data_sub
  
  # add columns with DISTRICT and ST_NM to main data 
  # same group ID => same grid/district/state 
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER") # to setup for the future left_join
  coordinates(temp) <- ~LONGITUDE + LATITUDE # convert to SPDF?
  proj4string(temp) <- "+proj=longlat +datum=WGS84"
  temp <- temp %>% over(districtmap) %>% # returns only ATTRIBUTES of districtmap (DISTRICT and ST_NM)
    select(1:2) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") # add column to join with the main data
  data_map <- left_join(temp, data_map)
  
  # add columns with GRID ATTRIBUTES to main data
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER")
  coordinates(temp) <- ~LONGITUDE + LATITUDE
  temp <- temp %>% over(gridmapg1) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") 
  data_map <- left_join(temp, data_map)
  names(data_map)[2] <- "GRIDG1"
  
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER")
  coordinates(temp) <- ~LONGITUDE + LATITUDE
  temp <- temp %>% over(gridmapg2) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") 
  data_map <- left_join(temp, data_map)
  names(data_map)[2] <- "GRIDG2"
  
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER")
  coordinates(temp) <- ~LONGITUDE + LATITUDE
  temp <- temp %>% over(gridmap3) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") 
  data_map <- left_join(temp, data_map)
  names(data_map)[2] <- "GRIDG3"
  
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER")
  coordinates(temp) <- ~LONGITUDE + LATITUDE
  temp <- temp %>% over(gridmapg4) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") 
  data_map <- left_join(temp, data_map)
  names(data_map)[2] <- "GRIDG4"
  
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER")
  coordinates(temp) <- ~LONGITUDE + LATITUDE
  temp <- temp %>% over(g2clip) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") 
  data_map <- left_join(temp, data_map)
  names(data_map)[2] <- "G2CLIP"
  
  temp <- data_map %>% group_by(SAMPLING.EVENT.IDENTIFIER) %>% slice(1)
  temp <- temp %>% column_to_rownames("SAMPLING.EVENT.IDENTIFIER")
  coordinates(temp) <- ~LONGITUDE + LATITUDE
  temp <- temp %>% over(g3clip) %>% rownames_to_column("SAMPLING.EVENT.IDENTIFIER") 
  data_map <- left_join(temp, data_map)
  names(data_map)[2] <- "G3CLIP"
  
  ## 
  
  assign("gridlevels", gridlevels, .GlobalEnv)
  
  assign("data_map", data_map, .GlobalEnv)
  rm(list = setdiff(ls(envir = .GlobalEnv), c("data_map","gridlevels")), pos = ".GlobalEnv")
  
  save.image("data_map.RData")
  rm(data_map, gridlevels, pos = ".GlobalEnv")
  
}


