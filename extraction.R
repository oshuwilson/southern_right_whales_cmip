#script for extracting environmental variables to tracks and pseudoabsences
setwd("~/OneDrive - University of Southampton/Documents/Southern Right Whales")


{
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyterra)
  library(ggplot2)
}

rm(list=ls())

#specify region
this.pop <- "NZ"

#read in tracks and buffers
tracks <- read.csv(paste0("data/", this.pop, "_SRW_SSM_track_data.csv"))
buffers <- read.csv(paste0("output/buffers/", this.pop, "_buffers.csv"))

#format for extraction - keep relevant columns
tracks <- tracks %>% select(id, date, lon, lat)
buffers <- buffers %>% select(-X)

#format
tracks <- tracks %>% rename(x = lon, y = lat)
tracks$date <- as_datetime(tracks$date)

#---------------
#Static Variables

###Depth###
depth <- rast("D:/Satellite_Data/static/depth/depth.nc")

#create SpatVector for tracks and buffers
tracks <- vect(tracks,
               geom=c("x", "y"),
               crs=crs(depth)) #this ensures crs are the same as rasters
buffers <- vect(buffers,
                geom=c("x", "y"),
                crs=crs(depth))

#extract
tracks$depth <- extract(depth, tracks, ID=F)
buffers$depth <- extract(depth, buffers, ID=F)

#remove rows where depth is NA - will be NA for every GLORYS variable
tracks <- tracks %>% drop_na(depth)
buffers <- buffers %>% drop_na(depth)


###Slope###
slope <- rast("D:/Satellite_Data/static/slope/slope.nc")
tracks$slope <- extract(slope, tracks, ID=F)
buffers$slope <- extract(slope, buffers, ID=F)

###dShelf### - does this work beyond 40 south?
dshelf <- rast("D:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
tracks$dshelf <- extract(dshelf, tracks, ID=F)
buffers$dshelf <- extract(dshelf, buffers, ID=F)

#cleanup static
rm(depth, slope, dshelf)


#---------------
#Dynamic Variables

#dynamic_extract function from 05a script
source("code/dynamic_extract_function.R")


###SST###
tracks <- dynamic_extract(predictor = "sst", tracks)
buffers <- dynamic_extract(predictor = "sst", buffers)

###MLD###
tracks <- dynamic_extract(predictor = "mld", tracks)
buffers <- dynamic_extract(predictor = "mld", buffers)

###SAL###
tracks <- dynamic_extract(predictor = "sal", tracks)
buffers <- dynamic_extract(predictor = "sal", buffers)

###SSH###
tracks <- dynamic_extract(predictor = "ssh", tracks)
buffers <- dynamic_extract(predictor = "ssh", buffers)

###SIC###
tracks <- dynamic_extract(predictor = "sic", tracks)
tracks$sic[is.na(tracks$sic)] <- 0 #SIC values of 0 print as NA in GLORYS
buffers <- dynamic_extract(predictor = "sic", buffers)
buffers$sic[is.na(buffers$sic)] <- 0

###CURR###
tracks <- dynamic_extract(predictor = "uo", tracks) #uo is eastward velocity
buffers <- dynamic_extract(predictor = "uo", buffers) 

tracks <- dynamic_extract(predictor = "vo", tracks) #vo is northwards velocity
buffers <- dynamic_extract(predictor = "vo", buffers)

tracks$curr <- sqrt((tracks$uo^2) + (tracks$vo^2)) #current speed
buffers$curr <- sqrt((buffers$uo^2) + (buffers$vo^2))

###CHL###
#uses different function for resampled files
source("code/dynamic_chlorophyll_function.R") #unique function for different file structure

tracks <- dynamic_chlorophyll(predictor = "chl", tracks)
buffers <- dynamic_chlorophyll(predictor = "chl", buffers)

###WIND###
#uses different function for monthly file structure
source("code/dynamic_wind_function.R")

tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "east")
tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "north")
tracks$wind <- sqrt(tracks$wind_east^2 + tracks$wind_north^2)

buffers <- dynamic_wind(predictor = "wind", buffers, direction = "east")
buffers <- dynamic_wind(predictor = "wind", buffers, direction = "north")
buffers$wind <- sqrt(buffers$wind_east^2 + buffers$wind_north^2)


#---------------
#Export
plot(tracks, pch=".")
tracks <- as.data.frame(tracks, geom="XY")

plot(buffers, pch=".")
buffers <- as.data.frame(buffers, geom="XY")

write.csv(tracks, 
          file=paste0("output/extraction/", this.pop, "/presences.csv"))

write.csv(buffers, 
          file=paste0("output/extraction/", this.pop, "/buffers.csv"))