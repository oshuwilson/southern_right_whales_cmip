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
this.pop <- "OZ"

#read in tracks and background
tracks <- read.csv(paste0("data/", this.pop, "_SRW_SSM_track_data.csv"))
background <- read.csv(paste0("output/background/", this.pop, "_background.csv"))

#format for extraction - keep relevant columns
tracks <- tracks %>% select(id, date, lon, lat)
background <- background %>% select(-X)

#format
tracks <- tracks %>% rename(x = lon, y = lat)
tracks$date <- as_datetime(tracks$date)

background$date <- as_datetime(background$date)

#---------------
#Static Variables

###Depth###
depth <- rast("D:/Satellite_Data/static/depth/depth.nc")

#create SpatVector for tracks and background
tracks <- vect(tracks,
               geom=c("x", "y"),
               crs=crs(depth)) #this ensures crs are the same as rasters
background <- vect(background,
                geom=c("x", "y"),
                crs=crs(depth))

#extract
tracks$depth <- extract(depth, tracks, ID=F)
background$depth <- extract(depth, background, ID=F)

#remove rows where depth is NA - will be NA for every GLORYS variable
tracks <- tracks %>% drop_na(depth)
background <- background %>% drop_na(depth)


###Slope###
slope <- rast("D:/Satellite_Data/static/slope/slope.nc")
tracks$slope <- extract(slope, tracks, ID=F)
background$slope <- extract(slope, background, ID=F)

###dShelf### - does this work beyond 40 south?
dshelf <- rast("D:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
tracks$dshelf <- extract(dshelf, tracks, ID=F)
background$dshelf <- extract(dshelf, background, ID=F)

#cleanup static
rm(depth, slope, dshelf)


#---------------
#Dynamic Variables

#dynamic_extract function from 05a script
source("code/dynamic_extract_function.R")


###SST###
tracks <- dynamic_extract(predictor = "sst", tracks)
background <- dynamic_extract(predictor = "sst", background)

###MLD###
tracks <- dynamic_extract(predictor = "mld", tracks)
background <- dynamic_extract(predictor = "mld", background)

###SAL###
tracks <- dynamic_extract(predictor = "sal", tracks)
background <- dynamic_extract(predictor = "sal", background)

###SSH###
tracks <- dynamic_extract(predictor = "ssh", tracks)
background <- dynamic_extract(predictor = "ssh", background)

###SIC###
tracks <- dynamic_extract(predictor = "sic", tracks)
tracks$sic[is.na(tracks$sic)] <- 0 #SIC values of 0 print as NA in GLORYS
background <- dynamic_extract(predictor = "sic", background)
background$sic[is.na(background$sic)] <- 0

###CURR###
tracks <- dynamic_extract(predictor = "uo", tracks) #uo is eastward velocity
background <- dynamic_extract(predictor = "uo", background) 

tracks <- dynamic_extract(predictor = "vo", tracks) #vo is northwards velocity
background <- dynamic_extract(predictor = "vo", background)

tracks$curr <- sqrt((tracks$uo^2) + (tracks$vo^2)) #current speed
background$curr <- sqrt((background$uo^2) + (background$vo^2))

###CHL###
#uses different function for resampled files
source("code/dynamic_chlorophyll_function.R") #unique function for different file structure

tracks <- dynamic_chlorophyll(predictor = "chl", tracks)
background <- dynamic_chlorophyll(predictor = "chl", background)

###WIND###
#uses different function for monthly file structure
source("code/dynamic_wind_function.R")

tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "east")
tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "north")
tracks$wind <- sqrt(tracks$wind_east^2 + tracks$wind_north^2)

background <- dynamic_wind(predictor = "wind", background, direction = "east")
background <- dynamic_wind(predictor = "wind", background, direction = "north")
background$wind <- sqrt(background$wind_east^2 + background$wind_north^2)


#---------------
#Export
plot(tracks, pch=".")
tracks <- as.data.frame(tracks, geom="XY")

plot(background, pch=".")
background <- as.data.frame(background, geom="XY")

write.csv(tracks, 
          file=paste0("output/extraction/", this.pop, "/", this.pop, "_presences_extracted.csv"))

write.csv(background, 
          file=paste0("output/extraction/", this.pop, "/", this.pop, "_background_extracted.csv"))

