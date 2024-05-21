rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Southern Right Whales")

{
  library(dynamicSDM)
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
  library(ggplot2)
  library(tictoc)
}

#define population
this.pop <- "ARG"

#read in tracks for this population
tracks <- read.csv(paste0("data/", this.pop, "_SRW_SSM_track_data.csv"))

#download ocean mask
oceans <- ne_download(scale = "medium", category = "physical", type = "ocean", returnclass = "sf")
oceans <- vect(oceans)

#plot tracks
tracks_terra <- vect(tracks,
                     geom = c("lon", "lat"),
                     crs = "epsg:4326")

#project oceans and tracks_terra
e <- ext(-180, 180, -86, -25)
oceans <- crop(oceans, e)
oceans <- project(oceans, "EPSG:6932")
plot(oceans)

tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch=".", add=T)

#only keep necessary variables
tracks <- tracks %>% select(id, date, lon, lat)

#format date and id column
tracks$id <- as.factor(tracks$id)
tracks$date <- as_datetime(tracks$date)

# 1. Create MCP Mask 
mch <- convHull(tracks_terra)
plot(mch)

mch_buff <- buffer(mch, width=84000)
plot(mch_buff, add=T)

mch_buff <- st_as_sf(mch_buff)
oceans <- st_as_sf(oceans)
oceans <- st_buffer(oceans,0)
mask <- st_intersection(mch_buff, oceans$geometry)
mask <- st_buffer(mask, 0)
plot(mask)

mask <- st_transform(mask, 4326)
plot(mask)
#alternative for non-projection method - can cause fatal error
# mask <- terra::intersect(mch_buff, oceans)
# plot(mask)
# mask <- st_as_sf(mask)
# plot(mask)

# 2. Create Background Samples

#format day/month/year column for dynamicSDM
tracks$day <- as.numeric(day(tracks$date))
tracks$month <- month(tracks$date)
tracks$year <- year(tracks$date)

tracks <- tracks %>% rename(x = lon, y = lat)

#filter to remove duplicates or invalid coordinates and dates
tracks <- spatiotemp_check(occ.data = tracks,
                           na.handle = "exclude",
                           date.handle = "exclude",
                           date.res = "day",
                           coord.handle = "exclude",
                           duplicate.handle = "exclude")



#run buffer creation script
allback <- NULL

tic()
for(i in 1:nrow(tracks)){
  
  tracks_test <- tracks[i,]
  
  suppressMessages( #suppress message - trust me you don't want messages
    background <- spatiotemp_pseudoabs(occ.data = tracks_test,
                                    spatial.method = "random",
                                    spatial.ext = mask, 
                                    temporal.method = "buffer",
                                    temporal.buffer = 0,
                                    n.pseudoabs = 1,
                                    prj = "+proj=longlat +datum=WGS84"))
  
  allback <- rbind(allback, background)  
  
  print(i)
}
toc()

allback$id <- tracks$id
allback$date <- tracks$date

#plot
plot(vect(allback[, c("x", "y")],
          geom = c("x", "y"),
          crs = "+proj=longlat +datum=WGS84"),
     pch = ".", col = "red") 
tracks_terra <- project(tracks_terra, "EPSG:4326")
plot(tracks_terra, col="black", pch=".", add=T)

# 5. Export allbuff
#only when happy export
write.csv(allback, paste0("output/background/", this.pop, "_background.csv"))

