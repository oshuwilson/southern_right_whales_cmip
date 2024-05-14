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
  library(amt)
  library(ggplot2)
  library(tictoc)
}

#define population
this.pop <- "NZ"

#read in tracks for this population
tracks <- read.csv(paste0("data/", this.pop, "_SRW_SSM_track_data.csv"))

#plot tracks
tracks_terra <- vect(tracks,
                     geom = c("x", "y"),
                     crs = "epsg:4326")
plot(tracks_terra, pch=".")

#only keep necessary variables
tracks <- tracks %>% select(id, date, lon, lat)

#format date and id column
tracks$id <- as.factor(tracks$id)
tracks$date <- as_datetime(tracks$date)



# 1. Calculate step lengths

#make sure all rows are complete with necessary info
tracks <- tracks[complete.cases(tracks),]

#make tracks
trks <- make_track(tracks, lon, lat, date, id = id,
                   crs=4326)

#transform coordinates to meters from lat lon
trks <- transform_coords(trks, 6932)

#nest by ID
trks <- trks |> nest(data = -"id")

#resample to daily rate and run steps_by_burst for all individuals for step lengths
trks2 <- trks |>
  mutate(steps = map(data, function(x)
    x |> track_resample(rate = days(1), tolerance = hours(1)) |> steps_by_burst()))

#unnest for calculations
trks3 <- trks2 |> amt::select(id, steps) |> unnest(cols = steps)

#visualise
ggplot(trks3, aes(x=sl_)) + geom_density()

#stats - 75th percentile
max.buffer <- as.numeric(quantile(trks3$sl_, probs=0.75))

#remove amt variables
rm(trks, trks2, trks3)



# 2. Create Buffer Samples

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

#download ocean mask
oceans <- ne_download(scale = "medium", category = "physical", type = "ocean", returnclass = "sf")

#run buffer creation script
allbuff <- NULL

tic()
for(i in 1:nrow(tracks)){
  
  tracks_test <- tracks[i,]
  
suppressMessages( #suppress message - trust me you don't want messages
  buffers <- spatiotemp_pseudoabs(occ.data = tracks_test,
                                  spatial.method = "buffer",
                                  spatial.ext = oceans, 
                                  temporal.method = "buffer",
                                  spatial.buffer = c(12000, max.buffer), #12000 ensures that pseudoabs falls in a different cell
                                  temporal.buffer = 0,
                                  n.pseudoabs = 1,
                                  prj = "+proj=longlat +datum=WGS84"))

allbuff <- rbind(allbuff, buffers)  

print(i)
}
toc()

allbuff$id <- tracks$id
allbuff$date <- tracks$date

#plot
plot(vect(allbuff[, c("x", "y")],
                  geom = c("x", "y"),
                  crs = "+proj=longlat +datum=WGS84"),
          pch = ".", col = "red") 
plot(tracks_terra, col="black", pch=".", add=T)

# 5. Export allbuff
#only when happy export
write.csv(allbuff, paste0("output/buffers/", this.pop, "_buffers.csv"))

