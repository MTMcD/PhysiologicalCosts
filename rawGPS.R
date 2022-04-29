#############################################################################
# Molly McDermott
# created 12/3/21
# adapted from code written by Sage Madden
# putting together GPS data from raw tag data
#############################################################################


#### SET UP ####
library(tidyverse)
library(dplyr)

setwd("~/SafranLab/BSM_Physiology/gps_data")

# create function for reading in csv and retaining file name
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}


#### READ IN GPS DATA ####
# read in all csv files in working directory using function above
tbl <-
  list.files(pattern = "*.csv") %>% 
  map_df(~read_plus(.))

#format date column
tbl$Date <- as.Date(tbl$`RTC-date`, "%y/%m/%d")

# parse file name to extract tag #, site, nest
tbl$filename <- str_replace(tbl$filename, ".csv", "")
tbl$tag <- as.numeric(str_sub(tbl$filename, 10, 14))
tbl$year <- as.numeric(str_sub(tbl$filename, 16, 19))
tbl$nest <- as.numeric(str_extract(str_sub(tbl$filename, 21, str_length(tbl$filename)), "\\d\\d\\d|\\d\\d|\\d"))
tbl$site <- as.factor(str_split_fixed(
  str_sub(tbl$filename, 21, str_length(tbl$filename)), " ", n = 2)[,1])

#filter GPS data points to only those with high accuracy
tbl <- tbl %>%
  filter(Status == "Valid" & HDOP < 5)

#reset working directory
setwd("~/SafranLab/PhysiologicalCosts")

#read in physiology dataset
df <- read.csv("BSM_FemalePhys_wide.csv")

#summarize dataset to find out how many days I have for each bird, and then join with phys dataset
date_sum <- tbl %>%
  group_by(nest, site, tag) %>%
  summarize(min_date = min(Date), max_date = max(Date))

date_range <- tbl %>%
  group_by(tag) %>%
  summarize(min_date = min(Date), max_date = max(Date)) %>%
  inner_join(df, by = c("tag" = "GPS_tag"))


#### FILTER DATA ####
#Filter out only the dates for D9, and D10
library(lubridate)
library(dplyr)

#format hatch date
df$Hatch_1 <- mdy(df$Hatch_1)
tbl$Date <- ymd(tbl$`FIX-date`)

#calculate date when chicks were 9 days old
df$Day9 <- df$Hatch_1 + 9

#calculate date when chicks were 10 days old
df$Day10 <- df$Hatch_1 + 10

#filter GPS dataset using chick ages from physiology dataset
GPS_age_9 <- df %>% 
  select(Band, GPS_tag, Day9, Day10) %>%
  filter(!is.na(GPS_tag)) %>%
  left_join(tbl, by = c('Day9' = 'Date',
                        'GPS_tag' = 'tag'))

GPS_age <- df %>% 
  select(Band, GPS_tag, Day9, Day10) %>%
  filter(!is.na(GPS_tag)) %>%
  left_join(tbl, by = c('Day10' = 'Date',
                        'GPS_tag' = 'tag')) %>%
  bind_rows(GPS_age_9)

GPS_age %>%
  group_by(Band) %>%
  summarize(n = length(Latitude))

#save data file with filtered data
#write.csv(GPS_age, "filteredGPS19-20.csv")


#### CALCULATE RANGE SIZE ####
GPS_age <- read.csv("filteredGPS19-20.csv")


#Select only relevant columns
GPS <- GPS_age %>%
  select(Band, GPS_tag, Day9, Day10, year, site, nest, 
         `RTC-date`, `RTC-time`, Latitude, Longitude, HDOP)


#Remove rows with NAs because otherwise the MCP function won't work
GPS <- GPS[!is.na(GPS$Latitude) & !is.na(GPS$Longitude),]

#Only include ID and coordinates (MCP can't take any other info)
library(dplyr)
GPS_mcp <- dplyr::select(GPS, Band, Latitude, Longitude)


#################### All birds ################
#Create spatial points data frame by defining coordinates
library(sp)
library(rgdal)

GPS_mcp <- data.frame(GPS_mcp)

#Indicate which columns give my coordinates
coordinates(GPS_mcp) <- c("Longitude", "Latitude")

#Indicate what format the coordinates are in
proj4string(GPS_mcp) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

#Transform to easting and northing coordinates
gpsEN <- spTransform(GPS_mcp, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))

#Calculate MCPs for each ID (each bird)
library(adehabitatHR)

#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
GPSMCPs <- mcp(GPS_mcp, percent = 100) 

gpsMCPs <- mcp(gpsEN, percent = 100, unout = "m2")

forranges100 <- cbind(gpsMCPs$id, gpsMCPs$area)
colnames(forranges100) <- c("band", "foraging_area100")
forranges100 <- data.frame(forranges100)


#Look at MCP by percentile--this shows you have the MCP calculation changes
#depending on the percentage of points included
par("mar")
par(mar=c(1,1,1,1))
areas <- mcp.area(gpsEN, percent = seq(50, 100, by = 5))

#Link up with physiology dataset


#### MAPS ####
#Make spatial polygon maps
library(ggmap)

# get basemap
mybasemap <- get_stamenmap(bbox = c(left = min(GPS@coords[,1])-0.005, 
                                    bottom = min(GPS@coords[,2])-0.005, 
                                    right = max(GPS@coords[,1])+0.005, 
                                    top = max(GPS@coords[,2])+0.005), 
                           maptype = "toner",
                           zoom = 12)

#Turn the spatial data frame of points into just an R data frame for plotting in ggmap
GPSdf <- data.frame(GPS@coords, 
                           id = GPS@data$band)

#Create the map with foraging range polygons! 
mymap.hr <- ggmap(mybasemap) + 
  geom_polygon(data = fortify(GPSMCPs),  
               # Polygon layer needs to be "fortified" to add geometry to the dataframe
               aes(long, lat, colour = id, fill = id),
               alpha = 0.3) + 
  geom_point(data = GPSdf, 
             aes(x = Longitude, y = Latitude, colour = id, shape = id),
             size = 4)  +
  scale_shape_manual(values = c(0,1,2,
                                0,1,2,
                                0,1,2,
                                0,1,2,
                                0,1,2,
                                0,1,2,
                                0))

pdf("RangeSize_After_bw.pdf", width = 6.5, height = 8)
mymap.hr
dev.off()

