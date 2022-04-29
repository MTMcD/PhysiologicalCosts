#############################################################################
# Molly McDermott
# created 12/3/21
# edited 4/27/22
# putting together GPS data from raw tag data
#############################################################################

library(tidyverse)
library(dplyr)

setwd("~/SafranLab/BSM_Physiology/gps_data")

# create function for reading in csv and retaining file name
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

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
  filter(HDOP < 5)

#reset working directory
setwd("~/SafranLab/PhysiologicalCosts")

#read in physiology dataset
df <- read.csv("BSM_FemalePhys_wide.csv")
date_range$Hatch_1 <- as.Date(date_range$Hatch_1, "%m/%d/%Y")

#summarize dataset to find out how many days I have for each bird, and then join with phys dataset
date_sum <- tbl %>%
  group_by(nest, site, tag) %>%
  summarize(min_date = min(Date), max_date = max(Date))



date_range <- tbl %>%
  group_by(tag) %>%
  summarize(min_date = min(Date), max_date = max(Date)) %>%
  inner_join(df, by = c("tag" = "GPS_tag"))

write.csv(date_range, "BSM_FemalePhys_withGPS.csv")


difftime(date_range$min_date, date_range$Hatch_1)


