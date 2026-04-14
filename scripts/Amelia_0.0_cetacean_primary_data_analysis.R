# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(suncalc)
#install.packages("deeptime")
#library(deeptime)
#update.packages("ggplot2")
library(ggplot2)
setwd(here())
#to find sunset sunrise times
#install.packages("suntools")
library(suntools)
#to find timezones based on coordinates
#install.packages("lutz")
library(lutz)
library(parsedate)
#install.packages("R.matlab")
library(R.matlab)

setwd(here())
source("scripts/fish_sleep_functions.R")

# # Section 1: Hyperoodon planifrons data ---------------------------------

#data from https://doi.org/10.1111/mms.12216 
hyper <- read_xlsx("C:/Users/ameli/OneDrive/Documents/cetacean_echo_data/Trickey_2015_Hyperoodon_planifrons.xlsx")

hyper <- hyper %>% separate("Encounter ID Date/time (GMT) Latitude Longitude Signal count", into = c("Encounter_ID", "Date", "Month", "Year", "Time", "Latitude", "Longitude", "Signal_count"), sep = " ")
hyper <- hyper %>% separate("Time", into = c("Start_time", "End_time"), sep = "–")

hyper$Signal_count <- as.numeric(hyper$Signal_count)
hyper$Start_time <- as.numeric(hyper$Start_time)
hyper$End_time <- as.numeric(hyper$End_time)

hyper <- hyper %>% mutate(Start_time = Start_time/100) %>% 
  mutate(Start_hour = as.integer(Start_time), Start_min = ((Start_time - as.integer(Start_time)) * 100/60)) %>%
  mutate(Start_time = Start_hour + Start_min)

hyper <- hyper %>% mutate(date = paste("2014", "02", Date, sep = "-")) 
hyper$date <- as.Date(parse_date_time(hyper$date, orders = "ymd"))

#LOESS (locally estimated scatterplot smoothing) which does not require you to describe a model
#It takes small subsets of the data along the independent variable and makes many models (usually first or second degree polynomials) and joins them together.
#from https://andrewirwin.github.io/data-visualization/working-models.html
ggplot(hyper, aes(y = Signal_count, x = Start_time)) + geom_point() + geom_smooth(method = "loess", formula = "y~x")

#since we are south and west, transform the coordinates to be negative 
hyper$lat <- substr(hyper$Latitude, start= 1, stop = 2)
hyper$lon <- substr(hyper$Longitude, start= 1, stop = 2)
hyper <- transform(hyper, lat = as.numeric(lat)*(-1))
hyper <- transform(hyper, lon = as.numeric(lon)*(-1))

#function to determine timezone from coordinates
hyper$timezone <- tz_lookup_coords(hyper$lat, hyper$lon, method = "accurate")

#test for first row. Sunrise at 9:21 pm (21:21), sunset at 12:18 (12:28) pm
#returns 4:23 sunrise and 19:19 sunset
getSunlightTimes(date = as.Date(parse_date_time("2014-2-19", orders = "ymd")), 
                 lat = -60, lon = -54, tz = tz_lookup_coords(-60, -54, method = "accurate"))

#find out the coordinates of each location and use to calculate sunrise and set times
#use nautical dawn as the start of dawn and sunrise end as the end (encompasses dawn) 
#use sunset start as the start of dusk and nautical dusk as the end (encompasses dusk) 

sun_times <- getSunlightTimes(data = hyper[, c("date", "lon", "lat")], 
                              keep = c("sunriseEnd", "sunsetStart", "nauticalDawn", "nauticalDusk"),
                              tz = "Atlantic/South_Georgia")

sun_times <- sun_times %>% separate_wider_delim(cols = c(4:7), delim = " ", names_sep = "_") %>% select(date, nauticalDawn_2, nauticalDusk_2, sunriseEnd_2 ,sunsetStart_2)
sun_times <- sun_times %>% separate_wider_delim(cols = c(2:5), delim = ":", names_sep = "_") %>% select(date, nauticalDawn_2_1, nauticalDawn_2_2, nauticalDusk_2_1, nauticalDusk_2_2, sunriseEnd_2_1, sunriseEnd_2_2, sunsetStart_2_1, sunsetStart_2_2)
sun_times <- sun_times %>% mutate(dawn_start = as.numeric(nauticalDawn_2_1) + as.numeric(nauticalDawn_2_2)/60) %>% 
  mutate(dawn_end = as.numeric(sunriseEnd_2_1) + as.numeric(sunriseEnd_2_2)/60) %>% 
  mutate(dusk_start = as.numeric(sunsetStart_2_1) + as.numeric(sunsetStart_2_2)/60) %>% 
  mutate(dusk_end = as.numeric(nauticalDusk_2_1) + as.numeric(nauticalDusk_2_2)/60) %>% select(date, dusk_start, dusk_end, dawn_start, dawn_end)

hyper <- cbind(hyper, sun_times) %>% select("Date", "Month", "Year", "Start_time", "Signal_count", "date", "dusk_start", "dusk_end", "dawn_start", "dawn_end")
#study took place in the anarctic (South Orkney Islands, South Shetland Islands, and Antarctic Peninsula)

#with suncalc times
ggplot(hyper, aes(y = Signal_count, x = Start_time)) + 
  theme_minimal() +
  annotate(geom = "rect", xmin = mean(hyper$dawn_start), xmax = mean(hyper$dawn_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = mean(hyper$dusk_start), xmax = mean(hyper$dusk_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = mean(hyper$dawn_start), ymin = -Inf, ymax = Inf, fill = "grey70") +
  annotate(geom = "rect", xmin = mean(hyper$dusk_end), xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_point() + 
  scale_x_continuous(limits = c(7, 24), breaks = c(0:24)) + 
  geom_smooth(method = "loess", formula = "y~x", colour = "black")

#with times put in manually
#during february 19-23 2014. Sunrise at 9:21 pm (21:21), sunset at 12:18 (12:18) pm
#nautical dawn - 7:28pm (19:28)
#nautical dusk- 2:12pm (14:12) 
#https://www.timeanddate.com/sun/@-60.446,-51.697?month=2&year=2014 

ggplot(hyper, aes(y = Signal_count, x = Start_time)) + 
  theme_minimal() +
  annotate(geom = "rect", xmin = 19.5, xmax = 22.3, ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 12.3, xmax = 14.2, ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 14.2, xmax = 19.5, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_point() + 
  scale_x_continuous(limits = c(7, 24), breaks = c(0:24)) + 
  geom_smooth(method = "loess", formula = "y~x", colour = "black")

hyper1 <- hyper

###new hyperoodon data from Barlow et al, 2021
#https://doi.org/10.1016/j.dsr2.2021.104973

hyper <- read_xlsx("C:/Users/ameli/OneDrive/Documents/cetacean_echo_data/Barlow_2021_Hyperoodon_planifrons.xlsx")
hyper$`Event sequential number Event ID Start date/time (UTC) Event type Number of echolocation signals South latitude West longitude` <- str_replace(hyper$`Event sequential number Event ID Start date/time (UTC) Event type Number of echolocation signals South latitude West longitude`, pattern = "Southern bottlenose whale", replacement = "Southern_bottlenose_whale")

hyper <- hyper %>% separate("Event sequential number Event ID Start date/time (UTC) Event type Number of echolocation signals South latitude West longitude", into = c("Event", "sequential_number", "Event_ID", "Start_date_time", "Event_type", "Signal_count", "South_latitude", "West_longitude"), sep = " ")

#remove anything that isn't a southern bottlenose whale 
hyper <- filter(hyper, Event_type == "Southern_bottlenose_whale")
hyper$Start_date_time <- str_replace(hyper$Start_date_time, pattern = ":", replacement = "")
hyper$Signal_count <- as.numeric(hyper$Signal_count)

#since we are south and west, transform the coordinates to be negative 
hyper <- mutate(hyper, lat = as.numeric(South_latitude)*(-1))
hyper <- mutate(hyper, lon = as.numeric(West_longitude)*(-1))

#make columns for the dates so it can be interpreted by the sunriset function
hyper <- hyper %>% separate(Event_ID, into = c("Month", "Day", "Year"), sep = "/")
hyper$date <- as.Date(as.POSIXct(paste0(hyper$Year, "-", hyper$Month, "-", hyper$Day)))
     
hyper$Start_date_time <- as.numeric(hyper$Start_date_time)      
hyper <- hyper %>% mutate(Start_time = Start_date_time/100) %>% 
  mutate(Start_hour = as.integer(Start_time), Start_min = ((Start_time - as.integer(Start_time)) * 100/60)) %>%
  mutate(Start_time = Start_hour + Start_min)

#function to determine timezone from coordinates
hyper$timezone <- tz_lookup_coords(hyper$lat, hyper$lon, method = "fast")

#test run
#off the Falkland islands (-55, -47) on January 1 sunrise was at 4:33am and sunset at 9:48pm (21:48)
#https://www.timeanddate.com/sun/@-55,-47?month=1&year=2020
#suncalc returns 3:34 sunrise and 20:50 sunset. Not perfect but close
getSunlightTimes(date = as.Date(parse_date_time("2020-01-01", orders = "ymd")), 
                 lat = -55, lon = -47, tz = tz_lookup_coords(-54, -47, method = "accurate"))

#find out the coordinates of each location and use to calculate sunrise and set times
#use nautical dawn as the start of dawn and sunrise end as the end (encompasses dawn) 
#use sunset start as the start of dusk and nautical dusk as the end (encompasses dusk) 

sun_times <- getSunlightTimes(data = hyper[, c("date", "lon", "lat")], 
                              keep = c("sunriseEnd", "sunsetStart", "nauticalDawn", "nauticalDusk"),
                              tz = "Atlantic/South_Georgia")

sun_times <- sun_times %>% separate_wider_delim(cols = c(4:7), delim = " ", names_sep = "_") %>% select(date, nauticalDawn_2, nauticalDusk_2, sunriseEnd_2 ,sunsetStart_2)
sun_times <- sun_times %>% separate_wider_delim(cols = c(2:5), delim = ":", names_sep = "_") %>% select(date, nauticalDawn_2_1, nauticalDawn_2_2, nauticalDusk_2_1, nauticalDusk_2_2, sunriseEnd_2_1, sunriseEnd_2_2, sunsetStart_2_1, sunsetStart_2_2)
sun_times <- sun_times %>% mutate(dawn_start = as.numeric(nauticalDawn_2_1) + as.numeric(nauticalDawn_2_2)/60) %>% 
  mutate(dawn_end = as.numeric(sunriseEnd_2_1) + as.numeric(sunriseEnd_2_2)/60) %>% 
  mutate(dusk_start = as.numeric(sunsetStart_2_1) + as.numeric(sunsetStart_2_2)/60) %>% 
  mutate(dusk_end = as.numeric(nauticalDusk_2_1) + as.numeric(nauticalDusk_2_2)/60) %>% select(dusk_start, dusk_end, dawn_start, dawn_end)

#replace 0 with 24 so mean dusk end is calculated properly later
sun_times$dusk_end[is.na(sun_times$dusk_end)] <- 1
sun_times$dusk_end[as.integer(sun_times$dusk_end) == 0] <- sun_times$dusk_end[as.integer(sun_times$dusk_end) == 0] + 24
sun_times$dusk_end[sun_times$dusk_end == 1] <- NA

hyper <- cbind(hyper, sun_times) %>% select(c(Signal_count, lat, lon, date, Start_time, timezone, dusk_start, dusk_end, dawn_start, dawn_end))

ggplot(hyper, aes(y = Signal_count, x = Start_time)) + 
  theme_minimal() +
  annotate(geom = "rect", xmin = mean(hyper$dawn_start, na.rm = TRUE), xmax = mean(hyper$dawn_end, na.rm = TRUE), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = mean(hyper$dusk_start, na.rm = TRUE), xmax = mean(hyper$dusk_end, na.rm = TRUE), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = mean(hyper$dawn_start, na.rm = TRUE), ymin = -Inf, ymax = Inf, fill = "grey70") +
  annotate(geom = "rect", xmin = mean(hyper$dusk_end, na.rm = TRUE), xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_point() + 
  scale_x_continuous(breaks = c(0:24)) +
  geom_smooth(method = "loess", formula = "y~x", colour = "black")

#study took place in  Falkland Islands to the South Sandwich Islands and South Georgia from December 30, 2019 to January 29, 2020 (Leg 1)
#only one detection occurred on the second leg of the trip from King George Island to Puerto Williams, Chile via the Antarctic Peninsula from February 11 to 27, 2020 (Leg 2)


hyper_merged <- rbind(hyper[, c("Start_time", "Signal_count", "dusk_start", "dusk_end", "dawn_start", "dawn_end")], hyper1[, c("Start_time", "Signal_count", "dusk_start", "dusk_end", "dawn_start", "dawn_end")])
ggplot(hyper_merged, aes(y = Signal_count, x = Start_time)) + 
  theme_minimal() +
  annotate(geom = "rect", xmin = mean(hyper_merged$dawn_start, na.rm = TRUE), xmax = mean(hyper_merged$dawn_end, na.rm = TRUE), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = mean(hyper_merged$dusk_start, na.rm = TRUE), xmax = mean(hyper_merged$dusk_end, na.rm = TRUE), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = mean(hyper_merged$dawn_start, na.rm = TRUE), ymin = -Inf, ymax = Inf, fill = "grey70") +
  annotate(geom = "rect", xmin = mean(hyper_merged$dusk_end, na.rm = TRUE), xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_point() +
  scale_x_continuous(breaks = c(0:24)) + 
  geom_smooth(method = "loess", formula = "y~x", colour = "black")


# Section 2: vaquita ---------------------------------------------------------
#study: https://iucn-csg.org/wp-content/uploads/2023/06/Vaquita-Survey-2023-Main-Report.pdf
#data from appendix: https://iucn-csg.org/wp-content/uploads/2023/06/Vaquita-Survey-2023-Appendices-FINAL.pdf

#vaquita moment
vaquita <- read.csv("C:/Users/ameli/Downloads/vaquita_PAM.csv")
vaquita$datetime <- parse_date(vaquita$Start)
vaquita$hour <- hour(vaquita$datetime)
vaquita <- vaquita[order(vaquita$hour, decreasing = TRUE), ]

#red lines indicate onset of dawn and dusk
ggplot(vaquita, aes(x = hour, fill = Date)) + geom_bar() + geom_vline(xintercept = 5.75, color = "red", size = 2) + geom_vline(xintercept = 19.5, color = "red", size = 2)

ggplot(vaquita, aes(x = hour)) + geom_density(size = 2) + geom_vline(xintercept = 5.75, color = "red", size = 2) + geom_vline(xintercept = 19.5, color = "red", size = 2)

#all data is taken from May in San Felipe, Mexico so sunrise times are about 5:45am and sunset is around 19:30pm
vaquita$diel <- "day"
for(i in 1:60)
  if(vaquita[i, "hour"] > 19.5){
    vaquita[i, "diel"] <- "night"
  } 

for(i in 1:60)
  if(vaquita[i, "hour"] < 5.75){
    vaquita[i, "diel"] <- "night"
  } 

#civil twilight lasts for about an 45 minutes in san felipe 
#so dawn is from 5am-5:45 and dusk is 7:30-8:15 (19:30-20:15)

for(i in 1:60)
  if(vaquita[i, "hour"] >= 19.5 & vaquita[i, "hour"] <= 20.25){
    vaquita[i, "diel"] <- "dusk"
  } 

for(i in 1:60)
  if(vaquita[i, "hour"] <= 5.75 & vaquita[i, "hour"] >= 5){
    vaquita[i, "diel"] <- "dawn"
  } 

ggplot(vaquita, aes(x = hour, fill = diel)) + geom_bar() + geom_vline(xintercept = 5.75, color = "red", size = 1) + geom_vline(xintercept = 19.5, color = "red", size = 1)

ggplot(vaquita, aes(x = hour)) +   
  theme_minimal() +
  annotate(geom = "rect", xmin = 5, xmax = 5.75, ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 19.5, xmax = 20.25, ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf, fill = "grey") +
  annotate(geom = "rect", xmin = 20.25, xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey") +
  geom_density(size = 2) 
   

# Section 3: Camera trap ungulates -------------------------------------------
#data from https://doi.org/10.1111/jzo.70062
#For this analysis, we focused on six camera trap grids in the savanna biome in 
#northern South Africa: the Associated Private Nature Reserves (around Kruger National Park),
#Kruger National Park, Madikwe Game Reserve, Pilanesberg National Park, Somkhanda Game Reserve and Venetia Limpopo Nature Reserve 

#Camera trap data for impala, kudu and wildebeest
camera_trap.df <- read.csv("C:/Users/ameli/Downloads/camtrapHawkes-v.2.0.0/camtrapHawkes/data/camtrap_data/data.csv")
camera_trap.df$eventDateTime <- paste(camera_trap.df$eventDate, camera_trap.df$eventTime, sep = " ")
camera_trap.df$eventDate <- parse_date_time(camera_trap.df$eventDate, orders = "ymd")
camera_trap.df$eventDateTime <- as_datetime(camera_trap.df$eventDateTime)
camera_trap.df$hour <- hour(camera_trap.df$eventDateTime)
#convert minute to fraction of hour (out of 60)
camera_trap.df$min <- minute(camera_trap.df$eventDateTime)/60
camera_trap.df$hourmin <- as.numeric(camera_trap.df$hour) + as.numeric(camera_trap.df$min)

#plot number of detections per hour and minute
ggplot(camera_trap.df, aes(x = hourmin)) + geom_density(size = 2)

#separate by location (six different parks) 
ggplot(camera_trap.df, aes(x = hourmin)) + geom_density(size = 2) + facet_wrap(~locationID)

#sites are labelled A-F so its undetermined which site is which park
#therefore we can't use the specific coordinates for each so use coordinates roughly equidistant from all (-25, 30)

#use suncalc to get sunrise and sunset times
#test run: sunrise should be 5:02, sunset at 6:32. Returns 5:03 sunrise and 6:33 sunset!
getSunlightTimes(date = as.Date(parse_date_time("2019-11-26", orders = "ymd")), lat = -25, lon = 30, tz = tz_lookup_coords(-25, 30))

#find out the coordinates of each location and use to calculate sunrise and set times
#use nautical dawn as the start of dawn and sunrise end as the end (encompasses dawn) -about an hour
#use sunset start as the start of dusk and nautical dusk as the end (encompasses dusk) -also about an hour
sun_times <- getSunlightTimes(date = as.Date(camera_trap.df$eventDate), 
                              lat = -25, lon = 30, 
                              keep = c("sunriseEnd", "sunsetStart", "nauticalDawn", "nauticalDusk"),
                              tz = tz_lookup_coords(-25, 30))

sun_times <- sun_times %>% separate_wider_delim(cols = c(4:7), delim = " ", names_sep = "_") %>% select(date, nauticalDawn_2, nauticalDusk_2, sunriseEnd_2 ,sunsetStart_2)
sun_times <- sun_times %>% separate_wider_delim(cols = c(2:5), delim = ":", names_sep = "_") %>% select(date, nauticalDawn_2_1, nauticalDawn_2_2, nauticalDusk_2_1, nauticalDusk_2_2, sunriseEnd_2_1, sunriseEnd_2_2, sunsetStart_2_1, sunsetStart_2_2)
sun_times <- sun_times %>% mutate(dawn_start = as.numeric(nauticalDawn_2_1) + as.numeric(nauticalDawn_2_2)/60) %>% 
  mutate(dawn_end = as.numeric(sunriseEnd_2_1) + as.numeric(sunriseEnd_2_2)/60) %>% 
  mutate(dusk_start = as.numeric(sunsetStart_2_1) + as.numeric(sunsetStart_2_2)/60) %>% 
  mutate(dusk_end = as.numeric(nauticalDusk_2_1) + as.numeric(nauticalDusk_2_2)/60) %>% select(date, dusk_start, dusk_end, dawn_start, dawn_end)
  
camera_trap.df <- cbind(camera_trap.df, sun_times)

#sunset and sunrise times seem to vary most from month/season?
camera_trap.df$month <- month(camera_trap.df$eventDate)

#sun times vary by only 1.7 hours across the year, so taking the mean should be a good approximation
camera_trap.df %>% summarize(dawn_end = max(dawn_end) - min(dawn_end),  dawn_start = max(dawn_start) - min(dawn_start),
                             dusk_end = max(dusk_start) - min(dusk_start), dusk_start = max(dusk_start) - min(dusk_start))

#Impala Aepyceros melampus

impala <- camera_trap.df %>% filter(snapshotName == "impala") 
ggplot(impala, aes(x = hourmin)) +
  theme_minimal() + 
  annotate(geom = "rect", xmin = mean(impala$dawn_start), xmax = mean(impala$dawn_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = mean(impala$dusk_start), xmax = mean(impala$dusk_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = mean(impala$dawn_start), ymin = -Inf, ymax = Inf, fill = "grey70") +
  annotate(geom = "rect", xmin = mean(impala$dusk_end), xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_density(size = 2) + ggtitle("Camera trap detections of Aepyceros melampus")

#Greater kudu Tragelaphus strepsiceros
kudu <- camera_trap.df %>% filter(snapshotName == "kudu") 
ggplot(kudu, aes(x = hourmin)) +
  theme_minimal() + 
  annotate(geom = "rect", xmin = mean(kudu$dawn_start), xmax = mean(kudu$dawn_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = mean(kudu$dusk_start), xmax = mean(kudu$dusk_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = mean(kudu$dawn_start), ymin = -Inf, ymax = Inf, fill = "grey70") +
  annotate(geom = "rect", xmin = mean(kudu$dusk_end), xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_density(size = 2) + ggtitle("Camera trap detections of Tragelaphus strepsiceros")

#Blue wildebeest Connochaetes taurinus
blue <- camera_trap.df %>% filter(snapshotName == "wildebeestblue") 
ggplot(blue, aes(x = hourmin)) +
  theme_minimal() + 
  annotate(geom = "rect", xmin = mean(blue$dawn_start), xmax = mean(blue$dawn_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = mean(blue$dusk_start), xmax = mean(blue$dusk_end), ymin = -Inf, ymax = Inf, fill = "pink") +
  annotate(geom = "rect", xmin = 0, xmax = mean(blue$dawn_start), ymin = -Inf, ymax = Inf, fill = "grey70") +
  annotate(geom = "rect", xmin = mean(blue$dusk_end), xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey70") +
  geom_density(size = 2) + ggtitle("Camera trap detections of Connochaetes taurinus")

# Section 4: Narwhal ---------------------------------------------------------

#### Narwhals https://www.science.org/doi/10.1126/sciadv.ade0440?adobe_mc=MCMID%3D53649406453315412110550155571971043555%7CMCORGID%3D242B6472541199F70A4C98A6%2540AdobeOrg%7CTS%3D1695155886#supplementary-materials

narwhale <- read.csv("C:/Users/ameli/Downloads/doi_10.5061_dryad.8gtht76tq__v5/Data_Buzz.txt", sep = "\t")

head(narwhale)
View(narwhale)
dim(narwhale)

## Parse time into datetime format
narwhale$datetime <- parse_date(narwhale$GPS_time)
narwhale$hour <- hour(narwhale$datetime)
narwhale$minute <- minute(narwhale$datetime)
narwhale$day <- day(narwhale$datetime)

## When do Buzz's occur
#this plots the number of seconds with and without a buzz in each hour
buzzes <- narwhale %>% group_by(hour, Buzz) %>% count()

#plotting out the number of seconds without a buzz (top plot) and with a buzz (bottom plot)
ggplot(buzzes, aes(x = hour, y = n, group = Buzz)) + geom_point() + geom_line() + facet_wrap(~Buzz, scales = "free", ncol = 1)

#plot the number of detection positive minutes (more standard)
returnMinPerHour <- function(hour = 1){
  buzzes <- narwhale %>% filter(hour == hour) %>% group_by(minute, Buzz) %>% count()
  return(buzzes)
}

buzz_list <- lapply(unique(narwhale$hour), function(x) returnMinPerHour(hour = x))
buzzes <- do.call(rbind.data.frame, buzz_list)
buzzes <- buzzes %>% filter(Buzz == 1)
hour_list <- lapply(0:23, function(x) rep(x, times = 60))
#buzzes$hour <- hour_list


buzzes$total_minute <- 0:(nrow(buzzes)-1)

ggplot(buzzes, aes(x = total_minute, y = n)) + geom_point() + geom_line()

#detection positive minutes per hour


## Plot depth versus time

## plot 1 Hz data (1 timepoint per second)
ggplot(narwhale[1:10000, ], aes(x = datetime, y = Depth)) + geom_point()

## Can we average by minute, or hour, or day?

#AVERAGE BY...
day_depth <- narwhale %>% group_by(day) %>% summarise(mean_depth = mean(Depth))

hour_depth <- narwhale %>% group_by(hour) %>% summarise(mean_depth = mean(Depth))
hour_depth

dayhour_depth <- narwhale %>% group_by(Ind, day, hour) %>% summarise(mean_depth = mean(Depth))

ggplot(dayhour_depth, aes(x = interaction(hour,day), y = mean_depth, colour = Ind, group = Ind)) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~Ind, scales = "free", ncol = 1)


minute_depth <- narwhale %>% group_by(day, hour, minute) %>% summarise(mean_depth = mean(Depth))


ggplot(minute_depth, aes(x = interaction(day,hour,minute), y = mean_depth)) + geom_point() + theme(axis.text.x = element_blank())

## Plot both buzzes and depth?
## take mean of "Buzz" for each day/hour/min, you will get a range from 0-1
## Find a way to ask if there was a buzz or not in each day/hour/min?

mean_buzz_day <- narwhale %>% group_by(day, Ind) %>% summarise(mean_buzz = mean(Buzz))
mean_buzz_day
#creates a tibble with three rows: the day, the mean frequency of a buzz occuring (0=no buzz, 1=buzz))
#and the individual whale making the buzz

ggplot(mean_buzz_day, aes(x = day, y = mean_buzz, colour = Ind, group = Ind)) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~Ind, scales = "free", ncol = 1)
#creates a graph showing the mean frequency of buzzes occurring per day by individual

mean_buzz_dayhour <- narwhale %>% group_by(day, hour, Ind) %>% summarise(mean_buzz = mean(Buzz))
mean_buzz_dayhour
#same as above by breaks down the mean number of buzzes by hour as well as day and individual

ggplot(mean_buzz_dayhour, aes(x = interaction(hour, day), y = mean_buzz, colour = Ind, group = Ind)) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~Ind, scales = "free", ncol = 1)

#Same as before but adding the minute as well
mean_buzz_dayhourminute <- narwhale %>% group_by(day, hour, minute, Ind) %>% summarise(mean_buzz = mean(Buzz))
mean_buzz_dayhourminute

#plot
ggplot(mean_buzz_dayhourminute, aes(x = interaction(hour, day, minute), y = mean_buzz, colour = Ind, group = Ind)) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~Ind, scales = "free", ncol = 1)


#plot depth over time, noting when the mean buzzes are occurring (ie hunting)
depth_buzz <- narwhale %>% group_by(hour, Ind) %>% summarise(mean_buzz = mean(Buzz), mean_depth = mean(Depth))
depth_buzz

ggplot(depth_buzz, aes(x = interaction(hour, Ind), y = mean_depth, colour = mean_buzz, size = mean_buzz, group = Ind)) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~Ind, scales = "free")



# Section 5: Beaked whale HARP data ----------------------------------------------------

#this is passive acoustic monitoring data on rare beaked whale species from dryad
#https://doi.org/10.5061/dryad.gf1vhhmw0

#MTT: Time of event as Matlab datenumber (days elapsed since January 0, 0000). Each row represents one detection.
#we need to convert this into a readable datetime from the matlab format into the r format
#then plot the number of detections per hour (or half hour)

#to get all observations, load in the other mat files 
#for mesoplodon mirus
# filenames <- list.files("C:/Users/ameli/Downloads/Mm", pattern = "*.mat", full.names = TRUE)
# whalename <- "mesoplodon_mirus"

# # #for mesoplodon europaeus
# filenames <- list.files("C:/Users/ameli/Downloads/Me", pattern = "*.mat", full.names = TRUE)
# whalename <- "mesoplodon_europaeus"

#for kogiia
filenames <- list.files("C:/Users/ameli/Downloads/Kogia", pattern = "*.mat", full.names = TRUE)
whalename <- "kogia"

#each file is named after the site it was recorded at and the disk number (tends to be multiple per site)

files <- lapply(filenames, readMat)
#make a dataframe out of just the time column from each mat file, all other information is irrelevant
files <- lapply(files, function(x) as.data.frame(x$MTT))

#extract the location name from the file name
filenames <- str_remove(filenames, "C:/Users/ameli/Downloads/")
#add in metadata from file names
names(files) <- filenames

#append all mat files together into one large dataframe of observations
test <- do.call(rbind,files)

#now we have a df with 184,081 rows, each with the timepoint of a detection in MATLAB datetime format
test$metadata <- row.names(test)
colnames(test) <- c("MATLAB_datetime", "metadata")


#MTT: Time of event as Matlab datenumber (days elapsed since January 0, 0000).
#refer to https://stackoverflow.com/questions/30072063/how-to-extract-the-time-using-r-from-a-matlab-serial-date-number
Matlab2Rdate <- function(val) as.Date(val - 1, origin = '0000-01-01') 
Matlab2Rdate(733038.6)
"2006-12-27"
Matlab2Rdate(735147.4)
"2012-10-05"

(735147.4 - 719529)*86400
#test datetime conversion, this doesn't seem right is this data really from 2012?
as.POSIXct(1349427015.16854, origin = "1970-01-01", tz = "UTC")
"2012-10-05 08:50:15 UTC"

#first convert the times into the R format 
#test$R_datetime <- lapply(test$MATLAB_datetime, function(x) (x - 719529)*86400)
test <- test %>% mutate(R_datetime = (MATLAB_datetime - 719529)*86400)

#extract the datetime
test <- test %>% mutate(times = as.POSIXct(R_datetime, origin = "1970-01-01", tz = "UTC"))

#separate out into component parts 
test$year <- year(test$times)
test$month <- month(test$times)
test$day <- day(test$times)
test$hour <- hour(test$times)
test$minute <- minute(test$times)
test$second <- second(test$times)

#plot the total detections in each hour bin
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))

ggplot(test, aes(x = hour)) + geom_density() 

#for each site add the latitude and longitude (manually)
#HH site: 25 (N), -85 (W)
#BS site: 30 (N), -77 (W)
#

#find sunrise and sunset times using suncalc


#use suncalc to get sunrise and sunset times
#site is the Golf of Mexico, Howell Hook. Lat = 25N, long = 85W
#test run: sunrise should be x, sunset at x. Returns 5:35 sunrise and 5:23pm (17:23) sunset!
getSunlightTimes(date = as.Date(parse_date_time("2012-10-05", orders = "ymd")), lat = 25, lon = -85, tz = tz_lookup_coords(25, -85))

#find out the coordinates of each location and use to calculate sunrise and set times
#use nautical dawn as the start of dawn and sunrise end as the end (encompasses dawn) -about an hour
#use sunset start as the start of dusk and nautical dusk as the end (encompasses dusk) -also about an hour
sun_times <- getSunlightTimes(date = as.Date(camera_trap.df$eventDate), 
                              lat = 25, lon = -85, 
                              keep = c("sunriseEnd", "sunsetStart", "nauticalDawn", "nauticalDusk"),
                              tz = tz_lookup_coords(25, -85))

sun_times <- sun_times %>% separate_wider_delim(cols = c(4:7), delim = " ", names_sep = "_") %>% select(date, nauticalDawn_2, nauticalDusk_2, sunriseEnd_2 ,sunsetStart_2)
sun_times <- sun_times %>% separate_wider_delim(cols = c(2:5), delim = ":", names_sep = "_") %>% select(date, nauticalDawn_2_1, nauticalDawn_2_2, nauticalDusk_2_1, nauticalDusk_2_2, sunriseEnd_2_1, sunriseEnd_2_2, sunsetStart_2_1, sunsetStart_2_2)
sun_times <- sun_times %>% mutate(dawn_start = as.numeric(nauticalDawn_2_1) + as.numeric(nauticalDawn_2_2)/60) %>% 
  mutate(dawn_end = as.numeric(sunriseEnd_2_1) + as.numeric(sunriseEnd_2_2)/60) %>% 
  mutate(dusk_start = as.numeric(sunsetStart_2_1) + as.numeric(sunsetStart_2_2)/60) %>% 
  mutate(dusk_end = as.numeric(nauticalDusk_2_1) + as.numeric(nauticalDusk_2_2)/60) %>% select(date, dusk_start, dusk_end, dawn_start, dawn_end)

camera_trap.df <- cbind(camera_trap.df, sun_times)
