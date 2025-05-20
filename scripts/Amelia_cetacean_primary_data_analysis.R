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
#install.packages("deeptime")
library(deeptime)
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

setwd(here())
source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_cetacean_primary_data_analysis.R")

# # Section 1: Hyperoodon planifrons data ---------------------------------

#data from https://doi.org/10.1111/mms.12216 
hyper <- read_xlsx("C:/Users/ameli/OneDrive/Documents/cetacean_echo_data/Trickey_2015_Hyperoodon_planifrons.xlsx")

hyper <- hyper %>% separate("Encounter ID Date/time (GMT) Latitude Longitude Signal count", into = c("Encounter_ID", "Date", "Month", "Year", "Time", "Latitude", "Longitude", "Signal_count"), sep = " ")
hyper <- hyper %>% separate("Time", into = c("Start_time", "End_time"), sep = "–")
hyper$Signal_count <- as.numeric(hyper$Signal_count)

#since we are south and west, transform the coordinates to be negative 
hyper$Latitude <- substr(hyper$Latitude, start= 1, stop = 2)
hyper$Longitude <- substr(hyper$Longitude, start= 1, stop = 2)
hyper <- transform(hyper, Latitude = as.numeric(Latitude)*(-1))
hyper <- transform(hyper, Longitude = as.numeric(Longitude)*(-1))

#function to determine timezone from coordinates
hyper$timezone <- tz_lookup_coords(hyper$Latitude, hyper$Longitude, method = "fast")

#this function takes the geographic coordinates in a matrix c(lon, lat) with each coordinate as its own row
#requires date and timezone (get from lutz function)
#won't let me use the timezone column so just pick most common timezone (check if gives different results)
sunset <- sunriset(matrix(c(hyper$Longitude, hyper$Latitude), nrow = nrow(hyper)), as.POSIXct(paste0("2014-02", "-", hyper$Date), tz = "Etc/GMT+3"), direction = 'sunset', POSIXct.out = TRUE)
hyper <- cbind(hyper, sunset)
hyper <- hyper %>% separate("time", into = c("sunset_Date", "sunset_Time"), sep = " ")
hyper$sunset_Time <- substr(hyper$sunset_Time, start= 1, stop = 5)
hyper$sunset_Time <- str_replace(hyper$sunset_Time, pattern = ":", replacement = "")
hyper$sunset_Time <- as.numeric(hyper$sunset_Time)

sunrise <- sunriset(matrix(c(hyper$Longitude, hyper$Latitude), nrow = nrow(hyper)), as.POSIXct(paste0("2014-02", "-", hyper$Date), tz = "Etc/GMT+3"), direction = 'sunrise', POSIXct.out = TRUE)
colnames(sunrise) <- c("sunrise_day_frac", "sunrise_time")
hyper <- cbind(hyper, sunrise)
hyper <- hyper %>% separate("sunrise_time", into = c("sunrise_Date", "sunrise_Time"), sep = " ")
hyper$sunrise_Time <- substr(hyper$sunrise_Time, start= 1, stop = 5)
hyper$sunrise_Time <- str_replace(hyper$sunrise_Time, pattern = ":", replacement = "")
hyper$sunrise_Time <- as.numeric(hyper$sunrise_Time)
hyper <- hyper[, -c(1, 11, 12, 14, 15)]

ggplot(hyper, aes(x = Start_time, y = Signal_count)) + geom_bar(stat = "identity")

#study took place in the anarctic (South Orkney Islands, South Shetland Islands, and Antarctic Peninsula)
#during february 19-23 2014. Sunrise was at about 21:00pm and sunset was at 12:00pm 
#https://www.timeanddate.com/sun/@-60.446,-51.697?month=2&year=2014 

hyper$diel_period <- "Night"

for(i in 1:nrow(hyper)){
  if(hyper[i, "Start_time"] %in% hyper[i, "sunrise_Time"]:hyper[i, "sunset_Time"]){
    hyper[i, "diel_period"] <- "Day"
  }
}

ggplot(hyper, aes(x = Start_time, y = Signal_count, fill = diel_period)) + geom_bar(stat = "identity")

ggplot(hyper, aes(x = diel_period, y = Signal_count, fill = diel_period)) + geom_boxplot() + geom_jitter()

#replace zeros with a very small value, doesn't change the result
hyper[hyper == 0] <- 0.001

hyper_ANOVA <- aov(Signal_count ~ diel_period, data = hyper)
summary(hyper_ANOVA)

#including dawn and dusk timepoints
#set dusk and dawn to be the hour before and after sunset/sunrise
for(i in 1:nrow(hyper)){
  if(hyper[i, "Start_time"] %in% (hyper[i, "sunset_Time"]-100):(hyper[i, "sunset_Time"]+100)){
    hyper[i, "diel_period"] <- "Dusk"
  }
  if(hyper[i, "Start_time"] %in% (hyper[i, "sunrise_Time"]-100):(hyper[i, "sunrise_Time"]+100)){
    hyper[i, "diel_period"] <- "Dawn"
  }
}

#one bout of calls at dusk
ggplot(hyper, aes(x = Start_time, y = Signal_count, fill = diel_period)) + geom_bar(stat = "identity")
ggplot(hyper, aes(x = diel_period, y = Signal_count, fill = diel_period)) + geom_boxplot() + geom_jitter()

#no significant associations between diel period and signal count
hyper_ANOVA <- aov(Signal_count ~ diel_period, data = hyper)
summary(hyper_ANOVA)

#new hyperoodon data from Barlow et al, 2021
#https://doi.org/10.1016/j.dsr2.2021.104973

hyper <- read_xlsx("C:/Users/ameli/OneDrive/Documents/cetacean_echo_data/Barlow_2021_Hyperoodon_planifrons.xlsx")
hyper$`Event sequential number Event ID Start date/time (UTC) Event type Number of echolocation signals South latitude West longitude` <- str_replace(hyper$`Event sequential number Event ID Start date/time (UTC) Event type Number of echolocation signals South latitude West longitude`, pattern = "Southern bottlenose whale", replacement = "Southern_bottlenose_whale")

hyper <- hyper %>% separate("Event sequential number Event ID Start date/time (UTC) Event type Number of echolocation signals South latitude West longitude", into = c("Event", "sequential_number", "Event_ID", "Start_date_time", "Event_type", "Signal_count", "South_latitude", "West_longitude"), sep = " ")

#remove anything that isn't a southern bottlenose whale 
hyper <- filter(hyper, Event_type == "Southern_bottlenose_whale")
hyper$Start_date_time <- str_replace(hyper$Start_date_time, pattern = ":", replacement = "")
hyper$Signal_count <- as.numeric(hyper$Signal_count)

#since we are south and west, transform the coordinates to be negative 
hyper <- transform(hyper, South_latitude = as.numeric(South_latitude)*(-1))
hyper <- transform(hyper, West_longitude = as.numeric(West_longitude)*(-1))

#make columns for the dates so it can be interpreted by the sunriset function
hyper <- hyper %>% separate(Event_ID, into = c("Month", "Day", "Year"), sep = "/")

#function to determine timezone from coordinates
hyper$timezone <- tz_lookup_coords(hyper$South_latitude, hyper$West_longitude, method = "fast")

#this function takes the geographic coordinates in a matrix c(lon, lat) with each coordinate as its own row
#requires date and timezone (get from lutz function)
sunset <- sunriset(matrix(c(hyper$West_longitude, hyper$South_latitude), nrow = nrow(hyper)), as.POSIXct(paste0(hyper$Year, "-", hyper$Month, "-", hyper$Day), tz = "Atlantic/South_Georgia"), direction = 'sunset', POSIXct.out = TRUE)
hyper <- cbind(hyper, sunset)
hyper <- hyper %>% separate("time", into = c("sunset_Date", "sunset_Time"), sep = " ")
hyper$sunset_Time <- substr(hyper$sunset_Time, start= 1, stop = 5)
hyper$sunset_Time <- str_replace(hyper$sunset_Time, pattern = ":", replacement = "")
hyper$sunset_Time <- as.numeric(hyper$sunset_Time)

sunrise <- sunriset(matrix(c(hyper$West_longitude, hyper$South_latitude), nrow = nrow(hyper)), as.POSIXct(paste0(hyper$Year, "-", hyper$Month, "-", hyper$Day), tz = "Atlantic/South_Georgia"), direction = 'sunrise', POSIXct.out = TRUE)
colnames(sunrise) <- c("sunrise_day_frac", "sunrise_time")
hyper <- cbind(hyper, sunrise)
hyper <- hyper %>% separate("sunrise_time", into = c("sunrise_Date", "sunrise_Time"), sep = " ")
hyper$sunrise_Time <- substr(hyper$sunrise_Time, start= 1, stop = 5)
hyper$sunrise_Time <- str_replace(hyper$sunrise_Time, pattern = ":", replacement = "")
hyper$sunrise_Time <- as.numeric(hyper$sunrise_Time)

hyper <- hyper[, -c(1, 2, 7, 12, 13, 15, 16)]

ggplot(hyper, aes(x = Start_date_time, y = Signal_count)) + geom_bar(stat = "identity")

#study took place in  Falkland Islands to the South Sandwich Islands and South Georgia from December 30, 2019 to January 29, 2020 (Leg 1)
#only one detection occurred on the second leg of the trip from King George Island to Puerto Williams, Chile via the Antarctic Peninsula from February 11 to 27, 2020 (Leg 2)

#off the Falkland islands in January sunrise was at 500h and sunset at 2130
#https://www.timeanddate.com/sun/@-55,-47?month=1&year=2020

hyper$diel_period <- "Night"

for(i in 1:nrow(hyper)){
  if(hyper[i, "Start_date_time"] %in% hyper[i, "sunrise_Time"]:hyper[i, "sunset_Time"]){
    hyper[i, "diel_period"] <- "Day"
  }
}

ggplot(hyper, aes(x = Start_date_time, y = Signal_count, fill = diel_period)) + geom_bar(stat = "identity")

ggplot(hyper, aes(x = diel_period, y = Signal_count, fill = diel_period)) + geom_boxplot() + geom_jitter()

#replace zeros with a very small value, doesn't change the result
hyper[hyper == 0] <- 0.001

hyper_ANOVA <- aov(Signal_count ~ diel_period, data = hyper)
summary(hyper_ANOVA)

#including dawn and dusk timepoints
#set dusk and dawn to be the hour before and after sunset/sunrise
for(i in 1:nrow(hyper)){
  if(hyper[i, "Start_date_time"] %in% (hyper[i, "sunset_Time"]-100):(hyper[i, "sunset_Time"]+100)){
    hyper[i, "diel_period"] <- "Dusk"
  }
  if(hyper[i, "Start_date_time"] %in% (hyper[i, "sunrise_Time"]-100):(hyper[i, "sunrise_Time"]+100)){
    hyper[i, "diel_period"] <- "Dawn"
  }
}

ggplot(hyper, aes(x = Start_date_time, y = Signal_count, fill = diel_period)) + geom_bar(stat = "identity")

ggplot(hyper, aes(x = diel_period, y = Signal_count, fill = diel_period)) + geom_boxplot() + geom_jitter()

#replace zeros with a very small value, doesn't change the result
hyper[hyper == 0] <- 0.001

hyper_ANOVA <- aov(Signal_count ~ diel_period, data = hyper)
summary(hyper_ANOVA)


# Section 2: vaquita ---------------------------------------------------------

#vaquita moment
vaquita <- read.csv("C:/Users/ameli/Downloads/vaquita_PAM.csv")
vaquita$datetime <- parse_date(vaquita$Start)
vaquita$hour <- hour(vaquita$datetime)
vaquita <- vaquita[order(vaquita$hour, decreasing = TRUE), ]
View(vaquita)
vaquita %>% group_by(hour, Date) %>%
  ggplot(., aes(x = Date, y = hour)) + geom_jitter(alpha = 0.5)  + geom_hline(yintercept = 5.75, color = "red", size = 1) + geom_hline(yintercept = 19.5, color = "red", size = 1)
ggplot(vaquita, aes(x = hour, fill = Date)) + geom_bar() + geom_vline(xintercept = 5.75, color = "red", size = 2) + geom_vline(xintercept = 19.5, color = "red", size = 2)
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
#red lines indicate onset of dawn and dusk
ggplot(vaquita, aes(x = hour, fill = diel)) + geom_bar() + geom_vline(xintercept = 5.75, color = "red", size = 1) + geom_vline(xintercept = 19.5, color = "red", size = 1)
ggplot(vaquita, aes(x = diel)) + geom_bar()

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


# Section 3: Camera trap ungulates -------------------------------------------


#Camera trap data for impala, kudu and wildebeest
impala <- read.csv("C:/Users/ameli/Downloads/camtrapHawkes-v.2.0.0/camtrapHawkes/data/camtrap_data/data.csv")
View(impala)

impala <- impala %>% filter(snapshotName == "impala")
impala <- impala %>% separate(eventTime, c("hour", "minute", "second"))

ggplot(impala, aes(x = hour)) + geom_bar()

kudu <- read.csv("C:/Users/ameli/Downloads/camtrapHawkes-v.2.0.0/camtrapHawkes/data/camtrap_data/data.csv")
View(kudu)

kudu <- kudu %>% filter(snapshotName == "kudu")
kudu <- kudu %>% separate(eventTime, c("hour", "minute", "second"))

ggplot(kudu, aes(x = hour)) + geom_bar()

blue <- read.csv("C:/Users/ameli/Downloads/camtrapHawkes-v.2.0.0/camtrapHawkes/data/camtrap_data/data.csv")
View(blue)

blue <- blue %>% filter(snapshotName == "wildebeestblue")
blue <- blue %>% separate(eventTime, c("hour", "minute", "second"))

ggplot(blue, aes(x = hour)) + geom_bar()


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

