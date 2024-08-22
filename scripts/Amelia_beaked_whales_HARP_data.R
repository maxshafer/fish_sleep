# #import packages --------------------------------------------------------
library(phangorn)
library(stringr)
library(here)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(dplyr)
library(readxl)
library(tidyr)
library(parsedate)
library(lubridate)
library(ggplot2)
library(forcats)
library(ggforce)
#install.packages("ggstream")
library(ggstream)
library(magrittr)
#library(tidyverse)
#install.packages("R.matlab")
library(R.matlab)

setwd(here())

#this is new passive acoustic monitoring data on rare beaked whale species from dryad
#https://doi.org/10.5061/dryad.gf1vhhmw0


# # Section 1: Mesoplodon mirus -------------------------------------------

mirus <- readMat("C:/Users/ameli/Downloads/Mm/WAT_NC_01_disk13_Mm_TPWS1.mat")

#MTT: Time of event as Matlab datenumber (days elapsed since January 0, 0000). Each row represents one detection.
#we need to convert this into a readable datetime from the matlab format into the r format
#then plot the number of detections per hour (or half hour)

#to get all observations, load in the other mat files 
filenames <- list.files("C:/Users/ameli/Downloads/Mm", pattern = "*.mat", full.names = TRUE)
files <- lapply(filenames, readMat)
#make a dataframe out of just the time column from each mat file, all other information is irrelevant
files <- lapply(files, function(x) as.data.frame(x$MTT))
#append all mat files together into one large dataframe of observations
test <- do.call(rbind,files)

#now we have a df with 184,081 rows, each with the timepoint of a detection in MATLAB datetime format
colnames(test) <- c("MATLAB_datetime")

#for some reason its very difficult to keep the values converted into R datetime values
#this is the only way that worked

#first convert the times into the R format 
#refer to https://stackoverflow.com/questions/30072063/how-to-extract-the-time-using-r-from-a-matlab-serial-date-number
test$R_datetime <- lapply(test$MATLAB_datetime, function(x) (x - 719529)*86400)

#extract the datetime
times <- lapply(test$R_datetime, function(x) as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))

#the data includes entries for every second there was a detection, so detections spanning multiple seconds have many entries
#should we filter these out so we only have detections for each minute? Instead of each second?
#or detections every 5 minutes? 

times <- as.data.frame(times)
#split into two dataframes so there's enough memory to pivot longer
times1 <- times[, 1:(length(times)%/% 2)]
times2 <- times[, (length(times1)+1):length(times)]
#pivot longer, we'll now have two dataframes of about equal length
test1 <- pivot_longer(times1, cols = 1:length(times1))
test2 <- pivot_longer(times2, cols = 1:length(times2))

#separate out into component parts (too difficult to remove duplicate datetimes since the seconds go on for 13 decimal points)
test1$year <- lapply(test1$value, year)
test1$month <- lapply(test1$value, month)
test1$day <- lapply(test1$value, day)
test1$hour <- lapply(test1$value, hour)
test1$minute <- lapply(test1$value, minute)
test1$second <- lapply(test1$value, second)

#need to round seconds so we can remove duplicate observations from the same second 
test1$second <- lapply(test1$second, as.integer)

#remove duplicates (same day, hour, minute, second)
#this will bring us from 92040 observations to 22,787 observations(yay!)
test1 <- test1 %>% distinct(year, month, day, hour, minute, second, .keep_all = TRUE)

#do the same for the other half of the dataframe
test2$year <- lapply(test2$value, year)
test2$month <- lapply(test2$value, month)
test2$day <- lapply(test2$value, day)
test2$hour <- lapply(test2$value, hour)
test2$minute <- lapply(test2$value, minute)
test2$second <- lapply(test2$value, second)
test2$second <- lapply(test2$second, as.integer)

#this will bring us from 92,041 observations to 20,412 observations(yay!)
test2 <- test2 %>% distinct(year, month, day, hour, minute, second, .keep_all = TRUE)

#now that the dataframe is less large, we should have enough memory to bind them back together
test <- rbind(test1, test2)
#because the format is weird lmao
test <- as.data.frame(test)

#each row in the hour column is actually a 1x1 list containing the datetime so use unlist
test$hour <- unlist(test$hour)

#plot the total detections in each hour bin
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))

#can colour by month to determine if there's any major differences according to season
test$month <- unlist(test$month)
test$month_factor <- factor(test$month)
ggplot(test, aes(x = hour, fill = month_factor, color = month_factor)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))

#save out
png("C:/Users/ameli/Downloads/beaked_whales_HARP_data/mesoplodon_mirus.png", width=20,height=10,units="cm",res=200)
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))
dev.off()

#we can statistically test if the detections are higher at night vs day using a T-test
#first, Mesoplodon mirus detections occurred mainly in the Western North Atlantic  (between North Carolina and Massachusetts)
#HARP stations were located at Heezen Canyon–HZ (271 days), Bear Seamount–BR (23 days), Nantucket Canyon–NC (145 days), Norfolk Canyon–NFC (424 days)
#At highest latitude (HZ) sunrise is at 7:15, sunset 4:30, at lowest latitude (NFC) sunrise is at 7:15, sunset at 5:00 (for January 1st)

#most detections occur in January
ggplot(test, aes(x = month)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$month), max(test$month), by = 1),1))

#at Bear Seamount in January sunrise occurs at approximately 700h and sunset at 1700h
#do we want to pick a time where the day and nighttime are approximately even?
# i should do the day/night by season thing :'(

#for now: set sunrise to 700h and sunset to 1900h
#subset to smaller df
test <- test[sample(1:nrow(test),100),]
test$daynight <- NA

#categorize each timepoint as belonging to the night or the day class
for(i in 1:nrow(test)){
  if(test[i, "hour"] > 7 & test[i, "hour"] < 20){
    test[i, "daynight"]  <- "Day"
  }
  else{ 
    test[i, "daynight"] <- "Night"
  }
}

ggplot(test, aes(x = hour, fill = daynight, color = daynight)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))

ggplot(test, aes(x = hour, colour = daynight)) + geom_histogram(fill = "white", alpha = 0.5, position = "identity")

ggplot(test, aes(y = hour, x = daynight)) + geom_point() + geom_boxplot(fill = "transparent")

#we want to vizualize the number of detections per hour in day vs night
tally <- test %>% count(hour)

for(i in 1:nrow(tally)){
  if(tally[i, "hour"] > 7 & tally[i, "hour"] < 20){
    tally[i, "daynight"]  <- "Day"
  }
  else{ 
    tally[i, "daynight"] <- "Night"
  }
}

#this plots the number of detections in a given hour, divided by day and night
ggplot(tally, aes(y = n, x = daynight, color = daynight)) + geom_point() + geom_boxplot(fill = "transparent")

#right now it appears that there are more detections during nighttime hours
#we can test if this difference is statistically significant
#null hypothesis: there is no difference between detections at night or at day
#alternative hypothesis: there is a significant difference between

ggplot(tally, aes(y = n, x = daynight, color = as.factor(hour))) + geom_point() 
#the highest detections are at 1400h, followed by 600h, then 1900h, 400h

#we want to know if the number of detections in nighttime hours is statistically different than the number of detections during daytime hours
#Welch Two Sample t-test: investigates if there is a significant difference between the mean of two independent groups that may have unequal variance. 
#The test compares the means of two groups while considering the variability within each group.

t.test(n ~ daynight, data = tally)
#says there is a significant difference in means but the p-value is 0.7881?

# # Section 2: Kogia -------------------------------------------
#to get all observations, load in the other mat files 
filenames <- list.files("C:/Users/ameli/Downloads/Kogia", pattern = "*.mat", full.names = TRUE)
files <- lapply(filenames, readMat)
#make a dataframe out of just the time column from each mat file, all other information is irrelevant
files <- lapply(files, function(x) as.data.frame(x$MTT))
#append all mat files together into one large dataframe of observations
test <- do.call(rbind, files)

#now we have a df with 726 rows, each with the timepoint of a detection in MATLAB datetime format
colnames(test) <- c("MATLAB_datetime")

#for some reason its very difficult to keep the values converted into R datetime values
#this is the only way that worked

#first convert the times into the R format 
#refer to https://stackoverflow.com/questions/30072063/how-to-extract-the-time-using-r-from-a-matlab-serial-date-number
test$R_datetime <- lapply(test$MATLAB_datetime, function(x) (x - 719529)*86400)

#extract the datetime
times <- lapply(test$R_datetime, function(x) as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))

#the data includes entries for every second there was a detection, so detections spanning multiple seconds have many entries
#should we filter these out so we only have detections for each minute? Instead of each second?
#or detections every 5 minutes? 

times <- as.data.frame(times)
#split into two dataframes so there's enough memory to pivot longer
times1 <- times[, 1:(length(times)%/% 2)]
times2 <- times[, (length(times1)+1):length(times)]
#pivot longer, we'll now have two dataframes of about equal length
test1 <- pivot_longer(times1, cols = 1:length(times1))
test2 <- pivot_longer(times2, cols = 1:length(times2))

#separate out into component parts (too difficult to remove duplicate datetimes since the seconds go on for 13 decimal points)
test1$year <- lapply(test1$value, year)
test1$month <- lapply(test1$value, month)
test1$day <- lapply(test1$value, day)
test1$hour <- lapply(test1$value, hour)
test1$minute <- lapply(test1$value, minute)
test1$second <- lapply(test1$value, second)

#need to round seconds so we can remove duplicate observations from the same second 
test1$second <- lapply(test1$second, as.integer)

#remove duplicates (same day, hour, minute, second)
test1 <- test1 %>% distinct(year, month, day, hour, minute, second, .keep_all = TRUE)

#do the same for the other half of the dataframe
test2$year <- lapply(test2$value, year)
test2$month <- lapply(test2$value, month)
test2$day <- lapply(test2$value, day)
test2$hour <- lapply(test2$value, hour)
test2$minute <- lapply(test2$value, minute)
test2$second <- lapply(test2$value, second)
test2$second <- lapply(test2$second, as.integer)

#remove duplicates
test2 <- test2 %>% distinct(year, month, day, hour, minute, second, .keep_all = TRUE)

#now that the dataframe is less large, we should have enough memory to bind them back together
test <- rbind(test1, test2)
#because the format is weird lmao
test <- as.data.frame(test)

#each row in the hour column is actually a 1x1 list containing the datetime so use unlist
test$hour <- unlist(test$hour)

#plot the total detections in each hour bin
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))

#save out
png("C:/Users/ameli/Downloads/beaked_whales_HARP_data/kogia.png", width=20,height=10,units="cm",res=200)
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))
dev.off()

# # Section 2: Mesoplodon europaeus -------------------------------------------
#to get all observations, load in the other mat files 
filenames <- list.files("C:/Users/ameli/Downloads/Me", pattern = "*.mat", full.names = TRUE)
files <- lapply(filenames, readMat)
#make a dataframe out of just the time column from each mat file, all other information is irrelevant
files <- lapply(files, function(x) as.data.frame(x$MTT))
#append all mat files together into one large dataframe of observations
test <- do.call(rbind, files)

#now we have a df with 1,293,807 rows, each with the timepoint of a detection in MATLAB datetime format
colnames(test) <- c("MATLAB_datetime")

#for some reason its very difficult to keep the values converted into R datetime values
#this is the only way that worked

#first convert the times into the R format 
#refer to https://stackoverflow.com/questions/30072063/how-to-extract-the-time-using-r-from-a-matlab-serial-date-number
test$R_datetime <- lapply(test$MATLAB_datetime, function(x) (x - 719529)*86400)

#extract the datetime
times <- lapply(test$R_datetime, function(x) as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))

#the data includes entries for every second there was a detection, so detections spanning multiple seconds have many entries
#should we filter these out so we only have detections for each minute? Instead of each second?
#or detections every 5 minutes? 

times <- as.data.frame(times)
#split into two dataframes so there's enough memory to pivot longer
times1 <- times[, 1:(length(times)%/% 2)]
times2 <- times[, (length(times1)+1):length(times)]
#pivot longer, we'll now have two dataframes of about equal length
test1 <- pivot_longer(times1, cols = 1:length(times1))
test2 <- pivot_longer(times2, cols = 1:length(times2))

#separate out into component parts (too difficult to remove duplicate datetimes since the seconds go on for 13 decimal points)
test1$year <- lapply(test1$value, year)
test1$month <- lapply(test1$value, month)
test1$day <- lapply(test1$value, day)
test1$hour <- lapply(test1$value, hour)
test1$minute <- lapply(test1$value, minute)
test1$second <- lapply(test1$value, second)

#need to round seconds so we can remove duplicate observations from the same second 
test1$second <- lapply(test1$second, as.integer)

#remove duplicates (same day, hour, minute, second), this takes us down to 243,598 rows
test1 <- test1 %>% distinct(year, month, day, hour, minute, second, .keep_all = TRUE)

#do the same for the other half of the dataframe
test2$year <- lapply(test2$value, year)
test2$month <- lapply(test2$value, month)
test2$day <- lapply(test2$value, day)
test2$hour <- lapply(test2$value, hour)
test2$minute <- lapply(test2$value, minute)
test2$second <- lapply(test2$value, second)
test2$second <- lapply(test2$second, as.integer)

#remove duplicates, this takes us down to 200,920 rows
test2 <- test2 %>% distinct(year, month, day, hour, minute, second, .keep_all = TRUE)

#now that the dataframe is less large, we should have enough memory to bind them back together
test <- rbind(test1, test2)
#because the format is weird lmao
test <- as.data.frame(test)

#each row in the hour column is actually a 1x1 list containing the datetime so use unlist
test$hour <- unlist(test$hour)

#plot the total detections in each hour bin
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))

#save out
png("C:/Users/ameli/Downloads/beaked_whales_HARP_data/Mesoplodon_europaeus.png", width=20,height=10,units="cm",res=200)
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))
dev.off()

