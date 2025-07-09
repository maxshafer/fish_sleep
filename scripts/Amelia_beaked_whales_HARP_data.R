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


# Load in and clean HARP data -------------------------------------------------------


#this is new passive acoustic monitoring data on rare beaked whale species from dryad
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

files <- lapply(filenames, readMat)
#make a dataframe out of just the time column from each mat file, all other information is irrelevant
files <- lapply(files, function(x) as.data.frame(x$MTT))
#append all mat files together into one large dataframe of observations
test <- do.call(rbind,files)

#now we have a df with 184,081 rows, each with the timepoint of a detection in MATLAB datetime format
colnames(test) <- c("MATLAB_datetime")

#first convert the times into the R format 
#refer to https://stackoverflow.com/questions/30072063/how-to-extract-the-time-using-r-from-a-matlab-serial-date-number
test$R_datetime <- lapply(test$MATLAB_datetime, function(x) (x - 719529)*86400)

#extract the datetime
times <- lapply(test$R_datetime, function(x) as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))

#the data includes entries for every second there was a detection, so detections spanning multiple seconds have many entries
#we will filter these out so we only have detections for each minute

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
#because the format is weird 
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
png(paste0("C:/Users/ameli/Downloads/beaked_whales_HARP_data/", whalename, "raw_data.png"), width=20,height=10,units="cm",res=200)
ggplot(test, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$hour), max(test$hour), by = 1),1))
dev.off()

#HARP stations were located at Heezen Canyon–HZ (271 days), Bear Seamount–BR (23 days), Nantucket Canyon–NC (145 days), Norfolk Canyon–NFC (424 days)
#At highest latitude (HZ) sunrise is at 7:15, sunset 4:30, at lowest latitude (NFC) sunrise is at 7:15, sunset at 5:00 (for January 1st)

#most detections occur in January
ggplot(test, aes(x = month)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(test$month), max(test$month), by = 1),1))


# Divide detections into their diel periods and plot ----------------------


#divide into day, night, dawn and dusk based on average timepoints
test$daynight <- NA

for(i in 1:nrow(test)){
  if(between(test[i, "hour"], 5, 7)){
    test[i, "daynight"] <- "Dawn"
  }
  else if(between(test[i, "hour"], 7, 18)){
    test[i, "daynight"] <- "Day"
  }
  else if(between(test[i, "hour"], 18, 20)){
    test[i, "daynight"] <- "Dusk"
  }
  else if(between(test[i, "hour"], 20, 24)){
    test[i, "daynight"] <- "Night"
  }
  else if(between(test[i, "hour"], 0, 5)){
    test[i, "daynight"] <- "Night"
  }
}

#we can plot the number of detections per hour, coloured by diel period (dawn, dusk, night, day), facet wrap by month (or season?)
ggplot(test, aes(x = hour, colour = daynight)) + geom_histogram(fill = "white", alpha = 0.5, position = "identity")

#how many detections are in each diel period
test2 <- test %>% count(daynight)
#I think we need to normalize for the number of hours in a diel period
#total calls in daytime/number of hours in daytime = average number of total calls per hour in daytime
#otherwise there will always be more total daytime calls than dusk or dawn calls because they're only about an hour long
test2$totalHours <- c(2, 11, 2, 9)
test2$normalized_count <- test2$n/test2$totalHours
norm_plot <- ggplot(test2, aes(y = normalized_count, x = daynight)) + geom_point() + geom_boxplot(fill = "transparent") 

#save out
png(paste0("C:/Users/ameli/Downloads/beaked_whales_HARP_data/", whalename, "normalized_diel_data.png"), width=20,height=10,units="cm",res=200)
norm_plot
dev.off()

#total number of detections per hour
test2 <- test %>% count(hour)
clean_detections <- ggplot(test2, aes(y = n, x = hour)) + geom_point() + geom_line(linewidth = 2) + geom_vline(xintercept = c(5, 7, 18, 20), linetype = "dotted", colour = c("blue", "deepskyblue", "hotpink", "red"), size = 1.5)

#save out
png(paste0("C:/Users/ameli/Downloads/beaked_whales_HARP_data/", whalename, "detections_clean.png"), width=20,height=10,units="cm",res=200)
clean_detections
dev.off()

#Methods from Munger et al, 2008  https://doi.org/10.1111/j.1748-7692.2008.00219.x
#We used a Kruskal-Wallis test (Zar 1999) to rank and compare mean-adjusted calling rates across diel periods,
#and a Tukey-Kramer multiple comparison test to determine if any of the diel periods showed significant difference in mean rank of mean-adjusted calling rates. 