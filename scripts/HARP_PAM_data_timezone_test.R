#install.packages("suncalc")
library(suncalc)

#first try with lowest latitude and longitude (37-09.991, 74-27.996),  NFC station
sun <- getSunlightTimes(as.Date(2016-06-26), lat = 37.1, lon = 74.3, keep = c("nauticalDawn", "sunriseEnd", "sunsetStart", "nauticalDusk"), tz = "UTC")

#Dawn = from when nautical dawn begins to when sunrise occurs, so from "nauticalDawn" to "sunriseEnd"
#Day = from when the sun is fully above the horizon to when it touches the horizon, so from "sunriseEnd" to "sunsetStart
#Dusk = from when sun drops below horizon to when nautical twilight ends, so from "sunsetStart" to "nauticalDusk"
#Night = from when astronomical twilight starts to when nautical dawn begins, so from "nauticalDusk" to "nauticalDawn"

for(i in 1:nrow(test10)){
  sun_times <- getSunlightTimes(test10[i, "date"], lat = 37.1, lon = 74.3, keep = c("nauticalDawn", "sunriseEnd", "sunsetStart", "nauticalDusk"), tz = "UTC")
  cbind(sun_times, test10, by = "date")
}

test10 <- test[sample(1:nrow(test),100),]
#test10 <- test

test10<- test10 %>% separate(value, into = c("date", "time"), sep = " ", remove = FALSE)

#extract for just one row to see if it works
# sun_times <- getSunlightTimes(as.Date(test10[1, "date"]$date), lat = 37.1, lon = 74.3, keep = c("nauticalDawn", "sunriseEnd", "sunsetStart", "nauticalDusk"), tz = "UTC")
# huh <- merge(test1, sun_times, by = "date")
# huh$daynight <- NA

#to confirm, on January 18th in HZ station, sunrise is at 6:44am, sunset at 16:30pm, 
#nautical dawn starts at 5:40, nautical dusk ends at 17:30pm
#according to the sunlightTimes function, sunrise ends at 21:18, sunset starts at  7:10
#nautical dawn starts at 20:15, nautical dusk is at 8:13
#it doesn't match up but this is in universal time zone so maybe that's why?

#find the time of night, day, dusk and dawn for all dates in the dataframe
#first make an empty list
datalist = vector("list", length = nrow(test10))

#for each row find the sun position times and then append them to the list
for(i in 1:nrow(test10)){
  dat <- getSunlightTimes(as.Date(test10[i, "date"]), lat = 37.1, lon = 74.3, keep = c("nauticalDawn", "sunriseEnd", "sunsetStart", "nauticalDusk"))
  datalist[[i]] <- dat
}

#bind the lists into a dataframe
times_final = do.call(rbind, datalist)

#merge the sun position times dataframe with the original call activity dataframe
#it doesn't matter if some rows have the same date since they will therefore have the same sunrise/sunset times regardless
times_final$identifier <- c(1:nrow(times_final))
test10$identifier <- c(1:nrow(test10))
huh <- merge(test10, times_final, by = "identifier")

#we want just the time to compare to the time of the calls
huh <- huh %>% separate(nauticalDawn, into = c("dawndate", "nauticalDawn"), sep = " ")
huh <- huh %>% separate(nauticalDusk, into = c("duskdate", "nauticalDusk"), sep = " ")
huh <- huh %>% separate(sunsetStart, into = c("sunsetdate", "sunsetStart"), sep = " ")
huh <- huh %>% separate(sunriseEnd, into = c("sunrisedate", "sunriseEnd"), sep = " ")

huh$daynight <- NA

#drop all unnecessary columns
huh <- huh[, c("date.x", "time", "year", "month", "day", "hour", "minute", "second", "nauticalDawn", "sunriseEnd", "sunsetStart", "nauticalDusk", "daynight")]

#then iterate through the dataframe and compare the time of the call to see if it falls between the timepoints for dawn, day, dusk or night

#this isn't working properly so troubleshoot

for(i in 1:nrow(huh)){
  if(huh[i, "time"] > huh[i, "nauticalDawn"] & huh[i, "time"] < huh[i,"sunriseEnd"]){
    huh[i, "daynight"] <- "Dawn"
  }
  if(huh[i, "time"] > huh[i,"sunriseEnd"] & huh[i, "time"] < huh[i,"sunsetStart"]){
    huh[i, "daynight"] <- "Day"
  }
  if(huh[i, "time"] > huh[i, "sunsetStart"] & huh[i, "time"] < huh[i, "nauticalDusk"]){
    huh[i, "daynight"] <- "Dusk"
  }
  if(huh[i, "time"] > huh[i, "nauticalDusk"] & huh[i, "time"] <= 24){
    huh[i, "daynight"] <- "Night"
  }
  if(huh[i, "time"] >= 0 & huh[i, "time"] < huh[i, "nauticalDusk"]){
    huh[i, "daynight"] <- "Night"
  }
}

#fake numbers while the other numbers are being strange
for(i in 1:nrow(huh)){
  if(huh[i, "hour"] >= 5 & huh[i, "hour"] < 6){
    huh[i, "daynight"] <- "Dawn"
  }
  else if(huh[i, "hour"] >= 6 & huh[i, "hour"] < 19){
    huh[i, "daynight"] <- "Day"
  }
  else if(huh[i, "hour"] >= 19 & huh[i, "hour"] < 20){
    huh[i, "daynight"] <- "Dusk"
  }
  else if(huh[i, "hour"] > huh[i, "nauticalDusk"] & huh[i, "hour"] <= 24){
    huh[i, "daynight"] <- "Night"
  }
  else if(huh[i, "hour"] >= 0 & huh[i, "hour"] < 5){
    huh[i, "daynight"] <- "Night"
  }
}

#now we have a dataframe with accurate information about the time of day each detection was made

#we can plot the number of detections per hour, coloured by diel period (dawn, dusk, night, day), facet wrap by month (or season?)
ggplot(huh, aes(x = hour, colour = daynight)) + geom_histogram(fill = "white", alpha = 0.5, position = "identity")

ggplot(huh, aes(y= hour, x = daynight)) + geom_point(alpha = 0.5) + geom_jitter() + geom_boxplot(fill = "transparent")

ggplot(huh, aes(x = as.factor(daynight))) + geom_histogram()

#how many detections are in each diel period
huh2 <- huh %>% count(daynight)
huh2 <- huh %>% count(hour)
ggplot(huh2, aes(y = n, x = daynight)) + geom_point() + geom_boxplot(fill = "transparent") 
ggplot(huh2, aes(y = n, x = hour)) + geom_point() 

#Methods from Munger et al, 2008  https://doi.org/10.1111/j.1748-7692.2008.00219.x
#We computed hourly calling rates within each of these diel periods and subtracted the overall calling rate (calls/hour) for that 24-h day to obtain a mean-adjusted calling rate in order to correct for variation in the total number of calls detected each day. 
#We used a Kruskal-Wallis test (Zar 1999) to rank and compare mean-adjusted calling rates across diel periods, and a Tukey-Kramer multiple comparison test to determine which, if any of the diel periods showed significant difference in mean rank of mean-adjusted right whale calling rates. 

#perform the one-way ANOVA

res.aov <- aov(n ~ daynight, data = huh2)
summary(res.aov)
