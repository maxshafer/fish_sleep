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

setwd(here())

# # # Phocoena sinus  -------------------------------------------------------

#import dataset from iNaturalist
#how to cite Phocoena sinus GBIF.org (02 August 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.wr59d3
#iNaturalist csv format doesn't open properly in excel so open from google sheets
url <- 'https://docs.google.com/spreadsheets/d/1o8fBjCxBdDgAxp9oSEa8GfjHotXQDUaWTiisN7s_3ew/edit?usp=sharing'
sinus <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
sinus <- sinus[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for only observations
sinus <- sinus %>% filter(basisOfRecord == "HUMAN_OBSERVATION") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
sinus$datetime <- parse_date(sinus$eventDate)
sinus$hour <- hour(sinus$datetime)

#filter out entries with a date but no time 
sinus <- sinus %>% filter(nchar(eventDate) > 10) 

#plot the occurrences by time
ggplot(sinus, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(sinus$hour), max(sinus$hour), by = 0.5),1))

#issue: these are mainly going to be daytime observations
#other issue: hard to determine the behaviour of the animal when the photo was taken (swimming, sleeping at surface etc)

#cite for dataset of all cetaceans
#GBIF.org (06 August 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.dap9kh


# # # Hyperoodontidae iNaturalist species ------------------------------------
# cite for hyperoodontidae GBIF.org (06 August 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.uenbs9
url <- 'https://docs.google.com/spreadsheets/d/1kzTjAyWu3OycmduUm0sX2wJw9e5GgDegRsRGPPsLB6g/edit?gid=351943248#gid=351943248'
hyper <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
hyper <- hyper[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for only observations
#hyper <- hyper %>% filter(basisOfRecord == "HUMAN_OBSERVATION") 
hyper <- hyper %>% filter(species != "")
hyper <- hyper %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
hyper$datetime <- parse_date(hyper$eventDate)
hyper$hour <- hour(hyper$datetime)

#need to remove NA values that are converted into 0 (these do not represent observations at midnight)
hyper <- hyper %>% filter(nchar(eventDate) > 10) 

#plot the occurrences by time, facet wrap by species
ggplot(hyper, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(hyper$hour), max(hyper$hour), by = 1),1)) + facet_wrap(~species)

#save out as individual plots
for(i in 1:21){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/hyperoodontidae", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(hyper, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(hyper$hour), max(hyper$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}

# # # Phocoenidae iNaturalist species ------------------------------------
# cite for hyperoodontidae GBIF.org (06 August 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.uenbs9
url <- ''
Phoco <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
Phoco <- Phoco[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for only observations
#hyper <- hyper %>% filter(basisOfRecord == "HUMAN_OBSERVATION") 
Phoco <- Phoco %>% filter(species != "")
Phoco <- Phoco %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
Phoco$datetime <- parse_date(Phoco$eventDate)
Phoco$hour <- hour(Phoco$datetime)

#need to remove NA values that are converted into 0 (these do not represent observations at midnight)
#find better method to not remove actual midnight observations
Phoco <- Phoco %>% filter(hour != "0") 

#plot the occurrences by time, facet wrap by species
ggplot(Phoco, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Phoco$hour), max(Phoco$hour), by = 1),1)) + facet_wrap(~species, ncol = 5, nrow = 4)

#save out as individual plots
for(i in 1:20){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/Phocoenidae", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(hyper, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(hyper$hour), max(hyper$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}

# # # Phocoena phocoena Duke Harbour iNaturalist species ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1EQ8Guf57wcYXQd0Th0xMxwBvJPCuRrjL8Mo-6UZ7lwo/edit?usp=sharing'
Phoco <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
Phoco <- Phoco[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for only observations
#hyper <- hyper %>% filter(basisOfRecord == "HUMAN_OBSERVATION") 
Phoco <- Phoco %>% filter(species != "")
Phoco <- Phoco %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
Phoco$datetime <- parse_date(Phoco$eventDate)
Phoco$hour <- hour(Phoco$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(Phoco, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Phoco$hour), max(Phoco$hour), by = 1),1)) 

#save out
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/Phocoena_Duke_Harbour.png")
ggplot(Phoco, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Phoco$hour), max(Phoco$hour), by = 1),1)) 
dev.off()

# # # Eschrichtius robustus Granite Canyon iNaturalist species ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/15zoFhIO3725B8rAQK0kVokX5hu_8J4Dn5yi5aus8pv4/edit?usp=sharing'
Robustus <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
Robustus <- Robustus[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter dataset to get rid of entries without a date or a species
Robustus <- Robustus %>% filter(species != "")
Robustus <- Robustus %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
Robustus$datetime <- parse_date(Robustus$eventDate)
Robustus$hour <- hour(Robustus$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(Robustus, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Robustus$hour), max(Robustus$hour), by = 1),1)) 

#save out
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/Eschrichtius_Robustus_Granite_Canyon.png")
ggplot(Robustus, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Robustus$hour), max(Robustus$hour), by = 1),1)) 
dev.off()

# # # Balaenoptera edeni Kauai iNaturalist species ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1FBAbO80TGHH5jtzS3O5Amf05sKG53y1GOPaDUFe7zJ0/edit?usp=sharing'
Bryde <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
Bryde <- Bryde[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter dataset to get rid of entries without a date or a species
Bryde <- Bryde %>% filter(species != "")
Bryde <- Bryde %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
Bryde$datetime <- parse_date(Bryde$eventDate)
Bryde$hour <- hour(Bryde$datetime)

#need to remove NA values that are converted into 0 (these do not represent observations at midnight)
#find better method to not remove actual midnight observations
Bryde <- Bryde %>% filter(hour != "0") 

#plot the occurrences by time, facet wrap by species
ggplot(Bryde, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Bryde$hour), max(Bryde$hour), by = 1),1)) 

#save out
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/Balaenoptera_edeni_kauai.png")
ggplot(Bryde, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(Bryde$hour), max(Bryde$hour), by = 1),1)) 
dev.off()

# # # DECAF iNaturalist data ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1-iIVPYGeLYB-kxSNaDqE-DpdaLhcDu2iOPjfkpk9gBo/edit?usp=sharing'
inat <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
inat <- inat[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for rows with species info and a date/time
inat <- inat %>% filter(species != "")
inat <- inat %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
inat$datetime <- parse_date(inat$eventDate)
inat$hour <- hour(inat$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 

#save out as individual plots
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/DECAF.png", width=20,height=10,units="cm",res=200)
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 
dev.off()

# # # minke kauai iNaturalist data ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1J_HctaHP6A4mKYaHm1DsO_LgJ7LEb2p46NsbekvgmIA/edit?usp=sharing'
inat1 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
inat1 <- inat1[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#because file is too large to import into google sheets I split into two sheets
#load in second csv file and combine
url <- 'https://docs.google.com/spreadsheets/d/1gPWXArfTjejhdtXTnmpMtazIjhtoLRZYxswgv6nJwPg/edit?usp=sharing'
inat2 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
inat2 <- inat2[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

inat <- rbind(inat1, inat2)

#filter for rows with species info and a date/time
inat <- inat %>% filter(species != "")
inat <- inat %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
inat$datetime <- parse_date(inat$eventDate)
inat$hour <- hour(inat$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 

#save out as individual plots
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/minke_kauai.png", width=20,height=10,units="cm",res=200)
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 
dev.off()
# # # antarctic iNaturalist data ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/10T1qLbi1Wwvr8QM-Idw8oWfVy5XLkulk6Thn9NlRKAo/edit?usp=sharing'
inat <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
inat <- inat[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#convert eventDate column into datetime format, extract hour
inat$datetime <- parse_date(inat$eventDate)
inat$hour <- hour(inat$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 

#save out as individual plots
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/antarctic.png", width=20,height=10,units="cm",res=200)
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 
dev.off()
# # # fin whale kauai iNaturalist data ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/129BviM5Ak0fBF4sPveXlKMkRqWEqzuSB20QuRTLT5Ps/edit?usp=sharing'
inat <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
inat <- inat[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for rows with species info and a date/time
inat <- inat %>% filter(species != "")
inat <- inat %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
inat$datetime <- parse_date(inat$eventDate)
inat$hour <- hour(inat$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 

#save out as individual plots
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/fin_whale_kauai.png", width=20,height=10,units="cm",res=200)
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 
dev.off()
# # # Lifewatch phocoena phocoena iNaturalist data ------------------------------------
url <- ''
inat1 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

url <- ''
inat2 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

url <- ''
inat3 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

inat <- rbind(inat1, inat2)
inat <- rbind(inat, inat3)

inat <- inat[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord")]

#filter for rows with species info and a date/time
inat <- inat %>% filter(species != "")
inat <- inat %>% filter(eventDate != "") 

#convert eventDate column into datetime format
#sinus$eventDate <- as.POSIXct(sinus$eventDate, format="%Y-%m-%d %H:%M:%S")
inat$datetime <- parse_date(inat$eventDate)
inat$hour <- hour(inat$datetime)

#plot the occurrences by time, facet wrap by species
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 

#save out as individual plots
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/lifewatch_phocoena_phocoena.png", width=20,height=10,units="cm",res=200)
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) 
dev.off()

# # # Human observations cetaceans iNaturalist species ------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1etslv-6q_lf1sw6d-ZKX__aAYCnBnPFsLnPwnLerT5w/edit?usp=sharing'
inat <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
inat <- inat[,c("datasetKey", "species", "countryCode", "eventDate", "basisOfRecord", "identifiedBy")]

#filter for only observations with species and date info
inat <- inat %>% filter(species != "")
inat <- inat %>% filter(eventDate != "") 

#convert eventDate column into datetime format
inat$datetime <- parse_date(inat$eventDate)
inat$hour <- hour(inat$datetime)

#filter out entries with a date but no time 
inat <- inat %>% filter(nchar(eventDate) > 10)

#view raw data for species with no or limited information
NA_species <- inat %>% filter(species %in% c("Berardius arnuxii", "Lagenorhynchus cruciger", "Mesoplodon bowdoini", "Mesoplodon carlhubbsi", "Mesoplodon eueu", "Mesoplodon grayi", "Mesoplodon hectori", "Mesoplodon layardii", "Mesoplodon mirus", "Mesoplodon perrini", "Mesoplodon peruvianus", "Mesoplodon traversii", "Phocoena dioptrica", "Phocoena spinipinnis", "Sousa sahulensis", "Tasmacetus shepherdi"))
ggplot(NA_species, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(NA_species$hour), max(NA_species$hour), by = 1),1)) + facet_wrap(~species)

#save out as individual plots
for(i in 1:11){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/NA_species_not_normalized", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(NA_species, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(NA_species$hour), max(NA_species$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}


low_data_species <- inat %>% filter(species %in% c("Balaenoptera omurai", "Caperea marginata", "Cephalorhynchus eutropia", "Hyperoodon planifrons", "Kogia breviceps", "Kogia sima", "Lagenodelphis hosei", "Lagenorhynchus australis", "Lissodelphis borealis", "Lissodelphis peronii", "Mesoplodon europaeus", "Mesoplodon ginkgodens", "Mesoplodon stejnegeri", "Orcaella heinsohni", "Phocoena sinus"))
ggplot(low_data_species, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(low_data_species$hour), max(low_data_species$hour), by = 1),1)) + facet_wrap(~species)


#filter for species with more than 50? occurrences 
inat_small <- inat %>% group_by(species) %>% filter(n() < 50)
inat_large <- inat %>% group_by(species) %>% filter(n() > 200)
inat_med <- inat %>% group_by(species) %>% filter(n() >50 & n() <200)

#plot the occurrences by time, facet wrap by species (more than 50 but less than 100 occurrences)
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/human_obs_over50.png", width=20,height=10,units="cm",res=200)
ggplot(inat_med, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat_med$hour), max(inat_med$hour), by = 1),1)) + facet_wrap(~species)
dev.off()

#plot the species with less than 50 occurrences 
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/human_obs_under50.png", width=20,height=10,units="cm",res=200)
ggplot(inat_small, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat_small$hour), max(inat_small$hour), by = 1),1)) + facet_wrap(~species)
dev.off()

#plot the species with over 200 occurrences 
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/human_obs_over200.png", width=20,height=10,units="cm",res=200)
ggplot(inat_large, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat_large$hour), max(inat_large$hour), by = 1),1)) + facet_wrap(~species)
dev.off()

#save out as individual plots
for(i in 1:21){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/human_obs_over50_species", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}

#save out as individual plots
for(i in 1:35){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/human_obs_under50_species", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(inat_small, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat_small$hour), max(inat_small$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}

#save out as individual plots
for(i in 1:27){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/human_obs_over100_species", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(inat_large, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat_large$hour), max(inat_large$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}


# Testing and correcting for diurnal bias ---------------------------------

#plot all occurrence data to see if there is a bias for diurnal observations
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/all_species_human_obs.png", width=20,height=10,units="cm",res=200)
ggplot(inat, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(inat$hour), max(inat$hour), by = 1),1))
dev.off()

#extract the number of occurrences per hour for all species
hour_bias <- inat %>% count(hour)
#extract the number of occurrences per hour from one species (Berardius arnuxii)
single_species_occurrence <- inat %>% filter(species == "Berardius arnuxii") %>% count(hour)
normalized_occurrence <- merge(hour_bias, single_species_occurrence, by = "hour")

#divide the single species data by the total species data (to account for more observations in daytime)
normalized_occurrence$normalized_count <- normalized_occurrence$n.y / normalized_occurrence$n.x

#plot the normalized data
ggplot(normalized_occurrence, aes(x = hour, y = normalized_count)) + geom_point() + geom_line()

#save out 
png("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/Berardius arnuxii_normalized.png", width=20,height=10,units="cm",res=200)
ggplot(normalized_occurrence, aes(x = hour, y = normalized_count)) + geom_point() + geom_line()
dev.off()

#repeat this for all data deficient species (there should be 11)
#create a dataframe with counts for all individual NA species, divided by total occurrence
hour_bias <- inat %>% count(hour)
#get total hour counts for each species
sps1 <- NA_species %>% filter(species == "Lagenorhynchus cruciger") %>% count(hour)
sps2 <- NA_species %>% filter(species == "Berardius arnuxii") %>% count(hour)
single_sps <- merge(sps1, sps2, by = "hour")

#combine with the total counts from all species
test <- merge(hour_bias, single_sps, by = "hour")
colnames(test) <- c("hour", "total occurrence", "Berardius occurrence")

hour_bias <- inat %>% count(hour)
colnames(hour_bias) <- c("hour", "total occurrence")
for(i in unique(NA_species$species)){
  sps <- NA_species %>% filter(species == "Mesoplodon traversii") %>% count(hour)
  hour_bias <- merge(hour_bias, sps, by = "hour")
}


#plotting total occurrences coloured by species
ggplot(inat, aes(x = hour, fill = species)) + geom_histogram()

#plotting species occurrences as a proportion of total occurrences
inat %>% group_by(hour, species) %>% summarise(total_occurrence = sum(hour, na.rm = TRUE)) %>% ggplot(aes(x = hour, y = total_occurrence, fill = species)) + ggstream::geom_stream(type = "proportion")

#same plot as before but filtered for species with no existing diel pattern data
NA_species %>% group_by(hour, species) %>% summarise(total_occurrence = sum(hour, na.rm = TRUE)) %>% ggplot(aes(x = hour, y = total_occurrence, fill = species)) + ggstream::geom_stream(type = "proportion")

#filtering for species that had limited data
low_data_species %>% group_by(hour, species) %>% summarise(total_occurrence = sum(hour, na.rm = TRUE)) %>% ggplot(aes(x = hour, y = total_occurrence, fill = species)) + ggstream::geom_stream(type = "proportion")

#filtering for species with under 50 occurrences
inat_small %>% group_by(hour, species) %>% summarise(total_occurrence = sum(hour, na.rm = TRUE)) %>% ggplot(aes(x = hour, y = total_occurrence, fill = species)) + ggstream::geom_stream(type = "proportion")

#filtering for species with 50-200 occurrences
inat_med %>% group_by(hour, species) %>% summarise(total_occurrence = sum(hour, na.rm = TRUE)) %>% ggplot(aes(x = hour, y = total_occurrence, fill = species)) + ggstream::geom_stream(type = "proportion")

#filtering for species with over 200 occurrences
inat_large %>% group_by(hour, species) %>% summarise(total_occurrence = sum(hour, na.rm = TRUE)) %>% ggplot(aes(x = hour, y = total_occurrence, fill = species)) + ggstream::geom_stream(type = "proportion")


#download NA_species as a csv to manually add information on if the observed specimen is dead or alive
#write.csv(NA_species, "C:/Users/ameli/Downloads/NA_species.csv")
url <- 'https://docs.google.com/spreadsheets/d/1kJypfUy5Nhih-_e8ujiPTT_iYaOX9j2JGRLrnc1IZ44/edit?usp=sharing'
NA_species <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#re-analyze the data excluding the dead specimens, this leaves 8 species with data
NA_species <- NA_species %>% filter(specimenStatus == "Alive")
ggplot(NA_species, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(NA_species$hour), max(NA_species$hour), by = 1),1)) + facet_wrap(~species)

#save out as individual plots
for(i in 1:8){
  png(paste0("C:/Users/ameli/Downloads/inaturalist_cetacean_activity/NA_species_not_norm_alive", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(NA_species, aes(x = hour)) + geom_histogram() + scale_x_continuous(breaks = round(seq(min(NA_species$hour), max(NA_species$hour), by = 1),1)) + facet_wrap_paginate(~species, scales = "free", ncol = 1, nrow = 1, page = i))
  dev.off()
}