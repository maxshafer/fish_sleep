#Packages
library(stringr)
library(here)
library(ggtree)
library(gsheet)
library(dplyr)
library(readxl)
library(tidyr)
install.packages("networkD3")
library(networkD3)

#how well does data from different sources (categories) agree with each other

#import fish database
fish_df <- read.csv(here("sleepy_fish_database_local.csv"))

#look for exact matches
df$match <- "Idk"
  
for(i in 1:length(fish_df$Species_name)){
  if(fish_df[i, "column_1"] == fish_df[i, "Column_2"]){
    fish_df[i, "match"] <- "Yes"
  } else if(fish_df[i, "column_1"] %in% c("crepuscular/diurnal", "crepuscular", "diurnal") & fish_df[i, "column_2"] %in% c("crepuscular/diurnal", "crepuscular", "diurnal")){
    fish_df[i, "match"] <- "Approximate"
  } else if(fish_df[i, "column_1"] %in% c("crepuscular/nocturnal", "crepuscular", "nocturnal") & fish_df[i, "column_2"] %in% c("crepuscular/nocturnal", "crepuscular", "nocturnal")){
    fish_df[i, "match"] <- "Approximate"
  } else {
    fish_df[i, "match"] <- "No"
  }
}

#allow for approximate matches
fish <- c("crepuscular/diurnal", "crepuscular", "diurnal")
a <- "crepuscular"
b <- "diurnal"

a %in% fish & b %in% fish

#How to visualize agreement across confidence levels

#load in sleepy fishies
url <- "https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit?gid=0#gid=0"
fish_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

View(fish_full)
#subset to only include species with multiple sources
fish_full <- fish_full %>% filter(Source.2 != "")

#if applicable, filter to confidence level 

#only keep species name and confidence data
fish_full <- fish_full[, c(1, 10:21)]

#for now also filter to the section I have done
fish_full <- fish_full[1182:1194,]

#rename columns 
colnames(fish_full) <- c("Species_name", "5", "5.1", "4", "4.1", "4.2", "3", "3.1", "3.2", "2", "2.1", "2.2", "2.3")



