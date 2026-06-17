##Packages we will use ---------------------------------------------------
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)
#open tree of life
library(rotl)
#adds timescale
#library(deeptime)
#colours
library(RColorBrewer)
#apply two separate colour palettes
library(ggnewscale)
#more colours
library(pals)
#useful
library(tidyr)
#package to work with shapefiles
#install.packages("sf")
library(sf)
#also helps with shapefiles
library(raster)
#colours
library(viridis)
#sankey diagram
# install.packages("ggsankey")
library(ggsankey)
library(ggpubr)
## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)
library(DescTools)
library(ggbeeswarm)
library(rstatix)

#extra packages
library(lubridate)
#install.packages("deeptime")
library(scales)
#install.packages("ggforce")
library(ggforce)
library(forcats)
#install.packages("ggdist")
library(ggdist)
library(knitr)
#install.packages("kableExtra")
library(kableExtra)
#install.packages("webshot")
library(webshot)
library(forcats)
library(ggpmisc)
library(stats)
library(phyloint)
library(patchwork)
library(ggridges)
library(rlang)
#pantheria database
#library(trait.data)
# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

# Section 1: Transform cetacean data into wide format -------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#take only the columns we're interested in
diel_full <- diel_full[, c(1, 3, 9:16)] 

diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")

#remove species without any information, drops from 98 species to 90
diel_full <- diel_full[!is.na(diel_full$Confidence),]

#convert strings to lowercase
diel_full$Conf_1 <- tolower(diel_full$Conf_1)
diel_full$Conf_2_daytime <- tolower(diel_full$Conf_2_daytime)
diel_full$Conf_2_nighttime <- tolower(diel_full$Conf_2_nighttime)
diel_full$Conf_3_PAM <- tolower(diel_full$Conf_3_PAM)
diel_full$Conf_3_Stomach_bycatch <- tolower(diel_full$Conf_3_Stomach_bycatch)
diel_full$Conf_4 <- tolower(diel_full$Conf_4)
diel_full$Conf_5 <- tolower(diel_full$Conf_5)

#separate all entries in a category of evidence (1-5) into separate columns
diel_full <- separate(data = diel_full, col = Conf_1, into = c("Conf1.1", "Conf1.2", "Conf1.3", "Conf1.4", "Conf1.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_daytime, into = c("Conf2.1", "Conf2.2", "Conf2.3", "Conf2.4", "Conf2.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_nighttime, into = c("Conf2N.1", "Conf2N.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_PAM, into = c("Conf3.1", "Conf3.2", "Conf3.3", "Conf3.4", "Conf3.5", "Conf3.6", "Conf3.7", "Conf3.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_Stomach_bycatch, into = c("Conf3ByS.1", "Conf3ByS.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_4, into = c("Conf4.1", "Conf4.2", "Conf4.3", "Conf4.4", "Conf4.5", "Conf4.6", "Conf4.7", "Conf4.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_5, into = c("Conf5.1", "Conf5.2", "Conf5.3", "Conf5.4"), sep = ",")

#replace strings 
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral-variable", "cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral-invariate", "cathemeral", x)}))

#may change to cathemeral/nocturnal/crepuscular in future, there are 5 weakly nocturnal crepuscular entries
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal/crepuscular", "nocturnal/crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal", "nocturnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-diurnal", "diurnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral-dvm", "cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("nocturnal-dvm", "nocturnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("-", "/", x)}))

#save out
write.csv(diel_full, file = here("cetacean_confidence_wide.csv"), row.names = FALSE)

# Section 2: Transform cetacean data into long format -------------
#read in wide format data
diel_full <- read.csv(here("cetacean_confidence_wide.csv"))
#convert into long format
diel_full_long <- diel_full %>% pivot_longer(cols = Conf1.1:Conf5.4, names_to = "column", values_to = "value")
#remove whitespace
diel_full_long$value <- str_trim(diel_full_long$value)

#remove unclear as an option since it gives no new information
#this removes 46 entries all from confidence level 1 and 2 from 25 species. Most of level 1 and 2 data isn't informative
diel_full_long[diel_full_long == "unclear"] <- NA

#remove rows with empty values
diel_full_long <- diel_full_long[diel_full_long$value != "",]
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

#check for any weird values
unique(diel_full_long$value)

#save out
write.csv(diel_full_long, file = here("cetacean_confidence_long_df.csv"), row.names = FALSE)

# Section 3: Objective method for making calls on activity patterns  ------------------------------------------------------------

#read in the confidence data in long format
test_diel_long <- read.csv(here("cetacean_confidence_long_df.csv"))
 
#take only the first part of the column name (ie conf1, conf2)
test_diel_long$column <- gsub("\\..*","",test_diel_long$column)
#we won't worry about different types of level 3 and level 2 evidence
test_diel_long$column <- str_replace(test_diel_long$column, pattern = "Conf3ByS", replacement = "Conf3")
test_diel_long$column <- str_replace(test_diel_long$column, pattern = "Conf2N", replacement = "Conf2")

#deal with x/cathemeral species first, makes no difference in final call whether you call them dinoc or cathemeral
unique(test_diel_long$value)
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "nocturnal/cathemeral", replacement = "cathemeral")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "diurnal/cathemeral", replacement = "cathemeral")

#deal with unclear values which are all in confidence level 2
#the 3 species it effects are not in the final mammal tree so we can keep as unclear and exclude later
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "unclear/diurnal", replacement = "unclear")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "unclear/cathemeral", replacement = "unclear")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "unclear/nocturnal", replacement = "unclear")

#separate out crepuscularity into its own column
#confidence level 2 data will be unclear/crepuscular since they don't show evidence of nocturnal or diurnal activity
test_diel_long <- separate(test_diel_long, col = value, into = c("new_diel", "crepuscular"), sep = "/")

#replace "unclear" with NA since it adds no new information
test_diel_long[test_diel_long == "unclear"] <- NA #84 species with some data

write.csv(test_diel_long, here("cetacean_confidence_long_final.csv"), row.names = FALSE)

#a custom function to take the mode where x is a list of numbers corresponding to the number of times each diel pattern appears for a given species
#if x = c(2,1,3) then for that species diurnal appeared twice, nocturnal appeared once and cathemeral appeared thrice
which.max.simple=function(x,na.rm=TRUE,tie_value="NA"){
  if(na.rm)
  {
    x=x[!is.na(x)] #removes NA values
  }
  if(length(x)==0) #if there is no activity pattern data return NA
  {
    return(NA)
  }
  if(length(x)==1) #if there is no activity pattern (will return an integer zero) data return NA
  { if(x == 0){
    return(NA)
  }
  }
  maxval=max(x) #takes the list of occurrences and picks whichever is the largest (ie cathemeral appears 3 times)
  if(is.na(maxval)) #if the highest value is NA then return NA
  {
    return(NA)
  }
  if(sum(x %in% maxval) > 1) #if the highest number appears twice or more (ie there's a tie)
  {
    # Ties exist, figure out what to do with them. Two options
    if(tie_value=="NA") #if there's a tie return NA
    {
      return(NA)
    }
    
    if(tie_value=="random") #if there's a tie, randomly chose one over the other
    {
      tie_postions=which(x==maxval)
      return(sample(tie_postions,size=1))
    }
    
    if(tie_value=="first") #if there's a tie, chose the first value
    {
      tie_postions=which(x==maxval)
      return(tie_postions[1])
    }
    
  } else
  {
    return(which.max(x)) #return the maximum value
  }
}

#the below function takes the all the entries for a given species in the specified confidence levels (ie for first pass its conf4)
#first checks is there is a mode for that confidence level or if its inconclusive (will return an NA)
#takes the diel patterns in each entry for that species and assigns it a number based by matching it to a position the list of unique values
#example 1 = diurnal, 2 = nocturnal, 3 = cathemeral, tabulate counts the number of times each of these numbers appears in the input
#example, diurnal appears twice, nocturnal appears once and cathemeral appears three times. So it returns cathemeral as the activity pattern

#example species
#x <- filter(test_diel_long, Species_name == "Pontoporia blainvillei")

# my method: take mode of mutliple level 4 sources if they exist, if not mode of Conf3 and Conf4
#if unclear add in Conf5 data, if unclear use single level 4 source (if it exists), if still unclear call cathemeral-variable
#if only level 1 and 2 data exist use these to make the call 
#can categorize each of these as A,B,C,D etc based on what evidence was used to make  the call

tabulateFuncCet <- function(x) {
  if(x %>% filter(column == "Conf4") %>% nrow() > 1 & !is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4")], unique(x$new_diel))), tie_value = "NA"))){
    activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4")], unique(x$new_diel))), tie_value = "NA")]
    activity_pattern <- paste(activity_pattern, "A")
  } else {
    if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA"))|(nrow(x[x$column %in% c("Conf3", "Conf4"),]) == 1)){
      if (nrow(x[x$column %in% c("Conf4"),]) < 1) {
        if (nrow(x[x$column %in% c("Conf3"),]) != 1){
          if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA"))) {
            if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3", "Conf1")], unique(x$new_diel))), tie_value = "NA")) | nrow(x[x$column %in% c("Conf3", "Conf4", "Conf5"),]) > 0){
              activity_pattern <- "cathemeral-variable J"
            } else{
              activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3", "Conf1")], unique(x$new_diel))), tie_value = "NA")]
              activity_pattern <- paste(activity_pattern, "I")
            }
          } else {
            activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA")]
            activity_pattern <- paste(activity_pattern, "H")
          }
        } else {
          activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf3")], unique(x$new_diel))), tie_value = "NA")]
          activity_pattern <- paste(activity_pattern, "G")
        }
      } else {
        activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4")], unique(x$new_diel))), tie_value = "NA")]
        activity_pattern <- paste(activity_pattern, "F")
      }
    } else {
      activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA")]
      activity_pattern <- paste(activity_pattern, "E")}
  }
  return(activity_pattern)
}

#run each species through this function, 82 species with activity pattern data (di, noc or cath)
activity_pattern_df <- test_diel_long[!is.na(test_diel_long$new_diel),] %>% group_by(Species_name) %>% do(tabulated_diel_pattern = tabulateFuncCet(.)) %>% unnest()

activity_pattern_df <- separate(activity_pattern_df, col = "tabulated_diel_pattern", into = c("tabulated_diel_pattern", "level"), sep = " ")
activity_pattern_df$level <- activity_pattern_df$level %>% replace_na("G")
table(activity_pattern_df$level)
table(activity_pattern_df$tabulated_diel_pattern)

unique(activity_pattern_df$tabulated_diel_pattern)

#replace cathemeral-variable with cathemeral since we aren't delineating between the two in the analysis
activity_pattern_df$tabulated_diel_pattern <- str_replace(activity_pattern_df$tabulated_diel_pattern, pattern = "cathemeral-variable", replacement = "cathemeral")

#function to determine if species is crepuscular
#for each species, how many sources say they are crepuscular? What is the confidence of these sources?
#percentage of each source
test_diel_long[is.na(test_diel_long)] <- "0" #replaces all the nas in crepuscular column with 0
test_diel_long$crepuscular <- str_replace(test_diel_long$crepuscular, pattern = "crepuscular", replacement = "1")
test_diel_long$crepuscular <- as.numeric(test_diel_long$crepuscular) #mark all crepuscular species with a value of 1
test_diel_long$total <- 1 #used to calculate the percentage of crepuscular sources out of the total sources

#use to test function below
x <- filter(test_diel_long, Species_name == "Monodon monoceros")

#if category 3 greater than 50%, category 4 higher than 0% and if category 2 higher than 75% then call this source crepuscular
#if majority of evidence categories call it crepuscular then designate it crepuscular
#category 5 has no crepuscular calls and category 1 has only one, so don't use these to make calls

tabulateCrep = function(x){
  df <- aggregate(x$crepuscular, by = list(Category = x$column), FUN = sum)
  #confidence level 1 and 5 don't have any information on crepuscularity so remove them, they will always be zero
  df <- df[df$Category %in% c("Conf2", "Conf3", "Conf4"), ]
  df2 <- aggregate(x$total, by = list(Category = x$column), FUN = sum)
  
  #if species only has level 1 and/or 5 data then total evidence will always be FALSE (we can't determine crepuscularity)
  if(all(!df2$Category %in% c("Conf2", "Conf3", "Conf4")) == TRUE){
    total_evidence <- "no sources"
    return(total_evidence)
  }
  
  df2 <- df2[df2$Category %in% c("Conf2", "Conf3", "Conf4"), ]
  df2 <- merge(df, df2, by = "Category")
  df2$percentage <- round((df2$x.x/df2$x.y) * 100, digits = 1)
  #define what percentage of each category we want to call a species crepuscular,
  df3 <- data.frame(Category = c("Conf2", "Conf3", "Conf4"), cutoff = c(50, 50, 20))
  df2 <- merge(df2, df3, by = "Category")  #species can have Conf2, Conf3 and/or Conf4 evidence, when merged any missing categories will be dropped
  df2$crep_evidence <- df2$percentage >= df2$cutoff
  df2[df2$crep_evidence == TRUE, "crep_evidence"] <- "crepuscular"
  df2[df2$crep_evidence == FALSE, "crep_evidence"] <- "non"
  unique_percents <- unique(df2$crep_evidence)
  #check if there's a tie, since we only use level 2,3 and 4, a tie will only occur when using two evidence levels 
  if(nrow(df2) == 2){
    if("crepuscular" %in% df2$crep_evidence & "non" %in% df2$crep_evidence){
      #to break ties could look at which has more sources (to avoid making calls off one source)
      #total_evidence <- df2[which.max(df2$x.y), "crep_evidence"]
      #instead could also break ties by looking at the higher confidence level (this will always be the second row but use which to be safe)
      total_evidence <- df2[which(df2$Category == max(df2$Category)), "crep_evidence"]
    } else{
      total_evidence <- unique_percents[which.max(tabulate(match(df2$crep_evidence, unique_percents)))]
    }
  } else {total_evidence <- unique_percents[which.max(tabulate(match(df2$crep_evidence, unique_percents)))]}
  #additional screen: if there is level 4 evidence then evaluate to true 
  if(nrow(filter(df2, Category == "Conf4"))== 1){
    if(df2[df2$Category == "Conf4", "crep_evidence"] == "crepuscular"){
      total_evidence <- "crepuscular"
    }
  } 
  return(total_evidence)
}

crep_df <- test_diel_long %>% group_by(Species_name) %>% do(tabulated_crep = tabulateCrep(.)) %>% unnest()

#values for the crepuscularity pipeline flowchart
test <- test_diel_long %>% filter(column %in% c("Conf2", "Conf3", "Conf4")) %>% group_by(Species_name, column) %>% 
  summarize(sum_crep = sum(crepuscular), sum_total = sum(total))  %>% mutate(percent_crep = (sum_crep/sum_total)*100) %>% 
  pivot_wider(id_cols = !c(sum_total, sum_crep), names_from = "column", values_from = percent_crep)

test <- test %>% mutate(evidence2 = as.numeric(Conf2 >50), evidence3 = as.numeric(Conf3 >50), evidence4 = as.numeric(Conf4 >20)) %>%
  mutate(total_evidence = evidence2 + evidence3 + evidence4)

table(test$evidence4)


final_df <- merge(crep_df, activity_pattern_df, by = "Species_name")
final_df$tabulated_diel <- final_df$tabulated_diel_pattern
for(i in 1:nrow(final_df)){
  if(final_df[i, "tabulated_crep"] == "crepuscular"){
    final_df[i, "tabulated_diel"] <- paste(test[i, "tabulated_diel_pattern"], "crepuscular", sep = "/")
  }
}

current_dataset  <- final_df[, c("Species_name", "tabulated_diel")]
previous_dataset <- read.csv(here("cetacean_tabulated_full.csv"))

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
previous_dataset <- previous_dataset[previous_dataset$tips %in% mam.tree$tip.label, ]
table(previous_dataset$tabulated_diel)
current_dataset$tips <- str_replace(current_dataset$Species_name, pattern = " ", replacement = "_")
current_dataset <- current_dataset[current_dataset$tips %in% mam.tree$tip.label, ]
table(current_dataset$tabulated_diel)

#check that nothing about the data has changed since running it last 
all(previous_dataset == current_dataset)
if(all(previous_dataset == current_dataset) == FALSE) stop("Dataset is not the same!")

#add a column for tips, formatted as the species names appear in the phylogenetic tree
final_df$tips <- str_replace(final_df$Species_name, pattern = " ", replacement = "_")

#save out the new tabulated activity pattern dataframe
write.csv(final_df[, c("Species_name", "tabulated_diel", "tips")], here("cetacean_tabulated_full.csv"), row.names = FALSE)

# Section 3.5 Save out cetacean data frame with additional details -------
#load in the dataframe with the tabulated activity patterns (objective calls based on source concordance)
cetaceans_tabulated_full <- read.csv(here("cetacean_tabulated_full.csv")) #84 species with data 82??

#load in full primary source dataframe, 98 species and subspecies
url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit?usp=sharing'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#add in the tabulated diel patterns
cetaceans_full <- cetaceans_full %>% select(-c(Diel_Pattern_1, Diel_Pattern_2, Diel_Pattern_3))
cetaceans_full <- merge(cetaceans_full, cetaceans_tabulated_full, by = "Species_name", all.x = TRUE)
#save out full version with sources
write.csv(cetaceans_full, here("cetaceans_full_with_sources.csv"), row.names = FALSE)

#remove unnecessary columns
cetaceans_full <- cetaceans_full[,c("Species_name", "Confidence", "Parvorder", "Family", "tabulated_diel", "tips")]
cetaceans_full$tips <- str_replace(cetaceans_full$Species_name, pattern = " ", replacement = "_")
colnames(cetaceans_full) <- c("Species_name", "Confidence", "Parvorder", "Family", "Diel_Pattern", "tips")

#add suborder taxonomic info for future reference
cetaceans_full$Suborder <- "Whippomorpha"
cetaceans_full$Order <- "Artiodactyla"

#create the three databases we will use 
#Diel_Pattern includes all 6 possible trait states: di, di/crep, noc, noc/crep, cath, cath/crep
#Max_crep will include 4 trait states and maximize crepuscularity: di, noc, cath and crep (di/crep, noc/crep, cath/crep)
cetaceans_full$max_crep <- cetaceans_full$Diel_Pattern
cetaceans_full$max_crep <- str_replace(cetaceans_full$max_crep, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
cetaceans_full$max_crep <- str_replace(cetaceans_full$max_crep, pattern = "diurnal/crepuscular", replacement = "crepuscular")
cetaceans_full$max_crep <- str_replace(cetaceans_full$max_crep, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
#Max_dinoc will include 4 trait states and maximize di and noc, di (including di/crep), noc (including noc/crep), cath, crep (cath/crep)
cetaceans_full$max_dinoc <- cetaceans_full$Diel_Pattern
cetaceans_full$max_dinoc <- str_replace(cetaceans_full$max_dinoc, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
cetaceans_full$max_dinoc <- str_replace(cetaceans_full$max_dinoc, pattern = "diurnal/crepuscular", replacement = "diurnal")
cetaceans_full$max_dinoc <- str_replace(cetaceans_full$max_dinoc, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#create a column with the max confidence level for that species (out of the confidence level for all sources)
#the confidence values are characters so convert to numerics and then take the maximum value
cetaceans_full$Confidence[is.na(cetaceans_full$Confidence)] <- 0
cetaceans_full$Confidence <- gsub(",", "\\1 ", cetaceans_full$Confidence)
cetaceans_full$Confidence <- strsplit(cetaceans_full$Confidence, " ")
cetaceans_full$Confidence <- lapply(cetaceans_full$Confidence, max)
cetaceans_full$Confidence <- unlist(cetaceans_full$Confidence)
cetaceans_full[cetaceans_full == 0.5] <- 1 #these are all the inaturalist observations

cetaceans_full <- cetaceans_full %>% select("Species_name", "Order", "Suborder", "Parvorder", "Family", "Diel_Pattern", "max_crep", "max_dinoc", "Confidence", "tips")

#save out a local copy in case google goes bankrupt
write.csv(cetaceans_full, file = here("cetaceans_full.csv"), row.names = FALSE)

# Section 4 Concordance table within confidence levels -------------------------------------------
diel_full_long <- read.csv(here("cetacean_confidence_long_df.csv"))
#read in the tabulated activity patterns
diel_full <- read.csv(here("cetaceans_full.csv"))
diel_full <- merge(diel_full[, c("Species_name", "Parvorder", "Diel_Pattern", "max_crep")], diel_full_long[c("Species_name", "column", "value")])

diel_full$column <- substr(diel_full$column, 1,5)

diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/crepuscular", "crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/diurnal", "diurnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/nocturnal", "nocturnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/cathemeral", "cathemeral", x)}))

diel_full <- data.frame(lapply(diel_full, function(x) {gsub("diurnal/cathemeral", "cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("nocturnal/cathemeral", "cathemeral", x)}))

#for max crep dataset
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral/crepuscular", "crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("diurnal/crepuscular", "crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("nocturnal/crepuscular", "crepuscular", x)}))

#filter
mulitple_sources <- diel_full %>% count(Species_name) %>% filter(n>1)
diel_full_filtered <- diel_full[diel_full$Species_name %in% mulitple_sources$Species_name,]
#this removes 12 species without an informative second source
#"Berardius arnuxii", "Caperea marginata", "Inia boliviensis","Inia humboldtiana", "Lagenorhynchus albirostris","Lissodelphis borealis", "Lissodelphis peronii","Mesoplodon hotaula"        
#"Mesoplodon mirus", "Orcaella heinsohni","Phocoena sinus", "Sousa teuszii"   

#check to see if there is a difference in accuracy in mysticetes vs odontocetes
#diel_full_filtered <- filter(diel_full_filtered, Parvorder == "Mysticeti")
#diel_full_filtered <- filter(diel_full_filtered, Parvorder == "Odontoceti")

concordance <- as.data.frame(table(diel_full_filtered$max_crep, diel_full_filtered$value))
colnames(concordance) <- c("actual", "predicted", "freq")
totals_df <- aggregate(concordance$freq, by=list(Category=concordance$actual), FUN=sum)
colnames(totals_df) <- c("actual", "total")
concordance <- merge(concordance, totals_df, by = "actual")
concordance$percent <- round(concordance$freq / concordance$total * 100, 1)

confusion_plot_cet <-
  ggplot(concordance, aes(actual, predicted, fill = percent)) + geom_tile() + geom_text(aes(label = paste0(percent, "%")), size = 3) +
  scale_fill_gradient(low = "#F5FBFF", high = "#0070D1") + 
  labs(x = "Actual (final activity pattern)", y = "Predicted (activity pattern of individual source)") + 
  theme_void() +
  scale_x_discrete(labels = c("Cathemeral", "Crepuscular", "Diurnal", "Nocturnal")) +
  scale_y_discrete(labels = c("Cathemeral", "Crepuscular", "Diurnal", "Nocturnal")) +
  theme(legend.position = "none", axis.text = element_text(size = 9), axis.title = element_text(size = 11))

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/cetacean_all_conf_levels_confusion_matrix.pdf", width = 4, height = 2)
confusion_plot_cet 
dev.off()

# Section 4.5: Concordance between confidence levels ----------------------
diel_full_long <- read.csv(here("cetacean_confidence_long_df.csv"))

#remove unclear as an option, since we aren't interested in if unclear matches with unclear
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/nocturnal", replacement = "nocturnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/diurnal", replacement = "diurnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/crepuscular", replacement = "crepuscular")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/cathemeral", replacement = "crepuscular")

diel_full_long[diel_full_long == "unclear"] <- NA
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

#check which diel patterns we are comparing
unique(diel_full_long$value)

#check to see if there is a difference in concordance for mysticeti vs odontoceti
#diel_full_long <- filter(diel_full_long, Parvorder == "Odontoceti")
#diel_full_long <- filter(diel_full_long, Parvorder == "Mysticeti")

#get a list of all the species with more than one source (should be most of them)
species_list <- table(diel_full_long$Species_name)
#should be 72 species with all cetaceans w multiple sources, 67 with only cetaceans in tree with mulitple sources
species_list <- names(species_list[species_list > 1])

#function Max wrote for comparing entries
compTwo <- function(comp1 = "comp1", comp2 = "comp2") {
  
  if(any(is.na(c(comp1, comp2)))) {
    return(NA)
  } else {
    #splits any entries with a backslash into two components (ie nocturnal/crepuscular into nocturnal and crepuscular)
    comp1 <- str_split(comp1, "/")[[1]]
    comp2 <- str_split(comp2, "/")[[1]]
    #then compares if any of the components match
    if(any(comp1 %in% comp2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
}

#apply this function across all species with multiple entries
output <- lapply(species_list, function(species) {
  
  #filter for one species at a time
  df <- diel_full_long[diel_full_long$Species_name == species,]
  #rename the column names to be unique for every entry for this species (ie for multiple column 2 entries column 2.1, 2.2 etc)
  df$column <- make.unique(df$column)
  
  #converts the dataframe so it compares every entry with each other (ie for A,B,C A-A, A-B, A-C, B-A, B-B, B-C, etc)
  df_lists_comb <- expand(df, nesting(var = column, vector = value), nesting(var2 = column, vector2 = value), .name_repair = "universal")
  
  df_lists_comb <- df_lists_comb %>% filter(var != var2) %>% arrange(var, var2) %>% mutate(vars = paste0(var, ".", var2)) %>% select(contains("var"), everything())
  
  #evaluates the activity patterns for each of these sources and returns if they agree or not (TRUE or FALSE)
  comparisons <- df_lists_comb %>% group_by(vars) %>% mutate(comp = compTwo(comp1 = vector, comp2 = vector2))
  #manipulate the strings for both variable names to revert them back to the original name (back to column 2 from col 2.1)
  comparisons$var <- str_sub(comparisons$var, start = 1, end = 5)
  comparisons$var2 <- str_sub(comparisons$var2, start = 1, end = 5)
  
  #create a column returning the comparison being made (ie col2-col2, col1-col2, etc)
  comparisons$var_final <- paste(comparisons$var, comparisons$var2, sep = "_")
  
  #return just the comparison result column (TRUE or FALSE match) and the comparison being made (ie col1 vs col1)
  return(comparisons[,c("comp","var_final")])
})

#combine this list of results 
output <- Reduce(rbind, output)

table <- table(output$var_final)
# prop.table(table, margin = 1)
table2 <- as.data.frame(prop.table(table(output$var_final, output$comp), margin = 1))
table2$Comp1 <- sapply(str_split(table2$Var1, "_"), `[`, 1)
table2$Comp2 <- sapply(str_split(table2$Var1, "_"), `[`, 2)
table2 <- table2[table2$Var2 == TRUE,]
table2$count <- table

#want to make a plot that has both the frequency and the counts
table2$freq_count <- paste0((round(table2$Freq, 2) * 100), "%", "\n", "(n=", table2$count, ")")
plot_countfreq_cet <- table2[c(1:5, 7:10, 13:15, 19:20, 25), ] %>% 
  ggplot(., aes(x = Comp1, y = Comp2, fill = Freq, label = freq_count)) +
  geom_tile() + geom_text(size = 3) + scale_fill_viridis(begin = 0.2, end = 1, limits = c(1,0)) + 
  theme_minimal() + ylab("Primary source category") + xlab("Secondary source category") +
  scale_x_discrete(labels = c("A", "B", "C", "D", "E")) +
  scale_y_discrete(labels = c("A", "B", "C", "D", "E")) +
  theme(legend.position = "none")

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/cetacaean_btw_source_concordance.pdf", width = 7, height = 7, bg = "transparent")
#pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/cetacaean_btw_source_concordance_odontoceti.pdf", width = 7, height = 7, bg = "transparent")
#pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/cetacaean_btw_source_concordance_mysticeti.pdf", width = 7, height = 7, bg = "transparent")
plot_countfreq_cet
dev.off()


# Section 5: Concordance sankey -------------------------------------------

#create dataframe of the number of species that had activity patterns determined at each step
df <- data.frame(
  step_6 = c(rep("A. Multiple category D \n sources in concordance?",82)),
  step_5 = c(rep("B. Return category D \n (n = 26)", 26), rep("C. Category D + C \n in concordance?", 56)),
  step_4 = c(rep(NA, 26), rep("D. Return category D + C \n (n = 25)", 25), rep("E. Single category D source?", 31)),
  step_3 = c(rep(NA, 51), rep("F. Return single cateory D  \n (n = 5)", 5), rep("G. Single category C source?", 26)),
  step_2 = c(rep(NA, 56), rep("H. Return single category C \n (n = 13)", 13), rep("I. Category E + D + C \n sources in concordance?", 13)),
  step_1 = c(rep(NA, 69), rep("J. Return category E + D \n + C (n = 1)", 1), rep("K. Multiple category C, D, or E \n sources in concordance?", 12)),
  step_0 = c(rep(NA, 70), rep("L. Return category A \n (n = 7)", 7), rep("M. Else return \n cathemeral (n = 5)", 5))
)

#convert to long format for geomsankey
df <- df %>% make_long(step_0, step_1, step_2, step_3, step_4, step_5, step_6)
df <- df[!is.na(df$node), ]

blues <- c("#010661", "#070E8A","#070E8A", "#0044A3","#0044A3", "#0070D1","#0070D1","#2E9DFF","#2E9DFF","#8AC8FF","#8AC8FF","#B8DEFF","#B8DEFF")

sankey_cet <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = substr(node, 4, 300))) +
  geom_sankey(flow.alpha= 0.5, node.color = 0.5) + geom_sankey_label(size = 3, color = 1, fill = "white")  + 
  theme_sankey(base_size = 11) + scale_fill_manual(values = blues) +
  theme(legend.position = "none", axis.text.x = element_blank(), panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent', colour = NA)) + labs(x = NULL) 

sankey_cet

#save out to figure folder
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "cetacean_flowchart.pdf"), height = 7, width = 14)
sankey_cet
dev.off()
# Section 6: Crepuscularity sankey -------------------------------------------

#create dataframe of the number of species that had activity patterns determined at each step
df <- data.frame(
  step_4 = c(rep("A. Total species (n = 84) ",84)),
  step_3 = c(rep("B. No sources, non-crepuscular (n = 76)", 8), rep("C. Level B, C, D sources", 76)),
  step_2 = c(rep(NA, 8), rep("D. Category B sources majority (n = 70)", 32), rep("E. Category C sources majority (n = 33)", 22), rep("F. Category D sources (n = 39)", 22)),
  step_1 = c(rep(NA, 8), rep("G. No evidence (n = 23)", 20), rep("H. Crepuscular evidence (n = 10)", 12), rep("I. No evidence (n = 67)", 10), rep("J. Crepuscular evidence (n = 13)", 12), rep("K. No evidence (n = 29)", 10), rep("L. Crepuscular evidence (n = 10)", 12))
  #step_0 = c(rep(NA, 8), rep("M. Non-crepuscular (n = 53)", 53), rep("N. Crepuscular (n = 23)", 23))
  )

#convert to long format for geomsankey
df <- df %>% make_long(step_1, step_2, step_3, step_4)
df <- df[!is.na(df$node), ]

blues <- c("#010661", "#070E8A","#070E8A", "#0044A3","#0044A3", "#0070D1","#0070D1","#2E9DFF","#2E9DFF","#8AC8FF","#8AC8FF","#B8DEFF","#B8DEFF")

sankey_crep_cet <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = substr(node, 4, 300))) +
  geom_sankey(flow.alpha= 0.5, node.color = 0.5) + geom_sankey_label(size = 3, color = 1, fill = "white")  + 
  theme_sankey(base_size = 11) +  scale_fill_manual(values = blues) +
  theme(legend.position = "none", axis.text.x = element_blank(), panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent', colour = NA)) + labs(x = NULL) 

sankey_crep_cet + coord_flip()

#save out to figure folder
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "cetacean_flowchart.pdf"), height = 7, width = 14)
sankey_cet
dev.off()


