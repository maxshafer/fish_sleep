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
library(deeptime)
#colours
library(RColorBrewer)
#apply two separate colour palettes
library(ggnewscale)
#more colours
library(pals)
#useful
library(tidyr)
#colours
library(viridis)
#sankey diagram
library(ggsankey)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")



# Section 1: Transform cetacean data into wide format -------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#take only the columns we're interested in
diel_full <- diel_full[, c(1, 3:6, 9:16)]

diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")
diel_full$Diel_Pattern <- tolower(diel_full$Diel_Pattern)
diel_full$New_Pattern <- tolower(diel_full$New_Pattern)

#remove species without any information, drops from 98 species to 83
diel_full <- diel_full[!is.na(diel_full$Diel_Pattern),]

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
#double check this primacy source data I analzed
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal-tbd", "weakly-nocturnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral-tbd", "cathemeral", x)}))
#may change to cathemeral/nocturnal/crepuscular in future, there are 5 weakly nocturnal crepuscular entries
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal/crepuscular", "nocturnal/crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal", "nocturnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-diurnal", "diurnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral-dvm", "cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("nocturnal-dvm", "nocturnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("tbd", "unclear", x)}))
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
#we won't worry about different types of level 3 and level 2 evidence for now
test_diel_long$column <- str_replace(test_diel_long$column, pattern = "Conf3ByS", replacement = "Conf3")
test_diel_long$column <- str_replace(test_diel_long$column, pattern = "Conf2N", replacement = "Conf2")

#deal with x/cathemeral species first, makes no difference in final call whether you call them dinoc or cathemeral
unique(test_diel_long$value)
# test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "nocturnal/cathemeral", replacement = "nocturnal")
# test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "diurnal/cathemeral", replacement = "diurnal")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "nocturnal/cathemeral", replacement = "cathemeral")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "diurnal/cathemeral", replacement = "cathemeral")

#deal with unclear values which are all in confidence level 2
#tested with including these as unclear or not and is does change the call on any activity pattern, it does not
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "unclear/diurnal", replacement = "diurnal")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "unclear/cathemeral", replacement = "cathemeral")
test_diel_long$value <- str_replace_all(test_diel_long$value, pattern = "unclear/nocturnal", replacement = "nocturnal")

#separate out crepuscularity into its own column
#confidence level 2 data will be unclear/crepuscular since they don't show evidence of nocturnal or diurnal activity
test_diel_long <- separate(test_diel_long, col = value, into = c("new_diel", "crepuscular"), sep = "/")
#replace "unclear" with NA since it adds no new information
test_diel_long[test_diel_long == "unclear"] <- NA

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
x <- filter(test_diel_long, Species_name == "Phocoena spinipinnis")

# my method: take mode of mutliple level 4 sources if they exist, if not mode of Conf3 and Conf4
#if unclear add in Conf5 data, if unclear use single level 4 source (if it exists), if still unclear call cathemeral-variable
#if only level 1 and 2 data exist use these to make the call 
#can categorize each of these as A,B,C,D etc based on what evidence was used to make  the call

tabulateFunc2 <- function(x) {
  #first check if there is multiple level 4 confidence entries
  if(x %>% filter(column == "Conf4") %>% nrow() > 1){
    if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4")], unique(x$new_diel))), tie_value = "NA"))) {
      if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA"))) {
        if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA"))) {
          activity_pattern <- "cathemeral-variable" #result if there is a tie between all confidence levels (case D)
        } else {
          activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA")]
          #result if level 4 and level 4+3 are inconclusive (case C)
        }
      } else {
        activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA")]
        #result if level 4 alone is inconclusive (case B)
      }
    } else {
      activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4")], unique(x$new_diel))), tie_value = "NA")]
      #result if there is multiple level 4 evidence and it is not inconclusive (case A)
    }
  } else {
    if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA"))) {
      if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA"))) {
        #if there is a tie use the single level 4 entry to make the call 
        activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4")], unique(x$new_diel))), tie_value = "NA")]
        #if there is no level 4 source it will return an NA and move on to here
        if(NA %in% activity_pattern){
          #only use level 1 and 2 to make the call if that's the only confidence levels (and check if there is a tie, if its tied to go else, call it cathemeral variable)
          if(nrow(x[x$column %in% c("Conf1", "Conf2"),]) == nrow(x) & !is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf2", "Conf1")], unique(x$new_diel))), tie_value = "NA"))){
            activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf2", "Conf1")], unique(x$new_diel))), tie_value = "NA")]
          } else {activity_pattern <- "cathemeral-variable"} #this is the case where there are multiple level 1,2,3 or 5 sources that are tied, so we can call it cathemeral-variable (case I)
          }
        } else {
          activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf5", "Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA")]
          #result if level 4 and level 4+3 are inconclusive (case F)
          }
      } else {
        activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("Conf4", "Conf3")], unique(x$new_diel))), tie_value = "NA")]
        #result if level 4 alone is inconclusive (case E)
        } 
  } 
  
  return(activity_pattern)  
}

#run each species through this function, x species with activity pattern data (di, noc or cath)
activity_pattern_df <- test_diel_long[!is.na(test_diel_long$new_diel),] %>% group_by(Species_name) %>% do(tabulated_diel_pattern = tabulateFunc2(.)) %>% unnest()

unique(activity_pattern_df$tabulated_diel_pattern)

#replace cathemeral-variable with cathemeral since we aren't delineating between the two in the analysis
activity_pattern_df$tabulated_diel_pattern <- str_replace(activity_pattern_df$tabulated_diel_pattern, pattern = "cathemeral-variable", replacement = "cathemeral")

#function to determine if species is crepuscular
#for each species, how many sources say they are crepuscular? What is the confidence of these sources?
#percentage of each source
test_diel_long[is.na(test_diel_long)] <- "0" #replaces all the nas in crepuscular column with 0
test_diel_long$crepuscular <- str_replace(test_diel_long$crepuscular, pattern = "crepuscular", replacement = "1")
test_diel_long$crepuscular <- as.numeric(test_diel_long$crepuscular) #mark all crepuscular species with a value of 1
test_diel_long$total <- 1 #used to calculate the percentage of crepuscular species out of the total species

#use to test function below
#x <- filter(test_diel_long, Species_name == "Eubalaena glacialis")

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
    total_evidence <- FALSE
    return(total_evidence)
  }
  
  df2 <- df2[df2$Category %in% c("Conf2", "Conf3", "Conf4"), ]
  df2 <- merge(df, df2, by = "Category")
  df2$percentage <- round((df2$x.x/df2$x.y) * 100, digits = 1)
  #define what percentage of each category we want to call a species crepuscular,
  df3 <- data.frame(Category = c("Conf2", "Conf3", "Conf4"), cutoff = c(50, 50, 20))
  df2 <- merge(df2, df3, by = "Category")  #species can have Conf2, Conf3 and/or Conf4 evidence, when merged any missing categories will be dropped
  df2$crep_evidence <- df2$percentage >= df2$cutoff
  #find if there are more instances of crepuscular evidence or not. This doesn't work I think TRUE always overrides
  #so every species with at least one level of crepuscular evidence will say TRUE
  unique_percents <- unique(df2$crep_evidence)
  #check if there's a tie, since we only use level 2,3 and 4, a tie will only occur when using two evidence levels 
  if(nrow(df2) == 2){
    if(TRUE & FALSE %in% df2$crep_evidence){
      #to break ties could look at which has more sources (to avoid making calls off one source)
      #total_evidence <- df2[which.max(df2$x.y), "crep_evidence"]
      #instead could also break ties by looking at the higher confidence level (this will always be the second row but use which to be safe)
      total_evidence <- df2[which(df2$Category == max(df2$Category)), "crep_evidence"]
    } else{
      total_evidence <- unique_percents[which.max(tabulate(match(df2$crep_evidence, unique_percents)))]
    }
  } else {total_evidence <- unique_percents[which.max(tabulate(match(df2$crep_evidence, unique_percents)))]}
  #additional screen: if there is level 4 evidence then evaluate to true 
  if(nrow(filter(df2, Category == "Conf4"))==1){
    total_evidence <- any(df2[which(df2$Category == "Conf4"), "crep_evidence"], total_evidence)
  } else {return(total_evidence)
  }
}

#can't determine crepuscular activity from level or level 5 data. So manually fix any species with only level 1 and/or 5 data

crep_df <- test_diel_long %>% group_by(Species_name) %>% do(tabulated_crep = tabulateCrep(.)) %>% unnest()

test <- merge(crep_df, activity_pattern_df, by = "Species_name")
test$tabulated_diel <- test$tabulated_diel_pattern
for(i in 1:nrow(test)){
  if(test[i, "tabulated_crep"] == TRUE){
    test[i, "tabulated_diel"] <- paste(test[i, "tabulated_diel_pattern"], "crepuscular", sep = "/")
  }
}

unique(test$tabulated_diel)

test <- test[, c("Species_name", "tabulated_diel")]
#save out the new tabulated activity pattern dataframe
write.csv(test, here("cetacean_tabulated_full.csv"), row.names = FALSE)

diel_full_new <- merge(diel_full, test, by = "Species_name")
#diel pattern 2 is all diel categories, diel pattern is maxcrep, new pattern is manually determined diel patterns, tabulated pattern is diel patterns determined from above function
diel_full_new <- diel_full_new[, c("Species_name", "Parvorder", "Diel_Pattern_2", "Diel_Pattern","New_Pattern", "Confidence", "tips", "tabulated_diel")]

# Section 3.5: Plot new vs old diel pattern comparison ------------------------

trait.data <- diel_full_new[diel_full_new$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue", "grey","pink")
custom.colours <- c("#dd8ae7","#FC8D62", "yellow", "#66C2A5", "green")
custom.colours.2 <- c("#dd8ae7", "peachpuff2" ,"#FC8D62", "yellow", "#66C2A5","green","black","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "tabulated_diel")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = tabulated_diel), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Tabulated activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

#with max crep designations
trait.data$tabulated_diel <- str_replace(trait.data$tabulated_diel, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$tabulated_diel <- str_replace(trait.data$tabulated_diel, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$tabulated_diel <- str_replace(trait.data$tabulated_diel, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
custom.colours.2 <- c("#dd8ae7", "peachpuff2", "#FC8D62","#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "tabulated_diel")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = tabulated_diel), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Tabulated activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3)
diel.plot


# Section 3.5 Save out cetacean data frame with additional details -------
#load in the dataframe with the tabulated activity patterns (objective calls based on source concordance)
cetaceans_tabulated_full <- read.csv(here("cetacean_tabulated_full.csv")) #83 species with data

#load in full primary source dataframe, 98 species and subspecies
url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit?usp=sharing'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#add in the tabulated diel patterns
cetaceans_full <- merge(cetaceans_full, cetaceans_tabulated_full, by = "Species_name", all.x = TRUE)
#save out full version with sources
write.csv(cetaceans_full, here("cetaceans_full_with_sources.csv"))

#remove unnecessary columns
cetaceans_full <- cetaceans_full[c("Species_name", "Confidence", "Parvorder", "tabulated_diel")]
#add a column for tips, formatted as the species names appear in the phylogenetic tree
cetaceans_full$tips <- str_replace(cetaceans_full$Species_name, pattern = " ", replacement = "_")
colnames(cetaceans_full) <- c("Species_name", "Confidence", "Parvorder", "Diel_Pattern", "tips")

#add suborder taxonomic info for future reference
cetaceans_full$Suborder <- "whippomorpha"

#rename the row names to be the tip names so it's easier to subset by the tree tip labels later
row.names(cetaceans_full) <- cetaceans_full$tips

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
cetaceans_full[cetaceans_full == 0.5] <- 0 #these are all the inaturalist observations

#save out a local copy in case google goes bankrupt
write.csv(cetaceans_full, file = here("cetaceans_full.csv"), row.names = FALSE)

#save out a version with hippos
whippomorpha <- read.csv(here("cetaceans_full.csv"))
whippomorpha <- rbind(whippomorpha, c("Hexaprotodon liberiensis", 3, "Hippopotamidae", "nocturnal/crepuscular", "Hexaprotodon_liberiensis", "whippomorpha", "crepusuclar", "nocturnal"))
whippomorpha <- rbind(whippomorpha, c("Hippopotamus amphibius", 4, "Hippopotamidae", "nocturnal/crepuscular", "Hippopotamus_amphibius", "whippomorpha", "crepusuclar", "nocturnal"))
rownames(whippomorpha) <- whippomorpha$tips
write.csv(whippomorpha, file = here("whippomorpha.csv"), row.names = FALSE)

# Section 4 Concordance table within confidence levels -------------------------------------------
diel_full_long <- read.csv(here("cetacean_confidence_long_df.csv"))
#this is cheating but read in the tabulated activity patterns
diel_full <- read.csv(here("cetaceans_full.csv"))
diel_full <- merge(diel_full[, c("Species_name", "Diel_Pattern", "max_crep", "max_dinoc")], diel_full_long[c("Species_name", "column", "value")])

diel_full$column <- substr(diel_full$column, 1,5)

diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/crepuscular", "crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/diurnal", "diurnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/nocturnal", "nocturnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/cathemeral", "cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("diurnal-low-confidence", "diurnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("nocturnal-low-confidence", "nocturnal", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("cathemeral-low-confidence", "cathemeral", x)}))

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
#"Berardius arnuxii"          "Caperea marginata"          "Inia boliviensis"           "Inia humboldtiana"         
#"Lagenorhynchus albirostris" "Lissodelphis borealis"      "Lissodelphis peronii"       "Mesoplodon hotaula"        
#"Mesoplodon mirus"           "Orcaella heinsohni"         "Phocoena sinus"             "Sousa teuszii"   
#of these only I humboldtiana is not the final tree

#to include only species in the final tree, is this necessary? If so I should do it consistently 
diel_full_filtered$tips <- str_replace(diel_full_filtered$Species_name, pattern = " ", replacement = "_")
diel_full_filtered <- diel_full_filtered[diel_full_filtered$tips %in% mam.tree$tip.label,]

concordance <- as.data.frame(table(diel_full_filtered$max_crep, diel_full_filtered$value))
colnames(concordance) <- c("actual", "predicted", "freq")
totals_df <- aggregate(concordance$freq, by=list(Category=concordance$actual), FUN=sum)
colnames(totals_df) <- c("actual", "total")
concordance <- merge(concordance, totals_df, by = "actual")
concordance$percent <- round(concordance$freq / concordance$total * 100, 1)

ggplot(concordance, aes(actual, predicted, fill = percent)) + geom_tile() + geom_text(aes(label = percent)) +
  scale_fill_gradient(low = "white", high = "dodgerblue") + labs(x = "Actual", y = "Predicted") 
#since there are four diel categories, anything above 25% is better than chance?

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/all_conf_levels_confusion_matrix.pdf")
ggplot(concordance, aes(actual, predicted, fill = percent)) + geom_tile() + geom_text(aes(label = percent)) +
  scale_fill_gradient(low = "white", high = "dodgerblue") + labs(x = "Actual", y = "Predicted") 
dev.off()

#function to plot the concordance for each of the confidence levels
plotConcordance = function(set_column = "Conf2"){
  diel_full_filtered <- diel_full %>% filter(column == set_column)
  #need to filter for species with more than one entry or else concordance will always be 100%
  mulitple_sources <- diel_full_filtered %>% count(Species_name) %>% filter(n>1)
  diel_full_filtered <- diel_full_filtered[diel_full_filtered$Species_name %in% mulitple_sources$Species_name,]
  concordance <- as.data.frame(table(diel_full_filtered$max_crep, diel_full_filtered$value))
  colnames(concordance) <- c("actual", "predicted", "freq")
  totals_df <- aggregate(concordance$freq, by=list(Category=concordance$actual), FUN=sum)
  colnames(totals_df) <- c("actual", "total")
  concordance <- merge(concordance, totals_df, by = "actual")
  concordance$percent <- round(concordance$freq / concordance$total * 100, 1)
  return(concordance)
}

concordance_list <- lapply(sort(unique(diel_full$column)), function(x){plotConcordance(x)})

for(i in seq_along(concordance_list)){
  pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "confidence", i, "_confusion_matrix.pdf"))
  print(ggplot(as.data.frame(concordance_list[i]), aes(actual, predicted, fill = percent)) + geom_tile() + geom_text(aes(label = percent)) +
    scale_fill_gradient(low = "white", high = "dodgerblue") + labs(x = "Actual", y = "Predicted") +
    ggtitle(paste("Confidence level ", i, " concordance"))) 
  dev.off()
}

# Section 4.5: Concordance between confidence levels ----------------------
diel_full_long <- read.csv(here("cetacean_confidence_long_df.csv"))
diel_full_long[diel_full_long == "unclear"] <- NA
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

#use below for four state data, will also need to make a call on partially cathemeral species
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "diurnal/crepuscular", replacement = "diurnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "nocturnal/cathemeral", replacement = "cathemeral")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "diurnal/cathemeral", replacement = "cathemeral")

unique(diel_full_long$value)

#remove unclear as an option, since we aren't interested in if unclear matches with unclear
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/nocturnal", replacement = "nocturnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/diurnal", replacement = "diurnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/crepuscular", replacement = "crepuscular")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear/cathemeral", replacement = "crepuscular")

#with only species included in the final tree. From 83 to 75.
diel_full_long <- diel_full_long[diel_full_long$tips %in% mam.tree$tip.label,]

#get a list of all the species with more than one source (should be most of them)
species_list <- table(diel_full_long$Species_name)
#should be 76 species with all cetaceans, x with cetaceans 
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

#species <- "Balaenoptera acutorostrata"

#apply this function across all species with multiple entries
output <- lapply(species_list, function(species) {
  
  #filter for one species at a time
  df <- diel_full_long[diel_full_long$Species_name == species,]
  #rename the column names to be unique for every entry for this species (ie for multiple column 2 entries column 2.1, 2.2 etc)
  df$column <- make.unique(df$column)
  
  #converts the dataframe so it compares every entry with each other (ie for A,B,C A-A, A-B, A-C, B-A, B-B, B-C, etc)
  df_lists_comb <- expand(df, nesting(var = column, vector = value), nesting(var2 = column, vector2 = value), .name_repair = "universal")
  
  #??? idk
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
plot_freq <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = round(Freq, digits = 2))) + 
  geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_freq
plot_count <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = count)) +
  geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_count

# Section 5: Cetacean confidence sankey ----------------------------------
diel_full <- read.csv(here("cetacean_confidence_wide.csv"))

#use below for all columns
df <- diel_full %>% make_long(6:39)

#use below for custom column selection
#df <- diel_full %>% make_long(Conf1.1, Conf2.1, Conf3.1, Conf4.1, Conf5.1)

#for all the columns with level 5 data, how does it agree with other data sources?
#df <- diel_full %>% make_long(Conf3.1, Conf3.2, Conf3.3, Conf3.4, Conf3.5)
#df <- diel_full %>% make_long(Conf5.1, Conf5.2, Conf5.3, Conf5.4)
#df <- diel_full %>% make_long(Conf4.1, Conf4.2, Conf4.3, Conf4.4, Conf4.5)
#df <- diel_full %>% make_long(Conf3.1, Conf4.1, Conf5.1)

df$node <- tolower(df$node)
df$node <- str_trim(df$node)
unique(df$node)
#calling these nocturnal, diurnal, cathemeral (removing unclear as an option)
df$node <- str_replace_all(df$node, pattern = "unclear/nocturnal", replacement = "nocturnal")
df$node <- str_replace_all(df$node, pattern = "unclear/diurnal", replacement = "diurnal")
df$node <- str_replace_all(df$node, pattern = "unclear/cathemeral", replacement = "cathemeral")
df$node <- str_replace_all(df$node, pattern = "unclear/crepuscular", replacement = "cathemeral")
df$node <- str_replace_all(df$node, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

df$next_node <- tolower(df$next_node)
df$next_node <- str_trim(df$next_node)
unique(df$next_node)

#calling these nocturnal, diurnal, cathemeral (removing unclear as an option)
df$next_node <- str_replace_all(df$next_node, pattern = "unclear/nocturnal", replacement = "nocturnal")
df$next_node <- str_replace_all(df$next_node, pattern = "unclear/diurnal", replacement = "diurnal")
df$next_node <- str_replace_all(df$next_node, pattern = "unclear/cathemeral", replacement = "cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "unclear/crepuscular", replacement = "cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

df[df == "unclear"] <- NA
df[df == ""] <- NA

#remove entries where both the starting and ending node are NA
df <- df %>% filter_at(vars(node, next_node), any_vars(!is.na(.)))

#optional: for four states
df$node <- str_replace_all(df$node, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
df$node <- str_replace_all(df$node, pattern = "diurnal/crepuscular", replacement = "diurnal")
df$next_node <- str_replace_all(df$next_node, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
df$next_node <- str_replace_all(df$next_node, pattern = "diurnal/crepuscular", replacement = "diurnal")
df$node <- str_replace_all(df$node, pattern = "diurnal/cathemeral", replacement = "cathemeral")
df$node <- str_replace_all(df$node, pattern = "nocturnal/cathemeral", replacement = "cathemeral")
#not sure whether to call these cathemeral or di/noc
df$next_node <- str_replace_all(df$next_node, pattern = "diurnal/cathemeral", replacement = "cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "nocturnal/cathemeral", replacement = "cathemeral")

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) + geom_sankey() + theme_sankey(base_size = 16)


# Section 5: Artiodactyl   ------------------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1JGC7NZE_S36-IgUWpXBYyl2sgnBHb40DGnwPg2_F40M/edit?gid=562902012#gid=562902012'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#filter for entries that have an activity pattern, I'm ignoring the 7th source because only two species have seven sources
diel_full <- diel_full[!is.na(diel_full$Confidence_primary_source), c(1, 2, 6:17)]
#filter for species with at least two sources so we can compare 
diel_full <- diel_full %>% filter(Diel_Pattern_secondary != "")
diel_full$Diel_Pattern_primary <- tolower(diel_full$Diel_Pattern_primary)
diel_full$Diel_Pattern_secondary <- tolower(diel_full$Diel_Pattern_secondary)
diel_full$Diel_Pattern_tertiary <- tolower(diel_full$Diel_Pattern_tertiary)
diel_full$Diel_pattern_4th_source <- tolower(diel_full$Diel_pattern_4th_source)
diel_full$Diel_Pattern_5th_source <- tolower(diel_full$Diel_Pattern_5th_source)
diel_full$Diel_Pattern_6th_source <- tolower(diel_full$Diel_Pattern_6th_source)

diel_full <- diel_full %>% pivot_longer(cols = c(Confidence_primary_source, Confidence_secondary_source, Confidence_tertiary_source, Confidence_4th_source, Confidence_5th_source, Confidence_6th_source), names_to = "column", values_to = "values")

#for each species we can see what confidence level the primary, secondary, tertiary sources are
diel_full$Confidence_level <- paste(diel_full$column, diel_full$values, sep = "_")

#filter for each source (primary -4th) to match each entry with its confidence level
diel_full_1 <- diel_full %>% filter(column == "Confidence_primary_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_primary) 
diel_full_1 <- diel_full_1[, c(1, 2, 10:14)]
diel_full_2 <- diel_full %>% filter(column == "Confidence_secondary_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_secondary)
diel_full_2 <- diel_full_2[, c(10:14)]
diel_full_3 <- diel_full %>% filter(column == "Confidence_tertiary_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_tertiary)
diel_full_3 <- diel_full_3[, c(10:14)]
diel_full_4 <- diel_full %>% filter(column == "Confidence_4th_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_pattern_4th_source)
diel_full_4 <- diel_full_4[, c(10:14)]
diel_full_5 <- diel_full %>% filter(column == "Confidence_5th_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_5th_source)
diel_full_5 <- diel_full_5[, c(10:13)]
diel_full_6 <- diel_full %>% filter(column == "Confidence_6th_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_6th_source)
diel_full_6 <- diel_full_6[, c(10:12)]
diel_full <- cbind(diel_full_1, diel_full_2, diel_full_3, diel_full_4, diel_full_5, diel_full_6)
#rename the columns to the confidence level, R will add .1 for every duplicate so they will end up with unique identifiers
colnames(diel_full) <- c("Species_name", "Family", "Conf5", "Conf4", "Conf3", "Conf2", "Conf1", "Conf3", "Conf4", "Conf2", "Conf1", "Conf5", "Conf3", "Conf4", "Conf1", "third_source_NA", "Conf2", "Conf3", "Conf4", "fourth_source_NA", "Conf1", "Conf2", "fifth_source_NA", "Conf3", "Conf5", "Conf4", "sixth_source_NA", "Conf3", "Conf4")
#drop the two columns with only NA values ("tertiary_source_NA", "fourth_source_NA", etc)
diel_full <- diel_full[, -c(16, 20, 23, 27)]

diel_full <- diel_full[,order(colnames(diel_full))]
diel_full <- diel_full %>% relocate("Species_name", .before = "Conf1")
diel_full <- diel_full %>% relocate("Family", .after = "Species_name")

diel_full[diel_full == ""] <- NA

#replace strings
#may change to cathemeral/nocturnal/crepuscular in future, 
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("diurnal/weakly-crepuscular", "diurnal/crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal/crepuscular", "nocturnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-nocturnal", "nocturnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-diurnal/crepuscular", "diurnal/crepuscular", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-diurnal", "diurnal/cathemeral", x)}))
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("weakly-cathemeral/crepuscular", "cathemeral/crepuscular", x)}))
#check what this is
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("diurnal-maybe", "diurnal", x)}))
#and this
diel_full <- data.frame(lapply(diel_full, function(x) {gsub("unclear/diurnal", "diurnal", x)}))

write.csv(diel_full, here("confidence_artio_wide.csv"), row.names = FALSE)


# Section 6: Artiodactyla long format -----------------------------------
#read in wide data
diel_full <- read.csv(here("confidence_artio_wide.csv"))

#convert into long format
diel_full_long <- diel_full %>% pivot_longer(cols = c(3:25), names_to = "column", values_to = "value")

diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]
unique(diel_full_long$value)

#replace strings, can check to see later if these change anything
diel_full_long$value <- str_replace(diel_full_long$value, pattern = "nocturnal/cathemeral", replacement = "cathemeral")
diel_full_long$value <- str_replace(diel_full_long$value, pattern = "diurnal/cathemeral", replacement = "cathemeral")

#check for any unconventional strings
unique(diel_full_long$value)

write.csv(diel_full_long, here("confidence_artio_long.csv"), row.names = FALSE)

#separate out crepuscularity into its own column
#confidence level 2 data will be unclear/crepuscular since they don't show evidence of nocturnal or diurnal activity
diel_full_long <- separate(diel_full_long, col = value, into = c("new_diel", "crepuscular"), sep = "/")

#set unclear to NA since it does not provide any additional information
diel_full_long[diel_full_long == "unclear"] <- NA

#run each species through this function, x species with activity pattern data (di, noc or cath)
activity_pattern_df <- diel_full_long[!is.na(diel_full_long$new_diel),] %>% group_by(Species_name) %>% do(tabulated_diel_pattern = tabulateFunc2(.)) %>% unnest()

#now determine whether or not each species is crepuscular
diel_full_long[is.na(diel_full_long)] <- "0"
diel_full_long$crepuscular <- str_replace(diel_full_long$crepuscular, pattern = "crepuscular", replacement = "1")
diel_full_long$crepuscular <- as.numeric(diel_full_long$crepuscular)
diel_full_long$total <- 1

#can't determine crepuscular activity from level or level 5 data. So manually fix any species with only level 1 and/or 5 data

crep_df <- diel_full_long %>% group_by(Species_name) %>% do(tabulated_crep = tabulateCrep(.)) %>% unnest()

test <- merge(crep_df, activity_pattern_df, by = "Species_name")
test$tabulated_diel <- test$tabulated_diel_pattern
for(i in 1:nrow(test)){
  if(test[i, "tabulated_crep"] == TRUE){
    test[i, "tabulated_diel"] <- paste(test[i, "tabulated_diel_pattern"], "crepuscular", sep = "/")
  }
}

unique(test$tabulated_diel)

test <- test[, c("Species_name", "tabulated_diel")]
#save out the new tabulated activity pattern dataframe
write.csv(test, here("artio_tabulated_full.csv"), row.names = FALSE)


# Section 7: plot to compare with existing data ------------------------
diel_full_new <- merge(diel_full, test, by = "Species_name")

#diel pattern 2 is all diel categories, diel pattern is maxcrep, new pattern is manually determined diel patterns, tabulated pattern is diel patterns determined from above function
diel_full_new <- diel_full_new[, c("Species_name", "Parvorder", "Diel_Pattern_2", "Diel_Pattern","New_Pattern", "Confidence", "tips", "tabulated_diel")]


# Section 8: Plot new vs old diel pattern comparison ------------------------

trait.data <- diel_full_new[diel_full_new$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue", "grey","pink")
custom.colours <- c("#dd8ae7","#FC8D62", "yellow", "#66C2A5", "green")
custom.colours.2 <- c("#dd8ae7", "peachpuff2", "#FC8D62","yellow", "#66C2A5", "green","black","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "tabulated_diel")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = tabulated_diel), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Tabulated activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

custom.colours <- c("white", "#dd8ae7",  "peachpuff2", "#FC8D62", "yellow", "#66C2A5", "green", "black", "blue", "red")
custom.colours.2 <- c("#dd8ae7", "peachpuff2", "#FC8D62","yellow", "#66C2A5", "green", "black", "grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "New_Pattern", "tabulated_diel")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "New activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = tabulated_diel), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Tabulated activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

#with max crep designations
trait.data$tabulated_diel <- str_replace(trait.data$tabulated_diel, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$tabulated_diel <- str_replace(trait.data$tabulated_diel, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$tabulated_diel <- str_replace(trait.data$tabulated_diel, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
custom.colours.2 <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "tabulated_diel")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = tabulated_diel), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Tabulated activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

# #Section 9: Concordance ---------------------------------------------
diel_full_long <- read.csv(here("confidence_artio_long.csv"))

unique(diel_full_long$value)

#remove unclear since it gives no new information
diel_full_long$value <- str_replace(diel_full_long$value, pattern = "unclear/crepuscular", replacement = "crepuscular")
diel_full_long[diel_full_long == "unclear"] <- NA
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

species_list <- table(diel_full_long$Species_name)
#should be 116 species
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
  
  #??? idk
  df_lists_comb <- df_lists_comb %>% filter(var != var2) %>% arrange(var, var2) %>% mutate(vars = paste0(var, ".", var2)) %>% select(contains("var"), everything())
  
  #evaluates the activity patterns for each of these sources and returns if they agree or not (TRUE or FALSE)
  comparisons <- df_lists_comb %>% group_by(vars) %>% mutate(comp = compTwo(comp1 = vector, comp2 = vector2))
  #manipulate the strings for both variable names to revert them back to the original name (back to column 2 from col 2.1)
  comparisons$var <- str_sub(comparisons$var, start = 1, end = 5)
  comparisons$var2 <- str_sub(comparisons$var2, start = 1, end = 5)
  
  #create a column returning the comparison being made (ie col2-col2, col1-col2, etc)
  comparisons$var_final <- paste(comparisons$var, comparisons$var2, sep = "-")
  
  #return just the comparison result column (TRUE or FALSE match) and the comparison being made (ie col1 vs col1)
  return(comparisons[,c("comp","var_final")])
})

#combine this list of results 
output <- Reduce(rbind, output)

table <- table(output$var_final)
# prop.table(table, margin = 1)
table2 <- as.data.frame(prop.table(table(output$var_final, output$comp), margin = 1))
table2$Comp1 <- sapply(str_split(table2$Var1, "-"), `[`, 1)
table2$Comp2 <- sapply(str_split(table2$Var1, "-"), `[`, 2)
table2 <- table2[table2$Var2 == TRUE,]
table2$count <- table
plot_freq <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = round(Freq, digits = 2))) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_freq
plot_count <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = count)) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_count


# Section 7: Plot manual data on tree --------------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#diel_full <- diel_full %>% filter(Parvorder == "Mysticeti")

diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")
diel_full$Diel_Pattern <- tolower(diel_full$Diel_Pattern)
diel_full$New_Pattern <- tolower(diel_full$New_Pattern)
diel_full <- diel_full[!is.na(diel_full$Diel_Pattern),]

diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "cathemeral-variable", replacement = "cathemeral")
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "cathemeral-invariate", replacement = "cathemeral")
#calling these nocturnal, diurnal, cathemeral (removing unclear as an option)
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "unclear-nocturnal", replacement = "nocturnal")
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "unclear-diurnal", replacement = "diurnal")
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "unclear-cathemeral", replacement = "cathemeral")
#max_crep dataset, optional
# diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
# diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
# diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")

#calling these nocturnal, diurnal, crepuscular, etc may switch to cathemeral later
#idk what to do about this one
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "weakly-nocturnal/cathemeral", replacement = "cathemeral")
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "weakly-diurnal", replacement = "cathemeral")
#diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "weakly-crepuscular", replacement = "crepuscular/cathemeral")
diel_full$New_Pattern <- str_replace_all(diel_full$New_Pattern, pattern = "weakly-nocturnal", replacement = "cathemeral")

trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label, ]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue", "grey","pink")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
custom.colours.2 <- c("grey", "#dd8ae7", "peachpuff2", "#FC8D62", "yellow", "#66C2A5", "green", "black")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "New_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "New activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue", "grey","pink")
custom.colours <- c("grey", "#dd8ae7", "peachpuff2", "#FC8D62", "#66C2A5", "black")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "New_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/new_cetacean_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()