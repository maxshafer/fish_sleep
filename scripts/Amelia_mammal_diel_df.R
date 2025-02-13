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

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: How well do the Bennie and Maor dfs agree with each --------

#read in the Bennie diel activity patterns
#from https://doi.org/10.1073/pnas.1216063110 
Bennie_mam_data <- read_excel(here("Bennie_diel_activity_data.xlsx"))
colnames(Bennie_mam_data) <- "SpeciesBehaviourReference"
Bennie_mam_data$SpeciesBehaviourReference <- str_replace(string = Bennie_mam_data$SpeciesBehaviourReference, pattern = " ", replacement  = "_")
Bennie_mam_data <- separate(Bennie_mam_data, col = SpeciesBehaviourReference, into = c("Species", "Activity_pattern", "Reference"), sep = " ")

#read in the Maor diel actiivty patterns
#from https://doi.org/10.1038/s41559-017-0366-5 
maor_mam_data <- read_excel(here("Maor_diel_activity_data.xlsx"))
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 3:5]
colnames(maor_mam_data) <- c("Species", "Activity_pattern", "Reference")

#the issue with this df is that if they have an alternative activity they add it in a new row
duplicates1 <- maor_mam_data[duplicated(maor_mam_data$Species),]
#make another dataframe since some sps are repeated twice
duplicates2 <- duplicates1[duplicated(duplicates1$Species),] #1080 species have at least one alt diel pattern
duplicates1 <- duplicates1[!duplicated(duplicates1$Species),] #208 species have 2 alt diel patterns

#remove duplicates from Maor dataframe for now, we'll add alternative diel patterns back next
maor_mam_data <-maor_mam_data[!duplicated(maor_mam_data$Species),]
maor_mam_data$Species <- str_replace(string = maor_mam_data$Species, pattern = " ", replacement  = "_")

#format some entries to be max_crep (di/crep -> crep, noc/crep -> crep)
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal/ Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal/Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal /Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal or Cathemeral", "Cathemeral")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal/Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal - EXTINCT", "Nocturnal")

#Don't add back alternative patterns for now
# #add all the extra diel patterns, then sort the columns after since they're in a random order anyway
# maor_mam_data <- merge(maor_mam_data, duplicates1, by='Species', all.x = TRUE, all.y = TRUE)
# maor_full <- merge(maor_mam_data, duplicates2, by='Species', all.x = TRUE, all.y = TRUE)
# maor_full <- maor_full[, c("Species", "Activity_pattern", "Activity_pattern.x", "Activity_pattern.y")]
# maor_full <- relocate(maor_full, "Activity_pattern.x", .after = "Activity_pattern.y")
# colnames(maor_full) <- c("Species", "alt_pattern_1", "alt_pattern_2", "Activity_pattern")


# Section 2: How well do these sources agree? -----------------------------
diel_merge <- merge(Bennie_mam_data,maor_mam_data,by="Species")
colnames(diel_merge) <- c("Species", "Bennie_diel", "Bennie_source", "Maor_diel", "Maor_source")
diel_merge$Bennie_diel <- tolower(diel_merge$Bennie_diel)
diel_merge$Maor_diel <- tolower(diel_merge$Maor_diel)
#set default to idk
diel_merge$match <- "Idk"

for(i in 1:length(diel_merge$Species)){
  if(diel_merge[i, "Bennie_diel"] == diel_merge[i, "Maor_diel"]){
    diel_merge[i, "match"] <- "Yes"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal") & diel_merge[i, "Maor_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")& diel_merge[i, "Maor_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("cathemeral", "cathemeral/crepuscular")& diel_merge[i, "Maor_diel"] %in% c("cathemeral", "cathemeral/crepuscular")){
    diel_merge[i, "match"] <- "Approximate"
  } else {
    diel_merge[i, "match"] <- "No"
  }
}

table(diel_merge$match)

#A surprising amount of consistency!
#170/2199 approximate matches, 244/2199 no matches, and 1785/2199 yes matches
#81% match!
diel_table <- diel_merge %>% count(match)
diel_table <- transform(diel_table, percent = (n/sum(diel_table$n)) * 100)

diel_table <- trait.data.all %>% count(Confidence)
diel_table <- transform(diel_table, percent = (n/sum(diel_table$n)) * 100)

ggplot(diel_table, aes(x="", y=n, fill=Confidence)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + theme_void()

#which species did not match
mismatch_species <- diel_merge %>% filter(match == "No")
table(mismatch_species$Bennie_diel, mismatch_species$Maor_diel)
#we should clean some of these entries up, could be a formatting error in why they don't match

#plot this on the mammal tree, which species tend to have conflicting data?
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- mismatch_species[,1:2]
trait.data <- trait.data[trait.data$Species %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$Species)
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)

#seem to be distributed throughout the tree

#save out the final 
colnames(diel_merge) <- c("tips", "Bennie_diel", "Bennie_source", "Maor_diel", "Maor_source", "match")
write.csv(diel_merge, here("sleepy_mammals.csv"))

# Section 3: Visualizing the mammal data -------------------------------
trait.data <- diel_merge
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

trait.data <- trait.data[trait.data$Species %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$Species)

diel.plot.all <- ggtree(trpy_n_mam, layout = "circular") %<+% trait.data[,c("Species", "Bennie_diel")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_mam$tip.label),], aes(x=x, y=y, fill = Bennie_diel), inherit.aes = FALSE, colour = "black", width = 5)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 1, offset = 3)
diel.plot.all