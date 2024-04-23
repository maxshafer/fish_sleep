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

setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")


# Section 1: Formatting the cetacean diel dataframe -----------------------

###### For cetaceans

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit#gid=0'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#set all diel pattern entries to lower case, helps keep consistency
cetaceans_full$Diel_Pattern_1 <- tolower(cetaceans_full$Diel_Pattern_1)
cetaceans_full$Diel_Pattern_2 <- tolower(cetaceans_full$Diel_Pattern_2)
cetaceans_full$Diel_Pattern_3 <- tolower(cetaceans_full$Diel_Pattern_3)

#they use a different name for sperm whales in the mam tree so we can change that
cetaceans_full$Species_name <- str_replace(cetaceans_full$Species_name, "Physeter catodon", "Physeter macrocephalus")
cetaceans_full <- cetaceans_full[,1:10]

#cycle through diel pattern 1, 2, 3 to check that all entries are correctly spelled and see spread of data
#20 cath sps, 20 crep, 15 di, 25 noc
table(cetaceans_full$Diel_Pattern_3)

#add another row "tips" with the species name formatted as they appear in the tree (mam tree and otl)
#saves you from remaking it for each of the trait.data dataframes
cetaceans_full$tips <- cetaceans_full$Species_name
cetaceans_full$tips <- str_replace(cetaceans_full$tips, pattern = " ", replacement = "_")

#rename the row names to be the tip names so it's easier to subset by the tree tip labels later
row.names(cetaceans_full) <- cetaceans_full$tips

## Probably should save out a local copy in case google goes bankrupt
write.csv(cetaceans_full, file = here("cetaceans_full.csv"))

# Section 3a: Formatting the artiodactyla diel dataframe ??? data (only contains di/noc)-----------------------

#we still need to load in the cetacean dataframe
#will combine this into one dataframe with the rest of artiodactyla

# Load in the cetacean dataframe
cetaceans_full <- read.csv("cetaceans_full.csv")
## Remove species without diel data
cetaceans_full <- cetaceans_full[!(is.na(cetaceans_full$Diel_Pattern_2)),]

#load in mammal trait data (diel activity patterns)
mammals_full <- readRDS("trait_data_mammals.rds")

#subset for just cetartiodactyla
artiodactyla <- subset(mammals_full, Order == "Cetartiodactyla")
# This data doesn't find any artio species to be cathemeral or crepuscular. 
#Which isn't accurate. This is just a binary dataset
# So all three diel columns are interchangeable

#edit cetacean dataframe so it has similar format to add it to artio dataframe
# diel 1 is only di or noc
# diel 2 is di, noc, di/crep, noc/crep, cathemeral
# diel 3 is di, noc, crep, cath
#since there are no crepuscular or cathemeral artiodactyla can match these up with any of their columns

cetaceans_full <- subset(cetaceans_full, select = c(Species_name, Diel_Pattern_1, Diel_Pattern_2, Diel_Pattern_3, tips))
artiodactyla <- artiodactyla[, c(1, 5, 6, 7)]
colnames(artiodactyla) <- c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3")
artiodactyla$tips <- str_replace(artiodactyla$Species_name, " ", "_")

#add cetaceans data to the rest of artiodactyla
artiodactyla_full <- rbind(artiodactyla, cetaceans_full)
row.names(artiodactyla_full) <- artiodactyla_full$Species_name

#save out as CSV file
write.csv(artiodactyla_full, here("artiodactyla_binary_df.csv"))


# Section 3b: Formatting the artiodactyla diel dataframe Cox data  --------



# Section 3c: Formatting the artiodactyla diel dataframe Maor data  --------
#from https://doi.org/10.1038/s41559-017-0366-5
maor_mam_data <- read_excel(here("Maor_diel_activity_data.xlsx"))
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 1:4]
colnames(maor_mam_data) <- c("Order", "Family", "Species", "Activity_pattern_1")

#subset just for artiodactyla
maor_mam_data <- subset(maor_mam_data, maor_mam_data$Order == "Artiodactyla")

#if a species has an alternative diel pattern its added in a new row :[
# find how many duplicated rows there are
dim(maor_mam_data)
#there are 219 species entries
length(unique(maor_mam_data$Species))
#only 165 are unique
#make a dataframe of all the duplicated species 
duplicates1 <- maor_mam_data[duplicated(maor_mam_data$Species),]
#make another dataframe since some sps are repeated twice
duplicates2 <- duplicates1[duplicated(duplicates1$Species),] #44 species have at least one alt diel pattern
duplicates1 <- duplicates1[!duplicated(duplicates1$Species),] #10 species have 2 alt diel patterns

#remove duplicates from Maor dataframe for now, we'll add alternative diel patterns back next
maor_mam_data <-maor_mam_data[!duplicated(maor_mam_data$Species),]

#diel pattern 1 = di/noc strictly
#diel pattern 2 = di, noc, di/crep, noc/crep, cath and crep
#diel pattern 3 = di, noc, cath, maximize for crep

#add all the extra diel patterns, then sort the columns after since they're in a random order anyway
maor_mam_data <- merge(maor_mam_data, duplicates1, by='Species', all.x = TRUE, all.y = TRUE)
maor_full <- merge(maor_mam_data, duplicates2, by='Species', all.x = TRUE, all.y = TRUE)
maor_full <- maor_full[, c("Species", "Activity_pattern_1", "Activity_pattern_1.x", "Activity_pattern_1.y")]
maor_full <- relocate(maor_full, "Activity_pattern_1.x", .after = "Activity_pattern_1.y")
#create column for activity pattern 2
#diel pattern 2 includes all variation, di, noc, di/crep, noc/crep, cath, crep as reported
colnames(maor_full) <- c("Species", "alt_pattern_1", "alt_pattern_2", "Diel_Pattern_2")

#create column for activity pattern 1
# diel pattern 1 maximizes for only diurnal and nocturnal
maor_full$Diel_Pattern_1 <- maor_full$Diel_Pattern_2
for(i in 1:nrow(maor_full)){
  if(maor_full[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_1"] <- "Diurnal"
  } else if(maor_full[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_1"] <- "Nocturnal"
  } else if(maor_full[i, "Diel_Pattern_2"] %in% c("Crepuscular", "Cathemeral")){
    maor_full[i, "Diel_Pattern_1"] <- NA
  }
}

#create activity pattern 3 column
# diel pattern 2 maximizes for crepuscularity
maor_full$Diel_Pattern_3 <- maor_full$Diel_Pattern_2
for(i in 1:nrow(maor_full)){
  if(maor_full[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if(maor_full[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_3"] <- "Crepuscular"
  } 
}

maor_full <- relocate(maor_full, "Diel_Pattern_2", .after = "Diel_Pattern_1")
#now add in the alternative diel patterns 

#method A: if alt pattern doesn't match current diel pattern check the literature
#method B: if alt pattern matches current pattern replace with NA (alt columns will be deleted after)
#if alternative pattern is crepuscular, add it as a secondary diel pattern

#Method A: 35 species with unresolved diel patterns
#now to decide on the 35 species with conflicting data
#Maor et al gets actiivty pattern data from four sources, 57, 4, 33 and 29
#Source 57: 57.Jones, K.E. et al., 2009. PanTHERIA: A species-level database of life history, ecology and geography of extant and recently extinct mammals. Ecology, 90, p.2648.
#Source 4: Aulagnier, S. & Thevenot, M., 1986. Catalogue des Mammiferes Sauvages du Maroc, Rabat, Morocco: Institute Scientifique.
#Source 33: Hufnagl, E., 1972. Lybian Mammals, Harrow, England: The Oleander Press.
# Source 29: 29. Gray, G.G. & Simpson, C.D., 1980. Ammotragus lervia. Mammalian Species, 144, pp.1–7.

# Sources 4, 29, 33 are all for  Ammotragus lervia and from the 1970s/1980s. 
# https://doi.org/10.25225/jvb.20055 more recent study from 2020 utilizing 24h camera traps found that A lervia is crepuscular and diurnal
#edit maor_full to reflect this
maor_full[6,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"
maor_full[6,"Diel_Pattern_1"] <- "Diurnal"
maor_full[6, "Diel_Pattern_3"] <- "Crepuscular"
maor_full[6,"alt_pattern_2"] <- NA

#All other 34 species have conflicting data from panTHERIA database (source 57)
#will examine these separatelyon a one by one basis and update the maor_full dataset
#import dataframe of updated information on each of these species and match with their current row
#uncomment lines below

# missing_sps <- read.csv("filepath/missing_sps.csv")
# row.names(missing_sps) <- missing_sps$Species
# 
# replacements <- match(missing_sps$Species, maor_full$Species)
# for(i in replacements){
#   maor_full[i, "Diel_Pattern_2"] <- missing_sps[maor_full[i, "Species"], "Diel_Pattern_2"]
# }

#Method B: Matches or near matches 
#Alcelaphus lichtensteinii
maor_full[4,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"
maor_full[4,"Diel_Pattern_3"] <- "Crepuscular"
maor_full[4,"alt_pattern_2"] <- NA

#Ammotragus lervia
maor_full[6,"alt_pattern_1"] <- NA

#Axis porcinus
maor_full[13,"alt_pattern_2"] <- NA

#Blastocerus dichotomus
maor_full[17,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"
maor_full[17,"Diel_Pattern_3"] <- "Crepuscular"
maor_full[17, "alt_pattern_2"] <- NA

#Bubalus bubalis
maor_full[23,"alt_pattern_2"] <- NA

#Bubalus mindorensis
maor_full[25,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[25,"Diel_Pattern_1"] <- "Nocturnal"
maor_full[25,"alt_pattern_2"] <- NA

#Cervus elaphus. Row 52 duplicated because there are four entries for it
#no new information in row 52 so drop it and renumber
maor_full <- maor_full[-c(52), ]
row.names(maor_full) <- 1:nrow(maor_full)

#Gazella subgutturosa
maor_full[63,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[63,"Diel_Pattern_1"] <- "Nocturnal"
maor_full[63,"alt_pattern_2"] <- NA

#Giraffa camelopardalis
maor_full[64,"alt_pattern_2"] <- NA

#Moschus moschiferus
maor_full[92,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[92,"Diel_Pattern_1"] <- "Nocturnal"
maor_full[92,"alt_pattern_2"] <- NA

#Philantomba maxwellii
maor_full[123,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"
maor_full[123,"Diel_Pattern_1"] <- "Diurnal"
maor_full[123,"alt_pattern_2"] <- NA

#Rangifer tarandus
maor_full[131,"alt_pattern_2"] <- NA

#Sus barbatus
maor_full[145,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[145,"Diel_Pattern_3"] <- "Crepuscular"
maor_full[145,"alt_pattern_2"] <- NA

#Tragelaphus angasii
maor_full[156,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"
maor_full[156,"Diel_Pattern_3"] <- "Crepuscular"
maor_full[156,"alt_pattern_2"] <- NA


#now we can format the artiodactyla data so it matches with the cetacean data and merge them
maor_full <- maor_full[, c("Species", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3")]
maor_full$Diel_Pattern_1 <- tolower(maor_full$Diel_Pattern_1)
maor_full$Diel_Pattern_2 <- tolower(maor_full$Diel_Pattern_2)
maor_full$Diel_Pattern_3 <- tolower(maor_full$Diel_Pattern_3)

maor_full$tips <- str_replace(maor_full$Species, " ", "_")
colnames(maor_full) <- c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "tips")

cetacaceans <- read.csv("cetaceans_full.csv")
cetaceans <- cetaceans_full[, c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "tips")]
maor_full <- rbind(maor_full, cetaceans)
row.names(maor_full) <- c(maor_full$Species_name)

write.csv(maor_full, here("artiodactyla_full.csv"))

# Section 4: Formatting the artiodactyla_without_cetaceans diel dataframe -----------------------

#easiest way to do this is to drop the last chunk of rows in maor_full (to keep all the formatting)
maor_full <- read.csv("artiodactyla_full.csv")
just_artio <- maor_full[1:151, ]
colnames(just_artio) <- c("tips", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3")
write.csv(just_artio, here("artiodactyla_without_cetaceans.csv"))


# Section 5: Load in and examine the mam tree -----------------------------------------

## Read in the mammalian phylogeny
#mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
#mam.tree <- maxCladeCred(mammal_trees, tree = TRUE) 
#this function takes a long time to run so save the result out

# Save out maxcladecred, so we don't have to recalculate it every time
#saveRDS(mam.tree,"maxCladeCred_mammal_tree.rds")

mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")

#check to see what cetaceans are in mammal tree by subsetting to just include cetaceans
#use cetaceans_full because it has all species
#trpy_n_test <- keep.tip(mam.tree, tip = cetaceans_full$tips) -uncomment to see missing species
#mam tree is missing 12 species
#described in 2014: Sousa_plumbea, Sousa_sahulensis, Mesoplodon_hotaula
#described in 2019: Berardius_minimus, Inia_araguaiaensis, 
#described in 2021: Balaenoptera_ricei, Mesoplodon_eueu, Platanista_minor
#subspecies, not a valid species: Inia_humboldtiana, Neophocaena_sunameri, Inia_boliviensis, Balaenoptera_brydei

#we need to drop the missing species to create the tree
cetaceans_full <- read.csv("cetaceans_full.csv")
cetaceans_full <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n_test <- keep.tip(mam.tree, tip = cetaceans_full$tips)
#check if time calibrated
test <- ggtree(trpy_n_test, layout = "circular", size = 0.5) + geom_tiplab(size = 1.5)
#test$data #29.36 million years ago to the root. Seems correct!

#look at the tree for whippomorpha (hippos + cetacea) and see if time calibrated
whippomorpha <- cetaceans_full
whippomorpha <- rbind(cetaceans_full, c("Hexaprotodon liberiensis", "", NA, NA, "Nocturnal", "Nocturnal", "Nocturnal", "", "3", "nocturnal", "Choeropsis_liberiensis"))
whippomorpha <- rbind(whippomorpha, c("Hippopotamus amphibius", "", NA, NA, "Nocturnal", "Nocturnal/crepuscular", "Crepuscular", "", "4", "crepuscular", "Hippopotamus_amphibius"))
write.csv(whippomorpha, file = here("whippomorpha.csv"))
whippomorpha <- whippomorpha[whippomorpha$tips %in% mam.tree$tip.label,]
trpy_n_whippo <- keep.tip(mam.tree, tip = whippomorpha$tips)
whippo <- ggtree(trpy_n_whippo, layout = "circular", size = 0.5) + geom_tiplab(size = 1.5)
#whippo$data shows the branch lengths (53.7 my to root). Also seems correct!
animals <- read_excel(here("Maor_diel_activity_data.xlsx"))
colnames(animals) <- c("1", "order", "tips", "4", "5", "6")
mustellidae <- animals[animals$order == "Mustelidae", ]
mustellidae <- mustellidae[!(is.na(mustellidae$order)),]
mustellidae$tips <- str_replace(mustellidae$tips, pattern = " ", replacement = "_")

mustellidae <- mustellidae[mustellidae$tips %in% mam.tree$tip.label, ]
trpy_must <- keep.tip(mam.tree, tip = mustellidae$tips)
tree <- ggtree(trpy_must, layout = "circular") + geom_tiplab(size = 2.5) #+ geom_text(aes(label=node))
tree

png('C:/Users/ameli/OneDrive/Documents/R_projects/mustellidae.png', width=21,height=20,units="cm",res=1200)
tree
dev.off()
