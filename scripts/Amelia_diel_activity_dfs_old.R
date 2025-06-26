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


# Section 1: Formatting the cetacean diel dataframe -----------------------

###### For cetaceans

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit?usp=sharing'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#save out full version with sources
write.csv(cetaceans_full, here("cetaceans_full_with_sources.csv"))

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

#drop na values
#cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$Diel_Pattern_2),]
#Alternatively, replace NA values
cetaceans_full <- cetaceans_full %>% replace_na(list(Diel_Pattern_2 = "unknown", Diel_Pattern_3 = "unknown", Confidence = "unknown"))


# transform diel pattern 3 to maximize for crepuscularity (nocturnal, diurnal and crepuscular trait states)
for(i in 1:nrow(cetaceans_full)){
  if(cetaceans_full[i, "Diel_Pattern_3"] == "cathemeral"){
    cetaceans_full[i, "Diel_Pattern_3"] <- "crepuscular"
  } 
}

#create new diel column for max_crep
# diel pattern 4 maximizes for diurnality and nocturnality (while keeping cathemerality as a trait state)
cetaceans_full$Diel_Pattern_4 <- cetaceans_full$Diel_Pattern_2
for(i in 1:nrow(cetaceans_full)){
  if(cetaceans_full[i, "Diel_Pattern_2"] == "diurnal/crepuscular"){
    cetaceans_full[i, "Diel_Pattern_4"] <- "diurnal"
  } else if(cetaceans_full[i, "Diel_Pattern_2"] == "nocturnal/crepuscular"){
    cetaceans_full[i, "Diel_Pattern_4"] <- "nocturnal"
  } else if (cetaceans_full[i, "Diel_Pattern_2"] == "cathemeral/crepuscular"){
    cetaceans_full[i, "Diel_Pattern_4"] <- "cathemeral"
  }
}

cetaceans_full <- cetaceans_full %>% relocate(Diel_Pattern_4, .after = Diel_Pattern_3)

#we want to take only the highest confidence level from all the sources
#create a column with the max confidence level for that species (out of the confidence level for all sources)
#the confidence values are characters so convert to numerics and then take the maximum value
cetaceans_full$Confidence <- lapply(cetaceans_full$Confidence, function(x) gsub(",", "\\1 ", x))
cetaceans_full$Confidence <- lapply(cetaceans_full$Confidence,function(x) strsplit(x, " ")[[1]])
cetaceans_full$max_conf <- lapply(cetaceans_full$Confidence, as.numeric)
cetaceans_full$max_conf <- lapply(cetaceans_full$Confidence, max)

#to fit with the artiodactyl dataset, rename max conf and drop the old confidence column
cetaceans_full$Confidence <- cetaceans_full$max_conf
cetaceans_full <- cetaceans_full[, -13]

#this allows you to save it out without an error 
cetaceans_full <- apply(cetaceans_full,2,as.character)

## Probably should save out a local copy in case google goes bankrupt
write.csv(cetaceans_full, file = here("cetaceans_full.csv"), row.names = FALSE)

#save out a version with hippos
whippomorpha <- read.csv(here("cetaceans_full.csv"))
whippomorpha <- rbind(whippomorpha, c("Hexaprotodon liberiensis", "", NA, NA, "nocturnal", "nocturnal/crepuscular", "crepuscular", "nocturnal", "NA", "3, 3, 3", "Hippopotimiade", "Choeropsis_liberiensis"))
whippomorpha <- rbind(whippomorpha, c("Hippopotamus amphibius", "", NA, NA, "nocturnal", "nocturnal/crepuscular", "crepuscular", "nocturnal", "NA", "4, 4, 4", "Hippopotimiade", "Hippopotamus_amphibius"))
rownames(whippomorpha) <- whippomorpha$tips
write.csv(whippomorpha, file = here("whippomorpha.csv"), row.names = FALSE)


# Section 1.5 High confidence cetacean data -------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit?usp=sharing'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#drop na values
cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$Diel_Pattern_2),]

cetaceans_full <- cetaceans_full[,1:10] 

#create a column with the max confidence level for that species (out of the confidence level for all sources)
#the confidence values are characters so convert to numerics and then take the maximum value
cetaceans_full$Confidence <- lapply(cetaceans_full$Confidence, function(x) gsub(",", "\\1 ", x))
cetaceans_full$Confidence <- lapply(cetaceans_full$Confidence,function(x) strsplit(x, " ")[[1]])
cetaceans_full$max_conf <- lapply(cetaceans_full$Confidence, as.numeric)
cetaceans_full$max_conf <- lapply(cetaceans_full$Confidence, max)

#filter for species that are in the mammal tree
#from 83 species to 75
cetaceans_full <- cetaceans_full %>% filter(!(Alt_name.notes %in% c("Not in mam tree", "Not in mam tree, Limited data", "Not in mam tree, Likely not a valid subspecies")))

#subset for the higher confidence levels (3-5)
#from 75 species to 70, may be worth dropping them
#for artiodactyla there's more, 50 level 1 only, 30 level 2
cetaceans_full <- cetaceans_full %>% filter(max_conf >= 3)

#save out

# Section 2: Formatting the binary artiodactyla diel dataframe ??? data (only contains di/noc)-----------------------

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


# Section 3: Formatting the artiodactyla diel dataframe Cox data from https://www.nature.com/articles/s41467-021-22023-4#data-availability --------

Cox_df <- read.csv(here("Cox_diel_activity_data.csv"))
Cox_df <- Cox_df %>% filter(Order == "Cetartiodactyla")

Cox_df <- Cox_df[,c("Binomial_iucn", "Activity_IM")]
Cox_df$tips <- str_replace(Cox_df$Binomial_iucn, " ", "_")
colnames(Cox_df) <- c("Species_name", "Diel_Pattern_2", "tips")
row.names(Cox_df) <- Cox_df$Species_name

#change any crepuscular entries to be cathemeral/crepuscular (lacks day-night preference but peaks activity at twilight)
#25 crepuscular only entries
for(i in 1:length(Cox_df$Diel_Pattern_2)){
  if(Cox_df[i, "Diel_Pattern_2"] %in% c("Crepuscular")){
    Cox_df[i, "Diel_Pattern_2"] <- "Cathemeral/Crepuscular"
  }
}

#create binary, max dinoc, max crep and 6-state columns

#diel pattern 1 is binary (di-noc only) -keeping this to preserve functioning of earlier models, not very biologically relevant
Cox_df$Diel_Pattern_1 <- Cox_df$Diel_Pattern_2
for(i in 1:nrow(Cox_df)){
  if(Cox_df[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    Cox_df[i, "Diel_Pattern_1"] <- "Diurnal"
  } else if(Cox_df[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    Cox_df[i, "Diel_Pattern_1"] <- "Nocturnal"
  } else if(Cox_df[i, "Diel_Pattern_2"] %in% c("Cathemeral", "Cathemeral/Crepuscular")){
    Cox_df[i, "Diel_Pattern_1"] <- NA
  }
}

# diel pattern 3 maximizes for crepuscularity
Cox_df$Diel_Pattern_3 <- Cox_df$Diel_Pattern_2
for(i in 1:nrow(Cox_df)){
  if(Cox_df[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    Cox_df[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if(Cox_df[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    Cox_df[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if (Cox_df[i, "Diel_Pattern_2"] == "Cathemeral/Crepuscular"){
    Cox_df[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if (Cox_df[i, "Diel_Pattern_2"] == "Cathemeral"){
    Cox_df[i, "Diel_Pattern_3"] <- "Crepuscular"}
}

# diel pattern 4 maximizes for diurnality and nocturnality (while keeping cathemerality as a trait state)
Cox_df$Diel_Pattern_4 <- Cox_df$Diel_Pattern_2
for(i in 1:nrow(Cox_df)){
  if(Cox_df[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    Cox_df[i, "Diel_Pattern_4"] <- "Diurnal"
  } else if(Cox_df[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    Cox_df[i, "Diel_Pattern_4"] <- "Nocturnal"
  } else if (Cox_df[i, "Diel_Pattern_2"] == "Cathemeral/Crepuscular"){
    Cox_df[i, "Diel_Pattern_4"] <- "Cathemeral"
  } 
}

Cox_df$Diel_Pattern_1 <- tolower(Cox_df$Diel_Pattern_1)
Cox_df$Diel_Pattern_2 <- tolower(Cox_df$Diel_Pattern_2)
Cox_df$Diel_Pattern_3 <- tolower(Cox_df$Diel_Pattern_3)
Cox_df$Diel_Pattern_4 <- tolower(Cox_df$Diel_Pattern_4)

#save out as a csv
write.csv(Cox_df, here("Cox_artiodactyla_without_cetaceans.csv"))

#now add in cetaceans
cetaceans <- read.csv("cetaceans_full.csv")
cetaceans <- cetaceans[, c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "Diel_Pattern_4", "tips")]
Cox_full <- rbind(Cox_df, cetaceans)
row.names(Cox_full) <- c(Cox_full$Species_name)

write.csv(Cox_full, here("Cox_artiodactyla_full.csv"))


# Section 4:Formatting the artiodactyla primary source data -------------------------------------

sheet <- 'https://docs.google.com/spreadsheets/d/1JGC7NZE_S36-IgUWpXBYyl2sgnBHb40DGnwPg2_F40M/edit?usp=sharing'
sleepy_artio <- read.csv(text=gsheet2text(sheet, format='csv'), stringsAsFactors=FALSE)
#save out a full version
write.csv(sleepy_artio, here("sleepy_artiodactyl_with_sources.csv"))

#subset the dataframe to exclude the columns that only contain source information
#this subset only takes the first confidence column, which is the highest confidence source
#when I made the spreadsheet I always put the highest confidence source first
sleepy_artio <- sleepy_artio[, 1:6]
sleepy_artio$tips <- sleepy_artio$Species_name
sleepy_artio$tips <- str_replace(sleepy_artio$tips, pattern = " ", replacement = "_")

#rename the row names to be the tip names so it's easier to subset by the tree tip labels later
row.names(sleepy_artio) <- sleepy_artio$tips

#drop na values
#sleepy_artio <- sleepy_artio[!is.na(sleepy_artio$Diel_Pattern_2),]
sleepy_artio$Confidence <- as.character(sleepy_artio$Confidence)
sleepy_artio <- sleepy_artio %>% replace_na(list(Diel_Pattern_2 = "unknown", Confidence = "unknown"))

#create binary, max dinoc, max crep and 6-state columns

#diel pattern 1 is binary (di-noc only) -keeping this to preserve functioning of earlier models, not very biologically relevant
sleepy_artio$Diel_Pattern_1 <- sleepy_artio$Diel_Pattern_2
for(i in 1:nrow(sleepy_artio)){
  if(sleepy_artio[i, "Diel_Pattern_2"] == "Diurnal/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_1"] <- "Diurnal"
  } else if(sleepy_artio[i, "Diel_Pattern_2"] == "Nocturnal/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_1"] <- "Nocturnal"
  } else if(sleepy_artio[i, "Diel_Pattern_2"] %in% c("Cathemeral", "Cathemeral/crepuscular")){
    sleepy_artio[i, "Diel_Pattern_1"] <- NA
  }
}

# diel pattern 3 maximizes for crepuscularity
sleepy_artio$Diel_Pattern_3 <- sleepy_artio$Diel_Pattern_2
for(i in 1:nrow(sleepy_artio)){
  if(sleepy_artio[i, "Diel_Pattern_2"] == "Diurnal/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if(sleepy_artio[i, "Diel_Pattern_2"] == "Nocturnal/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if (sleepy_artio[i, "Diel_Pattern_2"] == "Cathemeral/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_3"] <- "Crepuscular"
  } else if (sleepy_artio[i, "Diel_Pattern_2"] == "Cathemeral"){
    sleepy_artio[i, "Diel_Pattern_3"] <- "Crepuscular"}
}

# diel pattern 4 maximizes for diurnality and nocturnality (while keeping cathemerality as a trait state)
sleepy_artio$Diel_Pattern_4 <- sleepy_artio$Diel_Pattern_2
for(i in 1:nrow(sleepy_artio)){
  if(sleepy_artio[i, "Diel_Pattern_2"] == "Diurnal/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_4"] <- "Diurnal"
  } else if(sleepy_artio[i, "Diel_Pattern_2"] == "Nocturnal/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_4"] <- "Nocturnal"
  } else if (sleepy_artio[i, "Diel_Pattern_2"] == "Cathemeral/crepuscular"){
    sleepy_artio[i, "Diel_Pattern_4"] <- "Cathemeral"
  } 
}

sleepy_artio$Diel_Pattern_1 <- tolower(sleepy_artio$Diel_Pattern_1)
sleepy_artio$Diel_Pattern_2 <- tolower(sleepy_artio$Diel_Pattern_2)
sleepy_artio$Diel_Pattern_3 <- tolower(sleepy_artio$Diel_Pattern_3)
sleepy_artio$Diel_Pattern_4 <- tolower(sleepy_artio$Diel_Pattern_4)

sleepy_artio <- sleepy_artio %>% relocate(Diel_Pattern_1, .before = Diel_Pattern_2)
sleepy_artio <- sleepy_artio %>% relocate(Diel_Pattern_3, .after = Diel_Pattern_2)
sleepy_artio <- sleepy_artio %>% relocate(Diel_Pattern_4, .after = Diel_Pattern_3)

#drop na values
#sleepy_artio <- sleepy_artio[!is.na(sleepy_artio$Diel_Pattern_2),]

write.csv(sleepy_artio, here("sleepy_artiodactyla_minus_cetaceans.csv"))

cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- cetaceans_full[, c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "Diel_Pattern_4", "Confidence", "Parvorder", "tips")]
cetaceans_full <- cetaceans_full %>% relocate("Parvorder", .after = "Species_name")
colnames(cetaceans_full) <- c("Species_name", "Family", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "Diel_Pattern_4", "Confidence", "tips")
row.names(cetaceans_full) <- cetaceans_full$tips

sleepy_artio <- sleepy_artio[, c("Species_name", "Family", "Diel_Pattern_1", "Diel_Pattern_3", "Diel_Pattern_2", "Diel_Pattern_4", "Confidence", "tips")]
sleepy_artiodactyla_full <- rbind(sleepy_artio, cetaceans_full)

sleepy_artiodactyla_full$Order <- "idk"

for(i in 1:length(sleepy_artiodactyla_full$Species_name)){
  if(sleepy_artiodactyla_full[i, "Family"] %in% c("Camelidae")){
    sleepy_artiodactyla_full[i, "Order"] <- "Tylopoda"}
  else if(sleepy_artiodactyla_full[i, "Family"] %in% c("Suidae", "Tayassuidae")){
    sleepy_artiodactyla_full[i, "Order"] <- "Suina"}
  else if(sleepy_artiodactyla_full[i, "Family"] %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae")){
    sleepy_artiodactyla_full[i, "Order"] <- "Ruminantia"}
  else if(sleepy_artiodactyla_full[i, "Family"] %in% c("Mysticeti", "Odontoceti", "Hippopotamidae")){
    sleepy_artiodactyla_full[i, "Order"] <- "Whippomorpha"}
}

sleepy_artiodactyla_full <- sleepy_artiodactyla_full %>% relocate(Order, .after= Family)

#drop na values
#sleepy_artiodactyla_full <- sleepy_artiodactyla_full[!is.na(sleepy_artiodactyla_full$Diel_Pattern_2),]

write.csv(sleepy_artiodactyla_full, here("sleepy_artiodactyla_full.csv"))


# Section 5: Ruminant diel dataframe -------------------------------------------------
#read in artiodactyla minus cetaceans df
artio <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv"))

#filter for only ruminants, leaves 207 species
artio <- artio %>% filter(Family %in% c("Antilocapridae", "Bovidae", "Cervidae", "Giraffidae", "Moschidae", "Tragulidae"))

#203 of these are in the mammal tree
# artio <- artio[artio$tips %in% mam.tree$tip.label,]
# trpy <- keep.tip(mam.tree, tip = artio$tips)

write.csv(artio, here("ruminants_full.csv"))

#artio <- read.csv(here("ruminants_full.csv"))
artio <- read.csv(here("sleepy_artiodactyla_full.csv"))

#high confidence 
#create a column with the max confidence level for that species (out of the confidence level for all sources)
#the confidence values are characters so convert to numerics and then take the maximum value
artio$Confidence <- lapply(artio$Confidence, function(x) gsub(",", "\\1 ", x))
artio$Confidence <- lapply(artio$Confidence,function(x) strsplit(x, " ")[[1]])
artio$max_conf <- lapply(artio$Confidence, as.numeric)
artio$max_conf <- lapply(artio$Confidence, max)

artio$max_conf <- as.integer(artio$max_conf)

#of the 207 species, 145 have high confidence data
artio <- artio %>% filter(max_conf >= 3)

#gets rid of the list so we can save as a csv
artio <- apply(artio,2,as.character)

write.csv(artio, here("ruminants_full_high_conf.csv"))

# Section 6: Formatting the artiodactyla diel dataframe Maor data  --------
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

#now add in the alternative diel patterns 

#Conflicting diel patterns in Maor et al
#method A: if alt pattern doesn't match current diel pattern check the literature (currently 34 species, source = pantheria database)
#method B: if alt pattern matches current pattern replace with NA (alt columns will be deleted after)
#if alternative pattern is crepuscular, add it as a secondary diel pattern

#Method A: 35 species with unresolved diel patterns
#Maor et al gets activity pattern data from four sources, 57, 4, 33 and 29
#Source 57: 57.Jones, K.E. et al., 2009. PanTHERIA: A species-level database of life history, ecology and geography of extant and recently extinct mammals. Ecology, 90, p.2648.
#Source 4: Aulagnier, S. & Thevenot, M., 1986. Catalogue des Mammiferes Sauvages du Maroc, Rabat, Morocco: Institute Scientifique.
#Source 33: Hufnagl, E., 1972. Lybian Mammals, Harrow, England: The Oleander Press.
# Source 29: 29. Gray, G.G. & Simpson, C.D., 1980. Ammotragus lervia. Mammalian Species, 144, pp.1–7.
# Sources 4, 29, 33 are all for  Ammotragus lervia and from the 1970s/1980s. 


## Manually replace conflicting data
#Method B: Matches or near matches 
#Alcelaphus lichtensteinii
maor_full[4,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"

#Ammotragus lervia
maor_full[6,"alt_pattern_1"] <- NA

#Axis porcinus
maor_full[13,"alt_pattern_2"] <- NA

#Blastocerus dichotomus
maor_full[17,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"

#Bubalus bubalis
maor_full[23,"alt_pattern_2"] <- NA

#Bubalus mindorensis
maor_full[25,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"

#Cervus elaphus. Row 52 duplicated because there are four entries for it
#no new information in row 52 so drop it and renumber
maor_full <- maor_full[-c(52), ]
row.names(maor_full) <- 1:nrow(maor_full)

#Gazella subgutturosa
maor_full[63,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"

#Giraffa camelopardalis
maor_full[64,"alt_pattern_2"] <- NA

#Moschus moschiferus
maor_full[92,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"

#Philantomba maxwellii
maor_full[123,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"

#Rangifer tarandus
maor_full[131,"alt_pattern_2"] <- NA

#Sus barbatus
maor_full[145,"Diel_Pattern_2"] <- "Nocturnal/Crepuscular"

#Tragelaphus angasii
maor_full[156,"Diel_Pattern_2"] <- "Diurnal/Crepuscular"


#Resolve conflicting diel patterns with new information
#import dataframe of updated information on 34 conflicting species and match with their current row

missing_sps <- read.csv(here("Sleepy_artiodactyls_34_sps.csv"))
#Hexaprotodon liberiensis is not called Choeropsis liberiensis so no need to change (unlike in Cox et al)
row.names(missing_sps) <- missing_sps$Species_name


replacements <- match(missing_sps$Species_name, maor_full$Species)

for(i in replacements){
  maor_full[i, "Diel_Pattern_2"] <- missing_sps[maor_full[i, "Species"], "Diel_Pattern_2"]
}


# Rename and rearrange columns
#change any crepuscular entries to be cathemeral/crepuscular (lacks day-night preference but peaks activity at twilight)
for(i in 1:length(maor_full$Diel_Pattern_2)){
  if(maor_full[i, "Diel_Pattern_2"] %in% c("Crepuscular")){
    maor_full[i, "Diel_Pattern_2"] <- "Cathemeral/Crepuscular"
  }
}

#create column for activity pattern 1
# diel pattern 1 maximizes for only diurnal and nocturnal
maor_full$Diel_Pattern_1 <- maor_full$Diel_Pattern_2
for(i in 1:nrow(maor_full)){
  if(maor_full[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_1"] <- "Diurnal"
  } else if(maor_full[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_1"] <- "Nocturnal"
  } else if(maor_full[i, "Diel_Pattern_2"] %in% c("Crepuscular", "Cathemeral", "Cathemeral/Crepuscular")){
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
  } else if (maor_full[i, "Diel_Pattern_2"] == "Cathemeral/Crepuscular"){
    maor_full[i, "Diel_Pattern_3"] <- "Crepuscular"
  }
}

maor_full <- relocate(maor_full, "Diel_Pattern_2", .after = "Diel_Pattern_1")

# diel pattern 4 maximizes for diuranlity and nocturnality (while keeping cathemerality as a trait state)
maor_full$Diel_Pattern_4 <- maor_full$Diel_Pattern_2
for(i in 1:nrow(maor_full)){
  if(maor_full[i, "Diel_Pattern_2"] == "Diurnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_4"] <- "Diurnal"
  } else if(maor_full[i, "Diel_Pattern_2"] == "Nocturnal/Crepuscular"){
    maor_full[i, "Diel_Pattern_4"] <- "Nocturnal"
  } else if (maor_full[i, "Diel_Pattern_2"] == "Cathemeral/Crepuscular"){
    maor_full[i, "Diel_Pattern_4"] <- "Cathemeral"
  }
}


# Final formatting and saving out
#now we can format the artiodactyla data so it matches with the cetacean data and merge them
maor_full <- maor_full[, c("Species", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "Diel_Pattern_4")]
maor_full$Diel_Pattern_1 <- tolower(maor_full$Diel_Pattern_1)
maor_full$Diel_Pattern_2 <- tolower(maor_full$Diel_Pattern_2)
maor_full$Diel_Pattern_3 <- tolower(maor_full$Diel_Pattern_3)
maor_full$Diel_Pattern_4 <- tolower(maor_full$Diel_Pattern_4)

maor_full$tips <- str_replace(maor_full$Species, " ", "_")
colnames(maor_full) <- c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "Diel_Pattern_4", "tips")

cetaceans <- read.csv("cetaceans_full.csv")
cetaceans <- cetaceans[, c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "Diel_Pattern_4", "tips")]
maor_full <- rbind(maor_full, cetaceans)
row.names(maor_full) <- c(maor_full$Species_name)

write.csv(maor_full, here("Maor_artiodactyla_full.csv"))

# Section 7: Formatting the artiodactyla_without_cetaceans diel dataframe -----------------------

#easiest way to do this is to drop the last chunk of rows in maor_full (to keep all the formatting)
maor_full <- read.csv("Maor_artiodactyla_full.csv")
just_artio <- maor_full[1:151, ]
write.csv(just_artio, here("Maor_artiodactyla_without_cetaceans.csv"))


# Section 8: See how Cox and Maor data compare ---------------------------------------
maor_full <- read.csv(here("Maor_artiodactyla_without_cetaceans.csv"))
Cox_df <- read.csv(here("Cox_artiodactyla_without_cetaceans.csv"))

diel_merge <- merge(Cox_df,maor_full,by="Species_name")
diel_merge <- diel_merge[, c("Species_name", "Diel_Pattern_2.x", "Diel_Pattern_2.y")]
colnames(diel_merge) <- c("Species_name", "Cox_diel", "Maor_diel")
diel_merge$Cox_diel <- tolower(diel_merge$Cox_diel)
diel_merge$match <- "No"

for(i in 1:length(diel_merge$Species_name)){
  if(diel_merge[i, "Cox_diel"] == diel_merge[i, "Maor_diel"]){
    diel_merge[i, "match"] <- "Yes"
  } else if(diel_merge[i, "Cox_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal") & diel_merge[i, "Maor_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Cox_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")& diel_merge[i, "Maor_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Cox_diel"] %in% c("cathemeral", "cathemeral/crepuscular")& diel_merge[i, "Maor_diel"] %in% c("cathemeral", "cathemeral/crepuscular")){
    diel_merge[i, "match"] <- "Approximate"
  } else {
    diel_merge[i, "match"] <- "No"
  }
}

table(diel_merge$match)
#results: yes matches (73%), no matches (15%), approximate matches (13%) when including my cetacean data (will be the same for both)
#results: yes 86 (62%), no matches 28 (20%), approximate matches 24 (17%)
diel_merge <- diel_merge[order(diel_merge$match),]

write.csv(diel_merge, here("diel_merge.csv"))



# Section 9: Comparison of Cox and Maor data to primary literature -------------------
#see how well Cox and Maor datasets agree with the data I collected from the primary literature

#import missing species dataset
#use this for the preliminary 34 species dataset
#missing_sps <- read.csv(here("Sleepy_artiodactyls_34_sps.csv"))

#use this for the full 239 species dataset
sheet <- 'https://docs.google.com/spreadsheets/d/1JGC7NZE_S36-IgUWpXBYyl2sgnBHb40DGnwPg2_F40M/edit?usp=sharing'
missing_sps <- read.csv(text=gsheet2text(sheet, format='csv'), stringsAsFactors=FALSE)

#missing_sps <- read.csv(here("Sleepy_artiodactyls.csv"))

#missing_sps[14,"Species_name"] <- "Choeropsis liberiensis"
missing_sps[93,"Species_name"] <- "Choeropsis liberiensis"
missing_sps <- missing_sps[, c("Species_name", "Diel_Pattern_2", "Confidence")]
colnames(missing_sps) <- c("Species_name", "Amelia_diel", "Confidence")

#look at the distribution of data
table(missing_sps$Amelia_diel)
table(missing_sps$Amelia_diel, missing_sps$Confidence)

#format to match the Cox and Maor datasets (cathemeral/crepuscular species listed as crepuscular)
missing_sps$Amelia_diel <- str_replace_all(string = missing_sps$Amelia_diel, pattern = "Cathemeral/crepuscular", replacement  = "Crepuscular")
#remove NA species
missing_sps <- missing_sps[!(is.na(missing_sps$Amelia_diel)),]


#import unedited Cox dataset
Cox_df <- read.csv(here("Cox_diel_activity_data.csv"))
Cox_df <- Cox_df[, c("Binomial_iucn", "Activity_IM")]
colnames(Cox_df) <- c("Species_name", "Cox_diel")

Maor_df <- read_excel(here("Maor_diel_activity_data.xlsx"))
Maor_df <- Maor_df[17:nrow(maor_mam_data), 3:4]
colnames(Maor_df) <- c("Species_name", "Maor_diel")

#subset Cox and Maor df to only include species that are in the missing sps list
Cox_df <- Cox_df[Cox_df$Species_name %in% missing_sps$Species_name,]
Maor_df <- Maor_df[Maor_df$Species_name %in% missing_sps$Species_name,]

#merge and compare
diel_merge <- merge(Cox_df,missing_sps,by="Species_name")
diel_merge$Cox_diel <- tolower(diel_merge$Cox_diel)
diel_merge$Amelia_diel <- tolower(diel_merge$Amelia_diel)
diel_merge$match <- "No"

for(i in 1:length(diel_merge$Species_name)){
  if(diel_merge[i, "Cox_diel"] == diel_merge[i, "Amelia_diel"]){
    diel_merge[i, "match"] <- "Yes"
  } else if(diel_merge[i, "Cox_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal") & diel_merge[i, "Amelia_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Cox_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")& diel_merge[i, "Amelia_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Cox_diel"] %in% c("cathemeral", "cathemeral/crepuscular")& diel_merge[i, "Amelia_diel"] %in% c("cathemeral", "cathemeral/crepuscular")){
    diel_merge[i, "match"] <- "Approximate"
  } else {
    diel_merge[i, "match"] <- "No"
  }
}

table(diel_merge$match)
#results for 33 species: 3 yes matches, 11 no matches, 19 approximate matches
#results for all species: Approximate matches: 90, No: 76, Yes: 62
diel_merge <- diel_merge[order(diel_merge$match, diel_merge$Confidence),]
table(diel_merge$match, diel_merge$Confidence)
table(diel_merge$match, diel_merge$Amelia_diel)

#merge and compare
diel_merge <- merge(Maor_df,missing_sps,by="Species_name")
diel_merge$Maor_diel <- tolower(diel_merge$Maor_diel)
diel_merge$Amelia_diel <- tolower(diel_merge$Amelia_diel)
diel_merge$match <- "No"

for(i in 1:length(diel_merge$Species_name)){
  if(diel_merge[i, "Maor_diel"] == diel_merge[i, "Amelia_diel"]){
    diel_merge[i, "match"] <- "Yes"
  } else if(diel_merge[i, "Maor_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal") & diel_merge[i, "Amelia_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Maor_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")& diel_merge[i, "Amelia_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Maor_diel"] %in% c("cathemeral", "cathemeral/crepuscular")& diel_merge[i, "Amelia_diel"] %in% c("cathemeral", "cathemeral/crepuscular")){
    diel_merge[i, "match"] <- "Approximate"
  } else {
    diel_merge[i, "match"] <- "No"
  }
}

table(diel_merge$match)
#results for 50 species (repeat rows from Maor et al): 3 yes matches, 20 no matches, 27 approximate matches
#results for all species: Approximate matches: 52, No: 39, Yes: 22
diel_merge <- diel_merge[order(diel_merge$match, diel_merge$Confidence),]
table(diel_merge$match, diel_merge$Confidence)
table(diel_merge$match, diel_merge$Amelia_diel)


# Section 10: Load in and examine the mam tree, check time calibration -----------------------------------------

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
cetaceans_full <- cetaceans_full[,-(1)]
cetaceans_full <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n_test <- keep.tip(mam.tree, tip = cetaceans_full$tips)
#check if time calibrated
test <- ggtree(trpy_n_test, layout = "circular", size = 0.5) + geom_tiplab(size = 1.5)
#test$data #29.36 million years ago to the root. Seems correct!

#look at the tree for whippomorpha (hippos + cetacea) and see if time calibrated
whippomorpha <- read.csv(here("whippomorpha.csv"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
whippomorpha <- whippomorpha[whippomorpha$tips %in% mam.tree$tip.label,]
trpy <- keep.tip(mam.tree, tip = whippomorpha$tips)
#to add the timescale: geom_timescale() adds a scale for one region
#theme_tree2() adds a scale to the entire tree
whippo <- ggtree(trpy, layout = "rectangular", size = 0.5) + theme_tree2()
whippo <- revts(whippo)
#add label for epochs
whippo <- whippo + coord_geo(dat = epochs, xlim = c(-60, 0), ylim = c(-2, Ntip(trpy)), neg = TRUE, size =3, abbrv = FALSE)
whippo <- whippo + geom_tiplab(size = 2) 
whippo
#whippo$data shows the branch lengths (53.7 my to root). Also seems correct! 

#timescale for artiodactyla
artio <- read.csv(here("sleepy_artiodactyla_full.csv"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
artio <- artio[artio$tips %in% mam.tree$tip.label,]
trpy <- keep.tip(mam.tree, tip = artio$tips)
#to add the timescale: geom_timescale() adds a scale for one region
#theme_tree2() adds a scale to the entire tree
artio <- ggtree(trpy, layout = "rectangular", size = 0.5) + theme_tree2()
artio <- revts(artio)
#add label for epochs
artio <- artio + coord_geo(dat = epochs, xlim = c(-80, 0), ylim = c(-2, Ntip(trpy)), neg = TRUE, size =2, abbrv = FALSE)
artio <- artio + geom_tiplab(size = 2) 
artio



