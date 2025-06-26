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

#load in the dataframe with the tabulated activity patterns (objective calls based on source concordance)
#81 species with data, should be 83
cetaceans_tabulated_full <- read.csv(here("cetacean_tabulated_full.csv"))

#load in full primary source dataframe, 98 species
url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit?usp=sharing'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#add in the tabulated diel patterns
cetaceans_full <- merge(cetaceans_full, cetaceans_tabulated_full, by = "Species_name", all.x = TRUE)
#save out full version with sources
write.csv(cetaceans_full, here("cetaceans_full_with_sources.csv"))

#remove unnecessary columns
cetaceans_full <- cetaceans_full[c("Species_name", "Confidence", "Parvorder", "tabulated_diel")]
#add a column for tips, formatted as the species names appear in the phylogenetic tree
cetaceans_full$tips <- cetaceans_full$Species_name
cetaceans_full$tips <- str_replace(cetaceans_full$tips, pattern = " ", replacement = "_")
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

#save out a local copy in case google goes bankrupt
write.csv(cetaceans_full, file = here("cetaceans_full.csv"), row.names = FALSE)

#save out a version with hippos
whippomorpha <- read.csv(here("cetaceans_full.csv"))
whippomorpha <- rbind(whippomorpha, c("Hexaprotodon liberiensis", 3.0, "Hippopotamidae", "nocturnal/crepuscular", "Hexaprotodon_liberiensis", "whippomorpha", "crepusuclar", "nocturnal"))
whippomorpha <- rbind(whippomorpha, c("Hippopotamus amphibius", 4.0, "Hippopotamidae", "nocturnal/crepuscular", "Hippopotamus_amphibius", "whippomorpha", "crepusuclar", "nocturnal"))
rownames(whippomorpha) <- whippomorpha$tips
write.csv(whippomorpha, file = here("whippomorpha.csv"), row.names = FALSE)


# Section 2: Formatting the artiodactyl diel dataframe  -----------------


