library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)
library(phytools)
library(rfishbase)
library(xlsx)
library(geiger)
library(here)

setwd(here())

### LOAD IN OTHER DATA ### 

# Fetch the taxonomic levels from fishbase for all of the species and make it into a dataframe
# available_releases()
# [1] "23.01" "23.05" "21.06" "19.04"

fishbase_df <- load_taxa(collect = T, version = "21.06")
fishbase_df <- as.data.frame(fishbase_df)

### LOAD IN GOOGLE SHEET DATA FROM LOCAL ### 

sleepy_fish <- read.csv(file = here("sleepy_fish_database_local_biorxiv_20230519.csv"))


###############################################################################################################################################
### FETCH DATA FROM TREE OF LIFE ### 
################################################################################################################################################

## Remove fish without diel data, or real species names (necessary?)
sleepy_fish <- sleepy_fish[sleepy_fish$Diel_Pattern != "",]
sleepy_fish <- sleepy_fish[!(grepl("sp.|spp.|unidentified", sleepy_fish$Species)),]

## Fetch species from tree of life using rotl package
resolved_names <- tnrs_match_names(sleepy_fish$Species_name, context_name = "Vertebrates", do_approximate_matching = FALSE)

# ## Can use this to check those that don't have an exact match (in case there are new fish added)
# resolved_names_2 <- tnrs_match_names(resolved_names$search_string[is.na(resolved_names$unique_name)], context_name = "Vertebrates", do_approximate_matching = TRUE)

## Remove any that don't have exact matches
resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]

## Remove ambiguous matches (approximate matches, things with multiple matches)
resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]

## Print the number of matches
print(paste("rotl found matches for", nrow(resolved_names), "out of", nrow(sleepy_fish), "from the Sleepy fish database", sep = " "))

# Remove excess information, clean up, and add tip label ids that will match the tree
# resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
resolved_names <- resolved_names[!duplicated(resolved_names$tips),]

resolved_names$genus <- fishbase_df$Genus[match(resolved_names$unique_name, fishbase_df$Species)]
resolved_names$family <- fishbase_df$Family[match(resolved_names$unique_name, fishbase_df$Species)]
resolved_names$order <- fishbase_df$Order[match(resolved_names$unique_name, fishbase_df$Species)]

## Add data on traits (from sleepy_fish, and/or from fishbase)

resolved_names$diel <- sleepy_fish$Diel_Pattern[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]
resolved_names$diel <- factor(resolved_names$diel, levels = c("diurnal", "nocturnal", "unclear", "crepuscular", "crepuscular/diurnal", "crepuscular/nocturnal", "crepuscular/unclear", "unclear/diurnal", "unclear/nocturnal"))
resolved_names$diel2 <- ifelse(resolved_names$diel == "diurnal", "diurnal", ifelse(resolved_names$diel == "nocturnal", "nocturnal", ifelse(resolved_names$diel == "crepuscular", "crepuscular", ifelse(resolved_names$diel == "crepuscular/diurnal", "crepuscular", ifelse(resolved_names$diel == "crepuscular/nocturnal", "crepuscular", ifelse(resolved_names$diel == "crepuscular/unclear", "crepuscular", ifelse(resolved_names$diel == "unclear/diurnal", "unclear", ifelse(resolved_names$diel == "unclear/nocturnal", "nocturnal", "unknown"))))))))
resolved_names$diel_confidence <- sleepy_fish$Confidence[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]
resolved_names$genome <- sleepy_fish$Genome[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]


print(paste("Sleepy fish database covers ", round((length(unique(resolved_names$tips))/length(unique(fishbase_df$Species)))*100), "% of Species, ", round((length(unique(resolved_names$genus))/length(unique(fishbase_df$Genus)))*100), "% of Genuses, ", round((length(unique(resolved_names$family))/length(unique(fishbase_df$Family)))*100), "% of Families, and ", round((length(unique(resolved_names$order))/length(unique(fishbase_df$Order)))*100), "% of Orders", sep = ""))

################################################################################################################################################
### Determine searches I've already done ### 
################################################################################################################################################
## quick wrapper to make things simplier
'%out%' <- Negate('%in%')

url <- 'https://docs.google.com/spreadsheets/d/1ZFj589su_77lAKKQe8GHQftpFPD4808yrWNxSUrhwO4/edit?gid=0#gid=0'
random1500 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
url <- 'https://docs.google.com/spreadsheets/d/16p1vDdWxARV07aHkxo6Ub1RTNNviHKVoa5qYMOJ61hs/edit?gid=0#gid=0'
eco <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
eco$species <- str_split_fixed(eco$X, "%22",4)[,c(2)]
url <- 'https://docs.google.com/spreadsheets/d/1ACisO47K3LoMXJP6xCOcaKQvKWyUxcsNEDKB7sWlxKE/edit?gid=1842590696#gid=1842590696'
fam <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

not_searched <- fishbase_df[fishbase_df$Species %out% random1500$species,]
not_searched <- not_searched[not_searched$Species %out% eco$species,]
not_searched <- not_searched[not_searched$Genus %out% fam$genus,]

not_searched_not_found <- not_searched[not_searched$Species %out% resolved_names$unique_name,]

###############################################################################################################################################
### Determine missing clades ### 
################################################################################################################################################

species <- unique(not_searched$Species)[unique(not_searched$Species) %out% unique(resolved_names$unique_name)]
genus <- unique(not_searched$Genus)[unique(not_searched$Genus) %out% unique(resolved_names$genus)]
family <- unique(not_searched$Family)[unique(not_searched$Family) %out% unique(resolved_names$family)]
order <- unique(not_searched$Order)[unique(not_searched$Order) %out% unique(resolved_names$order)]

scholar_genus <- paste("https://scholar.google.com/scholar?start=0&q=%22", genus, " ", "%22+AND+(%22circadian%22+OR+%22diel%22+OR+%22diurnal%22+OR+%22nocturnal%22+OR+%22crepuscular%22)&hl=en&as_sdt=0,5", sep = "")

library(clipr)
write_clip(scholar_genus)

### This was to generate the random1500 above (already searched as of 19.07.2024)
# ## I think I should just take a reproducible subset of the missing species (500-1000), and generate searches for these
# 
# set.seed(12345678)
# sub.species <- sample(species, 1500)
# 
# 
# new_df <- data.frame(species = sub.species, genus = not_searched$Genus[match(sub.species, not_searched$Species)], family = not_searched$Family[match(sub.species, not_searched$Species)], order = not_searched$Order[match(sub.species, not_searched$Species)], data = sub.species %in% sleepy_fish$Species.name)
# 
# new_df$scholar <- paste("https://scholar.google.com/scholar?start=0&q=%22", new_df$species, "%22+AND+(%22circadian%22+OR+%22diel%22+OR+%22diurnal%22+OR+%22nocturnal%22+OR+%22crepuscular%22)&hl=en&as_sdt=0,5", sep = "")
# 
# library(clipr)
# write_clip(new_df)