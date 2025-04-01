# Section 0: Packages -----------------------------------------------------
library(corHMM)   
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(tidyr)
library(lubridate)
library(ggplot2)
library(geiger)
library(ape)

setwd(here())

source("scripts/Amelia_functions.R")

# Section 1: Formatting the arthropod diel dataframe -----------------------

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1eWc1JeGCXx9Gk8DtCQRGZaOPcRmT01pzlZfAcbccV-o/edit?gid=0#gid=0'
#use the url to load in the google sheet as a csv, specify anything you want encoded as na in na.strings
bugs_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE, na.strings = c("Couldn't find", "Generalist"))
#save out full version with sources
#write.csv(bugs_full, here("bugs_full_with_sources.csv"))

#drop everything except for species name and diel pattern

bugs_full <- bugs_full[, c("Species_name", "Diel_Pattern", "Order")]

#in the open tree of life, species names are written in the format Genus_species
#create a column called tips formatting the species names this way
#saves you from remaking it for each of the trait.data dataframes
bugs_full$tips <- bugs_full$Species_name
# This line below was causing the problem, it now can successfully match species to open tree of life
bugs_full$tips <- str_replace(bugs_full$tips, pattern = " ", replacement = "_") 

#there are duplicates of the same name so we have to resolve these
#Litocala sexsignata is in the dataframe twice so delete one 
#this will remove any duplicates
dupe_list <- bugs_full[duplicated(bugs_full$Species_name),]  #just added this line to check the dupes that got removed, so I can go back in my google sheet and update it
bugs_full <- bugs_full[!duplicated(bugs_full$Species_name),]

#rename the row names to be the tip names so it's easier to subset by the tree tip labels later
row.names(bugs_full) <- bugs_full$tips

#we can use three trait states for now: diurnal, nocturnal and crepuscular

#set all diel pattern entries to lower case, helps keep consistency
bugs_full$Diel_Pattern <- tolower(bugs_full$Diel_Pattern)

#drop na values  #david note to self: 
bugs_full <- bugs_full[!is.na(bugs_full$Diel_Pattern),]
bugs_full <- bugs_full[bugs_full$Diel_Pattern %in% c("nocturnal", "diurnal", "crepuscular"),]


#reclassify all species that aren't diurnal, nocturnal or crepuscular
#if later we want to use crepuscular dusk and crepuscular dawn/dusk as separate states we can revisit

#this is to get rid of the brackets, no longer needed
# bugs_full$Diel_Pattern <- str_replace(bugs_full$Diel_Pattern, pattern = "[(]", replacement = "")
# bugs_full$Diel_Pattern <- str_replace(bugs_full$Diel_Pattern, pattern = "[)]", replacement = "")
# 
# #use string replace to remove "dusk and dawn" and "dusk" to get just crepuscular
# bugs_full$Diel_Pattern <- str_replace(bugs_full$Diel_Pattern, pattern = "dawn and dusk", replacement = "")
# bugs_full$Diel_Pattern <- str_replace(bugs_full$Diel_Pattern, pattern = "dusk, males unspecified", replacement = "")
# bugs_full$Diel_Pattern <- str_replace(bugs_full$Diel_Pattern, pattern = "dusk", replacement = "")

#keep only rows that are crepuscular, diurnal or nocturnal  #david note to self: (removed 5)
diel_typos <- bugs_full %>% filter(!Diel_Pattern %in% c("diurnal", "nocturnal", "crepuscular")) #I added this line to catch the potential typos in diel pattern
bugs_full <- bugs_full %>% filter(Diel_Pattern %in% c("diurnal", "nocturnal", "crepuscular"))

# Section 2: Plotting on the open tree of life -----------------------
#the open tree of life is a synthetic tree that can be accessed online here https://tree.opentreeoflife.org/opentree/argus/opentree15.1@ott93302

#trait data has all the behavioural data for the 953 species with data 
trait.data <- bugs_full

#find the species names that are in the open tree of life, will return a warning with any missing species
#use tnrs_contexts() to see what possible categories we can search within, for you probably Arthropods
resolved_names <- tnrs_match_names(names = trait.data$Species_name, context_name = "Animals", do_approximate_matching = TRUE)

## Remove any that don't have exact matches or are synonyms (this is really just for finding ancestral state, so missing a few species won't matter)
#Remove species without open tree of life entries. Reduces 953 sps to 887
resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
#Removes species that are synonyms for other species. 
#no approximate matches, so this removes no species currently
#with approximate matching, this removes 53 species
approximate_matches <- resolved_names[resolved_names$approximate_match == TRUE,] #David should go through these and add alternative names to his database
resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]

#we've removed 53 species, so now we have 834 arthropods in our resolved_names dataframe to work with

# Remove excess information, clean up, and add tip label ids that will match the tree
resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
resolved_names <- resolved_names[!duplicated(resolved_names$search_string),]

# Add data on activity patterns using diel pattern column (crep, di, noc)
# match() finds the position of its first argument in its second
resolved_names$diel <- trait.data$Diel_Pattern[match(resolved_names$search_string, tolower(trait.data$Species_name))]
#check that the activity patterns look correct
table(resolved_names$diel)
#weeeeee
## Fetch the tree
#subset the open tree of life to only include the species in resolved_names (find them by their ott_id)
#David has other flags: "incertae_sedis_inherited", "sibling_higher", "incertae_sedis", "infraspecific"           
tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") 
# I need to use the id option here, and then use that to map the tip labels from resolved_names (that way I don't run into the issue with the difference in formatting between the two tools)

## Print the number of matches
print(paste("rotl found matches for", nrow(resolved_names), "out of", nrow(bugs_full), "from the arthropod database", sep = " "))
# currently 834 out of 953

# Resolve polytomies ~randomly using multi2dr
tr <- multi2di(tr)

# Make the reference file
# Ensure that the rownames and tip.labels in the target match the species names in the reference

#currently the tips are identified by their ott_id, can change to their species name (tip column)
tr$tip.label <- resolved_names$tips[match(tr$tip.label, paste("ott", resolved_names$ott_id, sep = ""))]

#this is the tree of all the cetacean species we are working with in the open tree of life
bug_tree <- ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
bug_tree

png(here("David_tree.png"), units = 'cm', height = 10, width = 10, res = 900)   #OG setting: 30x30, res = 900
bug_tree
dev.off()

resolved_names <- resolved_names[resolved_names$tips %in% tr$tip.label,]
tr <- keep.tip(tr, tip = resolved_names$tips)

#add diel pattern data
##add the diel pattern data as a geom tile
custom.colours <- c("aquamarine1", "khaki1", "goldenrod")
diel.plot <- ggtree(tr, layout = "circular") %<+% resolved_names[, c("tips", "diel")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = diel), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot.labelled <- diel.plot + geom_tiplab(color = "black", size = 1.5)
diel.plot.labelled

png(here("David_tree_diel_pattern.png"), units = 'cm', height = 20, width = 20, res = 900)   #OG setting: 30x30, res = 900
diel.plot.labelled
dev.off()

# # Section 3: Plotting the diel patterns labelled by major orders ---------------

#we want to create a tree labelled by order (column 3)
resolved_names$order <- trait.data$Order[match(resolved_names$tips, trait.data$tips)]
resolved_names <- resolved_names[!is.na(resolved_names$order), ]
table(resolved_names$order)

#this function will not label any order with one species
#findMRCANode2(phylo = tr, trait.data = resolved_names, taxonomic_level_col = 7, taxonomic_level_name = "Araneae")

resolved_names <- resolved_names %>% group_by(order) %>% filter(n()>1)

all_nodes <- lapply(unique(resolved_names$order), function(x) findMRCANode2(phylo = tr, trait.data = resolved_names, taxonomic_level_col = 7, taxonomic_level_name = x))
nodes_df <- do.call(rbind, all_nodes)

node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node), colour = "blue", hjust=-.2, size = 3) + geom_tiplab(size = 3, hjust = -0.1)
node_labels 

#can now easily label all clades within that taxonomic level on the tree using the nodes_df
#nodes_df$colour <- c("navy", "slateblue", "mediumpurple", "dodgerblue", "darkorchid1", "royalblue", "lightslateblue", "purple3", "steelblue", "red", "green")
order_tree <- diel.plot + geom_cladelab(node = nodes_df$node_number, label = nodes_df$clade_name, offset = 1, offset.text = 3, hjust = 0.3, align = TRUE, fontsize = 1.8, barsize = 1.5, barcolour = "grey45") #nodes_df$colour
order_tree

png(here("David_tree_bug_orders_diel_pattern.png"), units = 'cm', height = 25, width = 25, res = 900)   #OG setting: 30x30, res = 900
order_tree
dev.off()


# Section 4: Diel plots for each order ------------------------------------

#filter for families that have 15 or more species (can change this number)
trait.data.filtered <- trait.data %>% group_by_at(taxonomic_level_col) %>% filter(n()>15)

tree_list <- list()

for(i in trait.data.filtered$Order){
  to_drop <- trait.data.filtered[!(trait.data.filtered$Order == i), ]
  #use functions drop tip to remove these species
  trim_tr <- drop.tip(tr, to_drop$tips)
  tree_list[[i]] <- trim_tr
}

order_subtree <- ggtree(tree_list$Coleoptera, layout = "circular") %<+% trait.data[, c("tips", "Diel_Pattern")]
order_subtree <- order_subtree + geom_tile(data = order_subtree$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours)
order_subtree <- order_subtree + geom_tiplab(color = "black", size = 1.5)
order_subtree


# Section 5: Time calibrating the tree ---------------
# Time calibrate it using geiger and timetree.org 

#Retrieve the taxonomic info from GBIF
#load in the large dataframe of taxonomic info (around 200MB), 1,048,575 rows
GBIF_masterlist <- read.csv(here("scripts/arthropoda_master_list_gbif_20230828.csv"))

#remove rows with missing or blank names, leaves 911,981 rows
GBIF_masterlist <- GBIF_masterlist[!is.na(GBIF_masterlist$canonicalName),]
GBIF_masterlist <- filter(GBIF_masterlist, !(canonicalName == ""))

#filter to only entries with taxonRank species (and subspecies?) Ask David if he has subspecies. Leaves 728,357 rows
GBIF_masterlist <- filter(GBIF_masterlist, taxonRank == "species")

arthropod_df  <- filter(GBIF_masterlist, canonicalName %in% resolved_names$unique_name)

missing_bugs  <- resolved_names %>% filter(!(unique_name %in% GBIF_masterlist$canonicalName))

print(paste("gbif found matches for", nrow(arthropod_df), "out of", nrow(bugs_full), "from the arthropod database", sep = " "))
# currently 305 out of 940
# 397 out of 945

#save out the subset dataframe so we can push to the github (yay)
#write.csv(bugs, here("arthropoda_master_list_subset.csv"))
# arthropod_df <- read.csv(here("arthropoda_master_list_subset.csv"))
# arthropod_df <- GBIF_masterlist

resolved_names$genus <- arthropod_df$genus[match(resolved_names$unique_name, arthropod_df$canonicalName)]
resolved_names$family <- arthropod_df$family[match(resolved_names$unique_name, arthropod_df$canonicalName)]
resolved_names$order <- arthropod_df$order[match(resolved_names$unique_name, arthropod_df$canonicalName)]

#remove any species without genus, family, order information
resolved_names <- resolved_names[!(is.na(resolved_names$genus)), ]

reference.df <- resolved_names[resolved_names$ott_id %in% tr$tip.label,c("order", "family", "genus", "unique_name", "tips", "ott_id")]
colnames(reference.df) <- c("order", "family", "genus", "unique_name", "tips_species", "tips")

#Biston_betularia is in there twice so remove any duplicates
reference.df <- reference.df[!duplicated(reference.df$tips),]

rownames(reference.df) <- reference.df$tips

# Load the timetree tree (genus level data works, but not species)
# Have download timetree data for species, genus, family, and order
# Genus level data has the most calibration points

#replace with arthropod tree
timetree_order <- ape::read.tree(here("arthropods_order.nwk"))
timetree_family <- ape::read.tree(here("arthropods_family.nwk"))
timetree_genus <- ape::read.tree(here("arthropods_genus.nwk"))

# Use geiger to congruify the tree, works with treePL
setwd(here())
geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd(here())
geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd(here())
geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

tr.calibrated <- geiger.genus$phy

## Save out files

saveRDS(tr.calibrated, file = "calibrated_phylo_2024-09-07.rds")


