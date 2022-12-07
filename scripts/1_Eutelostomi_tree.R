library(rotl)
library(stringr)
library(ape)
library(geiger)
library(xlsx)

## OK, first part has to be run locally (because it needs internet access)

setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")

######################################################## 
###################### FIRST PART ###################### 
######################################################## 

# setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")
# 
# # Load in data files with names and activity patterns
# # Modify, subset, and combine
# 
# fish_data <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
# fish_data <- fish_data[,c("unique_name", "diel", "genus", "family", "order")]
# fish_data$group <- "fish"
# 
# mam_data <- read.xlsx("Cox_mammal_data/Supplementary Data 2.xlsx", 1)
# mam_data <- mam_data[,c("Binomial_iucn", "Activity_DD", "Genus", "Family", "Order")]
# colnames(mam_data) <- c("unique_name", "diel", "genus", "family", "order")
# mam_data$group <- "mammals"
# mam_data <- mam_data[mam_data$diel != "NA",]
# 
# tet_data <- read.csv("tetrapod_data/suppfile10appendix1914.csv")
# tet_data <- tet_data[tet_data$Species %in% tet_data$Species[tet_data$Class != "Mammalia"],]
# tet_data$group <- tet_data$Class
# tet_data$Genus <- gsub("\\_.*","",tet_data$Species)
# tet_data <- tet_data[,c("Species", "State", "Genus", "Family", "Order", "group")]
# colnames(tet_data) <- c("unique_name", "diel", "genus", "family", "order", "group")
# 
# all_data <- Reduce(rbind, list(fish_data, mam_data, tet_data))
# all_data$unique_name <- str_replace(all_data$unique_name, "_", " ")
# 
# all_data$diel <- str_replace(all_data$diel, "ARR", "unclear")
# all_data$diel <- str_replace(all_data$diel, "Cathemeral", "unclear")
# all_data$diel <- str_replace(all_data$diel, "NOC", "nocturnal")
# all_data$diel <- str_replace(all_data$diel, "DIU", "diurnal")
# all_data$diel <- str_replace(all_data$diel, "CRE", "crepuscular")
# all_data$diel <- tolower(all_data$diel)
# all_data <- all_data[!(is.na(all_data$diel)),]
# 
# # Fetch species from tree of life using rotl package
# 
# resolved_names <- tnrs_match_names(all_data$unique_name, context_name = "Vertebrates", do_approximate_matching = FALSE)
# 
# # Remove any that don't have exact matches or are synonyms (this is really just for finding ancestral state, so missing a few species won't matter)
# resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
# resolved_names <- resolved_names[resolved_names$is_synonym == FALSE,]
# resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]
# 
# # Remove excess information, clean up, and add tip label ids that will match the tree
# resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
# resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
# resolved_names <- resolved_names[!duplicated(resolved_names$tips),]
# 
# # Add data on activity
# 
# resolved_names$diel <- all_data$diel[match(resolved_names$search_string, tolower(all_data$unique_name))]
# 
# resolved_names$genus <- all_data$genus[match(resolved_names$unique_name, all_data$unique_name)]
# resolved_names$family <- all_data$family[match(resolved_names$unique_name, all_data$unique_name)]
# resolved_names$order <- all_data$order[match(resolved_names$unique_name, all_data$unique_name)]
# 
# ## Fetch the tree
# 
# tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") # I need to use the id option here, and then use that to map the tip labels from resolved_names (that way I don't run into the issue with the difference in formatting between the two tools)
# 
# # Time calibrate it using geiger and timetree.org
# # First resolve polytomies ~randomly using multi2dr
# 
# tr <- multi2di(tr)
# 
# # Save and re-load files
# saveRDS(tr, file = "tr_tree_AllGroups.rds")
# saveRDS(resolved_names, file = "resolved_names_AllGroups.rds")


######################################################## 
###################### SECOND PART #####################
######################################################## 

## This part is done on the cluster (calibrating with the timetrees)
## Run the 1st part, commit and then run this part

tr <- readRDS(file = "tr_tree_AllGroups.rds")
resolved_names <- readRDS(file = "resolved_names_AllGroups.rds")

# Make the reference file
# Ensure that the rownames and tip.labels in the target match the species names in the reference

resolved_names$ott_id <- paste("ott", resolved_names$ott_id, sep = "")

reference.df <- resolved_names[resolved_names$ott_id %in% tr$tip.label,c("order", "family", "genus", "unique_name", "tips", "ott_id")] 
colnames(reference.df) <- c("order", "family", "genus", "unique_name", "tips_species", "tips")
rownames(reference.df) <- reference.df$tips

# # two tips can't be found in the resolved_names df, but I cannot figure out why
# > tr$tip.label[!(tr$tip.label %in% resolved_names$ott_id)]
# [1] "mrcaott320143ott351725" "mrcaott106188ott185786"

# some are duplicated, or have missing data, remove them
reference.df <- reference.df[!duplicated(reference.df$unique_name),]
reference.df <- reference.df[!is.na(reference.df$unique_name),]

# Load the timetree tree (genus level data works, but not species)
# Have download timetree data for species, genus, family, and order
# Genus level data has the most calibration points
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

timetree_order <- ape::read.tree("timetree_data/euteleostomi_order.nwk")
timetree_family <- ape::read.tree("timetree_data/euteleostomi_family.nwk")
timetree_genus <- ape::read.tree("timetree_data/euteleostomi_genus.nwk")

# Use geiger to congruify the tree, works with treePL
# This seems to work up to genus, but not species (by replacing tip.labels with the same names)
setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")
geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")
geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")
geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

tr.calibrated <- geiger.genus$phy
tr.calibrated$tip.label <- resolved_names$tips[match(tr.calibrated$tip.label, resolved_names$ott_id)]

## Save out files

saveRDS(tr.calibrated, file = "tr_tree_calibrated_AllGroups.rds")

######################################################## 
###################### THIRD PART ######################
######################################################## 

## load back in to do by species manually locally

# tr.calibrated <- readRDS("calibrated_phylo.rds")

# Below works if you modify the heights.phylo function
# trace(geiger:::heights.phylo, edit = TRUE)
# depth = max(xx[!(is.na(xx))])
# Also have to do it manually - which is ugh!
# timetree_species <- ape::read.tree("timetree_data/actinopterygii_species.nwk")
# timetree_species <- multi2di(timetree_species)
# setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

# geiger.species <- congruify.phylo(reference = timetree_species, target = tr.calibrated, taxonomy = reference.df, tol = 0, scale = "treePL")
# tr.calibrated <- geiger.species$phy

# tr.calibrated$tip.label <- resolved_names$tips[match(tr.calibrated$tip.label, resolved_names$ott_id)]

# saveRDS(tr.calibrated, file = "calibrated_phylo.rds")


# ## Make trait.data file
# 
# trait.data <- data.frame(ott_id = tr.calibrated$tip.label, species = resolved_names$tips[match(tr.calibrated$tip.label, paste("ott", resolved_names$ott_id, sep = ""))], diel = resolved_names$diel[match(tr.calibrated$tip.label, paste("ott", resolved_names$ott_id, sep = ""))]) # OK, some species tip labels are more complicated and cause issues here
# 
# # Create vectors including crepuscular/unclear, or not
# trait.data$diel1 <- ifelse(trait.data$diel %in% c("diurnal", "crepuscular/diurnal"), "diurnal", ifelse(trait.data$diel %in% c("nocturnal", "crepuscular/nocturnal"), "nocturnal", "unclear/crepuscular"))
# levels(trait.data$diel1) <- c("diurnal", "nocturnal", "unclear/crepuscular")
# trait.data$diel2 <- ifelse(trait.data$diel %in% c("diurnal"), "diurnal", ifelse(trait.data$diel %in% c("nocturnal"), "nocturnal", ifelse(trait.data$diel %in% c("crepuscular", "crepuscular/diurnal", "crepuscular/nocturnal"), "crepuscular", "unclear")))
# levels(trait.data$diel2) <- c("diurnal", "nocturnal", "crepuscular", "unclear")
# 
# trait.data <- trait.data[!(is.na(trait.data$diel)),]
# rownames(trait.data) <- trait.data$species
# 
# trait.data$tips <- resolved_names$tips[match(trait.data$species, resolved_names$tips)]
# trait.data$order <- resolved_names$order[match(trait.data$species, resolved_names$tips)]
# 
# saveRDS(trait.data, file = "trait_data_AllGroups.rds")
