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

# ### LOAD IN OTHER DATA ### 
# 
# # Fetch the taxonomic levels from fishbase for all of the species and make it into a dataframe
# # available_releases()
# # [1] "23.01" "23.05" "21.06" "19.04"
# 
# fishbase_df <- load_taxa(collect = T, version = "21.06")
# fishbase_df <- as.data.frame(fishbase_df)
# 
# ### LOAD IN GOOGLE SHEET DATA ### 
# 
# # Load in the data from google sheets
# # and do some clean up on the data
# 
# url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit#gid=0'
# sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
# sleepy_fish$Diel_Pattern <- tolower(sleepy_fish$Diel_Pattern)
# 
# write.csv(sleepy_fish, file = here("sleepy_fish_database_local_2024-08-23.csv"))
# 
# # ################################################################################################################################################
# # ### OUTPUT DATA FOR ZUZANNA ### 
# # ################################################################################################################################################
# # 
# # species_for_nocturnal_database <- read.xlsx("~/Downloads/Species_for_nocturnal_database.xlsx", 1, header = F)
# # 
# # sleepy_fish_sub <- sleepy_fish[sleepy_fish$Species_name %in% species_for_nocturnal_database$X1 | sleepy_fish$Common_name %in% species_for_nocturnal_database$X1, c(1,6,8)] # 16:24 for references
# # sleepy_fish_sub$Diel_Pattern[grep("unclear", sleepy_fish_sub$Diel_Pattern)] <- "arrhythmic"
# # 
# # write.csv(sleepy_fish_sub, file = "species_for_nocturnal_database_filled.csv", row.names = F)
# 
# ###############################################################################################################################################
# ### FETCH DATA FROM TREE OF LIFE ### 
# ################################################################################################################################################
# 
# ## Remove fish without diel data, or real species names (necessary?)
# sleepy_fish <- sleepy_fish[sleepy_fish$Diel_Pattern != "",]
# sleepy_fish <- sleepy_fish[!(grepl("sp\\.|spp\\.|unidentified", sleepy_fish$Species)),]
# # sleepy_fish <- sleepy_fish[sleepy_fish$NEW == "NEW_SPECIES",]
# 
# ## Fetch species from tree of life using rotl package
# resolved_names <- tnrs_match_names(sleepy_fish$Species_name, context_name = "Vertebrates", do_approximate_matching = FALSE)
# 
# # ## Can use this to check those that don't have an exact match (in case there are new fish added)
# # resolved_names_2 <- tnrs_match_names(resolved_names$search_string[is.na(resolved_names$unique_name)], context_name = "Vertebrates", do_approximate_matching = TRUE)
# 
# ## Remove any that don't have exact matches
# resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
# 
# ## Remove ambiguous matches (approximate matches, things with multiple matches)
# resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]
# 
# ## Print the number of matches
# print(paste("rotl found matches for", nrow(resolved_names), "out of", nrow(sleepy_fish), "from the Sleepy fish database", sep = " "))
# 
# # Remove excess information, clean up, and add tip label ids that will match the tree
# # resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
# resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
# resolved_names <- resolved_names[!duplicated(resolved_names$tips),]
# 
# resolved_names$genus <- fishbase_df$Genus[match(resolved_names$unique_name, fishbase_df$Species)]
# resolved_names$family <- fishbase_df$Family[match(resolved_names$unique_name, fishbase_df$Species)]
# resolved_names$order <- fishbase_df$Order[match(resolved_names$unique_name, fishbase_df$Species)]
# 
# ## Add data on traits (from sleepy_fish, and/or from fishbase)
# 
# resolved_names$diel <- sleepy_fish$Diel_Pattern[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]
# resolved_names$diel <- factor(resolved_names$diel, levels = c("diurnal", "nocturnal", "unclear", "crepuscular", "crepuscular/diurnal", "crepuscular/nocturnal", "crepuscular/unclear", "unclear/diurnal", "unclear/nocturnal"))
# resolved_names$diel2 <- ifelse(resolved_names$diel == "diurnal", "diurnal", ifelse(resolved_names$diel == "nocturnal", "nocturnal", ifelse(resolved_names$diel == "crepuscular", "crepuscular", ifelse(resolved_names$diel == "crepuscular/diurnal", "crepuscular", ifelse(resolved_names$diel == "crepuscular/nocturnal", "crepuscular", ifelse(resolved_names$diel == "crepuscular/unclear", "crepuscular", ifelse(resolved_names$diel == "unclear/diurnal", "unclear", ifelse(resolved_names$diel == "unclear/nocturnal", "nocturnal", "unknown"))))))))
# resolved_names$diel_confidence <- sleepy_fish$Confidence[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]
# resolved_names$genome <- sleepy_fish$Genome[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]
# resolved_names$NEW <- sleepy_fish$NEW[match(resolved_names$search_string, tolower(sleepy_fish$Species_name))]
# 
# # These break tol_induced_subtree
# # c("incertae_sedis_inherited", "unplaced_inherited","incertae_sedis", "not_otu, incertae_sedis")
# # These do not break tol_induced_subtree, but do break congruify.phylo
# # c("infraspecific", "sibling_higher")
# # resolved_names <- resolved_names[resolved_names$flags %out% c("incertae_sedis_inherited", "unplaced_inherited", "incertae_sedis", "not_otu, incertae_sedis", "extinct_inherited, incertae_sedis"),]
# 
# setwd(here())
# 
# write.csv(resolved_names, file = here("resolved_names_local_2024-08-23.csv")) 
# 
# print(paste("Sleepy fish database covers ", round((length(unique(resolved_names$tips))/length(unique(fishbase_df$Species)))*100, 2), "% of Species, ", round((length(unique(resolved_names$genus))/length(unique(fishbase_df$Genus)))*100,2), "% of Genuses, ", round((length(unique(resolved_names$family))/length(unique(fishbase_df$Family)))*100,2), "% of Families, and ", round((length(unique(resolved_names$order))/length(unique(fishbase_df$Order)))*100,2), "% of Orders", sep = ""))

###############################################################################################################################################
### Determine missing clades ### 
################################################################################################################################################

# ## Figure out which Species, Genus, Family, Order, are missing from the datas
# '%out%' <- Negate('%in%')
# 
# species <- unique(fishbase_df$Species)[unique(fishbase_df$Species) %out% unique(resolved_names$unique_name)]
# genus <- unique(fishbase_df$Genus)[unique(fishbase_df$Genus) %out% unique(resolved_names$genus)]
# family <- unique(fishbase_df$Family)[unique(fishbase_df$Family) %out% unique(resolved_names$family)]
# order <- unique(fishbase_df$Order)[unique(fishbase_df$Order) %out% unique(resolved_names$order)]
# 
# # Find species which are from missing genuses, families, and orders
# # match the genus, then ask if it is in genus
# 
# species2 <- unique(species[fishbase_df$Genus[match(species, fishbase_df$Species)] %in% genus])
# species3 <- unique(species2[fishbase_df$Family[match(species2, fishbase_df$Species)] %in% family])
# species4 <- unique(species3[fishbase_df$Order[match(species3, fishbase_df$Species)] %in% order])
# 
# new_df <- data.frame(species = species4, genus = fishbase_df$Genus[match(species4, fishbase_df$Species)], family = fishbase_df$Family[match(species4, fishbase_df$Species)], order = fishbase_df$Order[match(species4, fishbase_df$Species)], data = species4 %in% sleepy_fish$Species.name)
# new_df_family <- unique(data.frame(genus = fishbase_df$Genus[match(species3, fishbase_df$Species)], family = fishbase_df$Family[match(species3, fishbase_df$Species)], order = fishbase_df$Order[match(species3, fishbase_df$Species)], data = species3 %in% sleepy_fish$Species.name))
# 
# new_df_family$scholar <- paste("https://scholar.google.com/scholar?start=0&q=%22", new_df_family$genus, "%22+AND+(%22circadian%22+OR+%22diel%22+OR+%22diurnal%22+OR+%22nocturnal%22+OR+%22crepuscular%22)&hl=en&as_sdt=0,5", sep = "")
# new_df_family$google <- paste("https://www.google.com/search?q=%22", new_df_family$genus, "%22+AND+%28circadian+OR+diel+OR+diurnal+OR+nocturnal+OR+crepuscular+%29&rlz=1C5CHFA_enCH822CH822&biw=1680&bih=831&sxsrf=AOaemvIvLVifugmAKccz0kD3lCJFQ9-laQ%3A1637151255625&ei=F_KUYYXQJYjjkgXZgb74BQ&oq=%22", new_df_family$genus, "%22+AND+%28circadian+OR+diel+OR+diurnal+OR+nocturnal+OR+crepuscular+%29&gs_lcp=Cgdnd3Mtd2l6EANKBAhBGAFQ5QJY5QJgpARoAnAAeACAAWqIAWqSAQMwLjGYAQCgAQKgAQHAAQE&sclient=gws-wiz&ved=0ahUKEwjFt6PYr5_0AhWIsaQKHdmAD18Q4dUDCA4&uact=5", sep = "")
# 
# library(clipr)
# write_clip(new_df_family)
# 
# # Checked all Gobiesociformes, Amiiformes

################################################################################################################################################
### FETCH AND TIME-CALIBRATE THE TREE ### 
################################################################################################################################################

resolved_names <- read.csv(file = here( here("resolved_names_local_2024-08-23.csv")), row.names = "X") 

# Fetch the combined tree from tree of life for the species ids found in resolved_names
# "sibling_higher" is the only flag that can be included where I can both fetch the tree and time calibrate it
'%out%' <- Negate('%in%')

# tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %out% c("incertae_sedis_inherited", "unplaced_inherited", "incertae_sedis", "not_otu", "not_otu, incertae_sedis", "extinct_inherited, incertae_sedis", "infraspecific")], label_format = "name")

tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") # I need to use the id option here, and then use that to map the tip labels from resolved_names (that way I don't run into the issue with the difference in formatting between the two tools)

# Time calibrate it using geiger and timetree.org
# First resolve polytomies ~randomly using multi2dr

tr <- multi2di(tr)

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

saveRDS(reference.df, file = here("reference_df_2024-08-23.rds"))

# Load the timetree tree (genus level data works, but not species)
# Have download timetree data for species, genus, family, and order
# Genus level data has the most calibration points
setwd(here())

timetree_order <- ape::read.tree(here("timetree_data/actinopterygii_order.nwk"))
timetree_family <- ape::read.tree(here("timetree_data/actinopterygii_family.nwk"))
timetree_genus <- ape::read.tree(here("timetree_data/actinopterygii_genus.nwk"))

# Use geiger to congruify the tree, works with treePL
# This seems to work up to genus, but not species (by replacing tip.labels with the same names)
setwd(here())
geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd(here())
geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd(here())
geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

tr.calibrated <- geiger.genus$phy

## Save out files

saveRDS(tr.calibrated, file = "calibrated_phylo_2024-08-23.rds")

# Add in a stop here

print("this is the last message")
stop()
print("you should not see this")

# ## load back in to do by species manually
# 
# reference.df <- readRDS(file = "reference_df_2024-08-23.rds")
# 
# tr.calibrated <- readRDS("calibrated_phylo_2024-08-23.rds")
# 
# # Below works if you modify the heights.phylo function
# trace(geiger:::heights.phylo, edit = TRUE)
# # depth = max(xx[!(is.na(xx))])
# # Also have to do it manually - which is ugh!
# timetree_species <- ape::read.tree(here("timetree_data/actinopterygii_species.nwk"))
# timetree_species <- multi2di(timetree_species)
# setwd(here())
# 
# geiger.species <- congruify.phylo(reference = timetree_species, target = tr.calibrated, taxonomy = reference.df, tol = 0, scale = "treePL")
# tr.calibrated <- geiger.species$phy
# 
# tr.calibrated$tip.label <- resolved_names$tips[match(tr.calibrated$tip.label, resolved_names$ott_id)]
# 
# saveRDS(tr.calibrated, file = "tr_tree_calibrated_fish_2024-08-23.rds")
# 
# 
# ## Generate subset with Nocturnal vs Diurnal
# 
# trait.data <- data.frame(species = tr.calibrated$tip.label, diel = resolved_names$diel[match(tr.calibrated$tip.label, resolved_names$tips)]) # OK, some species tip labels are more complicated and cause issues here
# 
# # Create vectors including crepuscular/unclear, or not
# trait.data$diel1 <- ifelse(trait.data$diel %in% c("diurnal", "crepuscular/diurnal", "unclear/diurnal"), "diurnal", ifelse(trait.data$diel %in% c("nocturnal", "crepuscular/nocturnal", "unclear/nocturnal"), "nocturnal", "unclear/crepuscular"))
# levels(trait.data$diel1) <- c("diurnal", "nocturnal", "unclear/crepuscular")
# trait.data$diel2 <- ifelse(trait.data$diel %in% c("diurnal"), "diurnal", ifelse(trait.data$diel %in% c("nocturnal"), "nocturnal", ifelse(trait.data$diel %in% c("crepuscular", "crepuscular/diurnal", "crepuscular/nocturnal"), "crepuscular", "unclear")))
# levels(trait.data$diel2) <- c("diurnal", "nocturnal", "crepuscular", "unclear")
# 
# trait.data <- trait.data[!(is.na(trait.data$diel)),]
# rownames(trait.data) <- trait.data$species
# 
# trait.data$tips <- resolved_names$tips[match(trait.data$species, resolved_names$tips)]
# trait.data$order <- resolved_names$order[match(trait.data$species, resolved_names$tips)]
# trait.data$confidence <- resolved_names$diel_confidence[match(trait.data$species, resolved_names$tips)]
# trait.data$NEW <- resolved_names$NEW[match(trait.data$species, resolved_names$tips)]
# 
# trait.data$crepuscular <- ifelse(trait.data$diel2 == "crepuscular", "crepuscular", ifelse(trait.data$confidence > 4, "non_crepuscular", NA))
# trait.data$diel_continuous <- ifelse(trait.data$diel1 == "diurnal", trait.data$confidence, ifelse(trait.data$diel1 == "nocturnal", trait.data$confidence*-1, NA))
# 
# # # Add fresh/marine to trait.data
# # fishbase_ecosystem <- ecosystem() # Salinity is here duh
# # trait.data$marine <- fishbase_ecosystem$Salinity[match(gsub("_", " ", trait.data$species), fishbase_ecosystem$Species)]
# 
# # Add acanthomorpha
# acanthomorpha <- extract.clade(tr.calibrated, node = getMRCA(tr.calibrated, tip = c("Saccogaster_melanomycter", "Apolemichthys_xanthopunctatus")))
# cartilagenous <- extract.clade(tr.calibrated, node = getMRCA(tr.calibrated, tip = c("Rhizoprionodon_terraenovae", "Rhynchobatus_djiddensis")))
# trait.data$acanthomorpha <- ifelse(trait.data$tips %in% acanthomorpha$tip.label, "acanthomorpha", ifelse(trait.data$tips %in% cartilagenous$tip.label, "cartilagenous", "outgroup"))
# 
# saveRDS(trait.data, file = "trait_data_fish_2024-08-23.rds")
# 
# 
# # trait.data$FeedingType <- fishbase_ecology$FeedingType[match(gsub("_", " ", trait.data$species), fishbase_ecology$Species)]









