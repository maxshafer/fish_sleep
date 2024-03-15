# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
#library(scales)
library(gsheet)
#library(patchwork)
#library(ggpubr)
library(dplyr)
library(phytools)
#library(rfishbase)
library(geiger)
library(ggplot2)
install.packages("rphylopic")
library(rphylopic)
library(RColorBrewer)
install.packages("ggimage")
library(ggimage)

setwd(here())

source("scripts/fish_sleep_functions.R")


# Section 1: Load in and clean up cetacean data ---------------------------

###### For cetaceans

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit#gid=0'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#View(cetaceans_full)

#set all diel pattern entries to lower case, helps keep consistency
cetaceans_full$Diel_Pattern_1 <- tolower(cetaceans_full$Diel_Pattern_1)
cetaceans_full$Diel_Pattern_2 <- tolower(cetaceans_full$Diel_Pattern_2)
cetaceans_full$Diel_Pattern_3 <- tolower(cetaceans_full$Diel_Pattern_3)

#cycle through diel pattern 1, 2, 3 to check that all entries are correctly spelled and see spread of data
#20 cath sps, 20 crep, 15 di, 25 noc
table(cetaceans_full$Diel_Pattern_3)

#add another row "tips" with the species name formatted as they appear in the tree (mam tree and otl)
#saves you from remaking it for each of the trait.data dataframes
cetaceans_full$tips <- cetaceans_full$Species_name
cetaceans_full$tips <- str_replace(cetaceans_full$tips, pattern = " ", replacement = "_")


# Section 2: Create the basic tree ----------------------------------------


#create trait.data dataframe to use to build tree and add trait data (in geom tile)
#check that all the entries have a value for diel pattern (I use diel pattern 3 instead of 1 because di/noc has NA for all cath sps)
#could also do !(is.na) for this
#for trait data only needs the diel patterns, species names, tips (unless we want to add more data to the tree like confidence etc)
trait.data <- cetaceans_full[cetaceans_full$Diel_Pattern_3 %in% c("diurnal", "nocturnal", "cathemeral", "crepuscular"),]
#for trait data only needs the diel patterns, species names, tips (unless we want to add more data to the tree like confidence etc)
trait.data <- trait.data[, c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "tips")]
#trait data has all the behavioural data for the 80 species with data 

#find the species names that are in the open tree of life
#looks for species in the same format as Species_names
resolved_names <- tnrs_match_names(names = trait.data$Species_name, context_name = "Vertebrates", do_approximate_matching = FALSE)

#Balaenoptera riceii, Inia humboldtiana and Neophocaena sunameri are not currently in the open tree of life
#20 species have no information and are also left out of the tree
#this leaves 75 species in our tree (and in resolved names) to work with

# Remove any that don't have exact matches or are synonyms (this is really just for finding ancestral state, so missing a few species won't matter)
resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
resolved_names <- resolved_names[resolved_names$is_synonym == FALSE,]
resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]

# Remove excess information, clean up, and add tip label ids that will match the tree
resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
resolved_names <- resolved_names[!duplicated(resolved_names$tips),]

# Add data on activity
# match() finds the position of its first argument in its second
resolved_names$diel <- trait.data$Diel_Pattern_3[match(resolved_names$search_string, tolower(trait.data$Species_name))]
table(resolved_names$diel)


## Fetch the tree

tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") # I need to use the id option here, and then use that to map the tip labels from resolved_names (that way I don't run into the issue with the difference in formatting between the two tools)

# Time calibrate it using geiger and timetree.org
# First resolve polytomies ~randomly using multi2dr

tr <- multi2di(tr)

tr$tip.label <- resolved_names$tips[match(tr$tip.label, paste("ott", resolved_names$ott_id, sep = ""))]
#drop tips with NA values leftover from the strictly diurnal/nocturnal column
to_drop <- trait.data[is.na(trait.data$Diel_Pattern_1), ]
to_drop$tips
tr2 <- drop.tip(tr, to_drop$tips)

ggtree(tr2) + geom_tiplab()

#pick colours! Using a custom colour palette for consistency between the plots
display.brewer.all(type = "qual", colorblindFriendly = TRUE)
display.brewer.pal(8, "Set2")
custom.colours1 <- c("#FC8D62","#66C2A5")

diel.plot1 <- ggtree(tr2, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_1")]
diel.plot1 <- diel.plot1 + geom_tile(data = diel.plot1$data[1:length(tr2$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_1), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours1)
diel.plot1.labelled <- diel.plot1 + geom_tiplab(color = "black", size = 1)
diel.plot1.labelled

#export plot as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel1.png", width=11,height=8,units="cm",res=1200)
print(diel.plot1.labelled)
dev.off()

custom.colours2 <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

diel.plot2 <- ggtree(tr, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2")]
diel.plot2 <- diel.plot2 + geom_tile(data = diel.plot2$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_2), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours2)
diel.plot2.labelled <- diel.plot2 + geom_tiplab(color = "black", size = 1)
diel.plot2.labelled

png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel2.png", width=14,height=10,units="cm",res=1200)
print(diel.plot2.labelled)
dev.off()

custom.colours3 <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
trait.data$random <- as.character(sample(1:4,80, replace=T))
diel.plot3 <- ggtree(tr, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_3")]
diel.plot3 <- diel.plot3 + geom_tile(data = diel.plot3$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_3, inherit.aes = FALSE, color = "transparent")) #+ scale_fill_manual(values = custom.colours3))
diel.plot3.labelled <- diel.plot3 + geom_tiplab(color = "black", size = 1)
diel.plot3.labelled


png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel3.png", width=14,height=10,units="cm",res=1200)
print(diel.plot3.labelled)
dev.off()

#diel with three states
trait.data3state <- cetaceans[!(is.na(cetaceans$Diel_Pattern_2)),]
trait.data3state$Diel_Pattern_2 <- str_replace_all(trait.data3state$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data3state$Diel_Pattern_2 <- str_replace_all(trait.data3state$Diel_Pattern_2, "diurnal/crepuscular", "nocturnal")
trait.data3state$tips <- trait.data3state$Species_name
trait.data3state$tips <- str_replace(trait.data3state$tips, pattern = " ", replacement = "_")
custom.colours4 <- c("#dd8ae7", "#FC8D62", "#66C2A5")

diel.plot3_state <- ggtree(tr, layout = "circular") %<+% trait.data3state[,c("tips", "Diel_Pattern_2")]
diel.plot3_state <- diel.plot3_state + geom_tile(data = diel.plot3_state$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_2), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours4)
diel.plot3_state <- diel.plot3_state + geom_tiplab(color = "black", size = 1.2)
diel.plot3_state

png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel4.png", width=15,height=12,units="cm",res=1200)
print(diel.plot3_state)
dev.off()

#add clade labels to find the nodes for each family
node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node), hjust=-.2, size = 1.8) + geom_tiplab(size = 1.8, hjust = -0.1)
node_labels 

#create a dataframe with each cetacean family and its node
family_names_right <- data.frame(names = c("Eschrichitiidae", "Neobalaenidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Platanistidae", "Kogiidae", "Physteridae"), nodes = c(72, 73, 151, 143, 88, 81, 139, 60))
family_names_left <- data.frame(names = c("Ziphiidae", "Lipotidae", "Phocoenidae", "Iniidae", "Monodontidae"), nodes = c(128, 3, 123, 85, 127))
#add the colour you want each family to be
family_names_right$colour <- c("slateblue3", "royalblue4", "darkslateblue", "slateblue1", "cyan1", "dodgerblue3","royalblue3", "blueviolet" )
family_names_left$colour <- c("royalblue1", "#1899f5", "#19e2ff","deepskyblue", "turquoise2")

diel_families <- diel.plot2 + geom_tiplab(color = "black", size = 1.5) +
geom_cladelab(family_names_right, node = family_names_right$nodes, label = family_names_right$names, offset=11, align=FALSE, angle=12, offset.text=1, barsize=1, fontsize=2.5, textcolor=family_names_right$colour, barcolor=family_names_right$colour, alpha = 0.75) +
geom_cladelab(family_names_left, node = family_names_left$nodes, label = family_names_left$names, offset=11, align=FALSE, angle=192, offset.text=1, barsize=1, fontsize=2.5, textcolor= family_names_left$colour, barcolor= family_names_left$colour, alpha = 0.75)
#does this actually make the plot smaller?
diel_families <- diel_families + theme(plot.margin = unit(c(14,8,14,8), "mm"))
#geom_cladelab doesn't give a bar to clades with single sps so add them manually
diel_families <- diel_families + geom_strip(60, 61, barsize=1, color='blueviolet', offset = 11)
diel_families <- diel_families + geom_strip(73, 72, barsize=1, color='royalblue4', offset = 11)
diel_families <- diel_families + geom_strip(63, 72, barsize=1, color='slateblue3', offset = 11)
diel_families <- diel_families + geom_strip(3, 6, barsize=1, color='#1899f5', offset = 11)
diel_families
                          
png("C:/Users/ameli/OneDrive/Documents/R_projects/diel3_with_families.png", width=30,height=20,units="cm",res=1200)
print(diel_families)
dev.off()

#better family plot
#create a dataframe with each cetacean family and its node
family_names_right <- data.frame(names = c("Eschrichitiidae", "Neobalaenidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Platanistidae", "Kogiidae", "Physteridae"), nodes = c(72, 73, 151, 143, 88, 81, 139, 60))
family_names_left <- data.frame(names = c("Ziphiidae", "Lipotidae", "Phocoenidae", "Iniidae", "Monodontidae"), nodes = c(128, 3, 123, 85, 127))
#add the colour you want each family to be
family_names_right$colour <- c("steelblue2", "steelblue3", "steelblue4", "steelblue1", "grey1", "grey50","grey60", "grey65" )
family_names_left$colour <- c("grey45", "grey40", "grey15","grey30", "grey20")

#odontocetes in grey, mysticetes in blue
family_names_right$parvorder <- "grey40"
family_names_right[c(1,2,3,4), "parvorder"] <- "steelblue2"
family_names_left$parvorder <- "grey40"

diel_families <- diel.plot3 +
  geom_cladelab(family_names_right, node = family_names_right$nodes, label = family_names_right$names, offset=1, align=FALSE, angle=16, offset.text=1, barsize=1, fontsize=5, textcolor= family_names_right$colour, barcolor=family_names_right$parvorder, alpha = 0.75) +
  geom_cladelab(family_names_left, node = family_names_left$nodes, label = family_names_left$names, offset=1, align=FALSE, angle=196, offset.text=1, barsize=1, fontsize=5, textcolor= family_names_left$colour, barcolor= family_names_left$parvorder, alpha = 0.75)
#does this actually make the plot smaller?
diel_families 
#geom_cladelab doesn't give a bar to clades with single sps so add them manually
diel_families <- diel_families + geom_strip(60, 61, barsize=1, color='grey50', offset = 1, align=FALSE)
diel_families <- diel_families + geom_strip(73, 72, barsize=1, color='steelblue1', offset = 1)
diel_families <- diel_families + geom_strip(63, 72, barsize=1, color='steelblue1', offset = 1)
diel_families <- diel_families + geom_strip(3, 6, barsize=1, color='grey50', offset = 1)
diel_families

png("C:/Users/ameli/OneDrive/Documents/R_projects/diel3_with_families.png", width=30,height=20,units="cm",res=1200)
print(diel_families)
dev.off()
#find images from phylopics collection, not necessary to rerun once you have the image info
# img <- pick_phylopic(name = "Mysticeti", n = 19, view = 19)
# img <- pick_phylopic(name = "Delphinida", n = 54, view = 54)
# img <- pick_phylopic(name = "Platanistoidea", n = 3, view = 3)
# img <- pick_phylopic(name = "pan-Physeteroidea", n = 10, view = 10)
# img <- pick_phylopic(name = "Ziphioidea", n = 20, view = 20)
# img <- pick_phylopic(name = "Lipotidae", n = 1, view = 1)
# img <- pick_phylopic(name = "Phocoena", n = 4, view = 4)
# img <- pick_phylopic(name = "Inia", n = 1, view = 1)
# img <- pick_phylopic(name = "Monodon", n = 2, view = 2)


#add image names into the dataframe
family_names_pics <- data.frame(names = c("Eschrichitiidae", "Neobalaenidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Platanistidae", "Kogiidae", "Physteridae", "Ziphiidae", "Lipotidae", "Phocoenidae", "Iniidae", "Monodontidae"), 
                                nodes = c(70, 71, 147, 138, 86, 79, 135, 58, 125, 3, 121, 83, 124),
                                images = c("8b73f54f-15e8-41b8-8c9c-46c86a185104", "941169f7-bc86-4030-bb28-4419b78214de", "fdff3c1b-dd0a-44d5-9fd3-2ecda7939846", "012afb33-55c3-4fc6-9ae3-3a91fda32fd5", "3caf4fbd-ca3a-48b4-925a-50fbe9acd887", "a55581b9-72c9-4ede-8fba-b908b08d94c9", "5bfb840e-071f-4a1a-b101-0a747a5453e7", "dc76cbdb-dba5-4d8f-8cf3-809515c30dbd", "2fad47bd-9a34-4bc1-b6cc-b7c4415b109d", "103dbb51-1bce-4df6-b1ac-ff0d38d188a3", "970f7e8d-a823-45b9-85b9-767704d0c13f", "f36c9daa-a102-42dd-88ac-a126753943d2", "bfe45de4-d8a8-423c-abf3-087a7c7d0d6c"))

#map these onto the phylogeny
diel3_fam_pics <- diel.plot3.labelled + geom_cladelab(data = family_names_pics, 
                        mapping = aes(node = nodes, label = names, image = images), 
                        geom = "phylopic", imagecolor = "grey30", 
                        offset=11, offset.text=0.5)

diel3_fam_pics


cetaceans_full <- cetaceans_full[,1:9]
View(cetaceans_full)
## Probably should save out a local copy in case google goes bankrupt
write.csv(cetaceans_full, file = here("sleepy_fish_database_local.csv"))


#doing the same as above but with Cox's mammal tree instead of open tree of life
#need to use the mamm tree for modelling and ancestral trait reconstruction

# ## Read in the mammalian phylogeny
# mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
# mam.tree <- maxCladeCred(mammal_trees, tree = TRUE)
# 
# trait.data <- cetaceans[cetaceans$Diel_Pattern_1 %in% c("diurnal", "nocturnal"),]
# trait.data$tips <- trait.data$Species_name
# trait.data$tips <- str_replace(trait.data$tips, pattern = " ", replacement = "_")
# trait.data2 <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
# 
# "Lagenorhynchus_obscurus" %in% mam.tree$tip.label
# 
# View(trait.data$tips)
# View(mam.tree$tip.label)
# mam.tree$tip.label[(grepl("Balaenoptera", mam.tree$tip.label))]
# 
# #find which species in trait data are missing from mam.tree
# trait.data$tips[!(trait.data$tips %in% mam.tree$tip.label)]
# 
# row.names(trait.data) <- trait.data$tips
# trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# #change this to work with the open tree of life
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# model<- corHMM(phy = trpy_n_mam, data = trait.data[trpy_n_mam$tip.label, c("tips", "diel")], rate.cat = 1, model = "ARD", node.states = "marginal")
# 
# #models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
# #model <- models$HMM_2state_2rate_marg
# 
# lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
# colnames(lik.anc) <- c("diurnal", "nocturnal")
# 
# lik.anc$node <- c(1:length(trpy_n_mam$tip.label), (length(trpy_n_mam$tip.label) + 1):(trpy_n_mam$Nnode + length(trpy_n_mam$tip.label)))
# 
# ancestral_plot <- ggtree(trpy_n_mam, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
# 
# 
# ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)


