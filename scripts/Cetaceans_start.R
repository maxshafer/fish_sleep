# For retreiving data from the open tree of life
library(rotl)
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)

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


###### For cetaceans

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit#gid=0'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
cetaceans_full$Diel_Pattern <- tolower(cetaceans_full$Diel_Pattern)

cetaceans_full <- cetaceans_full[,1:9]
cetaceans_full$diel <- tolower(cetaceans_full$Diel_Pattern_3)


## Probably should save out a local copy in case google goes bankrupt
# write.csv(cetaceans, file = here("sleepy_fish_database_local.csv"))


## Remove species without diel data
cetaceans <- cetaceans_full[!(is.na(cetaceans_full$diel)),]


## Read in the mammalian phylogeny
mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
mam.tree <- maxCladeCred(mammal_trees, tree = TRUE)

# Save out maxcladecred, so we don't have to recalculate it every time
saveRDS(mam.tree,"maxCladeCred_mammal_tree.rds")

mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")

## Quick change to diurnal / nocturnal only trait.data
## I use standard nomenclature for trees ('trpy_n') and trait data ('trait.data'), so that the script below can be run with any other tree etc

# This selects only data that is diurnal or nocturnal
trait.data <- cetaceans_full[cetaceans_full$diel %in% c("diurnal", "nocturnal"),]
# selects only data that is in the mammal tree
trait.data$tips <- str_replace(trait.data$Species_name, " ", "_")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
row.names(trait.data) <- trait.data$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)


## You should double check which species exist or not in your tree (I think we did this already)
## For the below, you can try all different combinations (diurnal/nocturnal only, or include crepuscular or cathemeral (or not), and see what difference it makes)

## Try simple models first, ER and ARD models
## This uses the ace package, but you can run these models with any of the libraries
## They have slightly different inputs, for example, ace wants a vector of traits
trait.vector <- trait.data$diel
ace_model_er <- ace(trait.vector, trpy_n, model = "ER", type = "discrete")
ace_model_sym <- ace(trait.vector, trpy_n, model = "SYM", type = "discrete")
ace_model_ard <- ace(trait.vector, trpy_n, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package

cor_model_er <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 1, model = "ARD", node.states = "marginal")


## Extract the data to plot
## Depending on which package you use, the ancestral data is in different places (or not determined at all)
## It typically only includes the internal node states, and not also the tip states
str(ace_model_er, max.level = 1)
head(ace_model_er$lik.anc)

str(cor_model_er, max.level = 1)
head(cor_model_er$states)

# In this case, I am putting together the internal nodes and tip states into one data frame
lik.anc <- as.data.frame(rbind(cor_model_ard$tip.states, cor_model_ard$states))
# dim of this should be equal to the tips and internal nodes
trpy_n
dim(lik.anc)

colnames(lik.anc) <- c("diurnal", "nocturnal")

# Tips are also nodes, and all nodes have a number, and they are number sequentially beginning with tips, so the first internal node is the number of tips +1

lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

ancestral_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")


ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
