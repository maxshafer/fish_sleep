library(ape) 
library(corHMM)
library(xlsx)
library(phangorn)
library(stringr)
library(here)

setwd(here())

## Read in data
tet_data <- read.csv("tetrapod_data/suppfile10appendix1914.csv")
# tet_data <- tet_data[tet_data$Species %in% tet_data$Species[tet_data$Class != "Mammalia"],]
tet_data$group <- tet_data$Class
tet_data$Genus <- gsub("\\_.*","",tet_data$Species)
tet_data <- tet_data[,c("Species", "State", "Genus", "Family", "Order", "group")]
colnames(tet_data) <- c("species", "diel", "genus", "family", "order", "group")

## Convert to my diel versions
tet_data$diel <- c("unclear", "diurnal", "nocturnal", "crepuscular")[match(tet_data$diel, c("ARR", "DIU", "NOC", "CRE"))]
tet_data$diel2 <- tet_data$diel
# Load tree
tet_tree <- read.nexus("tetrapod_data/Hackett1914.tree")


# Reciprocally trim

trait.data <- tet_data[tet_data$species %in% tet_tree$tip.label,]
colnames(trait.data) <- c("species", "diel1", "genus", "family", "order", "group", "diel2")
rownames(trait.data) <- trait.data$species
  
trpy_n_tet <- keep.tip(tet_tree, tip = trait.data$species)

### Save out trait data and tree
### This is the Maximum Clade Credibility tree, and associated Diurnal/Nocturnal trait data (mutually exclusive extant species)
dataset <- "tetrapods"

saveRDS(trpy_n_tet, here(paste("tr_tree_calibrated_", dataset, ".rds", sep = "")))

saveRDS(trait.data, here(paste("trait_data_", dataset, ".rds", sep = "")))





