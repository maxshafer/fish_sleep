library(ape) 
library(corHMM)
library(xlsx)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(patchwork)
library(ggpubr)
library(dplyr)
library(phytools)
library(rfishbase)
library(geiger)

setwd(here())

source("scripts/fish_sleep_functions.R")


###### For cetaceans

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit#gid=0'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
cetaceans_full$Diel_Pattern <- tolower(cetaceans_full$Diel_Pattern)

cetaceans_full <- cetaceans_full[,1:7]
cetaceans_full$diel <- ifelse(grepl("crepuscular", cetaceans_full$Diel_Pattern), "crepuscular", ifelse(grepl("nocturnal", cetaceans_full$Diel_Pattern), "nocturnal", ifelse(grepl("diurnal", cetaceans_full$Diel_Pattern), "diurnal", ifelse(grepl("cathemeral", cetaceans_full$Diel_Pattern), "cathemeral", NA))))

## Probably should save out a local copy in case google goes bankrupt
# write.csv(cetaceans, file = here("sleepy_fish_database_local.csv"))


## Remove species without diel data
cetaceans <- cetaceans_full[!(is.na(cetaceans_full$diel)),]


## Read in the mammalian phylogeny
mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
mam.tree <- maxCladeCred(mammal_trees, tree = TRUE)

trait.data <- resolved_names[resolved_names$diel %in% c("diurnal", "nocturnal"),]
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
row.names(trait.data) <- trait.data$tips
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)



model<- corHMM(phy = trpy_n_mam, data = trait.data[trpy_n_mam$tip.label, c("tips", "diel")], rate.cat = 1, model = "ARD", node.states = "marginal")

#models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
#model <- models$HMM_2state_2rate_marg

lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc) <- c("diurnal", "nocturnal")

lik.anc$node <- c(1:length(trpy_n_mam$tip.label), (length(trpy_n_mam$tip.label) + 1):(trpy_n_mam$Nnode + length(trpy_n_mam$tip.label)))

ancestral_plot <- ggtree(trpy_n_mam, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")


ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
