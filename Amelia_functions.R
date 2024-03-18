##Packages we will use ---------------------------------------------------
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)

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




## Read in the mammalian phylogeny
mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")

length(mammal_trees)
mammal_trees[[1]]


#first custom function SubsetTrees takes the 1k trees and the subset of mammals (cetaceans, artio)
#returns the 1k subsetted trees

subsetTrees <- function(tree = mammal_trees[[1]], subset_names = cetaceans_full$tips) {
  # will take a tree and keep only those tips that match subset_names
  
  subset_names <- subset_names[subset_names %in% tree$tip.label]
  
  out_tree <- keep.tip(tree, subset_names)
  
  return(out_tree)
  
}

cetacean_trees <- lapply(mammal_trees, function(x) subsetTrees(tree = x, subset_names = cetaceans_full$tips))



#second custom function returnModels that takes the subsetted 1k trees, the input data (diel patterns +species names) and the model type you want to run (ER, SYM, ARD)
#returns the model results for each of the 1000 trees (list of 1000 trees)

returnAceModels <- function(tree = cetacean_trees[[1]], trait.data = trait_data, column = "Diel_Pattern_1", model = "SYM") {
  trait.data <- trait.data[trait.data$tips %in% tree$tip.label,]
  trait.data <- trait.data[!(is.na(trait.data[,column])),]
  tree <- keep.tip(tree, trait.data$tips)
  ace_model <- ace(trait.data[,column], tree, model = model, type = "discrete")
  return(ace_model)
  
}

model_test <- returnAceModels(tree = cetacean_trees[[2]], trait.data = cetaceans_full, column = "Diel_Pattern_1", model = "SYM" )

cetacean_sim <- lapply(cetacean_trees, function(x) returnAceModels(tree = x, trait.data = cetaceans_full, column = "Diel_Pattern_1", model = "SYM"))

#third custom function returnModels that takes the subsetted 1k trees, the input data (diel patterns +species names) and the model type you want to run (ER, SYM, ARD)
#same as before but for the cor model
#returns results for each of the 1000 trees as a list

returnCorModels <- function(tree = cetacean_trees[[2]], trait.data = cetaceans_full, species_col = "species", diel_col = "Diel_Pattern_1", rate.cat = 1, model = "SYM", node.states = "marginal"){
  trait.data <- trait.data[trait.data$tips.label %in% tree$tip.label,]
  trait.data <- trait.data[!(is.na(trait.data[,diel_col])),]
  trait.data <- trait.data[tree$tip.label, c(species_col, diel_col)]
  tree <- keep.tip(tree, trait.data$species_col)
  cor_model <- corHMM(phy = tree, data = trait.data, rate.cat = rate.cat, model = model, node.states = node.states)
  return(cor_model)
}

test_time <- returnCorModels(cetacean_trees[[2]], trait.data = cetaceans_full, species_col = "species", diel_col = "Diel_Pattern_1", rate.cat = 1, model = "ER", node.states = "marginal")

#fourth custom function returnLikelihoods takes model results from returnModels and returns just the likelihoods of those models 

returnLikelihoods <- function(model = cetacean_sim[[1]], return = "loglik"){
  
  return(model[return])
}

cet_likelihoods <- unlist(lapply(cetacean_sim, function(x) returnLikelihoods(model = x)))


