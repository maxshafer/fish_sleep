# Section 0: Packages -----------------------------------------------------
#packages we will use
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
library(ggtree)
#manipulating dataframes
library(dplyr)
#install.packages("tictoc")
library(tictoc)

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
source("scripts/Amelia_functions.R")

# Section 1: Arguments ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(!(args[1] %in% c("max_crep", "max_dinoc", "four_state_max_crep", "four_state_max_crep", "six_state"))) {  
  stop("first argument must be states in the model")
}

if(!(args[2] %in% c("cetaceans", "whippomorpha", "artiodactyla", "artiodactyla_minus_cetaceans", "ruminants"))) {  
  stop("second argument must be the phylogenetic tree")
}

if(length(args) < 3) {  
  stop("must supply at least one model")
}  


# Section 2: Create trait data  --------
#we want to run select models on the 1k possible trees 
#allow us to better compare likelihoods for significant differences

#want to compare four models of the trinary cetacean data (cathemeral, diurnal, nocturnal)
#ER, SYM, ARD, cathemeral/crepuscular bridge

#first we'll import the trait data dataframe and subset the tree to include cetacean species only
if(args[2] == "cetaceans"){
  trait.data <- read.csv(here("cetaceans_full.csv"))
}

if(args[2] == "whippomorpha"){
  trait.data <- read.csv(here("whippomorpha.csv"))
}


if(args[2] == "artiodactyla"){
  trait.data <- read.csv(here("Sleepy_artiodactyla_full.csv"))
}

if(args[2] == "artiodactyla_minus_cetaceans"){
  trait.data <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv"))
}

if(args[2] == "ruminants"){
  trait.data <- read.csv(here("ruminants_full.csv"))
}


if(args[1] == "max_crep"){
  trait.data <- trait.data[,c("Diel_Pattern_2", "tips")]
  trait.data$Diel_Pattern_2 <- tolower(trait.data$Diel_Pattern_2)
  trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)),]
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "nocturnal/crepuscular", "crep_cath")
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "diurnal/crepuscular", "crep_cath")
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "cathemeral/crepuscular", "crep_cath")
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "cathemeral", "crep_cath")
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "crep_cath", "crepuscular_cathemeral")
  
}
if(args[1] == "max_dinoc"){
  trait.data <- trait.data[,c("Diel_Pattern_2", "tips")]
  trait.data$Diel_Pattern_2 <- tolower(trait.data$Diel_Pattern_2)
  trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)),]
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")
  trait.data$Diel_Pattern_2 <- str_replace_all(trait.data$Diel_Pattern_2, "cathemeral/crepuscular", "cathemeral")
  trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "nocturnal", "cathemeral"),]
}

if(args[1] == "four_state_max_crep"){
  trait.data <- trait.data[, c("tips", "Diel_Pattern_2")]
  trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, "nocturnal/crepuscular", "crepuscular")
  trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, "diurnal/crepuscular", "crepuscular")
  trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, "cathemeral/crepuscular", "crepuscular")
}

if(args[1] == "four_state_max_dinoc"){
  trait.data <- trait.data[, c("tips", "Diel_Pattern_2")]
  trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
  trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")
  trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, "cathemeral/crepuscular", "crepuscular")
}

if(args[1] == "six_state"){
  trait.data <- trait.data[, c("tips", "Diel_Pattern_2")]
  trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)),]
}

# Section 3: Subset the trees ---------------------------------------------
#we need to subset by species names in trait data (only species with behavioural data)
#this takes 15-25 seconds 
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))

#subset cetacean trees for now to see if it runs
mammal_trees <- mammal_trees[sample(1:length(mammal_trees), 100)]

phylo_trees <- lapply(mammal_trees, function(x) subsetTrees(tree = x, subset_names = trait.data$tips))
   
#subset trait data to only include species that are in the tree
trait.data <- trait.data[trait.data$tips %in% phylo_trees[[2]]$tip.label,]


# Section 4: Use returnCorModels to run corHMM models (ER, SYM, ARD, and/or bridge_only) on 1k possible trees --------

for(i in args[-(1:2)]){
  if(!(i %in% c("ER", "SYM", "ARD", "bridge_only", "CONSYM"))){
    stop("model not recognized")
  }
}

if("ER" %in% args){
  ER <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1, custom.rate.mat = "none", model = "ER", node.states = "marginal"))
}

if("SYM" %in% args){
  SYM <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1, custom.rate.mat = "none",  model = "SYM", node.states = "marginal"))
  
}

if("ARD" %in% args){
  ARD <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1,  custom.rate.mat = "none", model = "ARD", node.states = "marginal"))
  
}

# if("bridge_only" %in% args & !("six_state" %in% args)){
#   bridge_only <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1, custom.rate.mat = matrix(c(0,2,3,4,0,0,7,0,0), ncol = 3, nrow = 3, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)"), c("(1, R1)", "(2, R1)", "(3,R1)"))), model = "ARD", node.states = "marginal"))
#   
# }

# #don't allow transitions from noc -> di, noc -> di/crep, noc/crep -> di, noc/crep -> di/crep, di -> noc, di -> noc/crep, di/crep -> noc, di/crep -> noc/crep 
# generic_ratemat <- getStateMat4Dat(ARD$data)$rate.mat
# #need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> nocturnal/3 (6)
# custom_rate_matrix <- dropStateMatPars(generic_ratemat, c(14, 15, 19, 20, 23, 24, 28, 29))

#for three state
if("bridge_only" %in% args & "max_dinoc" %in% args){
  bridge_only <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1,  custom.rate.mat = matrix(c(0,2,3,4,0,0,7,0,0), ncol = 3, nrow = 3, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)"), c("(1, R1)", "(2, R1)", "(3,R1)"))), model = "ARD", node.states = "marginal"))
}

if("bridge_only" %in% args & "max_crep" %in% args){
  bridge_only <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1,  custom.rate.mat = matrix(c(0,2,3,4,0,0,7,0,0), ncol = 3, nrow = 3, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)"), c("(1, R1)", "(2, R1)", "(3,R1)"))), model = "ARD", node.states = "marginal"))
}

#for four state
if("bridge_only" %in% args & "four_state_max_dinoc" %in% args){
  bridge_only <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1,  custom.rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), model = "ARD", node.states = "marginal"))
}

if("bridge_only" %in% args & "four_state_max_crep" %in% args){
  bridge_only <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1, custom.rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), model = "ARD", node.states = "marginal"))
}

#try running a constrained bridge only model based on a symmetrical model
if("CONSYM" %in% args & "four_state_max_crep" %in% args){
  bridge_only <- lapply(phylo_trees, function(x) returnCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1, custom.rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), model = "SYM", node.states = "marginal"))
}

# Section 5: Save the results out and extract likelihoods  --------
#use paste() to create a filename with all of the arguments
result_list <- lapply(args[-(1:2)], function(x) eval(as.name(x)))
names(result_list) <- paste(args[-(1:2)], "_model", sep = "")

saveRDS(result_list, paste(args[2], args[1], "traits", paste0(args[-(1:2)], sep = "", collapse = "_"), "models", sep = "_"))
