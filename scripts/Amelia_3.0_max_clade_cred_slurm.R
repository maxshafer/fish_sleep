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
if(!(args[1] %in% c("max_crep", "max_dinoc", "six_state", "four_state_max_crep", "four_state_max_dinoc"))) {  
  stop("first argument must be states in the model")
}

if(!(args[2] %in% c("cetaceans", "whippomorpha", "artiodactyla", "artiodactyla_minus_cetaceans", "ruminants", "mammals"))) {  
  stop("second argument must be the phylogenetic tree")
}

if(length(args) < 3) {  
  stop("must supply at least one model")
}  


# Section 2: Create trait data  --------

#first we'll import the trait data dataframe and subset the tree to include certain clades only

if(args[2] == "cetaceans"){
  trait.data <- read.csv(here("cetaceans_full.csv")) #from step 1 
}

if(args[2] == "artiodactyla"){
  trait.data <- read.csv(here("Sleepy_artiodactyla_full.csv")) #from step 2
}

if(args[2] == "artiodactyla_minus_cetaceans"){
  trait.data <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv")) #from step 2
}

if(args[2] == "ruminants"){
  trait.data <- read.csv(here("ruminants_full.csv")) #from step 2
}

if(args[2] == "whippomorpha"){
  trait.data <- read.csv(here("whippomorpha.csv")) #from step 1
}

if(args[2] == "mammals"){
  trait.data <- read.csv(here("sleepy_mammals.csv"))
  trait.data$Diel_Pattern <- trait.data$Bennie_diel
}

#take the column required based on the number of trait states
if(args[1] == "six_state"){
  trait.data <- trait.data[!is.na(trait.data$Diel_Pattern), c("tips", "Diel_Pattern")]
}

if(args[1] == "four_state_max_crep"){
  trait.data <- trait.data[!is.na(trait.data$max_crep), c("tips", "max_crep")]
}

if(args[1] == "four_state_max_dinoc"){
  trait.data <- trait.data[!is.na(trait.data$max_dinoc), c("tips", "max_dinoc")]
}

# Section 3: Subset the trees ---------------------------------------------

#for these models we will only use the max clade credibility tree from Cox et al
phylo_trees <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#subset trait data to only include species that are in the tree
trait.data <- trait.data[trait.data$tips %in% phylo_trees$tip.label,]

# this selects a tree that is only the subset with data (mutual exclusive)
phylo_trees <- keep.tip(phylo_trees, tip = trait.data$tips)

#determine if the models are considering a hidden rate

if("hidden_rate" %in% args){
  hidden_rate <- 2
} else {
  hidden_rate <- 1
}

# Section 4: Use returnCorModels to run corHMM models (ER, SYM, ARD, and/or bridge_only) on 1k possible trees --------

for(i in args[-(1:2)]){
  if(!(i %in% c("ER", "SYM", "ARD", "bridge_only", "hidden_rate", "CONSYM"))){
    stop("model not recognized")
  }
}

if("ER" %in% args){
  ER <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, model = "ER", node.states = "marginal")
}

if("SYM" %in% args){
  SYM <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, model = "SYM", node.states = "marginal")
  
}

#a constrained bridge only model based on a symmetrical model (does not allow diurnal-nocturnal transitions)
if("CONSYM" %in% args & "four_state_max_crep" %in% args){
  CONSYM <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = matrix(c(0,1,2,3,1,0,5,6,2,5,0,0,3,6,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), model = "SYM", node.states = "marginal")
}

if("CONSYM" %in% args & "four_state_max_dinoc" %in% args){
  CONSYM <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = matrix(c(0,1,2,3,1,0,5,6,2,5,0,0,3,6,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), model = "SYM", node.states = "marginal")
}

#fix the matrix for this
# if("CONSYM" %in% args & "six_state" %in% args){
#   CONSYM <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = matrix(c(), ncol = 6, nrow = 6, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), model = "SYM", node.states = "marginal")
# }

if("ARD" %in% args){
  ARD <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, model = "ARD", node.states = "marginal")
}

#a constrained model based on an all rates different model (does not allow diurnal-nocturnal transitions)
if("bridge_only" %in% args & "four_state_max_dinoc" %in% args){
  bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), node.states = "marginal")
}

if("bridge_only" %in% args & "four_state_max_crep" %in% args){
  bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), node.states = "marginal")
}

#fix the matrix for this
# if("bridge_only" %in% args & "six_state" %in% args){
#   bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"), c("(1, R1)", "(2, R1)", "(3,R1)", "(4, R1)"))), node.states = "marginal")
# }

# Section 5: Save the results out and extract likelihoods  --------
#use paste() to create a filename with all of the arguments
result_list <- lapply(args[-(1:2)], function(x) eval(as.name(x)))
names(result_list) <- paste(args[-(1:2)], "_model", sep = "")

saveRDS(result_list, paste(args[2], "fixed_max_clade_cred", args[1], "traits", paste0(args[-(1:2)], sep = "", collapse = "_"), "models", sep = "_"))
