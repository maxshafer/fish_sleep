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

#args <- c("six_state", "artiodactyla", "ARD", "hidden_rate", "bridge_only")
# Section 1: Arguments ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(!(args[1] %in% c("max_crep", "max_dinoc", "six_state", "four_state_max_crep", "four_state_max_dinoc"))) {  
  stop("first argument must be states in the model")
}

if(!(args[2] %in% c("cetaceans", "artiodactyla", "artiodactyla_minus_cetaceans"))) {  
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

if(args[2] == "artiodactyla"){
  trait.data <- read.csv(here("Sleepy_artiodactyla_full.csv"))
}

if(args[2] == "artiodactyla_minus_cetaceans"){
  trait.data <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv"))
}

#the cetacean and artiodactyla dataframes are formatted to have Diel_Pattern_3 be max_crep (di/crep and noc/crep classified as crep, cath and crep pooled together)
#and Diel_Pattern_4 to be max_dinoc (di/crep classified as di, noc/crep as noc, cathemeral alone as third state)

if(args[1] == "max_crep"){
  trait.data <- trait.data[,c("tips", "Diel_Pattern_3")]
  
}
if(args[1] == "max_dinoc"){
  trait.data <- trait.data[,c("tips", "Diel_Pattern_4")]
}

if(args[1] == "six_state"){
  trait.data <- trait.data[, c("tips", "Diel_Pattern_2")]
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

# determine which custom rate matrix will be needed
#there are 3 options we would need different custom rate matrices for
#artiodactyla six state bridge only
if("six_state" %in% args & "artiodactyla" %in% args){
  #don't allow transitions from noc -> di, noc -> di/crep, noc/crep -> di, noc/crep -> di/crep, di -> noc, di -> noc/crep, di/crep -> noc, di/crep -> noc/crep 
  generic_ratemat <- getStateMat4Dat(trait.data)$rate.mat
  #need to disallow transition 14 (Noc -> Di), 15 (Noc/crep -> Di), 19 (Noc -> Di/crep), 20 (Noc/crep -> Di/crep), 23 (Di -> Noc), 24 (Di/crep -> Noc), 28 (Di -> Noc/crep), 29 (Di/crep -> Noc/crep)
  custom_rate_matrix <- dropStateMatPars(generic_ratemat, c(14, 15, 19, 20, 23, 24, 28, 29))
}

#cetacean or artiodactyla maxcrep or maxdinoc bridge only
if(!("six_state" %in% args)){
  #don't allow transitions from noc -> di and di -> noc
  generic_ratemat <- getStateMat4Dat(trait.data)$rate.mat
  #need to disallow transition 4 (noc -> di), and 6 (di -> noc)
  custom_rate_matrix <- dropStateMatPars(generic_ratemat, c(4,6))
}

#cetacean six state bridge only
# this is actually a four state since there are no cathemeral crepuscular or diurnal crepuscular species
if("six_state" %in% args & "cetaceans" %in% args){
  #don't allow transitions from noc -> di, noc -> di/crep, noc/crep -> di, di -> noc, di -> noc/crep, di/crep -> noc 
  generic_ratemat <- getStateMat4Dat(trait.data)$rate.mat
  #need to disallow transition 5 (Noc -> Di), 6 (Noc/crep -> Di),8 (Di -> Noc), 11 (Di -> Noc/crep)
  custom_rate_matrix <- dropStateMatPars(generic_ratemat, c(5,6,8,11))
}


# Section 4: Use returnCorModels to run corHMM models (ER, SYM, ARD, and/or bridge_only) on 1k possible trees --------

for(i in args[-(1:2)]){
  if(!(i %in% c("ER", "SYM", "ARD", "bridge_only", "hidden_rate"))){
    stop("model not recognized")
  }
}

if("ER" %in% args){
  ER <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, model = "ER", node.states = "marginal")
}

if("SYM" %in% args){
  SYM <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, model = "SYM", node.states = "marginal")

}

if("ARD" %in% args){
  ARD <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, model = "ARD", node.states = "marginal")

}

if("bridge_only" %in% args){
  bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = hidden_rate, rate.mat = custom_rate_matrix, node.states = "marginal")
}

# Section 5: Save the results out and extract likelihoods  --------
#use paste() to create a filename with all of the arguments
result_list <- lapply(args[-(1:2)], function(x) eval(as.name(x)))
names(result_list) <- paste(args[-(1:2)], "_model", sep = "")

saveRDS(result_list, paste(args[2], "max_clade_cred", args[1], "traits", paste0(args[-(1:2)], sep = "", collapse = "_"), "models", sep = "_"))
