#packages we will use
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
source("scripts/Amelia_functions.R")

# Section 1: Create trait data and subset the trees  --------


#we want to run select models on the 1k possible trees 
#allow us to better compare likelihoods for significant differences

#want to compare five models of the trinary cetacean data (cathemeral, diurnal, nocturnal)
#ER, SYM, ARD, cathemeral dead-end, cathemeral bridge

#first we'll import the trait data dataframe and subset the tree to include cetacean species only

trait.data3 <- cetaceans_full
trait.data3 <- trait.data3[!(is.na(trait.data3$Diel_Pattern_2)),]
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")

#we need to subset by species names in trait data (only species with behavioural data)
cetacean_trees <- lapply(mammal_trees, function(x) subsetTrees(tree = x, subset_names = trait.data3$tips))

#subset trait data to only include species that are in the tree
test_tree <- cetacean_trees[[3]]
trait.data3 <- trait.data3[trait.data3$tips %in% test_tree$tip.label,]

# Section 2: Use custom functions to re-run ace models on 1k possible  --------
cet_3state_ER <- lapply(cetacean_trees, function(x) returnAceModels(tree = x, trait.data = trait.data3, column = "Diel_Pattern_2", model = "ER"))
cet_3state_SYM <- lapply(cetacean_trees, function(x) returnAceModels(tree = x, trait.data = trait.data3, column = "Diel_Pattern_2", model = "SYM"))
cet_3state_ARD <- lapply(cetacean_trees, function(x) returnAceModels(tree = x, trait.data = trait.data3, column = "Diel_Pattern_2", model = "ARD"))


# Section 3: Save the results out and plot likelihoods  --------
result_list <- c(cet_3state_ER, cet_3state_SYM, cet_3state_ARD)
names(cet_ace_3state_results) <- c("cet_ace_3state_ER", "cet_ace_3state_SYM", "cet_ace_3tate_ARD")

saveRDS(result_list, "cet_ace_3state_results")

#get the likehoods for all model results 
cet_ER_likelihoods <- unlist(lapply(cet_3state_ER, function(x) returnLikelihoods(model = x)))
cet_SYM_likelihoods <- unlist(lapply(cet_3state_SYM, function(x) returnLikelihoods(model = x)))
cet_ARD_likelihoods <- unlist(lapply(cet_3state_ARD, function(x) returnLikelihoods(model = x)))

#organize this information into a dataframe
df1 <- as.data.frame(model = "ER", likelihoods = cet_ER_likelihoods)
df2 <- as.data.frame(model = "SYM", likelihoods = cet_SYM_likelihoods)
df3 <- as.data.frame(model = "ARD", likelihoods = cet_ARD_likelihoods)
df4 <- rbind(df1, df2)
df <- rbdind(df4, df3)

#plot
ggplot(df, aes(x = model, y = likelihoods)) + geom_point() + geom_boxplot()

#save out as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/ace_cet_3state_results.png", width=14,height=20,units="cm",res=1200)
ggplot(df, aes(x = model, y = likelihoods)) + geom_point() + geom_boxplot()
dev.off()
