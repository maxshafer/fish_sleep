# Section 0: Packages -----------------------------------------------------
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

trait.data3 <- cetaceans_full[,c("Diel_Pattern_2", "tips")]
trait.data3$Diel_Pattern_2 <- tolower(trait.data3$Diel_Pattern_2)
trait.data3 <- trait.data3[!(is.na(trait.data3$Diel_Pattern_2)),]
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")

#we need to subset by species names in trait data (only species with behavioural data)
cetacean_trees <- lapply(mammal_trees, function(x) subsetTrees(tree = x, subset_names = trait.data3$tips))

#subset trait data to only include species that are in the tree
test_tree <- cetacean_trees[[3]]
trait.data3 <- trait.data3[trait.data3$tips %in% test_tree$tip.label,]

#subset cetacean trees to 50 trees for now to see if it runs
cetacean_trees <- cetacean_trees[1:10]

# Section 2: Use custom functions to re-run ace models (ER, SYM, ARD) on 1k possible trees --------
cet_3state_ER <- lapply(cetacean_trees, function(x) returnCorModels(tree = x, trait.data = trait.data3, diel_col = "Diel_Pattern_2", rate.cat = 1, model = "ER", node.states = "marginal"))
cet_3state_SYM <- lapply(cetacean_trees, function(x) returnCorModels(tree = x, trait.data = trait.data3, diel_col = "Diel_Pattern_2", rate.cat = 1,model = "SYM", node.states = "marginal"))
cet_3state_ARD <- lapply(cetacean_trees, function(x) returnCorModels(tree = x, trait.data = trait.data3, diel_col = "Diel_Pattern_2", rate.cat = 1, model = "ARD", node.states = "marginal"))

#Section 2.5: Cathemeral dead-end and cathemeral bridge model
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
#this is a generic all rates different matrix that we can now edit
test_ratemat <- generic_ratemat$rate.mat
#to drop transitions from cathemerality to diurnality (1->2) and from cathemerality to nocturnality (1 -> 3)
#this drops 3 and 5 from the matrix 
deadend_matrix <- dropStateMatPars(test_ratemat, c(3, 5))
#easier way to do this below but isn't working
deadend_matrix2 <- matrix(c(0,1,2,0,0,3,0,4,0), ncol = 3, nrow = 3, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)"), c("(1, R1)", "(2, R1)", "(3,R1)")))
#cet_3state_deadend <- lapply(cetacean_trees, function(x) returnCustomCorModels(tree = x, trait.data = trait.data3, diel_col = "Diel_Pattern_2", rate.cat = 1, rate.mat = deadend_matrix, model = "ARD"))

bridge_matrix <- matrix(c(0,1,2,3,0,0,4,0,0), ncol = 3, nrow = 3, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)"), c("(1, R1)", "(2, R1)", "(3,R1)")))
cet_3state_bridge <- lapply(cetacean_trees, function(x) returnCustomCorModels(tree = x, trait.data = trait.data3, diel_col = "Diel_Pattern_2", rate.cat = 1, rate.mat = bridge_matrix, model = "ARD", node.states = "marginal"))


# Section 3: Save the results out and plot likelihoods  --------
result_list <- c(cet_3state_ER, cet_3state_SYM, cet_3state_ARD, cet_3state_bridge)
names(result_list) <- c("cet_ace_3state_ER", "cet_ace_3state_SYM", "cet_ace_3tate_ARD", "cet_3state_bridge")

saveRDS(result_list, "cet_ace_3state_results_10trees")

#get the likehoods for all model results 
cet_ER_likelihoods <- unlist(lapply(cet_3state_ER, function(x) returnLikelihoods(model = x)))
cet_SYM_likelihoods <- unlist(lapply(cet_3state_SYM, function(x) returnLikelihoods(model = x)))
cet_ARD_likelihoods <- unlist(lapply(cet_3state_ARD, function(x) returnLikelihoods(model = x)))
#cet_deadend_likelihoods <- unlist(lapply(cet_3state_deadend, function(x) returnLikelihoods(model =x)))
cet_bridge_likelihoods <- unlist(lapply(cet_3state_bridge, function(x) returnLikelihoods(model = x)))

#organize this information into a dataframe (use pivot_wider instead)
df1 <- data.frame(model = "ER", likelihoods = cet_ER_likelihoods)
df2 <- data.frame(model = "SYM", likelihoods = cet_SYM_likelihoods)
df3 <- data.frame(model = "ARD", likelihoods = cet_ARD_likelihoods)
df4 <- rbind(df1, df2)
df <- rbind(df4, df3)
#df <- rbind(df, data.frame(model = "cathemeral dead-end", likelihoods = cet_deadend_likelihoods))
df <- rbind(df, data.frame(model = "cathemeral bridge", likelihoods = cet_bridge_likelihoods))

#plot
ggplot(df, aes(x = model, y = likelihoods)) + geom_point() + geom_boxplot() + geom_jitter()

#save out as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_cet_3state_results.png", width=14,height=20,units="cm",res=1200)
ggplot(df, aes(x = model, y = likelihoods)) + geom_point() + geom_boxplot()
dev.off()

#to do: use tictoc to time how long it takes to run these on the full 1k trees https://www.jumpingrivers.com/blog/timing-in-r/
#look at the distribution of cetacean trees using ggdensitree https://yulab-smu.top/treedata-book/chapter4.html
#repeat with all of artiodactyla
#plot the distribution of the transition rates from these models similar to https://www.sciencedirect.com/science/article/pii/S0960982219301563#mmc1