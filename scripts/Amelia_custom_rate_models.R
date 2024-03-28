# Section 0: Import and run packages -------------------------------------------------
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

# Section 1.0 Cetacean cathemeral dead end-------------------------------

#create a custom matrix to test if cathemerality is a dead end (start with disallowing transitions back to diurnality)
# #ie rate from cathemeral into other diel patterns is zero

generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat

#to drop transitions from cathemerality to diurnality (1->2)
#this drops 3 from the matrix 
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(3))
custom_rate_matrix
 
cor_deadend_di <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

likelihoods <- rbind(likelihoods, c("cor_deadend_di", cor_deadend_di$loglik, "three states (di, noc, cath), ARD do not allow cath to di, cetacea"))

#try with just cathemeral to nocturnal
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat

#to drop transitions from cathemerality to nocturnality (1 -> 3)
#this drops 5 from the matrix 
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(5))
custom_rate_matrix

cor_deadend_noc <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

likelihoods <- rbind(likelihoods, c("cor_deadend_noc", cor_deadend_noc$loglik, "three states (di, noc, cath), ARD do not allow cath to noc, cetacea"))

#create a custom matrix to test if cathemerality is a dead end 
# #ie rate from cathemeral into other diel patterns is zero

generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat

#to drop transitions from cathemerality to diurnality (1->2) and from cathemerality to nocturnality (1 -> 3)
#this drops 3 and 5 from the matrix 
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(3, 5))
custom_rate_matrix


custom_rate_matrix2 <- matrix(c(0,1,2,0,0,3,0,4,0), ncol = 3, nrow = 3)
custom_rate_matrix2

cor_deadend <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix, model = "ARD", node.states = "marginal")
plotMKmodel(cor_deadend)

cor_dead_test <- corHMM(phy=cetacean_trees[[5]], data = trait.data3, rate.cat = 1, model = "ARD")

cor_deadend2 <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix2, model = "ARD")
plotMKmodel(cor_deadend2)

likelihoods <- rbind(likelihoods, c("cor_deadend", cor_deadend$loglik, "three states (di, noc, cath), ARD do not allow any transitions out of cath, cetacea"))

# Section X.0 Cathemeral dead-end artiodactyla ----------------------------
#can run total dead-end first (doesn't allow transitions out of cathemerality at all)
generic_ratemat <- getStateMat4Dat(cor_maor_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
test_ratemat
#to drop transitions from cathemerality to diurnality (1->2) and from cathemerality to nocturnality (1 -> 3)
#this drops 3 and 5 from the matrix 
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(3, 5))
custom_rate_matrix
cor_deadend_artio <- corHMM(phy = cor_maor_ard3$phy, data = cor_maor_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)
artio_likelihoods <- rbind(artio_likelihoods, c("cor_deadend_artio", cor_deadend_artio$loglik, "three states (di, noc, cath), ARD do not allow transitions out of cath, artio with cetcea"))


# Section X.0 Cathemeral dead-end artiodactyla without cetacea ------------
generic_ratemat <- getStateMat4Dat(cor_maorja_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
test_ratemat
#to drop transitions from cathemerality to diurnality (1->2) and from cathemerality to nocturnality (1 -> 3)
#this drops 3 and 5 from the matrix 
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(3, 5))
custom_rate_matrix
cor_deadend_artio <- corHMM(phy = cor_maorja_ard3$phy, data = cor_maorja_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)
artio_likelihoods <- rbind(artio_likelihoods, c("cor_deadend_ja", cor_deadend_ja$loglik, "three states (di, noc, cath), ARD do not allow transitions out of cath, artio without cetcea"))


# Section 2.0 Cetacean cathemeral bridge ----------------------------------
#don't allow transitions from nocturnality to diurnality
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat
#need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> nocturnal/3 (6)
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(4, 6))
custom_rate_matrix

cor_cath_bridge <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

likelihoods <- rbind(likelihoods, c("cor_cath_bridge", cor_cath_bridge$loglik, "three states (di, noc, cath), ARD, does not allow noc <-> di, cetacea"))

# Section 3.0 Artiodactyla cathemeral dead-end Custom rates with artiodactyla ------------------------------------------
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat
#need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> nocturnal/3 (6)
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(4, 6))
custom_rate_matrix

cor_cath_bridge <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

artio_likelihoods <- rbind(artio_likelihoods, c("cor_cath_bridge", cor_cath_bridge$loglik, "three states (di, noc, cath), ARD, does not allow noc <-> di, cetacea"))

# Section X.0 Cathemeral bridge, artio with cetaceans ----------------------------------
#don't allow transitions between nocturnality and diurnality
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat
#need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> noctunral/3 (6)
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(4, 6))
custom_rate_matrix


cor_cath_bridge_artio <- corHMM(phy = cor_maor_ard3$phy, data = cor_maor_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

artio_likelihoods <- rbind(artio_likelihoods, c("cor_cath_bridge_artio", cor_cath_bridge_artio$loglik, "three states (di, noc, cath), ARD, does not allow noc <-> di, artiofactyla with cetacea"))

# Section X.0 Cathemeral bridge, artio without cetaceans
#don't allow transitions between nocturnality and diurnality
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat
#need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> nocturnal/3 (6)
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(4, 6))
custom_rate_matrix

cor_cath_bridge_ja <- corHMM(phy = cor_maorja_ard3$phy, data = cor_maorja_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

artio_likelihoods <- rbind(artio_likelihoods, c("cor_cath_bridge_ja", cor_cath_bridge_ja$loglik, "three states (di, noc, cath), ARD, does not allow noc <-> di, artio without cetacea"))

# Cetacean discrete traits  -----------------------------------------------
###discrete traits to model with
discrete_cet <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")
more_cet_traits <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Coombs_et_al_2022.xlsx")
more_cet_traits <- more_cet_traits[!(is.na(more_cet_traits$No.)),]


# Section X.0 Hidden rates models, cetacea, maybe move this to another section --------
cor_model_er1_2r <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 2, model = "ER", node.states = "marginal")
cor_model_sym1_2r <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 2, model = "SYM", node.states = "marginal")
cor_model_ard1_2r <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 2, model = "ARD", node.states = "marginal")

likelihoods <- rbind(likelihoods, c("cor_model_er1_2r", cor_model_er1_2r$loglik, "diurnal, nocturnal, one hidden rate, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym1_2r", cor_model_sym1_2r$loglik, "diurnal, nocturnal, one hidden rate, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_ard1_2r", cor_model_ard1_2r$loglik, "diurnal, nocturnal, one hidden rate, cetacea"))

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_er_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_er1_2r)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_sym_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_sym1_2r)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_ard_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_ard1_2r)
dev.off()


# Section X.0 Hidden rates models, artiodactyla with cetacea, maybe move this to another section --------
cor_maor_er1_2r <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 2, model = "ER", node.states = "marginal")
cor_maor_sym1_2r <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 2, model = "SYM", node.states = "marginal")
cor_maor_ard1_2r <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 2, model = "ARD", node.states = "marginal")

artio_likelihoods <- rbind(artio_likelihoods, c("cor_model_er1_2r", cor_model_er1_2r$loglik, "diurnal, nocturnal, one hidden rate, cetacea"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_model_sym1_2r", cor_model_sym1_2r$loglik, "diurnal, nocturnal, one hidden rate, cetacea"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_model_ard1_2r", cor_model_ard1_2r$loglik, "diurnal, nocturnal, one hidden rate, cetacea"))

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_er_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_er1_2r)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_sym_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_sym1_2r)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_ard_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_ard1_2r)
dev.off()


# Cetaceans with crepuscularity as a secondary trait ----------------------

#We want to run some models on dataset of cetaceans with six categories of diel patterns
# Di-crep, Di-not-crep, Noc-crep, Noc-not-crep, Cath-crep, Cath-not-crep
# diel pattern 2 in cetaceans full already has all cetacean species categorized this way except for cath-crep (ie true crepuscular)

# We only need the tips (species names) and diel pattern 2. We start with 80 species with data
trait.data.crep <- cetaceans_full[, c(6,11)]
#change a species to be true crepuscular to see what it looks like
#trait.data.crep[1,1] <- "Cathemeral/crepuscular"
# selects only data that is in the mammal tree. 72 species are in the tree
trait.data.crep <- trait.data.crep[trait.data.crep$tips %in% mam.tree$tip.label,]
row.names(trait.data.crep) <- trait.data.crep$tips
# this subsets the tree to only include the species which we have data on (72 species)
trpy_n_crep <- keep.tip(mam.tree, tip = trait.data.crep$tips)

ggtree(trpy_n_crep) + geom_tiplab()

#first run the ace models
trait.vector.crep <- trait.data.crep$Diel_Pattern_2
ace_model_er_crep <- ace(trait.vector.crep, trpy_n_crep, model = "ER", type = "discrete")
ace_model_sym_crep <- ace(trait.vector.crep, trpy_n_crep, model = "SYM", type = "discrete")
ace_model_ard_crep <- ace(trait.vector.crep, trpy_n_crep, model = "ARD", type = "discrete")

#run the cor models
cor_model_er_crep <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym_crep <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard_crep <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

#add likelihoods and plot

#repeat but with a hidden rate model

#run the cor models
cor_model_er_crep_2r <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 2, model = "ER", node.states = "marginal")
cor_model_sym_crep_2r <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 2, model = "SYM", node.states = "marginal")
cor_model_ard_crep_2r <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 2, model = "ARD", node.states = "marginal")

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_ard_crep_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_ard_crep_2r)
dev.off()