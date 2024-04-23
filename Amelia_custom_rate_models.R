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

#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("cetaceans_full.csv"))

#run cor_maor_ard and cor_model ard
#these are the basis for the artiodactyla models and cetacean models respectively


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

# Section 3.0 Artiodactyla cathemeral dead-end with artiodactyla ------------------------------------------
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


# In progress: cetacean dead-end 1k trees ----------------------

# #Section 2.5 from R script: Cathemeral dead-end model
# generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
# #this is a generic all rates different matrix that we can now edit
# test_ratemat <- generic_ratemat$rate.mat
# #to drop transitions from cathemerality to diurnality (1->2) and from cathemerality to nocturnality (1 -> 3)
# #this drops 3 and 5 from the matrix 
# deadend_matrix <- dropStateMatPars(test_ratemat, c(3, 5))
# #easier way to do this below but isn't working
# deadend_matrix2 <- matrix(c(0,1,2,0,0,3,0,4,0), ncol = 3, nrow = 3, dimnames = list(c("(1, R1)", "(2, R1)", "(3,R1)"), c("(1, R1)", "(2, R1)", "(3,R1)")))
# #cet_3state_deadend <- lapply(cetacean_trees, function(x) returnCustomCorModels(tree = x, trait.data = trait.data, diel_col = "Diel_Pattern_2", rate.cat = 1, rate.mat = deadend_matrix, model = "ARD"))

#In progress: cetacean cortistatin association

# In progress: Prescence of cortistatin ------------------------------------------------

#create a phylogeny of just the cetaceans in the cortistatin paper
cort_cetaceans <- resolved_names[resolved_names$unique_name %in% c("Tursiops truncatus", "Tursiops aduncus", "Sousa chinesis", "Globicephala melas", "Peponocephala electra", "Lagenorhynchus obliquidens","Orcinus orca", 
                                                                   "Neophocaena asiaeorientalis", "Phocoena phocoena", "Phocoena sinus", "Delphinapterus leucas", "Monodon monoceros", "Inia geoffrensis", "Pontoporia blainvillei", "Lipotes vexillifer",
                                                                   "Platanista gangetica", "Mesoplodon bidens", "Ziphius cavirostris", "Kogia breviceps", "Physeter catodon", "Balaenoptera bonaerensis", "Balaenoptera acutorostrata", "Balaenoptera musculus",
                                                                   "Balaenoptera edeni", "Balaenoptera physalus", "Megaptera novaeangliae", "Eschrichtius robustus", "Eubalaena japonica", "Eubalaena glacialis"), ]
View(cort_cetaceans)

#add a column containing whether they have intact cort (missing Platanista gangetica and Sousa chinesis from the tree)
#just the type of mutation, should add info about where the mutation is
#do each of these have the same effect on the phenotype?
#cort_cetaceans$cort_144bp <- c("NA", "") need to think of how to do this
cort_cetaceans$cort_219bp <- c("1D", "1D", "NA", "NA", "NA", "1D", "NA", "NA", "NA", "NA", "1D", "1I", "NA", "PS", "1I", "NA", "NA", "NA", "1I", "NA", "PS", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

#full list of species and cort status 
#cort_list<- c("Tursiops truncatus", "Tursiops aduncus", "Sousa chinesis", "Globicephala melas", "Lagenorhynchus obliquidens","Orcinus orca", "Neophocaena asiaeorientalis", "Phocoena phocoena", "Phocoena sinus", "Delphinapterus leucas", "Monodon monoceros", "Inia geoffrensis", "Pontoporia blainvillei", "Lipotes vexillifer", "Mesoplodon bidens", "Ziphius cavirostris", "Kogia breviceps", "Physeter catodon", "Balaenoptera bonaerensis", "Balaenoptera acutorostrata", "Balaenoptera musculus", "Balaenoptera edeni", "Balaenoptera physalus", "Megaptera novaeangliae", "Eschrichtius robustus", "Eubalaena japonica", "Eubalaena glacialis")
#cort_219bp_list <- c("1D", "1D", "1D", "NA", "NA", "NA", "1D", "NA", "NA", "NA", "NA", "1D", "1I", "NA", "PS", "1I", "NA", "NA", "NA", "1I", "NA", "PS", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

cort_tr <- tol_induced_subtree(ott_ids = cort_cetaceans$ott_id[cort_cetaceans$flags %in% c("sibling_higher", "")], label_format = "id")
cort_tr$tip.label <- cort_cetaceans$tips[match(cort_tr$tip.label, paste("ott", cort_cetaceans$ott_id, sep = ""))]
cort_diel_tree <- ggtree(cort_tr, layout = "circular") %<+% cort_cetaceans[,c("tips", "diel", "cort_219bp")]
custom.colours3 <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
cort_diel_tree <- cort_diel_tree + geom_tile(data = cort_diel_tree$data[1:length(cort_tr$tip.label),], aes(y=y, x=x, fill = diel),  inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours3)
cort_diel_tree <- cort_diel_tree + geom_tiplab(color = "black", size = 2.5, hjust = -0.5)
cort_diel_tree

cort_diel_tree2 <- cort_diel_tree + geom_tile(data = cort_diel_tree$data[1:length(cort_tr$tip.label),], aes(y=y, x = x + 1, fill = cort_219bp), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = c('red', 'blue', "#d5cab2", "#dd8ae7", "#FC8D62", "black","#66C2A5", "green"))
cort_diel_tree2 

png("C:/Users/ameli/OneDrive/Documents/R_projects/cortistatin_cetaceans.png", width=18,height=15,units="cm",res=1200)
print(cort_diel_tree)
dev.off()


# In progress: other discrete traits --------------------------------------

#Load in the discrete traits 

##Data from Manger et al  2013, on cetacean body mass, brain mass, and social traits
#MBM – male body mass (g); FBM – female body mass (g); ABM – average body mass (g); BrM – brain mass (g); EQ – encephalization quotient; L – longevity (days); SMF – age at sexual maturity of females (days) 
#AGS – average group size in number of individuals; RGS – range of group size in number of individuals; GSD – group social dynamics; FS – feeding strategy. 
#Group Social Dynamics: Sol – mostly solitary lifestyle; SSG – small stable groups of less than 10 individuals; FF? – possible fission–fusion social dynamic; SG - stable groups of more than 10 individuals; FF – fission–fusion social dynamic. 
#Feeding strategies: SF – skim feeder; OCO – shows occasional co-operation in feeding; LF – lunge feeder; SwF – swallow feeder; CO – shows regular co-operation in feeding; SuF – suction feeder; R – raptorial feeder.
#Data for this table derived from the following sources: Nowak (1999), Lefebvre et al. (2006a), Manger (2006), Best et al. (2009), Perrin et al. (2009)

behav <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Manger_2013_discrete_cet.csv")

# Cetacean discrete traits  -----------------------------------------------
###other cetacean dataframes of discrete traits to model with
discrete_cet <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")
more_cet_traits <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Coombs_et_al_2022.xlsx")
more_cet_traits <- more_cet_traits[!(is.na(more_cet_traits$No.)),]


