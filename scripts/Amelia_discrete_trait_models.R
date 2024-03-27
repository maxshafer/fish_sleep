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
#open tree of life
library(rotl)

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


# Prescence of cortistatin ------------------------------------------------

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