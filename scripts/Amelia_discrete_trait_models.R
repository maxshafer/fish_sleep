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


# Prescence of cortistatin ------------------------------------------------

trait.data <- read.csv(here("cetaceans_full.csv"))

#create a phylogeny of just the cetaceans in the cortistatin paper
trait.data <- trait.data[trait.data$Species_name %in% c("Tursiops truncatus", "Tursiops aduncus", "Sousa chinesis", "Globicephala melas", "Peponocephala electra", "Lagenorhynchus obliquidens","Orcinus orca", 
                                                                   "Neophocaena asiaeorientalis", "Phocoena phocoena", "Phocoena sinus", "Delphinapterus leucas", "Monodon monoceros", "Inia geoffrensis", "Pontoporia blainvillei", "Lipotes vexillifer",
                                                                   "Platanista gangetica", "Mesoplodon bidens", "Ziphius cavirostris", "Kogia breviceps", "Physeter catodon", "Balaenoptera bonaerensis", "Balaenoptera acutorostrata", "Balaenoptera musculus",
                                                                   "Balaenoptera edeni", "Balaenoptera physalus", "Megaptera novaeangliae", "Eschrichtius robustus", "Eubalaena japonica", "Eubalaena glacialis"), ]

#add a column containing whether they have intact cort (missing Platanista gangetica and Sousa chinesis from the tree)
#just the type of mutation, should add info about where the mutation is
#do each of these have the same effect on the phenotype?
#cort_cetaceans$cort_144bp <- c("NA", "") need to think of how to do this
trait.data$cort_219bp <- c("1D", "1D", "NA", "NA", "NA", "1D", "NA", "NA", "NA", "NA", "1D", "1I", "NA", "PS", "1I", "NA", "NA", "NA", "1I", "NA", "PS", "NA", "NA", "NA", "NA", "NA", "NA")

#full list of species and cort status 
#cort_list<- c("Tursiops truncatus", "Tursiops aduncus", "Sousa chinesis", "Globicephala melas", "Lagenorhynchus obliquidens","Orcinus orca", "Neophocaena asiaeorientalis", "Phocoena phocoena", "Phocoena sinus", "Delphinapterus leucas", "Monodon monoceros", "Inia geoffrensis", "Pontoporia blainvillei", "Lipotes vexillifer", "Mesoplodon bidens", "Ziphius cavirostris", "Kogia breviceps", "Physeter catodon", "Balaenoptera bonaerensis", "Balaenoptera acutorostrata", "Balaenoptera musculus", "Balaenoptera edeni", "Balaenoptera physalus", "Megaptera novaeangliae", "Eschrichtius robustus", "Eubalaena japonica", "Eubalaena glacialis")
#cort_219bp_list <- c("1D", "1D", "1D", "NA", "NA", "NA", "1D", "NA", "NA", "NA", "NA", "1D", "1I", "NA", "PS", "1I", "NA", "NA", "NA", "1I", "NA", "PS", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

# cort_tr <- tol_induced_subtree(ott_ids = cort_cetaceans$ott_id[cort_cetaceans$flags %in% c("sibling_higher", "")], label_format = "id")
# cort_tr$tip.label <- cort_cetaceans$tips[match(cort_tr$tip.label, paste("ott", cort_cetaceans$ott_id, sep = ""))]
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), ]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)


#custom.colours3 <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
cort_diel_tree <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "cort_219bp")]
cort_diel_tree <- cort_diel_tree + geom_tile(data = cort_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_2),  inherit.aes = FALSE, color = "transparent")# + scale_fill_manual(values = custom.colours)
cort_diel_tree <- cort_diel_tree + geom_tiplab(color = "black", size = 2.5, offset = 3)
cort_diel_tree <- cort_diel_tree + geom_tile(data = cort_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x = x + 1.5, fill = cort_219bp), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = c('steelblue1', 'lightskyblue', "#dd8ae7", "#FC8D62", "grey","#66C2A5", "#A6D854", "royalblue"))
cort_diel_tree

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/cortistatin_cetaceans.png", width=18,height=15,units="cm",res=1200)
cort_diel_tree
dev.off()


# # Prescence of echolocation ---------------------------------------------
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]

echo <- read_xlsx("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_discrete_traits/Coombs_et_al_2022.xlsx")
echo <- echo[!(is.na(echo$Echo)), c("Museum ID", "Echo", "Diet", "Dentition", "FM", "Habitat")]
echo$tips <- str_replace(echo$`Museum ID`, " ", "_")

echo <- echo %>% filter(echo$tips %in% trait.data$tips)

trait.data <- merge(trait.data, echo, by = "tips")

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#should fix colours
custom.colours <- c("#dd8ae7", "pink", "black", "grey80", "#66C2A5", "#A6D854")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Echo")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Echo), inherit.aes = FALSE, colour = "transparent") 
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/echo_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#do the same for other discrete traits, like habitat
custom.colours <- c("#dd8ae7", "black", "grey70", "pink", "#66C2A5", "#A6D854", "grey30", "grey90")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Habitat")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Habitat), inherit.aes = FALSE, colour = "transparent") 
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/habitat_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#feeding mechanism
custom.colours <- c( "black","#dd8ae7","pink", "grey40",  "#66C2A5", "#A6D854", "grey80")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "FM")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = FM), inherit.aes = FALSE, colour = "transparent")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/feeding_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#feeding mechanism
custom.colours <- c( "black","#dd8ae7","pink", "grey40",  "#66C2A5", "#A6D854", "grey80")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "FM")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = FM), inherit.aes = FALSE, colour = "transparent")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/feeding_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#especially want to look at body size (see if it associates with cathemerality)

# Cetacean discrete traits  -----------------------------------------------
###other cetacean dataframes of discrete traits to model with
discrete_cet <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")

#dive depth
#body size


behav <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Manger_2013_discrete_cet.csv")

#body mass
#encephalization quotient 
#group size
#social dynamics
#longevity