# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(ggplot2)
library(RColorBrewer)

setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Load in Cox mammal tree + Amelia primary source data -----------------------

#need to use the mamm tree for modelling and ancestral trait reconstruction
#trees to make: trinary dataset -maxdinoc, trinary dataset -maxcrep, 5state

# ## Read in the mammalian phylogeny
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia primary source cetacean data
trait.data <- read.csv(here("cetaceans_full.csv"))
clade_name <- "cetcaeans"

#use below for primary source whippomorpha data
#trait.data <- read.csv(here("whippomorpha_full.csv"))
#clade_name <- "whippomorpha"

#use below for all artiodactyla
# trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
# clade_name <- "artiodactyla"

#use below for artiodactyla minus cetaceans
#trait.data <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv"))
#clade_name <- "artiodactyla_minus_cetaceans"

#trait.data <- read.csv(here("ruminants_full.csv"))
#clade_name <- "ruminants"

# trait.data <- read.csv(here("sleepy_mammals.csv"))
# clade_name <- "mammals"

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#does not include species with no activity pattern data (they are filtered out when processing the sleepy cetacean dataframe)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/all_cetacean_tree.png", width=40,height=40, units="cm",res=100)
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 1)
dev.off()

# Section 2: Plot all state trait data tree --------------------------------------

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

#plot the tree
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 1)
diel.plot.all 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", clade_name, "diel_plot_all_state.png"), width=24,height=18,units="cm",res=1200)
print(diel.plot.all)
dev.off()

#plot tree with confidence data
#trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2", "Confidence")]
trait.data.all <- trait.data[, c("tips", "Diel_Pattern_2", "Confidence")]
#change any sub 1 confidence to confidence 1
trait.data.all$Confidence[trait.data.all$Confidence == 0] <- 1
trait.data.all$Confidence[trait.data.all$Confidence == 0.5] <- 1
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

custom.colours.all <- c("grey90", "grey", "grey50", "grey30", "black","white", "#dd8ae7","pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "blue", "brown")
#custom.colours.all <- c("red", "red", "brown3", "grey50", "grey30", "black","#dd8ae7","pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

#plot the tree, separate geoms for diel pattern and confidence
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2", "Confidence")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "black", width = 2.5)
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x + 2.5, y=y, fill = as.factor(Confidence)), inherit.aes = FALSE, colour = "black", width = 2.5) + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all #+ geom_tiplab(size = 2, offset = 3)
diel.plot.all 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", "final", clade_name, "_confidence_diel_plot_all_state.png"), width=20,height=14,units="cm",res=1200)
print(diel.plot.all)
dev.off() 

# Section 3: Plot 3-state maxdinoc tree --------------------------------------

#creating a tree with three trait states (cathemeral, diurnal, nocturnal)

#remove species with no data
trait.data.maxDN <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trait.data.maxDN$Diel_Pattern_2 <- str_replace_all(trait.data.maxDN$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data.maxDN$Diel_Pattern_2 <- str_replace_all(trait.data.maxDN$Diel_Pattern_2, "diurnal/crepuscular", "nocturnal")

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxDN$tips)

#plot the tree
custom.colours.MDN <- c("#dd8ae7", "#FC8D62", "#66C2A5")
diel.plot.MDN <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxDN[,c("tips", "Diel_Pattern_2")]
diel.plot.MDN <- diel.plot.MDN + geom_tile(data = diel.plot.MDN$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MDN)
diel.plot.MDN <- diel.plot.MDN + geom_tiplab(size = 2)
diel.plot.MDN 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", clade_name, "diel_plot_max_dinoc.png"), width=24,height=18,units="cm",res=1200)
print(diel.plot.MDN)
dev.off()

# Section 4: Plot 3-state maxcrep tree --------------------------------------

#creating a tree with three trait states (crepuscular, diurnal, nocturnal)

#remove species with no data
trait.data.maxcrep <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "nocturnal/crepuscular", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "diurnal/crepuscular", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "cathemeral", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "crep/cath", "crepuscular/cathemeral")

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxcrep$tips)

#plot the tree
custom.colours.MC <- c("pink", "#FC8D62", "#66C2A5")
diel.plot.MC <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxcrep[,c("tips", "Diel_Pattern_2")]
diel.plot.MC <- diel.plot.MC + geom_tile(data = diel.plot.MC$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MC)
diel.plot.MC <- diel.plot.MC + geom_tiplab(size = 2)
diel.plot.MC 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", clade_name, "diel_plot_max_crep.png"), width=24,height=18,units="cm",res=1200)
print(diel.plot.MC)
dev.off()


# Section 5: Plot the ambiguity in the cetacean tree ----------------------

#look at the distribution of cetacean trees using ggdensitree https://yulab-smu.top/treedata-book/chapter4.html

#change this so we save out the result and load it in
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
cetacean_trees <- lapply(mammal_trees, function(x) subsetTrees(tree = x, subset_names = trait.data.all$tips))

#sample 100 random trees from the total 1k
#onehundred_trees <- cetacean_trees[sample(1:length(cetacean_trees), 100)]

#plot these using ggdensitree
# density_tree100 <- ggdensitree(onehundred_trees, alpha = 0.3, colour = "yellowgreen")
# density_tree100
# 
# #save out
# png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/density_tree100.png", width=24,height=18,units="cm",res=1200)
# print(density_tree100)
# dev.off()

#repeat the same process but with 50 trees
twenty_trees <- cetacean_trees[sample(1:length(cetacean_trees), 1000)]
ggdensitree(twenty_trees, alpha = 0.3, colour = "yellowgreen") + geom_tiplab(size=3, color = "black") 


png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/density_treebig.png", width=100,height=18,units="cm",res=1000)
density_20
dev.off()



# Section 6: Artiodactyla suborders ---------------------------------------

# ## Read in the mammalian phylogeny
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia data
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#filter for the suborder you want
#trait.data <- trait.data %>% filter(Order == "Ruminantia")
#trait.data <- trait.data %>% filter(Order == "Suina")
#trait.data <- trait.data %>% filter(Order == "Tylopoda")
trait.data <- trait.data %>% filter(Order == "Whippomorpha")

trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
#custom.colours.all <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "grey")
custom.colours.all <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "grey")
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x+1.1, y=y, fill = Diel_Pattern_2), width = 2, inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all, name = "Temporal activity pattern")
diel.plot.all <- diel.plot.all + geom_tiplab(size = 4, offset = 2.2)
diel.plot.all 

#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Ruminantia_6_state.png", width=46,height=42,units="cm",res=1200)
#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Suina_6_state.png", width=46,height=40,units="cm",res=1200)
#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Tylopoda_6_state.png", width=46,height=40,units="cm",res=200)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Whippomorpha_6_state.png", width=46,height=40,units="cm", res= 1200)

print(diel.plot.all)
dev.off()


# Section 7: Venn diagram of confidence levels ---------------------------------------

library(venneuler)
MyVenn <- venneuler(c(A=86,B=85,C=189,D=92, E = 21, "A&B"=13, 
                      "B&C"=42,"C&D"=53, "D&E"=5,"C&D&E"=5,"C&E"=12))
MyVenn$labels <- c("A\n86","B\n85","C\n189","D\n92", "E\n21")
plot(MyVenn)

png("C:/Users/ameli/OneDrive/Documents/R_projects/fish_sleep/confidence_venn.png", width=20,height=20,units="cm", res= 1200)
MyVenn <- venneuler(c(A=86,B=85,C=189,D=92, E = 21, "A&B"=13, "A&C" = 9, "A&D" = 5, "A&E" = 5, 
                      "B&C"=42,"C&D"=53, "D&E"=7,"C&D&E"=5,"C&E"=12, "A&D&E" = 2, "B&C&E" = 7))
MyVenn$labels <- c("A\n86","B\n85","C\n189","D\n92", "E\n21")
plot(MyVenn)
dev.off()


# Section 8: Crepuscular vs day-night diel plots -------------------------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia data
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#we want to play cath, di, noc preference separately than crep- non crep preference
#so separate the diel pattern 2 column, into its component parts
trait.data <- separate(trait.data, col = Diel_Pattern_2, into = c("Diel", "Crepuscularity"), sep = "/")

custom_colours <- c("#dd8ae7","steelblue1", "#FC8D62", "#66C2A5")

crep_diel_tree <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel", "Crepuscularity")]
crep_diel_tree <- crep_diel_tree + geom_tile(data = crep_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x=x, fill = Diel, width = 4),  inherit.aes = FALSE, color = "transparent")
crep_diel_tree <- crep_diel_tree + geom_tile(data = crep_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x = x + 5, fill = Crepuscularity, width = 4), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom_colours, na.value = "grey80")
crep_diel_tree

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/diel_v_crepuscular.png", width=46,height=40,units="cm", res= 1200)
crep_diel_tree
dev.off()


# # Section 9: Mammal diel plots ------------------------------------------

#chose either Bennie or Maor dataset
trait.data$Diel_Pattern_2 <- trait.data$Bennie_diel
trait.data$Diel_Pattern_2 <- trait.data$Maor_diel

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

custom.colours <- c("#dd8ae7", "pink", "#FC8D62", "#66C2A5")

#plot the tree
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent", width = 4) + scale_fill_manual(values = custom.colours)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 1, offset = 4)
diel.plot.all 
