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

# Section 1: Create diel plots with Cox mammal tree -----------------------

#need to use the mamm tree for modelling and ancestral trait reconstruction
#trees to make: trinary dataset -maxdinoc, trinary dataset -maxcrep, 5state

# ## Read in the mammalian phylogeny
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#use trait data below for Cox data
trait.data <- read.csv(here("cetaceans_full.csv"))

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#does not include species with no activity pattern data (they are filtered out when processing the sleepy cetacean dataframe)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/all_cetacean_tree.png", width=40,height=40, units="cm",res=100)
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)
dev.off()


# Section 2: 5-state trait data tree --------------------------------------
#use for cetaceans
trait.data.all <- read.csv(here("cetaceans_full.csv"))

#use for artiodactyla
trait.data.all <- read.csv(here("cetaceans_full.csv"))

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]

#rename nocturnal crepuscular to crepuscular to create 4 state max crep
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "nocturnal/crepuscular", "crepuscular")

#create a subtree with only the cetacean species
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
custom.colours.all <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
#use for a black and white plot
#custom.colours.all <- c("black", "grey70", "grey40", "grey90")

diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all #+ geom_tiplab(size = 2)
diel.plot.all 

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_4_state.png", width=24,height=18,units="cm",res=1200)
print(diel.plot.all)
dev.off()


# Section 3: maxdinoc tree --------------------------------------

#creating a tree with three trait states (cathemeral, diurnal, nocturnal)
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

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

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_maxdinoc.png", width=24,height=18,units="cm",res=1200)
print(diel.plot.MDN)
dev.off()

# Section 4: maxcrep tree --------------------------------------

#creating a tree with three trait states (crepuscular, diurnal, nocturnal)
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

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

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_maxcrep.png", width=24,height=18,units="cm",res=1200)
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



# Section ?: NOT DONE Create diel plots with Cox mammal tree -----------------------

#doing the same as above but with Cox's mammal tree instead of open tree of life
#need to use the mamm tree for modelling and ancestral trait reconstruction

# ## Read in the mammalian phylogeny
#mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
#mam.tree <- maxCladeCred(mammal_trees, tree = TRUE)

# trait.data <- cetaceans[cetaceans$Diel_Pattern_1 %in% c("diurnal", "nocturnal"),]
# trait.data$tips <- trait.data$Species_name
# trait.data$tips <- str_replace(trait.data$tips, pattern = " ", replacement = "_")
# trait.data2 <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
# 
# "Lagenorhynchus_obscurus" %in% mam.tree$tip.label
# 
# View(trait.data$tips)
# View(mam.tree$tip.label)
# mam.tree$tip.label[(grepl("Balaenoptera", mam.tree$tip.label))]
# 
# #find which species in trait data are missing from mam.tree
# trait.data$tips[!(trait.data$tips %in% mam.tree$tip.label)]
# 
# row.names(trait.data) <- trait.data$tips
# trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# #change this to work with the open tree of life
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# model<- corHMM(phy = trpy_n_mam, data = trait.data[trpy_n_mam$tip.label, c("tips", "diel")], rate.cat = 1, model = "ARD", node.states = "marginal")
# 
# #models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
# #model <- models$HMM_2state_2rate_marg
# 
# lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
# colnames(lik.anc) <- c("diurnal", "nocturnal")
# 
# lik.anc$node <- c(1:length(trpy_n_mam$tip.label), (length(trpy_n_mam$tip.label) + 1):(trpy_n_mam$Nnode + length(trpy_n_mam$tip.label)))
# 
# ancestral_plot <- ggtree(trpy_n_mam, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
# 
# 
# ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)


# Section 6: Create artiodactyla diel plots with Cox mammal tree or Amelia primary source data -----------------------

#need to use the mamm tree for modelling and ancestral trait reconstruction
#trees to make: trinary dataset -maxdinoc, trinary dataset -maxcrep, 5state

# ## Read in the mammalian phylogeny
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#use below for Cox data
#trait.data <- read.csv(here("Cox_artiodactyla_full.csv"))

#use trait data below for Amelia data
trait.data <- read.csv(here("sleepy_artiodactyla.csv"))


trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#currently includes species with no activity pattern data
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)


# Section 7: 6-state trait data tree --------------------------------------

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral, cath/crep)

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
custom.colours.all <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 2)
diel.plot.all 

#png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_artio_6_state.png", width=26,height=22,units="cm",res=1200)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/diel_artio_6_state.png", width=26,height=22,units="cm",res=1200)

print(diel.plot.all)
dev.off()


# Section 8: maxdinoc tree --------------------------------------

#creating a tree with only three trait states (nocturnal, diurnal, cathemeral)

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#remove species with no data
trait.data.maxDN <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_4")]

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxDN$tips)

#plot the tree
custom.colours.MDN <- c("#dd8ae7", "#FC8D62", "#66C2A5")
diel.plot.MDN <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxDN[,c("tips", "Diel_Pattern_4")]
diel.plot.MDN <- diel.plot.MDN + geom_tile(data = diel.plot.MDN$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_4), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MDN)
diel.plot.MDN <- diel.plot.MDN + geom_tiplab(size = 2)
diel.plot.MDN 

#png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_artio_maxdinoc.png", width=27,height=24,units="cm",res=1200)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/diel_artio_maxdinoc.png", width=27,height=24,units="cm",res=1200)

print(diel.plot.MDN)
dev.off()

# Section 9: maxcrep tree --------------------------------------

#creating a tree with three trait states (crepuscular, diurnal, nocturnal)

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#remove species with no data
trait.data.maxcrep <- trait.data[!(is.na(trait.data$Diel_Pattern_3)), c("tips", "Diel_Pattern_3")]

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxcrep$tips)

#plot the tree
custom.colours.MC <- c("pink", "#FC8D62", "#66C2A5")
diel.plot.MC <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxcrep[,c("tips", "Diel_Pattern_3")]
diel.plot.MC <- diel.plot.MC + geom_tile(data = diel.plot.MC$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_3), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MC)
diel.plot.MC <- diel.plot.MC + geom_tiplab(size = 2)
diel.plot.MC 

#png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_artio_maxcrep.png", width=27,height=24,units="cm",res=1200)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/diel_artio_maxcrep.png", width=27,height=24,units="cm",res=1200)

print(diel.plot.MC)
dev.off()



# Section 10: Create artiodactyla without cetaceans diel plots with Cox mammal tree -----------------------

#need to use the mamm tree for modelling and ancestral trait reconstruction
#trees to make: trinary dataset -maxdinoc, trinary dataset -maxcrep, 5state

# ## Read in the mammalian phylogeny
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("Cox_artiodactyla_without_cetaceans.csv"))

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#currently includes species with no activity pattern data
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)


# Section 11: 6-state trait data tree --------------------------------------
#use for cetaceans
trait.data <- read.csv(here("cetaceans_full.csv"))

#use for artiodactyla (including cetaceans)
trait.data <- read.csv(here("Sleepy_artiodactyla_full.csv"))

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral, cath/crep)

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#optional rename nocturnal/crepuscular and diurnal/crepuscular to crepuscular to create 4 state max crep
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "nocturnal/crepuscular", "crepuscular")
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "diurnal/crepuscular", "crepuscular")
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "cathemeral/crepuscular", "crepuscular")

#optional rename nocturnal/crepuscular and diurnal/crepuscular to crepuscular to create 4 state max dinoc
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, "cathemeral/crepuscular", "crepuscular")


#colours for 6 state
#custom.colours.all <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
#colours for 4 state
custom.colours.all <- c ("#934da3", "#faab8c", "#d5cab2", "#5bb095") 

#Plot
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2, width = 4), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all, name = "Temporal activity pattern")
diel.plot.all <- diel.plot.all #+ geom_tiplab(size = 2) 
diel.plot.all 

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_artio_4_state_max_crep.png", width=26,height=22,units="cm",res=1200)
print(diel.plot.all)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_cet_4_state.png", width=26,height=22,units="cm",res=1200)
print(diel.plot.all)
dev.off()


# Section 12: maxdinoc tree --------------------------------------

#creating a tree with only three trait states (nocturnal, diurnal, cathemeral)
trait.data <- read.csv(here("Cox_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#remove species with no data
trait.data.maxDN <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_4")]

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxDN$tips)

#plot the tree
custom.colours.MDN <- c("#dd8ae7", "#FC8D62", "#66C2A5")
diel.plot.MDN <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxDN[,c("tips", "Diel_Pattern_4")]
diel.plot.MDN <- diel.plot.MDN + geom_tile(data = diel.plot.MDN$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_4), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MDN)
diel.plot.MDN <- diel.plot.MDN + geom_tiplab(size = 2)
diel.plot.MDN 

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_artio_minus_cet_maxdinoc.png", width=27,height=24,units="cm",res=1200)
print(diel.plot.MDN)
dev.off()

# Section 13: maxcrep tree --------------------------------------

#creating a tree with three trait states (crepuscular, diurnal, nocturnal)
trait.data <- read.csv(here("Cox_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#remove species with no data
trait.data.maxcrep <- trait.data[!(is.na(trait.data$Diel_Pattern_3)), c("tips", "Diel_Pattern_3")]

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxcrep$tips)

#plot the tree
custom.colours.MC <- c("pink", "#FC8D62", "#66C2A5")
diel.plot.MC <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxcrep[,c("tips", "Diel_Pattern_3")]
diel.plot.MC <- diel.plot.MC + geom_tile(data = diel.plot.MC$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_3), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MC)
diel.plot.MC <- diel.plot.MC + geom_tiplab(size = 2)
diel.plot.MC 

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_artio_minus_cet_maxcrep.png", width=27,height=24,units="cm",res=1200)
print(diel.plot.MC)
dev.off()


# Section X: Artiodactyla suborders ---------------------------------------

# ## Read in the mammalian phylogeny
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia data
trait.data <- read.csv(here("sleepy_artiodactyla.csv"))
trait.data <- trait.data[!trait.data$tips %in% mam.tree$tip.label,]

#filter for the suborder you want
#trait.data <- trait.data %>% filter(Order == "Ruminantia")
#trait.data <- trait.data %>% filter(Order == "Suina")
#trait.data <- trait.data %>% filter(Order == "Tylopoda")
#trait.data <- trait.data %>% filter(Order == "Whippomorpha")

trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
custom.colours.all <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "grey")
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 2)
diel.plot.all 

#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Ruminantia_6_state.png", width=46,height=42,units="cm",res=1200)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Suina_6_state.png", width=46,height=40,units="cm",res=1200)
#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Tylopoda_6_state.png", width=46,height=40,units="cm",res=200)
#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Whippomorpha_6_state.png", width=46,height=40,units="cm", res= 1200)

print(diel.plot.all)
dev.off()

