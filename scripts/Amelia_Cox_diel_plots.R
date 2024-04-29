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
library(tidyr)

setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Create diel plots with Cox mammal tree -----------------------

#need to use the mamm tree for modelling and ancestral trait reconstruction
#trees to make: trinary dataset -maxdinoc, trinary dataset -maxcrep, 5state

# ## Read in the mammalian phylogeny
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("cetaceans_full.csv"))

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#currently includes species with no activity pattern data
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)


# Section 2: 5-state trait data tree --------------------------------------

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
custom.colours.all <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 2)
diel.plot.all 

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_5_state.png", width=24,height=18,units="cm",res=1200)
print(diel.plot.all)
dev.off()

#same plot with NA species included
#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

#remove species with no data
trait.data.all <- trait.data[, c("tips", "Diel_Pattern_2")]
trait.data.all$Diel_Pattern_2 <- trait.data.all$Diel_Pattern_2 %>% replace_na("data unavailable")

trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
custom.colours.all <- c("#dd8ae7", "lightgrey", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 3)
diel.plot.all 

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/diel_plot_5_state_missing_sps.png", width=40,height=30,units="cm",res=1200)
print(diel.plot.all)
dev.off()



# Section 3: maxdinoc tree --------------------------------------

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

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

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

#remove species with no data
trait.data.maxcrep <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "nocturnal/crepuscular", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "diurnal/crepuscular", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "cathemeral", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "crep/cath", "crepuscular/cathemeral")

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxcrep$tips)

#plot the tree
custom.colours.MC <- c("#dd8ae7", "#FC8D62", "#66C2A5")
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
density_tree50 <- ggdensitree(fifty_trees, alpha = 0.3, colour = "blue") +  geom_tiplab(color = "black")
density_tree50 

sample_trees <- cetacean_trees[sample(1:length(cetacean_trees), 10)]
densitree <- ggdensitree(sample_trees, alpha = 0.3, colour = "pink") +  geom_tiplab(color = "black")

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/density_tree100.png", width=40,height=30,units="cm",res=1200)
densitree
dev.off()


# Section 4: NOT DONE Create diel plots with Cox mammal tree -----------------------

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
custom.colours.all <- c("#dd8ae7", "yellow", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

trait.data <- read.csv("~/R_projects/fish_sleep/artiodactyla_full.csv")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_a <- keep.tip(mam.tree, tip = trait.data$tips)
artio_tree<- ggtree(trpy_a, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2")]
artio_tree <- artio_tree +  geom_tile(data = artio_tree$data[1:length(trpy_a$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all) +
  geom_tiplab(size = 2) + geom_strip('Sus_celebensis', 'Catagonus_wagneri', barsize=2, color='tomato1', label="Suiformes", offset.text=2, offset = 20) + 
  geom_strip('Camelus_dromedarius', 'Vicugna_vicugna', barsize=2, color='tomato2', label="Tylopoda", offset.text=2, offset = 20) +
  geom_strip('Hippopotamus_amphibius', 'Stenella_coeruleoalba', barsize=2, color='tomato3', label="Whippomorpha", offset.text=2, offset = 21) +
  geom_strip('Hyemoschus_aquaticus', 'Capra_caucasica', barsize=2, color='tomato4', label="Ruminantia", offset.text= 6, offset = 20) 

artio_tree

png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/artio_diel.png", width=40,height=30,units="cm",res=1200)
artio_tree
dev.off()




# Diel data table ---------------------------------------------------------

trait.data <- read.csv(here("cetaceans_full.csv"))
#add breakdown of families
test <- diel3_families$data
test <- test[1:72,c("label")]
test$family <- "blank"
test[1:2, "family"] <- "Platanistidae"
test[3, "family"] <- "Lipotidae"
test[4:6, "family"] <- "Iniidae"
test[7:41, "family"] <- "Delphinidae"
test[42:46, "family"] <- "Phocoenidae"
test[47:48, "family"] <- "Monodontidae"
test[49:59, "family"] <- "Ziphiidae"
test[60, "family"] <- "Physteridae"
test[61:62, "family"] <- "Kogiidae"
test[63:71, "family"] <- "Balaenopteridae"
test[72, "family"] <- "Eschrichitiidae"
test[73:77, "family"] <- "Balaenidae"
test <- col.names()

df_merge <- merge(test,cetaceans_full,by="")

cet_sum <- table(trait.data$Diel_Pattern_2)

