# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
#install.packages("deeptime")
library(deeptime)
#update.packages("ggplot2")
library(ggplot2)
setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Plotting the cetacean diel dataset -----------------------
cetaceans_full <- read.csv(here("cetaceans_full.csv"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#all diel patterns
custom.colours <- c("#dd8ae7", "peachpuff2" ,"#FC8D62", "yellow", "#66C2A5", "#A6D854","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 1)
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/six_state_plot_with_Na_labelled.pdf")
diel.plot
dev.off() 

#maximum crepuscular dataset, with NAs, species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 1)
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_with_Na_labelled.pdf")
diel.plot
dev.off() 

#maximum crepuscular dataset with NAs, no species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_with_Na_unlabelled.pdf")
diel.plot
dev.off() 

cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$Diel_Pattern), ]
trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#maximum crepuscular dataset, without NAs, species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 1)
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_labelled.pdf")
diel.plot
dev.off() 

#maximum crepuscular dataset without NAs, no species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_unlabelled.pdf")
diel.plot
dev.off()

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
custom.colours.2 <- c("grey90", "grey80","grey60", "grey30", "grey1")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "Confidence")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = as.factor(Confidence)), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Evidence category")
diel.plot <- diel.plot
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_unlabelled_confidence.pdf")
diel.plot
dev.off()

cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#add clade labels
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = "Mysticeti")
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = "Odontoceti")

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + 
  geom_cladelab(node = 137, label = "Mysticeti", align = TRUE, geom = "label", offset=1, align=TRUE, offset.text=1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot <- diel.plot + 
  geom_cladelab(node = 77, label = "Odontoceti", align = FALSE, geom = "label", offset=1, align=FALSE, offset.text=1, hjust = 1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_unlabelled_cladelabels.pdf")
diel.plot
dev.off()

findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = "Mysticeti")
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = "Odontoceti")

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + 
  geom_cladelab(node = 137, label = "Mysticeti", align = TRUE, geom = "label", offset=1, align=TRUE, offset.text=1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot <- diel.plot + 
  geom_cladelab(node = 77, label = "Odontoceti", align = FALSE, geom = "label", offset=1, align=FALSE, offset.text=1, hjust = 1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot

# Section 2: Phylogenetic signal -----------------

#can we do this for discrete traits as well? 
# some people use delta statistic https://github.com/mrborges23/delta_statistic 
source("scripts/Amelia_delta_code.R")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("cetaceans_full.csv"))

#keep only necessary columns
trait.data <- trait.data[, c("Species_name", "max_crep", "tips")]
trait.data <- trait.data[trait.data$max_crep %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
#this doesn't replace anything in the cetacean tree so comment it out for now
mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#vector of trait data, with species in same order as in tree (tree$tip.label)
sps_order <- as.data.frame(mam.tree$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
trait.data <- merge(trait.data, sps_order, by = "tips")
trait.data <- trait.data[order(trait.data$id), ]
trait <- trait.data$max_crep

#now we calculate delta using their custom function
delta_diel <- delta(trait, mam.tree, 0.1, 0.0589, 1000, 10, 100)
#returns a value of 0.7779191, will change slightly every time its calculated
#significance is only determined relative to a random simulation of the trait data

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,mam.tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>delta_diel)/length(random_delta)
boxplot(random_delta)
abline(h=delta_diel,col="red")

#if p-value less than 0.05 there is evidence of phylogenetic signal between the trait and character
#if its more than 0.05 there is not evidence for phylogenetic signal
#p value is 0.04, so there is a phylogenetic signal for diel activity patterns in cetaceans

