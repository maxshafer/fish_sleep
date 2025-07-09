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

#maximum crepuscular dataset, with NAs, species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 1)
diel.plot

#maximum crepuscular dataset with NAs, no species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot
diel.plot

cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$Diel_Pattern), ]
trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#maximum crepuscular dataset, without NAs, species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 1)
diel.plot

#maximum crepuscular dataset without NAs, no species labels
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot
diel.plot

# Section 2: Plotting the artiodactyla diel dataset  -----------------


