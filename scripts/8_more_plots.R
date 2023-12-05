library(phangorn)
library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(zoo)
library(ggtree)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)
library(xlsx)
library(gsheet)
library(RColorBrewer)
library(viridis)
library(scales)

#### This script to make some extra figures, based on extinction, fossil, and other data

setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

dataset_variable <- "fish"
name_variable <- "all"

## Load in the tree
trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c, only_model_data = TRUE)
trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable, only_model_data = TRUE)

models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
model <- models$HMM_2state_2rate_marg

########## ########## ########## ########## ########## 
########## Make radial plot ########## ########## 
########## ########## ########## ########## ########## 

lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2", "diurnal_R3", "nocturnal_R3", "diurnal_R4", "nocturnal_R4")[1:ncol(lik.anc)]

lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
lik.anc$noc.sum <- rowSums(lik.anc[,grep("nocturnal", colnames(lik.anc))])
lik.anc$di.sum <- rowSums(lik.anc[,grep("diurnal", colnames(lik.anc))])
lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])

ancestral_plot <- ggtree(trpy_n, layout = "fan", open.angle = 180, size = 0.5) %<+% lik.anc + aes(color = di.sum - noc.sum) + scale_color_distiller(palette = "RdBu") + geom_nodelab(size = 0)

orders <- c("Anguilliformes", "Clupeiformes", "Gymnotiformes", "Cypriniformes", "Characiformes", "Siluriformes",  
            "Salmoniformes", "Myctophiformes", "Gadiformes", "Holocentriformes", "Ophidiiformes",
            "Centrarchiformes", "Cottoidei", "Lophiiformes", "Gobiesociformes",
            "Kurtiformes", "Gobiiformes", "Syngnathiformes", "Pleuronectiformes", 
            "Cichliformes", "Perciformes/Scorpaenoidei", "Tetraodontiformes", 
            "Carcharhiniformes", "Holocentriformes", "Blenniiformes")

# orders <- unique(trait.data_n$order)
my_colors <- hue_pal()(length(orders))
ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral_plot, colours = my_colors, orders = orders, resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 5, alpha = 0.25)

pdf(file = here("outs/Figures/plot_XX_fish_tree.pdf"), width = 8, height = 8)
ancestral.plot
dev.off()

### For rates

ancestral_plot <- ggtree(trpy_n, layout = "fan", open.angle = 180, size = 0.5) %<+% lik.anc + aes(color = R1.sum) + scale_color_distiller(palette = "PRGn") + geom_nodelab(size = 0)
ancestral_plot <- ancestral_plot + guides(color=guide_legend(title="di.sum - noc.sum"))

orders <- c("Anguilliformes", "Clupeiformes", "Gymnotiformes", "Cypriniformes", "Characiformes", "Siluriformes",  
            "Salmoniformes", "Myctophiformes", "Gadiformes", "Holocentriformes", "Ophidiiformes",
            "Centrarchiformes", "Cottoidei", "Lophiiformes", "Gobiesociformes",
            "Kurtiformes", "Gobiiformes", "Syngnathiformes", "Pleuronectiformes", 
            "Cichliformes", "Perciformes/Scorpaenoidei", "Tetraodontiformes", 
            "Carcharhiniformes", "Holocentriformes", "Blenniiformes")

# orders <- unique(trait.data_n$order)
my_colors <- hue_pal()(length(orders))
ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral_plot, custom = "no", colours = my_colors, orders = orders, resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 5, alpha = 0.25)

pdf(file = here("outs/Figures/plot_XX_fish_tree_rates.pdf"), width = 8, height = 8)
ancestral.plot
dev.off()

########## ########## ########## ########## ########## 
########## Make circular plot ########## ########## 
########## ########## ########## ########## ########## 


dataset_variable <- "AllGroups"
name_variable <- "all"

## Load in the tree
trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c, only_model_data = TRUE)
trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable, only_model_data = TRUE)

models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
model <- models$HMM_2state_2rate_marg

########## ########## ########## ########## ########## 
########## Make radial plot ########## ########## 
########## ########## ########## ########## ########## 

lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2", "diurnal_R3", "nocturnal_R3", "diurnal_R4", "nocturnal_R4")[1:ncol(lik.anc)]

lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
lik.anc$noc.sum <- rowSums(lik.anc[,grep("nocturnal", colnames(lik.anc))])
lik.anc$di.sum <- rowSums(lik.anc[,grep("diurnal", colnames(lik.anc))])
lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])

ancestral_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + scale_color_distiller(palette = "RdBu") + geom_nodelab(size = 0)

## Add clade labels


ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral_plot, custom = c("Polypterus_delhezi", "Lutjanus_fulvus"), label = c("Osteichthyes"), colours = c("black"), resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 20, alpha = 0.25)
ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral.plot, custom = c("Pristis_pristis", "Prionace_glauca"), label = c("Chondrichthyes"), colours = c("black"), resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 20, alpha = 0.25)
ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral.plot, custom = c("Ornithorhynchus_anatinus", "Mus_mahomet"), label = c("Mammalia"), colours = c("black"), resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 20, alpha = 0.25)
ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral.plot, custom = c("Sphenodon_punctatus", "Melospiza_georgiana"), label = c("Sauropsida"), colours = c("black"), resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 20, alpha = 0.25)
ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral.plot, custom = c("Scolecomorphus_vittatus", "Smilisca_cyanosticta"), label = c("Amphibia"), colours = c("black"), resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 20, alpha = 0.25)



pdf(file = here("outs/Figures/plot_XX_AllGroups_tree.pdf"), width = 2.5, height = 2.5)
ancestral.plot + theme(legend.position = "none")
dev.off()


