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


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/Fish_sleep_functions.R")

## Load files

resolved_names <- readRDS("resolved_names_AllGroups.rds")
tr.calibrated <- readRDS("tr_tree_calibrated_AllGroups.rds")
trait.data <- readRDS("trait_data_AllGroups.rds")

name_variable <- "AllGroups"


# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]


# Load the best model, which is the HMM 2 state 2 rate model
models <- readRDS(file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
model <- models[[3]]

###############################################################################################################################################
### Create the object with lineages through time and diel switches data ### 
###############################################################################################################################################

# First, extract the ancestral states from the best fit model
anc_states <- returnAncestralStates(phylo_model = model, phylo_tree = trpy_n)

# Then, calculate transitions between states
anc_states <- calculateStateTranstitions(ancestral_states = anc_states, phylo_tree = trpy_n)

# Determine transition histories (types of lineages)
anc_states <- calculateLinTransHist(ancestral_states = anc_states, phylo_tree = trpy_n)

# Calculate cumsums through time (for ltt plots)
anc_states <- returnCumSums(ancestral_states = anc_states, phylo_tree = trpy_n)


###############################################################################################################################################
### Make Plots! ### 
###############################################################################################################################################

switch.histo <- switchHisto(ancestral_states = anc_states, replace_variable_names = T, backfill = T)

# So, this works now, and I think correctly
# One question is though, because all the tips are clustered at the same age, there's a lot of tips and transitions right near the end of the tree
# Transistions are associated with the tip or node that is reconstructed, which is why this appears like this
# I could associate them with the parental node, but I don't have evidence that that is when it was reconstructed (which would fix the artifact at the end)

switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02)

# This highlights which nodes have undergone the most transitions
numb_switch_tree <- switchTree(ancestral_states = anc_states, phylo_tree = trpy_n, layout = "circular", replace_variable_names = TRUE)

# Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void()



###############################################################################################################################################
### Save Plots! ### 
###############################################################################################################################################


# This really shows the difference between methods, and where to assign the switch
# Also super affected by the time calibration

xlims <- c(max(anc_states$node.age), min(anc_states$node.age))

pdf(file = paste("outs/Figures/AllGroups_phylogeny_diel_plot_transitions_histo", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
((switch.histo / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5))
dev.off()

pdf(file = paste("outs/Figures/AllGroups_phylogeny_diel_plot_transitions_swith", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
((switch.ratio / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5))
dev.off()


pdf(file = paste("outs/Figures/AllGroups_phylogeny_diel_highswitchlineages", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), height = 60, width = 60)
numb_switch_tree
dev.off()


