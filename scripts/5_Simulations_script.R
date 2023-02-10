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


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/fish_sleep_functions.R")

# Which tree?
name_variable <- "all"

## Load in the tree
trpy_n <- loadTree(return = "tree", dataset = "fish", subset = name_variable, custom_tips = NA)

## Load the best model, which is the HMM 2 state 2 rate model
# models <- readRDS(file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
model <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

## Load the ancestral states data
anc_states <- readRDS(file = paste("fish_diel_ancestral_states", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))

###############################################################################################################################################
### Simulate diel switchs based on derived model parameters ### 
###############################################################################################################################################

## This runs a simulation based on some arguments (wrapper for rTraitDisc)
## I can pull the rates directly from the models I've loaded above
model_type = "ER"
sim_numb <- 1000
simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = models$fitER$rates, states = c("nocturnal", "diurnal"), simulation_numb = sim_numb)

saveRDS(simulation, file = paste("fish_diel_switch_simulations", name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_"))

simulation <- readRDS(file = paste("fish_diel_switch_simulations", name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_"))
## This calculates for each node whether it is a switch from it's parental node
simulated_transitions <- calculateSimulatedTransitions(simulated_data = simulation, phylo_tree = trpy_n)

## Function that returns cumsums, rows are now ordered by node.age, whether or it is included as a column in the output
cumsums <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions, phylo_tree = trpy_n, include_node_age = TRUE)


###############################################################################################################################################
### Plot the results of the simulation, along with the actual data ### 
###############################################################################################################################################

## Plot the results of the simulation, either the summary (mean +/- stdev), the simulations, or both ('overlay')
plot <- simulatedSwitchRatio(simulated_cumsums = cumsums, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "black") #+ ylim(c(0,0.15))

## Can I add the computed value? Actual transitions
## Easy way is to make the plot, then use the data from it to add to the above
switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02)

plot <- plot + geom_line(data = switch.ratio$data, aes(x=node.age,y=ratio), colour = "red")
plot <- plot + ggtitle(paste(model_type, "model simulation,", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", ncol(simulation), ", Avg. transitions:", mean(colSums(simulated_transitions)), "+/-", round(sd(colSums(simulated_transitions)),1), sep = " "))

# This 'zooms-in', which preserves the ribbon (setting limits removes part of the ribbon stdev)
recon_plot <- plot + coord_cartesian(ylim = c(0, layer_scales(plot)$y$range$range[[2]]), xlim = abs(layer_scales(plot)$x$range$range)) 

# Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot)$x$range$range)) 


## Save out the plot

pdf(file = paste("outs/Figures/fish_phylogeny_diel_plot_transitions_simulation", name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simulation), "zoom.pdf", sep = "_"), width = 10, height = 10)
recon_plot / geo_scale + plot_layout(nrow = 2, heights = c(5,0.5))
dev.off()


