library(phangorn)
library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(zoo)
library(ggtree)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/Fish_sleep_functions.R")

index_list <- list()
index_list[[1]] <- c("all", "only_highqual", "only_cartilaginous", "only_ingroup")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

recon <- "marg"

for (i in 1:length(index_list)) {
  dataset_variable <- names(index_list)[[i]]
  
  for (j in 1:length(index_list[[i]])) {
    name_variable <- index_list[[i]][[j]]
    
    ## Load in the tree
    trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
    
    ## Load the best model, which is the HMM 2 state 2 rate model
    models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
    
    ## I think I will always take the 2 rate model, 3 is too hard to comprehend
    
    ## If I want to use the joint reconstruction then the structure of the data changes (marginal gives % for tips + nodes, joint gives states as vectors)
    
    if(recon == "joint") {
      model <- models$HMM_2state_2rate_joint
    } else {
      model <- models$HMM_2state_2rate_marg
    }
    
    # model <- models$MK_2state_marg
    
    
    # View(unlist(lapply(models, function(x) x[names(x[grep("loglik", names(x))])])))
    
    ###############################################################################################################################################
    ### Create the object with lineages through time and diel switches data ### 
    ###############################################################################################################################################
    
    # First, extract the ancestral states from the best fit model
    anc_states <- returnAncestralStates(phylo_model = model, phylo_tree = trpy_n, rate.cat = 2, recon = recon)
    
    # Then, calculate transitions between states (or rate categories)
    anc_states <- calculateStateTransitions(ancestral_states = anc_states, phylo_tree = trpy_n, rate.cat = F)
    
    # Determine transition histories (types of lineages)
    anc_states <- calculateLinTransHist(ancestral_states = anc_states, phylo_tree = trpy_n)
    
    # Calculate cumsums through time (for ltt plots)
    anc_states <- returnCumSums(ancestral_states = anc_states, phylo_tree = trpy_n)
    
    ## Save out the file (to be re-used in #5)
    saveRDS(anc_states, file = paste(dataset_variable, "diel_ancestral_states", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
    
    anc_states <- readRDS(file = paste(dataset_variable, "diel_ancestral_states", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
    
    
    ###############################################################################################################################################
    ### Make Plots! ### 
    ###############################################################################################################################################
    
    switch.histo <- switchHisto(ancestral_states = anc_states, replace_variable_names = F, backfill = T)
    
    # So, this works now, and I think correctly
    # One question is though, because all the tips are clustered at the same age, there's a lot of tips and transitions right near the end of the tree
    # Transistions are associated with the tip or node that is reconstructed, which is why this appears like this
    # I could associate them with the parental node, but I don't have evidence that that is when it was reconstructed (which would fix the artifact at the end)
    
    switch.ratio.types <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02, use_types = T)
    switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02, use_types = F)
    
    
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
    
    pdf(file = paste("outs/Figures/phylogeny_diel_plot_transitions_histo", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
    print(((switch.histo / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5)))
    dev.off()
    
    pdf(file = paste("outs/Figures/phylogeny_diel_plot_transitions_swithTypes", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
    print(((switch.ratio.types / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5)))
    dev.off()
    
    pdf(file = paste("outs/Figures/phylogeny_diel_plot_transitions_swith", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
    print(((switch.ratio / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5)))
    dev.off()
    
    pdf(file = paste("outs/Figures/phylogeny_diel_highswitchlineages", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), height = 60, width = 60)
    print(numb_switch_tree)
    dev.off()
    
    
    ############# Do it for rates (not states)
    
    
    ###############################################################################################################################################
    ### Create the object with lineages through time and diel switches data ### 
    ###############################################################################################################################################
    
    anc_rates <- returnAncestralStates(phylo_model = model, phylo_tree = trpy_n)
    anc_rates <- calculateStateTransitions(ancestral_states = anc_rates, phylo_tree = trpy_n, rate.cat = T)
    anc_rates <- calculateLinTransHist(ancestral_states = anc_rates, phylo_tree = trpy_n)
    anc_rates <- returnCumSums(ancestral_states = anc_rates, phylo_tree = trpy_n)
    
    ## Save out the file (to be re-used in #5)
    saveRDS(anc_rates, file = paste(dataset_variable, "diel_ancestral_rates", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
    anc_rates <- readRDS(file = paste(dataset_variable, "diel_ancestral_rates", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
    
    
    ###############################################################################################################################################
    ### Make Plots! ### 
    ###############################################################################################################################################
    
    switch.ratio <- switchRatio(ancestral_states = anc_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, use_types = F)
    geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void()
    
    ###############################################################################################################################################
    ### Save Plots! ### 
    ###############################################################################################################################################
    
    xlims <- c(max(anc_rates$node.age), min(anc_rates$node.age))
    
    pdf(file = paste("outs/Figures/phylogeny_diel_plot_transitions_rates_switch", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
    print(((switch.ratio / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5)))
    dev.off()
    
  }
  
}










