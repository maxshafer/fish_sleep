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


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/Fish_sleep_functions.R")

index_list <- list()
index_list[[1]] <- c("all", "only_highqual", "only_cartilaginous", "only_ingroup")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

# Set simulation parameters (type and #)
model_types <- c("ARD", "HR")
sim_numb <- 500

joint <- FALSE

for (y in 2:length(model_types)) {
  model_type <- model_types[[y]]
  
  for (i in 1:length(index_list)) {
    dataset_variable <- names(index_list)[[i]]
    
    for (j in 1:length(index_list[[i]])) {
      name_variable <- index_list[[i]][[j]]
      
      setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")
      
      ## Load in the tree
      trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
      
      ## Load the best model, which is the HMM 2 state 2 rate model
      models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
      
      ## I think I will always take the 2 rate model, 3 is too hard to comprehend
      
      ## If I want to use the joint reconstruction then the structure of the data changes (marginal gives % for tips + nodes, joint gives states as vectors)
      if (joint == TRUE) {
        model <- models$HMM_2state_2rate_joint
        model_ER <- models$MK_2state_joint
      } else {
        model <- models$HMM_2state_2rate_marg
        model_ER <- models$MK_2state_marg
      }
      
      ## Load the ancestral states data
      anc_states <- readRDS(file = paste(dataset_variable, "diel_ancestral_states", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      anc_rates <- readRDS(file = paste(dataset_variable, "diel_ancestral_rates", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      
      ###############################################################################################################################################
      ### Simulate diel switchs based on derived model parameters ### 
      ###############################################################################################################################################
      
      ## This runs a simulation based on some arguments (wrapper for rTraitDisc)
      ## I can pull the rates directly from the models I've loaded above
      
      # if (model_type == "ER") {
      #   model_rates <- models$MK_2state_marg$solution[2]
      #   simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = model_rates, states = c("diurnal", "nocturnal"), simulation_numb = sim_numb)
      # }
      # 
      # if (model_type == "ARD") {
      #   model_rates <- models$MK_2state_marg$solution[2:3]
      #   simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = model_rates, states = c("diurnal", "nocturnal"), simulation_numb = sim_numb, root_state = "nocturnal")
      # }
      # 
      # if (model_type == "HR") {
      #   model_rates <- model$solution
      #   simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = model_rates, states = c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2"), simulation_numb = sim_numb, root_state = "nocturnal_R2")
      # }
      # 
      # saveRDS(simulation, file = paste(dataset_variable, "diel_switch_simulations", name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_"))
      
      simulation <- readRDS(file = paste(dataset_variable, "diel_switch_simulations", name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_"))
      
      ## This calculates for each node whether it is a switch from it's parental node
      
      if (model_type == "HR") {
        # Make different for rates vs states
        simulation_rates <- simulation
        
        # Replace numbers with states
        simulation[simulation == "diurnal_R1"] <- "diurnal"
        simulation[simulation == "diurnal_R2"] <- "diurnal"
        simulation[simulation == "nocturnal_R1"] <- "nocturnal"
        simulation[simulation == "nocturnal_R2"] <- "nocturnal"
        
        simulation_rates[simulation_rates == "diurnal_R1"] <- "R1"
        simulation_rates[simulation_rates == "nocturnal_R1"] <- "R1"
        simulation_rates[simulation_rates == "diurnal_R2"] <- "R2"
        simulation_rates[simulation_rates == "nocturnal_R2"] <- "R2"
        
      } 
      
      simulated_transitions <- calculateSimulatedTransitions(simulated_data = simulation, phylo_tree = trpy_n)
      
      ## Function that returns cumsums, rows are now ordered by node.age, whether or it is included as a column in the output
      cumsums <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions, phylo_tree = trpy_n, include_node_age = TRUE)
      
      if (model_type == "HR") {
        simulated_transitions_rates <- calculateSimulatedTransitions(simulated_data = simulation_rates, phylo_tree = trpy_n)
        cumsums_rates <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_rates, phylo_tree = trpy_n, include_node_age = TRUE)
      }
      
      ###############################################################################################################################################
      ### Plot the results of the simulation, along with the actual data ### 
      ###############################################################################################################################################
      
      ## Plot the results of the simulation, either the summary (mean +/- stdev), the simulations, or both ('overlay')
      plot <- simulatedSwitchRatio(simulated_cumsums = cumsums, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "red", error = "both") #+ ylim(c(0,0.15))
      
      ## Can I add the computed value? Actual transitions
      ## Easy way is to make the plot, then use the data from it to add to the above
      switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02)
      
      plot <- plot + geom_line(data = switch.ratio$data, aes(x=node.age,y=ratio), colour = "black")
      plot <- plot + ggtitle(paste(model_type, "model simulation,", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", ncol(simulation), ", Avg. transitions:", mean(colSums(simulated_transitions)), "+/-", round(sd(colSums(simulated_transitions)),1), sep = " "))
      
      # This 'zooms-in', which preserves the ribbon (setting limits removes part of the ribbon stdev)
      recon_plot <- plot + coord_cartesian(ylim = c(0, layer_scales(plot)$y$range$range[[2]]), xlim = abs(layer_scales(plot)$x$range$range)) 
      
      # Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
      geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot)$x$range$range)) 
      
      
      ## Save out the plot
      setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/outs/Figures")
      
      pdf(file = paste(dataset_variable, "phylogeny_diel_plot_transitions_simulation", name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simulation), "zoom.pdf", sep = "_"), width = 10, height = 10)
      print(recon_plot / geo_scale + plot_layout(nrow = 2, heights = c(5,0.5)))
      dev.off()
      
      if (model_type == "HR") {
        plot_simulation_rates <- simulatedSwitchRatio(simulated_cumsums = cumsums_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "red", error = "both") #+ ylim(c(0,0.15))
        switch.ratio.rates <- switchRatio(ancestral_states = anc_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02)
        plot.rates <- plot_simulation_rates + geom_line(data = switch.ratio.rates$data, aes(x=node.age,y=ratio), colour = "black")
        plot.rates <- plot.rates + ggtitle(paste(model_type, "model simulation (rates),", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", ncol(simulation), ", Avg. transitions:", mean(colSums(simulated_transitions_rates)), "+/-", round(sd(colSums(simulated_transitions_rates)),1), sep = " "))
        recon_plot_rates <- plot.rates + coord_cartesian(ylim = c(0, layer_scales(plot.rates)$y$range$range[[2]]), xlim = abs(layer_scales(plot.rates)$x$range$range)) 
        geo_scale_rates <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot.rates)$x$range$range))
        
        pdf(file = paste(dataset_variable, "phylogeny_diel_plot_transitions_simulation_rates", name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simulation), "zoom.pdf", sep = "_"), width = 10, height = 10)
        print(recon_plot_rates / geo_scale_rates + plot_layout(nrow = 2, heights = c(5,0.5)))
        dev.off()
      }
      
    }
  }
  
}





