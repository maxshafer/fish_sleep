library(phytools)
library(corHMM)
library(ggtree)
library(tidyr)
library(dplyr)
library(tibble)
library(phangorn)
library(ggplot2)
library(patchwork)


##### This script runs makeSimmap from the corrHMM package, to generate stochastic character maps from the distribution of
##### states from either an ARD model, or a HR model (best fitting model in most cases)
##### This then saves out the simmaps (typically 500 or 1000), which can be used downstream in plotting
##### This also creates some basic figures, which show the mean +/- SD of cumulative transitions

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/fish_sleep_functions.R")


index_list <- list()
index_list[[1]] <- c("all", "only_highqual", "only_cartilaginous", "only_ingroup")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

# Set simulation parameters (type and #)
model_types <- c("ARD", "HR")
sim_numb <- 500

## This doesn't actually change anything, because both reconstruction methods use the same rates
## It would just change which reconstruction to plot along side the SIMMAPs
joint <- FALSE


for (y in 1:length(model_types)) {
  model_type <- model_types[[y]]
  
  for (i in 1:length(index_list)) {
    dataset_variable <- names(index_list)[[i]]
    
    for (j in 1:length(index_list[[i]])) {
      name_variable <- index_list[[i]][[j]]
      
      setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")
      
      ## Load in the tree
      trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
      
      trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)
      
      ## Load the best model, which is the HMM 2 state 2 rate model
      models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
      
      ## I think I will always take the 2 rate model, 3 is too hard to comprehend
      
      ## If I want to use the joint reconstruction then the structure of the data changes (marginal gives % for tips + nodes, joint gives states as vectors)
      if (joint == TRUE) {
        model <- models$HMM_2state_2rate_joint
        model_ARD <- models$MK_2state_joint
      } else {
        model <- models$HMM_2state_2rate_marg
        model_ARD <- models$MK_2state_marg
      }
      
      ## Load the ancestral states and rates data
      anc_states <- readRDS(file = paste(dataset_variable, "diel_ancestral_states", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      anc_rates <- readRDS(file = paste(dataset_variable, "diel_ancestral_rates", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      
      ### Run simmap

      if (model_type == "ARD") {
        simmaps <- corHMM::makeSimmap(tree = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = model_ARD$solution, nSim = sim_numb, nCores = 5)
      }

      if (model_type == "HR") {
        simmaps <- corHMM::makeSimmap(tree = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = model$solution, nSim = sim_numb, nCores = 5)
      }

      ## Save SIMMAP

      saveRDS(simmaps, file = paste(dataset_variable, "diel_switch_SIMMAP", name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "corrHMM.rds", sep = "_"))

      simmaps <- readRDS(file = paste(dataset_variable, "diel_switch_SIMMAP", name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "corrHMM.rds", sep = "_"))
      
      
      ### Extract the node states from the simmaps
      ### This was my code, but phytools already has an implementation
      
      # simmap_simulations <- lapply(seq_along(simmaps), function(x) {
      # 
      #   simtree <- simmaps[[x]]
      #   # This identifies the end node for each edge (the second column in simtree$edge, which is ordered by edge #)
      #   # The root doesn't have an edge leading to it, so we have to add it here as well
      # 
      #   node_states <- data.frame(edge = 1:nrow(simtree$edge), tips = simtree$edge[,1], state = unlist(lapply(1:length(simtree$maps), function(x) names(simtree$maps[x][[1]])[length(simtree$maps[x])])))
      #   # node_states <- node_states[unique(node_states$tips),]
      # 
      #   tip_states <- data.frame(edge = 1:nrow(simtree$edge), tips = simtree$edge[,2])
      #   tip_states$state <- unlist(lapply(tip_states$edge, function(x) names(simtree$maps[x][[1]])[length(simtree$maps[x])]))
      # 
      #   # Add root
      #   root_edge <- grep(Ntip(trpy_n)+1, simtree$edge[1,])
      #   tip_states <- rbind(tip_states, data.frame(edge = NA, tips = Ntip(trpy_n)+1, state = names(simtree$maps[root_edge][[1]])[1]))
      # 
      #   # These are wrong for the tips (I think it's the last simulated data, not necessarily what is at the tip)
      #   # Duh! This is also the same for all edges (node values are the beginning of the edge?)
      #   # So replace tip values with actual tip values (this is what the plots do)
      #   tip_states$state[match(1:Ntip(trpy_n), tip_states$tips)] <- trait.data_n$diel1
      # 
      #   ## Replace node state values in tip_states with those from node_states (should be different)
      #   ## Match by edges (which should match)
      #   tip_states$state[match(node_states$tips, tip_states$tips)] <- node_states$state
      # 
      # 
      #   ## I need to convert tips to tip_names and node_x
      #   tip_states <- tip_states[order(tip_states$tips),]
      #   tip_states$node <- c(trpy_n$tip.label, (Ntip(trpy_n)+1):(Ntip(trpy_n)+Nnode(trpy_n)))
      # 
      #   ## Convert to character names
      #   # tip_states$state <- ifelse(tip_states$state == 1, "diurnal", "nocturnal")
      # 
      #   colnames(tip_states) <- c("edge", "tips", paste("state", x, sep = "_"), "node")
      # 
      #   return(tip_states[,c(4,3)])
      # 
      # })
      
      ## This uses phytools function to extract tip and node states (should have looked for this before)
      ## However, it gives me different downstream results (maybe because of the order!)
      simmap_simulations <- lapply(seq_along(simmaps), function(x) {
        df <- getStates(simmaps[[x]], type ="both")
        df <- data.frame(nodes = names(df), state = df)
        colnames(df) <- c("node", paste("state", x, sep = "_"))
        return(df)
        }
        )
      
      # Add the marginal reconstruction states, to see if it is calculated the same way
      marg_states <- apply(anc_states$lik.anc, 1, function(x) which.max(x))
      marg_states <- data.frame(node = names(marg_states), marg_states = marg_states)
      marg_states <- marg_states[row.names(simmap_simulations[[1]]),]
      
      
      # simmap_simulations[[501]] <- marg_states
      
      # Then full_join them together!
      simmap_simulations <- Reduce(full_join, simmap_simulations)
      rownames(simmap_simulations) <- simmap_simulations$node
      simmap_simulations <- simmap_simulations[,!(names(simmap_simulations) %in% c("node"))]
      
      if (model_type == "HR") {
        # Make different for rates vs states
        simmap_simmulations_rates <- simmap_simulations
        
        # Replace numbers with states
        simmap_simulations[simmap_simulations == 1] <- "diurnal"
        simmap_simulations[simmap_simulations == 2] <- "nocturnal"
        simmap_simulations[simmap_simulations == 3] <- "diurnal"
        simmap_simulations[simmap_simulations == 4] <- "nocturnal"
        
        # Replace with rate names
        simmap_simmulations_rates[simmap_simmulations_rates == "diurnal"] <- "R1"
        simmap_simmulations_rates[simmap_simmulations_rates == "nocturnal"] <- "R1"
        
        simmap_simmulations_rates[simmap_simmulations_rates == 1] <- "R1"
        simmap_simmulations_rates[simmap_simmulations_rates == 2] <- "R1"
        simmap_simmulations_rates[simmap_simmulations_rates == 3] <- "R2"
        simmap_simmulations_rates[simmap_simmulations_rates == 4] <- "R2"
        
      } else {
        simmap_simulations[simmap_simulations == 1] <- "diurnal"
        simmap_simulations[simmap_simulations == 2] <- "nocturnal"
      }
      
      
      ## Compare the output of simmaps to the marginal ancestral reconstruction, they should correlate well
      
      data <- simmap_simulations
      data[data == "diurnal"] <- 1
      data[data == "nocturnal"] <- 0
      data <- data %>% mutate_if(is.character,as.numeric)
      simmap_avg <- data.frame(names = row.names(data), lik.anc = rowMeans(data))
      
      lik.anc <- anc_states$lik.anc
      lik.anc$diurnal <- lik.anc$diurnal_R1+lik.anc$diurnal_R2
      lik.anc$simmap_avg_names <- simmap_avg$names[match(row.names(lik.anc), simmap_avg$names)]
      lik.anc$simmap_avg <- simmap_avg$lik.anc[match(row.names(lik.anc), simmap_avg$names)]
      lik.anc$diff <- abs(lik.anc$diurnal - lik.anc$simmap_avg)
      
      comp_anc_methods <- ggplot(lik.anc[], aes(x = diurnal, y = simmap_avg)) + geom_point()
      pdf(file = paste(dataset_variable, "phylogeny_diel_plot_AncRec_SIMMAPvsSimulation", name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simmap_simulations), "zoom.pdf", sep = "_"), width = 10, height = 10)
      print(comp_anc_methods)
      dev.off()
      
      ## Need to make sure the tips are in the correct order, then phytools::getStates matches what I was doing
      simmap_simulations <- simmap_simulations[row.names(lik.anc),]
      
      ### Calculate transitions and then cummulative transitions through time
      simulated_transitions_simmap <- calculateSimulatedTransitions(simulated_data = simmap_simulations, phylo_tree = trpy_n)
      
      ## Function that returns cumsums, rows are now ordered by node.age, whether or it is included as a column in the output
      cumsums_simmap <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_simmap, phylo_tree = trpy_n, include_node_age = TRUE)
      
      if (model_type == "HR") {
        simulated_transitions_simmap_rates <- calculateSimulatedTransitions(simulated_data = simmap_simmulations_rates, phylo_tree = trpy_n)
        cumsums_simmap_rates <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_simmap_rates, phylo_tree = trpy_n, include_node_age = TRUE)
      }
      
      ###############################################################################################################################################
      ### Plot the results of the simulation, along with the actual data ### 
      ###############################################################################################################################################
      
      ## Plot the results of the simulation, either the summary (mean +/- stdev), the simulations, or both ('overlay')
      plot_simmap <- simulatedSwitchRatio(simulated_cumsums = cumsums_simmap, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "blue") #+ ylim(c(0,0.15))
      
      ## Can I add the computed value? Actual transitions
      ## Easy way is to make the plot, then use the data from it to add to the above
      switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02)
      
      plot <- plot_simmap + geom_line(data = switch.ratio$data, aes(x=node.age,y=ratio), colour = "black")
      plot <- plot + ggtitle(paste(model_type, "model SIMMAP,", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", ncol(simmap_simulations), ", Avg. transitions:", mean(colSums(simulated_transitions_simmap)), "+/-", round(sd(colSums(simulated_transitions_simmap)),1), sep = " "))
      
      # This 'zooms-in', which preserves the ribbon (setting limits removes part of the ribbon stdev)
      recon_plot <- plot + coord_cartesian(ylim = c(0, layer_scales(plot)$y$range$range[[2]]), xlim = abs(layer_scales(plot)$x$range$range)) 
      
      # Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
      geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot)$x$range$range)) 
      
      ## Save out the plot
      setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/outs/Figures")
      
      pdf(file = paste(dataset_variable, "phylogeny_diel_plot_transitions_SIMMAP", name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simmap_simulations), "zoom.pdf", sep = "_"), width = 10, height = 10)
      print(recon_plot / geo_scale + plot_layout(nrow = 2, heights = c(5,0.5)))
      dev.off()
      
      if (model_type == "HR") {
        plot_simmap_rates <- simulatedSwitchRatio(simulated_cumsums = cumsums_simmap_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "blue") #+ ylim(c(0,0.15))
        switch.ratio.rates <- switchRatio(ancestral_states = anc_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02)
        plot.rates <- plot_simmap_rates + geom_line(data = switch.ratio.rates$data, aes(x=node.age,y=ratio), colour = "black")
        plot.rates <- plot.rates + ggtitle(paste(model_type, "model SIMMAP (rates),", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", ncol(simmap_simulations), ", Avg. transitions:", mean(colSums(simulated_transitions_simmap_rates)), "+/-", round(sd(colSums(simulated_transitions_simmap_rates)),1), sep = " "))
        recon_plot_rates <- plot.rates + coord_cartesian(ylim = c(0, layer_scales(plot.rates)$y$range$range[[2]]), xlim = abs(layer_scales(plot.rates)$x$range$range)) 
        geo_scale_rates <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot.rates)$x$range$range))
        
        pdf(file = paste(dataset_variable, "phylogeny_diel_plot_transitions_SIMMAP_rates", name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simmap_simulations), "zoom.pdf", sep = "_"), width = 10, height = 10)
        print(recon_plot_rates / geo_scale_rates + plot_layout(nrow = 2, heights = c(5,0.5)))
        dev.off()
      }
      
    }
  }
}











# # Which tree?
# name_variable <- "all" # all, only_highqual, only_cartilaginous, or only_ingroup
# dataset_variable <- "fish" # fish or AllGroups
# 
# ## Load in the tree
# trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
# trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)
# 
# models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
# model <- models$HMM_2state_2rate_marg
# 
# 
# ### Run simmap
# 
# trait.data_n$diel_numb <- ifelse(trait.data_n$diel1 == "diurnal", 0, 1)
# 
# simmaps <- corHMM::makeSimmap(tree = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = model$solution, nSim = 5, nCores = 7)
# 
# # data.vector <- trait.data_n$diel1
# # names(data.vector) <- rownames(trait.data_n)
# # simmaps <- make.simmap(tree = trpy_n, x = data.vector, model = "ARD", nsim = 1)
# 
# 
# # # This plots one of them
# # plotSimmap(simmaps[[1]], node.numbers = TRUE)
# 
# ## I can use the first value of maps and match it with the node from edges
# ## Or, maybe use descendant nodes (which aren't duplicated, and includes tips)
# ## for each row of edge, match the last value of maps with the second node value
# ## Ahh, this doesn't include the root, because it is not a descendent
# 
# 
# ## This is weird though, the tip states should be constant across simulations, right?
# ## I think it doesn't matter, because I dont' plot them below.
# ## I don't know, but it somehow does not seem to matter, the make.simmap function always returns plots where the tips are the same
# ## But the distribution of transitions is not even close to the reconstructed states????
# 
# 
# simmap_simulations <- lapply(seq_along(simmaps), function(x) {
#   
#   simtree <- simmaps[[x]]
#   tip_states <- data.frame(edge = 1:nrow(simtree$edge), tips = simtree$edge[,2])
#   tip_states$state <- unlist(lapply(tip_states$edge, function(x) names(simtree$maps[x][[1]])[length(simtree$maps[x])]))
#   # Add root
#   root_edge <- grep(Ntip(trpy_n)+1, simtree$edge[1,])
#   tip_states <- rbind(tip_states, data.frame(edge = NA, tips = Ntip(trpy_n)+1, state = names(simtree$maps[root_edge][[1]])[1]))
#   
#   # These are wrong for the tips (I think it's the last simulated data, not necessarily what is at the tip)
#   # So replace tip values with actual tip values (this is what the plots do)
#   tip_states$state[match(1:Ntip(trpy_n), tip_states$tips)] <- trait.data_n$diel1
#   
#   ## I need to convert tips to tip_names and node_x
#   tip_states <- tip_states[order(tip_states$tips),]
#   tip_states$node <- c(trpy_n$tip.label, paste("Node", (Ntip(trpy_n)+1):(Ntip(trpy_n)+Nnode(trpy_n)), sep = ""))
#   
#   ## Convert to character names
#   # tip_states$state <- ifelse(tip_states$state == 1, "diurnal", "nocturnal")
#   
#   colnames(tip_states) <- c("edge", "tips", paste("state", x, sep = "_"), "node")
#   
#   return(tip_states[,c(4,3)])
#   
# })
# 
# # Then full_join them together!
# simmap_simulations <- Reduce(full_join, simmap_simulations)
# rownames(simmap_simulations) <- simmap_simulations$node
# simmap_simulations <- simmap_simulations[,!(names(simmap_simulations) %in% c("node"))]
# 
# simmap_simmulations_rates <- simmap_simulations
# # Replace numbers with states
# 
# simmap_simulations[simmap_simulations == 1] <- "diurnal"
# simmap_simulations[simmap_simulations == 2] <- "nocturnal"
# simmap_simulations[simmap_simulations == 3] <- "diurnal"
# simmap_simulations[simmap_simulations == 4] <- "nocturnal"
# 
# 
# simmap_simmulations_rates[simmap_simmulations_rates == "diurnal"] <- "R1"
# simmap_simmulations_rates[simmap_simmulations_rates == "nocturnal"] <- "R1"
# 
# simmap_simmulations_rates[simmap_simmulations_rates == 1] <- "R1"
# simmap_simmulations_rates[simmap_simmulations_rates == 3] <- "R2"
# simmap_simmulations_rates[simmap_simmulations_rates == 2] <- "R1"
# simmap_simmulations_rates[simmap_simmulations_rates == 4] <- "R2"
# 
# 
# simulated_transitions_simmap <- calculateSimulatedTransitions(simulated_data = simmap_simulations, phylo_tree = trpy_n)
# simulated_transitions_simmap_rates <- calculateSimulatedTransitions(simulated_data = simmap_simmulations_rates, phylo_tree = trpy_n)
# 
# ## Function that returns cumsums, rows are now ordered by node.age, whether or it is included as a column in the output
# cumsums_simmap <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_simmap, phylo_tree = trpy_n, include_node_age = TRUE)
# cumsums_simmap_rates <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_simmap_rates, phylo_tree = trpy_n, include_node_age = TRUE)
# 
# ###############################################################################################################################################
# ### Plot the results of the simulation, along with the actual data ### 
# ###############################################################################################################################################
# 
# anc_states <- readRDS(file = paste(dataset_variable, "diel_ancestral_states", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
# anc_rates <- readRDS(file = paste(dataset_variable, "diel_ancestral_rates", name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
# 
# 
# ## Plot the results of the simulation, either the summary (mean +/- stdev), the simulations, or both ('overlay')
# plot_simmap <- simulatedSwitchRatio(simulated_cumsums = cumsums_simmap, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "black") #+ ylim(c(0,0.15))
# plot_simmap_rates <- simulatedSwitchRatio(simulated_cumsums = cumsums_simmap_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "black") #+ ylim(c(0,0.15))
# 
# ## Can I add the computed value? Actual transitions
# ## Easy way is to make the plot, then use the data from it to add to the above
# switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02)
# switch.ratio.rates <- switchRatio(ancestral_states = anc_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02)
# 
# 
# 
# plot_simmap_3 <- plot_simmap_2 + geom_line(data = switch.ratio$data, aes(x=node.age,y=ratio), colour = "red")
# 
# switch.ratio + plot_simmap_2 + plot_layout(ncol = 1)
# 
# 
# ################################################################################################################
# ################################### Other non-used ways to access the edges ####################################
# ################################################################################################################
# 
# ## Both of these next two methods try to access the edges to determine both when the transitions occur (phytools method)
# ## or the whether one occurs on an edge (from paper). Phytools method extracts info from the plot, but could also be used
# ## to extract info from the node ages. However, both methods struggle with the underlying effect that the fast transition rates
# ## go back and forth a lot, and it looks messed up (way more transitions eg. 1343434343434343434342). I would have to come up with some
# ## solution to how to treat them (fast rates as one type (neither diurnal or nocturnal)).
# ## I think the most reasonable method will be to just ID the node states from the simmaps, and use those. This is also then the same as
# ## when I simulate the data using the other method.
# 
# 
# ### This is a solution from phytools guy
# ### It works for simple cases, but not the back-and-forths I see
# ### I will have to first transform things into each other (both fast rate states as 1)
# ## Then run something like the following
# 
# tree <- simmaps[[5]]
# 
# states <- sort(unique(getStates(tree)))
# 
# colors <- setNames(palette()[1:length(states)],states)
# 
# obj <- get("last_plot.phylo", envir=.PlotPhyloEnv)
# 
# nc <- sapply(tree$maps,length)-1 # This is whenever there is more than 1 state per edge
# 
# ii <- which(nc>0) # which edges have transitions
# 
# nc <- nc[ii] # the number of transitions (if one, then 2 states, j and j+1)
# 
# h <- vector()
# 
# for(i in 1:length(ii)){
#   for(j in 1:nc[i]){
#     ss <- names(tree$maps[[ii[i]]])[j+1] # this returns the second state
#     mm <- tree$edge[ii[i],1] # This returns the parental node for the edge
#     dd <- tree$edge[ii[i],2] # descendant node
#     
#     x <- rep(obj$xx[mm]+cumsum(tree$maps[[ii[i]]])[j],2) # This is the x position of the transition (the age)
#     # Instead of querying the plot ('obj$xx[mm]'), I can just query the node ages
#     
#     y<-c(obj$yy[dd]-0.5*mean(strheight(LETTERS)*cex),
#          obj$yy[dd]+0.5*mean(strheight(LETTERS)*cex)) # I suppose this is the y position
#     
#     lines(x,y,lwd=lwd,col=colors[ss],lend=2)
#     
#     h<-c(h,x[1])
#     
#   }
# }
# invisible(h)
# 
# 
# 
# #### This is stuff from the paper https://www.nature.com/articles/s41467-023-36501-4.epdf?sharing_token=QfMaSEdzTqF5T6YrXpW1pNRgN0jAjWel9jnR3ZoTv0OUoPP6nlC3aeT4B3x29l9M9UoeMzYPnHHHYj8ZAMisGmVLEM_TKDrD1V848azrWE5wnd64ljbzUMAuqmWIDS5MneI61VcU-gi5yqmh72WY59Y-WFX82AP_OoeNsEIEQyY%3D
# 
# ## https://github.com/stfriedman/Depth-transitions-paper/blob/main/Rcode/02_transition_analysis.R
# 
# ### What this does is take the 'mapped_edge' slot of the simmap, which is the sum of the time spent in each state for each edge, and ask how many times
# ### for each descendant node of an mrca, there is more than 1 state on that edge. If more than 1 state, that means there was a transition on that edge
# ### This should work, but could fail if you expect transitions from A -> B -> C on a branch, as this would only be counted as 1 transition.
# ### At this point, it avoids double counting nodes, by only searching for each node in the "node_1" column
# 
# ### So this works for my data, and I could adapt it so that it combines X1 & X3, and X2 & X4, to be clear it is correct  for example, only X1 and X3 might appear 
# ### on an edge, and be counted as a transition when it actually isn't). Easy to modify this to count transitions between rates as well
# ### However, it does not know when the transition took place, just that it took place on that branch (nothing in 'mapped_edge' lets you know which state was first, etc)
# ### For this, I would need to access the maps slot, which has the ordered list of states
# ### I could use the same framework, but then instead of looking at the 'mapped_edges', I use that to index 'maps', then the strings of fast rate transitions can be 'averaged'
# ### and I can take the midpoint + the age of node_1 as the transition time.
# 
# 
# ### grep the states in the 'maps' vector, take the mean(), and use that to order the df?
# 
# ## These are some functions, they use purrr, not lapply
# 
# num_trans <- function(mrca, simtree) {
#   trans_count <- function(desc) {
#     edges <- get_mapped_edge(simtree) %>%
#       mutate_at(vars(node_1:node_2), as.numeric)
#     edges <- edges %>% mutate(index = row_number(), Di = X1, Noc = X2)
# 
#     
#     df <- edges %>%
#             filter(node_1 %in% desc) %>%
#             select(-c(X1, X2)) %>%
#             pivot_longer(-c(node_1, node_2, index), names_to = "state", values_to = "value") %>%
#             filter(value != 0) %>%
#             group_by(node_1, node_2) %>%
#             add_count() %>%
#             select(node_1, node_2, value, state, index, n) %>%
#             filter(n != 1)
#     
#     ## At this point, I know that there was a transition, but not when (9.14 after 153, or 16.1 after 153?)
#     ## If there is a transition, I can just ask if the first and last states differ,
#     ## if so, then the time is node + the time spent in the original state. If they are the same, but this says there is a transition
#     ## then there were two transitions?? How to deal with this?
#     ## in this case, it starts at 4, ends as 4, but spends 9 million years as diurnal 3?
#     
#     grepl("1|3", names(simtree$maps[unique(df$index)][[1]])[1]) # not diurnal, therefore starts as noc
#     
#     
#   }
#   tibble(mrca) %>%
#     mutate(
#       desc = map(mrca, ~ c(.x, getDescendants(simtree, .x))),
#       trans_num = map_dbl(desc, trans_count)
#     )
# }
# 
# get_mapped_edge <- function(x) {
#   data.frame(x$mapped.edge) %>%
#     rownames_to_column("edges") %>%
#     as_tibble() %>%
#     separate(edges, c("node_1", "node_2"), sep = ",")
# }
# 
# 
# ## This is the actual command that returns the averaged transitions
# 
# clade_transitions <- tibble(tree_id = 1:length(simmaps)) %>%
#   mutate(simmap = simmaps,
#          trans = map(simmaps, ~num_trans(mrcas$mrca, .x))) %>%
#   unnest(trans) %>%
#   select(-simmap, -desc) %>%
#   group_by(mrca) %>%
#   summarize(trans_num = median(trans_num)) 
# 
# 
# 
# 
# 
# ## OK, so 1 and 2 are the fast rate transitions
# 
# test <- simmaps[[1]]$maps[6][[1]]
# 
# 
# test2 <- c(test)
# 
# 
# 
# 
# 
# 
# ## keep
# num_trans <- function(mrca, simtree) {
#   trans_count <- function(desc) {
#     edges <- get_mapped_edge(simtree) %>%
#       mutate_at(vars(node_1:node_2), as.numeric)
#     
#     edges %>%
#       filter(node_1 %in% desc) %>%
#       pivot_longer(-c(node_1, node_2), names_to = "state", values_to = "value") %>%
#       filter(value != 0) %>%
#       group_by(node_1, node_2) %>%
#       add_count() %>%
#       filter(n > 1) %>%
#       select(node_1, node_2) %>%
#       unique() %>%
#       nrow(.)
#   }
#   tibble(mrca) %>%
#     mutate(
#       desc = map(mrca, ~ c(.x, getDescendants(simtree, .x))),
#       trans_num = map_dbl(desc, trans_count)
#     )
# }
