### FUNCTIONS FOR LINEAGE THROUGH TIME PLOTS
### Associated with Fish_sleep project
### Copyright Maxwell Shafer 2022
### The following functions are associated with scripts for analysing diurnal and nocturnal states through time across a phylogeny
### and can be used with both data from fish, but also data from all vertebrates
### Updated 07.12.2022

## Source the script that allows for new scales on a ggplot
source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")

## Source the script with the geom_gg scale
source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/gggeo_scale.R")

### Need to make these work for if there is another trait involved (for example, marine/fresh). Would still be useful to plot the Di/Noc for these models

### This function loads either the tree or the trait data (only the tree and the model are required for the ancestral reconstruction)
loadTree <- function(return = "tree", dataset = c("fish", "AllGroups", "tetrapods", "mammals"), subset = c("no", "all", "only_highqual", "only_ingroup", "only_cartilaginous", "not_mammals", "custom"), custom_tips = NA) {
  require(ape)
  
  ## Load tree and trait data
  tr.calibrated <- readRDS(paste("tr_tree_calibrated_", dataset, ".rds", sep = ""))
  trait.data <- readRDS(paste("trait_data_", dataset, ".rds", sep = ""))
  
  ## Subset
  
  if (subset %in% c("no", "all")) {
    name_variable <- "all"
  }
  
  if (subset == "only_highqual") {
    trait.data <- trait.data[trait.data$confidence > 1,]
    name_variable <- "only_highqual"
  }
  
  if (subset == "only_ingroup") {
    node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Lepidosiren_paradoxa", "Lutjanus_fulvus"))
    tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
    trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
    name_variable <- "only_ingroup"
  }
  
  if (subset == "only_cartilaginous") {
    node_of_interest <- getMRCA(tr.calibrated, tip = c("Rhizoprionodon_terraenovae", "Rhynchobatus_djiddensis"))
    tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
    trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
    name_variable <- "only_cartilaginous"
  }
  
  if (subset == "not_mammals") {
    trait.data <- trait.data[!(trait.data$group %in% "Mammalia"),]
    tr.calibrated <- keep.tip(phy = tr.calibrated, tip = trait.data$unique_name)
    name_variable <- "not_mammals"
  }
  
  # Just diurnal/nocturnal
  trait.vector_n <- trait.data$diel1
  names(trait.vector_n) <- trait.data$species
  trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
  trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
  trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
  trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]
  
  if (subset == "custom") {
    if (is.na(custom_tips) | length(custom_tips) != 2) {
      stop("Custom subsets require 2 tip taxa names (tip labels)")
    }
    node_of_interest <- getMRCA(tr.calibrated, tip = custom_tips)
    tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
    trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
    name_variable <- "custom"
  }
  
  
  ## I need either just the phylo_tree (trpy_n) and/or the trait.data_n
  if (return == "tree") {
    return(trpy_n)
  }
  if (return == "trait_data") {
    return(trait.data_n)
  }
  if (return == "both") {
    return(list(trpy_n, trait_data_n))
  }
}


### Write a function to extract ancestral likelihoods from a model
returnAncestralStates <- function(phylo_model = model, phylo_tree = trpy_n, rate.cat = FALSE) {
  
  # Create a data frame with trait values from reconstruction
  lik.anc <- as.data.frame(rbind(phylo_model$tip.states, phylo_model$states))
  row.names(lik.anc) <- c(row.names(lik.anc)[1:length(phylo_tree$tip.label)], (Ntip(phylo_tree) + 1):(Ntip(phylo_tree) + Nnode(phylo_tree)))
  
  if (phylo_model$rate.cat == 1) {
    print("model has single rate category")
    states <- unique(phylo_model$data[,2])
    names(states) <- as.numeric(unique(phylo_model$data.legend[,2]))
    rate_states <- sort(states)
    states <- sort(states)
  }
  if (phylo_model$rate.cat == 2) {
    print("model has 2 rate categories")
    states <- unique(phylo_model$data[,2])
    names(states) <- as.numeric(unique(phylo_model$data.legend[,2]))
    states <- sort(states)
    rate_states <- c(paste(states, "R1", sep = "_"), paste(states, "R2", sep = "_"))
  }
  if (phylo_model$rate.cat == 3) {
    print("model has 3 rate categories")
    states <- unique(phylo_model$data[,2])
    names(states) <- as.numeric(unique(phylo_model$data.legend[,2]))
    states <- sort(states)
    rate_states <- c(paste(states, "R1", sep = "_"), paste(states, "R2", sep = "_"), paste(states, "R3", sep = "_"))
  }
  
  colnames(lik.anc) <- rate_states
  
  ancestral_states <- list()
  ancestral_states$lik.anc <- lik.anc
  ancestral_states$node <- 1:(length(phylo_tree$tip.label) + phylo_tree$Nnode) # The only thing that isn't correct is the row.names on lik.anc, but I don't use those!
  ancestral_states$states <- states
  ancestral_states$rate_states <- rate_states
  print(paste("returning ancestral states for", phylo_tree$Nnode, "internal nodes corresponding to", length(phylo_tree$tip.label), "tips in tree provided", sep = " "))
  return(ancestral_states)
}

### Function for running simulations

simulateCustom <- function(phylo_tree = trpy_n, models_list = models, model_type = "ER", rates = c(0.1,0.1), states = c("nocturnal", "diurnal"), simulation_numb = 100) {
  require(ape)
  
  if (model_type == "ER") {
    if (length(rates) != 1) {
      stop("ER models must have a single rate")
    }
    # Use the ER model (actually the ARD model but with symmetric rates, b/c I'm not sure how to code in ER model)
    simulation <- replicate(simulation_numb, rTraitDisc(phy = phylo_tree, model = "SYM", k = 2, rate = rates, states = states, ancestor = T))
  }
  
  if (model_type == "ARD") {
    if (length(rates) != length(states)*(length(states)-1)) {
      stop("ARD models must have rates equal to n*(n-1) states")
    }
    ## Simulate data based on the ARD model
    simulation <- replicate(simulation_numb, rTraitDisc(phylo_tree, model = "ARD", k = 2, rate = rates, states = states, ancestor = T))
  }
  
  if (model_type == "HR") {
    # Use the best fit model
    if (!(is.matrix(rates))) {
      stop("for HR models, rates must be a solution matrix")
    }
    rates[is.na(rates)] <- 0
    
    simulation <- replicate(simulation_numb, rTraitDisc(trpy_n, model = rates, states = states, ancestor = T, root.value = 4))
  }
  
  return(simulation)
}


### Function to generate object from simulated data

calculateSimulatedTransitions <- function(simulated_data = simulation, phylo_tree = trpy_n) {
  # Determine ages of each node
  node_heights <- nodeHeights(phylo_tree)
  node <- c(1:(Ntip(phylo_tree) + Nnode(phylo_tree)))
  
  node.age <- node_heights[match(node, phylo_tree$edge[,2]),2]
  node.age[node[is.na(node.age)]] <- 0
  node.age <- (node.age - max(node.age))*-1
  
  simulated_data[simulated_data == "diurnal"] <- 1  
  simulated_data[simulated_data == "nocturnal"] <- 0 
  
  ## I want to have this run it for all simulations, and return a df the same as simulation
  # ID parental nodes and parental states
  parental.node <- unlist(lapply(node, function(x) Ancestors(phylo_tree, x, type = "parent")))
  
  # This now becomes a df with each column representing the parent.diel for each simulation
  parent.diel <- apply(simulated_data, 2, function(x) unlist(lapply(parental.node, function(y) x[match(y, node)])) )
  
  # parent of the root is NA
  parent.diel[is.na(parent.diel)] <- 0
  
  # calculated whether each is a transition, Reduce + cbind it
  transition <- Reduce(cbind, lapply(seq_along(1:ncol(simulated_data)), function(x) ifelse(ifelse(parent.diel[,x] > 0.5, 1, 0) != ifelse(simulated_data[,x] > 0.5, 1, 0), 1, 0) ))
  
  rownames(transition) <- node
  colnames(transition) <- c(paste("simulation", 1:ncol(transition), sep = "_"))

  print(paste("Identified an average of", mean(colSums(transition)), "transitions across", ncol(simulated_data), "simulations", sep = " "))
  
  return(transition)
}


### Function to determine the number of transitions between states
calculateStateTransitions <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, ancestor = 0, rate.cat = FALSE) {
  
  states <- ancestral_states$states
  rate_states <- ancestral_states$rate_states
  
  # Determine the number of transitions or diurnal/nocturnal taxa by age
  if (length(ancestral_states$states) == length(ancestral_states$rate_states)){
    ancestral_states$recon_states <- ancestral_states$lik.anc[,states[[1]]]
  } else {
    ancestral_states$recon_states <- as.numeric(rowSums(ancestral_states$lik.anc[,rate_states[grep(states[[1]], rate_states)]]))
  }
  
  ## Add functionality to do this for rate categories?
  if (rate.cat) {
    print("returning rate categories instead of states")
    ancestral_states$recon_states <- as.numeric(rowSums(ancestral_states$lik.anc[,rate_states[grep("R1", rate_states)]]))
  }
  
  # Determine ages of each node
  node_heights <- nodeHeights(phylo_tree)
  
  ancestral_states$node.age <- node_heights[match(ancestral_states$node, phylo_tree$edge[,2]),2] # This extracts all but one that is missing from phylo_tree$edge[,2], and generates an NA
  ancestral_states$node.age[ancestral_states$node[is.na(ancestral_states$node.age)]] <- 0 # This is the root
  # ancestral_states$node.age[is.na(ancestral_states$node.age)] <- 0
  ancestral_states$node.age <- (ancestral_states$node.age - max(ancestral_states$node.age))*-1
  
  # ancestral_states <- ancestral_states[order(ancestral_states$Time),]
  
  # ID parental nodes and parental states
  ancestral_states$parental.node <- unlist(lapply(ancestral_states$node, function(x) Ancestors(phylo_tree, x, type = "parent")))
  ancestral_states$parent.diel <- unlist(lapply(ancestral_states$parental.node, function(x) ancestral_states$recon_states[match(x, ancestral_states$node)]))
  # parent of the root is NA
  ancestral_states$parent.diel[is.na(ancestral_states$parent.diel)] <- ancestor
  
  ancestral_states$transition <- ifelse(ifelse(ancestral_states$parent.diel > 0.5, 1, 0) != ifelse(ancestral_states$recon_states > 0.5, 1, 0), 1, 0)
  # For those with transitions, I can just ask what they are, and that's the switch type!
  ancestral_states$trans.ND <- ifelse(ancestral_states$transition == 1, ifelse(ancestral_states$recon_states > 0.5, 1, 0),0)
  ancestral_states$trans.DN <- ifelse(ancestral_states$transition == 1, ifelse(ancestral_states$recon_states < 0.5, 1, 0),0)
  # ancestral_states$transition[is.na(ancestral_states$transition)] <- 0
  
  print(paste("Identified", table(ancestral_states$transition)[2], "transitions between", ancestral_states$states[1], "and", ancestral_states$states[2], sep = " "))
  return(ancestral_states)
}


### Function to calculate transition history on tree for the state
calculateLinTransHist <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n) {
  
  ancestors <- lapply(c(1:length(phylo_tree$tip.label)), function(x) Ancestors(phylo_tree, x, type = "all"))
  ancestors <- lapply(seq_along(ancestors), function(x) append(c(1:length(phylo_tree$tip.label))[[x]], ancestors[[x]]))
  print("calculating ancestral states for all nodes")
  recon_states <- ifelse(ancestral_states$recon_states > 0.5, 1, 0)
  ancestors.diel <- lapply(ancestors, function(x) lapply(x, function(y) recon_states[match(y, ancestral_states$node)]))
  
  print("identifying switch types")
  # This works, can I simplify it so it works on a vector?
  switch.type <- unlist(lapply(ancestors.diel, function(x) {
    df <- data.frame(test = unlist(x))
    df <- df[with(df, c(test[-1]!= test[-nrow(df)], TRUE)),]
    return(paste(df, collapse = ""))
  }))
  
  ancestral_states$switch.type <- as.character(switch.type)
  
  print("Identified the following switch types")
  print(table(switch.type))
  
  return(ancestral_states)
}


### Function to calculate the cummulative sums (of switches and switch types)
returnCumSums <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, use.height = TRUE) {
  # This is all nodes (for transitions only)
  node_order <- order(ancestral_states$node.age, decreasing = T) # Oldest node first (root)
  
  ancestral_states$transition_cumsum <- cumsum(ancestral_states$transition[node_order])
  ancestral_states$trans.DN_cumsum <- cumsum(ancestral_states$trans.DN[node_order])
  ancestral_states$trans.ND_cumsum <- cumsum(ancestral_states$trans.ND[node_order])
  
  # This now no longer works, because it is the age of the branch points?
  
  ## OK, for this graph, I want to index the first half of node_heights!
  
  # First column from 'edge' is the higher node, second is lower
  # both columns from 'node_heights' should match by the index ('edge')
  # If I find the tip in edge[,2], but pull node_heights[,1], I get the age of the node connected to the tip (the branch point of the tip)

  node_heights <- nodeHeights(phylo_tree)
  node.age2 <- node_heights[match(c(1:length(phylo_tree$tip.label)), phylo_tree$edge[,2]),1] # This finds the position of each node # in the 2nd column of 'edge', and returns the first column of 'edge', which is the age of the split

  node.age2[is.na(node.age2)] <- 0
  node.age2 <- (node.age2 - max(node.age2))*-1
  
  cumsums <- data.frame(row.names = phylo_tree$tip.label[order(node.age2, decreasing = T)])
  cumsums$node.age <- node.age2[order(node.age2, decreasing = T)]
  # cumsums$node.age <- ancestral_states$node.age[node_order_2]
  cumsums$Lineage_Cumsum <- cumsum(rep(1, length(phylo_tree$tip.label)))
  
  for (i in unique(ancestral_states$switch.type)) {
    cumsums[,i] <- cumsum(ifelse(ancestral_states$switch.type[order(node.age2, decreasing = T)] == i, 1, 0))
  }
  
  ancestral_states$cumsums <- cumsums
  return(ancestral_states)
}


### Function to make the Switch history histogram

switchHisto <- function(ancestral_states = ancestral_states, replace_variable_names = TRUE, backfill = TRUE) {
  require(tidyr)
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  require(RColorBrewer)
  
  ## Prepare the data by gathering
  node.data <- gather(ancestral_states$cumsums, "variable", "Cummulative_ratio", -node.age, -Lineage_Cumsum)
  node.data <- node.data %>% group_by(node.age, variable) %>% summarise(n = sum(Cummulative_ratio)) %>% mutate(percentage = n / sum(n))
  
  ## Set factor levels for switch history types
  levels <- sort(unique(node.data$variable))
  
  # This is the correct order of the levels
  factor_order <- match(levels[order(nchar(levels), levels)], levels)
  
  node.data$variable <- factor(node.data$variable, levels = levels[factor_order])
  
  ## Determine the colour scales
  levels_str <- str_sub(levels, start = 1, end = 1)
  noc_levels <- length(levels_str[levels_str == "0"])
  di_levels <- length(levels_str[levels_str == "1"])
  
  # This doesn't work if there are only 2 levels each
  if (noc_levels < 3) {
    noc_colours <- brewer.pal(3, "Blues")[(3-length(noc_levels)):3]
  } else {
    noc_colours <- brewer.pal(noc_levels, "Blues") 
  }
  if (di_levels < 3) {
    di_colours <- brewer.pal(3, "Reds")[(3-length(di_levels)):3]
  } else {
    di_colours <- brewer.pal(di_levels, "Reds") 
  }
  colours <- c(noc_colours, di_colours)
  colours <- colours[factor_order]
  
  if(replace_variable_names) {
    variable_names <- levels
    variable_names <- gsub("0", "Noc", variable_names)
    variable_names <- gsub("1", "-Di-", variable_names)
    node.data$variable <- variable_names[match(node.data$variable, levels)]
    node.data$variable <- factor(node.data$variable, levels = variable_names[factor_order])
  }
  if (backfill) {
    node.data[nrow(node.data)+1,] <- list(max(ancestral_states$node.age), levels(node.data$variable)[1], 1, 1)
  }
  ## Should come up with function to make this? Yeah, and auto-find colors maybe from the red and blue scales
  plot <- ggplot(node.data, aes(x = node.age, y = percentage, fill = variable)) + theme_classic() + geom_area() + scale_x_reverse()
  plot <- plot + xlab("Millions of years ago") + ylab("Fraction of lineages") #+ scale_fill_manual(values = c("blue4", "red4", "blue3", "red3", "blue1", "red1"))
  

  
  plot <- plot + scale_fill_manual(values = colours)
  
  return(plot)
}


### Function to make the Switch ratio plot (line)

switchRatio <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, node.age.cutoff = 0.1, use_types = FALSE) {
  require(ape)
  
  node_order <- order(ancestral_states$node.age, decreasing = T)
  df <- data.frame(node = ancestral_states$node[node_order], node.age = ancestral_states$node.age[node_order], transition_cumsum = ancestral_states$transition_cumsum, trans.ND.cumsum = ancestral_states$trans.ND_cumsum, trans.DN.cumsum = ancestral_states$trans.DN_cumsum)
  
  ## calculate it the way the ape function does, but with my node.ages?
  ## This works close to the function, but neglects the node at the root
  df$ltt_source <- ifelse(df$node %in% c(1:Ntip(phylo_tree)), -1, 1)
  df$ltt_cumsum <- cumsum(df$ltt_source)

  
  ## OK this seems to be really close. I can calculate ltt's properly, but it also makes a weird thing were ltt decreases near present day
  df$ratio <- df$transition_cumsum/df$ltt_cumsum
  df$ratio_DN <- df$trans.DN.cumsum/df$ltt_cumsum
  df$ratio_ND <- df$trans.ND.cumsum/df$ltt_cumsum
  
  df$ratio[is.nan(df$ratio)] <- 0
  df$ratio[is.infinite(df$ratio)] <- 0
  
  if(use_types) {
    df2 <- gather(df, key = "trans_type", value = "ratio", -node, -node.age, -transition_cumsum, -trans.ND.cumsum, -trans.DN.cumsum, -ltt_source, -ltt_cumsum, -ratio)
    plot <- ggplot(df2[df2$node.age > node.age.cutoff,], aes(x = node.age, y = ratio, colour = trans_type)) + geom_line() + theme_classic() + scale_x_reverse() #+ scale_x_reverse(limits = c(413,0.1)) + ylim(c(0,0.15)) # ltt would be the number of times it was possible to switch? Or something like that?
    plot <- plot + xlab("Millions of years ago") + ylab("Fraction of lineages transitioning")
  } else {
    plot <- ggplot(df[df$node.age > node.age.cutoff,], aes(x = node.age, y = ratio)) + geom_line() + theme_classic() + scale_x_reverse() #+ scale_x_reverse(limits = c(413,0.1)) + ylim(c(0,0.15)) # ltt would be the number of times it was possible to switch? Or something like that?
    plot <- plot + xlab("Millions of years ago") + ylab("Fraction of lineages transitioning")
  }
  
  return(plot)
}

### This function returns cumsums for simulated transition data

calculateSimualtedTransitionCumsums <- function(simulated_transitions = simulated_transitions, phylo_tree = trpy_n, include_node_age = TRUE) {
  
  node_heights <- nodeHeights(phylo_tree)
  
  extracted_data <- as.data.frame(simulated_transitions)
  
  extracted_data$node.age <- node_heights[match(rownames(extracted_data), phylo_tree$edge[,2]),2]
  extracted_data$node.age[is.na(extracted_data$node.age)] <- 0 # This is the root (always node equal to Ntip(tree) + 1)
  extracted_data$node.age <- (extracted_data$node.age - max(extracted_data$node.age))*-1
  extracted_data <- extracted_data[order(extracted_data$node.age, decreasing = T),]
  
  extracted_data_2 <- as.data.frame(apply(extracted_data[,1:100], 2, function(x) cumsum(x)))
  if(include_node_age) {
    extracted_data_2$node.age <- extracted_data$node.age
  }
  
  return(extracted_data_2)
}


### Function to plot the lineages with X# of switches

switchTree <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, layout = "circular", replace_variable_names = TRUE) {
  
  possible.colors <- rev(c("grey75", "grey50", "black", "blue", "green", "red"))
  
  switch.types <- unique(ancestral_states$switch.type)
  switch.type <- c(as.character(ancestral_states$switch.type), rep("internal_node", Nnode(phylo_tree)))
  switch.type <- factor(switch.type, levels = c(switch.types, "internal_node"))

  colours <- possible.colors[1:length(switch.types)]
  colours <- c(rev(colours), possible.colors[6])
  
  if(replace_variable_names) {
    variable_names <- switch.types
    variable_names <- gsub("0", "Noc", variable_names)
    variable_names <- gsub("1", "-Di-", variable_names)
    variable_names <- c(variable_names, "internal_node")
    switch.type <- variable_names[match(c(ancestral_states$switch.type, rep("internal_node", Nnode(phylo_tree))), switch.types)]
    switch.type <- factor(switch.type, levels = variable_names)
  }
  
  tree <- ggtree(phylo_tree, layout = layout, size = 0.25) %<+% as.data.frame(switch.type) + aes(color = switch.type) + geom_tippoint(aes(color = switch.type), shape = 16) + scale_color_manual(values = colours)
  
  return(tree)
}

## Function that plots switch ratio

simulatedSwitchRatio <- function(simulated_cumsums = cumsums, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "red") {
  require(ape)
  
  # If there is a supplied node.age, use it to order, otherwise assume rows are already ordered by node.age
  if ('node.age' %in% colnames(simulated_cumsums)) {
    node_order <- order(simulated_cumsums$node.age, decreasing = T)
  } else {
    node_order <- rownames(simulated_cumsums)
  }
  
  simulated_cumsums <- simulated_cumsums[node_order,]
  
  df <- data.frame(node = rownames(simulated_cumsums)[node_order], node.age = simulated_cumsums$node.age[node_order])
  
  df$ltt_source <- ifelse(df$node %in% c(1:Ntip(phylo_tree)), -1, 1)
  df$ltt_cumsum <- cumsum(df$ltt_source)
  
  ## Calculate the ratio
  simulated_ratios <- apply(simulated_cumsums[,1:(ncol(simulated_cumsums)-1)], 2, function(x) x/df$ltt_cumsum )
  
  simulated_ratios[is.nan(simulated_ratios)] <- 0
  simulated_ratios[is.infinite(simulated_ratios)] <- 0
  
  simulated_ratios <- as.data.frame(simulated_ratios)
  
  
  if (plot_type == "summary") {
    df$mean <- rowMeans(simulated_ratios)
    df$stdev <- apply(simulated_ratios, 1, function(x) sd(x))
    
    plot <- ggplot(df[df$node.age > node.age.cutoff,], aes(x = node.age, y = mean)) + geom_line(colour = highlight_colour) + geom_ribbon(aes(ymin = mean-stdev, ymax = mean+stdev), fill = "blue", alpha = 0.1) + theme_classic() + scale_x_reverse() + theme(legend.position = "none")
    plot <- plot + xlab("Millions of years ago") + ylab("Fraction of lineages transitioning")
  }  
  
  if (plot_type == "simulation") {
    simulated_ratios$node.age <- df$node.age
    simulated_ratios_2 <- gather(as.data.frame(simulated_ratios), "simulation", "Ratio", -node.age)
    plot <- ggplot(simulated_ratios_2[simulated_ratios_2$node.age > node.age.cutoff,], aes(x = node.age, y = Ratio, color = simulation)) + geom_line() + theme_classic() + scale_x_reverse() + theme(legend.position = "none") #+ scale_x_reverse(limits = c(413,0.1)) + ylim(c(0,0.15)) # ltt would be the number of times it was possible to switch? Or something like that?
  }
  
  if (plot_type == "overlay") {
    simulated_ratios_2 <- simulated_ratios
    simulated_ratios_2$mean <- rowMeans(simulated_ratios)
    simulated_ratios_2$stdev <- apply(simulated_ratios, 1, function(x) sd(x))
    
    simulated_ratios_2$node.age <- df$node.age
    simulated_ratios_2$colour <- "colour"
    simulated_ratios_2 <- gather(as.data.frame(simulated_ratios_2), "simulation", "Ratio", -node.age, -mean, -stdev, -colour)
    
    # Plot the simulations
    # then plot on top with the subsetted df
    
    plot <- ggplot(simulated_ratios_2[simulated_ratios_2$node.age > node.age.cutoff,], aes(x = node.age, y = Ratio, group = simulation)) + geom_line(colour = "grey85")
    plot <- plot + new_scale("colour") + geom_line(data = simulated_ratios_2[simulated_ratios_2$node.age > node.age.cutoff & simulated_ratios_2$simulation == "simulation_1",], aes(x = node.age, y = mean), colour = highlight_colour) 
    plot <- plot + new_scale("alpha") + geom_ribbon(data = simulated_ratios_2[simulated_ratios_2$node.age > node.age.cutoff & simulated_ratios_2$simulation == "simulation_1",], aes(ymin = mean-stdev, ymax = mean+stdev), fill = "blue", alpha = 0.2) + theme_classic() + scale_x_reverse() + theme(legend.position = "none")
  }
  
  return(plot)
}


## Write a function that takes the rate index from computed models of trat evolution and make a plot which shows the rates

plotRateIndex <- function(model = model, HMM = FALSE, levels = NA) {
  require(tidyr)
  require(zoo)
  require(schoolmath)
  require(ggplot2)
  
  if(HMM == TRUE) {
    rates <- as.vector(model$solution)
    rates[is.na(rates)] <- 0
    
    matrix <- as.data.frame(x = model$index.mat, row.names = levels)
    colnames(matrix) <- levels
    matrix$levels <- levels
    
  } else {
    if (length(colnames(model$lik.anc)) != 4) {
      stop("model has more than or less than 4 categories")
    }
    rates <- model$rates
    matrix <- as.data.frame(x = model$index.matrix, 
                            row.names = colnames(model$lik.anc))
    colnames(matrix) <- colnames(model$lik.anc)
    
    matrix$levels <- factor(colnames(model$lik.anc), 
                            levels = colnames(model$lik.anc))
  }
  
  
  
  
  index.df <- data.frame(traits = colnames(matrix)[1:4], 
                         x = c(1,1,2,2), 
                         y = c(1,2,1,2))
  
  out <- Reduce(rbind, lapply(rownames(matrix), function(x) {
    if(grep(x, rownames(matrix)) == 1) {
      groups <- c(1,2,3,4)
    }
    if(grep(x, rownames(matrix)) == 2) {
      groups <- c(5,6,7,8)
    }
    if(grep(x, rownames(matrix)) == 3) {
      groups <- c(9,10,11,12)
    }
    if(grep(x, rownames(matrix)) == 4) {
      groups <- c(13,14,15,16)
    }
    df <- data.frame(from = rownames(matrix[x,1:4]), to = colnames(matrix)[1:4], value = rates[match(matrix[x,1:4], index(rates))])
    df1 <- cbind(df, data.frame(x = index.df$x[match(df$from, index.df$trait)], y = index.df$y[match(df$from, index.df$trait)], group = groups))
    #df2 <- data.frame(from = colnames(matrix)[1:4], to = rownames(matrix[x,1:4]), value = rates[match(matrix[x,1:4], index(rates))])
    df2 <- cbind(df, data.frame(x = index.df$x[match(df$to, index.df$trait)], y = index.df$y[match(df$to, index.df$trait)], group = groups))
    df <- rbind(df1, df2)
    df
  }))
  
  out$from <- factor(out$from, levels = matrix$levels)
  out$to <- factor(out$to, levels = matrix$levels)
  out$from_numb <- as.numeric(out$from)
  out$to_numb <- as.numeric(out$to)
  
  out$offset <- (is.odd(out$from_numb) & is.odd(out$to_numb)) |(is.even(out$from_numb) & is.even(out$to_numb))
  out$offset2 <- ifelse(as.numeric(out$from) > as.numeric(out$to), 
                        0.1, 
                        ifelse(as.numeric(out$from) < as.numeric(out$to), 
                               -0.1, 
                               0
                               )
                        )
  
  out$x2 <- ifelse(out$offset, out$x, out$x + out$offset2)
  out$y2 <- ifelse(out$offset, out$y + out$offset2, out$y)
  
  out$x3 <- ifelse(out$from_numb == 1 | out$to_numb == 1, out$x2 - 0.1, out$x2)
  out$y3 <- ifelse(out$from_numb == 1 | out$to_numb == 1, out$y2 - 0.1, out$y2)
  
  out$x3 <- ifelse(out$from_numb == 2 | out$to_numb == 2, out$x3 - 0.1, out$x3)
  out$y3 <- ifelse(out$from_numb == 2 | out$to_numb == 2, out$y3 + 0.1, out$y3)
  
  out$x3 <- ifelse(out$from_numb == 3 | out$to_numb == 3, out$x3 + 0.1, out$x3)
  out$y3 <- ifelse(out$from_numb == 3 | out$to_numb == 3, out$y3 - 0.1, out$y3)
  
  out$x3 <- ifelse(out$from_numb == 4 | out$to_numb == 4, out$x3 + 0.1, out$x3)
  out$y3 <- ifelse(out$from_numb == 4 | out$to_numb == 4, out$y3 + 0.1, out$y3)
  
  out$label <- NA
  label_df <- data.frame(from = rownames(matrix), to = rownames(matrix), value = NA, x = NA, y = NA, group = c(1,5,9,13), from_numb = NA, to_numb = NA, offset = NA, offset2 = NA, x2 = NA, y2 = NA, x3 = c(0.8, 0.8,2.2,2.2), y3 = c(0.8,2.2,0.8,2.2), label = rownames(matrix))
  out2 <- rbind(out,label_df)
  
  plot <- ggplot(out2, aes(x = x3, y = y3, group = group, size = value, colour = value, label = label)) + geom_point() + geom_path(arrow = arrow(angle = 30, length = unit(0.25, "inches"), ends = "last", type = "open")) + geom_label(colour = "black")
  plot <- plot + theme_classic() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank()) #+ ylim(0.5,2.5) + xlim(0.5,2.5) 
  plot <- plot + annotate('text', x = 1.2, y = 2.62, label = model$call, colour = "red", size = 6)
  plot <- plot + annotate('text', x = 1.2, y = 2.5, label = paste("Log-likelihood:", model$loglik, sep = " "), colour = "red", size = 6)
  return(plot)
}



gggeo_scale <- function(gg, fill = NULL, color = "black", alpha = 1, height = .05, size = 5, quat = FALSE, pos = "bottom", abbrv = TRUE, periods = NULL, neg = FALSE, blank.gg = FALSE) {
  #This is a function to add a geologic time scale to a ggplot object.
  #gg: the ggplot object
  #fill: the fill color of the boxes; default are the colors from the Commission for the Geological Map of the World (CGMW);
  #   custom fill colors can be provided and will be recycled if necessary
  #   if a custom dataset is provided with periods without color and without fill, a greyscale will be used
  #color: the outline color of the boxes
  #alpha: transparency of the fill colors
  #height: the proportional height of the plot to use for the scale 
  #size: the size of the text in the scale
  #quat: specifies whether the Quaternary should be labelled
  #pos: which side to add the scale to (left, right, top, or bottom)
  #abbrv: whether to use abbreviations instead of full period names
  #periods: a custom data set of time interval boundaries, with the following columns:
  #   period: name of each period (will be used as labels if no abbreviations are provided)
  #   max_age: the oldest boundary of each time interval
  #   min_age: the youngest boundary of each time interval
  #   abbr: (optional) abbreviations that will be used as labels
  #   color: (optional) a hex color code (which can be obtained with rgb()) for each time interval
  #neg: set this to true if your x-axis is actually negative values
  
  require(ggplot2)
  if(is.null(periods)){
    periods <- data.frame(period = c("Quaternary", "Neogene", "Paleogene", "Cretaceous", "Jurassic", "Triassic", "Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician", "Cambrian", "Ediacaran", "Cryogenian", "Tonian"),
                          max_age = c(2.588, 23.03, 66, 145, 201.3, 252.2, 298.9, 358.9, 419.2, 443.4, 485.4, 541, 635, 720, 1000),
                          min_age = c(0, 2.588, 23.03, 66, 145, 201.3, 252.2, 298.9, 358.9, 419.2, 443.4, 485.4, 541, 635, 720),
                          abbr = c("Q", "N", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Cm","E","Cr","To"),
                          color = c(rgb(249, 249, 127, maxColorValue = 255),rgb(255, 230, 25, maxColorValue = 255),rgb(253, 154, 82, maxColorValue = 255),rgb(127, 198, 78, maxColorValue = 255),rgb(52, 178, 201, maxColorValue = 255),rgb(129, 43, 146, maxColorValue = 255),rgb(240, 64, 40, maxColorValue = 255),rgb(103, 165, 153, maxColorValue = 255),rgb(203, 140, 55, maxColorValue = 255),rgb(179, 225, 182, maxColorValue = 255),rgb(0, 146, 112, maxColorValue = 255),rgb(127, 160, 86, maxColorValue = 255),rgb(254, 217, 106, maxColorValue = 255),rgb(254, 204, 92, maxColorValue = 255),rgb(254, 191, 78, maxColorValue = 255)),
                          stringsAsFactors = FALSE)
  }
  if(neg){
    periods$max_age <- -1 * (periods$max_age)
    periods$min_age <- -1 * (periods$min_age)
  }
  periods$mid_age <- (periods$max_age + periods$min_age)/2
  if(!is.null(fill)){
    periods$color <- rep(fill, length.out = nrow(periods))
  }else if(!("color" %in% colnames(periods))){
    periods$color <- rep(c("grey60","grey80"), length.out = nrow(periods))
  }
  lims <- ggplot_build(gg)$layout$panel_params[[1]]
  if(abbrv & "abbr" %in% colnames(periods)){
    periods$names <- periods$abbr
  }else{
    periods$names <- periods$period
  }
  if(!quat){
    periods$names[periods$abbr=="Q"] <- ""
  }
  if(pos %in% c("bottom", "top", "b", "t")){
    if(pos %in% c("top","t")){
      ymax <- max(lims$y.range)*1.1
      ymin <- max(lims$y.range) - height * (max(lims$y.range) - min(lims$y.range)) 
    }else{
      ymin <- min(lims$y.range)*1.1
      ymax <- min(lims$y.range) + height * (max(lims$y.range) - min(lims$y.range)) 
    }
    gg <- gg +
      annotate("rect", xmin = periods$min_age, xmax = periods$max_age, ymin = ymin, ymax = ymax,
               fill = periods$color, color = color, alpha = alpha) +
      annotate("text", x = periods$mid_age, label = periods$names, y = (ymin+ymax)/2,
               vjust = "middle", hjust = "middle", size = size)
    if (blank.gg) {
      gg <- gg + ylim(c(ymin, ymax))
    }
  }else if(pos %in% c("left", "right","l","r")){
    if(pos %in% c("right","r")){
      xmax <- max(lims$x.range)
      xmin <- max(lims$x.range) - height * (max(lims$x.range) - min(lims$x.range))
    }else{
      xmin <- min(lims$x.range)
      xmax <- min(lims$x.range) + height * (max(lims$x.range) - min(lims$x.range))
    }
    gg <- gg +
      annotate("rect", ymin = periods$min_age, ymax = periods$max_age, xmin = xmin, xmax = xmax,
               fill = periods$color, color = color, alpha = alpha) +
      annotate("text", y = periods$mid_age, label = periods$names, x = (xmin+xmax)/2,
               vjust = "middle", hjust = "middle", size = size, angle = 90)
    if (blank.gg) {
      gg <- gg + ylim(c(ymin, ymax))
    }
  }
  gg
}