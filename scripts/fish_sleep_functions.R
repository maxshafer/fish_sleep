### FUNCTIONS FOR LINEAGE THROUGH TIME PLOTS
### Associated with Fish_sleep project
### Copyright Maxwell Shafer 2022
### The following functions are associated with scripts for analysing diurnal and nocturnal states through time across a phylogeny
### and can be used with both data from fish, but also data from all vertebrates
### Updated 07.12.2022

### Need to make these work for if there is another trait involved (for example, marine/fresh). Would still be useful to plot the Di/Noc for these models

### Write a function to extract ancestral likelihoods from a model
returnAncestralStates <- function(phylo_model = model, phylo_tree = trpy_n) {
  
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
  colnames(transition) <- c(1:ncol(transition))

  print(paste("Identified an average of", mean(colSums(transition)), "transitions between", ancestral_states$states[1], "and", ancestral_states$states[2], "across", ncol(simulated_data), "simulations", sep = " "))
  
  return(transition)
}


### Function to determine the number of transitions between states
calculateStateTransitions <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n) {
  
  states <- ancestral_states$states
  rate_states <- ancestral_states$rate_states
  
  # Determine the number of transitions or diurnal/nocturnal taxa by age
  if (length(ancestral_states$states) == length(ancestral_states$rate_states)){
    ancestral_states$recon_states <- ancestral_states$lik.anc[,states[[1]]]
  } else {
    ancestral_states$recon_states <- as.numeric(rowSums(ancestral_states$lik.anc[,rate_states[grep(states[[1]], rate_states)]]))
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
  ancestral_states$parent.diel[is.na(ancestral_states$parent.diel)] <- 0
  
  ancestral_states$transition <- ifelse(ifelse(ancestral_states$parent.diel > 0.5, 1, 0) != ifelse(ancestral_states$recon_states > 0.5, 1, 0), 1, 0)
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
  
  noc_colours <- brewer.pal(noc_levels, "Blues") 
  di_colours <- brewer.pal(noc_levels, "Reds")
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
    node.data[nrow(node.data)+1,] <- list(max(ancestral_states$node.age), "Noc", 1, 1)
  }
  ## Should come up with function to make this? Yeah, and auto-find colors maybe from the red and blue scales
  plot <- ggplot(node.data, aes(x = node.age, y = percentage, fill = variable)) + theme_classic() + geom_area() + scale_x_reverse()
  plot <- plot + xlab("Millions of years ago") + ylab("Fraction of lineages") #+ scale_fill_manual(values = c("blue4", "red4", "blue3", "red3", "blue1", "red1"))
  

  
  plot <- plot + scale_fill_manual(values = colours)
  
  return(plot)
}


### Function to make the Switch ratio plot (line)

switchRatio <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, node.age.cutoff = 0.1) {
  require(ape)
  
  node_order <- order(ancestral_states$node.age, decreasing = T)
  df <- data.frame(node = ancestral_states$node[node_order], node.age = ancestral_states$node.age[node_order], transition_cumsum = ancestral_states$transition_cumsum)
  
  # ltt.data <- as.data.frame(ltt.plot.coords(phy = phylo_tree, backward = TRUE, tol = 0, type = "S"))
  # ltt.data$node.age <- (ltt.data$time - max(ltt.data$time))*-1
  
  ## calculate it the way the ape function does, but with my node.ages?
  ## This works close to the function, but neglects the node at the root
  df$ltt_source <- ifelse(df$node %in% c(1:Ntip(phylo_tree)), -1, 1)
  df$ltt_cumsum <- cumsum(df$ltt_source)

  ## OK this seems to be really close. I can calculate ltt's properly, but it also makes a weird thing were ltt decreases near present day
  df$ratio <- df$transition_cumsum/df$ltt_cumsum
  
  df$ratio[is.nan(df$ratio)] <- 0
  df$ratio[is.infinite(df$ratio)] <- 0
  
  plot <- ggplot(df[df$node.age > node.age.cutoff,], aes(x = node.age, y = ratio)) + geom_line() + theme_classic() + scale_x_reverse() #+ scale_x_reverse(limits = c(413,0.1)) + ylim(c(0,0.15)) # ltt would be the number of times it was possible to switch? Or something like that?
  plot <- plot + xlab("Millions of years ago") + ylab("Fraction of lineages transitioning")
  
  return(plot)
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