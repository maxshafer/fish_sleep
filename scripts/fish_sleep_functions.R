


### Write a function to extract ancestral likelihoods from a model

returnAncestralStates <- function(phylo_model = model, phylo_tree = trpy_n) {
  
  # Create a data frame with trait values from reconstruction
  lik.anc <- as.data.frame(rbind(phylo_model$tip.states, phylo_model$states))
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
  
  colnames(lik.anc) <- rate_states
  
  ancestral_states <- list()
  ancestral_states$lik.anc <- lik.anc
  ancestral_states$node <- 1:(length(phylo_tree$tip.label) + phylo_tree$Nnode)
  ancestral_states$states <- states
  ancestral_states$rate_states <- rate_states
  return(ancestral_states)
}



# Function to determine the number of transitions between states

calculateStateTranstitions <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n) {
  
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
  ancestral_states$node.age <- node_heights[match(ancestral_states$node, phylo_tree$edge[,2]),1]
  ancestral_states$node.age[is.na(ancestral_states$node.age)] <- 0
  ancestral_states$node.age <- (ancestral_states$node.age - max(ancestral_states$node.age))*-1
  
  # ancestral_states <- ancestral_states[order(ancestral_states$Time),]
  
  # ID parental nodes and parental states
  ancestral_states$parental.node <- unlist(lapply(ancestral_states$node, function(x) Ancestors(phylo_tree, x, type = "parent")))
  ancestral_states$parent.diel <- unlist(lapply(ancestral_states$parental.node, function(x) ancestral_states$recon_states[match(x, ancestral_states$node)]))
  # parent of the root is NA
  ancestral_states$parent.diel[is.na(ancestral_states$parent.diel)] <- 0
  
  ancestral_states$transition <- ifelse(ifelse(ancestral_states$parent.diel > 0.5, 1, 0) != ifelse(ancestral_states$recon_states > 0.5, 1, 0), 1, 0)
  # ancestral_states$transition[is.na(ancestral_states$transition)] <- 0
  return(ancestral_states)
}



calculateLinTransHist <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n) {

  ancestors <- lapply(c(1:length(phylo_tree$tip.label)), function(x) Ancestors(phylo_tree, x, type = "all"))
  ancestors <- lapply(seq_along(ancestors), function(x) append(c(1:length(phylo_tree$tip.label))[[x]], ancestors[[x]]))
  recon_states <- ifelse(ancestral_states$recon_states > 0.5, 1, 0)
  ancestors.diel <- lapply(ancestors, function(x) lapply(x, function(y) recon_states[match(y, ancestral_states$node)]))
  
  # This works, can I simplify it so it works on a vector?
  switch.type <- unlist(lapply(ancestors.diel, function(x) {
    df <- data.frame(test = unlist(x))
    df <- df[with(df, c(test[-1]!= test[-nrow(df)], TRUE)),]
    return(paste(df, collapse = ""))
  }))
  
  ancestral_states$switch.type <- as.character(switch.type)
  
  return(ancestral_states)
}



returnCumSums <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n) {
  
  node_order <- order(ancestral_states$node.age, decreasing = T)
  
  ancestral_states$transition_cumsum <- cumsum(ancestral_states$transition[node_order])
  
  node_order_2 <- order(ancestral_states$node.age[ancestral_states$node %in% c(1:length(phylo_tree$tip.label))], decreasing = T)
    
  cumsums <- data.frame(row.names = phylo_tree$tip.label[node_order_2])
  cumsums$node.age <- ancestral_states$node.age[node_order_2]
  cumsums$Lineage_Cumsum <- cumsum(rep(1, length(phylo_tree$tip.label)))
  
  for (i in unique(ancestral_states$switch.type)) {
    cumsums[,i] <- cumsum(ifelse(ancestral_states$switch.type == i, 1, 0))
  }
  
}


node.data3 <- data.table::melt(cumsums, id.vars = c("node.age", "Lineage_Cumsum"), value.name = "Cummulative_ratio")

node.data3$CRLC <- node.data3$Cummulative_ratio
node.data4 <- node.data3 %>% group_by(node.age, variable) %>% summarise(n = sum(CRLC)) %>% mutate(percentage = n / sum(n))

switch_histo <- ggplot(node.data4, aes(x = node.age, y = percentage, fill = variable)) + theme_classic() + geom_area() + scale_x_reverse() + xlab("Millions of years ago") + scale_fill_manual(values = c("blue4", "red4", "blue3", "red3", "blue1", "red1"))

  
  
  
  











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

