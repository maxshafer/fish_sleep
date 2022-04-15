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

