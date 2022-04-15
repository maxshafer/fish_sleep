library(phangorn)
library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(zoo)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/scripts/Trait_rate_matrix_figure_script.R")

## Load files

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
trpy <- readRDS("calibrated_phylo_tree_ancestral.rds")

# Load the best model, which is the HMM 2 state 2 rate model
model <- readRDS("best_fit_model.rds")

lik.anc <- as.data.frame(rbind(HMM_2state_2rate$tip.states, HMM_2state_2rate$states))
colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

node.data2 <- node.data <- lik.anc
### Plot the number of transitions or diurnal/nocturnal taxa by year
#node.data2 <- node.data[node.data$node %in% 1:length(trpy_n$tip.label),]
node.data2$diurnal <- node.data2$diurnal_R1 + node.data2$diurnal_R2
node.data2$Time <- branching.times(trpy_n)[node.data2$node]*-1
node.data2 <- node.data2[order(node.data2$Time),]
node.data2$Diurnal_Cumsum <- cumsum(ifelse(node.data2$diurnal > 0.5, 1, 0))
node.data2$Nocturnal_Cumsum <- cumsum(ifelse(node.data2$diurnal < 0.5, 1, 0))

node.data3 <- data.table::melt(node.data2[,c("Time", "Diurnal_Cumsum", "Nocturnal_Cumsum")], id.vars = "Time", value.name = "Cummulative_sum")


ltt.diel.plot <- ggplot(node.data3, aes(x = Time, y = Cummulative_sum, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_Cumsum" = "#d6604d", "Nocturnal_Cumsum" = "#4393c3"))
ltt.diel.plot <- ltt.diel.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10))

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ltt_", length(trpy$tip.label), "_species.pdf", sep = ""), width = 6, height = 4)
ltt.diel.plot
dev.off()




## 
# T
node_heights <- nodeHeights(trpy_n)

data.df <- node.data2
data.df$node.age <- node_heights[match(data.df$node, trpy_n$edge[,2]),1]
data.df$node.age[is.na(data.df$node.age)] <- 0
data.df$node.age <- (data.df$node.age - max(data.df$node.age))*-1

data.df$parental.node <- apply(data.df, 1, function(x) Ancestors(trpy_n, x[5], type = "parent"))
data.df$parent.diel <- apply(data.df, 1, function(x) ifelse(data.df$diurnal[x[5]] > 0.5, 1, 0))
data.df$diurnal <- ifelse(data.df$diurnal > 0.5, 1, 0)
data.df$switch <- ifelse(data.df$parent.diel != data.df$diurnal, "yes", "no")
data.df <- data.df[order(data.df$node.age, decreasing = T),]
data.df$switch2 <- cumsum(ifelse(data.df$parent.diel != data.df$diurnal, 1, 0))

histo <- ggplot(data.df, aes(x = node.age, fill = switch)) + geom_histogram(position = "fill", bins = 100)

histo.data <- ggplot_build(histo)$data[[1]][101:200,]

df.test <- data.frame(Time = histo.data$x, rollingmean = histo.data$y)
df.test$rollingmean[is.na(df.test$rollingmean)] <- 0
df.test <- df.test[order(df.test$Time),]
df.test$rollingmean2 <- rollmean(df.test$rollingmean, 2, align = "center", fill = 0)


histo / ggplot(df.test, aes(x = Time, y = rollingmean2)) + geom_point() + geom_path() + theme_classic()

ltt.switch.plot <- ggplot(data.df, aes(x = node.age, y = switch2)) + geom_line(size = 1.5) + theme_classic()

pdf(file = paste("outs/Figures/fish_phylogeny_switch_ltt_", length(trpy$tip.label), "_species.pdf", sep = ""), width = 4.25, height = 4)
ltt.switch.plot
dev.off()


## Plots the ratio between the cumulative switchs (diurnal to nocturnal and vice versa), and cumulative lineages

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/scripts/gggeo_scale.R")

data.df$switch3 <- data.df$switch2/(c(1:nrow(data.df)))

cum.ratio.plot <- ggplot(data.df, aes(x = node.age, y = (switch3))) + geom_line(size = 1.5) + theme_classic() #+ ylim(0, 0.6)
cum.ratio.plot <- gggeo_scale(cum.ratio.plot, pos = "top") + xlim(0,405.2853)



# Extract the data from the plot from the PNAS paper figur 1 (https://www.pnas.org/content/114/22/5653) using https://automeris.io/WebPlotDigitizer/
geo_frag <- read.csv("~/Downloads/Default Dataset.csv")
colnames(geo_frag) <- c("Time", "frag_index") # Time is millions of years ago
geo_frag$root <- (geo_frag$Time - max(geo_frag$Time))*-1
geo.frag.plot <- ggplot(geo_frag, aes(x = Time, y = frag_index)) + geom_line(size = 1.5) + theme_classic() + xlim(0,405.2853)

# Plot them together to see correlation

pdf(file = paste("outs/Figures/diel_switchs_cumsum_", length(trpy$tip.label), "_species.pdf", sep = ""), width = 10, height = 10)
cum.ratio.plot + plot_layout(nrow = 2) + geo.frag.plot
dev.off()


































## Now I need to add these to the phylo, then generate the lineages over time data and summarise that
## The below plots the ltt and the ratio of ltt and switches for a single dataset, I have 500 or more
## Write a function that does the below, then apply it over simulation

## OK trpy$edge is the same as nodeHeights(trpy)
## can use the indexing from trpy$edge to extract nodeHeights, which is the distance from the root
node_heights <- nodeHeights(trpy_n)

extractLTTData <- function(phylo = trpy, trait.data = simulation[,1], nodeHeights = node_heights) {
  nd <- data.frame(diurnal = ifelse(trait.data == "diurnal_R1" | trait.data == "diurnal_R2", 1, 0), node = 1:length(trait.data))
  # nd$node.age <- branching.times(phylo)[nd$node]*-1
  
  nd$node.age <- node_heights[match(nd$node, phylo$edge)]
  
  nd <- nd[order(nd$node.age),]
  #nd$Diurnal_Cumsum <- cumsum(ifelse(nd$diurnal > 0.5, 1, 0))
  #nd$Nocturnal_Cumsum <- cumsum(ifelse(nd$nocturnal > 0.5, 1, 0))
  
  #data.df <- nd
  #data.df$diurnal <- ifelse(data.df$diurnal >= 0.5, 1, 0)
  #data.df$node.age <- data.df$Time
  
  nd$parental.node <- apply(nd, 1, function(x) Ancestors(phylo, x[3], type = "parent"))
  nd$parent.diel[2:nrow(nd)] <- apply(nd[2:nrow(nd),], 1, function(x) nd$diurnal[x[[4]]])
  nd$parent.diel[1] <- nd$diurnal[1]
  nd$switch <- ifelse(nd$parent.diel != nd$diurnal, 1, 0)
  #nd <- nd[order(nd$node.age),]
  nd$switch2 <- cumsum(nd$switch)
  
  #ltt.switch.plot <- ggplot(data.df, aes(x = Time, y = log(switch2))) + geom_line(size = 1.5) + theme_classic()
  nd$switch3 <- nd$switch2/(c(1:nrow(nd)))
  #return(nd$switch3)
  return(nd$switch2)
}

simulation <- replicate(10000, rTraitDisc(trpy_n, model = "ARD", k = 2, rate = c(0.0074, 0.0059), states = c("nocturnal", "diurnal"), ancestor = T))

# Use the HMM model that fits best
simulation <- replicate(100, rTraitDisc(trpy_n, model = matrix(c(0,50.089200437,0.006767715,0, 99.907505815,0,0,0.006767715, 0.015389737,0,0,0.000000001, 0,0.01538974,0.00009360387,0), 4), states = c("diurnal_R1", "diurnal_R2", "nocturnal_R1", "nocturnal_R2"), ancestor = T, root.value = 4))


# The below is ordered by node.age! Not node number?
test <- apply(simulation, 2, function(x) extractLTTData(phylo = trpy_n, trait.data = x))
#test2 <- apply(simulation[,1:100], 2, function(x) extractLTTData(phylo = trpy, trait.data = x))

# Need to correct this

nodes <- c((length(trpy_n$tip.label)+1):(nrow(trpy_n$edge)+1), 1:length(trpy_n$tip.label))
node.age <- node_heights[match(nodes, trpy_n$edge)]

df <- data.frame(nodes = c((length(trpy_n$tip.label)+1):(nrow(trpy_n$edge)+1), 1:length(trpy_n$tip.label)), node.age = node_heights[match(nodes, trpy_n$edge)])
test.df <- df[order(df$node.age),]

test.df$means <- rowMeans(test)
stdev <- (unlist(apply(test, 1, function(x) sd(x))))
test.df$upper <- test.df$means + stdev
test.df$lower <- test.df$means - stdev
test.df$max <- apply(test, 1, function(x) max(x))
test.df$min <- apply(test, 1, function(x) min(x))
test.df$reconstruction <- data.df$switch3
#ci <- t(apply(test, 1, function(x) confidence_interval(vector = x, interval = 0.95)))
#test.df <- cbind(test.df, ci)

test.df$lineages <- cumsum(as.numeric(rownames(test.df)))


test.df.2 <- pivot_longer(test.df, cols = c("means", "upper", "lower", "reconstruction", "max", "min"))

ggplot(test.df.2, aes(x = as.numeric(node.age), y = value, group = name, colour = name)) + geom_line(size = 0.5) + theme_classic() #+ scale_colour_manual(values = c("lightblue", "black", "lightblue")) # + ylim(c(0,1))

df <- data.frame(diurnal = as.numeric(as.factor(simulation[,555])))

p1 <- ggtree(trpy_n, layout = "circular") %<+% data.frame(diurnal = simulation[,2]) + aes(color = simulation[,55]) + scale_color_manual(values = c("red", "blue")) #+ geom_tippoint(aes(color = simulation[,1]), shape = 16, size = 1.5)
p2 <- ggtree(trpy_n, layout = "circular") %<+% data.frame(diurnal = simulation[,2]) + aes(color = simulation[,99]) + scale_color_manual(values = c("red", "blue")) #+ geom_tippoint(aes(color = simulation[,1]), shape = 16, size = 1.5)

p1 + p2

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}


                      

plot <- ggplot(test.df, aes(x = node.age*-1, y = (means))) + geom_line(size = 1.5) + theme_classic() + ylim(c(0,1))

plot <- plot + ggplot(test.df, aes(x = node.age*-1, y = upper)) + geom_line(colour = "red", size = 1.5) + theme_classic() + ylim(c(0,1))

plot <- plot + ggplot(test.df, aes(x = node.age*-1, y = lower)) + geom_line(colour = "blue", size = 1.5) + theme_classic() + ylim(c(0,1))



ggplot(test.df, aes(x=node.age*-1, y = value)) + stat_summary(geom="ribbon", fun.ymin="min", fun.ymax="max", aes(fill=type), alpha=0.3) + theme_classic()
