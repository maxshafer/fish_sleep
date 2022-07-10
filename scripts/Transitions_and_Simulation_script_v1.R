library(phangorn)
library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(zoo)
library(ggtree)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/scripts/Trait_rate_matrix_figure_script.R")

## Load files

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
tr.calibrated <- readRDS("calibrated_phylo.rds")
trait.data <- readRDS("trait_data.rds")

name_variable <- "all"

# Remove low quality data
# trait.data <- trait.data[trait.data$confidence > 1,]

# Remove cartilaginous fish (just actinopterygii)
# Common ancestor of Lepisosteus_osseus (Order: Lepisosteiformes) and Lutjanus_fulvus (Order: Perciformes)
node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Lepisosteus_osseus", "Lutjanus_fulvus"))
tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
name_variable <- "only_ingroup"

# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]

# Make a tree for extracting node ages
tree <- ggtree(trpy_n)

# Load the best model, which is the HMM 2 state 2 rate model
# model <- readRDS(file = "best_fit_model_all_species.rds")
# model <- readRDS(file = paste("best_fit_model_", length(trpy_n$tip.label), "_species.rds", sep = ""))
model <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

###############################################################################################################################################
### Create the dataframe with lineages through time and diel switches data ### 
###############################################################################################################################################

# Create data frame with trait values
lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

# Determine the number of transitions or diurnal/nocturnal taxa by age
lik.anc$diurnal <- ifelse(lik.anc$diurnal_R1 + lik.anc$diurnal_R2 > 0.5, 1, 0)
# lik.anc$Time <- tree$data$x[match(lik.anc$node, tree$data$node)]
# lik.anc$Time2 <- (lik.anc$Time - max(lik.anc$Time))*-1
lik.anc <- lik.anc[order(lik.anc$Time),]

##
lik.anc$parental.node <- apply(lik.anc, 1, function(x) Ancestors(trpy_n, x["node"], type = "parent"))
lik.anc$parent.diel <- apply(lik.anc, 1, function(x) lik.anc$diurnal[match(x["parental.node"], lik.anc$node)]) # I should match, not index
lik.anc$parent.diel[is.na(lik.anc$parent.diel)] <- 0
lik.anc$switch <- ifelse(lik.anc$parent.diel != lik.anc$diurnal, 1, 0)
lik.anc$switch[is.na(lik.anc$switch)] <- 0

# Determine if switch is N->D (1) or D->N (0)
lik.anc$switch.ND <- ifelse(lik.anc$switch == 1 & lik.anc$diurnal == 1, 1, 0)
lik.anc$switch.DN <- ifelse(lik.anc$switch == 1 & lik.anc$diurnal == 0, 1, 0)

# Include or not
node_heights <- nodeHeights(trpy_n)
lik.anc$node.age <- node_heights[match(lik.anc$node, trpy_n$edge[,2]),1]
lik.anc$node.age[is.na(lik.anc$node.age)] <- 0
lik.anc$node.age <- (lik.anc$node.age - max(lik.anc$node.age))*-1

# Sort by node.age or Time2
lik.anc <- lik.anc[order(lik.anc$node.age, decreasing = T),]

# Determine cumulative sums
lik.anc$Diurnal_Cumsum <- cumsum(ifelse(lik.anc$diurnal > 0.5, 1, 0))
lik.anc$Nocturnal_Cumsum <- cumsum(ifelse(lik.anc$diurnal < 0.5, 1, 0))
lik.anc$Lineage_Cumsum <- cumsum(c(rep(1, nrow(lik.anc))))
lik.anc$Diurnal_ratio <- lik.anc$Diurnal_Cumsum/lik.anc$Lineage_Cumsum
lik.anc$Nocturnal_ratio <- lik.anc$Nocturnal_Cumsum/lik.anc$Lineage_Cumsum
lik.anc$switch2 <- cumsum(lik.anc$switch)
lik.anc$switch.NDCum <- cumsum(lik.anc$switch.ND)/lik.anc$Lineage_Cumsum
lik.anc$switch.DNCum <- cumsum(lik.anc$switch.DN)/lik.anc$Lineage_Cumsum

###############################################################################################################################################
### Make plots ### 
###############################################################################################################################################

# Make LTT plots
node.data3 <- data.table::melt(lik.anc[,c("node.age", "Diurnal_Cumsum", "Nocturnal_Cumsum", "Lineage_Cumsum")], id.vars = "node.age", value.name = "Cummulative_sum")
ltt.diel.plot <- ggplot(node.data3, aes(x = node.age, y = log(Cummulative_sum), colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_Cumsum" = "#d6604d", "Nocturnal_Cumsum" = "#4393c3"))
ltt.diel.plot <- ltt.diel.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

node.data3 <- data.table::melt(lik.anc[,c("node.age", "Diurnal_ratio", "Nocturnal_ratio")], id.vars = "node.age", value.name = "Cummulative_ratio")
ltt.diel.ratio.plot <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_ratio" = "#d6604d", "Nocturnal_ratio" = "#4393c3"))
ltt.diel.ratio.plot <- ltt.diel.ratio.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

node.data3 <- data.table::melt(lik.anc[,c("node.age", "switch.NDCum", "switch.DNCum")], id.vars = "node.age", value.name = "Cummulative_ratio")
ltt.diel.switch.dir.plot <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("switch.NDCum" = "#d6604d", "switch.DNCum" = "#4393c3"))
ltt.diel.switch.dir.plot <- ltt.diel.switch.dir.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

# Plots
switch_histo <- ggplot(lik.anc, aes(x = node.age, fill = as.character(switch))) + geom_histogram(position = "fill", bins = 100) + scale_x_reverse()
switch_ratio <- ggplot(lik.anc, aes(x = node.age, y = switch2/Lineage_Cumsum)) + geom_line() + theme_classic() + scale_x_reverse()
ltt.switch.plot <- ggplot(lik.anc, aes(x = node.age, y = log(switch2))) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + ylab("log(switches-through-time)") + xlab("Millions of years ago")

ggplot(lik.anc, aes(y = log(switch2), x = log(Lineage_Cumsum))) + geom_point() + theme_bw()

# Save plots

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ltt_switch_ratio", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 6, height = 4)
switch_histo / switch_ratio / ltt.diel.switch.dir.plot
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ltt", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 6, height = 4)
ltt.diel.plot + ltt.diel.ratio.plot + plot_layout(nrow = 2)
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_switch_stt", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 4.25, height = 4)
ltt.switch.plot
dev.off()


## Plots the ratio between the cumulative switchs (diurnal to nocturnal and vice versa), and cumulative lineages

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/scripts/gggeo_scale.R")


cum.ratio.plot <- ggplot(lik.anc, aes(x = node.age, y = switch2/Lineage_Cumsum)) + geom_line(size = 1.5) + theme_classic()
cum.ratio.plot <- gggeo_scale(cum.ratio.plot, pos = "top") + scale_x_reverse() + xlim(405.2853,0)


# Extract the data from the plot from the PNAS paper figur 1 (https://www.pnas.org/content/114/22/5653) using https://automeris.io/WebPlotDigitizer/
geo_frag <- read.csv("~/Downloads/Default Dataset.csv")
colnames(geo_frag) <- c("Time", "frag_index") # Time is millions of years ago
# geo_frag$root <- (geo_frag$Time - max(geo_frag$Time))*-1
geo.frag.plot <- ggplot(geo_frag, aes(x = Time, y = frag_index)) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + xlim(405.2853,0)
geo.frag.plot <- geo.frag.plot + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20))
# Plot them together to see correlation

pdf(file = paste("outs/Figures/fish_phylogeny_stt_geofrag_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 10, height = 10)
cum.ratio.plot + plot_layout(nrow = 2) + geo.frag.plot
dev.off()



pdf(file = paste("outs/Figures/fish_phylogeny_ltt_all_plots_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 15, height = 7)
cum.ratio.plot + ltt.diel.plot + ltt.diel.switch.dir.plot + ltt.diel.ratio.plot + plot_layout(nrow = 2)
dev.off()



















###############################################################################################################################################
### Simulate diel switchs based on derived model parameters ### 
###############################################################################################################################################

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

## Load files
resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
tr.calibrated <- readRDS("calibrated_phylo.rds")
trait.data <- readRDS("trait_data.rds")

name_variable <- "all"

# Remove low quality data
# trait.data <- trait.data[trait.data$confidence > 1,]
# name_variable <- "only_highqual"

# Remove cartilaginous fish (just actinopterygii)
# Common ancestor of Lepisosteus_osseus (Order: Lepisosteiformes) and Lutjanus_fulvus (Order: Perciformes)
node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Lepisosteus_osseus", "Lutjanus_fulvus"))
tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
name_variable <- "only_ingroup"

# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]

## Define function and the common variables

node_heights <- nodeHeights(trpy_n)
ancestors <- data.frame(node = c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label))), parental.node = Ancestors(trpy_n, c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label))), type = "parent"))

extractLTTData <- function(phylo = trpy_n, trait.data = simulation[,100], nodeHeights = node_heights, Ancestors = ancestors) {
  # Create data frame with trait values
  lik.anc.fnc <- as.data.frame(trait.data)
  lik.anc.fnc$diurnal <- ifelse(lik.anc.fnc$trait.data %in% c("diurnal_R1", "diurnal_R2", "diurnal"), 1, 0)
  lik.anc.fnc$nocturnal <- ifelse(lik.anc.fnc$trait.data %in% c("nocturnal_R1", "nocturnal_R2", "nocturnal"), 1, 0)
  
  lik.anc.fnc$node <- c(1:length(phylo$tip.label), (length(phylo$tip.label) + 1):(phylo$Nnode + length(phylo$tip.label)))
  
  lik.anc.fnc$parental.node <- Ancestors$parental.node
  lik.anc.fnc$parent.diel <- lik.anc.fnc$diurnal[match(lik.anc.fnc$parental.node, lik.anc.fnc$node)]
  #lik.anc.fnc$parent.diel <- apply(lik.anc.fnc, 1, function(x) lik.anc.fnc$diurnal[match(x["parental.node"], lik.anc.fnc$node)]) # I should match, not index
  lik.anc.fnc$parent.diel[is.na(lik.anc.fnc$parent.diel)] <- 0
  lik.anc.fnc$switch <- ifelse(lik.anc.fnc$parent.diel != lik.anc.fnc$diurnal, 1, 0)
  lik.anc.fnc$switch[is.na(lik.anc.fnc$switch)] <- 0
  
  # Determine if switch is N->D (1) or D->N (0)
  lik.anc.fnc$switch.ND <- ifelse(lik.anc.fnc$switch == 1 & lik.anc.fnc$diurnal == 1, 1, 0)
  lik.anc.fnc$switch.DN <- ifelse(lik.anc.fnc$switch == 1 & lik.anc.fnc$diurnal == 0, 1, 0)
  
  # Include or not
  
  lik.anc.fnc$node.age <- nodeHeights[match(lik.anc.fnc$node, phylo$edge[,2]),1]
  lik.anc.fnc$node.age[is.na(lik.anc.fnc$node.age)] <- 0
  lik.anc.fnc$node.age <- (lik.anc.fnc$node.age - max(lik.anc.fnc$node.age))*-1
  
  # Sort by node.age or Time2
  lik.anc.fnc <- lik.anc.fnc[order(lik.anc.fnc$node.age, decreasing = T),]
  
  # Determine cumulative sums
  lik.anc.fnc$Diurnal_Cumsum <- cumsum(ifelse(lik.anc.fnc$diurnal > 0.5, 1, 0))
  lik.anc.fnc$Nocturnal_Cumsum <- cumsum(ifelse(lik.anc.fnc$diurnal < 0.5, 1, 0))
  lik.anc.fnc$Lineage_Cumsum <- cumsum(c(rep(1, nrow(lik.anc.fnc))))
  lik.anc.fnc$Diurnal_ratio <- lik.anc.fnc$Diurnal_Cumsum/lik.anc.fnc$Lineage_Cumsum
  lik.anc.fnc$Nocturnal_ratio <- lik.anc.fnc$Nocturnal_Cumsum/lik.anc.fnc$Lineage_Cumsum
  lik.anc.fnc$switch2 <- cumsum(lik.anc.fnc$switch)
  lik.anc.fnc$switch.NDCum <- cumsum(lik.anc.fnc$switch.ND)/lik.anc.fnc$Lineage_Cumsum
  lik.anc.fnc$switch.DNCum <- cumsum(lik.anc.fnc$switch.DN)/lik.anc.fnc$Lineage_Cumsum
  lik.anc.fnc$return <- lik.anc.fnc$switch2/lik.anc.fnc$Lineage_Cumsum
  return(lik.anc.fnc)
}

## Choose which model to simulate:
model_type <- "ER"
simulation_numb <- 500

if (model_type == "ARD") {
  model <- ace(trait.vector_n, trpy_n, model = "ARD", type = "discrete")
  ## Simulate data based on the ARD model
  simulation <- replicate(simulation_numb, rTraitDisc(trpy_n, model = "ARD", k = 2, rate = c(0.0074, 0.0059), states = c("nocturnal", "diurnal"), ancestor = T))
}

# OK so there is a bug in the rTraitDisc code, it replaces

if (model_type == "ER") {
  model_d <- ace(trait.vector_n, trpy_n, model = "ER", type = "discrete")
  # Use the ER model (actually the ARD model but with symmetric rates, b/c I'm not sure how to code in ER model)
  simulation <- replicate(simulation_numb, rTraitDisc(phy = trpy_n, model = "SYM", k = 2, rate = model_d$rates, states = c("nocturnal", "diurnal"), ancestor = T))
}

if (model_type == "HR") {
  model_d <- readRDS(file = "best_fit_model_all_species.rds")
  # Use the best fit model
  solution <- model_d$solution
  solution[is.na(solution)] <- 0
  simulation <- replicate(simulation_numb, rTraitDisc(trpy_n, model = solution, states = c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2"), ancestor = T, root.value = 4))
}

saveRDS(simulation, file = paste("diel_switch_simulations", name_variable, model_type, ncol(simulation), "rTraitDisc.rds", sep = "_"))

# Extract the data from the simulations (ancestral reconstructions and switches)
# This is the full data.frame
extracted_data <- apply(simulation, 2, function(x) extractLTTData(phylo = trpy_n, trait.data = x))
extracted.df <- extracted_data[[1]][,c("node", "node.age")]

# extract a single column from simluated data (for example, ratio of switchs to lineages)
extracted_data_col <- Reduce(cbind, lapply(extracted_data, function(x) {
  return(x["return"])
}))

# Determine means, stdev, of cumulative switch numbers
extracted.df$means <- rowMeans(extracted_data_col)
stdev <- (unlist(apply(extracted_data_col, 1, function(x) sd(x))))
extracted.df$upper <- extracted.df$means + stdev
extracted.df$lower <- extracted.df$means - stdev
extracted.df$max <- apply(extracted_data_col, 1, function(x) max(x))
extracted.df$min <- apply(extracted_data_col, 1, function(x) min(x))
# This comes from the actual reconstruction model, and is made above
extracted.df$reconstruction <- lik.anc$switch2/lik.anc$Lineage_Cumsum
extracted.df$lineages <- cumsum(rep(1, nrow(extracted.df)))

# Convert for plotting
extracted.df <- pivot_longer(extracted.df, cols = c("means", "upper", "lower", "reconstruction"))
extracted.df$name <- factor(extracted.df$name, levels = c("upper", "means", "lower", "reconstruction"))

recon_plot <- ggplot(extracted.df, aes(x = as.numeric(node.age), y = value, group = name, colour = name)) + geom_line(size = 1) + theme_classic() + scale_x_reverse() + scale_color_manual(values = c("lightblue", "blue", "lightblue", "red")) + ggtitle(paste("Simulations vs Reconstruction - ", model_type, " model - ", ncol(simulation), " simulations", sep = "")) #+ ylim(c(0,0.385))


pdf(file = paste("outs/Figures/fish_phylogeny_stt", name_variable, model_type, ncol(simulation), length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 7.5, height = 5)
recon_plot
dev.off()
















ggplot(test.df, aes(x = as.numeric(node.age), y = log(reconstruction) - log(means))) + geom_line() + theme_bw()


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
