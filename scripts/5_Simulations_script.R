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

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/fish_sleep_functions.R")

## Load files

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
tr.calibrated <- readRDS("calibrated_phylo.rds")
trait.data <- readRDS("trait_data.rds")

name_variable <- "all"

# # Remove low quality data
# trait.data <- trait.data[trait.data$confidence > 1,]
# name_variable <- "only_highqual"

# # Remove cartilaginous fish (just actinopterygii)
# # Common ancestor of Lepisosteus_osseus (Order: Lepisosteiformes) and Lutjanus_fulvus (Order: Perciformes)
# node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Lepisosteus_osseus", "Lutjanus_fulvus"))
# tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
# trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
# name_variable <- "only_ingroup"

# # Keep only cartilaginous fish (no actinopterygii)
# node_of_interest <- getMRCA(tr.calibrated, tip = c("Rhizoprionodon_terraenovae", "Rhynchobatus_djiddensis"))
# tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
# trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
# name_variable <- "only_cartilaginous"


# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]


# Load the best model, which is the HMM 2 state 2 rate model
models <- readRDS(file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
model <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))



###############################################################################################################################################
### Simulate diel switchs based on derived model parameters ### 
###############################################################################################################################################

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
simulation_numb <- 100

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





## This works now

test <- calculateSimulatedTransitions(simulated_data = simulation, phylo_tree = trpy_n)

  
  
  
  

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