library(ape) 
library(corHMM)
library(xlsx)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")


# load data for mammals from Cox et al
# First is the tree
mam.table <- read.xlsx("Cox_mammal_data/Supplementary Data 2.xlsx", 1)

mammal_trees <- read.nexus("Cox_mammal_phylo/Complete_phylogeny.nex")

## Determine the consensus tree
mam.tree <- maxCladeCred(mammal_trees, tree = TRUE)


## There are two sources for diel activity, not sure which is which? Can I combine them?
mam.table$diel.comb <- paste(mam.table$Activity_DD, mam.table$Activity_IM, sep = "_")


mam.data <- mam.table[,c(1,7)]
mam.data$Binomial_iucn <- str_replace(mam.data$Binomial_iucn, " ", "_")
mam.data <- mam.data[!(is.na(mam.data$Activity_DD)),]
mam.data <- mam.data[mam.data$Activity_DD %in% c("Diurnal", "Nocturnal"),]

# Reciprically determine matching names
mam.data <- mam.data[mam.data$Binomial_iucn %in% mam.tree$tip.label,]
mam.data <- mam.data[!(duplicated(mam.data)),]
row.names(mam.data) <- mam.data$Binomial_iucn

trpy_n_mam <- keep.tip(mam.tree, tip = mam.data$Binomial_iucn)






###############################################################################################################################################
### Run Hidden Rates Models ### 
###############################################################################################################################################

# Run corHMM for just diel with 1 rat category (Markov model)
MK_2state <- corHMM(phy = trpy_n_mam, data = mam.data[trpy_n_mam$tip.label, c("Binomial_iucn", "Activity_DD")], rate.cat = 1, model = "ARD")

# Run corHMM for just diel, with 2 or 3 rate categories (large improvement for 2, diminishing returns for 3)
HMM_2state_2rate <- corHMM(phy = trpy_n_mam, data = mam.data[trpy_n_mam$tip.label, c("Binomial_iucn", "Activity_DD")], rate.cat = 2, model = "ARD", get.tip.states = TRUE)

# Save out 2 state and HMM with 2 rates
standard_tests <- list(c(MK_2state, HMM_2state_2rate))
saveRDS(standard_tests, file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

HMM_2state_3rate <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 3, model = "ARD", get.tip.states = TRUE)
standard_tests[[3]] <- HMM_2state_3rate
saveRDS(standard_tests, file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

















#### Plot switches over time ####

###############################################################################################################################################
### Create the dataframe with lineages through time and diel switches data ### 
###############################################################################################################################################

# Create data frame with trait values
lik.anc.mam <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc.mam) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
lik.anc.mam$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

# Determine the number of transitions or diurnal/nocturnal taxa by age
lik.anc.mam$diurnal <- ifelse(lik.anc.mam$diurnal_R1 + lik.anc.mam$diurnal_R2 > 0.5, 1, 0)
# lik.anc.mam$Time <- tree$data$x[match(lik.anc.mam$node, tree$data$node)]
# lik.anc.mam$Time2 <- (lik.anc.mam$Time - max(lik.anc.mam$Time))*-1
# lik.anc.mam <- lik.anc.mam[order(lik.anc.mam$Time),]

##
lik.anc.mam$parental.node <- apply(lik.anc.mam, 1, function(x) Ancestors(trpy_n, x["node"], type = "parent"))
lik.anc.mam$parent.diel <- apply(lik.anc.mam, 1, function(x) lik.anc.mam$diurnal[match(x["parental.node"], lik.anc.mam$node)]) # I should match, not index
lik.anc.mam$parent.diel[is.na(lik.anc.mam$parent.diel)] <- 0
lik.anc.mam$switch <- ifelse(lik.anc.mam$parent.diel != lik.anc.mam$diurnal, 1, 0)
lik.anc.mam$switch[is.na(lik.anc.mam$switch)] <- 0

# Determine if switch is N->D (1) or D->N (0)
lik.anc.mam$switch.ND <- ifelse(lik.anc.mam$switch == 1 & lik.anc.mam$diurnal == 1, 1, 0)
lik.anc.mam$switch.DN <- ifelse(lik.anc.mam$switch == 1 & lik.anc.mam$diurnal == 0, 1, 0)

# Include or not
node_heights <- nodeHeights(trpy_n)
lik.anc.mam$node.age <- node_heights[match(lik.anc.mam$node, trpy_n$edge[,2]),1]
lik.anc.mam$node.age[is.na(lik.anc.mam$node.age)] <- 0
lik.anc.mam$node.age <- (lik.anc.mam$node.age - max(lik.anc.mam$node.age))*-1

# Sort by node.age or Time2
lik.anc.mam <- lik.anc.mam[order(lik.anc.mam$node.age, decreasing = T),]

# Determine cumulative sums
lik.anc.mam$Diurnal_Cumsum <- cumsum(ifelse(lik.anc.mam$diurnal > 0.5, 1, 0))
lik.anc.mam$Nocturnal_Cumsum <- cumsum(ifelse(lik.anc.mam$diurnal < 0.5, 1, 0))
lik.anc.mam$Lineage_Cumsum <- cumsum(c(rep(1, nrow(lik.anc.mam))))
lik.anc.mam$Diurnal_ratio <- lik.anc.mam$Diurnal_Cumsum/lik.anc.mam$Lineage_Cumsum
lik.anc.mam$Nocturnal_ratio <- lik.anc.mam$Nocturnal_Cumsum/lik.anc.mam$Lineage_Cumsum
lik.anc.mam$switch2 <- cumsum(lik.anc.mam$switch)
lik.anc.mam$switch.NDCum <- cumsum(lik.anc.mam$switch.ND)/lik.anc.mam$Lineage_Cumsum
lik.anc.mam$switch.DNCum <- cumsum(lik.anc.mam$switch.DN)/lik.anc.mam$Lineage_Cumsum


###############################################################################################################################################
### Make plots ### 
###############################################################################################################################################

# Make LTT plots
node.data3 <- data.table::melt(lik.anc.mam[,c("node.age", "Diurnal_Cumsum", "Nocturnal_Cumsum", "Lineage_Cumsum")], id.vars = "node.age", value.name = "Cummulative_sum")
ltt.diel.plot.mam <- ggplot(node.data3, aes(x = node.age, y = log(Cummulative_sum), colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_Cumsum" = "#d6604d", "Nocturnal_Cumsum" = "#4393c3"))
ltt.diel.plot.mam <- ltt.diel.plot.mam + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

node.data3 <- data.table::melt(lik.anc.mam[,c("node.age", "Diurnal_ratio", "Nocturnal_ratio")], id.vars = "node.age", value.name = "Cummulative_ratio")
ltt.diel.ratio.plot.mam <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_ratio" = "#d6604d", "Nocturnal_ratio" = "#4393c3"))
ltt.diel.ratio.plot.mam <- ltt.diel.ratio.plot.mam + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

node.data3 <- data.table::melt(lik.anc.mam[,c("node.age", "switch.NDCum", "switch.DNCum")], id.vars = "node.age", value.name = "Cummulative_ratio")
ltt.diel.switch.dir.plot.mam <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("switch.NDCum" = "#d6604d", "switch.DNCum" = "#4393c3"))
ltt.diel.switch.dir.plot.mam <- ltt.diel.switch.dir.plot.mam + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

# Plots
switch_histo.mam <- ggplot(lik.anc.mam, aes(x = node.age, fill = as.character(switch))) + geom_histogram(position = "fill", bins = 100) + scale_x_reverse()
switch_ratio.mam <- ggplot(lik.anc.mam, aes(x = node.age, y = switch2/Lineage_Cumsum)) + geom_line() + theme_classic() + scale_x_reverse()
ltt.switch.plot.mam <- ggplot(lik.anc.mam, aes(x = node.age, y = log(switch2))) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + ylab("log(switches-through-time)") + xlab("Millions of years ago")



###############################################################################################################################################
### Plot the Hidden rate ancestral reconstructions (of rate and state) ### 
###############################################################################################################################################

lik.anc.mam <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc.mam) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
lik.anc.mam$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

lik.anc.mam$noc.sum <- rowSums(lik.anc.mam[,c(2,4)])
lik.anc.mam$di.sum <- rowSums(lik.anc.mam[,c(1,3)])

lik.anc.mam$R1.sum <- rowSums(lik.anc.mam[,c(1,2)])
lik.anc.mam$R2.sum <- rowSums(lik.anc.mam[,c(3,4)])

# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
diel_2state_2rate <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
diel_2state_2rate_rect <- ggtree(trpy_n) %<+% lik.anc.mam + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")

rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")


dir1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = diurnal_R1) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1)
dir2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = diurnal_R2) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
nocr1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = nocturnal_R1) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)
nocr2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = nocturnal_R2) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)
dir3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = diurnal_R3) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "YlOrBr", direction = 1)
nocr3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.mam + aes(color = nocturnal_R3) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "BuPu", direction = 1)

two_rate_plots <- dir1 + dir2 + nocr1 + nocr2 + plot_layout(nrow = 2)

# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_diel_rates_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width =20, height = 20)
# two_rate_plots
# dev.off()
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_rates_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 20, height = 20)
# rate_plot
# dev.off()




