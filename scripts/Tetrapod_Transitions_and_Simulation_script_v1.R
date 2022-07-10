## load packages, trees, and data. For this package only I parsed down the data file to only two columns (species and activity state) and created separate files for diurnal and nocturnal. This made referencing the data simpler in the corHMM function. 
library(ape) 
library(corHMM)
library(ggplot2)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

ericson <- read.nexus("tetrapod_data/Ericson1914.tree")
hackett <- read.nexus("tetrapod_data/Hackett1914.tree") 
# Full data with Fa
dataTet <- read.csv("tetrapod_data/suppfile10appendix1914.csv")
dataD <- read.csv("tetrapod_data/suppfile11diu1914.csv")
dataN <- read.csv("tetrapod_data/suppfile12noc1914.csv")

# Remove mammals
dataN <- dataN[dataN$Species %in% dataTet$Species[dataTet$Class != "Mammalia"],]
hackett <- keep.tip(hackett, tip = dataN$Species)

# ## run corHMM hidden rates model to determine fit, using default root state. Compare different numbers of hidden rates by running the model with all possible rate categories, from 1 to 5. The example listed below is for a 2 rate model. To run all additional corHMM models, simply change the rate.cat function. More rate catergories take more time to run. 
# ## ericson nocturnal, 2 rates 
# EN2.hmm<-corHMM(ericson, dataN[,c(1,2)], rate.cat = 2, node.states = c("marginal"))
# ## ericson diurnal, 2 rates
# ED2.hmm<-corHMM(ericson, dataD[,c(1,2)], rate.cat = 2, node.states = c("marginal")) 
# ## hackett nocturnal, 2 rates 
# HN2.hmm<-corHMM(hackett, dataN[,c(1,2)], rate.cat = 2, node.states = c("marginal"))
# ## hackett diurnal, 2 rates
# HD2.hmm<-corHMM(hackett, dataD[,c(1,2)], rate.cat = 2, node.states = c("marginal"))
# 
# ## compare corHMM models with different rates using the AICc and loglik values within the same tree and coding scheme only
# EN2.hmm$loglik
# EN2.hmm$AICc
# ED2.hmm$loglik 
# ED2.hmm$AICc
# 
# HN2.hmm$loglik
# HN2.hmm$AICc
# HD2.hmm$loglik 
# HD2.hmm$AICc

## test the affect of different root treatments on the fit of the model. The default in corHMM is root.p=NULL, assuming an equal weighting of all states at the root. The methods "yang" and "madfitz" can also be applied to weight the likelihoods of each state according to transition rates on the tree. This was only done on models using the best fitting number of hidden rate classes.
## testing "yang" root treatment
## ericson nocturnal
# EN2.yang <-corHMM(ericson, dataN[,c(1,2)], rate.cat = 2, node.states = c("marginal"), root.p="yang")
## ericson diurnal
# ED2.yang <-corHMM(ericson, dataD[,c(1,2)], rate.cat = 2, node.states = c("marginal"), root.p="yang") 
## hackett nocturnal
# HN2.yang <-corHMM(hackett, dataN[,c(1,2)], rate.cat = 2, node.states = c("marginal"), root.p="yang")
## hackett diurnal
# HD3.yang <-corHMM(hackett, dataD[,c(1,2)], rate.cat = 3, node.states = c("marginal"), root.p="yang")
## testing "madfitz" root treatment
## ericson nocturnal
# EN2.madfitz <-corHMM(ericson, dataN[,c(1,2)], rate.cat = 2, node.states = c("marginal"), root.p="madfitz")
## ericson diurnal
# ED2.madfitz <-corHMM(ericson, dataD[,c(1,2)], rate.cat = 2, node.states = c("marginal"), root.p="madfitz") 
## hackett nocturnal 
HN2.madfitz <-corHMM(hackett, dataN[,c(1,2)], rate.cat = 2, node.states = c("marginal"), root.p="madfitz")
## hackett diurnal
# HD3.madfitz <-corHMM(hackett, dataD[,c(1,2)], rate.cat = 3, node.states = c("marginal"), root.p="madfitz")

## for the best fitting number of rate classes and the best fitting root treatment, plot ancestral state reconstructions as pie charts at the nodes of the phylogeny use plotRECON
## colors are ordered (0_R1, 1_R1, 0_R2, 1_R2) with R1 (rate 1) being slower than R2
## ericson nocturnal 
plotRECON(ericson, EN2.madfitz$states, piecolors = c("yellow", "cornflowerblue", "orange", "dark blue"), cex = 0.1, pie.cex = 0.15, file = "EricsonNoc2.corHMM.pdf", height = 50, width = 12, show.tip.label = TRUE, label.offset = 0.4, title = "Ericson Nocturnal Two Rates")
## ericson diurnal 
plotRECON(ericson, ED2.madfitz$states, piecolors = c("yellow", "cornflowerblue", "orange", "dark blue"), cex = 0.1, pie.cex = 0.15, file = "EricsonDiu2.corHMM.pdf", height = 50, width = 12, show.tip.label = TRUE, label.offset = 0.4, title = "Ericson Diurnal Two Rates")
## hackett nocturnal 
plotRECON(hackett, HN2.madfitz$states, piecolors = c("yellow", "cornflowerblue", "orange", "dark blue"), cex = 0.1, pie.cex = 0.15, file = "HackettNoc2.corHMM.pdf", height = 50, width = 12, show.tip.label = TRUE, label.offset = 0.4, title = "Hackett Nocturnal Two Rates")
## hackett diurnal 
plotRECON(hackett, HD3.madfitz$states, piecolors = c("yellow", "cornflowerblue", "goldenrod1", "blue2", "orange", "dark blue"), cex = 0.1, pie.cex = 0.15, file = "HackettDiu3.corHMM.pdf", height = 50, width = 12, show.tip.label = TRUE, label.offset = 0.4, title = "Hackett Diurnal Three Rates")

## end




# 0 is Diurnal, 1 is Nocturnal for states
model <- HN2.madfitz
trpy_n <- keep.tip(hackett, tip = dataN$Species)

###############################################################################################################################################
### Plot the Hidden rate ancestral reconstructions (of rate and state) ### 
###############################################################################################################################################

lik.anc.tet <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc.tet) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
lik.anc.tet$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

lik.anc.tet$noc.sum <- rowSums(lik.anc.tet[,c(2,4)])
lik.anc.tet$di.sum <- rowSums(lik.anc.tet[,c(1,3)])

lik.anc.tet$R1.sum <- rowSums(lik.anc.tet[,c(1,2)])
lik.anc.tet$R2.sum <- rowSums(lik.anc.tet[,c(3,4)])

# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
diel_2state_2rate <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
diel_2state_2rate_rect <- ggtree(trpy_n) %<+% lik.anc.tet + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")

rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")


dir1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = diurnal_R1) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1)
dir2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = diurnal_R2) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
nocr1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = nocturnal_R1) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)
nocr2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = nocturnal_R2) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)
dir3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = diurnal_R3) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "YlOrBr", direction = 1)
nocr3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.tet + aes(color = nocturnal_R3) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "BuPu", direction = 1)

two_rate_plots <- dir1 + dir2 + nocr1 + nocr2 + plot_layout(nrow = 2)

# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_diel_rates_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width =20, height = 20)
# two_rate_plots
# dev.off()
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_rates_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 20, height = 20)
# rate_plot
# dev.off()




#### Plot switches over time ####

###############################################################################################################################################
### Create the dataframe with lineages through time and diel switches data ### 
###############################################################################################################################################

# Create data frame with trait values
lik.anc.tet <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc.tet) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
lik.anc.tet$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

# Determine the number of transitions or diurnal/nocturnal taxa by age
lik.anc.tet$diurnal <- ifelse(lik.anc.tet$diurnal_R1 + lik.anc.tet$diurnal_R2 > 0.5, 1, 0)
# lik.anc.tet$Time <- tree$data$x[match(lik.anc.tet$node, tree$data$node)]
# lik.anc.tet$Time2 <- (lik.anc.tet$Time - max(lik.anc.tet$Time))*-1
lik.anc.tet <- lik.anc.tet[order(lik.anc.tet$Time),]

##
lik.anc.tet$parental.node <- apply(lik.anc.tet, 1, function(x) Ancestors(trpy_n, x["node"], type = "parent"))
lik.anc.tet$parent.diel <- apply(lik.anc.tet, 1, function(x) lik.anc.tet$diurnal[match(x["parental.node"], lik.anc.tet$node)]) # I should match, not index
lik.anc.tet$parent.diel[is.na(lik.anc.tet$parent.diel)] <- 0
lik.anc.tet$switch <- ifelse(lik.anc.tet$parent.diel != lik.anc.tet$diurnal, 1, 0)
lik.anc.tet$switch[is.na(lik.anc.tet$switch)] <- 0

# Determine if switch is N->D (1) or D->N (0)
lik.anc.tet$switch.ND <- ifelse(lik.anc.tet$switch == 1 & lik.anc.tet$diurnal == 1, 1, 0)
lik.anc.tet$switch.DN <- ifelse(lik.anc.tet$switch == 1 & lik.anc.tet$diurnal == 0, 1, 0)

# Include or not
node_heights <- nodeHeights(trpy_n)
lik.anc.tet$node.age <- node_heights[match(lik.anc.tet$node, trpy_n$edge[,2]),1]
lik.anc.tet$node.age[is.na(lik.anc.tet$node.age)] <- 0
lik.anc.tet$node.age <- (lik.anc.tet$node.age - max(lik.anc.tet$node.age))*-1

# Sort by node.age or Time2
lik.anc.tet <- lik.anc.tet[order(lik.anc.tet$node.age, decreasing = T),]

# Determine cumulative sums
lik.anc.tet$Diurnal_Cumsum <- cumsum(ifelse(lik.anc.tet$diurnal > 0.5, 1, 0))
lik.anc.tet$Nocturnal_Cumsum <- cumsum(ifelse(lik.anc.tet$diurnal < 0.5, 1, 0))
lik.anc.tet$Lineage_Cumsum <- cumsum(c(rep(1, nrow(lik.anc.tet))))
lik.anc.tet$Diurnal_ratio <- lik.anc.tet$Diurnal_Cumsum/lik.anc.tet$Lineage_Cumsum
lik.anc.tet$Nocturnal_ratio <- lik.anc.tet$Nocturnal_Cumsum/lik.anc.tet$Lineage_Cumsum
lik.anc.tet$switch2 <- cumsum(lik.anc.tet$switch)
lik.anc.tet$switch.NDCum <- cumsum(lik.anc.tet$switch.ND)/lik.anc.tet$Lineage_Cumsum
lik.anc.tet$switch.DNCum <- cumsum(lik.anc.tet$switch.DN)/lik.anc.tet$Lineage_Cumsum


###############################################################################################################################################
### Make plots ### 
###############################################################################################################################################

# Make LTT plots
node.data3 <- data.table::melt(lik.anc.tet[,c("node.age", "Diurnal_Cumsum", "Nocturnal_Cumsum", "Lineage_Cumsum")], id.vars = "node.age", value.name = "Cummulative_sum")
ltt.diel.plot.tet <- ggplot(node.data3, aes(x = node.age, y = log(Cummulative_sum), colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_Cumsum" = "#d6604d", "Nocturnal_Cumsum" = "#4393c3"))
ltt.diel.plot.tet <- ltt.diel.plot.tet + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

node.data3 <- data.table::melt(lik.anc.tet[,c("node.age", "Diurnal_ratio", "Nocturnal_ratio")], id.vars = "node.age", value.name = "Cummulative_ratio")
ltt.diel.ratio.plot.tet <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_ratio" = "#d6604d", "Nocturnal_ratio" = "#4393c3"))
ltt.diel.ratio.plot.tet <- ltt.diel.ratio.plot.tet + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

node.data3 <- data.table::melt(lik.anc.tet[,c("node.age", "switch.NDCum", "switch.DNCum")], id.vars = "node.age", value.name = "Cummulative_ratio")
ltt.diel.switch.dir.plot.tet <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("switch.NDCum" = "#d6604d", "switch.DNCum" = "#4393c3"))
ltt.diel.switch.dir.plot.tet <- ltt.diel.switch.dir.plot.tet + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")

# Plots
switch_histo.tet <- ggplot(lik.anc.tet, aes(x = node.age, fill = as.character(switch))) + geom_histogram(position = "fill", bins = 100) + scale_x_reverse()
switch_ratio.tet <- ggplot(lik.anc.tet, aes(x = node.age, y = switch2/Lineage_Cumsum)) + geom_line() + theme_classic() + scale_x_reverse()
ltt.switch.plot.tet <- ggplot(lik.anc.tet, aes(x = node.age, y = log(switch2))) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + ylab("log(switches-through-time)") + xlab("Millions of years ago")



## Put the two together (get lik.anc from other script)

lik.anc$group <- "fish"
lik.anc.tet$group <- "tetrapods"
lik.anc.mam$group <- "mammals"

test <- rbind(lik.anc, lik.anc.tet, lik.anc.mam)
node.data3 <- data.table::melt(test[,c("node.age", "group", "switch.NDCum", "switch.DNCum")], id.vars = c("node.age","group"), value.name = "Cummulative_ratio")
node.data3$plot.group <- paste(node.data3$group, node.data3$variable, sep = "_")

switch_type_plot <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = plot.group)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("fish_switch.NDCum" = "red", "fish_switch.DNCum" = "blue","tetrapods_switch.NDCum" = "indianred1", "tetrapods_switch.DNCum" = "dodgerblue1","mammals_switch.NDCum" = "orange1", "mammals_switch.DNCum" = "slateblue1")) + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")
switch_type_plot <- switch_type_plot + scale_x_reverse(breaks = scales::pretty_breaks(n = 20)) + theme_bw()
switch_type_plot

switch_plot <- ggplot(test, aes(x = node.age, y = switch2/Lineage_Cumsum, colour = group)) + geom_line(size = 1.5) + scale_x_reverse(breaks = scales::pretty_breaks(n = 40)) + theme_bw()
switch_plot
ggplot(test, aes(x = node.age, y = log(switch2), colour = group)) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + ylab("log(switches-through-time)") + xlab("Millions of years ago")


pdf(file = paste("outs/Figures/fish_phylogeny_stt_mammals_tetrapods", length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 25, height = 10)
switch_plot + switch_type_plot
dev.off()




gggeo_scale(cum.ratio.plot, pos = "top") + scale_x_reverse() + xlim(405.2853,0)

co2 <- read.csv("CO2_levels_paleo.csv", header = F)
colnames(co2) <- c("Time", "CO2_levels")

co2_plot <- ggplot(co2, aes(x = Time, y = CO2_levels)) + geom_point() + geom_line() + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20)) + scale_y_reverse()

temp <- read.csv("Temp_levels_paleo.csv", header = F)
colnames(temp) <- c("Time", "temp")

temp_plot <- ggplot(temp, aes(x = Time, y = temp)) + geom_point() + geom_line() + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20))

speciation <- read.csv("Speciation_paleo.csv", header = F)
colnames(speciation) <- c("Time", "speciation")

speciation_plot <- ggplot(speciation, aes(x = Time, y = speciation)) + geom_point() + geom_line() + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20))


occurances_conodonta <- read.csv("paleobiodb/pbdb_data_conodonta_summary.csv")
occurances_conodonta$group <- "Conodonta"
occurances_actino <- read.csv("paleobiodb/pbdb_data_actinopterygii_summary.csv")
occurances_actino$group <- "Actinopterygii"
occurances_tetrapoda <- read.csv("paleobiodb/pbdb_data_tetrapoda_summary.csv")
occurances_tetrapoda$group <- "Tetrapoda"
occurances_ammo <- read.csv("paleobiodb/pbdb_data_ammonoidea_summary.csv")
occurances_ammo$group <- "Ammonoidea"
occurances_ichthyosaurs <- read.csv("paleobiodb/pbdb_data_ichthyosaurs_summary.csv")
occurances_ichthyosaurs$group <- "Ichthyosaurs"
occurances <- rbind(occurances_conodonta, occurances_actino, occurances_tetrapoda, occurances_ammo, occurances_ichthyosaurs)

occurances_plot <- ggplot(occurances[occurances$group %in% c("Conodonta", "Actinopterygii", "Ammonoidea"),], aes(x = min_ma, y = sampled_in_bin, colour = group)) + geom_line(size = 1) + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20))

switch_plot / occurances_plot

pdf(file = paste("outs/Figures/fish_phylogeny_stt_mammals_tetrapods_PaleoData", length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 15, height = 25)
co2_plot / temp_plot / switch_plot / speciation_plot / geo.frag.plot
dev.off()

switch_type_plot / temp_plot
co2_plot / plot / temp_plot





