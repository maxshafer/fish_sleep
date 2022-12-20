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

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/Fish_sleep_functions.R")

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

# View(unlist(lapply(models, function(x) x[names(x[grep("loglik", names(x))])])))

###############################################################################################################################################
### Create the object with lineages through time and diel switches data ### 
###############################################################################################################################################

# First, extract the ancestral states from the best fit model
anc_states <- returnAncestralStates(phylo_model = model, phylo_tree = trpy_n)

# Then, calculate transitions between states
anc_states <- calculateStateTransitions(ancestral_states = anc_states, phylo_tree = trpy_n)

# Determine transition histories (types of lineages)
anc_states <- calculateLinTransHist(ancestral_states = anc_states, phylo_tree = trpy_n)

# Calculate cumsums through time (for ltt plots)
anc_states <- returnCumSums(ancestral_states = anc_states, phylo_tree = trpy_n)


###############################################################################################################################################
### Make Plots! ### 
###############################################################################################################################################

switch.histo <- switchHisto(ancestral_states = anc_states, replace_variable_names = T, backfill = T)

# So, this works now, and I think correctly
# One question is though, because all the tips are clustered at the same age, there's a lot of tips and transitions right near the end of the tree
# Transistions are associated with the tip or node that is reconstructed, which is why this appears like this
# I could associate them with the parental node, but I don't have evidence that that is when it was reconstructed (which would fix the artifact at the end)

switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02)

# Add geological timescale
switch.ratio <- gggeo_scale(switch.ratio, pos = "bottom") + scale_x_reverse() + xlim(430,0)

# This highlights which nodes have undergone the most transitions
numb_switch_tree <- switchTree(ancestral_states = anc_states, phylo_tree = trpy_n, layout = "circular", replace_variable_names = TRUE)

# Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void()



###############################################################################################################################################
### Save Plots! ### 
###############################################################################################################################################


# This really shows the difference between methods, and where to assign the switch
# Also super affected by the time calibration

xlims <- c(max(anc_states$node.age), min(anc_states$node.age))

pdf(file = paste("outs/Figures/fish_phylogeny_diel_plot_transitions_histo", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
((switch.histo / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5))
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_plot_transitions_swith", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 10, height = 5)
((switch.ratio / geo_scale) & xlim(xlims)) + plot_layout(nrow = 2, heights = c(5,0.5))
dev.off()


pdf(file = paste("outs/Figures/fish_phylogeny_diel_highswitchlineages", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), height = 60, width = 60)
numb_switch_tree
dev.off()









# ## Maybe put this in a different script
# 
# # Extract the data from the plot from the PNAS paper figur 1 (https://www.pnas.org/content/114/22/5653) using https://automeris.io/WebPlotDigitizer/
# geo_frag <- read.csv("~/Downloads/Default Dataset.csv")
# colnames(geo_frag) <- c("Time", "frag_index") # Time is millions of years ago
# # geo_frag$root <- (geo_frag$Time - max(geo_frag$Time))*-1
# geo.frag.plot <- ggplot(geo_frag, aes(x = Time, y = frag_index)) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + xlim(405.2853,0)
# geo.frag.plot <- geo.frag.plot + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20))
# # Plot them together to see correlation
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_stt_geofrag_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 10, height = 10)
# cum.ratio.plot + plot_layout(nrow = 2) + geo.frag.plot
# dev.off()











# # Create data frame with trait values
# lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
# colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
# lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
# 
# # Determine the number of transitions or diurnal/nocturnal taxa by age
# lik.anc$diurnal <- ifelse(lik.anc$diurnal_R1 + lik.anc$diurnal_R2 > 0.5, 1, 0)
# lik.anc$Time <- tree$data$x[match(lik.anc$node, tree$data$node)]
# # lik.anc$Time2 <- (lik.anc$Time - max(lik.anc$Time))*-1
# lik.anc <- lik.anc[order(lik.anc$Time),]
# 
# ##
# lik.anc$parental.node <- apply(lik.anc, 1, function(x) Ancestors(trpy_n, x["node"], type = "parent"))
# lik.anc$parent.diel <- apply(lik.anc, 1, function(x) lik.anc$diurnal[match(x["parental.node"], lik.anc$node)]) # I should match, not index
# lik.anc$parent.diel[is.na(lik.anc$parent.diel)] <- 0
# lik.anc$switch <- ifelse(lik.anc$parent.diel != lik.anc$diurnal, 1, 0)
# lik.anc$switch[is.na(lik.anc$switch)] <- 0
# 
# # Determine if switch is N->D (1) or D->N (0)
# # lik.anc$switch.ND <- ifelse(lik.anc$switch == 1 & lik.anc$diurnal == 1, 1, 0)
# # lik.anc$switch.DN <- ifelse(lik.anc$switch == 1 & lik.anc$diurnal == 0, 1, 0)
# 
# ### Maybe its better to assign ND, NDN, NDND, NDNDN etc, then covert that to numeric, or whatever
# 
# 
# ## This gives me all the ancestors of a node, in order of descendance?
# # ancestors <- apply(lik.anc, 1, function(x) Ancestors(trpy_n, x["node"], type = "all"))
# ancestors <- apply(lik.anc[lik.anc$node %in% c(1:length(trpy_n$tip.label)),], 1, function(x) Ancestors(trpy_n, x["node"], type = "all"))
# ancestors <- lapply(seq_along(ancestors), function(x) append(c(1:length(trpy_n$tip.label))[[x]], ancestors[[x]]))
# 
# ancestors.diel <- lapply(ancestors, function(x) lapply(x, function(y) lik.anc$diurnal[match(y, lik.anc$node)]))
# 
# # This works, can I simplify it so it works on a vector?
# switch.type <- unlist(lapply(ancestors.diel, function(x) {
#   df <- data.frame(test = unlist(x))
#   df <- df[with(df, c(test[-1]!= test[-nrow(df)], TRUE)),]
#   return(paste(df, collapse = ""))
# }))
# 
# names(switch.type) <- lik.anc[lik.anc$node %in% c(1:length(trpy_n$tip.label)), "node"]
# 
# # lik.anc2 <- lik.anc[lik.anc$node %in% c(1:length(trpy_n$tip.label)),]
# lik.anc$switch.type <- as.character(switch.type[match(lik.anc$node, names(switch.type), nomatch = NA)])
# ## This doesn't make sense, because it counts lineages double.
# ## I think my original way of counting switches, and normalizing them is still Correct?
# ## Or should everything be normalized by the 'normal' cumsum of lineages through time?
# 
# lik.anc$switch.N <- ifelse(is.na(lik.anc$switch.type), 0, ifelse(lik.anc$switch.type == "0", 1, 0))
# lik.anc$switch.ND <- ifelse(is.na(lik.anc$switch.type), 0, ifelse(lik.anc$switch.type == "10", 1, 0))
# lik.anc$switch.NDN <- ifelse(is.na(lik.anc$switch.type), 0, ifelse(lik.anc$switch.type == "010", 1, 0))
# lik.anc$switch.NDND <- ifelse(is.na(lik.anc$switch.type), 0, ifelse(lik.anc$switch.type == "1010", 1, 0))
# lik.anc$switch.NDNDN <- ifelse(is.na(lik.anc$switch.type), 0, ifelse(lik.anc$switch.type == "01010", 1, 0))
# lik.anc$switch.NDNDND <- ifelse(is.na(lik.anc$switch.type), 0, ifelse(lik.anc$switch.type == "101010", 1, 0))
# 
# # Include or not
# node_heights <- nodeHeights(trpy_n)
# lik.anc$node.age <- node_heights[match(lik.anc$node, trpy_n$edge[,2]),1]
# lik.anc$node.age[is.na(lik.anc$node.age)] <- 0
# lik.anc$node.age <- (lik.anc$node.age - max(lik.anc$node.age))*-1
# 
# # Sort by node.age or Time2
# lik.anc <- lik.anc[order(lik.anc$node.age, decreasing = T),]
# 
# # Determine cumulative sums
# lik.anc$Diurnal_Cumsum <- cumsum(ifelse(lik.anc$diurnal > 0.5, 1, 0))
# lik.anc$Nocturnal_Cumsum <- cumsum(ifelse(lik.anc$diurnal < 0.5, 1, 0))
# lik.anc$Lineage_Cumsum <- cumsum(ifelse(lik.anc$node %in% c(1:length(trpy_n$tip.label)), 1, 0))
# lik.anc$Lineage_Cumsum2 <- cumsum(c(rep(1, nrow(lik.anc))))
# lik.anc$Diurnal_ratio <- lik.anc$Diurnal_Cumsum/lik.anc$Lineage_Cumsum # This doesn't work now, because I'm not normalizing by the right thing. If I just use the tips, it ignores lineags that were NDN I think?
# lik.anc$Nocturnal_ratio <- lik.anc$Nocturnal_Cumsum/lik.anc$Lineage_Cumsum # This doesn't work now, because I'm not normalizing by the right thing. If I just use the tips, it ignores lineags that were NDN I think?
# lik.anc$switch2 <- cumsum(lik.anc$switch)
# # lik.anc$switch.ND.Cum <- cumsum(lik.anc$switch.ND)/lik.anc$Lineage_Cumsum
# # lik.anc$switch.DN.Cum <- cumsum(lik.anc$switch.DN)/lik.anc$Lineage_Cumsum
# 
# lik.anc$switch.N.Cum <- cumsum(ifelse(is.na(lik.anc$switch.N), 0, lik.anc$switch.N))
# lik.anc$switch.ND.Cum <- cumsum(ifelse(is.na(lik.anc$switch.ND), 0, lik.anc$switch.ND))
# lik.anc$switch.NDN.Cum <- cumsum(ifelse(is.na(lik.anc$switch.NDN), 0, lik.anc$switch.NDN))
# lik.anc$switch.NDND.Cum <- cumsum(ifelse(is.na(lik.anc$switch.NDND), 0, lik.anc$switch.NDND))
# lik.anc$switch.NDNDN.Cum <- cumsum(ifelse(is.na(lik.anc$switch.NDNDN), 0, lik.anc$switch.NDNDN))
# lik.anc$switch.NDNDND.Cum <- cumsum(ifelse(is.na(lik.anc$switch.NDNDND), 0, lik.anc$switch.NDNDND))

# ###############################################################################################################################################
# ### Make plots ### 
# ###############################################################################################################################################
# 
# # Make LTT plots
# node.data3 <- data.table::melt(lik.anc[,c("node.age", "switch2", "Lineage_Cumsum")], id.vars = "node.age", value.name = "Cummulative_sum")
# ltt.diel.plot <- ggplot(node.data3, aes(x = node.age, y = log(Cummulative_sum), colour = variable)) + geom_line(size = 1.5) + theme_classic() #+ scale_colour_manual(values = c("Diurnal_Cumsum" = "#d6604d", "Nocturnal_Cumsum" = "#4393c3"))
# ltt.diel.plot <- ltt.diel.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + ylab("log(lineages-through-time)") + xlab("Millions of years ago")
# 
# # This is a no good plot, I'm not sure what it even shows? It's a misscounting of 'diurnal' and nocturnal lineages (overcounting), but I don't think the normal LTT plots do this data justice
# node.data3 <- data.table::melt(lik.anc[,c("node.age", "Diurnal_ratio", "Nocturnal_ratio")], id.vars = c("node.age"), value.name = "Cummulative_ratio")
# ltt.diel.ratio.plot <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio, colour = variable)) + geom_line(size = 1.5) + theme_classic() + scale_colour_manual(values = c("Diurnal_ratio" = "#d6604d", "Nocturnal_ratio" = "#4393c3"))
# ltt.diel.ratio.plot <- ltt.diel.ratio.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + xlab("Millions of years ago") #+ ylab("log(lineages-through-time)") 
# 
# # This is now accurate, but it looks strange
# # It's a classic LTT plot, where each tip is also assigned to a switch.type (N, ND, NDN, etc)
# # However, the earliest lineages are not N, but are NDN ???
# node.data3 <- data.table::melt(lik.anc[lik.anc$node %in% c(1:length(trpy_n$tip.label)),c("node.age", "Lineage_Cumsum", "switch.N.Cum", "switch.ND.Cum", "switch.NDN.Cum", "switch.NDND.Cum", "switch.NDNDN.Cum", "switch.NDNDND.Cum")], id.vars = c("node.age", "Lineage_Cumsum"), value.name = "Cummulative_ratio")
# ltt.diel.switch.dir.plot <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio/Lineage_Cumsum, colour = variable)) + geom_line(size = 1.5) + theme_classic() #+ scale_colour_manual(values = c("switch.ND.Cum" = "#d6604d", "switch.DN.Cum" = "#4393c3"))
# ltt.diel.switch.dir.plot <- ltt.diel.switch.dir.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + xlab("Millions of years ago") # + ylab("log(lineages-through-time)") 
# 
# node.data3$CRLC <- node.data3$Cummulative_ratio#/node.data3$Lineage_Cumsum
# node.data4 <- node.data3 %>% group_by(node.age, variable) %>% summarise(n = sum(CRLC)) %>% mutate(percentage = n / sum(n))
# # node.data4$variable <- factor(node.data4$variable, levels = c("switch.NDN.Cum", "switch.N.Cum", "switch.ND.Cum", "switch.NDND.Cum", "switch.NDNDN.Cum"))
# switch_histo <- ggplot(node.data4, aes(x = node.age, y = percentage, fill = variable)) + theme_classic() + geom_area() + scale_x_reverse() + xlab("Millions of years ago") + scale_fill_manual(values = c("blue4", "red4", "blue3", "red3", "blue1", "red1"))
# 
# 
# # Plots
# # switch_histo <- ggplot(lik.anc, aes(x = node.age, fill = as.character(switch.type))) + geom_histogram(position = "fill", bins = 100) + scale_x_reverse()
# switch_ratio <- ggplot(lik.anc, aes(x = node.age, y = switch2/Lineage_Cumsum2)) + geom_line() + theme_classic() + scale_x_reverse() # Lineage_Cumsum2 would be the number of times it was possible to switch? Or something like that?
# ltt.switch.plot <- ggplot(lik.anc, aes(x = node.age, y = log(switch2))) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + ylab("log(switches-through-time)") + xlab("Millions of years ago")
# 
# 
# # ggplot(lik.anc, aes(y = log(switch2), x = log(Lineage_Cumsum))) + geom_point() + theme_bw()
# 
# # Save plots
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ltt_multiswitch_ratio", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 6, height = 9)
# (ltt.diel.plot / switch_ratio / ltt.diel.switch.dir.plot / switch_histo) & xlim(c(430,0))
# dev.off()
# 
# source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")
# 
# numb_switch_tree <- ggtree(trpy_n, layout = "circular", size = 0.25) %<+% lik.anc3 + aes(color = switch.type) + geom_tippoint(aes(color = switch.type), shape = 16) + scale_color_manual(values = c("grey75", "black", "blue", "grey50", "green", "red"))

# pdf(file = paste("outs/Figures/fish_phylogeny_diel_highswitchlineages", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), height = 60, width = 60)
# numb_switch_tree
# dev.off()
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ltt_switch_ratio", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 6, height = 4)
# switch_histo / switch_ratio / ltt.diel.switch.dir.plot
# dev.off()
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ltt", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 6, height = 4)
# ltt.diel.plot + ltt.diel.ratio.plot + plot_layout(nrow = 2)
# dev.off()
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_switch_stt", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 4.25, height = 4)
# ltt.switch.plot
# dev.off()
# 
# ## Plots the ratio between the cumulative switchs (diurnal to nocturnal and vice versa), and cumulative lineages
# 
# source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/scripts/gggeo_scale.R")
# 
# 
# cum.ratio.plot <- ggplot(lik.anc, aes(x = node.age, y = switch2/Lineage_Cumsum)) + geom_line(size = 1.5) + theme_classic()
# cum.ratio.plot <- gggeo_scale(cum.ratio.plot, pos = "top") + scale_x_reverse() + xlim(405.2853,0)


# # Extract the data from the plot from the PNAS paper figur 1 (https://www.pnas.org/content/114/22/5653) using https://automeris.io/WebPlotDigitizer/
# geo_frag <- read.csv("~/Downloads/Default Dataset.csv")
# colnames(geo_frag) <- c("Time", "frag_index") # Time is millions of years ago
# # geo_frag$root <- (geo_frag$Time - max(geo_frag$Time))*-1
# geo.frag.plot <- ggplot(geo_frag, aes(x = Time, y = frag_index)) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + xlim(405.2853,0)
# geo.frag.plot <- geo.frag.plot + theme_bw() + scale_x_reverse(limits = c(408.341, 0), breaks = scales::pretty_breaks(n = 20))
# # Plot them together to see correlation
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_stt_geofrag_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 10, height = 10)
# cum.ratio.plot + plot_layout(nrow = 2) + geo.frag.plot
# dev.off()
# 
# 
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_ltt_all_plots_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 15, height = 7)
# cum.ratio.plot + ltt.diel.plot + ltt.diel.switch.dir.plot + ltt.diel.ratio.plot + plot_layout(nrow = 2)
# dev.off()































# ggplot(test.df, aes(x = as.numeric(node.age), y = log(reconstruction) - log(means))) + geom_line() + theme_bw()
# 
# 
# df <- data.frame(diurnal = as.numeric(as.factor(simulation[,555])))
# 
# p1 <- ggtree(trpy_n, layout = "circular") %<+% data.frame(diurnal = simulation[,2]) + aes(color = simulation[,55]) + scale_color_manual(values = c("red", "blue")) #+ geom_tippoint(aes(color = simulation[,1]), shape = 16, size = 1.5)
# p2 <- ggtree(trpy_n, layout = "circular") %<+% data.frame(diurnal = simulation[,2]) + aes(color = simulation[,99]) + scale_color_manual(values = c("red", "blue")) #+ geom_tippoint(aes(color = simulation[,1]), shape = 16, size = 1.5)
# 
# p1 + p2
# 
# confidence_interval <- function(vector, interval) {
#   # Standard deviation of sample
#   vec_sd <- sd(vector)
#   # Sample size
#   n <- length(vector)
#   # Mean of sample
#   vec_mean <- mean(vector)
#   # Error according to t distribution
#   error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
#   # Confidence interval as a vector
#   result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
#   return(result)
# }
# 
# 
#                       
# 
# plot <- ggplot(test.df, aes(x = node.age*-1, y = (means))) + geom_line(size = 1.5) + theme_classic() + ylim(c(0,1))
# 
# plot <- plot + ggplot(test.df, aes(x = node.age*-1, y = upper)) + geom_line(colour = "red", size = 1.5) + theme_classic() + ylim(c(0,1))
# 
# plot <- plot + ggplot(test.df, aes(x = node.age*-1, y = lower)) + geom_line(colour = "blue", size = 1.5) + theme_classic() + ylim(c(0,1))
# 
# 
# 
# ggplot(test.df, aes(x=node.age*-1, y = value)) + stat_summary(geom="ribbon", fun.ymin="min", fun.ymax="max", aes(fill=type), alpha=0.3) + theme_classic()
