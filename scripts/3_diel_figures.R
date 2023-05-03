library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)
library(geiger)
library(phytools)
library(rfishbase)
library(xlsx)
library(RColorBrewer)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/Fish_sleep_functions.R")

index_list <- list()
index_list[[1]] <- c("all", "only_highqual", "only_cartilaginous", "only_ingroup")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

for (i in 1:length(index_list)) {
  dataset_variable <- names(index_list)[[i]]
  
  for (j in 1:length(index_list[[i]])) {
    name_variable <- index_list[[i]][[j]]
    
    ## Load in the tree
    trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
    trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)
    
    ## Load the best model, which is the HMM 2 state 2 rate model
    models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
    
    ## I think I will always take the 2 rate model, 3 is too hard to comphrehend
    
    model <- models$HMM_2state_2rate_marg
    
    
  }
}
    
    
    
    
    
    
    
## Load files for making figures

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
trait.data <- readRDS("trait_data.rds")
tr.calibrated <- readRDS("calibrated_phylo.rds")

name_variable <- "all"

# # Remove low quality data or not
# trait.data <- trait.data[trait.data$confidence > 1,]
# tr.calibrated <- keep.tip(tr.calibrated, tip = trait.data$species)
# name_variable <- "only_highqual"

# # Remove cartilaginous fish (just actinopterygii)
# # Common ancestor of Lepisosteus_osseus (Order: Lepisosteiformes) Lepidosiren_paradoxa (lepidosiren) and Lutjanus_fulvus (Order: Perciformes)
# node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Lepidosiren_paradoxa", "Lutjanus_fulvus"))
# tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
# trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
# name_variable <- "only_ingroup"

# Keep only cartilaginous fish (no actinopterygii)
node_of_interest <- getMRCA(tr.calibrated, tip = c("Rhizoprionodon_terraenovae", "Rhynchobatus_djiddensis"))
tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
name_variable <- "only_cartilaginous"

################################################################################################################################################
### Make phylogenetic tree plots showing activity patterns ###
################################################################################################################################################

# Can I make a phylo tree, but add bars to indicate diel activity pattern (so that I can make them longer and skinnier) - like in Fig 2 here (https://www.pnas.org/doi/10.1073/pnas.1216063110)
# Can then use the same for adding orders, maybe

diel.plot <- ggtree(tr.calibrated, layout = "circular") %<+% trait.data[,c("tips", "diel2", "order")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(tr.calibrated$tip.label),], aes(y=y, x=x+15, fill = diel2), width = 30, inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = c("#abdda4", "#d7191c", "#2c7bb6", "#fdae61"))
diel.plot.labelled <- diel.plot + geom_tiplab(color = "black", size = 1.5, offset = 30)
# # Make a plot showing species names and diel activity
# diel.plot <- ggtree(tr.calibrated, layout = "circular") %<+% trait.data[,c("tips", "diel2", "order")] + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diel2), shape = 16, size = 0.5) + scale_color_manual(values = c("yellow", "red", "blue", "green"))

pdf(file = paste("outs/Figures/fish_phylogeny_diel_", length(tr.calibrated$tip.label), "_species.pdf", sep = ""), width = 40, height = 40)
diel.plot.labelled
dev.off()

# Make a plot showing diel activity and orders
# diel.plot <- ggtree(tr.calibrated, layout = "circular") %<+% resolved_names[,c("tips", "diel", "order")] + geom_tippoint(aes(color = diel), shape = 16, size = 1.5) + scale_color_manual(values = c("red", "blue", "yellow", "green3", "green1", "green4", "green2", "white"))

# Annotate clades?
orders <- unique(resolved_names$order)[!(is.na(unique(resolved_names$order)))]
# orders <- orders[!(grepl("Perciformes", orders))]
# orders <- orders[!(orders %in% c("Scorpaeniformes", "Lampriformes", "Osmeriformes", "Beryciformes", "Scombriformes", "Eupercaria/misc", "Acanthuriformes", "Acropomatiformes"))]
# orders <- orders[!(orders %in% c("Perciformes/Percoidei", "Perciformes/Notothenioidei", "Perciformes", "Scorpaeniformes", "Lampriformes", "Osmeriformes", "Beryciformes", "Scombriformes", "Eupercaria/misc", "Acanthuriformes", "Perciformes/Serranoidei", "Acropomatiformes", "Perciformes/Uranoscopoidei", "Perciformes/Zoarcoidei", "Perciformes/Gasterosteoidei", "Perciformes/Cottoidei", "Perciformes/Notothenoidei"))]

# Or just a few orders
# orders <- c("Anguilliformes", "Cypriniformes", "Characiformes", "Siluriformes", "Kurtiformes", "Gobiiformes", "Syngnathiformes", "Pleuronectiformes", "Cichliformes", "Perciformes/Scorpaenoidei", "Tetraodontiformes", "Carcharhiniformes", "Holocentriformes", "Blenniiformes")

my_colors <- hue_pal()(length(orders))

# Add highlights and labels using a for loop for each order

diel.plot.orders <- diel.plot

for (i in 1:length(orders)) {
  tips <- resolved_names$tips[resolved_names$order == orders[[i]] & !(is.na(resolved_names$genus))]
  tips <- tips[!(is.na(tips))]
  anchor.point <- round(median(1:length(tips)))
  anchor.point <- tips[anchor.point]
  tips <- tips[!(is.na(tips))]
  tips <- match(tips, tr.calibrated$tip.label)
  tips <- tips[!(is.na(tips))]
  
  offset = 25
  alpha = 0.25
  
  if (orders[[i]] %in% ("test")) {
    diel.plot.orders <- diel.plot.orders
  } else {
    if(length(tips) > 1) {
      node <- getMRCA(tr.calibrated, tips)
      diel.plot.orders <- diel.plot.orders + geom_hilight(node = node, fill = my_colors[[i]], alpha = alpha)
      # angle = diel.plot.orders$data$angle[match(node, diel.plot.orders$data$node)]
      angle <- mean(diel.plot.orders$data$angle[diel.plot.orders$data$node %in% tips])
      if(between(angle, 90, 270)) {
        angle <- angle + 180
        hjust <- 1
      } else {
        angle <- angle
        hjust <- 0
      }
      diel.plot.orders <- diel.plot.orders + geom_cladelabel(node = node, label = orders[[i]], color = my_colors[[i]], offset = offset, align = F, angle = angle, hjust = hjust, barsize = 1.5)
    } 
    
    if(length(tips) == 1) {
      diel.plot.orders <- diel.plot.orders + geom_hilight(node = tips, fill = my_colors[[i]], alpha = alpha)
      # angle = diel.plot.orders$data$angle[match(tips, diel.plot.orders$data$node)]
      angle <- mean(diel.plot.orders$data$angle[diel.plot.orders$data$node %in% tips])
      if(between(angle, 90, 270)) {
        angle <- angle + 180
        hjust <- 1
      } else {
        angle <- angle
        hjust <- 0
      }
      diel.plot.orders <- diel.plot.orders + geom_cladelab(node = tips, label = orders[[i]], color = my_colors[[i]], offset = offset, align = F, angle = angle, hjust = hjust, barsize = 1.5)
    }
  }
}
rm(tips, anchor.point, node, angle, hjust)

pdf(file = paste("outs/Figures/fish_phylogeny_diel_orders_", length(tr.calibrated$tip.label), "_species.pdf", sep = ""), width = 40, height = 40)
diel.plot.orders + xlim(0,600)
dev.off()

# To save the subset
pdf(file = paste("outs/Figures/fish_phylogeny_diel_orders_subsetLabelled_", length(tr.calibrated$tip.label), "_species.pdf", sep = ""), width = 15, height = 15)
diel.plot.orders + xlim(0,600)
dev.off()

################################################################################################################################################
### Plot Reconstructed ancestral states for Diurnal / Nocturnal ###
################################################################################################################################################


# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]

# Load most likely model from script #2
# best_fit_model <- readRDS("best_fit_model_all_species.rds")
best_fit_model <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

# Extract the liklihood for each node and make it into a ggplotable data frame

lik.anc <- as.data.frame(rbind(best_fit_model$tip.states, best_fit_model$states))

if (ncol(lik.anc) == 4) {
  colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
  lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
  
  lik.anc$noc.sum <- rowSums(lik.anc[,c(2,4)])
  lik.anc$di.sum <- rowSums(lik.anc[,c(1,3)])
  
  lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
  lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])
} else {
  colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2", "diurnal_R3", "nocturnal_R3")
  lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
  
  lik.anc$noc.sum <- rowSums(lik.anc[,c(2,4,6)])
  lik.anc$di.sum <- rowSums(lik.anc[,c(1,3,5)])
  
  lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
  lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])
  lik.anc$R3.sum <- rowSums(lik.anc[,c(5,6)])
}


# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
ancestral_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
# ancestral_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_gradient(low = "#2c7bb6", high = "#ffffbf")

ancestral_plot_rect <- ggtree(trpy_n) %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")

rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")

# Because it is a 2 state 2 rate plot
dir1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R1) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1)
dir2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R2) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
nocr1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R1) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)
nocr2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R2) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)

rate_plots <- dir1 + dir2 + nocr1 + nocr2 + plot_layout(nrow = 2)

if (ncol(lik.anc) == 4) {
  dir3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R3) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "YlOrBr", direction = 1)
  nocr3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R3) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "BuPu", direction = 1)
  rate_plots <- dir1 + dir2 + dir3 + nocr1 + nocr2 + nocr3 + plot_layout(nrow = 2)
}




# Add markers for those nodes where the reconstruction is confident (80% or higher)
lik.anc$cutoff <- ifelse(lik.anc$noc.sum > 0.8 | lik.anc$di.sum > 0.8, TRUE, FALSE)

ancestral_plot$data$cutoff <- ifelse(!(ancestral_plot$data$isTip) & lik.anc$cutoff, TRUE, FALSE)

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_", length(trpy_n$tip.label), "_species_node_cutoff.pdf", sep = ""), width = 40, height = 40)
ancestral_plot + geom_text2(aes(subset=cutoff, label = node), colour = "black")
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_", length(trpy_n$tip.label), "_species_cutoff.pdf", sep = ""), width = 20, height = 20)
ancestral_plot + geom_point(data = subset(ancestral_plot$data, subset = cutoff), aes(x = x, y = y), size = 1.5, colour = "black")
dev.off()

#### Save the plots

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 40, height = 40)
ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5)
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_", length(trpy_n$tip.label), "_species_nolab.pdf", sep = ""), width = 20, height = 20)
ancestral_plot
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_rec_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 25, height = 80)
ancestral_plot_rect + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5)
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_rec_", length(trpy_n$tip.label), "_species_nolab.pdf", sep = ""), width = 25, height = 80)
ancestral_plot_rect
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_", length(trpy_n$tip.label), "_species_nolab_small.pdf", sep = ""), width = 10, height = 10)
ancestral_plot
dev.off()

# Add order highlights to ancestral.plot

# ancestral.plot <- ggtree(trpy, layout = "circular") %<+% node.data + aes(color = diurnal) + scale_color_distiller(palette = "RdBu") + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu")
ancestral.plot <- ancestral_plot

orders <- c("Anguilliformes", "Cypriniformes", "Characiformes", "Siluriformes", "Kurtiformes", "Gobiiformes", "Syngnathiformes", "Pleuronectiformes", "Cichliformes", "Perciformes/Scorpaenoidei", "Tetraodontiformes", "Carcharhiniformes", "Holocentriformes", "Blenniiformes")

my_colors <- hue_pal()(length(orders))

# Add labels using a for loop for each order

for (i in 1:length(orders)) {
  tips <- resolved_names$tips[resolved_names$order == orders[[i]] & !(is.na(resolved_names$genus))]
  tips <- tips[!(is.na(tips))]
  anchor.point <- round(median(1:length(tips)))
  anchor.point <- tips[anchor.point]
  tips <- tips[!(is.na(tips))]
  tips <- match(tips, trpy_n$tip.label)
  tips <- tips[!(is.na(tips))]
  
  offset = 5
  alpha = 0.25
  
  if (orders[[i]] %in% ("test")) {
    ancestral.plot <- ancestral.plot
  } else {
    if(length(tips) > 1) {
      node <- getMRCA(trpy_n, tips)
      ancestral.plot <- ancestral.plot + geom_hilight(node = node, fill = my_colors[[i]], alpha = alpha)
      # angle = ancestral.plot$data$angle[match(node, ancestral.plot$data$node)]
      angle <- mean(ancestral.plot$data$angle[ancestral.plot$data$node %in% tips])
      if(between(angle, 90, 270)) {
        angle <- angle + 180
        hjust <- 1
      } else {
        angle <- angle
        hjust <- 0
      }
      ancestral.plot <- ancestral.plot + geom_cladelabel(node = node, label = orders[[i]], color = my_colors[[i]], offset = offset, align = F, angle = angle, hjust = hjust, barsize = 1.5)
    } 
    
    if(length(tips) == 1) {
      ancestral.plot <- ancestral.plot + geom_hilight(node = tips, fill = my_colors[[i]], alpha = alpha)
      # angle = ancestral.plot$data$angle[match(tips, ancestral.plot$data$node)]
      angle <- mean(ancestral.plot$data$angle[ancestral.plot$data$node %in% tips])
      if(between(angle, 90, 270)) {
        angle <- angle + 180
        hjust <- 1
      } else {
        angle <- angle
        hjust <- 0
      }
      ancestral.plot <- ancestral.plot + geom_cladelab(node = tips, label = orders[[i]], color = my_colors[[i]], offset = offset, align = F, angle = angle, hjust = hjust, barsize = 1.5)
    }
  }
}
rm(tips, anchor.point, node, angle, hjust)

pdf(file = paste("outs/Figures/fish_phylogeny_diel_orders_ancestral_", length(trpy_n$tip.label), "_species.pdf", sep = ""), width = 20, height = 20)
ancestral.plot + xlim(0,500)
dev.off()


