library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)


# setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")
setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")

# source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/scripts/Trait_rate_matrix_figure_script.R")

## Load files

resolved_names <- readRDS("resolved_names_AllGroups.rds")
tr.calibrated <- readRDS("tr_tree_calibrated_AllGroups.rds")
trait.data <- readRDS("trait_data_AllGroups.rds")

name_variable <- "AllGroups"

###############################################################################################################################################
### Now subset the tree based on which traits you want to test ### 
###############################################################################################################################################

# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]

trpy_n$node.label <- Ntip(trpy_n):Nnode(trpy_n)

###############################################################################################################################################
### Run Hidden Rates Models ### 
###############################################################################################################################################

# Run corHMM for just diel with 1 rat category (Markov model)
MK_2state <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = "ARD")

# Run corHMM for just diel, with 2 or 3 rate categories (large improvement for 2, diminishing returns for 3)
HMM_2state_2rate <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = "ARD", get.tip.states = TRUE)

# Save out 2 state and HMM with 2 rates
standard_tests <- list(c(MK_2state, HMM_2state_2rate))
saveRDS(standard_tests, file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

HMM_2state_3rate <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 3, model = "ARD", get.tip.states = TRUE)
standard_tests[[3]] <- HMM_2state_3rate
saveRDS(standard_tests, file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))


# # Look at models
# plotMKmodel(MK_2state) # High Diurnal -> Nocturnal
# plotMKmodel(HMM_2state_2rate) # 1 rate is high transition rate (both ways), other is low both ways, with slightly higher rate to go from high to low transition rate
# 
# ###############################################################################################################################################
# 
# # Re-load the best fit model
# # HMM_2state_2rate <- readRDS(file = paste("best_fit_model_", length(trpy_n$tip.label), "_species.rds", sep = ""))
# HMM_2state_2rate <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
# MK_2state <- readRDS(file = paste("MK_2state_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
# 
# ###############################################################################################################################################
# ### Plot the Hidden rate ancestral reconstructions (of rate and state) ### 
# ###############################################################################################################################################
# model <- MK_2state
# 
# lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
# colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1") #, "diurnal_R2", "nocturnal_R2")
# lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
# 
# # lik.anc$noc.sum <- rowSums(lik.anc[,c(2,4)])
# # lik.anc$di.sum <- rowSums(lik.anc[,c(1,3)])
# 
# # lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
# # lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])
# 
# # Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
# diel_allGroups <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
# diel_allGroups_rect <- ggtree(trpy_n) %<+% lik.anc + aes(color = diurnal_R1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
# 
# 
# ## OK, now collapse some nodes
# 
# mrca_amphibians <- getMRCA(phy = trpy_n, tip = c("Allophryne_ruthveni", "Dermophis_mexicanus")) # Allophryne ruthveni (Anura), Dermophis mexicanus (Gymnophiona)
# mrca_mammals <- getMRCA(phy = trpy_n, tip = c("Ornithorhynchus_anatinus", "Acinonyx_jubatus")) # Ornithorhynchus anatinus (Monotremata),	Acinonyx jubatus (Carnivora)
# mrca_lepidosaurs <- getMRCA(phy = trpy_n, tip = c("Anguis_fragilis", "Gonatodes_eladioi")) # Anguis fragilis (Anguimorpha), Gonatodes eladioi (Gekkota)
# mrca_testudines <- getMRCA(phy = trpy_n, tip = c("Chelodina_longicollis", "Cyclanorbis_elegans")) # 	Chelodina_longicollis	 (Chelidae), Cyclanorbis elegans (Trionychidae)
# mrca_crocodylia <- getMRCA(phy = trpy_n, tip = c("Gavialis_gangeticus", "Alligator_mississippiensis")) # Gavialis gangeticus (Gavialis), Alligator mississippiensis	(Alligatorinae)
# mrca_aves <- getMRCA(phy = trpy_n, tip = c("Nothocrax_urumutum", "Acrocephalus_arundinaceus")) # Nothocrax urumutum (Galliformes), Acrocephalus arundinaceus (Passeriformes)
# mrca_cartilagenous <- getMRCA(phy = trpy_n, tip = c("Tetronarce_californica", "Squatina_squatina"))
# mrca_bonyfish <- getMRCA(phy = trpy_n, tip = c("Lepisosteus_osseus", "Lutjanus_fulvus"))
# 
# # For tip sizes
# log(table(all_data$group))*3
# 
# p2 <- diel_allGroups_rect
# p2 <- p2 %>% collapse(node = mrca_bonyfish) + geom_point2(aes(subset = (node==mrca_bonyfish)), shape=0, size=20)
# p2 <- p2 %>% collapse(node = mrca_mammals) + geom_point2(aes(subset = (node==mrca_mammals)), shape=1, size=25)
# p2 <- p2 %>% collapse(node = mrca_aves) + geom_point2(aes(subset = (node==mrca_aves)), shape=2, size=19)
# p2 <- p2 %>% collapse(node = mrca_lepidosaurs) + geom_point2(aes(subset = (node==mrca_lepidosaurs)), shape=3, size=18)
# p2 <- p2 %>% collapse(node = mrca_amphibians) + geom_point2(aes(subset = (node==mrca_amphibians)), shape=4, size=18)
# p2 <- p2 %>% collapse(node = mrca_crocodylia) + geom_point2(aes(subset = (node==mrca_crocodylia)), shape=5, size=9)
# p2 <- p2 %>% collapse(node = mrca_cartilagenous) + geom_point2(aes(subset = (node==mrca_cartilagenous)), shape=6, size=5)
# p2 <- p2 %>% collapse(node = mrca_testudines) + geom_point2(aes(subset = (node==mrca_testudines)), shape=7, size=13)
# 
# pdf("outs/Figures/collapsed_all_groups_reconstruction2.pdf", height = 7.5, width = 7.5)
# p2
# dev.off()
# 
# major_groups_plot <- p2 + geom_tiplab() + geom_nodelab(size = 5, colour = "black")
# 
# # diel_allGroups + geom_nodelab(size = 2, colour = "black")
# 
# ## Not sure why the below doesn't work
# 
# order_numb <- data.frame(names = unique(trait.data_n$order)[!is.na(unique(trait.data_n$order))], numb = as.vector(table(trait.data_n$order)))
# order_numb <- order_numb[order_numb$numb > 15,]
# mrca <- list()
# for (i in 1:nrow(order_numb)) {
#   if (order_numb$numb[i] > 15) {
#     mrca[[i]] <- getMRCA(phy = trpy_n, tip = trait.data_n$tips[trait.data_n$order %in% c(order_numb$names[i])])
#   } else {
#     mrca[[i]] <- "none"
#   }
# }
# names(mrca) <- order_numb$names
# mrca_nodes <- unlist(mrca[mrca != "none"])
# 
# mrca_nodes <- sort(mrca_nodes, decreasing = TRUE)
# p2 <- diel_allGroups_rect
# 
# for (i in 1:length(mrca_nodes[1:25])) {
#   p2 <- p2 %>% collapse(node = mrca_nodes[[i]]) + geom_point2(aes(subset = (node==mrca_nodes[[i]])), shape=20, size=20)
# }
# 
# p2
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")
# 
# 
# dir1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R1) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1)
# dir2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R2) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
# nocr1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R1) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)
# nocr2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R2) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)
# dir3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R3) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "YlOrBr", direction = 1)
# nocr3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R3) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "BuPu", direction = 1)
# 
# two_rate_plots <- dir1 + dir2 + nocr1 + nocr2 + plot_layout(nrow = 2)
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_diel_rates", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width =20, height = 20)
# two_rate_plots
# dev.off()
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_rates", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 20, height = 20)
# rate_plot
# dev.off()



