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

# Load trees
trpy_n_fish <- loadTree(return = "tree", dataset = "fish", subset = "all")
trpy_n_all <- loadTree(return = "tree", dataset = "AllGroups", subset = "all")


anc_states_fish <- readRDS(file = paste("fish", "diel_ancestral_states", "all", Ntip(trpy_n_fish), "species.rds", sep = "_"))
anc_states_all <- readRDS(file = paste("AllGroups", "diel_ancestral_states", "all", Ntip(trpy_n_all), "species.rds", sep = "_"))

switch.plots.full <- list()
switch.plots.full[[1]] <- switchRatio(ancestral_states = anc_states_all, phylo_tree = trpy_n_all, node.age.cutoff = 0.02)
switch.plots.full[[2]] <- switchRatio(ancestral_states = anc_states_fish, phylo_tree = trpy_n_fish, node.age.cutoff = 0.02)

## Make the new_states
# Copy node info
trpy_n_all$node.label <- paste0("n", (1:trpy_n_all$Nnode) + Ntip(trpy_n_all))

## MRCAs for different groups
mrca <- list()
mrca$amphibians <- getMRCA(phy = trpy_n_all, tip = c("Allophryne_ruthveni", "Dermophis_mexicanus")) # Allophryne ruthveni (Anura), Dermophis mexicanus (Gymnophiona)
mrca$mammals <- getMRCA(phy = trpy_n_all, tip = c("Ornithorhynchus_anatinus", "Acinonyx_jubatus")) # Ornithorhynchus anatinus (Monotremata),	Acinonyx jubatus (Carnivora)
# mrca$lepidosaurs <- getMRCA(phy = trpy_n_all, tip = c("Anguis_fragilis", "Gonatodes_eladioi")) # Anguis fragilis (Anguimorpha), Gonatodes eladioi (Gekkota)
# mrca$testudines <- getMRCA(phy = trpy_n_all, tip = c("Chelodina_longicollis", "Cyclanorbis_elegans")) # 	Chelodina_longicollis	 (Chelidae), Cyclanorbis elegans (Trionychidae)
# mrca$crocodylia <- getMRCA(phy = trpy_n_all, tip = c("Gavialis_gangeticus", "Alligator_mississippiensis")) # Gavialis gangeticus (Gavialis), Alligator mississippiensis	(Alligatorinae)
# mrca$aves <- getMRCA(phy = trpy_n_all, tip = c("Nothocrax_urumutum", "Acrocephalus_arundinaceus")) # Nothocrax urumutum (Galliformes), Acrocephalus arundinaceus (Passeriformes)
mrca$sauropsida <- getMRCA(phy = trpy_n_all, tip = c("Phelsuma_madagascariensis", "Coracina_novaehollandiae"))

droptrees_all <- lapply(mrca, function(x) extract.clade(phy = trpy_n_all, node = x))
names(droptrees_all) <- names(mrca)

# Make the node index
node_index <- list()
node_index <- lapply(droptrees_all, function(x) data.frame(old.node = c( c(1:Ntip(trpy_n_all))[trpy_n_all$tip.label %in% x$tip.label], as.numeric(str_sub(x$node.label, start = 2, end = 10)) ), new.node = c(1:(Ntip(x) + Nnode(x)))) ) 

# subset the anc_states data
new_states <- list()
for (i in 1:length(node_index)) {
  print(paste("running new_states for", names(mrca)[[i]], sep = " "))
  new_states[[i]] <- anc_states_all[1:4]
  # This doesn't work correctly yet, need to also include the tips
  new_states[[i]]$lik.anc <- new_states[[i]]$lik.anc[node_index[[i]]$old.node,]
  new_states[[i]][[2]] <- node_index[[i]]$new.node
  
  # re-make new_states data
  new_states[[i]] <- calculateStateTransitions(ancestral_states = new_states[[i]], phylo_tree = droptrees_all[[i]], ancestor = round(anc_states_all$recon_states[mrca[[i]]]))
  new_states[[i]] <- calculateLinTransHist(ancestral_states = new_states[[i]], phylo_tree = droptrees_all[[i]])
  new_states[[i]] <- returnCumSums(ancestral_states = new_states[[i]], phylo_tree = droptrees_all[[i]])
}


switch.plots.all <- list()
switch.plots.all <- lapply(seq_along(new_states), function(x) switchRatio(ancestral_states = new_states[[x]], phylo_tree = droptrees_all[[x]], node.age.cutoff = 0.02))


## For fish
# Copy node info
trpy_n_fish$node.label <- paste0("n", (1:trpy_n_fish$Nnode) + Ntip(trpy_n_fish))

mrca <- list()
mrca$cartilagenous <- getMRCA(phy = trpy_n_fish, tip = c("Tetronarce_californica", "Squatina_squatina"))
mrca$bonyfish <- getMRCA(phy = trpy_n_fish, tip = c("Lepisosteus_osseus", "Lutjanus_fulvus"))

droptrees_fish <- lapply(mrca, function(x) extract.clade(phy = trpy_n_fish, node = x))
names(droptrees_fish) <- names(mrca)

# Make the node index
node_index <- list()
node_index <- lapply(droptrees_fish, function(x) data.frame(old.node = c( c(1:Ntip(trpy_n_fish))[trpy_n_fish$tip.label %in% x$tip.label], as.numeric(str_sub(x$node.label, start = 2, end = 10)) ), new.node = c(1:(Ntip(x) + Nnode(x)))) ) 

# subset the anc_states data
new_states <- list()
for (i in 1:length(node_index)) {
  print(paste("running new_states for", names(mrca)[[i]], sep = " "))
  new_states[[i]] <- anc_states_fish[1:4]
  # This doesn't work correctly yet, need to also include the tips
  new_states[[i]]$lik.anc <- new_states[[i]]$lik.anc[node_index[[i]]$old.node,]
  new_states[[i]][[2]] <- node_index[[i]]$new.node
  
  # re-make new_states data
  new_states[[i]] <- calculateStateTransitions(ancestral_states = new_states[[i]], phylo_tree = droptrees_fish[[i]], ancestor = round(anc_states_fish$recon_states[mrca[[i]]]))
  new_states[[i]] <- calculateLinTransHist(ancestral_states = new_states[[i]], phylo_tree = droptrees_fish[[i]])
  new_states[[i]] <- returnCumSums(ancestral_states = new_states[[i]], phylo_tree = droptrees_fish[[i]])
}


switch.plots.fish <- list()
switch.plots.fish <- lapply(seq_along(new_states), function(x) switchRatio(ancestral_states = new_states[[x]], phylo_tree = droptrees_fish[[x]], node.age.cutoff = 0.02))



## OK lets put them together
## This extracts the plot data, adds grouping info, then rbinds it together before plotting

switch.plots <- append(switch.plots.all, switch.plots.fish)
switch.plots <- append(switch.plots, switch.plots.full)

switch.data <- lapply(switch.plots, function(x) as.data.frame(x$data))

switch.data <- lapply(seq_along(switch.data), function(x) {
  switch.data[[x]]$group <- c(names(droptrees_all), names(droptrees_fish), "all_tetrapods", "all_fish" )[x]
  return(switch.data[[x]])
})

switch.data <- Reduce(rbind, switch.data)


ggplot(switch.data, aes(x = node.age, y = transition_cumsum, group = group, colour = group)) + geom_line(size = 2) + scale_x_reverse() + theme_classic()


##### This is the code to add them together

