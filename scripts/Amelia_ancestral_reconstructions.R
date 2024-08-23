# Section 0: Import and run packages -------------------------------------------------
# For retreiving data from the open tree of life
library(rotl)
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

# Set the working directory and source the functions (not used yet)
setwd(here())

#with this script we will add cetaceans to the rest of artiodactyla and do ancestral trait reconstruction
#as well as modelling the evolution of diel patterns, to see if it's similar or different
#to the results seen with only cetaceans

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")


# Section 1: Import data  -----------------------------------------------

#we want this script to be flexible enough to take any model input
#only need to change the model_results and file_name to create plots of new data

#currently I have max_clade_crep data for: artio max_crep, artio max_dinoc
#to do: cetacean max_crep, cetacean max_dinoc, artio w/out cetaceans max_crep, artio w/out cetaceans max_dinoc
all_model_results <- readRDS(here("artiodactyla_max_clade_cred_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))
#copy and paste first half of filename here (leave out the models)
file_name <- "artiodactyla_max_clade_cred_max_crep_traits"

#separate the results by the model types we want to use (ER, SYM, ARD, bridge_only)
#uncomment the model you want to plot

#model_results <- all_model_results$ER_model
#model_name <- "ER"

#model_results <- all_model_results$SYM_model
#model_name <- "SYM"

#model_results <- all_model_results$ARD_model
#model_name <- "ARD"

model_results <- all_model_results$bridge_only
model_name <- "bridge_only"

# Section 1: Plotting ancestral reconstruction from corHMM model  --------

#from the model results file, tip states describes the trait states at the tips, states describes the trait states at the nodes
lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("cathemeral", "diurnal", "nocturnal")
phylo_tree <- model_results$phy
#associate each of these species and their trait states with its node
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

#plot the ancestral reconstruction, displaying each of the three trait states (cathemeral, diurnal, nocturnal)
ancestral_plot_di <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
ancestral_plot_di
ancestral_plot_noc <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
ancestral_plot_noc
ancestral_plot_cath <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
ancestral_plot_cath

#create the name of the file by pasting together ancestral recon, the diel state and the file_name 
png(paste("C:/Users/ameli/OneDrive/Documents/R_projects/New_ancestral_recon/", "ancestral_recon_diurnal_", file_name, "_", model_name, ".png", sep = ""), width=17,height=16, units="cm",res=1200)
ancestral_plot_di
dev.off()

png(paste("C:/Users/ameli/OneDrive/Documents/R_projects/New_ancestral_recon/", "ancestral_recon_nocturnal_", file_name, "_", model_name, ".png", sep = ""), width=17,height=16,units="cm",res=1200)
ancestral_plot_noc
dev.off()

png(paste("C:/Users/ameli/OneDrive/Documents/R_projects/New_ancestral_recon/", "ancestral_recon_cathemeral_", file_name, "_", model_name,  ".png", sep = ""), width=17,height=16,units="cm",res=1200)
ancestral_plot_cath
dev.off()


# Alternative: Using ancRECON to retrieve the internal node states -------


#index from 1 to the max value from the rate index matrix (ignoring NA values with na.rm)
#so for ER the max value is 1 (all rates are equal), for ARD this is 6 (six different rates)
#na.omit removes na values from the vector containing the model solutions (the transition rates)
#then we index this vector of transition rates, by the vector of rate indices (ie 1,1,1,1,1,1 for ER, 1,2,3,4,5,6 for ARD)
#this creates a numeric vector of transition rates in the correct order (ie for ARD matrix index 1 = transition rate 1, matrix index 2 = transition rate 2, etc)
p <- sapply(1:max(ER_results$index.mat, na.rm = TRUE), function(x) na.omit(c(ER_results$solution))[na.omit(c(ER_results$index.mat) == x)][1])
  
#ER ancestral reconstruction
ER_recon <- ancRECON(phy = ER_results$phy, data = ER_results$data, p = p, method = "marginal", rate.cat = ER_results$rate.cat, rate.mat = ER_results$index.mat, model = "ER", root.p = "yang", get.likelihood = TRUE, get.tip.states = TRUE, collapse = TRUE)
ER_recon <- ancRECON(phy = ER_results$phy, data = ER_results$data, p = p, method = "marginal", rate.cat = ER_results$rate.cat, model = "ER", root.p = "yang")

#ARD ancestral reconstruction
p <- sapply(1:max(ARD_results$index.mat, na.rm = TRUE), function(x) na.omit(c(ARD_results$solution))[na.omit(c(ARD_results$index.mat) == x)][1])
ARD_recon <- ancRECON(phy = ARD_results$phy, data = ARD_results$data, p = p, method = "marginal", rate.cat = ARD_results$rate.cat, model = "ARD", root.p = ARD_results$root.p)

#new metho (using  anc recon)

#lik.tip.states gives the trait states at  the tips, lik.anc.states gives the trait states at the nodes
#rbind into one dataframe
lik.anc <- as.data.frame(rbind(ARD_results$lik.tip.states, ARD_results$lik.anc.states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("cathemeral", "diurnal", "nocturnal")
#add row labels based on the tip and node positions
phylo_tree <- ARD_results$phy
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

#this is exactly the same as using the tip.states and states directly from the corHMM results 
ancestral_plot_di <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
ancestral_plot_di
ancestral_plot_noc <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
ancestral_plot_noc
ancestral_plot_cath <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
ancestral_plot_cath



