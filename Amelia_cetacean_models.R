# Import and run packages -------------------------------------------------
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
library(dplyr)
## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)
library(readxl)

# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")

# Section 0 Load in cetacean diel data and mam tree -------------------------------

#read in the cetacean behavioural data 
cetaceans_full <- read.csv(here("cetaceans_full.csv"))

## Read in the mammalian phylogeny

mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")


# Section 1: Binary, cetacean model: di, noc --------------------------------------

## Quick change to diurnal / nocturnal only trait.data
## I use standard nomenclature for trees ('trpy_n') and trait data ('trait.data'), so that the script below can be run with any other tree etc

# This selects only data that is diurnal or nocturnal
trait.data <- cetaceans_full[cetaceans_full$Diel_Pattern_1 %in% c("diurnal", "nocturnal"),]
# selects only data that is in the mammal tree
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
row.names(trait.data) <- trait.data$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

## You should double check which species exist or not in your tree (I think we did this already)
## For the below, you can try all different combinations (diurnal/nocturnal only, or include crepuscular or cathemeral (or not), and see what difference it makes)

## Try simple models first, ER and ARD models
## This uses the ace package, but you can run these models with any of the libraries
## They have slightly different inputs, for example, ace wants a vector of traits
trait.vector <- trait.data$Diel_Pattern_1
ace_model_er1 <- ace(trait.vector, trpy_n, model = "ER", type = "discrete")
ace_model_sym1 <- ace(trait.vector, trpy_n, model = "SYM", type = "discrete")
ace_model_ard1 <- ace(trait.vector, trpy_n, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package

cor_model_er1 <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym1 <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard1 <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 1b Add likelihoods to table -----------------------------------

#create a table comparing all models
likelihoods <- data.frame(model = c("ace_model_er1", "ace_model_sym1", "ace_model_ard1"), likelihood = c(ace_model_er1$loglik, ace_model_sym1$loglik, ace_model_ard1$loglik), description = c("diurnal, nocturnal, cetacea", "diurnal, nocturnal, cetacea", "diurnal, nocturnal, cetacea"))
View(likelihoods)

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("cor_model_er1", cor_model_er$loglik, "diurnal, nocturnal, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym1", cor_model_sym$loglik, "diurnal, nocturnal, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard1", cor_model_ard$loglik, "diurnal, nocturnal, cetacea"))

saveRDS(likelihoods, here("likelihoods.rds"))

# Section 1c Plot ancestral recon and transition rates ------------------

## Extract the data to plot
## Depending on which package you use, the ancestral data is in different places (or not determined at all)
## It typically only includes the internal node states, and not also the tip states
str(ace_model_er1, max.level = 1)
head(ace_model_er1$lik.anc)

str(cor_model_er1, max.level = 1)
head(cor_model_er1$states)

# In this case, I am putting together the internal nodes and tip states into one data frame
lik.anc <- as.data.frame(rbind(cor_model_ard1$tip.states, cor_model_ard1$states))
# dim of this should be equal to the tips and internal nodes 
trpy_n
dim(lik.anc)

colnames(lik.anc) <- c("diurnal", "nocturnal")

# Tips are also nodes, and all nodes have a number, and they are number sequentially beginning with tips, so the first internal node is the number of tips +1
# attach each species (and ancestral sps/node) with its node number
lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

#colour the internal nodes and tips by the diurnal column (ie whether or not they are diurnal),
ancestral_plot_dinoc <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_dinoc

#label the tips (and save as png)
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_dinoc_cet.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_dinoc + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cor_model_rates_2states.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_model_ard1)
dev.off()

# Section 2a Four-state, cetacean model: di, noc, crep, cath --------------------------------

#repeating the same for crepuscular, cathemeral, diurnal nocturnal
#create dataframe to use. Starts with 98 species
trait.data2 <- cetaceans_full
#remove sps with no diel data. Should now have 80 species 
trait.data2 <- cetaceans_full[!(is.na(trait.data2$Diel_Pattern_3)),]
# selects only data that is in the mammal tree (start with 80 sps, only 71 in tree)
trait.data2 <- trait.data2[trait.data2$tips %in% mam.tree$tip.label,]
row.names(trait.data2) <- trait.data2$tips

# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n2 <- keep.tip(mam.tree, tip = trait.data2$tips)

trait.vector2 <- trait.data2$Diel_Pattern_3
ace_model_er2 <- ace(trait.vector2, trpy_n2, model = "ER", type = "discrete")
ace_model_sym2 <- ace(trait.vector2, trpy_n2, model = "SYM", type = "discrete")
ace_model_ard2 <- ace(trait.vector2, trpy_n2, model = "ARD", type = "discrete")

cor_model_er2 <- corHMM(phy = trpy_n2, data = trait.data2[trpy_n2$tip.label, c("tips", "Diel_Pattern_3")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym2 <- corHMM(phy = trpy_n2, data = trait.data2[trpy_n2$tip.label, c("tips", "Diel_Pattern_3")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard2 <- corHMM(phy = trpy_n2, data = trait.data2[trpy_n2$tip.label, c("tips", "Diel_Pattern_3")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 2b Add to table ------------------

#load in the likelihoods dataframe
likelihoods <- readRDS("likelihoods.rds")

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("ace_model_er2", ace_model_er2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_sym2", ace_model_sym2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_ard2", ace_model_ard2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("cor_model_er2", cor_model_er2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym2", cor_model_sym2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard2", cor_model_ard2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))

saveRDS(likelihoods, here("likelihoods.rds"))

# Section 2c Plot ancestral recon and transition rates ------------------

#create a dataframe of each species, its node and whether or not it is either of the diel patterns
cet.lik.anc <- as.data.frame(rbind(cor_model_ard2$tip.states, cor_model_ard2$states))
colnames(cet.lik.anc) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
#associate each row with specific nodes
cet.lik.anc$node <- c(1:length(trpy_n2$tip.label), (length(trpy_n2$tip.label) + 1):(trpy_n2$Nnode + length(trpy_n2$tip.label)))
custom.colours3 <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
#cath crep di noc
cet_ancestral_plot_di <- ggtree(trpy_n2, layout = "circular") %<+% cet.lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
cet_ancestral_plot_noc <- ggtree(trpy_n2, layout = "circular") %<+% cet.lik.anc + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
cet_ancestral_plot_crep <- ggtree(trpy_n2, layout = "circular") %<+% cet.lik.anc + aes(color = crepuscular) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5) + scale_color_distiller(palette = "Greys", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5)
cet_ancestral_plot_cath <- ggtree(trpy_n2, layout = "circular") %<+% cet.lik.anc + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)

png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard2_1.png", width=17,height=14,units="cm",res=1200)
cet_ancestral_plot_di
dev.off()
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard2_2.png", width=17,height=14,units="cm",res=1200)
cet_ancestral_plot_noc
dev.off()
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard2_3.png", width=17,height=14,units="cm",res=1200)
cet_ancestral_plot_cath
dev.off()
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard2_4.png", width=17,height=14,units="cm",res=1200)
cet_ancestral_plot_crep
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cor_rates_4states.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_model_ard2)
dev.off()

# Section 3a Three-state cetacean models: max-dinoc --------------------------------

#test with nocturnal, diurnal, cathemeral, with crepuscular species categorized as either di or noc
trait.data3 <- cetaceans_full
trait.data3 <- trait.data3[!(is.na(trait.data3$Diel_Pattern_2)),]
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")

trait.data3 <- trait.data3[trait.data3$tips %in% mam.tree$tip.label,]
row.names(trait.data3) <- trait.data3$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n3 <- keep.tip(mam.tree, tip = trait.data3$tips)

trait.vector3 <- trait.data3$Diel_Pattern_2
ace_model_er3 <- ace(trait.vector3, trpy_n3, model = "ER", type = "discrete")
ace_model_sym3 <- ace(trait.vector3, trpy_n3, model = "SYM", type = "discrete")
ace_model_ard3 <- ace(trait.vector3, trpy_n3, model = "ARD", type = "discrete")

cor_model_er3 <- corHMM(phy = trpy_n3, data = trait.data3[trpy_n3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym3 <- corHMM(phy = trpy_n3, data = trait.data3[trpy_n3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard3 <- corHMM(phy = trpy_n3, data = trait.data3[trpy_n3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 3b Add to table -----------------------------------------------

#load in the likelihoods dataframe
likelihoods <- readRDS("likelihoods.rds")

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("ace_model_er3", ace_model_er3$loglik, "maxdinoc: diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_sym3", ace_model_sym3$loglik, "maxdinoc: diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_ard3", ace_model_ard3$loglik, "maxdinoc: diurnal, nocturnal, cathemeral, cetacea"))

likelihoods <- rbind(likelihoods, c("cor_model_er3", cor_model_er3$loglik, "maxdinoc: diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym3", cor_model_sym3$loglik, "maxdinoc: diurnal, nocturnal, cathemeral, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard3", cor_model_ard3$loglik, "maxdinoc: diurnal, nocturnal, cathemeral, cetacea"))

saveRDS(likelihoods, here("likelihoods.rds"))

# Section 3c Plot ancestral recon and transition rates ------------------

#plot ancestral three states
cet.lik.anc3 <- as.data.frame(rbind(cor_model_ard3$tip.states, cor_model_ard3$states))
colnames(cet.lik.anc3) <- c("cathemeral", "diurnal", "nocturnal")
#associate each row with specific nodes
cet.lik.anc3$node <- c(1:length(trpy_n3$tip.label), (length(trpy_n3$tip.label) + 1):(trpy_n3$Nnode + length(trpy_n3$tip.label)))

cet_ancestral_plot_di3 <- ggtree(trpy_n3, layout = "circular") %<+% cet.lik.anc3 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
cet_ancestral_plot_noc3 <- ggtree(trpy_n3, layout = "circular") %<+% cet.lik.anc3 + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
cet_ancestral_plot_cath3 <- ggtree(trpy_n3, layout = "circular") %<+% cet.lik.anc3 + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)

#save as pngs
png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_maxdinoc_p1.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_di3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_maxdinoc_p2.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_noc3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_maxdinoc_p3.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_cath3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cor_model_rates_maxdinoc.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_model_ard3)
dev.off()

# Section 4a Three-state cetacean models: max-crep --------------------------------

#test with nocturnal, diurnal, cathemeral, with partially crepuscular species categorized as fully crepuscular
trait.data4 <- cetaceans_full
trait.data4 <- trait.data4[!(is.na(trait.data4$Diel_Pattern_2)),]
trait.data4$Diel_Pattern_2 <- str_replace_all(trait.data4$Diel_Pattern_2, "nocturnal/crepuscular", "crep/cath")
trait.data4$Diel_Pattern_2 <- str_replace_all(trait.data4$Diel_Pattern_2, "diurnal/crepuscular", "crep/cath")
trait.data4$Diel_Pattern_2 <- str_replace_all(trait.data4$Diel_Pattern_2, "cathemeral", "crep/cath")
trait.data4$Diel_Pattern_2 <- str_replace_all(trait.data4$Diel_Pattern_2, "crep/cath", "crepuscular/cathemeral")

trait.data4 <- trait.data4[trait.data4$tips %in% mam.tree$tip.label,]
row.names(trait.data4) <- trait.data4$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n4 <- keep.tip(mam.tree, tip = trait.data4$tips)

trait.vector4 <- trait.data4$Diel_Pattern_2
ace_model_er4 <- ace(trait.vector4, trpy_n4, model = "ER", type = "discrete")
ace_model_sym4 <- ace(trait.vector4, trpy_n4, model = "SYM", type = "discrete")
ace_model_ard4 <- ace(trait.vector4, trpy_n4, model = "ARD", type = "discrete")

cor_model_er4 <- corHMM(phy = trpy_n4, data = trait.data4[trpy_n4$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym4 <- corHMM(phy = trpy_n4, data = trait.data4[trpy_n4$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard4 <- corHMM(phy = trpy_n4, data = trait.data4[trpy_n4$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 4b Add to table -----------------------------------------------

#load in the likelihoods dataframe
likelihoods <- readRDS("likelihoods.rds")

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("ace_model_er4", ace_model_er4$loglik, "maxcrep: diurnal, nocturnal, cathemeral/crepuscular, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_sym4", ace_model_sym4$loglik, "maxcrep: diurnal, nocturnal, cathemeral/crepuscular, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_ard4", ace_model_ard4$loglik, "maxcrep: diurnal, nocturnal, cathemeral/crepuscular, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_er4", cor_model_er4$loglik, "maxcrep: diurnal, nocturnal, cathemeral/crepuscular, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym4", cor_model_sym4$loglik, "maxcrep: diurnal, nocturnal, cathemeral/crepuscular, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard4", cor_model_ard4$loglik, "maxcrep: diurnal, nocturnal, cathemeral/crepuscular, cetacea"))

saveRDS(likelihoods, here("likelihoods.rds"))

# Section 4c Plot ancestral recon and transition rates ------------------

#plot ancestral three states
cet.lik.anc4 <- as.data.frame(rbind(cor_model_ard4$tip.states, cor_model_ard4$states))
colnames(cet.lik.anc4) <- c("cathemeral_crepuscular", "diurnal", "nocturnal")
#associate each row with specific nodes
cet.lik.anc4$node <- c(1:length(trpy_n4$tip.label), (length(trpy_n4$tip.label) + 1):(trpy_n4$Nnode + length(trpy_n4$tip.label)))

cet_ancestral_plot_di4 <- ggtree(trpy_n4, layout = "circular") %<+% cet.lik.anc4 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
cet_ancestral_plot_noc4 <- ggtree(trpy_n4, layout = "circular") %<+% cet.lik.anc4 + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
cet_ancestral_plot_cath4 <- ggtree(trpy_n4, layout = "circular") %<+% cet.lik.anc4 + aes(color = cathemeral_crepuscular) + geom_tippoint(aes(color = cathemeral_crepuscular), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral_crepuscular), shape = 16, size = 1.5)

#save as pngs
png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_maxcrep_p1.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_di4
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_maxcrep_p2.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_noc4
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_maxcrep_p3.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_cath4
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cor_model_rates_maxcrep.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_model_ard4)
dev.off()

# Section 5 Bridge-only maxdinoc model --------------------------------------------
#want to model an evolutionary pattern where transitions between noc-di must occur through cath/crep intermediate
#don't allow transitions from nocturnality to diurnality
generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat
#need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> nocturnal/3 (6)
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(4, 6))
custom_rate_matrix

cor_cath_bridge <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)


#load in the likelihoods dataframe
likelihoods <- readRDS("likelihoods.rds")
#add how well this model did to the likelihood table
likelihoods <- rbind(likelihoods, c("cor_cath_bridge_maxdinoc", cor_cath_bridge$loglik, "three states (di, noc, cath), ARD, does not allow noc <-> di, cetacea"))
saveRDS(likelihoods, here("likelihoods.rds"))

# Section 6 Bridge-only maxcrep model --------------------------------------------
#want to model an evolutionary pattern where transitions between noc-di must occur through cath/crep intermediate
#don't allow transitions from nocturnality to diurnality
generic_ratemat <- getStateMat4Dat(cor_model_ard4$data)
test_ratemat <- generic_ratemat$rate.mat
#this is a generic all rates different matrix that we can now edit
test_ratemat
#need to disallow transitions from nocturnal/3 -> diurnal/2 (4) and diurnal/2 -> nocturnal/3 (6)
custom_rate_matrix <- dropStateMatPars(test_ratemat, c(4, 6))
custom_rate_matrix

cor_cath_bridge_maxcrep <- corHMM(phy = cor_model_ard4$phy, data = cor_model_ard4$data, rate.cat = 1, rate.mat = custom_rate_matrix)

#load in the likelihoods dataframe
likelihoods <- readRDS("likelihoods.rds")
#add how well this model did to the likelihood table
likelihoods <- rbind(likelihoods, c("cor_cath_bridge_maxcrep", cor_cath_bridge_maxcrep$loglik, "three states (di, noc, cath/crep), ARD, does not allow noc <-> di, cetacea"))
saveRDS(likelihoods, here("likelihoods.rds"))


# #Section 7: 5-state (di, noc, di/crep, noc/crep, cath) cetacean (1r and 2r) --------
#We want to run some models on dataset of cetaceans with the most detailed categories of diel patterns
# Di-crep, Di-not-crep, Noc-crep, Noc-not-crep, Cath-crep, Cath-not-crep
# diel pattern 2 in cetaceans full already has all cetacean species categorized this way except for cath-crep (ie true crepuscular)
#there are no purely crepuscular cetacean species, so we will have five categories

# We only need the tips (species names) and diel pattern 2. We start with 80 species with data
trait.data.crep <- cetaceans_full
trait.data.crep <- trait.data.crep[!(is.na(trait.data.crep$Diel_Pattern_2)),]
# selects only data that is in the mammal tree. 72 species are in the tree
trait.data.crep <- trait.data.crep[trait.data.crep$tips %in% mam.tree$tip.label,]
row.names(trait.data.crep) <- trait.data.crep$tips
# this subsets the tree to only include the species which we have data on (72 species)
trpy_n_crep <- keep.tip(mam.tree, tip = trait.data.crep$tips)

#first run the ace models
trait.vector.crep <- trait.data.crep$Diel_Pattern_2
ace_model_er_crep <- ace(trait.vector.crep, trpy_n_crep, model = "ER", type = "discrete")
ace_model_sym_crep <- ace(trait.vector.crep, trpy_n_crep, model = "SYM", type = "discrete")
ace_model_ard_crep <- ace(trait.vector.crep, trpy_n_crep, model = "ARD", type = "discrete")

#run the cor models
cor_model_er_crep <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym_crep <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard_crep <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

#repeat but with a hidden rate model
cor_model_er_crep_2r <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 2, model = "ER", node.states = "marginal")
cor_model_sym_crep_2r <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 2, model = "SYM", node.states = "marginal")
cor_model_ard_crep_2r <- corHMM(phy = trpy_n_crep, data = trait.data.crep[trpy_n_crep$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 2, model = "ARD", node.states = "marginal")

png("C:/Users/ameli/OneDrive/Documents/R_projects/cor_ard_crep_hiddenRate.png", width=14,height=20,units="cm",res=1200)
plotMKmodel(cor_model_ard_crep_2r)
dev.off()

#load in the likelihoods dataframe
likelihoods <- readRDS("likelihoods.rds")

#add likelihoods to the table
likelihoods <- rbind(likelihoods, c("ace_model_er_crep", ace_model_er_crep$loglik, "di/crep, noc.crep, diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_sym_crep", ace_model_sym_crep$loglik, "di/crep, noc/crep, diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_ard_crep", ace_model_ard_crep$loglik, "di/crep, noc/crep, diurnal, nocturnal, cathemeral, cetacea"))

likelihoods <- rbind(likelihoods, c("cor_model_er_crep", cor_model_er_crep$loglik, "di/crep, noc/crep, diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym_crep", cor_model_sym_crep$loglik, "di/crep, noc/crep, diurnal, nocturnal, cathemeral, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard_crep", cor_model_ard_crep$loglik, "di/crep, noc/crep, diurnal, nocturnal, cathemeral, cetacea"))

#save out updated table
saveRDS(likelihoods, here("likelihoods.rds"))

# Simmaps -----------------------------------------------------------------

#we can use functions from phytools to create character histories and plot them on the phylogeny
simmap_test <- makeSimmap(tree=cor_model_ard3$phy, data=cor_model_ard3$data, model=cor_model_ard3$solution, rate.cat=1, nSim=1, nCores=1)
# use plotSimmap to visualize (teal = nocturnal, orange = diurnal, purple = cathemeral)
simmap.colors <- setNames(c("#dd8ae7", "#FC8D62", "#66C2A5"), c("Cathemeral","Diurnal", "Nocturnal"))
plotSimmap(simmap_test[[1]], colors = simmap.colors, fsize=.5)
png("C:/Users/ameli/OneDrive/Documents/R_projects/diel.simmap.png", width=14,height=13,units="cm",res=1200)
plotSimmap(simmap_test[[1]], colors = simmap.colors, fsize=.5)
dev.off()
#this is stochastic mapping, creates a new tree every time based on the possible tree
#can we create an average of many simmmaps? 


# Export likelihood table -------------------------------------------------

#first add better model names
likelihoods <- readRDS("likelihoods.rds")
likelihoods$cetacean_model <- c("binary ER", "binary SYM", "binary ARD", "binary ER", "binary SYM", "binary ARD", "4-state ER", "4-state SYM", "4-state ARD", "4-state ER", "4-state SYM", "4-state ARD", "3-state max_dinoc ER", "3-state max_dinoc SYM", "3-state max_dinoc ARD", "3-state max_dinoc ER", "3-state max_dinoc SYM", "3-state max_dinoc ARD", "3-state max_crep ER", "3-state max_crep SYM", "3-state max_crep ARD", "3-state max_crep ER", "3-state max_crep SYM", "3-state max_crep ARD", "bridge-only maxdinoc", "bridge-only maxcrep", "5-state ER", "5-state SYM", "5-state ARD", "5-state ER", "5-state SYM", "5-state ARD")
acecor <- c(rep("ace", times = 3), rep("cor", times = 3))
likelihoods$package <- c(rep(acecor, times = 4), "cor", "cor", acecor)

#export dataframe of likelihood as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/likelihood_table.png", height = 30*nrow(likelihoods), width = 200*ncol(likelihoods), res = 90)
grid.table(likelihoods)
dev.off()

#make a plot of the model likelihoods

likelihood_plot <- ggplot(likelihoods, aes(x = cetacean_model, y = as.numeric(likelihood), colour = package)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png("C:/Users/ameli/OneDrive/Documents/R_projects/likelihood_plot.png")
likelihood_plot
dev.off()

#plot only the three state models we will be comparing
#row 13-26 of the current likelihoods dataframe
trinary_likelihoods <- likelihoods[13:26,]
trinary_likelihoods <- trinary_likelihoods[order(trinary_likelihoods$likelihood),]

trinary_likelihoods_plot <- ggplot(trinary_likelihoods, aes(x=factor(cetacean_model, level=unique(trinary_likelihoods$cetacean_model)), y = likelihood, colour = package)) + geom_point(size = 3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("model")

png("C:/Users/ameli/OneDrive/Documents/R_projects/trinary_likelihoods_plot.png")
trinary_likelihoods_plot
dev.off()