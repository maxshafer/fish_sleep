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

# Section 1a Binary dataframe, Cox???, Artiodactyla model ------------------------------------

## Read in the mammalian phylogeny
mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")
# read in the Cox binary df
artiodactyla_full <- read.csv("artiodactyla_binary_df.csv")

#we need to drop the missing species to create the tree
artiodactyla_full <- artiodactyla_full[artiodactyla_full$tips %in% mam.tree$tip.label,]
trpy_artio <- keep.tip(mam.tree, tip = artiodactyla_full$tips)
#tree of all the artio species with data in the mammal tree (shown without trait data)
ggtree(trpy_artio, layout = "fan", size = 0.5) + geom_tiplab(size = 1.5)

#isolate just the diurnal and nocturnal species (198 species) -all species in this tree are listed as only di or noc
trait.data <- artiodactyla_full[!(is.na(artiodactyla_full$Diel_Pattern_1)),]
# selects only data that is in the mammal tree (190 species)
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
row.names(trait.data) <- trait.data$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#now that we have a tree with the di/noc exclusive data we can run markov models on this data
trait.vector <- trait.data$Diel_Pattern_1
ace_binary_er <- ace(trait.vector, trpy_n, model = "ER", type = "discrete")
ace_binary_sym <- ace(trait.vector, trpy_n, model = "SYM", type = "discrete")
ace_binary_ard <- ace(trait.vector, trpy_n, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_binary_er <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_binary_sym <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_binary_ard <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 1b Add to table ----------------------------------

#create a df comparing all models
artio_likelihoods <- data.frame(model = c("ace_binary_er", "ace_binary_sym", "ace_binary_ard"), likelihood = c(ace_binary_er$loglik, ace_binary_sym$loglik, ace_binary_ard$loglik), description = c("diurnal, nocturnal, artiodactyla, ??? dataset", "diurnal, nocturnal, artiodactyla, ??? dataset", "diurnal, nocturnal, artiodactyla, ??? dataset"))

#add cor results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_binary_er", cor_binary_er$loglik, "diurnal, nocturnal, artiodactyla, ??? dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_binary_sym", cor_binary_sym$loglik, "diurnal, nocturnal, artiodactyla, ??? dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_binary_ard", cor_binary_ard$loglik, "diurnal, nocturnal, artiodactyla, ??? dataset"))

write.csv(artio_likelihoods, here("artio_likelihoods.csv"))

# Section 1c Plot ancestral recon and transition rates ------------------
## plot the transition rates for di/noc artiodactyla
png("C:/Users/ameli/OneDrive/Documents/R_projects/artio_binary_model.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_binary_ard)
dev.off()

## Extract the data to plot
## Depending on which package you use, the ancestral data is in different places (or not determined at all)
## It typically only includes the internal node states, and not also the tip states
str(ace_binary_er, max.level = 1)
head(ace_binary_er$lik.anc)

str(cor_binary_er, max.level = 1)
head(cor_binary_er$states)

# In this case, I am putting together the internal nodes and tip states into one data frame
lik.anc <- as.data.frame(rbind(cor_binary_ard$tip.states, cor_binary_ard$states))
# dim of this should be equal to the tips and internal nodes (190 tips and 189 nodes for 379 total)
trpy_n
dim(lik.anc)

colnames(lik.anc) <- c("diurnal", "nocturnal")

# Tips are also nodes, and all nodes have a number, and they are number sequentially beginning with tips, so the first internal node is the number of tips +1

lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

ancestral_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard_binary.png", width=17,height=16,units="cm",res=1200)
ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
dev.off()

# Section 2a Binary, Maor, Artiodactyla with cetacea model ----------------------------------------

#load in the tree
mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")
#load in the artiodactyla data
maor_final <- read.csv("artiodactyla_full.csv")

#update df to drop the artiodactyla species missing from the mam tree
#start with 245 species, 223 in the tree
maor_final <- maor_final[maor_final$tips %in% mam.tree$tip.label,]
trpy_n_maor1 <- keep.tip(mam.tree, tip = maor_final$tips)
#tree of all the artio species with data in the mammal tree (shown without trait data)
ggtree(trpy_n_maor1, layout = "fan", size = 0.5) + geom_tiplab(size = 1.5)

#First we run the model for just the diurnal and nocturnal species
maor.trait.data1 <- maor_final[!(is.na(maor_final$Diel_Pattern_1)),]
# selects only data that is in the mammal tree
maor.trait.data1 <- maor.trait.data1[maor.trait.data1$tips %in% mam.tree$tip.label,]
row.names(maor.trait.data1) <- maor.trait.data1$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_maor1 <- keep.tip(mam.tree, tip = maor.trait.data1$tips)

#now that we have the trait data for these 155 species and the subtree we can model
trait.vector.maor1 <- maor.trait.data1$Diel_Pattern_1
ace_maor_er1 <- ace(trait.vector.maor1, trpy_n_maor1, model = "ER", type = "discrete")
ace_maor_sym1 <- ace(trait.vector.maor1, trpy_n_maor1, model = "SYM", type = "discrete")
ace_maor_ard1 <- ace(trait.vector.maor1, trpy_n_maor1, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maor_er1 <- corHMM(phy = trpy_n_maor1, data = maor.trait.data1[trpy_n_maor1$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maor_sym1 <- corHMM(phy = trpy_n_maor1, data = maor.trait.data1[trpy_n_maor1$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maor_ard1 <- corHMM(phy = trpy_n_maor1, data = maor.trait.data1[trpy_n_maor1$tip.label, c("tips", "Diel_Pattern_1")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 2b Add to table -------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_er1", ace_maor_er1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_sym1", ace_maor_sym1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_ard1", ace_maor_ard1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_er1", cor_maor_er1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_sym1", cor_maor_sym1$loglik, "diurnal, nocturnal, artoidactyla, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_ard1", cor_maor_ard1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))


# Section 2c Plot ancestral recon and transition rates ------------------

#plot the most likely result
## Extract the data to plot
## Depending on which package you use, the ancestral data is in different places (or not determined at all)
## It typically only includes the internal node states, and not also the tip states
str(ace_maor_ard1, max.level = 1)
head(ace_maor_ard1$lik.anc)

str(cor_maor_ard1, max.level = 1)
head(cor_maor_ard1$states)

# In this case, I am putting together the internal nodes and tip states into one data frame
#still not sure how to plot ancestral reconstructions with ace 

#lik.anc.maor <- as.data.frame(rbind(dinoc.trait.data[,"Diel_Pattern_1"], ace_dinoc_ard1$lik.anc))
# dim of this should be equal to the tips and internal nodes (190 tips and 189 nodes for 379 total)
#trpy_n_maor1
#dim(lik.anc.maor)

#colnames(lik.anc) <- c("diurnal", "nocturnal")

# Tips are also nodes, and all nodes have a number, and they are number sequentially beginning with tips, so the first internal node is the number of tips +1

#lik.anc$node <- c(1:length(trpy_n_dinoc$tip.label), (length(trpy_n_dinoc$tip.label) + 1):(trpy_n_dinoc$Nnode + length(trpy_n_dinoc$tip.label)))

#ancestral_plot_maor_ace <- ggtree(trpy_n_maor1, layout = "circular") %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")

#ancestral_plot_maor_ace + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)

#cor model plots
lik.anc.maor1 <- as.data.frame(rbind(cor_maor_ard1$tip.states, cor_maor_ard1$states))
colnames(lik.anc.maor1) <- c("diurnal", "nocturnal")
lik.anc.maor1$node <- c(1:length(trpy_n_maor1$tip.label), (length(trpy_n_maor1$tip.label) + 1):(trpy_n_maor1$Nnode + length(trpy_n_maor1$tip.label)))

ancestral_plot_maor1 <- ggtree(trpy_n_maor1, layout = "circular") %<+% lik.anc.maor1 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard_maor1.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_maor1 + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
dev.off()

##plot the transition rates for maor, ard, dinoc
png("C:/Users/ameli/OneDrive/Documents/R_projects/maor_artio_model_dinoc.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maor_ard1)
dev.off()

# Section 3a Four-state, Maor, Artiodactyla with cetacea --------------------------------------------

#do the same (ancestral state reconstruction) but including crepuscular and cathemeral states
maor.trait.data2 <- maor_final[!(is.na(maor_final$Diel_Pattern_3)),]
# selects only data that is in the mammal tree
maor.trait.data2 <- maor.trait.data2[maor.trait.data2$tips %in% mam.tree$tip.label,]
row.names(maor.trait.data2) <- maor.trait.data2$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_maor2 <- keep.tip(mam.tree, tip = maor.trait.data2$tips)

#now that we have the trait data for these 155 Species_name and the subtree we can model
trait.vector.maor2 <- maor.trait.data2$Diel_Pattern_3
ace_maor_er2 <- ace(trait.vector.maor2, trpy_n_maor2, model = "ER", type = "discrete")
ace_maor_sym2 <- ace(trait.vector.maor2, trpy_n_maor2, model = "SYM", type = "discrete")
ace_maor_ard2 <- ace(trait.vector.maor2, trpy_n_maor2, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maor_er2 <- corHMM(phy = trpy_n_maor2, data = maor.trait.data2[trpy_n_maor2$tip.label, c("tips", "Diel_Pattern_3")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maor_sym2 <- corHMM(phy = trpy_n_maor2, data = maor.trait.data2[trpy_n_maor2$tip.label, c("tips", "Diel_Pattern_3")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maor_ard2 <- corHMM(phy = trpy_n_maor2, data = maor.trait.data2[trpy_n_maor2$tip.label, c("tips", "Diel_Pattern_3")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 3b Add to table ---------------------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_er2", ace_maor_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_sym2", ace_maor_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_ard2", ace_maor_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla,maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_er2", cor_maor_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_sym2", cor_maor_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artoidactyla, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_ard2", cor_maor_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))


# Section 3c Plot ancestral recon and transition rates ------------------

#cor model plots
lik.anc.maor2 <- as.data.frame(rbind(cor_maor_ard2$tip.states, cor_maor_ard2$states))
colnames(lik.anc.maor2) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
lik.anc.maor2$node <- c(1:length(trpy_n_maor2$tip.label), (length(trpy_n_maor2$tip.label) + 1):(trpy_n_maor2$Nnode + length(trpy_n_maor2$tip.label)))

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_4state_p1.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_maor2 <- ggtree(trpy_n_maor2, layout = "circular") %<+% lik.anc.maor2 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_maor2 + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_4state_p2.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_maor2 <- ggtree(trpy_n_maor2, layout = "circular") %<+% lik.anc.maor2 + aes(color = nocturnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_maor2 + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_4state_p3.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_maor2 <- ggtree(trpy_n_maor2, layout = "circular") %<+% lik.anc.maor2 + aes(color = cathemeral) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_maor2 + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_4state_p4.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_maor2 <- ggtree(trpy_n_maor2, layout = "circular") %<+% lik.anc.maor2 + aes(color = crepuscular) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_maor2 + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5)
dev.off()


# Section 4a Three state, Maor, Artiodactyla with cetacea --------------------------------------

#test models of artiodactyla with three states (cathemeral, diurnal, nocturnal)
maor.trait.data3 <- maor_final[!(is.na(maor_final$Diel_Pattern_2)),]
maor.trait.data3$Diel_Pattern_2 <- str_replace_all(maor.trait.data3$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
maor.trait.data3$Diel_Pattern_2 <- str_replace_all(maor.trait.data3$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")
#unlike cetaceans some artio are purely crepuscular, so switch into cathemeral
maor.trait.data3$Diel_Pattern_2 <- str_replace_all(maor.trait.data3$Diel_Pattern_2, "crepuscular", "cathemeral")

maor.trait.data3 <- maor.trait.data3[maor.trait.data3$tips %in% mam.tree$tip.label,]
row.names(maor.trait.data3) <- maor.trait.data3$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_maor3 <- keep.tip(mam.tree, tip = maor.trait.data3$tips)

maor.trait.vector3 <- maor.trait.data3$Diel_Pattern_2
ace_maor_er3 <- ace(maor.trait.vector3, trpy_n_maor3, model = "ER", type = "discrete")
ace_maor_sym3 <- ace(maor.trait.vector3, trpy_n_maor3, model = "SYM", type = "discrete")
ace_maor_ard3 <- ace(maor.trait.vector3, trpy_n_maor3, model = "ARD", type = "discrete")

#repeat with corHMM package
cor_maor_er3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maor_sym3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maor_ard3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 3.1: Setting the root state of the artiodactyla tree -------------

root1 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal", root.p = c(1,0,0))
root2 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal", root.p = c(0,1,0))
root3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("tips", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal", root.p = c(0,0,1))

rooted_comparisons <- data.frame(model = c("any_root", "root_1", "root_2", "root_3"), likelihood = c(cor_maor_ard3$loglik, root1$loglik, root2$loglik, root3$loglik))

# Section 4b Add to table -----------------------------------------

#add how well these models did to the table
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_er3", ace_maor_er3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_sym3", ace_maor_sym3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_ard3", ace_maor_ard3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))

artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_er3", cor_maor_er3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_sym3", cor_maor_sym3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_ard3", cor_maor_ard3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))


# Section 4c Plot ancestral recon and transition rates ------------------

#anc recon based on this model
#plot ancestral three states
art.lik.anc3 <- as.data.frame(rbind(cor_maor_ard3$tip.states, cor_maor_ard3$states))
colnames(art.lik.anc3) <- c("cathemeral", "diurnal", "nocturnal")
#associate each row with specific nodes
art.lik.anc3$node <- c(1:length(trpy_n_maor3$tip.label), (length(trpy_n_maor3$tip.label) + 1):(trpy_n_maor3$Nnode + length(trpy_n_maor3$tip.label)))

art_ancestral_plot_di3 <- ggtree(trpy_n_maor3, layout = "circular") %<+% art.lik.anc3 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
art_ancestral_plot_noc3 <- ggtree(trpy_n_maor3, layout = "circular") %<+% art.lik.anc3 + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
art_ancestral_plot_cath3 <- ggtree(trpy_n_maor3, layout = "circular") %<+% art.lik.anc3 + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
art_ancestral_plot_di3
art_ancestral_plot_noc3
art_ancestral_plot_cath3

#save as pngs
png("C:/Users/ameli/OneDrive/Documents/R_projects/best_artancrecon_3state_p1.png", width=16,height=15,units="cm",res=1200)
art_ancestral_plot_di3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_artancrecon_3state_p2.png", width=16,height=15,units="cm",res=1200)
art_ancestral_plot_noc3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_artancrecon_3state_p3.png", width=16,height=15,units="cm",res=1200)
art_ancestral_plot_cath3
dev.off()

#plot the rates for the artio 3 state models
png("C:/Users/ameli/OneDrive/Documents/R_projects/best_artiocor_model_rates_3states.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maor_ard3)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_artiocor_rates_4states.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maor_ard2)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_artiocor_model_rates_2states.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_dinoc_ard1)
dev.off()



# Section 5a Four-state, Maor, Artio without cetacea ------------------------
#load in the tree 

#load in the artiodactyla without cetacea df
mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")
#load in the artiodactyla data
just_artio <- read.csv("artiodactyla_without_cetaceans.csv")

#ancestral state reconstruction for 4 states, diurnal, nocturnal, crepuscular and cathemeral states
trait.data.ja2 <- just_artio[!(is.na(just_artio$Diel_Pattern_3)),]
# selects only data that is in the mammal tree 
trait.data.ja2 <- trait.data.ja2[trait.data.ja2$Species_name %in% mam.tree$tip.label,]
row.names(trait.data.ja2) <- trait.data.ja2$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_ja2 <- keep.tip(mam.tree, tip = trait.data.ja2$Species)

#now that we have the trait data for these 151 species and the subtree we can model
trait.vector.ja2 <- trait.data.ja2$Diel_Pattern_3
ace_maorja_er2 <- ace(trait.vector.ja2, trpy_n_ja2, model = "ER", type = "discrete")
ace_maorja_sym2 <- ace(trait.vector.ja2, trpy_n_ja2, model = "SYM", type = "discrete")
ace_maorja_ard2 <- ace(trait.vector.ja2, trpy_n_ja2, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maorja_er2 <- corHMM(phy = trpy_n_ja2, data = trait.data.ja2[trpy_n_ja2$tip.label, c("Species", "Diel_Pattern_3")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maorja_sym2 <- corHMM(phy = trpy_n_ja2, data = trait.data.ja2[trpy_n_ja2$tip.label, c("Species", "Diel_Pattern_3")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maorja_ard2 <- corHMM(phy = trpy_n_ja2, data = trait.data.ja2[trpy_n_ja2$tip.label, c("Species", "Diel_Pattern_3")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 5b Add to table -------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_er2", ace_maorja_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_sym2", ace_maorja_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_ard2", ace_maorja_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_er2", cor_maorja_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_sym2", cor_maorja_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artoidactyla without cetacea, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_ard2", cor_maorja_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))

# Section 5c Plot ancestral recon and transition rates -----------------

##plot the transition rates
png("C:/Users/ameli/OneDrive/Documents/R_projects/just_artio_ard_4state.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maorja_ard2)
dev.off()


#plot the best performing model
lik.anc.ja2 <- as.data.frame(rbind(cor_maorja_ard2$tip.states, cor_maorja_ard2$states))
colnames(lik.anc.ja2) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
lik.anc.ja2$node <- c(1:length(trpy_n_ja2$tip.label), (length(trpy_n_ja2$tip.label) + 1):(trpy_n_ja$Nnode + length(trpy_n_ja2$tip.label)))

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_ja4state_p1.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_ja_4state <- ggtree(trpy_n_ja2, layout = "circular") %<+% lik.anc.ja2 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_ja_4state + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_ja4state_p2.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_ja4state <- ggtree(trpy_n_ja2, layout = "circular") %<+% lik.anc.ja2 + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_ja4state + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_ja4state_p3.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_ja4state <- ggtree(trpy_n_ja2, layout = "circular") %<+% lik.anc.ja2 + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_ja4state + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/anc_recon_ja4state_p4.png", width=17,height=16,units="cm",res=1200)
ancestral_plot_ja4state <- ggtree(trpy_n_ja2, layout = "circular") %<+% lik.anc.ja2 + aes(color = crepuscular) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ancestral_plot_ja4state + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5)
dev.off()


# Section 6a Binary, Maor, Artio without cetacea ----------------
#repeat but with di/noc
trait.data.ja1 <- just_artio[!(is.na(just_artio$Diel_Pattern_1)),]
# selects only data that is in the mammal tree 
trait.data.ja1 <- trait.data.ja1[trait.data.ja1$Species %in% mam.tree$tip.label,]
row.names(trait.data.ja1) <- trait.data.ja1$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_ja1 <- keep.tip(mam.tree, tip = trait.data.ja1$Species)

trait.vector.ja1 <- trait.data.ja1$Diel_Pattern_1
ace_maorja_er1 <- ace(trait.vector.ja1, trpy_n_ja1, model = "ER", type = "discrete")
ace_maorja_sym1 <- ace(trait.vector.ja1, trpy_n_ja1, model = "SYM", type = "discrete")
ace_maorja_ard1 <- ace(trait.vector.ja1, trpy_n_ja1, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maorja_er1 <- corHMM(phy = trpy_n_ja1, data = trait.data.ja1[trpy_n_ja1$tip.label, c("Species", "Diel_Pattern_1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maorja_sym1 <- corHMM(phy = trpy_n_ja1, data = trait.data.ja1[trpy_n_ja1$tip.label, c("Species", "Diel_Pattern_1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maorja_ard1 <- corHMM(phy = trpy_n_ja1, data = trait.data.ja1[trpy_n_ja1$tip.label, c("Species", "Diel_Pattern_1")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 6b Add to table --------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_er1", ace_maorja_er1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_sym1", ace_maorja_sym1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_ard1", ace_maorja_ard1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_er1", cor_maorja_er1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_sym1", cor_maorja_sym1$loglik, "diurnal, nocturnal,  artoidactyla without cetacea, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_ard1", cor_maorja_ard1$loglik, "diurnal, nocturnal,  artiodactyla without cetacea, maor dataset"))


# Section 6c Plot ancestral recon and transition rates ------------------

##plot the transition rates
png("C:/Users/ameli/OneDrive/Documents/R_projects/just_artio_ard_dinoc.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maorja_ard1)
dev.off()

#plot the ancestral reconstruction
#need to do


# Section 7a Three-state, Maor, Artio without cetacea -------------------

#artio without cetacea three states models
trait.data.ja3 <- just_artio[!(is.na(just_artio$Diel_Pattern_2)),]
trait.data.ja3$Diel_Pattern_2 <- str_replace_all(trait.data.ja3$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data.ja3$Diel_Pattern_2 <- str_replace_all(trait.data.ja3$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")
#unlike cetaceans some artio are purely crepuscular, so switch into cathemeral
trait.data.ja3$Diel_Pattern_2 <- str_replace_all(trait.data.ja3$Diel_Pattern_2, "crepuscular", "cathemeral")

trait.data.ja3 <- trait.data.ja3[trait.data.ja3$Species %in% mam.tree$tip.label,]
row.names(trait.data.ja3) <- trait.data.ja3$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_ja3 <- keep.tip(mam.tree, tip = trait.data.ja3$Species)

trait.vector.ja3 <- trait.data.ja3$Diel_Pattern_2
ace_maorja_er3 <- ace(trait.vector.ja3, trpy_n_ja3, model = "ER", type = "discrete")
ace_maorja_sym3 <- ace(trait.vector.ja3, trpy_n_ja3, model = "SYM", type = "discrete")
ace_maorja_ard3 <- ace(trait.vector.ja3, trpy_n_ja3, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maorja_er3 <- corHMM(phy = trpy_n_ja3, data = trait.data.ja3[trpy_n_ja3$tip.label, c("Species", "Diel_Pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maorja_sym3 <- corHMM(phy = trpy_n_ja3, data = trait.data.ja3[trpy_n_ja3$tip.label, c("Species", "Diel_Pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maorja_ard3 <- corHMM(phy = trpy_n_ja3, data = trait.data.ja3[trpy_n_ja3$tip.label, c("Species", "Diel_Pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 7b Add to table and plot ------------------------------------------------------------
#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_er3", ace_maorja_er3$loglik, "diurnal, nocturnal, cathemeral. artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_sym3", ace_maorja_sym3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_ard3", ace_maorja_ard3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_er3", cor_maorja_er3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_sym3", cor_maorja_sym3$loglik, "diurnal, nocturnal, cathemeral, artoidactyla without cetacea, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_ard3", cor_maorja_ard3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))


# Section 7c Plot ancestral reconstruction and transition rates ---------

#plot the transition rates
png("C:/Users/ameli/OneDrive/Documents/R_projects/just_artio_ard_3state.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maorja_ard3)
dev.off()


# Export likelihood table -------------------------------------------------
#export dataframe of likelihood as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/artio_likelihood_table.png", height = 30*nrow(artio_likelihoods), width = 300*ncol(artio_likelihoods), res = 90)
grid.table(artio_likelihoods)
dev.off()

likelihood_plot <- ggplot(artio_likelihoods, aes(x = model, y = as.numeric(likelihood), colour = package)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
