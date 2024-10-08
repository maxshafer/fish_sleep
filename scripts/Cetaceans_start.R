# Import and run packages -------------------------------------------------
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


# Section 0 Load in and clean cetacean data -------------------------------

###### For cetaceans

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit?usp=sharing'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
table(cetaceans_full$Diel_Pattern_1)
View(cetaceans_full)
#they use a different name for sperm whales in the mam tree so we can change that
cetaceans_full$Species_name <- str_replace(cetaceans_full$Species_name, "Physeter catodon", "Physeter macrocephalus")
cetaceans_full <- cetaceans_full[,1:9]
#diel column is based on diel pattern 3 so that cathemeral species aren't dropped when we drop na values
cetaceans_full$diel <- tolower(cetaceans_full$Diel_Pattern_3)


## Probably should save out a local copy in case google goes bankrupt
write.csv(cetaceans, file = here("sleepy_fish_database_local.csv"))

## Remove species without diel data
cetaceans_full <- cetaceans_full[!(is.na(cetaceans_full$diel)),]


## Read in the mammalian phylogeny
#mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
#mam.tree <- maxCladeCred(mammal_trees, tree = TRUE) #this function takes a long time to run so save the result out
getwd()

# Save out maxcladecred, so we don't have to recalculate it every time
#saveRDS(mam.tree,"maxCladeCred_mammal_tree.rds")

mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")
#check to see what cetaceans are in mammal tree by subsetting to just include cetaceans
#use cetaceans_full because it has all species
cetaceans_full$tips <- str_replace(cetaceans_full$Species_name, " ", "_")
#trpy_n_test <- keep.tip(mam.tree, tip = cetaceans_full$tips) -uncomment to see missing species
#mam tree is missing 12 species
#described in 2014: Sousa_plumbea, Sousa_sahulensis, Mesoplodon_hotaula
#described in 2019: Berardius_minimus, Inia_araguaiaensis, 
#described in 2021: Balaenoptera_ricei, Mesoplodon_eueu, Platanista_minor
#subspecies, not a valid species: Inia_humboldtiana, Neophocaena_sunameri, Inia_boliviensis, Balaenoptera_brydei

#we need to drop the missing species to create the tree
all_cetaceans <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n_test <- keep.tip(mam.tree, tip = all_cetaceans$tips)
#check if time calibrated
test <- ggtree(trpy_n_test, layout = "circular", size = 0.5) + geom_tiplab(size = 1.5)
#test$data #29.36 million years ago to the root. Seems correct!

#look at the tree for whippomorpha (hippos + cetacea) and see if time calibrated
whippomorpha <- cetaceans_full
whippomorpha <- rbind(cetaceans_full, c("Hexaprotodon liberiensis", "", NA, NA, "Nocturnal", "Nocturnal", "Nocturnal", "", "3", "nocturnal", "Choeropsis_liberiensis"))
whippomorpha <- rbind(whippomorpha, c("Hippopotamus amphibius", "", NA, NA, "Nocturnal", "Nocturnal/crepuscular", "Crepuscular", "", "4", "crepuscular", "Hippopotamus_amphibius"))
whippomorpha <- whippomorpha[whippomorpha$tips %in% mam.tree$tip.label,]
trpy_n_whippo <- keep.tip(mam.tree, tip = whippomorpha$tips)
whippo <- ggtree(trpy_n_whippo, layout = "circular", size = 0.5) + geom_tiplab(size = 1.5)
#whippo$data shows the branch lengths (53.7 my to root). Also seems correct!
                     
write.csv(cetaceans_full, file = here("sleepy_fish_database_local.csv"))

# Section 1.1 Binary, cetacean model --------------------------------------

## Quick change to diurnal / nocturnal only trait.data
## I use standard nomenclature for trees ('trpy_n') and trait data ('trait.data'), so that the script below can be run with any other tree etc

# This selects only data that is diurnal or nocturnal
trait.data <- cetaceans_full[cetaceans_full$diel %in% c("diurnal", "nocturnal"),]
# selects only data that is in the mammal tree
trait.data$tips <- str_replace(trait.data$Species_name, " ", "_")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
row.names(trait.data) <- trait.data$tips
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

## You should double check which species exist or not in your tree (I think we did this already)
## For the below, you can try all different combinations (diurnal/nocturnal only, or include crepuscular or cathemeral (or not), and see what difference it makes)

## Try simple models first, ER and ARD models
## This uses the ace package, but you can run these models with any of the libraries
## They have slightly different inputs, for example, ace wants a vector of traits
trait.vector <- trait.data$diel
ace_model_er1 <- ace(trait.vector, trpy_n, model = "ER", type = "discrete")
ace_model_sym1 <- ace(trait.vector, trpy_n, model = "SYM", type = "discrete")
ace_model_ard1 <- ace(trait.vector, trpy_n, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package

cor_model_er1 <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_model_sym1 <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_model_ard1 <- corHMM(phy = trpy_n, data = trait.data[trpy_n$tip.label, c("tips", "diel")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 1.1b Add likelihoods to table -----------------------------------

#create a table comparing all models
likelihoods <- data.frame(model = c("ace_model_er1", "ace_model_sym1", "ace_model_ard1"), likelihood = c(ace_model_er1$loglik, ace_model_sym1$loglik, ace_model_ard1$loglik), description = c("diurnal, nocturnal, cetacea", "diurnal, nocturnal, cetacea", "diurnal, nocturnal, cetacea"))
View(likelihoods)

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("cor_model_er1", cor_model_er$loglik, "diurnal, nocturnal, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym1", cor_model_sym$loglik, "diurnal, nocturnal, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard1", cor_model_ard$loglik, "diurnal, nocturnal, cetacea"))


# Section 1.1c Plot ancestral recon and transition rates ------------------

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

# Section 1.2a Four-state, cetacean model --------------------------------

#repeating the same for crepuscular, cathemeral, diurnal nocturnal
View(cetaceans_full)
#create dataframe to use. Starts with 98 species
trait.data2 <- cetaceans_full[, 1:10]
#remove sps with no diel data. Should now have 80 species 
trait.data2 <- cetaceans_full[!(is.na(trait.data2$Diel_Pattern_3)),]
# formats the species names to be the same as in the tree
trait.data2$tips <- str_replace(trait.data2$Species_name, " ", "_")
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


# Section 1.2b Add to table ------------------

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("ace_model_er2", ace_model_er2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_sym2", ace_model_sym2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_ard2", ace_model_ard2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("cor_model_er2", cor_model_er2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym2", cor_model_sym2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard2", cor_model_ard2$loglik, "diurnal, nocturnal, crepuscular, cathemeral, cetacea"))

# Section 1.2c Plot ancestral recon and transition rates ------------------

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

# Section 1.3a Three-state cetacean models --------------------------------

#test with nocturnal, diurnal, cathemeral, with crepuscular species categorized as either di or noc
trait.data3 <- cetaceans_full
trait.data3 <- trait.data3[!(is.na(trait.data3$Diel_Pattern_2)),]
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data3$Diel_Pattern_2 <- str_replace_all(trait.data3$Diel_Pattern_2, "diurnal/crepuscular", "diurnal")

trait.data3$tips <- str_replace_all(trait.data3$Species_name, " ", "_")
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

# Section 1.3b Add to table -----------------------------------------------

#add how well these models did to the table
likelihoods <- rbind(likelihoods, c("ace_model_er3", ace_model_er3$loglik, "diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_sym3", ace_model_sym3$loglik, "diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("ace_model_ard3", ace_model_ard3$loglik, "diurnal, nocturnal, cathemeral, cetacea"))

likelihoods <- rbind(likelihoods, c("cor_model_er3", cor_model_er3$loglik, "diurnal, nocturnal, cathemeral, cetacea"))
likelihoods <- rbind(likelihoods, c("cor_model_sym3", cor_model_sym3$loglik, "diurnal, nocturnal, cathemeral, cetacea"))  
likelihoods <- rbind(likelihoods, c("cor_model_ard3", cor_model_ard3$loglik, "diurnal, nocturnal, cathemeral, cetacea"))


# Section 1.3c Plot ancestral recon and transition rates ------------------

#plot ancestral three states
cet.lik.anc3 <- as.data.frame(rbind(cor_model_ard3$tip.states, cor_model_ard3$states))
colnames(cet.lik.anc3) <- c("cathemeral", "diurnal", "nocturnal")
#associate each row with specific nodes
cet.lik.anc3$node <- c(1:length(trpy_n3$tip.label), (length(trpy_n3$tip.label) + 1):(trpy_n3$Nnode + length(trpy_n3$tip.label)))

cet_ancestral_plot_di3 <- ggtree(trpy_n3, layout = "circular") %<+% cet.lik.anc3 + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
cet_ancestral_plot_noc3 <- ggtree(trpy_n3, layout = "circular") %<+% cet.lik.anc3 + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_distiller(palette = "GnBu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
cet_ancestral_plot_cath3 <- ggtree(trpy_n3, layout = "circular") %<+% cet.lik.anc3 + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdPu", direction = 1) + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)

#save as pngs
png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_3state_p1.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_di3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_3state_p2.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_noc3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cetancrecon_3state_p3.png", width=16,height=15,units="cm",res=1200)
cet_ancestral_plot_cath3
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/best_cor_model_rates_3states.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_model_ard3)
dev.off()

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
#can we create a max credibility simmmap?


# Export likelihood table -------------------------------------------------


#export dataframe of likelihood as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/likelihood_table.png", height = 30*nrow(likelihoods), width = 200*ncol(likelihoods), res = 90)
grid.table(likelihoods)
dev.off()

#hello

