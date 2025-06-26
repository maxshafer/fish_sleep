# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(ggplot2)
library(RColorBrewer)

setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Load in Cox mammal tree + Amelia primary source data -----------------------

#need to use the mamm tree for modelling and ancestral trait reconstruction
#trees to make: trinary dataset -maxdinoc, trinary dataset -maxcrep, 5state

# ## Read in the mammalian phylogeny
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia primary source cetacean data
trait.data <- read.csv(here("cetaceans_full.csv"))
clade_name <- "cetcaeans"

#use below for primary source whippomorpha data
#trait.data <- read.csv(here("whippomorpha_full.csv"))
#clade_name <- "whippomorpha"

#use below for all artiodactyla
# trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
# clade_name <- "artiodactyla"

#use below for artiodactyla minus cetaceans
#trait.data <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv"))
#clade_name <- "artiodactyla_minus_cetaceans"

#trait.data <- read.csv(here("ruminants_full.csv"))
#clade_name <- "ruminants"

# trait.data <- read.csv(here("sleepy_mammals.csv"))
# clade_name <- "mammals"

#remove species with unknown diel patterns, optional
trait.data <- trait.data[!(trait.data$Diel_Pattern_2 %in% c("unknown")), ]

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

#does not include species with no activity pattern data (they are filtered out when processing the sleepy cetacean dataframe)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/all_cetacean_tree.png", width=40,height=40, units="cm",res=100)
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 1)
dev.off()

# Section 2: Plot all state trait data tree --------------------------------------

#creating a tree with all trait data categories (diurnal, nocturnal, di/crep, noc/crep, cathemeral)

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
#custom.colours <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

#plot the tree
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 1)
diel.plot.all 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", clade_name, "diel_plot_all_state.png"), width=24,height=18,units="cm",res=1200)
print(diel.plot.all)
dev.off()

trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
#trait.data.all$Diel_Pattern_2 <- str_replace_all(trait.data.all$Diel_Pattern_2, pattern = "cathemeral/crepuscular", replacement = "crepusuclar")

custom.colours <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#66C2A5", "red")
#custom.colours <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

#plot the tree
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all +
  geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent", width = 3) + 
  scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot.all <- diel.plot.all + geom_tiplab(size = 2, offset = 1.5) + theme(legend.position="none")
diel.plot.all

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "diel_plot_4_state_max_crep_no_label.pdf"))
diel.plot.all
dev.off()

#plot tree with confidence data
#trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2", "Confidence")]
trait.data.all <- trait.data[, c("tips", "Diel_Pattern_2", "Confidence")]
#change any sub 1 confidence to confidence 1
trait.data.all$Confidence[trait.data.all$Confidence == 0] <- 1
trait.data.all$Confidence[trait.data.all$Confidence == 0.5] <- 1
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

custom.colours.all <- c("grey90", "grey", "grey50", "grey30", "black","white", "#dd8ae7","pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "blue", "brown")
#custom.colours.all <- c("red", "red", "brown3", "grey50", "grey30", "black","#dd8ae7","pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")

#plot the tree, separate geoms for diel pattern and confidence
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2", "Confidence")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "black", width = 2.5)
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x + 2.5, y=y, fill = as.factor(Confidence)), inherit.aes = FALSE, colour = "black", width = 2.5) + scale_fill_manual(values = custom.colours.all)
diel.plot.all <- diel.plot.all #+ geom_tiplab(size = 2, offset = 3)
diel.plot.all 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", "final", clade_name, "_confidence_diel_plot_all_state.png"), width=20,height=14,units="cm",res=1200)
print(diel.plot.all)
dev.off() 

# Section 3: Plot 3-state maxdinoc tree --------------------------------------

#creating a tree with three trait states (cathemeral, diurnal, nocturnal)

#remove species with no data
trait.data.maxDN <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trait.data.maxDN$Diel_Pattern_2 <- str_replace_all(trait.data.maxDN$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data.maxDN$Diel_Pattern_2 <- str_replace_all(trait.data.maxDN$Diel_Pattern_2, "diurnal/crepuscular", "nocturnal")

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxDN$tips)

#plot the tree
custom.colours.MDN <- c("#dd8ae7", "#FC8D62", "#66C2A5")
diel.plot.MDN <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxDN[,c("tips", "Diel_Pattern_2")]
diel.plot.MDN <- diel.plot.MDN + geom_tile(data = diel.plot.MDN$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MDN)
diel.plot.MDN <- diel.plot.MDN + geom_tiplab(size = 2)
diel.plot.MDN 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", clade_name, "diel_plot_max_dinoc.png"), width=24,height=18,units="cm",res=1200)
print(diel.plot.MDN)
dev.off()

# Section 4: Plot 3-state maxcrep tree --------------------------------------

#creating a tree with three trait states (crepuscular, diurnal, nocturnal)

#remove species with no data
trait.data.maxcrep <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "nocturnal/crepuscular", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "diurnal/crepuscular", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "cathemeral", "crep/cath")
trait.data.maxcrep$Diel_Pattern_2 <- str_replace_all(trait.data.maxcrep$Diel_Pattern_2, "crep/cath", "crepuscular/cathemeral")

trpy_n <- keep.tip(mam.tree, tip = trait.data.maxcrep$tips)

#plot the tree
custom.colours.MC <- c("pink", "#FC8D62", "#66C2A5")
diel.plot.MC <- ggtree(trpy_n, layout = "circular") %<+% trait.data.maxcrep[,c("tips", "Diel_Pattern_2")]
diel.plot.MC <- diel.plot.MC + geom_tile(data = diel.plot.MC$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.MC)
diel.plot.MC <- diel.plot.MC + geom_tiplab(size = 2)
diel.plot.MC 

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", clade_name, "diel_plot_max_crep.png"), width=24,height=18,units="cm",res=1200)
print(diel.plot.MC)
dev.off()


# Section 5: Plot the ambiguity in the cetacean tree ----------------------

#look at the distribution of cetacean trees using ggdensitree https://yulab-smu.top/treedata-book/chapter4.html

#change this so we save out the result and load it in
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
cetacean_trees <- lapply(mammal_trees, function(x) subsetTrees(tree = x, subset_names = trait.data.all$tips))

#sample 100 random trees from the total 1k
#onehundred_trees <- cetacean_trees[sample(1:length(cetacean_trees), 100)]

#plot these using ggdensitree
# density_tree100 <- ggdensitree(onehundred_trees, alpha = 0.3, colour = "yellowgreen")
# density_tree100
# 
# #save out
# png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/density_tree100.png", width=24,height=18,units="cm",res=1200)
# print(density_tree100)
# dev.off()

#repeat the same process but with 50 trees
twenty_trees <- cetacean_trees[sample(1:length(cetacean_trees), 1000)]
ggdensitree(twenty_trees, alpha = 0.3, colour = "yellowgreen") + geom_tiplab(size=3, color = "black") 


png("C:/Users/ameli/OneDrive/Documents/R_projects/Cox_diel_plots/density_treebig.png", width=100,height=18,units="cm",res=1000)
density_20
dev.off()



# Section 6: Artiodactyla suborders ---------------------------------------

# ## Read in the mammalian phylogeny
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia data
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#filter for the suborder you want
#trait.data <- trait.data %>% filter(Order == "Ruminantia")
#trait.data <- trait.data %>% filter(Order == "Suina")
#trait.data <- trait.data %>% filter(Order == "Tylopoda")
trait.data <- trait.data %>% filter(Order == "Whippomorpha")

trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

#plot the tree
#custom.colours.all <- c("#dd8ae7", "pink", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "grey")
custom.colours.all <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "grey")
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x+1.1, y=y, fill = Diel_Pattern_2), width = 2, inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.all, name = "Temporal activity pattern")
diel.plot.all <- diel.plot.all + geom_tiplab(size = 4, offset = 2.2)
diel.plot.all 

#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Ruminantia_6_state.png", width=46,height=42,units="cm",res=1200)
#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Suina_6_state.png", width=46,height=40,units="cm",res=1200)
#png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Tylopoda_6_state.png", width=46,height=40,units="cm",res=200)
png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/Whippomorpha_6_state.png", width=46,height=40,units="cm", res= 1200)

print(diel.plot.all)
dev.off()


# Section 7: Venn diagram of confidence levels ---------------------------------------

library(venneuler)
MyVenn <- venneuler(c(A=86,B=85,C=189,D=92, E = 21, "A&B"=13, 
                      "B&C"=42,"C&D"=53, "D&E"=5,"C&D&E"=5,"C&E"=12))
MyVenn$labels <- c("A\n86","B\n85","C\n189","D\n92", "E\n21")
plot(MyVenn)

png("C:/Users/ameli/OneDrive/Documents/R_projects/fish_sleep/confidence_venn.png", width=20,height=20,units="cm", res= 1200)
MyVenn <- venneuler(c(A=86,B=85,C=189,D=92, E = 21, "A&B"=13, "A&C" = 9, "A&D" = 5, "A&E" = 5, 
                      "B&C"=42,"C&D"=53, "D&E"=7,"C&D&E"=5,"C&E"=12, "A&D&E" = 2, "B&C&E" = 7))
MyVenn$labels <- c("A\n86","B\n85","C\n189","D\n92", "E\n21")
plot(MyVenn)
dev.off()


# Section 8: Crepuscular vs day-night diel plots -------------------------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#use trait data below for Amelia data
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#we want to play cath, di, noc preference separately than crep- non crep preference
#so separate the diel pattern 2 column, into its component parts
trait.data <- separate(trait.data, col = Diel_Pattern_2, into = c("Diel", "Crepuscularity"), sep = "/")

custom_colours <- c("#dd8ae7","steelblue1", "#FC8D62", "#66C2A5")

crep_diel_tree <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel", "Crepuscularity")]
crep_diel_tree <- crep_diel_tree + geom_tile(data = crep_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x=x, fill = Diel, width = 4),  inherit.aes = FALSE, color = "transparent")
crep_diel_tree <- crep_diel_tree + geom_tile(data = crep_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x = x + 5, fill = Crepuscularity, width = 4), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom_colours, na.value = "grey80")
crep_diel_tree

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/diel_v_crepuscular.png", width=46,height=40,units="cm", res= 1200)
crep_diel_tree
dev.off()


# # Section 9: Mammal diel plots ------------------------------------------

#chose either Bennie or Maor dataset
trait.data$Diel_Pattern_2 <- trait.data$Bennie_diel
trait.data$Diel_Pattern_2 <- trait.data$Maor_diel

#remove species with no data
trait.data.all <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
trpy_n_all <- keep.tip(mam.tree, tip = trait.data.all$tips)

custom.colours <- c("#dd8ae7", "pink", "#FC8D62", "#66C2A5")

#plot the tree
diel.plot.all <- ggtree(trpy_n_all, layout = "circular") %<+% trait.data.all[,c("tips", "Diel_Pattern_2")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_all$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent", width = 4) + scale_fill_manual(values = custom.colours)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 1, offset = 4)
diel.plot.all 


# Section 10: Phylogenetic signal ----------------------------

#tutorial from http://blog.phytools.org/2012/11/testing-for-pagels-10.html

#First, let's simulate a tree & data for the demo, using the λ=0.854. We will use geiger::lambdaTree 
# simulate tree
tree <-pbtree(n=50)
# simulate data with lambda=0.854
x<-fastBM(lambdaTree(tree,0.854))

# fit lambda model
lambda <- phylosig(tree, x, method = "lambda")
lambda

#function phylosig tests for phylogenetic signal using two methods, K or lambda

#phylosig(tree = phylogenetic_tree, x = vector containing metrics for continuous trait, test = "K" or "lambda")
#tests to see if there is a strong phylogenetic signal
#null hypothesis that the trait is evolving purely through brownian motion under pagels lambda
#alternative hypothesis is that more closely related species tend to have more similar traits
#pagel's lambda is the scaling factor you have to adjust the phylogenetic tree by to get the evolution under brownian motion
#lambda can be 0 (no correlation btw species) or up to 1.00 (correlation between species equal to Brownian expectation)

#bloomberg's k is slightly different in that it is the scaled ratio of the variance among species 
#as determined by the variance between all the contrasts

#can we do this for discrete traits as well? 
# some people use delta statistic https://github.com/mrborges23/delta_statistic 
source("scripts/Amelia_delta_code.R")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("cetaceans_full.csv"))

#keep only necessary columns
trait.data <- trait.data[, c("Species_name", "Diel_Pattern", "tips")]
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data <- trait.data[trait.data$Diel_Pattern %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
#this doesn't replace anything in the cetacean tree so comment it out for now
mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#vector of trait data, with species in same order as in tree (tree$tip.label)
sps_order <- as.data.frame(mam.tree$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
trait.data <- merge(trait.data, sps_order, by = "tips")
trait.data <- trait.data[order(trait.data$id), ]
trait <- trait.data$Diel_Pattern

#now we calculate delta using their custom function
delta_diel <- delta(trait, mam.tree, 0.1, 0.0589, 1000, 10, 100)
#returns a value of 0.7779191, will change slightly every time its calculated
#significance is only determined relative to a random simulation of the trait data

#calculate the p-values
#works by randomly shuffling the order of the trait data and finding the delta of this new randomly ordered trait data
#repeats this 100 times to find a null distribution
#in this null distribution the traits don't associate by phylogeny any more than chance
#we can then see if there is a sufficiently small probability of getting our actual delta given the null distribution

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,mam.tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>delta_diel)/length(random_delta)
boxplot(random_delta)
abline(h=delta_diel,col="red")

#if p-value less than 0.05 there is evidence of phylogenetic signal between the trait and character
#if its more than 0.05 there is not evidence for phylogenetic signal
#p value is 0.04, so there is a phylogenetic signal for diel activity patterns in cetaceans


# Comparing delta statistics  ---------------------------------------------

#use a custom function that given a dataframe with species names and diel activity data returns the delta statistic for that group
#requires diel activity data to be in a column called Diel_Pattern
#requires species names to be in a column called tips with the format Genus_species

mammals_df <- read.csv(here("sleepy_mammals.csv"))
#start with species with activity patterns that agree based on both Maor and Bennie databases
mammals_df <- filter(mammals_df, match == "Yes") #leaves 1794 species
mammals_df <- mammals_df[, c("tips", "Bennie_diel", "Order", "Family")]
#mammals_df <- mammals_df[, c("tips", "Maor_diel", "Order", "Family")]
colnames(mammals_df) <- c("tips", "Diel_Pattern", "Order", "Family")
mammals_df <- relocate(mammals_df, "tips", .after = "Diel_Pattern")

#test run
source("scripts/Amelia_functions.R")

test <- calculateDelta(trait.data = mammals_df, taxonomic_group_name = "Afrosoricida")

#remove orders with less than 20 species
mammals_df_filtered <- mammals_df %>% group_by(Order) %>% filter(n()>20)
mammals_df_filtered <- mammals_df_filtered %>% group_by(Order) %>% filter(n()<50)


all_delta <- lapply(unique(mammals_df_filtered$Order), function(X) calculateDelta(trait.data = mammals_df_filtered, taxonomic_group_name = X))
delta_df <- do.call(rbind, all_delta)

ggplot(delta_df, aes(x= delta_diel, y= taxonomic_group_name)) + geom_bar(stat = "identity")

#add in my primary source data 
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[, c("Diel_Pattern_2", "tips", "Family", "Order")]
artio_full$Diel_Pattern <- artio_full$Diel_Pattern_2
artio_full$Diel_Pattern <- str_replace_all(artio_full$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
artio_full$Diel_Pattern <- str_replace_all(artio_full$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
artio_full$Diel_Pattern <- str_replace_all(artio_full$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
artio_full <- artio_full[artio_full$Diel_Pattern %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), c("tips", "Family", "Order", "Diel_Pattern")]
colnames(artio_full) <- c("tips", "Family", "Suborder", "Diel_Pattern")
artio_full$Order <- "Artiodactyla"
ggplot(artio_full, aes(x = Suborder, fill = Diel_Pattern)) + geom_bar(position = "fill")


artio_full$Order <- "Amelia_Artiodactyla"
artio_full <- relocate(artio_full, "Family", .after = "Order")
artio_full <- relocate(artio_full, "Diel_Pattern", .before = "tips")
artio_full <- artio_full[,c("Diel_Pattern", "tips", "Order", "Family")]

new_mammals <- rbind(artio_full, mammals_df_filtered)

#plot the proportion of species in each diel category for each order
#filter to orders with at least ??? 5 species?
mammals_df_filtered <- mammals_df %>% group_by(Order) %>% filter(n()>=5)

ggplot(mammals_df_filtered, aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill")

ggplot(new_mammals, aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill")

#for just my artiodactyla data
new_mammals %>% filter(Order == "Amelia_Artiodactyla") %>% ggplot(., aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill")

#for just my cetacean data
artio_full %>% filter(Suborder == "Whippomorpha") %>% ggplot(., aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill")

new_mammals %>% filter(Family == "Odontoceti") %>% ggplot(., aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill")
new_mammals %>% filter(Family == "Mysticeti") %>% ggplot(., aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill")


#for just my ruminant data

#for all mammals (data from bennie and maor?)
new_mammals$mammal <- "mammal"
ggplot(new_mammals, aes(x = mammal, fill = Diel_Pattern)) + geom_bar(position = "fill")
