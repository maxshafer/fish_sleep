# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
#install.packages("deeptime")
library(deeptime)
#update.packages("ggplot2")
library(ggplot2)
library(gridExtra)
setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Select the diel dataset -----------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#uncomment whichever clade you want to plot
# clade_name <- "cetaceans_full"
# clade_name <- "sleepy_artiodactyla_full"
# clade_name <- "ruminants_full"
clade_name <- "whippomorpha"
# clade_name <- "sleepy_artiodactyla_minus_cetaceans"

diel_full <- read.csv(here(paste0(clade_name, ".csv")))

trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#use below to remove NA species
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern), ]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

# Section 2: Create all the desired plots ---------------------------------

#all diel patterns
custom.colours <- c("#dd8ae7", "peachpuff2" ,"#FC8D62", "gold", "#66C2A5", "#A6D854","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3.2) 
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_six_state_plot_labelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 

custom.colours <- c("#dd8ae7", "peachpuff2" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3.2) 
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_max_crep_plot_labelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 

custom.colours <- c("#dd8ae7", "peachpuff2" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_max_crep_plot_unlabelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 


# Section 3: Clade labels -------------------------------------------------
cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#add clade labels
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Mysticeti")
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Odontoceti")

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular", fill = "transparent") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + 
  geom_cladelab(node = 137, label = "Mysticeti", align = TRUE, geom = "label", offset=1, align=TRUE, offset.text=1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot <- diel.plot + 
  geom_cladelab(node = 77, label = "Odontoceti", align = FALSE, geom = "label", offset=1, align=FALSE, offset.text=1, hjust = 1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_unlabelled_cladelabels.pdf", bg = "transparent")
diel.plot
dev.off()

#label major families
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Mysticeti")
cetacean_families <- trait.data %>% count(Family) %>% filter(n>1) #filter for clades with more than one species or it can't find the MRCA
#lapply(list(cetacean_families$Family), findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = x))

#I don't know why the lapply wasn't working so do a for loop
for(i in cetacean_families$Family){
  node_df <- findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = i)
  nodes_list[[i]] <- node_df
  nodes_df <- do.call(rbind, nodes_list)
}

nodes_df$hjust <- c(0,0, 1, 0, 1, 1, 1)
nodes_df$colour <- c("red", "pink", "grey", "blue", "slateblue", "black", "lightblue")
  
custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_cladelab(node = nodes_df$node_number, label = nodes_df$clade_name, align = TRUE, geom = "label", offset=1, align=TRUE, offset.text=1, barsize=2, fontsize=3, fill = "grey", barcolour = nodes_df$colour, textcolour = "black")
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_unlabelled_famlabels.pdf", width = 10, height = 10, bg = "transparent")
diel.plot
dev.off()


# Section 4 breakdown of % activity patterns ----------------------------
mammals_df <- read.csv(here("sleepy_mammals.csv")) #data from Bennie et al, 2014

#add in my primary source data 
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ]

#replace their artiodactyla data with my artiodactyla data?
# mammals_df <- mammals_df %>% filter(Order != "Artiodactyla")
# new_mammals <- rbind(mammals_df, artio_full)

#we want to compare all mammals, vs all artiodactyla vs cetaceans/ruminants
new_mammals$mammals <- "Mammals"

custom.colours <- c("#dd8ae7","peachpuff2", "#FC8D62", "#66C2A5")
mammals_plot <- ggplot(new_mammals, aes(x = mammals, fill = max_crep)) + geom_bar(position = "fill", width = 0.75) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.x = element_blank())
artiodactyla_plot <- artio_full %>% filter(Order == "Artiodactyla") %>% ggplot(., aes(x = Order, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
ruminantia_plot <- artio_full %>% filter(Suborder == "Ruminantia") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
whippomorpha_plot <- artio_full %>% filter(Suborder == "Whippomorpha") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_classic() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/barplot_percentages.pdf", width = 10, height = 10, bg = "transparent")
grid.arrange(mammals_plot, artiodactyla_plot, ruminantia_plot, whippomorpha_plot, nrow = 1)
dev.off()

# Section 2: Phylogenetic signal -----------------

#can we do this for discrete traits as well? 
# some people use delta statistic https://github.com/mrborges23/delta_statistic 
source("scripts/Amelia_delta_code.R")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("cetaceans_full.csv"))

#keep only necessary columns
trait.data <- trait.data[, c("Species_name", "max_crep", "tips")]
trait.data <- trait.data[trait.data$max_crep %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]

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
trait <- trait.data$max_crep

#now we calculate delta using their custom function
delta_diel <- delta(trait, mam.tree, 0.1, 0.0589, 1000, 10, 100)
#returns a value of 0.7779191, will change slightly every time its calculated
#significance is only determined relative to a random simulation of the trait data

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

