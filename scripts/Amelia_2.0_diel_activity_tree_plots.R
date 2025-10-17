setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Select the diel dataset -----------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#uncomment whichever clade you want to plot
clade_name <- "cetaceans_full"
#clade_name <- "sleepy_artiodactyla_full"
#clade_name <- "ruminants_full"
#clade_name <- "whippomorpha"
#clade_name <- "whippomorpha_high_conf"
# clade_name <- "sleepy_artiodactyla_minus_cetaceans"

diel_full <- read.csv(here(paste0(clade_name, ".csv")))

trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#use below to remove NA species
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern), ]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

# Section 2: Create all the desired plots ---------------------------------

#all diel patterns
custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "gold", "#66C2A5", "#A6D854","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3.2) 
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_six_state_plot_labelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 

custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3.2) 
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_max_crep_plot_labelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 

custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_max_crep_plot_unlabelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 


# Section 3: Clade labels -------------------------------------------------
#cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- read.csv(here("whippomorpha.csv"))
cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#add clade labels
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Mysticeti")
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Odontoceti")

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
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



# Whippomorpha with family labels ----------------------------------------
diel_full <- read.csv(here("whippomorpha.csv"))
diel_full <- diel_full[!is.na(diel_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#label major families
families <- trait.data %>% count(Family) %>% filter(n>1) #filter for clades with more than one species or it can't find the MRCA

#Use a for loop to find the nodes for all the families
nodes_list <- list()
for(i in families$Family){
  node_df <- findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 5, taxonomic_level_name = i)
  nodes_list[[i]] <- node_df
  nodes_df <- do.call(rbind, nodes_list)
}

nodes_left <- nodes_df[c("Delphinidae", "Phocoenidae", "Monodontidae", "Ziphiidae"),]
nodes_right <- nodes_df[c("Kogiidae", "Balaenopteridae", "Balaenidae", "Hippopotamidae"),]

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_cladelab(node = nodes_right$node_number, label = nodes_right$clade_name, offset=1.5, offset.text=2, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + geom_cladelab(node = nodes_left$node_number, label = nodes_left$clade_name, offset=1.5, offset.text=2, hjust = 1, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_nodelab()
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/whippo_with_families.pdf", width = 8, height = 7, bg = "transparent")
diel.plot
dev.off()

#use geom tile as a clade label instead
custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
custom.colours.2 <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5", "red", "grey","blue", "pink", "chocolate", "black", "skyblue", "magenta", "orchid", "orange")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "Family")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Family), inherit.aes = FALSE) #+ scale_fill_manual(values = custom.colours.2)
diel.plot

# Section 4: Artiodactyla with suborder labels ----------------------------
diel_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
diel_full <- diel_full[!is.na(diel_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#label major families
suborders <- trait.data %>% count(Suborder) %>% filter(n>1) #filter for clades with more than one species or it can't find the MRCA

#For loop to find the labels for each suborder
nodes_list <- list()
for(i in suborders$Suborder){
  node_df <- findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = i)
  nodes_list[[i]] <- node_df
  nodes_df <- do.call(rbind, nodes_list)
}

nodes_left <- nodes_df[c("Ruminantia"),]
nodes_right <- nodes_df[c("Whippomorpha", "Tylopoda", "Suina"),]

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_cladelab(node = nodes_right$node_number, label = nodes_right$clade_name, offset=1.5, offset.text=2, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + geom_cladelab(node = nodes_left$node_number, label = nodes_left$clade_name, offset=1.5, offset.text=2, hjust = 1, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/artio_with_suborders.pdf", width = 7, height = 7, bg = "transparent")
diel.plot
dev.off()


# Section 5: breakdown of % activity patterns proportion plots----------------------------
new_mammals <- read.csv(here("sleepy_mammals.csv")) #data from Bennie et al, 2014

#add in my primary source data 
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ]


#we want to compare all mammals, vs all artiodactyla vs cetaceans/ruminants
new_mammals$mammals <- "Mammals"

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
mammals_plot <- ggplot(new_mammals, aes(x = mammals, fill = max_crep)) + geom_bar(position = "fill", width = 0.75) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.x = element_blank())
artiodactyla_plot <- artio_full %>% filter(Order == "Artiodactyla") %>% ggplot(., aes(x = Order, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
ruminantia_plot <- artio_full %>% filter(Suborder == "Ruminantia") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
whippomorpha_plot <- artio_full %>% filter(Suborder == "Whippomorpha") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/barplot_percentages.pdf", width = 10, height = 5, bg = "transparent")
grid.arrange(mammals_plot, artiodactyla_plot, ruminantia_plot, whippomorpha_plot, nrow = 1)
dev.off()

#replace their artiodactyla data with my artiodactyla data?
artio_full$Order <- "Amelia_artiodactyla"
artio_full$mammals <- "Mammals"
mammals_df <- rbind(new_mammals, artio_full)

mammals_df %>% filter(Order %in% c("Artiodactyla", "Amelia_artiodactyla")) %>% ggplot(., aes(x = Order, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
custom.colours <- c("#dd8ae7","#EECBAD", "#EECBAD" ,"#FC8D62", "gold", "#66C2A5", "green")
mammals_df %>% filter(Order %in% c("Artiodactyla", "Amelia_artiodactyla")) %>% ggplot(., aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())


# Section 6: Mammal tree ----------------------------------------------------

#make plot of bennie et al data with the artiodactyla data replaced with my own
new_mammals <- read.csv(here("sleepy_mammals.csv")) #data from Bennie et al, 2014
new_mammals <- new_mammals %>% filter(Order != "Artiodactyla")
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ] #318 species (83 cetaceans, 235 non cetaceans)
diel_full <- rbind(new_mammals, artio_full)
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,] #should be 3,102 species
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+3, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 6) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "clade_name", "sleepy_mammals", "_max_crep_plot_unlabelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 

order_list <- trait.data %>% group_by(Order) %>% filter(n()>10)
order_list <- unique(order_list$Order)
node_labels <- lapply(order_list, function(x){findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 2, taxonomic_level_name = x)})
node_labels <- do.call(rbind.data.frame, node_labels)
node_labels$barsize <- 2
node_labels$vjust <- 0.5
node_labels_left <- node_labels[node_labels$clade_name %in% c("Lagomorpha", "Scandentia", "Primates", "Artiodactyla", "Carnivora", "Perissodactyla", "Dermoptera", "Pholidota"),]
#node_labels_left[node_labels_left$clade_name %in% c("Dermoptera", "Pholidota"), "vjust"] <- -2
#node_labels_left[node_labels_left$clade_name %in% c("Perissodactyla"), "vjust"] <- -2
node_labels_right <- node_labels[!node_labels$clade_name %in% c("Lagomorpha", "Scandentia", "Primates", "Artiodactyla", "Carnivora", "Perissodactyla", "Dermoptera", "Pholidota"),]
#node_labels_right[node_labels_right$clade_name %in% c("Proboscidea", "Monotremata"), "vjust"] <- -2
node_labels_right[node_labels_right$clade_name %in% c("Cingulata"), "vjust"] <- 0

diel.plot <- ggtree(trpy_n, layout = "circular", fill = "transparent") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 6) + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot + geom_cladelab(barsize = 1.5, barcolor = "grey50", node = node_labels_left$node_number, label = node_labels_left$clade_name, hjust = 1, offset = 3, vjust = node_labels_left$vjust, offset.text = 2)
diel.plot <- diel.plot + geom_cladelab(barsize = 1.5, barcolor = "grey50", node = node_labels_right$node_number, label = node_labels_right$clade_name, offset = 3, vjust = node_labels_right$vjust, offset.text = 2)
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/mammals_max_crep_plot_cladelabels.pdf", bg = "transparent", width = 10, height = 10)
diel.plot
dev.off()

# Section 7: Phylogenetic signal -----------------
#to calculate the phylogenetic signal of discrete traits
# some people use delta statistic https://github.com/mrborges23/delta_statistic 
setwd(here())
source("scripts/Amelia_delta_code.R")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#trait.data <- read.csv(here("cetaceans_full.csv"))

trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- filter(trait.data, Suborder == "Ruminantia")

#keep only necessary columns
trait.data <- trait.data[, c("Species_name", "max_crep", "tips")]
trait.data <- trait.data[trait.data$max_crep %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
#this doesn't replace anything in the cetacean tree so comment it out for now
mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#vector of trait data, with species in same order as in tree (mam.tree$tip.label)
sps_order <- as.data.frame(mam.tree$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
trait.data <- merge(trait.data, sps_order, by = "tips")
trait.data <- trait.data[order(trait.data$id), ]
trait <- trait.data$max_crep

#convert to a vector
trait <- str_replace(trait, pattern = "cathemeral", replacement = "1")
trait <- str_replace(trait, pattern = "crepuscular", replacement = "2")
trait <- str_replace(trait, pattern = "diurnal", replacement = "3")
trait <- str_replace(trait, pattern = "nocturnal", replacement = "4")
trait <- as.numeric(trait)

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


#alternative method for calculating phylogenetic signal of a discrete trait
#install.packages("remotes")
#remotes::install_github("stoufferlab/phyloint", force = TRUE)
library(phyloint)
install.packages('motmot.2.0')
library(motmot.2.0)

#requires a vector of the trait data in the same order as phy$tip.label and the tree

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#trait.data <- read.csv(here("cetaceans_full.csv"))

trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- filter(trait.data, Suborder == "Whippomorpha")
#trait.data <- filter(trait.data, Suborder == "Ruminantia")

trait.data <- read.csv(here("Bennie_mam_data.csv"))
#trait.data <- filter(trait.data, Order == "Eulipotyphla")

#keep only necessary columns
trait.data <- trait.data[!duplicated(trait.data$tips), c("max_crep", "tips")]
trait.data <- trait.data[trait.data$max_crep %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #78 cetaceans in final tree
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
#this doesn't replace anything in the cetacean tree so comment it out for now
#mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#vector of trait data, with species in same order as in tree (mam.tree$tip.label)
sps_order <- as.data.frame(mam.tree$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
trait.data <- merge(trait.data, sps_order, by = "tips")
trait.data <- trait.data[order(trait.data$id), ]
trait <- trait.data$max_crep
names(trait) <- trait.data$tips

#function phylo.signal isn't working so load in the base code
phylo.signal <- function(trait, phy, rep = 999) {
  if (length(attributes(factor(trait))$levels) == length(trait)) 
    stop("Are you sure this variable is categorical?")
  
  phy <- keep.tip(phy, tip = names(trait))
  
  # calculate likelihood corresponding to maximum likelihood value of lambda
  obs <- fitDiscrete(phy, trait, transform="lambda")
  
  # calculate likelihood of model with no phylogenetic signal
  # this wasn't working so I tested the likelihood of the lambda against a tree with traits randomly distributed
  trait.random <- sample(trait)
  names(trait.random) <- names(trait)
  null <- fitDiscrete(phy, trait.random, transform = "lambda")
  
  # calculate the likelihood ratio between the two models
  LLR <- -2*(null$opt$lnL - obs$opt$lnL)
  
  # what is the p value of this likelihood ratio?
  p <- pchisq(LLR, df=1, lower.tail=FALSE)
  
  return(data.frame(row.names=NULL, lambda=obs$opt$lambda, obs=obs$opt$lnL, null=null$opt$lnL, LLR=LLR, p=p))
}

signal <- phylo.signal(trait = trait, phy = mam.tree, rep = 999)

#get the phylogenetic signal for all families with more than 100 species
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data <- trait.data[!duplicated(trait.data$tips), c("max_crep", "tips", "Order")]
trait.data <- trait.data[trait.data$max_crep %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #4228 mammals in final tree
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)
table(trait.data$Order)

#filter for orders that have over 100 species
trait.data <- trait.data %>% group_by(Order) %>% filter(n() > 100)

makeTraitVector <- function(trait.data = trait.data, Order_name = "Primates"){
  trait.data <- trait.data[trait.data$Order == Order_name,]
  sps_order <- as.data.frame(mam.tree$tip.label)
  colnames(sps_order) <- "tips"
  sps_order$id <- 1:nrow(sps_order)
  trait.data <- merge(trait.data, sps_order, by = "tips")
  trait.data <- trait.data[order(trait.data$id), ]
  trait <- trait.data$max_crep
  names(trait) <- trait.data$tips
  return(trait)
}
  
trait.vector.list <- lapply(unique(trait.data$Order), function(x) makeTraitVector(trait = trait.data, Order_name = x))
names(trait.vector.list) <- unique(trait.data$Order)

phylo.sig.list <- lapply(trait.vector.list, function(x) phylo.signal(trait = x, phy = mam.tree, rep = 3))

phylo.sig.df <- do.call(rbind.data.frame, phylo.sig.list)
phylo.sig.df$Species <- row.names(phylo.sig.df)

#i think chiroptera only has nocturnal species so the random reordering isn't any different than the actual data on the tree
ggplot(phylo.sig.df, aes(x = Species, y = lambda, fill = log(p))) + geom_bar(stat = "identity") + geom_text(aes(label = round(lambda, digits = 3)), vjust = -0.2)

ggplot(trait.data, aes(x = Order, fill = max_crep)) + geom_bar(position = "fill")
