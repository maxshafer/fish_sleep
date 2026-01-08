setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Select the diel dataset -----------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#uncomment whichever clade you want to plot
#clade_name <- "cetaceans_full"
#clade_name <- "sleepy_artiodactyla_full"
#clade_name <- "ruminants_full"
clade_name <- "whippomorpha"
#clade_name <- "whippomorpha_high_conf"
# clade_name <- "sleepy_artiodactyla_minus_cetaceans"

diel_full <- read.csv(here(paste0(clade_name, ".csv")))
#use below to remove NA species
diel_full <- diel_full[!is.na(diel_full$Diel_Pattern), ]

trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,]
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



# Section 4: Whippomorpha with family labels ----------------------------------------
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
nodes_right <- nodes_df[c("Balaenopteridae", "Balaenidae", "Hippopotamidae"),]

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "Family")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_cladelab(node = nodes_right$node_number, label = nodes_right$clade_name, offset=1.5, offset.text=2, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + geom_cladelab(node = nodes_left$node_number, label = nodes_left$clade_name, offset=1.5, offset.text=2, hjust = 1, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_cladelab(node = 71, label = "Eschrichtiidae", offset = 3, fontsize = 3)
diel.plot <- diel.plot + geom_cladelab(node = 45, label = "Iniiae", hjust = 1, offset = 3, fontsize = 3)
diel.plot <- diel.plot + geom_cladelab(node = 46, label = "Lipotidae", offset = 3, hjust = 1, fontsize = 3)
diel.plot <- diel.plot + geom_cladelab(node = 44, label = "Pontoporiidae", offset = 3, hjust = 1, fontsize = 3)
diel.plot <- diel.plot + geom_cladelab(node = 59, label = "Platanistidae", offset =3, hjust = 0.9, fontsize = 3)
diel.plot <- diel.plot + geom_cladelab(node = 141, label = "Kogiidae", offset=1.5, offset.text=2, hjust = 0.4, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + geom_cladelab(node = 60, label = "Physeteridae", offset = 3, fontsize = 3)
diel.plot <- diel.plot + geom_cladelab(node = 72, label = "Neobalaenidae", offset = 3,fontsize = 3)
diel.plot <- diel.plot + new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),] %>% filter(Family %in% c("Eschrichtiidae", "Iniiae", "Lipotidae", "Pontoporiidae", "Platanistidae", "Physeteridae", "Neobalaenidae")), aes(x=x+2.5, y=y, fill = Family), inherit.aes = FALSE, colour = "transparent", width = 1.3, height = 0.6) + scale_fill_manual(values = c(rep("grey50", 19)))
diel.plot 

ggtree(trpy_n, layout = "circular") + geom_tiplab() + geom_label(aes(label = node))

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

# Section 5: Artiodactyla with suborder labels ----------------------------
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
diel.plot <- diel.plot + geom_cladelab(node = nodes_right$node_number, label = nodes_right$clade_name, offset=1.5, offset.text=2, barsize=2, fontsize=4, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + geom_cladelab(node = nodes_left$node_number, label = nodes_left$clade_name, offset=1.5, offset.text=2, hjust = 1, barsize=2, fontsize=4, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/artio_with_suborders.pdf", width = 7, height = 7, bg = "transparent")
diel.plot
dev.off()


# Section 6: breakdown of % activity patterns proportion plots----------------------------
new_mammals <- read.csv(here("Bennie_mam_data.csv")) #data from Bennie et al, 2014
new_mammals <- new_mammals[!is.na(new_mammals$max_crep), ]

#add in my primary source data 
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ]

#we want to compare all mammals, vs all artiodactyla vs cetaceans/ruminants
new_mammals$mammals <- "Mammals"

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
mammals_plot <- ggplot(new_mammals, aes(x = mammals, fill = max_crep)) + geom_bar(position = "fill", width = 0.75) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank())
artiodactyla_plot <- artio_full %>% filter(Order == "Artiodactyla") %>% ggplot(., aes(x = Order, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())
#ruminantia_plot <- artio_full %>% filter(Suborder == "Ruminantia") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())
whippomorpha_plot <- artio_full %>% filter(Suborder == "Whippomorpha") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/barplot_percentages.pdf", width = 5, height = 2, bg = "transparent")
grid.arrange(mammals_plot, artiodactyla_plot, whippomorpha_plot, nrow = 1)
dev.off()

#replace their artiodactyla data with my artiodactyla data?
artio_full$Order <- "Amelia_artiodactyla"
artio_full$mammals <- "Mammals"
artio_full <- artio_full %>% select(Species_name, Order, Family, Diel_Pattern, max_crep)
new_mammals$Diel_Pattern <- new_mammals$max_crep
new_mammals <- new_mammals %>% select(Species_name, Order, Family, Diel_Pattern, max_crep)
mammals_df <- rbind(new_mammals, artio_full)

mammals_df %>% filter(Order %in% c("Artiodactyla", "Amelia_artiodactyla")) %>% ggplot(., aes(x = Order, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
custom.colours.3 <- c("#dd8ae7","#EECBAD", "#EECBAD" ,"#FC8D62", "gold", "#66C2A5", "green")
mammals_df %>% filter(Order %in% c("Artiodactyla", "Amelia_artiodactyla")) %>% ggplot(., aes(x = Order, fill = Diel_Pattern)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours.3) + theme_minimal() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())


# Section 7: Mammal tree ----------------------------------------------------

#make plot of bennie et al data with the artiodactyla data replaced with my own
new_mammals <- read.csv(here("Bennie_mam_data.csv")) #data from Bennie et al, 2014, 4732 species
new_mammals <- new_mammals %>% filter(Order != "Artiodactyla") #4492 species (removes 240 artios)
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv")) 
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ] #317 species (82 cetaceans, 235 non cetaceans)
artio_full <- artio_full %>% select(Species_name, Order, Family, max_crep)
new_mammals <- new_mammals %>% select(Species_name, Order, Family, max_crep)
diel_full <- rbind(new_mammals, artio_full) #4809 species
diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,] #should be 4,400 species in tree (other 400 misnamed)
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+3, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 6) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/clade_name/sleepy_mammals_max_crep_plot_unlabelled.pdf", width = 9, height = 8, bg = "transparent")
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

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/mammals_max_crep_plot_cladelabels.pdf", bg = "transparent", width = 10, height = 10)
# diel.plot
# dev.off()

# Section 8: Phylogenetic signal Delta statistic -----------------
#to calculate the phylogenetic signal of discrete traits
# some people use delta statistic https://github.com/mrborges23/delta_statistic 
setwd(here())
source("scripts/Amelia_delta_code.R")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
#trait.data <- filter(trait.data, Suborder == "Whippomorpha")

#keep only necessary columns
trait.data <- trait.data[!is.na(trait.data$max_crep), c("Species_name", "max_crep", "tips")]

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #78 species for cetaceans
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#vector of trait data, with species in same order as in tree (mam.tree$tip.label)
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

random_delta <- rep(NA,100) #create a list of 100 NA values to fill in later
for (i in 1:100){
  rtrait <- sample(trait) #randomly place samples across the tree
  random_delta[i] <- delta(rtrait,mam.tree,0.1,0.0589,10000,10,100) #calculate the delta for each of these simulated trees
}

#compute the probability p(random_delta>deltaA) in the null distribution, which returns the p-value.
p_value <- sum(random_delta>delta_diel)/length(random_delta) #is there a significant difference between the real delta and the simulated ones?
p_value

#visualize how the calculated delta compares to the simulated ones
random_delta_df <- data.frame(random_delta)
ggplot(random_delta_df, aes(y = random_delta)) + geom_boxplot() + geom_hline(yintercept = delta_diel, color = "red", size = 2)

#if p-value less than 0.05 there is evidence of phylogenetic signal between the trait and character
#p value is 0.04, so there is a phylogenetic signal for diel activity patterns in cetaceans


#calculate delta statistic for all major mammal orders
calculateDelta2 <- function(trait.data = trait.data, order_name = "Artiodactyla"){
  mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
  trait.data <- filter(trait.data, Order == order_name)
  trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
  mam.tree.trimmed <- keep.tip(mam.tree, tip = trait.data$tips)
  sps_order <- as.data.frame(mam.tree.trimmed$tip.label)
  colnames(sps_order) <- "tips"
  sps_order$id <- 1:nrow(sps_order)
  trait.data <- merge(trait.data, sps_order, by = "tips")
  trait.data <- trait.data[order(trait.data$id), ]
  trait <- trait.data$max_crep
  delta_diel <- delta(trait, mam.tree.trimmed, 0.1, 0.0589, 10000, 10, 100)
  random_delta <- rep(NA,10) #create a list of 100 NA values to fill in later
  for (i in 1:10){
    rtrait <- sample(trait) #randomly place samples across the tree
    random_delta[i] <- delta(rtrait,mam.tree.trimmed,0.1,0.0589,10000,10,100) #calculate the delta for each of these simulated trees
  }
  p_value <- sum(random_delta>delta_diel)/length(random_delta)
  return(list(delta_diel, p_value))
}

trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data <- trait.data[, c("max_crep", "tips", "Order")]
trait.data.1 <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep), c("max_crep", "tips", "Order")]
trait.data.1$Order <- str_replace(trait.data.1$Order, pattern = "Artiodactyla", replacement = "Amelia_artiodactyla")
trait.data <- rbind(trait.data, trait.data.1) #4459 species

#minimum number they use is 20
trait.data <- trait.data %>% group_by(Order) %>% filter(length(unique(max_crep)) ==4)# %>% filter(n() > 20)
table(trait.data$Order)

#lapply not working so do each separately
a <- calculateDelta2(trait.data, order_name = "Cingulata")
b <- calculateDelta2(trait.data, order_name = "Amelia_artiodactyla")
c <- calculateDelta2(trait.data, order_name = "Artiodactyla")
d <- calculateDelta2(trait.data, order_name = "Didelphimorphia")
e <- calculateDelta2(trait.data, order_name = "Lagomorpha")
f <- calculateDelta2(trait.data, order_name = "Primates")
g <- calculateDelta2(trait.data, order_name = "Rodentia")
h <- calculateDelta2(trait.data, order_name = "Perissodactyla")
i <- calculateDelta2(trait.data, order_name = "Carnivora")

delta_list <- list(a,b,c,d,e,f,g,h,i)
delta_df <- do.call(rbind.data.frame, delta_list)
colnames(delta_df) <- c("delta_statistic", "p_value")
delta_df$clade <- c("Cingulata", "Amelia_artiodactyla", "Artiodactyla", "Didelphimorphia", "Lagomorpha", "Primates", "Rodentia", "Perissodactyla")

write.csv(delta_df, here("delta_df.csv"), row.names = FALSE)

trait.data <- trait.data %>% group_by(Order) %>% filter(length(unique(max_crep)) ==4) %>% filter(n() > 20) %>% filter()
delta_list <- lapply(unique(trait.data$Order), function(x) calculateDelta2(trait.data, order_name = x))
delta_df <- do.call(rbind.data.frame, delta_list)
delta_df$clade <- c("Cingulata", "Amelia_artiodactyla", "Artiodcatyla", "Didelphimorphia", "Lagomorpha", "Primates", "Rodentia")

# Section 9: Phylogenetic signal lambda -wrong  ----------

#requires a vector of the trait data in the same order as phy$tip.label
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#function phylo.signal isn't working so load in the base code
#rep doesn't do anything?
phylo.signal <- function(trait, phy, rep = 999) {
  if (length(attributes(factor(trait))$levels) == length(trait)) 
    stop("Are you sure this variable is categorical?")
  
  phy <- keep.tip(phy, tip = names(trait))
  
  # calculate likelihood corresponding to maximum likelihood value of lambda
  obs <- fitDiscrete(phy, trait, transform="lambda")
  
  # calculate likelihood of model with no phylogenetic signal
  #null <- fitDiscrete(transform(phylo, "lambda", 0), trait)
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

#function to create a vector of trait data, with species in same order as in tree (mam.tree$tip.label)
makeTraitVector <- function(trait.data = trait.data, taxonomic_level = "Order", clade_name = "Primates"){
  trait.data <- trait.data[trait.data[,taxonomic_level] == clade_name, ]
  sps_order <- as.data.frame(mam.tree$tip.label)
  colnames(sps_order) <- "tips"
  sps_order$id <- 1:nrow(sps_order)
  trait.data <- merge(trait.data, sps_order, by = "tips")
  trait.data <- trait.data[order(trait.data$id), ]
  trait <- trait.data$max_crep
  names(trait) <- trait.data$tips
  return(trait)
}

#calculate signal for all mammals
trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data$Kingdom <- "Mammals"
trait <- makeTraitVector(trait.data = trait.data, taxonomic_level = "Kingdom", clade_name = "Mammals")
mam_signal <- phylo.signal(trait = trait, phy = mam.tree, rep = 999)
mam_signal$clade <- "Mammals"
row.names(mam_signal) <- "Mammals"

#calculate signal for artiodactyla suborders
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.vector.list <- lapply(unique(trait.data$Suborder), function(x) makeTraitVector(trait = trait.data, taxonomic_level = "Suborder", clade_name = x))
names(trait.vector.list) <- unique(trait.data$Suborder)
phylo.sig.list <- lapply(trait.vector.list, function(x) phylo.signal(trait = x, phy = mam.tree, rep = 999))
suborder_df <- do.call(rbind.data.frame, phylo.sig.list)
suborder_df$clade <- row.names(suborder_df)

#write.csv(suborder_df, here("phylogenetic_signal_artio_suborder.csv"), row.names = FALSE)

#calculate phylogenetic signal for all families with more than 100 species
trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data <- trait.data[, c("max_crep", "tips", "Order")]

#to use my artio data instead (or in addition to)
#trait.data <- filter(trait.data, Order != "Artiodactyla")
trait.data.1 <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep), c("max_crep", "tips", "Order")]
trait.data.1$Order <- str_replace(trait.data.1$Order, pattern = "Artiodactyla", replacement = "Amelia_artiodactyla")

trait.data <- rbind(trait.data, trait.data.1) #4459 species

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #4303 mammals in final tree
table(trait.data$Order, trait.data$max_crep)

#filter for orders that have over x number of species
#trait.data <- trait.data %>% group_by(Order) %>% filter(n() > 70)

#filter for species with all four diel categories
trait.data <- trait.data %>% group_by(Order) %>% filter(length(unique(max_crep)) ==4)

trait.vector.list <- lapply(unique(trait.data$Order), function(x) makeTraitVector(trait = trait.data, taxonomic_level = , clade_name = x))
names(trait.vector.list) <- unique(trait.data$Order)

phylo.sig.list <- lapply(trait.vector.list, function(x) phylo.signal(trait = x, phy = mam.tree))
phylo.sig.df <- do.call(rbind.data.frame, phylo.sig.list)
phylo.sig.df$clade <- row.names(phylo.sig.df)

#save out 
#write.csv(phylo.sig.df, here("phylogenetic_signal_mammals.csv"), row.names = FALSE)

#final comparison
phylo_signal_final <- rbind(phylo.sig.df, suborder_df, mam_signal)
write.csv(phylo_signal_final, here("phylogenetic_signal.csv"), row.names = FALSE)

#i think chiroptera only has nocturnal species so the random reordering isn't any different than the actual data on the tree
ggplot(phylo_signal_final, aes(x = clade, y = lambda, fill = log(p))) + geom_bar(stat = "identity") + geom_text(aes(label = round(lambda, digits = 3)), vjust = -0.2)

knitr::kable(phylo_signal_final, format = "html", digits = 3, caption = "Table X") %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("phylosig_table_longer.html")
webshot("phylosig_table_longer.html", file = "C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/phylogenetic_signal.pdf")


# Section 9.5: Phylogenetic signal M statistic -----------------
#install.packages("phylosignalDB")
library(phylosignalDB)

#calculate for mammals
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("Bennie_mam_data.csv"))

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] 
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)
trait_df <- data.frame(B1 = as.factor(trait.data$max_crep), row.names = trait.data$tips)
trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
mam_signal <- phylosignal_M(trait_dist, mam.tree, reps = 99) # reps=999 better

#calculate for each major order
calculatePhylosig <- function(trait.data = trait.data, order_name = "Artiodactyla"){
  mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
  trait.data.filtered <- filter(trait.data, Order == order_name)
  trait.data.filtered <- trait.data.filtered[trait.data.filtered$tips %in% mam.tree$tip.label,]
  mam.tree.trimmed <- keep.tip(mam.tree, tip = trait.data.filtered$tips)
  trait_df <- data.frame(B1 = as.factor(trait.data.filtered$max_crep), row.names = trait.data.filtered$tips)
  trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
  signal <- phylosignal_M(trait_dist, mam.tree.trimmed, reps = 99)
  return(signal)
}

trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data <- trait.data[, c("tips", "max_crep", "Order")]
trait.data.art <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), c("tips", "max_crep", "Order")]
trait.data.art$Order <- str_replace(trait.data.art$Order, pattern = "Artiodactyla", replacement = "Amelia_artiodactyla")
trait.data <- rbind(trait.data, trait.data.art)
#filter for orders with more than one activity pattern or 4 activity patterns
trait.data <- trait.data %>% group_by(Order) %>% filter(length(unique(max_crep)) ==4)
table(trait.data$Order)
#filter for orders with a cutoff number of species
#trait.data <- trait.data %>% group_by(Order) %>% filter(n() > 20) 
signal_list <- lapply(unique(trait.data$Order), function(x) calculatePhylosig(trait.data, order_name = x))
names(signal_list) <- unique(trait.data$Order)
signal_list_df <- do.call(rbind.data.frame, signal_list)
signal_list_df$clade <- row.names(signal_list_df)

write.csv(signal_list_df, here("M_statistic_df.csv"), row.names = FALSE)

#calculate for artio
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #78 species for cetaceans
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)
trait_df <- data.frame(B1 = as.factor(trait.data$max_crep), row.names = trait.data$tips)
trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
artio_signal <- phylosignal_M(trait_dist, mam.tree, reps = 99) # reps=999 better

#calculate for whippo
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- filter(trait.data, Suborder == "Whippomorpha")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #78 species for cetaceans
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)
trait_df <- data.frame(B1 = as.factor(trait.data$max_crep), row.names = trait.data$tips)
trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
whippo_signal <- phylosignal_M(trait_dist, mam.tree, reps = 99) # reps=999 better

#calculate for ruminants
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data <- filter(trait.data, Suborder == "Ruminantia")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #78 species for cetaceans
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)
trait_df <- data.frame(B1 = as.factor(trait.data$max_crep), row.names = trait.data$tips)
trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
ruminant_signal <- phylosignal_M(trait_dist, mam.tree, reps = 99) # reps=999 better

#save out table of phylogenetic signal
delta_df <- read.csv(here("delta_df.csv"))
M_stat_df <- read.csv(here("M_statistic_df.csv"))

final_signal_df <- merge(M_stat_df, delta_df, by = "clade", all = TRUE) 
colnames(final_signal_df) <- c("Order", "M statistic", "p-value", "Delta statistic", "p-value ")

knitr::kable(final_signal_df, format = "html", digits = 3, caption = "Table X") %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("phylosig_table_final.html")
webshot("phylosig_table_final.html", file = "C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/phylogenetic_signal_final.pdf")

# Section 10: artio tree with whippo collapsed -----------------------------

diel_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
diel_full <- diel_full[!is.na(diel_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#label major families
suborders <- trait.data %>% filter(Suborder != "Whippomorpha") 

#For loop to find the labels for each suborder
nodes_list <- list()
for(i in unique(suborders$Suborder)){
  node_df <- findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = i)
  nodes_list[[i]] <- node_df
  nodes_df <- do.call(rbind, nodes_list)
}

nodes_left <- nodes_df[c("Ruminantia"),]
nodes_right <- nodes_df[c("Tylopoda", "Suina"),]

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- scaleClade(diel.plot, 511, .05) %>% collapse(node = 511, 'mixed', fill = "grey40") #whippo is at node 511
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_cladelab(node = nodes_right$node_number, label = nodes_right$clade_name, offset=1.5, offset.text=2, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + geom_cladelab(node = nodes_left$node_number, label = nodes_left$clade_name, offset=1.5, offset.text=2, hjust = 1, barsize=2, fontsize=3, barcolour = "grey50", textcolour = "black")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_cladelabel(node = 511, label = "Whippomorpha", offset = 60, hjust = 0.35, size = 0.5, barsize = 2, fontsize = 3)
diel.plot

diel.plot <- ggtree(trpy_n, layout = "circular") 
diel.plot <- diel.plot + geom_text(aes(label = node), colour= "blue")
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/artio_collapsed.pdf", width = 7, height = 7, bg = "transparent")
diel.plot
dev.off()


# Section 11: mammal with collapsed node ----------------------------------

new_mammals <- read.csv(here("Bennie_mam_data.csv")) #data from Bennie et al, 2014, 4477 species
new_mammals <- new_mammals %>% filter(Order != "Artiodactyla") #4242 species (removes 240 artios)
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv")) 
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ] #317 species (82 cetaceans, 235 non cetaceans)
artio_full <- artio_full %>% select(Species_name, Order, Family, max_crep)
new_mammals <- new_mammals %>% select(Species_name, Order, Family, max_crep)
diel_full <- rbind(new_mammals, artio_full) #4559 species
diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
missing_species <- diel_full[!diel_full$tips %in% mam.tree$tip.label,] #256 missing
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,] #should be 4,303 species in tree
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

# trait.data <- trait.data %>% filter(Order %in% c("Chiroptera", "Soricomorpha", "Cingulata", "Erinaceomorpha", "Pilosa", "Macroscelidae", "Afrosoricida")) 
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

order_list <- trait.data %>% group_by(Order) %>% filter(n()> 10)

node_labels <- lapply(unique(order_list$Order), function(x){findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 2, taxonomic_level_name = x)})
node_labels <- do.call(rbind.data.frame, node_labels)
node_labels_left <- node_labels[node_labels$clade_name %in% c("Lagomorpha", "Scandentia", "Primates", "Carnivora", "Dermoptera", "Pholidota", "Chiroptera"),]
node_labels_right <- node_labels[!node_labels$clade_name %in% c("Lagomorpha", "Scandentia", "Primates", "Carnivora", "Dermoptera", "Pholidota","Chiroptera", "Artiodactyla"),]

diel.plot <- ggtree(trpy_n, layout = "circular", size = 0.5, colour = "grey30") 
diel.plot <- scaleClade(diel.plot, 4877, .05) %>% collapse(node = 4877, 'mixed', fill = "grey40") #artio is at node 2717
diel.plot <- diel.plot %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 6) + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot + geom_cladelab(barsize = 1.5, barcolor = "grey50", node = node_labels_left$node_number, label = node_labels_left$clade_name, hjust = 1, offset =3, offset.text = 4, fontsize = 4)
diel.plot <- diel.plot + geom_cladelab(barsize = 1.5, barcolor = "grey50", node = node_labels_right$node_number, label = node_labels_right$clade_name, offset = 3, offset.text = 4, fontsize= 4)
#diel.plot <- diel.plot + geom_cladelab(node = , label = "Eschrichtiidae", offset = 3, fontsize = 3)
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = 'transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_cladelabel(node = 4877, label = "Artiodactyla", offset.text = 80, vjust = 1, barsize = 2, fontsize = 4)
diel.plot

# node_tree <- ggtree(trpy_n, layout = "circular") + geom_tiplab(size = 1)
# node_tree <- node_tree + geom_label(aes(label = node), size = 1)
# node_tree

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/mammals_collapsed.pdf", bg = "transparent", width = 10, height = 10)
diel.plot
dev.off()


# Section 12: Full table of activity patterns (supplemental) --------------

artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"),)
artio_full <- artio_full[, c("Species_name", "Suborder", "Diel_Pattern")]

#add in NA artio species since they got removed in the previous step
url <- 'https://docs.google.com/spreadsheets/d/1JGC7NZE_S36-IgUWpXBYyl2sgnBHb40DGnwPg2_F40M/edit?gid=562902012#gid=562902012'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
diel_full <- diel_full[is.na(diel_full$Confidence_primary_source), c("Species_name", "Diel_Pattern_primary")]
diel_full[diel_full == ""] <- NA
diel_full$Suborder <- "Ruminantia" #they are all ruminants
colnames(diel_full) <- c("Species_name", "Diel_Pattern", "Suborder")
diel_full <- diel_full %>% relocate(Diel_Pattern, .after = Suborder)

artio_full <- rbind(artio_full, diel_full)
artio_full$Diel_Pattern <- str_to_title(artio_full$Diel_Pattern)

test <- artio_full %>% count(Diel_Pattern)
test1 <- artio_full %>% filter(Suborder == "Whippomorpha") %>% count(Diel_Pattern)
test2 <- artio_full %>% filter(Suborder == "Ruminantia") %>% count(Diel_Pattern)
test3 <- artio_full %>% filter(Suborder == "Suina") %>% count(Diel_Pattern) #there are only 20sps
test4 <- artio_full %>% filter(Suborder == "Tylopoda") %>% count(Diel_Pattern) #there are only 6sps

test <- rbind(test, test1, test2, test3, test4)
test$Clade <- c(rep("All artiodactyla", 7), rep("Whippomorpha", 7), rep("Ruminantia", 7), rep("Suina", 6), rep("Tylopoda", 2))
test[is.na(test)] <- "Unknown"
test <- test %>% pivot_wider(names_from = Diel_Pattern, values_from = n)
test[is.na(test)] <- 0

knitr::kable(test, format = "html", digits = 3, caption = "Table 1") %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("table_final.html")
webshot("table_final.html", file = "C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Supplemental_diel_pattern_chart.pdf")

#with max crep categories
artio_full$Diel_Pattern <- str_replace(artio_full$Diel_Pattern, pattern = "Diurnal/Crepuscular", replacement = "Crepuscular")
artio_full$Diel_Pattern <- str_replace(artio_full$Diel_Pattern, pattern = "Cathemeral/Crepuscular", replacement = "Crepuscular")
artio_full$Diel_Pattern <- str_replace(artio_full$Diel_Pattern, pattern = "Nocturnal/Crepuscular", replacement = "Crepuscular")

test <- artio_full %>% count(Diel_Pattern)
test1 <- artio_full %>% filter(Suborder == "Whippomorpha") %>% count(Diel_Pattern)
test2 <- artio_full %>% filter(Suborder == "Ruminantia") %>% count(Diel_Pattern)
test3 <- artio_full %>% filter(Suborder == "Suina") %>% count(Diel_Pattern) #there are only 20sps
test4 <- artio_full %>% filter(Suborder == "Tylopoda") %>% count(Diel_Pattern) #there are only 6sps

test <- rbind(test, test1, test2, test3, test4)
test$Clade <- c(rep("All artiodactyla", 5), rep("Whippomorpha", 5), rep("Ruminantia", 5), rep("Suina", 4), rep("Tylopoda", 2))
test[is.na(test)] <- "Unknown"
test <- test %>% pivot_wider(names_from = Diel_Pattern, values_from = n)
test[is.na(test)] <- 0

knitr::kable(test, format = "html", digits = 3, caption = "Table 1") %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("table_final.html")
webshot("table_final.html", file = "C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Diel_pattern_chart.pdf")
