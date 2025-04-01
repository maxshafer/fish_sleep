##Packages we will use ---------------------------------------------------
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
#open tree of life
library(rotl)
#adds timescale
library(deeptime)
#colours
library(RColorBrewer)
#apply two separate colour palettes
library(ggnewscale)
#more colours
library(pals)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#Load in the discrete traits 

##Data from Manger et al  2013, on cetacean body mass, brain mass, and social traits
#https://doi.org/10.1016/j.neuroscience.2013.07.041
#MBM – male body mass (g); FBM – female body mass (g); ABM – average body mass (g); BrM – brain mass (g); EQ – encephalization quotient; L – longevity (days); SMF – age at sexual maturity of females (days) 
#AGS – average group size in number of individuals; RGS – range of group size in number of individuals; GSD – group social dynamics; FS – feeding strategy. 
#Group Social Dynamics: Sol – mostly solitary lifestyle; SSG – small stable groups of less than 10 individuals; FF? – possible fission–fusion social dynamic; SG - stable groups of more than 10 individuals; FF – fission–fusion social dynamic. 
#Feeding strategies: SF – skim feeder; OCO – shows occasional co-operation in feeding; LF – lunge feeder; SwF – swallow feeder; CO – shows regular co-operation in feeding; SuF – suction feeder; R – raptorial feeder.
#Data for this table derived from the following sources: Nowak (1999), Lefebvre et al. (2006a), Manger (2006), Best et al. (2009), Perrin et al. (2009)

#Data from Coombs, 2021 (PhD thesis)
#https://discovery.ucl.ac.uk/id/eprint/10135933/7/Coombs_10135933_thesis_revised.pdf
#includes fossil species
#data on dentition, diet, echolocation, feeding mechanism, habitat

#Data from Churchill & Baltz, 2021
#https://doi.org/10.1111/joa.13522
# Includes museum specimens
# Data on orbit length (proxy for eye size), bizygomatic width (skull width, a proxy for body size)
#dive depth, dive duration, mass, body length
#first occurrence data (FAD) and last occurrence data (LAD) for fossil species

#Data from ??? et al, 20??
#Data on body size, diet, divetype, feeding behaviour, habitat, regime

#Data from Groot et al, 2023
#Data on lifespan, length, mass, brain mass, encephalization quotient,
#female and male reproductive age, lifespan, group size, sociality, group foraging,
#acoustic communication and more 


# Phylogenetic ANOVA tutorial ---------------------------------------------
#from http://blog.phytools.org/2017/11/print-method-for-phylanova-phylogenetic.html praise Levell

#simulate tree and data
tree<-pbtree(n=100)
Q<-matrix(c(-2,1,1,
            1,-2,1,
            1,1,-2),3,3)
rownames(Q)<-colnames(Q)<-letters[1:3]
# x is a discrete trait that can have one of three states (A, B, C)
x<-as.factor(sim.history(tree,Q)$states)
# y is a continuous trait, we generate values for it under fast brownian motion (y1)
#and under brownian motion with half the variance plus the value of x (creates a stronger phylogenetic signal between x and y)
y1<-fastBM(tree)
y2<-fastBM(tree,sig2=0.5)+as.numeric(x)

# case one we're looking at the effect of the discrete trait on continuous trait y1
#the nodes are labelled by the x state (A, B or C) and their y trait is shown on the y axis
#we can see that for case 1, trait x doesn't appear to associate with trait y (there seem to be As spread across all y values etc)
phenogram(tree,y1,ftype="off",col=make.transparent("blue",0.5))
tiplabels(pie=to.matrix(x,letters[1:3]),
          piecol=colorRampPalette(c("blue", "yellow"))(3),
          cex=0.4)

#we confirm this lack of association with the phylogenetic ANOVA
#F-value is 0.55 and Pr(>F) is 0.687, paired p-values are high (above 0.05) for all comparisons
#meaning none of the categories for x have an effect on y
phylANOVA(tree,x,y1)

#we can plot this as well
#even without phylogenetically correcting we can see the means and distributions are similar between groups
example_trait <- as.data.frame(y1 = y1, x = x)
ggplot(example_trait, aes(x = x, y = y1)) + geom_boxplot() + geom_point()

#case two, we're looking at continuous trait y1
#we can see more of an association in this phenogram
#species with a yellow trait x tend to have a higher trait y, blues have lower y and grey are intermediate
phenogram(tree,y2,ftype="off",col=make.transparent("blue",0.5))
tiplabels(pie=to.matrix(x,letters[1:3]),
          piecol=colorRampPalette(c("blue", "yellow"))(3),
          cex=0.4) 

#the phylogenetic ANOVA supports this
#the corrected p-values are low between all comparisons
phylANOVA(tree,x,y2)
#we can plot this as well
#even without phylogenetically correcting we can see the means and distributions are dissimilar between groups
example_trait <- as.data.frame(y1 = y2, x = x)
ggplot(example_trait, aes(x = x, y = y2)) + geom_boxplot() + geom_point()


phylANOVA(tree, x, y, nsim=1000, posthoc=TRUE, p.adj="holm")
#where tree is a phylo object, x is a vector containing the groups. y is a vector containing the response variable (continuously valued).
#nsim = an integer specifying the number of simulations (including the observed data).
#posthoc = a logical value indicating whether or not to conduct posthoc tests to compare the mean among groups.
#p.adj = method to adjust P-values for the posthoc tests to account for multiple testing. Options same as p.adjust.

# Prescence of cortistatin ------------------------------------------------

trait.data <- read.csv(here("cetaceans_full.csv"))

#create a phylogeny of just the cetaceans in the cortistatin paper
#Valente et al, 2021. https://doi.org/10.1016/j.ygeno.2020.11.002 
trait.data <- trait.data[trait.data$Species_name %in% c("Tursiops truncatus", "Tursiops aduncus", "Sousa chinesis", "Globicephala melas", "Peponocephala electra", "Lagenorhynchus obliquidens","Orcinus orca", 
                                                                   "Neophocaena asiaeorientalis", "Phocoena phocoena", "Phocoena sinus", "Delphinapterus leucas", "Monodon monoceros", "Inia geoffrensis", "Pontoporia blainvillei", "Lipotes vexillifer",
                                                                   "Platanista gangetica", "Mesoplodon bidens", "Ziphius cavirostris", "Kogia breviceps", "Physeter catodon", "Balaenoptera bonaerensis", "Balaenoptera acutorostrata", "Balaenoptera musculus",
                                                                   "Balaenoptera edeni", "Balaenoptera physalus", "Megaptera novaeangliae", "Eschrichtius robustus", "Eubalaena japonica", "Eubalaena glacialis"), ]

#add a column containing whether they have intact cort (missing Platanista gangetica and Sousa chinesis from the tree)
#just the type of mutation, should add info about where the mutation is
#do each of these have the same effect on the phenotype?
#cort_cetaceans$cort_144bp <- c("NA", "") need to think of how to do this
trait.data$cort_219bp <- c("1D", "1D", "NA", "NA", "NA", "1D", "NA", "NA", "NA", "NA", "1D", "1I", "NA", "PS", "1I", "NA", "NA", "NA", "1I", "NA", "PS", "NA", "NA", "NA", "NA", "NA", "NA")

#full list of species and cort status 
#cort_list<- c("Tursiops truncatus", "Tursiops aduncus", "Sousa chinesis", "Globicephala melas", "Lagenorhynchus obliquidens","Orcinus orca", "Neophocaena asiaeorientalis", "Phocoena phocoena", "Phocoena sinus", "Delphinapterus leucas", "Monodon monoceros", "Inia geoffrensis", "Pontoporia blainvillei", "Lipotes vexillifer", "Mesoplodon bidens", "Ziphius cavirostris", "Kogia breviceps", "Physeter catodon", "Balaenoptera bonaerensis", "Balaenoptera acutorostrata", "Balaenoptera musculus", "Balaenoptera edeni", "Balaenoptera physalus", "Megaptera novaeangliae", "Eschrichtius robustus", "Eubalaena japonica", "Eubalaena glacialis")
#cort_219bp_list <- c("1D", "1D", "1D", "NA", "NA", "NA", "1D", "NA", "NA", "NA", "NA", "1D", "1I", "NA", "PS", "1I", "NA", "NA", "NA", "1I", "NA", "PS", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

# cort_tr <- tol_induced_subtree(ott_ids = cort_cetaceans$ott_id[cort_cetaceans$flags %in% c("sibling_higher", "")], label_format = "id")
# cort_tr$tip.label <- cort_cetaceans$tips[match(cort_tr$tip.label, paste("ott", cort_cetaceans$ott_id, sep = ""))]
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), ]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)


#custom.colours3 <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
cort_diel_tree <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "cort_219bp")]
cort_diel_tree <- cort_diel_tree + geom_tile(data = cort_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_2),  inherit.aes = FALSE, color = "transparent")# + scale_fill_manual(values = custom.colours)
cort_diel_tree <- cort_diel_tree + geom_tiplab(color = "black", size = 2.5, offset = 3)
cort_diel_tree <- cort_diel_tree + geom_tile(data = cort_diel_tree$data[1:length(trpy_n$tip.label),], aes(y=y, x = x + 1.5, fill = cort_219bp), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = c('steelblue1', 'lightskyblue', "#dd8ae7", "pink", "#FC8D62", "grey","#66C2A5", "#A6D854", "royalblue"))
cort_diel_tree

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/cortistatin_cetaceans.png", width=18,height=15,units="cm",res=1200)
cort_diel_tree
dev.off()


# Coombs et al traits (echo, habitat, feeding mechanism) ---------------------------------------------
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]

#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), c("tips", "Diel_Pattern_2")]

echo <- read_xlsx("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_discrete_traits/Coombs_et_al_2021.xlsx")
echo <- echo[!(is.na(echo$Echo)), c("Museum ID", "Age", "Echo", "Diet", "Dentition", "FM", "Habitat")]
echo$tips <- str_replace(echo$`Museum ID`, " ", "_")

echo <- echo %>% filter(echo$tips %in% trait.data$tips)

trait.data <- merge(trait.data, echo, by = "tips")

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#should fix colours
custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours.2 <- c("grey35", "grey66")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Echo")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Echo), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Prescence of echolocation")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/echo_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#do the same for other discrete traits, like habitat
custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours.2 <- c( "grey90", "grey40", "grey66", "black")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Habitat")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Habitat), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Habitat")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/habitat_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#feeding mechanism
custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours.2 <- c("black", "grey40", "grey75")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "FM")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = FM), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Feeding mechanism")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/feeding_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

# Coombs et al, echolocation ----------------------------------------------
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_4)), c("tips", "Diel_Pattern_4")]

echo <- read_xlsx("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_discrete_traits/Coombs_et_al_2022.xlsx")
echo <- echo[!(is.na(echo$Echo)), c("Museum ID", "Echo", "Diet", "Dentition", "FM", "Habitat")]
echo$tips <- str_replace(echo$`Museum ID`, " ", "_")

echo <- echo %>% filter(echo$tips %in% trait.data$tips)

trait.data <- merge(trait.data, echo, by = "tips")

#add hippo data
trait.data <- rbind(trait.data, c("Hexaprotodon_liberiensis", "nocturnal","Hexaprotodon liberiensis", "no echo", "Misc", "Unknown", "Biting", "Semi-aquatic"))
trait.data <- rbind(trait.data, c("Hippopotamus_amphibius",  "nocturnal","Hippopotamus amphibius", "no echo","Misc", "Unknown", "Biting", "Semi-aqautic"))

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#should fix colours
#custom.colours <- c("#dd8ae7", "pink", "orange", "black", "grey80", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7", "#FC8D62", "black", "grey80", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_4", "Echo")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_4), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Echo), inherit.aes = FALSE, colour = "transparent") 
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3)
diel.plot


# Ancestral reconstruction of echolocation ------------------------------
#run the markov models on the echolocation data
#since 1 = echo and 2 = no echo we can use the following to set the state at the root to no echo
#"If analyzing a binary or multistate character, the order of root.p is the same order as the traits – e.g., for states 1, 2, 3, a root.p=c(0,1,0) would fix the root to be in state 2"
ER <- corHMM(phy = trpy_n, data = trait.data[, c("tips", "Echo")], rate.cat = 1, model = "ER", node.states = "marginal", root.p = c(0, 1))
SYM <- corHMM(phy = trpy_n, data = trait.data[, c("tips", "Echo")], rate.cat = 1, model = "SYM", node.states = "marginal", root.p = c(0, 1))
ARD <- corHMM(phy = trpy_n, data = trait.data[, c("tips", "Echo")], rate.cat = 1, model = "ARD", node.states = "marginal", root.p = c(0, 1))

echo_model_results <- lapply(c("ER", "SYM", "ARD"), function(x) eval(as.name(x)))
names(echo_model_results) <- paste(c("ER", "SYM", "ARD"), "_model", sep = "")
#can now plot the likelihood metrics using the Amelia all model plots new functions

model_results <- ARD
file_name <- "cetacean_echolocation"

phylo_tree <- model_results$phy

#rename column names for consistency in the next steps
colnames(model_results$data) <- c("tips", "Echo")

#to make more clear we can colour the tips separately using geom_tipppoint 
#may have to adjust what trait data column is called in each
base_tree <- ggtree(phylo_tree, layout = "rectangular") + geom_tiplab(size = 3, hjust = -0.1)
base_tree <- base_tree %<+% model_results$data[, c("tips", "Echo")]
base_tree <- base_tree + geom_tippoint(aes(color = Echo), size = 4) 
base_tree

#make the dataframe of likelihoods at the internal nodes without the tips
lik.anc <- as.data.frame(model_results$states)

#for cetaceans we have to add 72 because we are skipping the tips (nodes 1-72)
#the internal nodes start at 73 and end at node 143
#for artiodactyla we add 300 because we are skipping the tips (nodes 1-300)
#the internal nodes start at 301 and end at node 599
lik.anc$node <- c(1:nrow(lik.anc)) + nrow(model_results$data)

#get the pie charts from this database using nodepie
#the number of columns changes depending on how many trait states
pie <- nodepie(lik.anc, 1:(length(lik.anc)-1))

pie_tree <- base_tree + geom_inset(pie, width = .03, height = .03) 

#adds a timescale
pie_tree <- pie_tree + theme_tree2()
#reverses the timescale so it starts at 0mya at the tips and extends back to 50mya at ancestor
pie_tree <- revts(pie_tree)
pie_tree

#save out
png(paste("C:/Users/ameli/OneDrive/Documents/R_projects/New_ancestral_recon/pie_chart/", "pie_chart_anc_recon_", "echolocation", ".png", sep = ""), width=40,height=20, units="cm",res=600)
pie_tree
dev.off()

# is the proportion of nocturnal species higher in the echolocating species vs non-echolocating?
#no it seems about half are 
table(trait.data$Diel_Pattern_4, trait.data$Echo)

# Unknown (dive depth, body size)  -----------------------------------------------
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]

dive <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")
dive$tips <- dive$Taxon

dive <- dive %>% filter(dive$tips %in% trait.data$tips)

trait.data <- merge(trait.data, dive, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#dive depth
custom.colours <- c("#dd8ae7", "grey20", "#FC8D62", "grey40",  "#66C2A5", "#A6D854", "grey80", "black")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Divetype")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Divetype), inherit.aes = FALSE, colour = "transparent")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/dive_depth_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#use cutoffs to make a continuous trait (body size) into a discrete trait
trait.data$body_size <- "NA"

for(i in 1:length(trait.data$Body.size)){
  if(trait.data[i, 'Body.size'] >= 100){
    trait.data[i, 'body_size'] <- "small"
  }
  if(trait.data[i, 'Body.size'] >= 100 & trait.data[i, 'Body.size'] <= 500){
    trait.data[i, 'body_size'] <- "medium"
  }
  if(trait.data[i, 'Body.size'] >= 500 & trait.data[i, 'Body.size'] <= 1000){
    trait.data[i, 'body_size'] <- "large"
  }
  if(trait.data[i, 'Body.size'] >= 1000){
    trait.data[i, 'body_size'] <- "very large"
  }
}

#especially want to look at body size (see if it associates with cathemerality)
#body size
custom.colours <- c("#dd8ae7", "#FC8D62", "grey30", "grey70",  "#66C2A5", "#A6D854", "grey90", "black")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "body_size")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = body_size), inherit.aes = FALSE, colour = "transparent")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/body_size_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

custom.colours <- c("#dd8ae7","grey20", "#FC8D62", "grey30", "grey70", "grey40", "grey50", "#66C2A5", "#A6D854", "grey90", "black", "grey60", 'grey80', "grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Family")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Family), inherit.aes = FALSE, colour = "transparent")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/family_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()


# Manger et al traits  -------------------------------------------------
Manger <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Manger_2013_discrete_cet.csv")
#body mass, encephalization quotient, group size, social dynamics, longevity, feeding strategy
#of these, body mass is probably the most relevant to activity patterns (ie are whales with large body size cathemeral)
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), c("tips", "Diel_Pattern_2")]

colnames(Manger) <- c("Genus_species", "Male_mass", "Female_mass", "Average_body_mass", "Brain_mass", "Encephalization_quotient", "Longevity_days", "Sexual_maturity_days", "Group_size", "Group_size_range", "Group_social_dynamics", "Feeding_strategy")
Manger$tips <- str_replace(Manger$Genus_species, " ", "_")
Manger <- Manger %>% filter(Manger$tips %in% trait.data$tips)

trait.data <- merge(trait.data, Manger, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#remove commas and convert to numbers
trait.data$Average_body_mass <- str_replace_all(trait.data$Average_body_mass, pattern = ",", replacement = "")
trait.data$Average_body_mass <- as.integer(trait.data$Average_body_mass)
#remove empty rows
trait.data <- trait.data[!apply(trait.data == "", 1, all),]
trait.data <- trait.data[!is.na(trait.data$Average_body_mass),]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Average_body_mass")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label), ], aes(x=x +2, y=y, fill = Average_body_mass), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = "Average body mass")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/orbit_ratio_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()
# Groot et al traits () ---------------------------------------------

groot <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Groot_et_al_2023.xlsx")
#this data is coded so need to decode it with the original paper
#contains trait data that is in the other dataframes
#lifespan, length, mass, brain mass, EQ, age to reproduction, group size, gestation, sociality, group foraging, learned foraging, communication

# Fisher's exact test (not phylogenetically corrected) -----------------------------------------------------
#requires dataframe with both discrete traits of interest
#habitat vs diel pattern (Church data)
test_df <- trait.data[trait.data$Diel_Pattern_2 %in% c("cathemeral", "diurnal", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), c("Diel_Pattern_2", "Habitat1")]
test <- fisher.test(table(test_df))

mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 

#echo vs diel pattern (Coombs data)
test_df <- trait.data[trait.data$Diel_Pattern_2 %in% c("cathemeral", "diurnal", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), c("Diel_Pattern_2", "Echo")]
test <- fisher.test(table(test_df))
#this has slight more than 10 entries per unique combination (cathemeral echo has 15 entires, noc echo has 11 etc)
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 

#feeding mechanism vs diel pattern (Coombs data)
test_df <- trait.data[trait.data$Diel_Pattern_2 %in% c("cathemeral", "diurnal", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), c("Diel_Pattern_2", "FM")]
test <- fisher.test(table(test_df))
#this has slight more than 10 entries per unique combination (cathemeral echo has 15 entires, noc echo has 11 etc)
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 

# Churchill et al, traits: orbit ratio ---------------------------------------
#orbit information
church_pt1 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt1.xlsx")
#we are interested in orbit ratio, length of orbit relative to the bigozymatic width (proxy for body size)
#this is a continuous trait

trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2", "Parvorder")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), ]

#use below to rerun with max-crep four state model
# trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, pattern = "diurnal/crepuscular", replacement = "crepuscular")
# trait.data$Diel_Pattern_2 <- str_replace(trait.data$Diel_Pattern_2, pattern = "nocturnal/crepuscular", replacement = "crepuscular")

colnames(church_pt1) <- c("Species", "Family", "Specimen_number", "Left_orbit_length", "Right_orbit_length", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio")

church_pt1 <- church_pt1[!(is.na(church_pt1$Orbit_ratio)), c("Species", "Family", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio")]
church_pt1$tips <- str_replace(church_pt1$Species, " ", "_")

#filter the orbit size data by the species in our trait data 
church_pt1 <- church_pt1 %>% filter(church_pt1$tips %in% trait.data$tips)
#remove duplicates, leaves 46 species
church_pt1 <- church_pt1[!duplicated(church_pt1$Species),]

#could also average values for each species duplicates


trait.data <- merge(trait.data, church_pt1, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Orbit_ratio")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Orbit_ratio), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = "Orbit Ratio")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/orbit_ratio_vs_diel_.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
ggplot(trait.data, aes(x = Diel_Pattern, y = Axial_length)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Axial length") + scale_fill_manual(values=as.vector(polychrome(26))) #+ facet_wrap(~ Parvorder)

ggplot(trait.data, aes(x = Diel_Pattern, y = Corneal_diameter)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Corneal diameter") + scale_fill_manual(values=as.vector(polychrome(26))) #+ facet_wrap(~ Parvorder)

ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) #+ facet_wrap(~ Parvorder)

# Baker et al artiodactyla axial length -----------------------------------
artio_eyes <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Baker_2019.xlsx")
artio_eyes <- artio_eyes[2: nrow(artio_eyes),]
colnames(artio_eyes) <- c("tips", "Order", "Corneal_diameter", "Axial_length", "Activity_pattern", "Source")
artio_eyes <- filter(artio_eyes, Order == "Artiodactyla")
artio_eyes <- artio_eyes[artio_eyes$Corneal_diameter != "n/a", c("tips", "Corneal_diameter", "Axial_length", "Activity_pattern") ]
sleepy_artio <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_eyes <- artio_eyes %>% filter(artio_eyes$tips %in% sleepy_artio$tips)
artio_eyes$Transverse_diameter <- "NA"
artio_eyes <- artio_eyes %>% relocate(Transverse_diameter, .before= Activity_pattern)

artio_eyes <- merge(sleepy_artio, artio_eyes, by = "tips")
artio_eyes$Axial_length <- as.numeric(artio_eyes$Axial_length)
artio_eyes$Corneal_diameter <- as.numeric(artio_eyes$Corneal_diameter)
artio_eyes <- artio_eyes[, c("tips", "Species_name", "Family", "Order", "Diel_Pattern_2", "Corneal_diameter", "Axial_length", "Activity_pattern")]

#we can look at the ratio of corneal diameter to axial length
#From Hall et al: https://doi.org/10.1098/rspb.2012.2258
#The ratio of corneal diameter to axial length of the eye is a useful measure of relative sensitivity 
#and relative visual acuity that has been used in previous studies as a way to compare animals of disparate size
artio_eyes$Orbit_ratio <- artio_eyes$Corneal_diameter/artio_eyes$Axial_length
                        
#use below to rerun with max-crep four state model
artio_eyes$Diel_Pattern <- str_replace_all(artio_eyes$Diel_Pattern_2, pattern = "diurnal/crepuscular", replacement = "crepuscular")
artio_eyes$Diel_Pattern <- str_replace_all(artio_eyes$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
artio_eyes$Diel_Pattern <- str_replace_all(artio_eyes$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
ggplot(artio_eyes, aes(x = Diel_Pattern, y = Axial_length)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Axial length") + scale_fill_manual(values=as.vector(polychrome(26))) 

ggplot(artio_eyes, aes(x = Diel_Pattern, y = Corneal_diameter)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Corneal diameter") + scale_fill_manual(values=as.vector(polychrome(26))) 

ggplot(artio_eyes, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Corneal diameter: axial length ratio") + scale_fill_manual(values=as.vector(polychrome(26))) 


trait.data <- artio_eyes[artio_eyes$tips %in% mam.tree$tip.label,]
#two artio species dropped from the tree Camelus bactrianus and Lama glama
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2", "Orbit_ratio")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_2), inherit.aes = FALSE, colour = "transparent", width = 2) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2.5, y=y, fill = Orbit_ratio), inherit.aes = FALSE, colour = "transparent", width = 2) + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = "Orbit Ratio")
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3.5)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/artio_orbit_ratio_vs_diel_.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

# Phylogenetic signal of orbit size ---------------------------------------

#calculate pagels lambda and bloombergs k for cetacean orbit ratio using phylosig function
#not sure if the vector has to be in the same species order as the tr structure
#will order it this way just in case
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
trait.vector <- test$Orbit_ratio
names(trait.vector) <- test$tips

orbit_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list(), niter=10)
#k = 1.8925 for cetaceans
#p-value = 0.001 for cetaceans

#k = 0.348275 for non-cetacean artiodactyls
#p-value = 0.038
orbit_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                         control=list(), niter=10)
#lambda = 0.993593 for cetaceans
#p-value = 1.18 e-11 for cetaceans

#lambda = 0.964172
# p-value = 0.00269743

#phylosig gives the same result whether it is ordered or not but no harm in ordering it 


# Phylogenetic ANOVA for orbit size vs activity pattern ---------------------------------

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
ggplot(trait.data, aes(x = Diel_Pattern_2, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) #+ facet_wrap(~ Parvorder)

#we can make a phenogram of eye size evolution
#requires the trait data to be a named vector in the same species order as appears in trpy_n$tip.labels
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
#y trait is the continuous (response) variable 
trait.y <- test$Orbit_ratio
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

#use function phenogram from phytools
phenogram(trpy_n, trait.y,ftype="reg", spread.labels = TRUE, spread.cost=c(1,0))

#colour the tips by the activity pattern
tiplabels(pie = to.matrix(trait.x, unique(trait.data$Diel_Pattern)), piecol = c("peachpuff2", "#dd8ae7", "#66C2A5", "#FC8D62"), cex = 0.25)

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
orbit_ANOVA <- aov(Orbit_ratio ~ Diel_Pattern, data = trait.data)
summary(orbit_ANOVA) #not statistically significant for cetaceans, p value of F statistic = 0.297
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.00213

#phylogenetically corrected ANOVA
orbit_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant for cetaceans, p value of F statistic = 0.441
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.002

#plot by families (more than one species)
trait.data %>% group_by(Family) %>% filter(n()>2) %>% ggplot(., aes(x = Diel_Pattern_2, y = Orbit_ratio))+ geom_boxplot() + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + facet_wrap(~ Family) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26)))
#plot by families with only one species
trait.data %>% group_by(Family) %>% filter(n()<2) %>% ggplot(., aes(x = Diel_Pattern_2, y = Orbit_ratio))+ geom_boxplot() + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + facet_wrap(~ Family) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26)))

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/orbit_ratio_vs_diel_boxplot.png", width=20,height=15,units="cm",res=1200)
ggplot(trait.data, aes(x = Diel_Pattern_2, y = Orbit_ratio)) + geom_boxplot() + geom_point()
dev.off()

# # Churchill et al: dive depth -------------------------------------------
#dive depth, dive duration, body length, mass
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")

trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2", "Parvorder")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), ]

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern_2, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")

colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

#filter for species with dive data, leaves 28 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Dive_depth)), c("Species", "Dive_depth", "Dive_duration", "Mass", "Total_body_length", "Assymetry_index", "Peak_frequency")]

#correct species spelling to match
church_pt2[8, "Species" ] = "Berardius bairdii"
church_pt2[17, "Species" ] = "Lagenorhynchus obliquidens"

#read in the additional dive data I collected
url <- 'https://docs.google.com/spreadsheets/d/1_0ZS_tbddOCckkcKn9H5HpVRDZty4jhkUU20Nc0YYQY/edit?usp=sharing'
dive_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#average dive depth is very inconsistent (deep vs shallow dives, day vs night dives)
#and measured inconsistently (average of means across individuals, average of medians etc) so exclude for now
dive_full <- dive_full[!is.na(dive_full$Max_dive_depth_m),1:9]

#filter to rows with an alternative value
dive_full_alt <- dive_full[!is.na(dive_full$Alt_Max_1), ]
dive_full <- dive_full[is.na(dive_full$Alt_Max_1), ]

#pick the bigger value between max dive depth and alt dive depth
dive_full_alt[is.na(dive_full_alt)] <- 0
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_1)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_2)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_3)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_4)

#combine with original values
dive_full <- rbind(dive_full, dive_full_alt)
dive_full <- dive_full[, c("Species_name", "Max_dive_depth_m")]
colnames(dive_full) <- c("Species", "Max_dive_depth_m")
church_pt2 <- merge(church_pt2, dive_full, all.x = TRUE, all.y = TRUE)
church_pt2[is.na(church_pt2)] <- 0
church_pt2$Dive_depth <- pmax(church_pt2$Max_dive_depth_m, church_pt2$Dive_depth)
church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")
church_pt2[church_pt2 == 0] <- NA

#filter the dive data by the species in our trait data 
church_pt2 <- church_pt2 %>% filter(church_pt2$tips %in% trait.data$tips)

trait.data <- merge(trait.data, church_pt2, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trait.data <- trait.data[!is.na(trait.data$Dive_depth),]
trait.data$Dive_depth <- as.numeric(trait.data$Dive_depth)

#use below to look just at odontocetes or mysticetes
#trait.data <- trait.data %>% filter(Parvorder == "Odontoceti")
trait.data <- trait.data %>% filter(Parvorder == "Mysticeti")

#all 28 species with dive data are in the tree
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Dive_depth")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Dive_depth), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = "Dive depth")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/improved_dive_depth_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#plot out dive depth vs activity pattern (NOT phylogenetically corrected)
ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Parvorder), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) 

#identify the phylogenetic signal of dive depth
trait.vector <- trait.data$Dive_depth
names(trait.vector) <- trait.data$tips

dive_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                    control=list(), niter=10)
#K: 1.3837
#p-value = 0.001

#with additional species
#K: 0.766911 
#p-value = 0.001

dive_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                         control=list(), niter=10)
#lambda: 0.910314 
#p-value = 3.2622 e-6

#with additional species
#lambda: 0.910314 
#p-value = 1.03408e-05  

#optional: remove sperm whales because it skews everything
# trait.data <- trait.data[-(20),]
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#we can make a phenogram of dive depth evolution
#requires the trait data to be a named vector in the same species order as appears in trpy_n$tip.labels
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
#y trait is the continuous (response) variable 
trait.y <- test$Dive_depth
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

#use function phenogram from phytools
phenogram(trpy_n, trait.y, ftype="reg", spread.labels = TRUE, spread.cost=c(1,0))

#colour the tips by the activity pattern
tiplabels(pie = to.matrix(trait.x, unique(trait.data$Diel_Pattern)), piecol = c("#dd8ae7", "#FC8D62", "#66C2A5","peachpuff2"), cex = 0.25)

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
dive_ANOVA <- aov(Dive_depth ~ Diel_Pattern, data = trait.data)
summary(dive_ANOVA) #statistically significant, p value of F statistic = 0.0003
#statistically significant, p value of F statistic = 0.02

#phylogenetically corrected ANOVA
dive_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#statistically significant, p value of F statistic = 0.01
#not statistically significant, p value of F statistic = 0.227

# #Body mass from Churchill et all ----------------------------------------
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")

trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2", "Parvorder")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), ]

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern_2, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")

colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

#filter for species with body mass data, leaves 62 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Mass)), c("Species", "Dive_depth", "Dive_duration", "Mass", "Total_body_length", "Assymetry_index", "Peak_frequency")]

#correct species spelling to match
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus commersoni", replacement = "Cephalorhynchus commersonii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Berardius bairdi", replacement = "Berardius bairdii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias australis", replacement = "Lagenorhynchus australis")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias cruciger", replacement = "Lagenorhynchus cruciger")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias obliquidens", replacement = "Lagenorhynchus obliquidens")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus heavisidi", replacement = "Cephalorhynchus heavisidii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Leucopleurus acutus", replacement = "Lagenorhynchus_acutus")

church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")

#filter the body mass by the species in our trait data 
#leaves 56 species, removes phocoena dioptrica, lagenorhynchus cruciger and some mesoplodon sps
church_pt2 <- church_pt2 %>% filter(church_pt2$tips %in% trait.data$tips)

trait.data <- merge(trait.data, church_pt2, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#55 species with body mass are in the tree, Sousa_plumbea is not in the tree
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Mass")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Mass), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = "Body mass")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/body_size_vs_diel_.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#plot out body size vs activity pattern (NOT phylogenetically corrected)
ggplot(trait.data, aes(x = Diel_Pattern, y = log(Mass))) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Parvorder), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Log body mass") + scale_fill_manual(values=as.vector(polychrome(26))) + facet_wrap(~Parvorder)

#identify the phylogenetic signal of dive depth
trait.vector <- trait.data$Mass
names(trait.vector) <- trait.data$tips

mass_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                    control=list(), niter=10)
#K: 1.96692
#p-value = 0.001

mass_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                         control=list(), niter=10)
#lambda: 1.01668
#p-value = 3.62866e-18 

#we can make a phenogram of body mass evolution
#requires the trait data to be a named vector in the same species order as appears in trpy_n$tip.labels
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]

#filter to examine mysticetes and odontocetes separately
#test <- filter(test, Parvorder == "Mysticeti")  #only 8 species
test <- filter(test, Parvorder == "Odontoceti")
trpy_n <- keep.tip(mam.tree, tip = test$tips)

#y trait is the continuous (response) variable 
trait.y <- test$Mass
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

#use function phenogram from phytools
phenogram(trpy_n, log(trait.y), ftype="reg", spread.labels = TRUE, spread.cost=c(1,0))

#colour the tips by the activity pattern
tiplabels(pie = to.matrix(trait.x, unique(trait.data$Diel_Pattern)), piecol = c("#dd8ae7", "#FC8D62", "#66C2A5","peachpuff2"), cex = 0.25)

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
orbit_ANOVA <- aov(Mass ~ Diel_Pattern, data = test)
summary(orbit_ANOVA) #not statistically significant, p value of F statistic = 0.142

#phylogenetically corrected ANOVA
dive_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant, p value of F statistic = 0.393


# Section: Organize dive data ---------------------------------------------

#read in the Churchill et al data
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")
colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

#filter for species with dive data, leaves 62 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Dive_depth)), c("Species", "Dive_depth", "Dive_duration")]

#correct species spelling to match
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus commersoni", replacement = "Cephalorhynchus commersonii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Berardius bairdi", replacement = "Berardius bairdii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias australis", replacement = "Lagenorhynchus australis")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias cruciger", replacement = "Lagenorhynchus cruciger")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias obliquidens", replacement = "Lagenorhynchus obliquidens")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus heavisidi", replacement = "Cephalorhynchus heavisidii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Leucopleurus acutus", replacement = "Lagenorhynchus_acutus")

church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")

#read in the additional dive data I collected
url <- 'https://docs.google.com/spreadsheets/d/1_0ZS_tbddOCckkcKn9H5HpVRDZty4jhkUU20Nc0YYQY/edit?usp=sharing'
dive_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#average dive depth is very inconsistent (deep vs shallow dives, day vs night dives)
#and measured inconsistently (average of means across individuals, average of medians etc) so exclude for now
dive_full <- dive_full[!is.na(dive_full$Diel_Pattern_2),1:9]

#filter to rows with an alternative value
dive_full_alt <- dive_full[!is.na(dive_full$Alt_Max_1), ]
dive_full <- dive_full[is.na(dive_full$Alt_Max_1), ]

#pick the bigger value between max dive depth and alt dive depth
dive_full_alt[is.na(dive_full_alt)] <- 0
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_1)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_2)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_3)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_4)

#combine with original values
dive_full <- rbind(dive_full, dive_full_alt)
dive_full <- dive_full[, c("Species_name", "Max_dive_depth_m")]
colnames(dive_full) <- c("Species", "Max_dive_depth_m")
church_pt2 <- merge(church_pt2, dive_full, all.x = TRUE, all.y = TRUE)
church_pt2[is.na(church_pt2)] <- 0
church_pt2$Final_dive_depth <- pmax(church_pt2$Max_dive_depth_m, church_pt2$Dive_depth)
church_pt2$tips <- str_replace(church_pt2$Species, pattern = " ", replacement = "_")

#more data
Laeta_data <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Laeta_2020_dive_data.xlsx")
Laeta_data$tips <- str_replace(Laeta_data$Species, " ", "_")
Laeta_data <- Laeta_data[!is.na(Laeta_data$Species), c("Depth (m)", "tips")]
colnames(Laeta_data) <- c("Dive_depth_Laeta", "tips")
Laeta_data <- filter(Laeta_data, Dive_depth_Laeta != "-" )
Laeta_data$Dive_depth_Laeta <- as.integer(Laeta_data$Dive_depth_Laeta)

Dive_depth <- merge(church_pt2, Laeta_data, all.x = TRUE, all.y = TRUE)
Dive_depth[is.na(Dive_depth)] <- 0
Dive_depth$Dive_depth_final <- pmax(Dive_depth$Final_dive_depth, Dive_depth$Dive_depth_Laeta)

#save out only the relevant columns, leaves 64 species
Dive_depth <- Dive_depth %>% filter(Dive_depth_final != 0)
Dive_depth <- Dive_depth[, c(1,2,8)]
colnames(Dive_depth) <- c("tips", "Species", 'Dive_depth')
#four entries only have genus (Delphinus_sp, Lagenorhynchus_cruciger, Lissodelphis_sp, Neophocaena_sp,)

trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern_2)), c("tips", "Diel_Pattern_2", "Parvorder")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern_2 %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), ]

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern_2, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")

trait.data <- merge(trait.data, Dive_depth, by = "tips")
