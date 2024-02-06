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


# Section 1.0 Format and combine cetacean and Cox artiodactyla data ------------------------------

#we still need to load in the cetacean sleep google sheets
#will combine this into one dataframe with the rest of artiodactyla

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit#gid=0'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
#check if all the correct diel patterns are there, cycle through pattern 2 and 3 as well
table(cetaceans_full$Diel_Pattern_1)

#they use a different name for sperm whales in the mam tree so we can change that
cetaceans_full$Species_name <- str_replace(cetaceans_full$Species_name, "Physeter catodon", "Physeter macrocephalus")
#only need these first rows (notes and sources aren't important)
cetaceans_full <- cetaceans_full[,1:10]
cetaceans_full$diel <- tolower(cetaceans_full$Diel_Pattern_3)

## Remove species without diel data
cetaceans_full <- cetaceans_full[!(is.na(cetaceans_full$diel)),]

#load in Cox mammal trait data (diel activity patterns)
mammals_full <- readRDS("trait_data_mammals.rds")

#subset for just cetartiodactyla
artiodactyla <- subset(mammals_full, Order == "Cetartiodactyla")
# This data doesn't find any artio species to be cathemeral or crepuscular. Is this accurate??
# So all three diel columns are interchangeable

#edit cetacean dataframe so it has similar format to add it to artio dataframe
# diel 1 is only di or noc
# diel 2 is di, noc, di/crep, noc/crep, cathemeral
# diel 3 is di, noc, crep, cath
#since there are no crepuscular or cathemeral artiodactyla can match these up with any of their columns

cetaceans_shorter <- subset(cetaceans_full, select = c(Species_name, Diel_Pattern_1, Diel_Pattern_2, Diel_Pattern_3, Parvorder))
cetaceans_shorter$Order <- "Cetartiodactyla"
cetaceans_shorter$Genus <- cetaceans_shorter$Species_name
#difficult to select just the first word of the name for the genus column
cetaceans_shorter <- relocate(cetaceans_shorter, Order, .after = Species_name)
cetaceans_shorter <- relocate(cetaceans_shorter, Parvorder, .after = Order)
cetaceans_shorter <- relocate(cetaceans_shorter, Genus, .after = Parvorder)
cetaceans_shorter$Species_name <- str_replace(cetaceans_shorter$Species_name, " ", "_")
cetaceans_shorter$Diel_Pattern_1 <- tolower(cetaceans_shorter$Diel_Pattern_1)
cetaceans_shorter$Diel_Pattern_2 <- tolower(cetaceans_shorter$Diel_Pattern_2)
cetaceans_shorter$Diel_Pattern_3 <- tolower(cetaceans_shorter$Diel_Pattern_3)

colnames(cetaceans_shorter) <- c("species", "Order", "Family", "Genus", "diel", "diel1", "diel2")

#add cetaceans data to the rest of artiodactyla
artiodactyla_full <- rbind(artiodactyla, cetaceans_shorter)
row.names(artiodactyla_full) <- artiodactyla_full$species
colnames(artiodactyla_full) <- c("species", "Order", "Family", "Genus", "diel1", "diel2", "diel3")

# Section 1.1a Binary, Cox, Artiodactyla model ------------------------------------

## Read in the mammalian phylogeny
#mammal_trees <- read.nexus("Cox_mammal_data/Complete_phylogeny.nex")
#mam.tree <- maxCladeCred(mammal_trees, tree = TRUE) #this function takes a long time to run
# Save out maxcladecred, so we don't have to recalculate it every time
#saveRDS(mam.tree,"maxCladeCred_mammal_tree.rds")
#keep commented out unless you need to recalculate it, will be saved otherwise
mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")

#we need to drop the missing species to create the tree
all_artiodactyla <- artiodactyla_full[artiodactyla_full$species %in% mam.tree$tip.label,]
trpy_n_cox1 <- keep.tip(mam.tree, tip = all_artiodactyla$species)
#tree of all the artio species with data in the mammal tree (shown without trait data)
ggtree(trpy_n_cox1, layout = "fan", size = 0.5) + geom_tiplab(size = 1.5)

#isolate just the diurnal and nocturnal species (198 species)
cox.trait.data1 <- artiodactyla_full[!(is.na(artiodactyla_full$diel1)),]
# selects only data that is in the mammal tree (190 species)
cox.trait.data1 <- cox.trait.data1[cox.trait.data1$species %in% mam.tree$tip.label,]
row.names(cox.trait.data1) <- cox.trait.data1$species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_cox1 <- keep.tip(mam.tree, tip = cox.trait.data1$species)

#now that we have a tree with the di/noc exclusive data we can run markov models on this data
trait.vector.cox1 <- cox.trait.data1$diel1
ace_cox_er1 <- ace(trait.vector.cox1, trpy_n_cox1, model = "ER", type = "discrete")
ace_cox_sym1 <- ace(trait.vector.cox1, trpy_n_cox1, model = "SYM", type = "discrete")
ace_cox_ard1 <- ace(trait.vector.cox1, trpy_n_cox1, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_cox_er1 <- corHMM(phy = trpy_n_cox1, data = cox.trait.data1[trpy_n_cox1$tip.label, c("species", "diel1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_cox_sym1 <- corHMM(phy = trpy_n_cox1, data = cox.trait.data1[trpy_n_cox1$tip.label, c("species", "diel1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_cox_ard1 <- corHMM(phy = trpy_n_cox1, data = cox.trait.data1[trpy_n_cox1$tip.label, c("species", "diel1")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 1.1b Add to table and plot ----------------------------------

#create a df comparing all models
artio_likelihoods <- data.frame(model = c("ace_cox_er1", "ace_cox_sym1", "ace_cox_ard1"), likelihood = c(ace_cox_er1$loglik, ace_cox_sym1$loglik, ace_cox_ard1$loglik), description = c("diurnal, nocturnal, artiodactyla, Cox dataset", "diurnal, nocturnal, artiodactyla, Cox dataset", "diurnal, nocturnal, artiodactyla, Cox dataset"))
View(artio_likelihoods)

#add cor results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_cox_er1", cor_cox_er1$loglik, "diurnal, nocturnal, artiodactyla, Cox dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_cox_sym1", cor_cox_sym1$loglik, "diurnal, nocturnal, artiodactyla, Cox dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_cox_ard1", cor_cox_ard1$loglik, "diurnal, nocturnal, artiodactyla, Cox dataset"))

# Section 1.1c Plot ancestral recon and transition rates ------------------
## plot the transition rates for di/noc artiodactyla
png("C:/Users/ameli/OneDrive/Documents/R_projects/cox_artio_model_dinoc.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_artio_ard1)
dev.off()

## Extract the data to plot
## Depending on which package you use, the ancestral data is in different places (or not determined at all)
## It typically only includes the internal node states, and not also the tip states
str(ace_cox_er1, max.level = 1)
head(ace_cox_er1$lik.anc)

str(cor_cox_er1, max.level = 1)
head(cor_cox_er1$states)

# In this case, I am putting together the internal nodes and tip states into one data frame
lik.anc.cox <- as.data.frame(rbind(cor_cox_ard1$tip.states, cor_cox_ard1$states))
# dim of this should be equal to the tips and internal nodes (190 tips and 189 nodes for 379 total)
trpy_n_cox1
dim(lik.anc.cox)

colnames(lik.anc.cox) <- c("diurnal", "nocturnal")

# Tips are also nodes, and all nodes have a number, and they are number sequentially beginning with tips, so the first internal node is the number of tips +1

lik.anc.cox$node <- c(1:length(trpy_n_cox1$tip.label), (length(trpy_n_cox1$tip.label) + 1):(trpy_n_cox1$Nnode + length(trpy_n_cox1$tip.label)))

ancestral_plot <- ggtree(trpy_n_cox1, layout = "circular") %<+% lik.anc.cox + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
png("C:/Users/ameli/OneDrive/Documents/R_projects/ancestral_recon_ard_cox1.png", width=17,height=16,units="cm",res=1200)
ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
dev.off()

# Section 2.0 Format the Maor et al artiodactyla data and combine with cetacean data ---------------------------------


#repeat the same analysis but from Maor et al dataset
#includes cathemeral as a state (missing from Cox et al)
maor_mam_data <- read_excel("Maor_diel_activity_data.xlsx")
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 1:4]
colnames(maor_mam_data) <- c("Order", "Family", "Species", "Activity_pattern_1")

#subset just for artiodactyla
maor_mam_data <- subset(maor_mam_data, maor_mam_data$Order == "Artiodactyla")

#if a species has an alternative diel pattern its added in a new row :[
# find how many duplicated rows there are
dim(maor_mam_data)
#there are 219 species entries
unique(maor_mam_data)
#216 are unique, so only 3 repeats?
duplicated(maor_mam_data$Species)
maor_mam_data$Species[duplicated(maor_mam_data$Species)]
duplicates1 <- maor_mam_data[duplicated(maor_mam_data$Species),]
maor_mam_data_final <-maor_mam_data[!duplicated(maor_mam_data$Species),]
#make another dataframe since some sps are repeated twice
duplicates2 <- duplicates1[duplicated(duplicates1$Species),]
duplicates1 <- duplicates1[!duplicated(duplicates1$Species),]
#activity pattern 1 = di/noc strictly
#act patt 2 = di, noc, di/crep, noc/crep, cath
#act patt 3 = di, noc, maximize for cath, maximize for crep??

#add all the extra diel patterns, then sort the columns after since they're in a random order anyway
newdf <- merge(maor_mam_data_final, duplicates1, by='Species', all.x = TRUE, all.y = TRUE)
maor_full <- merge(newdf, duplicates2, by='Species', all.x = TRUE, all.y = TRUE)
maor_full <- maor_full[, c("Species", "Activity_pattern_1", "Activity_pattern_1.x", "Activity_pattern_1.y")]
View(maor_full)

maor_full <- relocate(maor_full, "Activity_pattern_1.x", .after = "Activity_pattern_1.y")
colnames(maor_full) <- c("Species", "alt_pattern_1", "alt_pattern_2", "Activity_pattern_1")

#create column for activity pattern 2
#diel pattern 2 includes all variation, di, noc, di/crep, noc/crep, cath as reported
maor_full$Activity_pattern_2 <- maor_full$Activity_pattern_1

#create column for activity pattern 1
# diel pattern 1 maximizes for only diurnal and nocturnal
maor_full$Activity_pattern_2 <- maor_full$Activity_pattern_1
for(i in 1:nrow(maor_full))
  if(maor_full[i, "Activity_pattern_1"] == "Diurnal/Crepuscular"){
    maor_full[i, "Activity_pattern_1"] <- "Diurnal"
  } else if(maor_full[i, "Activity_pattern_1"] == "Nocturnal/Crepuscular"){
    maor_full[i, "Activity_pattern_1"] <- "Nocturnal"
  } else if(maor_full[i, "Activity_pattern_1"] %in% c("Crepuscular", "Cathemeral")){
    maor_full[i, "Activity_pattern_1"] <- NA
  }

#create activity pattern 3 column
# diel pattern 2 maximizes for cathemeral and crepuscular
maor_full$Activity_pattern_3 <- maor_full$Activity_pattern_2
for(i in 1:nrow(maor_full))
  if(maor_full[i, "Activity_pattern_2"] == "Diurnal/Crepuscular"){
    maor_full[i, "Activity_pattern_3"] <- "Crepuscular"
  } else if(maor_full[i, "Activity_pattern_2"] == "Nocturnal/Crepuscular"){
    maor_full[i, "Activity_pattern_3"] <- "Crepuscular"
  } 

#now add in the alternative diel patterns that fit with the given data :/
#First Alcelaphus lichtensteinii
maor_full[4,"Activity_pattern_2"] <- "Diurnal/Crepuscular"
maor_full[4,"Activity_pattern_3"] <- "Crepuscular"
maor_full[4,"alt_pattern_2"] <- NA

#Alces alces

#Ammotragus lervia
maor_full[6,"alt_pattern_1"] <- NA

#Axis porcinus
maor_full[13,"alt_pattern_2"] <- NA

#Bubalus bubalis
maor_full[23,"alt_pattern_2"] <- NA

#Cervus elaphus. Row 52 duplicated because there are four entries for it
#no new information in row 52 so drop it and renumber
maor_full <- maor_full[-c(52), ]
row.names(maor_full) <- 1:nrow(maor_full)

#Gazella subgutturosa
maor_full[63,"Activity_pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[63,"Activity_pattern_1"] <- "Nocturnal"
maor_full[63,"alt_pattern_2"] <- NA

#Giraffa camelopardalis
maor_full[64,"alt_pattern_2"] <- NA

#Moschus moschiferus
maor_full[92,"Activity_pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[92,"Activity_pattern_1"] <- "Nocturnal"
maor_full[92,"alt_pattern_2"] <- NA

#Philantomba maxwellii
maor_full[123,"Activity_pattern_2"] <- "Diurnal/Crepuscular"
maor_full[123,"Activity_pattern_1"] <- "Diurnal"
maor_full[123,"alt_pattern_2"] <- NA

#Rangifer tarandus
maor_full[131,"alt_pattern_2"] <- NA

#Sus barbatus
maor_full[145,"Activity_pattern_2"] <- "Nocturnal/Crepuscular"
maor_full[145,"Activity_pattern_3"] <- "Crepuscular"
maor_full[145,"alt_pattern_2"] <- NA

#Tragelaphus angasii
maor_full[156,"Activity_pattern_2"] <- "Diurnal/Crepuscular"
maor_full[156,"Activity_pattern_3"] <- "Crepuscular"
maor_full[156,"alt_pattern_2"] <- NA

#now to decide on the species with conflicting data
maor_sources <- read_excel("Maor_diel_activity_data.xlsx")
maor_sources <- maor_sources[17:3511, 1:5]
colnames(maor_sources) <- c("Order", "Family", "Species", "Activity_pattern_1", "Source")
#subset just for artiodactyla
maor_sources <- subset(maor_sources, maor_sources$Order == "Artiodactyla")

#first look at alt pattern 1
alt_df <- data.frame(maor_full[, c("Species", "alt_pattern_1", "alt_pattern_2")])
alt_df <- alt_df[!(is.na(alt_df$alt_pattern_2)),]
row.names(alt_df) <- 1:nrow(alt_df)

#gives just the sources for duplicates
alt_sources <- merge(maor_sources, alt_df, by = "Species")
alt_sources <- alt_sources[, c(1, 4, 5, 6, 7)]

unique(alt_sources$Source)
#only takes from four sources, 57, 4, 33 and 29
#Source 57: 57.Jones, K.E. et al., 2009. PanTHERIA: A species-level database of life history, ecology and geography of extant and recently extinct mammals. Ecology, 90, p.2648.
#Source 4: Aulagnier, S. & Thevenot, M., 1986. Catalogue des Mammiferes Sauvages du Maroc, Rabat, Morocco: Institute Scientifique.
#Source 33: Hufnagl, E., 1972. Lybian Mammals, Harrow, England: The Oleander Press.
# Source 29: 29. Gray, G.G. & Simpson, C.D., 1980. Ammotragus lervia. Mammalian Species, 144, pp.1–7.

# Sources 4, 29, 33 are all for  Ammotragus lervia and from the 1970s/1980s. 
# https://doi.org/10.25225/jvb.20055 more recent study from 2020 utilizing 24h camera traps found that A lervia is crepuscular and diurnal
#edit maor_full to reflect this
maor_full[6,"Activity_pattern_2"] <- "Diurnal/Crepuscular"
maor_full[6,"Activity_pattern_1"] <- "Diurnal"
maor_full[6, "Activity_pattern_3"] <- "Crepuscular"
maor_full[6,"alt_pattern_2"] <- NA

# All other 34 species have conflicting data from panTHERIA database (source 57)


#rerun models with this updated data

#create a dataframe with both maor data and cetacean data
maor_final <- maor_full[, -c(2,3)]
maor_final$Activity_pattern_1 <- tolower(maor_final$Activity_pattern_1)
maor_final$Activity_pattern_2 <- tolower(maor_final$Activity_pattern_2)
maor_final$Activity_pattern_3 <- tolower(maor_final$Activity_pattern_3)
maor_final$Species <- str_replace(maor_final$Species, " ", "_")
cetaceans_shorter2 <- cetaceans_shorter[, -c(2,3,4)]
colnames(cetaceans_shorter2) <- c("Species", "Activity_pattern_1", "Activity_pattern_2", "Activity_pattern_3")
maor_final <- rbind(maor_final, cetaceans_shorter2)
row.names(maor_final) <- c(1:nrow(maor_final))


# Section 2.1a Binary, Maor, Artiodactyla with cetacea model ----------------------------------------

#load in the tree
mam.tree <- readRDS("maxCladeCred_mammal_tree.rds")

#update df to drop the artiodactyla species missing from the mam tree
#start with 245 species, 223 in the tree
maor_final <- maor_final[maor_final$Species %in% mam.tree$tip.label,]
trpy_n_maor1 <- keep.tip(mam.tree, tip = maor_final$Species)
#tree of all the artio species with data in the mammal tree (shown without trait data)
ggtree(trpy_n_maor1, layout = "fan", size = 0.5) + geom_tiplab(size = 1.5)

#First we run the model for just the diurnal and nocturnal species
maor.trait.data1 <- maor_final[!(is.na(maor_final$Activity_pattern_1)),]
# selects only data that is in the mammal tree
maor.trait.data1 <- maor.trait.data1[maor.trait.data1$Species %in% mam.tree$tip.label,]
row.names(maor.trait.data1) <- maor.trait.data1$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_maor1 <- keep.tip(mam.tree, tip = maor.trait.data1$Species)

#now that we have the trait data for these 155 species and the subtree we can model
trait.vector.maor1 <- maor.trait.data1$Activity_pattern_1
ace_maor_er1 <- ace(trait.vector.maor1, trpy_n_maor1, model = "ER", type = "discrete")
ace_maor_sym1 <- ace(trait.vector.maor1, trpy_n_maor1, model = "SYM", type = "discrete")
ace_maor_ard1 <- ace(trait.vector.maor1, trpy_n_maor1, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maor_er1 <- corHMM(phy = trpy_n_maor1, data = maor.trait.data1[trpy_n_maor1$tip.label, c("Species", "Activity_pattern_1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maor_sym1 <- corHMM(phy = trpy_n_maor1, data = maor.trait.data1[trpy_n_maor1$tip.label, c("Species", "Activity_pattern_1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maor_ard1 <- corHMM(phy = trpy_n_maor1, data = maor.trait.data1[trpy_n_maor1$tip.label, c("Species", "Activity_pattern_1")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 2.1b Add to table -------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_er1", ace_maor_er1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_sym1", ace_maor_sym1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_ard1", ace_maor_ard1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_er1", cor_maor_er1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_sym1", cor_maor_sym1$loglik, "diurnal, nocturnal, artoidactyla, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_ard1", cor_maor_ard1$loglik, "diurnal, nocturnal, artiodactyla, maor dataset"))


# Section 2.1c Plot ancestral recon and transition rates ------------------

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
#lik.anc.maor <- as.data.frame(rbind(dinoc.trait.data[,"Activity_pattern_1"], ace_dinoc_ard1$lik.anc))
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

# Section 2.2a Four-state, Maor, Artiodactyla with cetacea --------------------------------------------

#do the same (ancestral state reconstruction) but including crepuscular and cathemeral states
maor.trait.data2 <- maor_final[!(is.na(maor_final$Activity_pattern_3)),]
# selects only data that is in the mammal tree
maor.trait.data2 <- maor.trait.data2[maor.trait.data2$Species %in% mam.tree$tip.label,]
row.names(maor.trait.data2) <- maor.trait.data2$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_maor2 <- keep.tip(mam.tree, tip = maor.trait.data2$Species)

#now that we have the trait data for these 155 species and the subtree we can model
trait.vector.maor2 <- maor.trait.data2$Activity_pattern_3
ace_maor_er2 <- ace(trait.vector.maor2, trpy_n_maor2, model = "ER", type = "discrete")
ace_maor_sym2 <- ace(trait.vector.maor2, trpy_n_maor2, model = "SYM", type = "discrete")
ace_maor_ard2 <- ace(trait.vector.maor2, trpy_n_maor2, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maor_er2 <- corHMM(phy = trpy_n_maor2, data = maor.trait.data2[trpy_n_maor2$tip.label, c("Species", "Activity_pattern_3")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maor_sym2 <- corHMM(phy = trpy_n_maor2, data = maor.trait.data2[trpy_n_maor2$tip.label, c("Species", "Activity_pattern_3")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maor_ard2 <- corHMM(phy = trpy_n_maor2, data = maor.trait.data2[trpy_n_maor2$tip.label, c("Species", "Activity_pattern_3")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 2.2b Add to table ---------------------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_er2", ace_maor_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_sym2", ace_maor_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_ard2", ace_maor_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla,maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_er2", cor_maor_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_sym2", cor_maor_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artoidactyla, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_ard2", cor_maor_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla, maor dataset"))


# Section 2.2c Plot ancestral recon and transition rates ------------------

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


# Section 2.3a Three state, Maor, Artiodactyla with cetacea --------------------------------------

#test models of artiodactyla with three states (cathemeral, diurnal, nocturnal)
maor.trait.data3 <- maor_final[!(is.na(maor_final$Activity_pattern_2)),]
maor.trait.data3$Activity_pattern_2 <- str_replace_all(maor.trait.data3$Activity_pattern_2, "nocturnal/crepuscular", "nocturnal")
maor.trait.data3$Activity_pattern_2 <- str_replace_all(maor.trait.data3$Activity_pattern_2, "diurnal/crepuscular", "diurnal")
#unlike cetaceans some artio are purely crepuscular, so switch into cathemeral
maor.trait.data3$Activity_pattern_2 <- str_replace_all(maor.trait.data3$Activity_pattern_2, "crepuscular", "cathemeral")

maor.trait.data3 <- maor.trait.data3[maor.trait.data3$Species %in% mam.tree$tip.label,]
row.names(maor.trait.data3) <- maor.trait.data3$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_maor3 <- keep.tip(mam.tree, tip = maor.trait.data3$Species)

maor.trait.vector3 <- maor.trait.data3$Activity_pattern_2
ace_maor_er3 <- ace(maor.trait.vector3, trpy_n_maor3, model = "ER", type = "discrete")
ace_maor_sym3 <- ace(maor.trait.vector3, trpy_n_maor3, model = "SYM", type = "discrete")
ace_maor_ard3 <- ace(maor.trait.vector3, trpy_n_maor3, model = "ARD", type = "discrete")

#repeat with corHMM package
cor_maor_er3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("Species", "Activity_pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maor_sym3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("Species", "Activity_pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maor_ard3 <- corHMM(phy = trpy_n_maor3, data = maor.trait.data3[trpy_n_maor3$tip.label, c("Species", "Activity_pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 2.3b Add to table -----------------------------------------

#add how well these models did to the table
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_er3", ace_maor_er3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_sym3", ace_maor_sym3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maor_ard3", ace_maor_ard3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))

artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_er3", cor_maor_er3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_sym3", cor_maor_sym3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maor_ard3", cor_maor_ard3$loglik, "diurnal, nocturnal, cathemeral, maor dataset"))


# Section 2.3c Plot ancestral recon and transition rates ------------------

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


# Section 3.0 Create dataframe of Artio without cetacea -------------------------------------

#repeat but with artiodactyla excluding cetaceans
#easiest way to do this is to drop the last chunk of rows in maor_final (to keep all the formatting)
just_artio <- maor_final[1:151, ]

just_artio <- just_artio[just_artio$Species %in% mam.tree$tip.label,]
trpy_n_ja <- keep.tip(mam.tree, tip = just_artio$Species)
#tree of all the artio species with data in the mammal tree (shown without trait data)
ggtree(trpy_n_ja, layout = "fan", size = 0.5) + geom_tiplab(size = 1.5)


# Section 3.1a Four-state, Maor, Artio without cetacea ------------------------

#ancestral state reconstruction for 4 states, diurnal, nocturnal, crepuscular and cathemeral states
trait.data.ja2 <- just_artio[!(is.na(just_artio$Activity_pattern_3)),]
# selects only data that is in the mammal tree 
trait.data.ja2 <- trait.data.ja2[trait.data.ja2$Species %in% mam.tree$tip.label,]
row.names(trait.data.ja2) <- trait.data.ja2$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_ja <- keep.tip(mam.tree, tip = trait.data.ja2$Species)

#now that we have the trait data for these 151 species and the subtree we can model
trait.vector.ja2 <- trait.data.ja2$Activity_pattern_3
ace_maorja_er2 <- ace(trait.vector.ja2, trpy_n_ja2, model = "ER", type = "discrete")
ace_maorja_sym2 <- ace(trait.vector.ja2, trpy_n_ja2, model = "SYM", type = "discrete")
ace_maorja_ard2 <- ace(trait.vector.ja2, trpy_n_ja2, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maorja_er2 <- corHMM(phy = trpy_n_ja2, data = trait.data.ja2[trpy_n_ja2$tip.label, c("Species", "Activity_pattern_3")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maorja_sym2 <- corHMM(phy = trpy_n_ja2, data = trait.data.ja2[trpy_n_ja2$tip.label, c("Species", "Activity_pattern_3")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maorja_ard2 <- corHMM(phy = trpy_n_ja2, data = trait.data.ja2[trpy_n_ja2$tip.label, c("Species", "Activity_pattern_3")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 3.1b Add to table -------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_er2", ace_maorja_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_sym2", ace_maorja_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_ard2", ace_maorja_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_er2", cor_maorja_er2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_sym2", cor_maorja_sym2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artoidactyla without cetacea, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_ard2", cor_maorja_ard2$loglik, "diurnal, nocturnal, cathemeral, crepuscular, artiodactyla without cetacea, maor dataset"))


# Section 3.1 c Plot ancestral recon and transition rates -----------------

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


# Section 3.2a Binary, Maor, Artio without cetacea ----------------
#repeat but with di/noc
trait.data.ja1 <- just_artio[!(is.na(just_artio$Activity_pattern_1)),]
# selects only data that is in the mammal tree 
trait.data.ja1 <- trait.data.ja1[trait.data.ja1$Species %in% mam.tree$tip.label,]
row.names(trait.data.ja1) <- trait.data.ja1$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_ja1 <- keep.tip(mam.tree, tip = trait.data.ja1$Species)

trait.vector.ja1 <- trait.data.ja1$Activity_pattern_1
ace_maorja_er1 <- ace(trait.vector.ja1, trpy_n_ja1, model = "ER", type = "discrete")
ace_maorja_sym1 <- ace(trait.vector.ja1, trpy_n_ja1, model = "SYM", type = "discrete")
ace_maorja_ard1 <- ace(trait.vector.ja1, trpy_n_ja1, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maorja_er1 <- corHMM(phy = trpy_n_ja1, data = trait.data.ja1[trpy_n_ja1$tip.label, c("Species", "Activity_pattern_1")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maorja_sym1 <- corHMM(phy = trpy_n_ja1, data = trait.data.ja1[trpy_n_ja1$tip.label, c("Species", "Activity_pattern_1")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maorja_ard1 <- corHMM(phy = trpy_n_ja1, data = trait.data.ja1[trpy_n_ja1$tip.label, c("Species", "Activity_pattern_1")], rate.cat = 1, model = "ARD", node.states = "marginal")


# Section 3.2b Add to table --------------------------------------

#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_er1", ace_maorja_er1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_sym1", ace_maorja_sym1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_ard1", ace_maorja_ard1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_er1", cor_maorja_er1$loglik, "diurnal, nocturnal, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_sym1", cor_maorja_sym1$loglik, "diurnal, nocturnal,  artoidactyla without cetacea, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_ard1", cor_maorja_ard1$loglik, "diurnal, nocturnal,  artiodactyla without cetacea, maor dataset"))


# Section 3.2c Plot ancestral recon and transition rates ------------------

##plot the transition rates
png("C:/Users/ameli/OneDrive/Documents/R_projects/just_artio_ard_dinoc.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maorja_ard1)
dev.off()

#plot the ancestral reconstruction
#need to do


# Section 3.3a Three-state, Maor, Artio without cetacea -------------------

#artio without cetacea three states models
trait.data.ja3 <- just_artio[!(is.na(just_artio$Activity_pattern_2)),]
trait.data.ja3$Activity_pattern_2 <- str_replace_all(trait.data.ja3$Activity_pattern_2, "nocturnal/crepuscular", "nocturnal")
trait.data.ja3$Activity_pattern_2 <- str_replace_all(trait.data.ja3$Activity_pattern_2, "diurnal/crepuscular", "diurnal")
#unlike cetaceans some artio are purely crepuscular, so switch into cathemeral
trait.data.ja3$Activity_pattern_2 <- str_replace_all(trait.data.ja3$Activity_pattern_2, "crepuscular", "cathemeral")

trait.data.ja3 <- trait.data.ja3[trait.data.ja3$Species %in% mam.tree$tip.label,]
row.names(trait.data.ja3) <- trait.data.ja3$Species
# this selects a tree that is only the subset with data (mutual exclusive)
trpy_n_ja3 <- keep.tip(mam.tree, tip = trait.data.ja3$Species)

trait.vector.ja3 <- trait.data.ja3$Activity_pattern_2
ace_maorja_er3 <- ace(trait.vector.ja3, trpy_n_ja3, model = "ER", type = "discrete")
ace_maorja_sym3 <- ace(trait.vector.ja3, trpy_n_ja3, model = "SYM", type = "discrete")
ace_maorja_ard3 <- ace(trait.vector.ja3, trpy_n_ja3, model = "ARD", type = "discrete")

## e.g. running the same as the above, but with the corHMM package
cor_maorja_er3 <- corHMM(phy = trpy_n_ja3, data = trait.data.ja3[trpy_n_ja3$tip.label, c("Species", "Activity_pattern_2")], rate.cat = 1, model = "ER", node.states = "marginal")
cor_maorja_sym3 <- corHMM(phy = trpy_n_ja3, data = trait.data.ja3[trpy_n_ja3$tip.label, c("Species", "Activity_pattern_2")], rate.cat = 1, model = "SYM", node.states = "marginal")
cor_maorja_ard3 <- corHMM(phy = trpy_n_ja3, data = trait.data.ja3[trpy_n_ja3$tip.label, c("Species", "Activity_pattern_2")], rate.cat = 1, model = "ARD", node.states = "marginal")

# Section 3.3b Add to table and plot ------------------------------------------------------------
#add results to dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_er3", ace_maorja_er3$loglik, "diurnal, nocturnal, cathemeral. artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_sym3", ace_maorja_sym3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("ace_maorja_ard3", ace_maorja_ard3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))

#add results to the dataframe
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_er3", cor_maorja_er3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_sym3", cor_maorja_sym3$loglik, "diurnal, nocturnal, cathemeral, artoidactyla without cetacea, maor dataset"))  
artio_likelihoods <- rbind(artio_likelihoods, c("cor_maorja_ard3", cor_maorja_ard3$loglik, "diurnal, nocturnal, cathemeral, artiodactyla without cetacea, maor dataset"))


# Section 3.3c Plot ancestral reconstruction and transition rates ---------

#plot the transition rates
png("C:/Users/ameli/OneDrive/Documents/R_projects/just_artio_ard_3state.png", width=14,height=13,units="cm",res=1200)
plotMKmodel(cor_maorja_ard3)
dev.off()


# Export likelihood table -------------------------------------------------
#export dataframe of likelihood as a png
png("C:/Users/ameli/OneDrive/Documents/R_projects/artio_likelihood_table.png", height = 30*nrow(artio_likelihoods), width = 300*ncol(artio_likelihoods), res = 90)
grid.table(artio_likelihoods)
dev.off()