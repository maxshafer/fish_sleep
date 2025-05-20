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
library(ggplot2)
#install.packages("rphylopic")
library(rphylopic)
library(RColorBrewer)
#install.packages("ggimage")
library(ggimage)

setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: Load in fish data
url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit?gid=0#gid=0' 
#url <- 'https://docs.google.com/spreadsheets/d/14vd7fH7Ra64nNFi5mq34iEUfIIfE2MGHqLEfxJ8C71c/edit?usp=sharing' 
fish_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

fish_full <- filter(fish_full, Source == "Schier lab Cichlid screen")

trait.data <- fish_full[, c("Species_name", "X5", "Source")]
colnames(trait.data) <- c("Species_name", "Temporal_activity_pattern")

#replace species with non-standard names with proxy names so they can be placed in the tree (won't keep these species names)
#trait.data$Species_name <- str_replace(trait.data$Species_name, pattern = "Haplochromis 'red back scraper'", replacement = )
#trait.data$Species_name <- str_replace(trait.data$Species_name, pattern = "Altolamprologus_shell", replacement = "Altolamprologus_shell")
# Astatotilapia_burtoni should be in tree


trait.data$tips <- str_replace(trait.data$Species_name, pattern = " ", replacement = "_")

fish.tree <- readRDS(here("tr_tree_calibrated_fish.rds"))

#trait.data <- trait.data[!trait.data$tips %in% fish.tree$tip.label,]
trait.data <- trait.data[trait.data$tips %in% fish.tree$tip.label,]
trpy_n <- keep.tip(fish.tree, tip = trait.data$tips)
rownames(trait.data) <- trait.data$tips

#only 50 species in this tree
ggtree(trpy_n, layout = "circular") + geom_tiplab()

custom.colours <- c("pink","#fbbe30", "#A6D854", "#FC8D62", "#66C2A5", "#dd8ae7")

#plot the tree
diel.plot.all <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Temporal_activity_pattern")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Temporal_activity_pattern), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours)
diel.plot.all <- diel.plot.all + geom_tiplab(size = 4)
diel.plot.all 

saveRDS(trait.data, file = here("cichlid.trait.data.rds"))
saveRDS(trpy_n, file = here("cichlid_trpy_n.rds"))
# OTL ---------------------------------------------------------------------


#same species still missing

#make tree with open tree of life
resolved_names <- tnrs_match_names(names = trait.data$Species_name, context_name = "Vertebrates", do_approximate_matching = FALSE)

## Remove any that don't have exact matches or are synonyms (this is really just for finding ancestral state, so missing a few species won't matter)
#Remove species without open tree of life entries. Removes Balaenoptera riceii (newly discovered, not in the tree)
resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
#Removes species that are synonyms for other species. Removes Inia humboldtiana (is a synonym/subspecies for Inia geoffrensis) and Neophocaena sunameri (subspecies of Neophocaena asiaeorientalis)
resolved_names <- resolved_names[resolved_names$is_synonym == FALSE,]
#no approximate matches, so this removes no species currently
resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]

#we've removed these three species, so now we have 77 cetaceans in our resolved_names dataframe to work with

# Remove excess information, clean up, and add tip label ids that will match the tree
resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
resolved_names <- resolved_names[!duplicated(resolved_names$tips),]

# Add data on activity patterns, starting with diel pattern 3 (crep, cath, di, noc)
# match() finds the position of its first argument in its second
resolved_names$diel <- trait.data$Temporal_activity_pattern[match(resolved_names$search_string, tolower(trait.data$Species_name))]
#check that the activity patterns look correct
table(resolved_names$diel)


## Fetch the tree
#subset the open tree of life to only include the species in resolved_names (find them by their ott_id)
tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") 
# I need to use the id option here, and then use that to map the tip labels from resolved_names (that way I don't run into the issue with the difference in formatting between the two tools)

# Time calibrate it using geiger and timetree.org - the open tree of life is already time calibrated, so we don't need to do this step
# First resolve polytomies ~randomly using multi2dr

tr <- multi2di(tr)

#currently the tips are identified by their ott_id, can change to their species name (tip column)
tr$tip.label <- resolved_names$tips[match(tr$tip.label, paste("ott", resolved_names$ott_id, sep = ""))]

#this is the tree of all the cetacean species we are working with in the open tree of life
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 2)


# Four state cichlid modelling --------------------------------------------

#read in the trait data and tree
trait.data <- readRDS(here("cichlid.trait.data.rds"))
trpy_n <- readRDS(here("cichlid_trpy_n.rds"))

#

#
ER <- corHMM(phy = trpy_n, data = trait.data[, c("tips", "Temporal_activity_pattern")], rate.cat = 1, model = "ER", node.states = "marginal")
SYM <- corHMM(phy = trpy_n, data = trait.data[, c("tips", "Temporal_activity_pattern")], rate.cat = 1, model = "SYM", node.states = "marginal")
ARD <- corHMM(phy = trpy_n, data = trait.data[, c("tips", "Temporal_activity_pattern")], rate.cat = 1, model = "ARD", node.states = "marginal")

model_results <- data.frame(MK_model = c("ER", "SYM", "ARD"), AIC_score = c(ER$AIC, SYM$AIC, ARD$AIC), AICc_score = c(ER$AICc, SYM$AICc, ARD$AICc), log_lik = c(ER$loglik, SYM$loglik, ARD$loglik))

ggplot(model_results, aes(x = fct_inorder(MK_model), y = AIC_score)) + geom_bar(stat = "identity") #+ annotate("text", x = 1, y = min(model_results$AIC_score), label = "*")
ggplot(model_results, aes(x = fct_inorder(MK_model), y = AICc_score)) + geom_bar(stat = "identity")
ggplot(model_results, aes(x = fct_inorder(MK_model), y = log_lik)) + geom_bar(stat = "identity")  #+ annotate("text", x = 1, y = max(model_results$AIC_score), label = "*")



