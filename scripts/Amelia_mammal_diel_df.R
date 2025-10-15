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
setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Section 1: How well do the Bennie and Maor dfs agree with each --------

#read in the Bennie diel activity patterns
#from https://doi.org/10.1073/pnas.1216063110 
Bennie_mam_data <- read_excel(here("Bennie_diel_activity_data.xlsx"))
colnames(Bennie_mam_data) <- "SpeciesBehaviourReference"
Bennie_mam_data$SpeciesBehaviourReference <- str_replace(string = Bennie_mam_data$SpeciesBehaviourReference, pattern = " ", replacement  = "_")
Bennie_mam_data <- separate(Bennie_mam_data, col = SpeciesBehaviourReference, into = c("Species", "Activity_pattern", "Reference"), sep = " ")

#check how many artio species are in bennie dataset -224 species
artio <- read.csv(here("sleepy_artiodactyla_minus_cetaceans.csv"))
nrow(Bennie_mam_data[Bennie_mam_data$Species %in% artio$tips, ])

#read in the Maor diel activity patterns
#from https://doi.org/10.1038/s41559-017-0366-5 
maor_mam_data <- read_excel(here("Maor_diel_activity_data.xlsx"))
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 3:5]
colnames(maor_mam_data) <- c("Species", "Activity_pattern", "Reference")

#the issue with this df is that if they have an alternative activity they add it in a new row
duplicates1 <- maor_mam_data[duplicated(maor_mam_data$Species),]
#make another dataframe since some sps are repeated twice
duplicates2 <- duplicates1[duplicated(duplicates1$Species),] #1080 species have at least one alt diel pattern
duplicates1 <- duplicates1[!duplicated(duplicates1$Species),] #208 species have 2 alt diel patterns

#remove duplicates from Maor dataframe for now, we'll add alternative diel patterns back next
maor_mam_data <-maor_mam_data[!duplicated(maor_mam_data$Species),]
maor_mam_data$Species <- str_replace(string = maor_mam_data$Species, pattern = " ", replacement  = "_")

#format some entries to be max_crep (di/crep -> crep, noc/crep -> crep)
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal/ Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal/Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal /Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal or Cathemeral", "Cathemeral")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal/Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal - EXTINCT", "Nocturnal")

#Don't add back alternative patterns for now
# #add all the extra diel patterns, then sort the columns after since they're in a random order anyway
# maor_mam_data <- merge(maor_mam_data, duplicates1, by='Species', all.x = TRUE, all.y = TRUE)
# maor_full <- merge(maor_mam_data, duplicates2, by='Species', all.x = TRUE, all.y = TRUE)
# maor_full <- maor_full[, c("Species", "Activity_pattern", "Activity_pattern.x", "Activity_pattern.y")]
# maor_full <- relocate(maor_full, "Activity_pattern.x", .after = "Activity_pattern.y")
# colnames(maor_full) <- c("Species", "alt_pattern_1", "alt_pattern_2", "Activity_pattern")


# Section 2: How well do these sources agree? -----------------------------
diel_merge <- merge(Bennie_mam_data,maor_mam_data,by="Species")
colnames(diel_merge) <- c("Species", "Bennie_diel", "Bennie_source", "Maor_diel", "Maor_source")
diel_merge$Bennie_diel <- tolower(diel_merge$Bennie_diel)
diel_merge$Maor_diel <- tolower(diel_merge$Maor_diel)
#set default to idk
diel_merge$match <- "Idk"

for(i in 1:length(diel_merge$Species)){
  if(diel_merge[i, "Bennie_diel"] == diel_merge[i, "Maor_diel"]){
    diel_merge[i, "match"] <- "Yes"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal") & diel_merge[i, "Maor_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")& diel_merge[i, "Maor_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("cathemeral", "cathemeral/crepuscular")& diel_merge[i, "Maor_diel"] %in% c("cathemeral", "cathemeral/crepuscular")){
    diel_merge[i, "match"] <- "Approximate"
  } else {
    diel_merge[i, "match"] <- "No"
  }
}

table(diel_merge$match)

#A surprising amount of consistency!
#170/2199 approximate matches, 244/2199 no matches, and 1785/2199 yes matches
#81% match!
diel_table <- diel_merge %>% count(match)
diel_table <- transform(diel_table, percent = (n/sum(diel_table$n)) * 100)

diel_table <- trait.data.all %>% count(Confidence)
diel_table <- transform(diel_table, percent = (n/sum(diel_table$n)) * 100)

ggplot(diel_table, aes(x="", y=n, fill=Confidence)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + theme_void()

#which species did not match
mismatch_species <- diel_merge %>% filter(match == "No")
table(mismatch_species$Bennie_diel, mismatch_species$Maor_diel)
#we should clean some of these entries up, could be a formatting error in why they don't match

#plot this on the mammal tree, which species tend to have conflicting data?
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- mismatch_species[,1:2]
trait.data <- trait.data[trait.data$Species %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$Species)
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)

#seem to be distributed throughout the tree

#get taxonomic information from Maor et al
maor_mam_data <- read_excel(here("Maor_diel_activity_data.xlsx"))
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 1:3]
colnames(maor_mam_data) <- c("Order", "Family", "tips")
maor_mam_data$tips <- str_replace(maor_mam_data$tips, pattern = " ", replacement = "_")

#diel_merge currently has 2,199 species
diel_merge <- merge(diel_merge, maor_mam_data, by = "tips") #now is 3,129 species
diel_merge <- diel_merge[!duplicated(diel_merge$tips),] #remove duplicates, back to 2,199

#save out the final 
colnames(diel_merge) <- c("tips", "Bennie_diel", "Bennie_source", "Maor_diel", "Maor_source", "match", "Order", "Family")
write.csv(diel_merge, here("sleepy_mammals_old.csv"))

# Section 3: Visualizing the mammal data -------------------------------
trait.data <- read.csv(here("sleepy_mammals.csv"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

trait.data <- filter(trait.data, match == "Yes") 
trait.data <- trait.data[, c("tips", "Bennie_diel", "Order", "Family")]
#mammals_df <- mammals_df[, c("tips", "Maor_diel", "Order", "Family")]
colnames(trait.data) <- c("tips", "Diel_Pattern", "Order", "Family")
#trait.data <- relocate(trait.data, "tips", .after = "Diel_Pattern")
 
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$tips)

diel.plot.all <- ggtree(trpy_n_mam, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern")]
diel.plot.all <- diel.plot.all + geom_tile(data = diel.plot.all$data[1:length(trpy_n_mam$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, width = 5)
diel.plot.all <- diel.plot.all #+ geom_tiplab(size = 1, offset = 3)
diel.plot.all

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", "mammals_bennie_maor_agree_unlabelled", "diel_plot.png"), width=24,height=18,units="cm",res=1200)
diel.plot.all
dev.off()

#label by major orders
source("scripts/Amelia_functions.R")
trait.data <- trait.data %>% group_by(Order) %>% filter(n()>1)
all_nodes <- lapply(unique(trait.data$Order), function(x) findMRCANode2(phylo = trpy_n_mam, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = x))
nodes_df <- do.call(rbind, all_nodes)

#needs some tweaking still
#nodes_df$colour <- c("navy", "slateblue", "mediumpurple", "dodgerblue", "darkorchid1", "royalblue", "lightslateblue", "purple3", "steelblue",)
#nodes_df$colour <- c("white", "grey95", "grey90", "grey85", "grey80", "grey75", "grey70", "grey65", "grey60", "grey55", "grey50", "grey45", "grey40", "grey35", "grey30", "grey25", "grey20", "grey15", "grey10", "grey5", "black")
order_tree <- diel.plot.all + geom_cladelab(node = nodes_df$node_number, label = nodes_df$clade_name, offset = 1.8, offset.text = 1, hjust = 0, vjust = 0, align = FALSE, fontsize = 2.5, barsize = 2.5, barcolour = "grey")
order_tree

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/", "mammals_bennie_maor_agree_order_labelled", "diel_plot.png"), width=24,height=18,units="cm",res=1200)
order_tree
dev.off()