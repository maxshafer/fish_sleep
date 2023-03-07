library(ape) 
library(corHMM)
library(xlsx)
library(phangorn)
library(stringr)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

# load data for mammals from Cox et al

mam.table <- read.xlsx("Cox_mammal_data/Supplementary Data 2.xlsx", 1)
mammal_trees <- read.nexus("Cox_mammal_phylo/Complete_phylogeny.nex")

## Determine the consensus tree
mam.tree <- maxCladeCred(mammal_trees, tree = TRUE)

## There are two sources for diel activity, not sure which is which? Can I combine them?
## I would say that the Activity_IM has more data (165 more entries), and there is only 1 where they disagree in a manner that isn't consistant with my analyses (allowed to be Crepuscular AND Diurnal/Nocturnal)
mam.table$diel <- tolower(apply(mam.table, 1, function(x) {
  sub <- x[c("Activity_IM", "Activity_DD")]
  out <- ifelse("NA" %in% sub, 
                sub[1], 
                ifelse(sub[1] == sub[2], 
                  sub[1],
                  ifelse("Crepuscular" %in% sub, 
                        paste(sub[grepl("Crepuscular", sub)], sub[!(grepl("Crepuscular", sub))], sep = "/"),
                        sub[1])))
  return(out)
}))

mam.table$diel1 <- ifelse(mam.table$diel == "diurnal", "diurnal", ifelse(mam.table$diel == "nocturnal", "nocturnal", ifelse(mam.table$diel == "crepuscular", "crepuscular", ifelse(mam.table$diel == "crepuscular/diurnal", "diurnal", ifelse(mam.table$diel == "crepuscular/nocturnal", "nocturnal", ifelse(mam.table$diel == "crepuscular/unclear", "crepuscular", "unknown"))))))
mam.table$diel2 <- ifelse(mam.table$diel == "diurnal", "diurnal", ifelse(mam.table$diel == "nocturnal", "nocturnal", ifelse(mam.table$diel == "crepuscular", "crepuscular", ifelse(mam.table$diel == "crepuscular/diurnal", "crepuscular", ifelse(mam.table$diel == "crepuscular/nocturnal", "crepuscular", ifelse(mam.table$diel == "crepuscular/unclear", "crepuscular", "unknown"))))))

## diel is raw (crep/di, crep/noc, di, noc)
## diel1 is only di or noc
## diel2 is crep, di or noc

trait.data <- mam.table[,c("Binomial_iucn", "Order", "Family", "Genus", "diel", "diel1", "diel2")]
colnames(trait.data) <- c("species", "Order", "Family", "Genus", "diel", "diel1", "diel2")
trait.data$species <- str_replace(trait.data$species, " ", "_")
trait.data <- trait.data[!(is.na(trait.data$diel)),]
trait.data <- trait.data[trait.data$diel1 %in% c("diurnal", "nocturnal"),]

# Reciprically determine matching names between trait data and phylogenetic tree
trait.data <- trait.data[trait.data$species %in% mam.tree$tip.label,]
trait.data <- trait.data[!(duplicated(trait.data)),]
row.names(trait.data) <- trait.data$species

trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$species)

### Save out trait data and tree
### This is the Maximum Clade Credibility tree, and associated Diurnal/Nocturnal trait data (mutually exclusive extant species)
dataset <- "mammals"

saveRDS(trpy_n_mam, paste("tr_tree_calibrated_", dataset, ".rds", sep = ""))

saveRDS(trait.data, paste("trait_data_", dataset, ".rds", sep = ""))





