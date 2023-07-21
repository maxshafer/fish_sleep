library(ape)
library(corHMM)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("must have 2 arguments for this Rscript command, 1) the name_variable,  2) the dataset_variable")
}

if (!(args[[1]] %in% c("all", "only_highqual", "only_ingroup", "only_cartilaginous", "not_mammals", "amphibians", "sauropsids", "lepidosauria", "testudines", "crocodylia", "aves"))) {
  stop("wrong name_variable")
}

if (!(args[[2]] %in% c("fish", "AllGroups", "mammals", "tetrapods"))) {
  stop("wrong dataset_variable")
}

# setwd("/Users/maxwellshafer/Documents/R_Projects/fish_sleep")
setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")

source("scripts/fish_sleep_functions.R")

## Load files
# Which tree?
name_variable <- args[[1]] # all, only_highqual, only_cartilaginous, or only_ingroup, or not_mammals
dataset_variable <- args[[2]] # fish or AllGroups, or mammals, or tetrapods

trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable)
trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)


###############################################################################################################################################
### Run Hidden Rates Models ### 
###############################################################################################################################################

standard_tests <- list()

## First run them with marginal reconstruction, then joint
## 1, 2, or 3 rate categories

standard_tests[[1]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = "ARD", node.states = "marginal")
standard_tests[[2]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = "ARD", get.tip.states = TRUE, node.states = "marginal")
standard_tests[[3]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 3, model = "ARD", get.tip.states = TRUE, node.states = "marginal")


standard_tests[[4]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = "ARD", node.states = "joint")
standard_tests[[5]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = "ARD", get.tip.states = TRUE, node.states = "joint")
standard_tests[[6]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 3, model = "ARD", get.tip.states = TRUE, node.states = "joint")

# Name the models in the list and save
names(standard_tests) <- c("MK_2state_marg", "HMM_2state_2rate_marg", "HMM_2state_3rate_marg", "MK_2state_joint", "HMM_2state_2rate_joint", "HMM_2state_3rate_joint")
saveRDS(standard_tests, file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))




