library(ape)
library(corHMM)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("must have 2 arguments for this Rscript command, 1) the name_variable,  2) the dataset_variable")
}

if (!(args[[1]] %in% c("all", "only_highqual", "only_cartilaginous", "only_ingroup", "not_mammals"))) {
  stop("wrong name_variable")
}

if (!(args[[2]] %in% c("fish", "AllGroups", "mammals", "tetrapods"))) {
  stop("wrong dataset_variable")
}


# setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")
setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/fish_sleep_functions.R")

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


## Output the simmap of the best model with 100 simulations?
## Phytools has to run the model itself, but the simulation should just need the Q matrix and tip states
## Can I hack make.simmap to use the Q matrix from the best model?
## Also, I'm not sure how to pull node values, since these are all edge lengths
library(phytools)

trait.data_v <- trait.data_n[,c("diel1")]
names(trait.data_v) <- trait.data_n$species

test <- make.simmap(trpy_n, trait.data_v, model = "ARD", nsim = 10)

# plotSimmap(test[[1]])


## Can also run it with an already computed Q matrix, somehow it takes a long time with the corHMM matrices... / never finishes?
## Seems that the rates for diurnal -> diurnal should be a slight negative????
test1 <- make.simmap(trpy_n, trait.data_v, model = "ARD", nsim = 1, Q = test[[1]]$Q)
Q <- standard_tests[[1]]$solution
rownames(Q) <- colnames(Q) <- c("diurnal", "nocturnal")
Q[is.na(Q)] <- -0.3
test2 <- make.simmap(trpy_n, trait.data_v, model = "ARD", nsim = 1, Q = Q)

plotSimmap(test2)





