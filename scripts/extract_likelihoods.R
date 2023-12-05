library(phangorn)
library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(zoo)
library(ggtree)
library(here)

setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

index_list <- list()
index_list[[1]] <- c("all", "only_highqual", "only_cartilaginous", "only_ingroup")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals", "amphibians", "sauropsids", "lepidosauria", "testudines", "aves")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

# signal_K <- index_list
# signal_lambda <- index_list
log.likelihood <- list()
AIC <- list()

recon <- "marg"

for (i in 1:length(index_list)) {
  dataset_variable <- names(index_list)[[i]]
  log.likelihood[[i]] <- list()
  AIC[[i]] <- list()
  
  for (j in 1:length(index_list[[i]])) {
    name_variable <- index_list[[i]][[j]]
    
    ## Load in the tree
    trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
    trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)
    
    ## Load the best model, which is the HMM 2 state 2 rate model
    models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
    
    plotMKmodel(models$HMM_2state_2rate_marg)
    models$HMM_2state_2rate_marg
    log.likelihood[[i]][[j]] <- unlist(lapply(models, function(x) x$loglik))[c(2:3,7)]
    AIC[[i]][[j]] <- unlist(lapply(models, function(x) x$AIC))[c(2:3,7)]
    
    
    # signal_K[[i]][[j]] <- phylosig(trpy_n, as.numeric(as.factor(trait.data_n$diel1)), method="K", test=FALSE, nsim=1000, se=NULL, start=NULL, control=list())
    # signal_lambda[[i]][[j]] <- phylosig(trpy_n, as.numeric(as.factor(trait.data_n$diel1)), method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL, control=list())[[1]]
    
    # data <- as.numeric(as.factor(trait.data_n$diel1))
    # names(data) <- trpy_n$tip.label
    # test <- fitDiscrete(phy = trpy_n, dat = data, model = "ER", transform = "lambda")
  }
}
    
    
    
    