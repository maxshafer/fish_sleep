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
index_list[[3]] <- c("amphibians", "sauropsids")
names(index_list) <- c("fish", "mammals", "tetrapods")

standard_tests <- list()

for (i in 1:length(index_list)) {
  dataset_variable <- names(index_list)[[i]]
  
  standard_tests[[i]] <- list()
  
  for (j in 1:length(index_list[[i]])) {
    name_variable <- index_list[[i]][[j]]
    
    trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable)
    trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)
    
    models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
    
    ###############################################################################################################################################
    ### Run Hidden Rates Models ### 
    ###############################################################################################################################################
    
    data <- trait.data_n[trpy_n$tip.label, c("species", "diel1")]
    
    RateCat1 <- getStateMat4Dat(data)$rate.mat # R1
    RateCat1 <- equateStateMatPars(RateCat1, c(1:4))
    
    RateCat2 <- getStateMat4Dat(data)$rate.mat # R2
    RateCat2 <- dropStateMatPars(RateCat2, 3)
    RateCat2[2] <- 0
    RateCat2[3] <- 0
    
    RateClassMat <- getRateCatMat(2) #
    
    StateMats <- list(RateCat1, RateCat2)
    
    FullMat <- getFullMat(StateMats, RateClassMat)
    
    standard_tests[[i]][[j]] <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, rate.mat = FullMat, get.tip.states = FALSE, node.states = "marginal")
    
    models[[7]] <- standard_tests[[i]][[j]]
    saveRDS(models, file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
    
  }
}