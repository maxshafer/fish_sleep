#first custom function SubsetTrees takes the 1k trees and the subset of mammals (cetaceans, artio)
#returns the 1k subsetted trees

subsetTrees <- function(tree = mammal_trees[[1]], subset_names = cetaceans_full$tips) {
  # will take a tree and keep only those tips that match subset_names
  
  subset_names <- subset_names[subset_names %in% tree$tip.label]
  
  out_tree <- keep.tip(tree, subset_names)
  
  return(out_tree)
  
}

#second custom function returnModels that takes the subsetted tree, the input data (diel patterns +species names) and the model type you want to run (ER, SYM, ARD)
#returns the model results for each of the trees provided. 

returnAceModels <- function(tree = cetacean_trees[[1]], trait.data = trait_data, column = "Diel_Pattern_1", model = "SYM") {
  trait.data <- trait.data[trait.data$tips %in% tree$tip.label,]
  trait.data <- trait.data[!(is.na(trait.data[,column])),]
  tree <- keep.tip(tree, trait.data$tips)
  ace_model <- ace(trait.data[,column], tree, model = model, type = "discrete")
  return(ace_model)
  
}

#third custom function returnAceModels that takes the subsetted tree, the input data (diel patterns +species names) and the model type you want to run (ER, SYM, ARD)
#same as before but for the cor model
#added a conditional so it can take custom rate matrices

returnCorModels <- function(tree = cetacean_trees[[2]], trait.data = cetaceans_full, diel_col = "Diel_Pattern_1", rate.cat = 1, custom.rate.mat = "none", model = "SYM", node.states = "marginal"){
  trait.data <- trait.data[trait.data$tips %in% tree$tip.label,]
  row.names(trait.data) <- trait.data$tips
  trait.data <- trait.data[tree$tip.label, c("tips", diel_col)]
  trait.data <- trait.data[!(is.na(trait.data[, diel_col])),]
  tree <- keep.tip(tree, trait.data$tips)
  
  if(!("none" %in% custom.rate.mat)) {
    cor_model <- corHMM(phy = tree, data = trait.data, rate.cat = rate.cat, rate.mat = custom.rate.mat, model = model, node.states = node.states)
  } else {
    cor_model <- corHMM(phy = tree, data = trait.data, rate.cat = rate.cat, model = model, node.states = node.states)
  }
  
  return(cor_model)
}

#testing to see if this works with and without custom rate matrices
# returnCorModels(phylo_trees[[2]], trait.data, "Diel_Pattern_1", rate.cat = 1, custom.rate.mat = matrix(c(0,0,2,0), nrow = 2, ncol = 2, dimnames = list(c("(1)", "(2)"), c("(1)", "(2)"))), model = "SYM", node.states = "marginal")
# returnCorModels(phylo_trees[[2]], trait.data, "Diel_Pattern_1", rate.cat = 1, custom.rate.mat = "none", model = "SYM", node.states = "marginal")
# returnCorModels(phylo_trees[[2]], trait.data, "Diel_Pattern_1", rate.cat = 1, model = "SYM", node.states = "marginal")

#fourth custom function returnLikelihoods takes model results from returnModels and returns just the likelihoods of those models 

returnLikelihoods <- function(model = cetacean_sim_ace[[1]], return = "loglik"){
  
  return(model[return])
}

#fifth custom  function
#function to extract the transition rates from all models, similar to likelihoods

returnRates <- function(model = cetacean_sim_ace[[1]], return = "solution"){
  
  return(model[return])
}

#sixth function
corhmm_model <- readRDS("~/R_projects/fish_sleep/cetaceans_max_clade_cred_six_state_traits_ER_SYM_ARD_models.rds")
corhmm_model <- corhmm_model$ARD_model
returnAICc <- function(model = corhmm_model, return = "AICc"){
  
  return(model[return])
}

#seventh function
corhmm_model <- readRDS("~/R_projects/fish_sleep/cetaceans_max_clade_cred_six_state_traits_ER_SYM_ARD_models.rds")
corhmm_model <- corhmm_model$ARD_model
returnAIC <- function(model = corhmm_model, return = "AIC"){
  
  return(model[return])
}

#eight function
#create a function that returns the node label for the mrca of a set of species/taxonomic group

#requires a dataframe with columns for species names and their taxonomic level names
#returns a dataframe with the clade name and the mrca node number

findMRCANode <- function(phylo = tr, trait.data = trait.data, taxonomic_level_col = 3){
  nodes_list <- list()
  for(x in 1:nrow(unique(trait.data[,taxonomic_level_col]))){
    i <- unique(trait.data[x,taxonomic_level_col])
    #ensure the species are in the tree you're working with
    trait.data <- trait.data[trait.data$tips %in% phylo$tip.label,]
    #remove any taxonomic levels with only one species (cannot find MRCA for one species)
    trait.data <- trait.data %>% group_by_at(taxonomic_level_col) %>% filter(n()>1)
    #subset the trait data into species belonging to the same taxonomic group
    trait.data.filtered <- trait.data[trait.data[,taxonomic_level_col] == i, ]
    #take the vector of these species names
    tip_vector <- as.vector(trait.data.filtered[, "tips"])
    #find the node number of the MRCA of these species
    MRCA_node <- findMRCA(phylo, tip_vector$tips)
    #create a dataframe
    loop_df <- data.frame(clade_name = i, node_number = MRCA_node)
    #add to list, to extract later (out of for loop)
    nodes_list[[i]] <- loop_df}
  
  nodes_df <- do.call(rbind, nodes_list)
  return(nodes_df)
}
  

findMRCANode2 <- function(phylo = tr, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = "Araneae"){
  nodes_list <- list()
  trait.data <- trait.data %>% group_by_at(taxonomic_level_col) %>% filter(n()>1)
  trait.data.filtered <- trait.data[trait.data[,taxonomic_level_col] == taxonomic_level_name, ]
  tip_vector <- as.vector(trait.data.filtered[, "tips"])
  #find the node number of the MRCA of these species
  MRCA_node <- findMRCA(phylo, tip_vector$tips)
  node_df <- data.frame(clade_name = taxonomic_level_name, node_number = MRCA_node)
  return(node_df)
}

#delta statistic function
#a function that given a df with species names and diel activity data returns the delta statistic for that group

calculateDelta <- function(trait.data = trait.data, taxonomic_group_name = "Cetacea"){
  mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
  trait.data <- filter(trait.data, Order == taxonomic_group_name)
  trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
  mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)
  mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(tree$edge.length, 0.1)*0.1
  
  sps_order <- as.data.frame(mam.tree$tip.label)
  colnames(sps_order) <- "tips"
  sps_order$id <- 1:nrow(sps_order)
  trait.data <- merge(trait.data, sps_order, by = "tips")
  trait.data <- trait.data[order(trait.data$id), ]
  trait <- trait.data$Diel_Pattern
  
  source("scripts/Amelia_delta_code.R")
  delta_diel <- delta(trait, mam.tree, 0.1, 0.0589, 1000, 10, 100)
  
  random_delta <- rep(NA,100)
  for (i in 1:100){
    rtrait <- sample(trait)
    random_delta[i] <- delta(rtrait,mam.tree,0.1,0.0589,10000,10,100)
  }
  p_value <- sum(random_delta>delta_diel)/length(random_delta)
  
  result <- data.frame(delta_diel, taxonomic_group_name, p_value)
  return(result)
}
