library(phangorn)
library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(zoo)
library(ggtree)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)
library(dplyr)
library(tibble)


setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

run_simulations <- TRUE
run_simulations_cumsums <- TRUE
run_simmap <- TRUE
run_simmap_cumsums <- TRUE

index_list <- list()
index_list[[1]] <- c("all", "only_cartilaginous", "only_ingroup", "only_highqual")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals", "amphibians", "sauropsids", "lepidosauria", "testudines", "aves")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

# Set simulation parameters (type and #)
model_types <- c("ARD", "HR")
sim_numb <- 500

recon <- "marg"

for (y in 1:length(model_types)) {
  model_type <- model_types[[y]]
  
  for (i in 1:length(index_list)) {
    dataset_variable <- names(index_list)[[i]]
    
    for (j in 1:length(index_list[[i]])) {
      name_variable <- index_list[[i]][[j]]
      
      setwd(here())
      
      ###############################################################################################################################################
      ### Load files (trees/data/models) ### 
      ###############################################################################################################################################
      
      ## Load in the tree + data
      trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
      trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)
      
      ## Load the best model, which is the HMM 2 state 2 rate model
      models <- readRDS(file = here(paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_")))
      
      ## I think I will always take the 2 rate model, 3 is too hard to comprehend
      
      ## If I want to use the joint reconstruction then the structure of the data changes (marginal gives % for tips + nodes, joint gives states as vectors)
      if (recon == "joint") {
        model <- models$HMM_2state_2rate_joint
        model_ER <- models$MK_2state_joint
        model_ARD <- models$MK_2state_joint
      } 
      if (recon == "marg") {
        model <- models$HMM_2state_2rate_marg
        model_ER <- models$MK_2state_marg
        model_ARD <- models$MK_2state_marg
      }
      
      ## Load the ancestral states data
      anc_states <- readRDS(file = paste("diel_ancestral_states", dataset_variable, name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      anc_rates <- readRDS(file = paste("diel_ancestral_rates", dataset_variable, name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      
      ###############################################################################################################################################
      ### Simulate diel switchs based on derived model parameters ### 
      ###############################################################################################################################################
      
      ## This runs a simulation based on some arguments (wrapper for rTraitDisc)
      ## I can pull the rates directly from the models I've loaded above
      
      if (run_simulations) {
        ## This runs a simulation based on some arguments (wrapper for rTraitDisc)
        ## I can pull the rates directly from the models I've loaded above
        
        if (model_type == "ER") {
          model_rates <- models$MK_2state_marg$solution[2]
          simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = model_rates, states = c("diurnal", "nocturnal"), simulation_numb = sim_numb)
        }
        
        if (model_type == "ARD") {
          model_rates <- models$MK_2state_marg$solution[2:3]
          simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = model_rates, states = c("diurnal", "nocturnal"), simulation_numb = sim_numb, root_state = "nocturnal")
        }
        
        if (model_type == "HR") {
          model_rates <- model$solution
          simulation <- simulateCustom(phylo_tree = trpy_n, model_type = model_type, rates = model_rates, states = c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2"), simulation_numb = sim_numb, root_state = "nocturnal_R2")
        }
        
        saveRDS(simulation, file = here(paste("diel_switch_simulations", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        
      }
      
      if (run_simulations_cumsums) {
        simulation <- readRDS(file = here(paste("diel_switch_simulations", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        
        ## This calculates for each node whether it is a switch from it's parental node
        
        if (model_type == "HR") {
          # Make different for rates vs states
          simulation_rates <- simulation
          
          # Replace numbers with states
          simulation[simulation == "diurnal_R1"] <- "diurnal"
          simulation[simulation == "diurnal_R2"] <- "diurnal"
          simulation[simulation == "nocturnal_R1"] <- "nocturnal"
          simulation[simulation == "nocturnal_R2"] <- "nocturnal"
          
          simulation_rates[simulation_rates == "diurnal_R1"] <- "R1"
          simulation_rates[simulation_rates == "nocturnal_R1"] <- "R1"
          simulation_rates[simulation_rates == "diurnal_R2"] <- "R2"
          simulation_rates[simulation_rates == "nocturnal_R2"] <- "R2"
          
        } 
        
        simulated_transitions <- calculateSimulatedTransitions(simulated_data = simulation, phylo_tree = trpy_n)
        
        ## Function that returns cumsums, rows are now ordered by node.age, whether or it is included as a column in the output
        cumsums <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions, phylo_tree = trpy_n, include_node_age = TRUE)
        
        saveRDS(cumsums, file = here(paste("diel_switch_simulations_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        
        if (model_type == "HR") {
          simulated_transitions_rates <- calculateSimulatedTransitions(simulated_data = simulation_rates, phylo_tree = trpy_n)
          cumsums_rates <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_rates, phylo_tree = trpy_n, include_node_age = TRUE)
          saveRDS(cumsums_rates, file = here(paste("diel_switch_simulations_rates_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        }
      }
      
      ###############################################################################################################################################
      ### Run SIMMAP based on derived model parameters ### 
      ###############################################################################################################################################
      
      if (run_simmap) {
        if (model_type == "ARD") {
          simmaps <- corHMM::makeSimmap(tree = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = model_ARD$solution, nSim = sim_numb, nCores = 5)
        }
        
        if (model_type == "HR") {
          simmaps <- corHMM::makeSimmap(tree = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = model$solution, nSim = sim_numb, nCores = 5)
        }
        
        ## Save SIMMAP
        
        saveRDS(simmaps, file = here(paste("diel_switch_SIMMAP", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "corrHMM.rds", sep = "_")))
      }
      
      if (run_simmap_cumsums) {
        simmaps <- readRDS(file = here(paste("diel_switch_SIMMAP", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "corrHMM.rds", sep = "_")))
        
        ## This uses phytools function to extract tip and node states (should have looked for this before)
        ## However, it gives me different downstream results (maybe because of the order!)
        simmap_simulations <- lapply(seq_along(simmaps), function(x) {
          df <- getStates(simmaps[[x]], type ="both")
          df <- data.frame(nodes = names(df), state = df)
          colnames(df) <- c("node", paste("state", x, sep = "_"))
          return(df)
        }
        )
        
        # Then full_join them together!
        simmap_simulations <- Reduce(full_join, simmap_simulations)
        rownames(simmap_simulations) <- simmap_simulations$node
        simmap_simulations <- simmap_simulations[,!(names(simmap_simulations) %in% c("node"))]
        
        if (model_type == "HR") {
          # Make different for rates vs states
          simmap_simmulations_rates <- simmap_simulations
          
          # Replace numbers with states
          simmap_simulations[simmap_simulations == 1] <- "diurnal"
          simmap_simulations[simmap_simulations == 2] <- "nocturnal"
          simmap_simulations[simmap_simulations == 3] <- "diurnal"
          simmap_simulations[simmap_simulations == 4] <- "nocturnal"
          
          # Replace with rate names
          simmap_simmulations_rates[simmap_simmulations_rates == "diurnal"] <- "R1"
          simmap_simmulations_rates[simmap_simmulations_rates == "nocturnal"] <- "R1"
          
          simmap_simmulations_rates[simmap_simmulations_rates == 1] <- "R1"
          simmap_simmulations_rates[simmap_simmulations_rates == 2] <- "R1"
          simmap_simmulations_rates[simmap_simmulations_rates == 3] <- "R2"
          simmap_simmulations_rates[simmap_simmulations_rates == 4] <- "R2"
          
        } else {
          simmap_simulations[simmap_simulations == 1] <- "diurnal"
          simmap_simulations[simmap_simulations == 2] <- "nocturnal"
        }
        
        
        ## Compare the output of simmaps to the marginal ancestral reconstruction, they should correlate well
        
        data <- simmap_simulations
        data[data == "diurnal"] <- 1
        data[data == "nocturnal"] <- 0
        data <- data %>% mutate_if(is.character,as.numeric)
        simmap_avg <- data.frame(names = row.names(data), lik.anc = rowMeans(data))
        
        lik.anc <- anc_states$lik.anc
        lik.anc$diurnal <- lik.anc$diurnal_R1+lik.anc$diurnal_R2
        lik.anc$simmap_avg_names <- simmap_avg$names[match(row.names(lik.anc), simmap_avg$names)]
        lik.anc$simmap_avg <- simmap_avg$lik.anc[match(row.names(lik.anc), simmap_avg$names)]
        lik.anc$diff <- abs(lik.anc$diurnal - lik.anc$simmap_avg)
        
        comp_anc_methods <- ggplot(lik.anc[], aes(x = diurnal, y = simmap_avg)) + geom_point()
        
        pdf(file = here(paste("outs/Figures/plot_15_SIMMAPvsSimulation", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, ncol(simmap_simulations), "zoom.pdf", sep = "_")), width = 10, height = 10)
        print(comp_anc_methods)
        dev.off()
        
        ## Need to make sure the tips are in the correct order, then phytools::getStates matches what I was doing
        simmap_simulations <- simmap_simulations[row.names(lik.anc),]
        
        ### Calculate transitions and then cummulative transitions through time
        simulated_transitions_simmap <- calculateSimulatedTransitions(simulated_data = simmap_simulations, phylo_tree = trpy_n)
        
        ## Function that returns cumsums, rows are now ordered by node.age, whether or it is included as a column in the output
        cumsums_simmap <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_simmap, phylo_tree = trpy_n, include_node_age = TRUE)
        
        saveRDS(cumsums_simmap, file = here(paste("diel_switch_SIMMAP_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        
        if (model_type == "HR") {
          simulated_transitions_simmap_rates <- calculateSimulatedTransitions(simulated_data = simmap_simmulations_rates, phylo_tree = trpy_n)
          cumsums_simmap_rates <- calculateSimualtedTransitionCumsums(simulated_transitions = simulated_transitions_simmap_rates, phylo_tree = trpy_n, include_node_age = TRUE)
          
          saveRDS(cumsums_simmap_rates, file = here(paste("diel_switch_SIMMAP_rates_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        }
      }
      
    }
  }
  
}





