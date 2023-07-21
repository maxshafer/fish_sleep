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
library(xlsx)
library(gsheet)

#### This script to make some extra figures, based on extinction, fossil, and other data

setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

index_list <- list()
index_list[[1]] <- c("all")
index_list[[2]] <- c("all")
index_list[[3]] <- c("not_mammals")
names(index_list) <- c("fish", "mammals", "tetrapods")

# Set simulation parameters (type and #)
model_types <- "HR"
sim_numb <- 500

recon <- "marg"

anc_rates <- list()
switch.ratio <- list()

for (i in 1:length(index_list)) {
      
      dataset_variable <- names(index_list)[[i]]
      name_variable <- index_list[[i]]
      
      setwd(here())
      
      ###############################################################################################################################################
      ### Load files ### 
      ###############################################################################################################################################
      
      ## Load in the tree
      trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
      
      ## Load the best model, which is the HMM 2 state 2 rate model
      models <- readRDS(file = here(paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_")))
      
      ## I think I will always take the 2 rate model, 3 is too hard to comprehend
      
      ## If I want to use the joint reconstruction then the structure of the data changes (marginal gives % for tips + nodes, joint gives states as vectors)
      if (recon == "joint") {
        model <- models$HMM_2state_2rate_joint
        model_ER <- models$MK_2state_joint
      } 
      if (recon == "marg") {
        model <- models$HMM_2state_2rate_marg
        model_ER <- models$MK_2state_marg
      }
      
      ## Load the ancestral states data
      anc_states[[i]] <- readRDS(file = paste("diel_ancestral_states", dataset_variable, name_variable, Ntip(trpy_n), "species.rds", sep = "_"))

      switch.ratio[[i]] <- switchRatio(ancestral_states = anc_states[[i]], phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = 5)
      
}


## Load data from Rohde et al 2005
url <- 'https://docs.google.com/spreadsheets/d/1mdd1tQ46U6MFIoFkqilHgpQXf24ue93RFoD6M6iZODo/edit#gid=0'
rohde_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
rohde_data <- rohde_data[,c("Time_Ma", "Extinction_Intensity", "Origination_Intensity")]
colnames(rohde_data) <- c("node.age", "extinction", "origination")

## Smooth the data by 5 my (same as mine)
rohde_data$node.age <- cut_width(rohde_data$node.age, 5, label = FALSE)
rohde_data <- rohde_data %>% group_by(node.age) %>% summarise(extinction = mean(extinction), origination = mean(origination))
rohde_data$node.age <- rohde_data$node.age*5


## Make a combined plot
combined.data <- Reduce(rbind, list(switch.ratio[[1]]$data, switch.ratio[[2]]$data, switch.ratio[[3]]$data))
combined.data$group <- c(rep("fish", 49), rep("mammals", 21), rep("tetrapods", 46))
combined.plot <- ggplot(combined.data, aes(x = node.age, y = ratio, group = group, colour = group)) + geom_line(size = 2, alpha = 0.75) + scale_color_viridis(discrete = TRUE, option = "D") + theme_classic() + scale_x_reverse()


extinction <- ggplot(rohde_data, aes(x = node.age, y = extinction)) + geom_line() + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + theme_classic()

origination <- ggplot(rohde_data, aes(x = node.age, y = origination)) + geom_line(colour = "red") + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + ylim(0,30) + theme_classic()

geo_scale <- gggeo_scale(combined.plot, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(combined.plot)$x$range$range)) 

combined.plot / geo_scale / extinction  + plot_layout(nrow = 3, heights = c(5,0.5,5))

combined.plot / geo_scale / origination + plot_layout(nrow = 3, heights = c(5,0.5,5))


### What about fossil data?
file_list <- list.files(here("paleobiodb/"), pattern = "summary")

db_data <- lapply(file_list, function(x) {
  data <- read.csv(file = here(paste("paleobiodb/", x, sep = "")))
  data$node.age <- (data$max_ma + data$min_ma)/2
  
  # bin by age
  data$node.age <- cut_width(data$node.age, 5, label = FALSE)
  data <- data %>% group_by(node.age) %>% summarise(fossil = mean(sampled_in_bin))
  data$node.age <- data$node.age*5
  data$group <- str_split(x, pattern = "_")[[1]][[3]]
  
  return(data)
})

db_data <- Reduce(rbind, db_data)

fossil_plot <- ggplot(db_data[db_data$group != "tetrapoda",], aes(x = node.age, y = fossil, group = group, colour = group)) + geom_line() + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + theme_classic()


combined.plot / geo_scale / fossil_plot + plot_layout(nrow = 3, heights = c(5,0.5,5))

