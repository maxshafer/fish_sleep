# Section 0: Import and run packages -------------------------------------------------
# For retreiving data from the open tree of life
library(rotl)
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

#extra
library(ggplot2)
library(scales)
library(RColorBrewer)
#install.packages("ggforce")
library(ggforce)
library(forcats)
library(tidyr)
# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")
#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))



# Section 1: Cetaceans, max_crep dataset#########
ER_results <- readRDS(here("cetaceans_tree_max_crep_traits_ER_models.rds"))
SYM_results <- readRDS(here("cetaceans_tree_max_crep_traits_SYM_models.rds"))
ARD_results <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))


ER_likelihoods <- unlist(lapply(ER_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(SYM_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_crep", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM max_crep", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD max_crep", likelihoods = ARD_likelihoods)
df4 <- data.frame(model = "bridge_only max_crep", likelihoods = bridge_only_likelihoods)
df_crep <- rbind(df1, df2, df3, df4)

all_maxcrep_plot <- ggplot(df_crep, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_maxcrep_all_models_1k_trees_plot.png", units = "cm", height = 15, width = 15, res = 225)
all_maxcrep_plot
dev.off()

# ###Cetaceans_max_crep_ER_rates### ----------------------------------------------------------------

model_results <- readRDS(here("cetaceans_tree_max_crep_traits_ER_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Cath/crep -> Di", "Cath/crep -> Noc", "Di -> Cath/crep", "Di -> Noc", "Noc -> Cath/crep", "Noc -> Di"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_crep_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Cetaceans_max_crep_SYM_rates### ----------------------------------------------------------------
model_results <- readRDS(here("cetaceans_tree_max_crep_traits_SYM_models.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_crep_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Cetaceans_max_crep_ARD_rates### ####
ARD_results <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models.rds"))

rates <- unlist(lapply(ARD_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which state transition they represent
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc, Cath -> Cath)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)


rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)
#ggplot(rates_df, aes(x= solution, y = rates, colour = solution)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) #+ scale_y_continuous(trans='log10') #+  scale_y_sqrt()
rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_crep_ARD/rates_plots_MAXCREP_ARD", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i) +  theme(strip.text = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")))
  dev.off()
}


# ###Cetacean_max_crep_bridge_only_rates### -----------------------------
bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))
rates <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc, Cath -> Cath)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"), 1000)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")

#ggplot(rates_df, aes(x= solution, y = rates, colour = solution)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) #+ scale_y_continuous(trans='log10') #+  scale_y_sqrt()
rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:4){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_crep_bridge_only/rates_plots_MAXCREP_BRIDGE", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i) +  theme(strip.text = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")))
  dev.off()
}



# Section 2: Cetaceans, max_dinoc dataset ###########
ER_SYM_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models.rds"))
ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models.rds"))


ER_likelihoods <- unlist(lapply(ER_SYM_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(ER_SYM_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_dinoc", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM max_dinoc", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD max_dinoc", likelihoods = ARD_likelihoods)
df4 <- data.frame(model = "bridge_only max_dinoc", likelihoods = bridge_only_likelihoods)
df_dinoc <- rbind(df1, df2, df3, df4)

all_maxdinoc_plot <- ggplot(df_dinoc, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#00BFC4") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_maxdinoc_all_models_1k_trees_plots.png", units = "cm", height = 15, width = 15, res = 225)
all_maxdinoc_plot
dev.off()

# ###Cetaceans_max_dinoc_ER_rates### ----------------------------------------------------------------

model_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_dinoc_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Cetaceans_max_dinoc_SYM_rates### ----------------------------------------------------------------
model_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_dinoc_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ####Cetaceans_max_dinoc_ARD_rates### #####

ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models.rds"))
rates <- unlist(lapply(ARD_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
  
#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_dinoc_ARD/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i) +  theme(strip.text = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")))
  dev.off()
  }

# ####Cetaceans_max_dinoc_bridge_only_rates###########
bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models.rds"))
rates <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Cath -> Noc"), 1000)

#ggplot(rates_df, aes(x= solution, y = rates, colour = solution)) + geom_boxplot() + geom_point() + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(trans = "log10")
#ggplot(rates_df, aes(x = rates, fill = solution)) + geom_histogram(position = "dodge") + scale_x_continuous(trans = "log2")
#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:4){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/cetacean_max_dinoc_bridge_only/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i)+  theme(strip.text = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")))
  dev.off()
}



# Section 3: Artiodactyla, max_crep--------------------------------------

#use below for Cox et al data
# ER_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_ER_SYM_models.rds"))
# SYM_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_ER_SYM_models.rds"))
# ARD_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_ARD_models.rds"))
# bridge_only_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_bridge_only_models.rds"))

#use below for primary source data
ER_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))
SYM_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))
ARD_results <- readRDS(here("artiodactyla_new_max_crep_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

ER_likelihoods <- unlist(lapply(ER_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(SYM_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_crep", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM max_crep", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD max_crep", likelihoods = ARD_likelihoods)
df4 <- data.frame(model = "bridge_only max_crep", likelihoods = bridge_only_likelihoods)
df_crep <- rbind(df1, df2, df3, df4)

all_maxcrep_artio_plot <- ggplot(df_crep, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/new_artio_maxcrep_all_models_plot.png", units = "cm", height = 15, width = 15, res = 225)
all_maxcrep_artio_plot
dev.off()

# ###Artio_max_crep_ER_rates### ----------------------------------------------------------------

#use below for Cox et al data
#model_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_ER_SYM_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_crep_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Artio_max_crep_SYM_rates### ----------------------------------------------------------------

#use below for Cox data
#model_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_ER_SYM_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM__bridge_only_models.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_crep_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Artio_max_crep_ARD_rates### ----------------------------------------------------------------

#use below for Cox data
#model_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_ARD_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_crep_traits_ARD_models.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution)

for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_crep_ARD/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

##To colour each of the rates by their model likelihood
#there are six rates for each model (except bridge_only which has 4), so duplicate the likelihood of each model 6 times
likelihood<- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 6, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)


# ###Artio_max_crep_bridge_only_rates### ----------------------------------------------------------------

#use below for Cox data
#model_results <- readRDS(here("artiodactyla_Cox_max_crep_traits_bridge_only_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")


#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

for(i in 1:4){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_crep_bridge_only/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

##To colour each of the rates by their model likelihood
#there are four rates for each bridge_only model, so duplicate the likelihood of each model 4 times
likelihood<- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 4, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + facet_wrap(~solution) #+ theme_bw() + theme(plot.title = element_blank(), legend.position = "none")



# Section 4: Artiodactyla, max_dinoc--------------------------------------
#use below for Cox et al, data
# ER_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ER_SYM_models.rds"))
# SYM_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ER_SYM_models.rds"))
# ARD_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ARD_models.rds"))
# bridge_only_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_bridge_only_models.rds"))

#use below for new primary source data
ER_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ER_SYM_bridge_only_models.rds"))
SYM_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ER_SYM_bridge_only_models.rds"))
ARD_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ER_SYM_bridge_only_models.rds"))


ER_likelihoods <- unlist(lapply(ER_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(SYM_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_dinoc", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM max_dinoc", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD max_dinoc", likelihoods = ARD_likelihoods)
df4 <- data.frame(model = "bridge_only max_dinoc", likelihoods = bridge_only_likelihoods)
df_dinoc <- rbind(df1, df2, df3, df4)

all_maxdinoc_artio_plot <- ggplot(df_dinoc, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#00BFC4") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/new_artio_maxdinoc_all_models_plot.png",  units = "cm", height = 15, width = 15, res = 225)
all_maxdinoc_artio_plot
dev.off()



# ###Artio_max_dinoc_ER_rates### ----------------------------------------------------------------

#use below for Cox data
#model_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ER_SYM_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i)) 
  dev.off()
}

# ###Artio_max_dinoc_SYM_rates### ----------------------------------------------------------------

#use below for Cox data
#model_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ER_SYM_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Artio_max_dinoc_ARD_rates### ----------------------------------------------------------------

#use below for Cox
#model_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ARD_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ARD_models.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution)

#save out as png, remember to change file name
for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_ARD/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

##To colour each of the rates by their model likelihood
#there are six rates for each model (except bridge_only which has 4), so duplicate the likelihood of each model 6 times
likelihood<- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 6, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)


# ###Artio_max_dinoc_bridge_only_rates### ----------------------------------------------------------------

#use below for Cox data
#model_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_bridge_only_models.rds"))
#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Cath -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution)

for(i in 1:4){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_bridge_only/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

##To colour each of the rates by their model likelihood
#there are four rates for each bridge_only model, so duplicate the likelihood of each model 4 times
likelihood<- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 4, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)


# Section 5: Artiodactyla without cetaceans, max_crep--------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("artiodactyla_Amelia_max_crep_traits_ER_SYM_models.rds"))

ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_crep", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM max_crep", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD max_crep", likelihoods = ARD_likelihoods)
df4 <- data.frame(model = "bridge_only max_crep", likelihoods = bridge_only_likelihoods)
df_crep <- rbind(df1, df2, df3, df4)

all_maxcrep_artio_minus_cet_plot <- ggplot(df_crep, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cox_artio_minus_cet_maxcrep_all_models_plot.png",  units = "cm", height = 15, width = 15, res = 225)
all_maxcrep_artio_minus_cet_plot
dev.off()

# ###Artio_minus_cetaceans_max_crep_ER_rates### ----------------------------------------------------------------

model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("artiodactyla_Amelia_max_crep_traits_ER_SYM_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_crep_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Artio_minus_cetaceans_max_crep_SYM_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_crep_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}


# ###Artio_minus_cetaceans_max_crep_ARD_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)


for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_crep_ARD/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}


# ###Artio_minus_cetaceans_max_crep_bridge_only_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_brdige_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:4){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_crep_bridge_only/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}



# Section 6: Artiodactyla without cetaceans, max_dinoc--------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_dinoc", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM max_dinoc", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD max_dinoc", likelihoods = ARD_likelihoods)
df4 <- data.frame(model = "bridge_only max_dinoc", likelihoods = bridge_only_likelihoods)
df_dinoc <- rbind(df1, df2, df3, df4)

all_maxdinoc_artio_minus_cet_plot <- ggplot(df_dinoc, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#00BFC4") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cox_artio_minus_cet_maxdinoc_all_models_plot.png",  units = "cm", height = 15, width = 15, res = 225)
all_maxdinoc_artio_minus_cet_plot
dev.off()

# ###Artio_minus_cetaceans_max_dinoc_ER_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_dinoc_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Artio_minus_cetaceans_max_dinoc_SYM_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_dinoc_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}


# ###Artio_minus_cetaceans_max_dinoc_ARD_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Noc -> Di", "Cath -> Noc", "Di -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:6){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_dinoc_ARD/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ###Artio_minus_cetaceans_max_dinoc_bridge_only_rates### ----------------------------------------------------------------

#use below for Cox data
model_results <- readRDS(here("artiodactyla_minus_cetaceans_Cox_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))
#use below for Amelia primary source data
#model_results <- readRDS(here("filename.rds"))

rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath", "Noc -> Cath", "Cath -> Di", "Cath -> Noc"), 1000)
rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")

rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 3, page = 2)

for(i in 1:4){
  png(paste0("C:/Users/ameli/Downloads/rates_dump/artiodactyla_minus_cet_max_dinoc_bridge_only/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}




# Section 7: Final comparisons --------------------------------------------
# ####Likelihoods of all 8 cetacean models#################

#combine dfs
df_cet <- rbind(df_dinoc, df_crep)
df_cet$Dataset <- c(rep("max_dinoc", 4000), rep("max_crep", 4000))
all_models_plot <- ggplot(df_cet, aes(x = fct_inorder(model), y = likelihoods, colour = Dataset)) + geom_jitter(alpha = 0.6) + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") +
  scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only","Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) +
  scale_colour_discrete(labels = c("Max crep", "Max di/noc"))
all_models_plot

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_all_8_models_1k_trees_plot.png", width = 15, height = 15, units = "cm", res = 900)
all_models_plot
dev.off()



#function to extract the transition rates from all models


# Section 8: Rates scatter plots -----------------------------------------
#we want to know how the distribution of rates to and from a state compare to each other
#for example, when rates from noc <- cath are low are cath -> noc also low?


# ### Rates Scatterplot Max_crep ARD ###----------------------------
#use below for cetaceans
#ARD_results <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models.rds"))
#use below for artiodactyla - primary source data
ARD_results <- readRDS(here("artiodactyla_new_max_crep_traits_ARD_models.rds"))

rates <- unlist(lapply(ARD_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which state transition they represent
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc, Cath -> Cath)
rates_df <- rates_df[!(is.na(rates_df$rates)),]
rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)


#we want to plot all the di to noc and noc to di points
rates_pivot <- rates_df %>% pivot_wider(names_from = solution, values_from = rates)
rates_pivot <- unnest(rates_pivot, cols = c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"))

#add in likelihoods of the models for each set of rates
likelihoods <- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
likelihoods_df <- as.data.frame(likelihoods)
rates_pivot <- cbind( rates_pivot, likelihoods_df)

colnames(rates_pivot) <- c("model", "DiCath", "NocCath", "CathDi", "NocDi", "CathNoc", "DiNoc", "likelihoods")

scatter_dinoc <- ggplot(rates_pivot, aes(x = NocDi, y = DiNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "seagreen4", high = "yellow") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)))

scatter_dicath <- ggplot(rates_pivot, aes(x = CathDi, y = DiCath, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "deeppink4", high = "plum") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

scatter_noccath <- ggplot(rates_pivot, aes(x = NocCath, y = CathNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "dodgerblue4", high = "cyan") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 



#save out as a png, remember to change filename
png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_crep_ARD_scatter_dinoc.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dinoc
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_crep_ARD_scatter_dicath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dicath
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_crep_ARD_scatter_noccath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_noccath
dev.off()

# ### Rates Scatterplot Max_crep bridge_only ###----------------------------
#use below for cetacean
#bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))
#use below for artiodactyla -primary source data 
bridge_only_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which state transition they represent
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc, Cath -> Cath)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"), 1000)


#we want to plot all the di to noc and noc to di points
rates_pivot <- rates_df %>% pivot_wider(names_from = solution, values_from = rates)
rates_pivot <- unnest(rates_pivot, cols = c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"))

#add in likelihoods of the models for each set of rates
likelihoods <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
likelihoods_df <- as.data.frame(likelihoods)
rates_pivot <- cbind( rates_pivot, likelihoods_df)

colnames(rates_pivot) <- c("model", "DiCath", "NocCath", "CathDi", "CathNoc", "likelihoods")

scatter_dicath <- ggplot(rates_pivot, aes(x = CathDi, y = DiCath, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "deeppink4", high = "plum") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

scatter_noccath <- ggplot(rates_pivot, aes(x = NocCath, y = CathNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "dodgerblue4", high = "cyan") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

#to be continued
# scatter_1 <- ggplot(rates_pivot, aes(x = NocCath, y = DiCath)) + geom_point(colour = "grey10") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))
# scatter_2 <- ggplot(rates_pivot, aes(x = NocCath, y = CathDi)) + geom_point(colour = "grey20") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))
# scatter_3 <- ggplot(rates_pivot, aes(x = CathNoc, y = DiCath)) + geom_point(colour = "grey30") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))
# scatter_4 <- ggplot(rates_pivot, aes(x = CathNoc, y = CathDi)) + geom_point(colour = "grey40") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))


#save out as a png

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_crep_bridge_only_scatter_dicath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dicath
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_crep_bridge_only_scatter_noccath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_noccath
dev.off()

# ### Rates Scatterplot Cetacean Max_dinoc ARD ###----------------------------
#use below for cetacean data
#ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models.rds"))
#use below for artiodactyla - primary source data
ARD_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ARD_models.rds"))

rates <- unlist(lapply(ARD_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which state transition they represent
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc, Cath -> Cath)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)


#we want to plot all the di to noc and noc to di points
rates_pivot <- rates_df %>% pivot_wider(names_from = solution, values_from = rates)
rates_pivot <- unnest(rates_pivot, cols = c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"))

#add in likelihoods of the models for each set of rates
likelihoods <- unlist(lapply(ARD_results$ARD_model, function(x) returnLikelihoods(model = x)))
likelihoods_df <- as.data.frame(likelihoods)
rates_pivot <- cbind( rates_pivot, likelihoods_df)

colnames(rates_pivot) <- c("model", "DiCath", "NocCath", "CathDi", "NocDi", "CathNoc", "DiNoc", "likelihoods")

scatter_dinoc <- ggplot(rates_pivot, aes(x = NocDi, y = DiNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "seagreen4", high = "yellow") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)))

scatter_dicath <- ggplot(rates_pivot, aes(x = CathDi, y = DiCath, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "deeppink4", high = "plum") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

scatter_noccath <- ggplot(rates_pivot, aes(x = NocCath, y = CathNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "dodgerblue4", high = "cyan") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

#save out as a png, remember to change file name
png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/new_artio_max_dinoc_ARD_scatter_dinoc.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dinoc
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/new_artio_max_dinoc_ARD_scatter_dicath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dicath
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/new_artio_max_dinoc_ARD_scatter_noccath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_noccath
dev.off()

# ### Rates Scatterplot Max_dinoc bridge_only ###----------------------------
#use below for cetacean data
#bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models.rds"))
#use below for artiodactyla -primary source data
bridge_only_results <- readRDS(here("artiodactyla_new_max_dinoc_traits_ER_SYM_bridge_only_models.rds"))

rates <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which state transition they represent
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc, Cath -> Cath)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"), 1000)


#we want to plot all the di to noc and noc to di points
rates_pivot <- rates_df %>% pivot_wider(names_from = solution, values_from = rates)
rates_pivot <- unnest(rates_pivot, cols = c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"))

#add in likelihoods of the models for each set of rates
likelihoods <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
likelihoods_df <- as.data.frame(likelihoods)
rates_pivot <- cbind( rates_pivot, likelihoods_df)

colnames(rates_pivot) <- c("model", "DiCath", "NocCath", "CathDi", "CathNoc", "likelihoods")

scatter_dicath <- ggplot(rates_pivot, aes(x = CathDi, y = DiCath, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "deeppink4", high = "plum") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

scatter_noccath <- ggplot(rates_pivot, aes(x = NocCath, y = CathNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "dodgerblue4", high = "cyan") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 


#save out as a png, remember to change file name

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_dinoc_bridge_only_scatter_dicath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dicath
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/artio_max_dinoc_bridge_only_scatter_noccath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_noccath
dev.off()





# #Section Y: Plot max_clade_cred tree results------------------------------------
model_results <- readRDS(here("artiodactyla_max_clade_cred_six_state_traits_ER_SYM_ARD_models.rds"))
filename <- "artiodactyla_max_clade_cred_six_state_traits_ER_SYM_ARD_models.rds"

# ER_results <- model_results$ER_model
# SYM_results <- model_results$SYM_model
# ARD_results <- model_results$ARD_model
# bridge_results <- model_results$bridge_only
# hidden_rates <- model_results$hidden_rates

likelihoods <- unlist(lapply(model_results, function(x) returnLikelihoods(model = x)))
likelihoods <- as.data.frame(likelihoods)
likelihoods$model <- rownames(likelihoods)

ggplot(likelihoods, aes(x = model, y = likelihoods)) + geom_point()
