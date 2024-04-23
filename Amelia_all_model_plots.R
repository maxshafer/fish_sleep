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
# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")

#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("cetaceans_full.csv"))

#run cor_maor_ard and cor_model ard
#these are the basis for the artiodactyla models and cetacean models respectively



# #Section X: Evaluate and plot the 1k model results -------------------------------
source("scripts/Amelia_functions.R")

#####max crep dataset, cetaceans, all models#########
ER_results <- readRDS(here("cetaceans_tree_max_crep_traits_ER_models"))
SYM_results <- readRDS(here("cetaceans_tree_max_crep_traits_SYM_models"))
ARD_results <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models"))


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
all_maxcrep_plot <- ggplot(df_crep, aes(x = model, y = likelihoods)) + geom_point() + geom_boxplot() + geom_jitter() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_maxcrep_all_models_1k_trees_plot.png")
all_maxcrep_plot
dev.off()

#########max_dinoc, cetaceans, all models ###########
ER_SYM_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models"))
ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models"))


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

all_maxdinoc_plot <- ggplot(df_dinoc, aes(x = model, y = likelihoods)) + geom_point() + geom_boxplot() + geom_jitter() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_maxdinoc_all_models_1k_trees_plot.png")
all_maxdinoc_plot
dev.off()

#####plot all 8 cetacean models##############

#combine dfs
df_cet <- rbind(df_dinoc, df_crep)
df_cet$dataset <- c(rep("max_crep", 4000), rep("max_dinoc", 4000))
all_models_plot <- ggplot(df_cet, aes(x = model, y = likelihoods, colour = dataset)) + geom_point(alpha = 0.2) + geom_boxplot() + geom_jitter() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
all_models_plot

#save plot as a PNG
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_all_8_models_1k_trees_plot.png")
all_models_plot
dev.off()


# Section 2: Plot transition rates ----------------------------------------

#function to extract the transition rates from all models

#####max_crep_ER#####
rates <- unlist(lapply(ER_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_crep_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Crep/Cath -> Di", "Crep/Cath -> Noc", "Di -> Crep/Cath", "Di -> Noc", "Noc -> Crep/Cath", "Noc -> Di"), 1000)

ggplot(rates_df, aes(x= solution, y = rates, colour = model)) + geom_point() + theme(axis.text.x = element_text(angle = 90))

#####max_dinoc_ARD#####
ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models"))
rates <- unlist(lapply(ARD_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Cath -> Di", "Cath -> Noc", "Di -> Cath", "Di -> Noc", "Noc -> Cath", "Noc -> Di"), 1000)

ggplot(rates_df, aes(x= solution, y = rates, colour = solution)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) +  scale_y_continuous(trans='log10')

#many values are close to zero, should we make a cutoff? ie anything less than 0.001? 

####max_dinoc_bridge_only####
bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models"))
rates <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Cath -> Di", "Cath -> Noc", "Di -> Cath", "Di -> Noc", "Noc -> Cath", "Noc -> Di"), 1000)

ggplot(rates_df, aes(x= solution, y = rates, colour = solution)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(trans='log10') #+  scale_y_sqrt()

bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models"))

