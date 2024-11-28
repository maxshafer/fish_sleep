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

setwd(here())
source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Set the working directory and source the functions (not used yet)
setwd("C:/Users/ameli/OneDrive/Documents/R_projects/fish_sleep/1k_model_results")


#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))



# Section 1: Cetaceans, max_crep dataset#########
#extract and plot the likelihoods
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

#extract and plot the AICc scores
ER_results <- readRDS(here("cetaceans_tree_max_crep_traits_ER_models.rds"))
SYM_results <- readRDS(here("cetaceans_tree_max_crep_traits_SYM_models.rds"))
ARD_results <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))


ER_AICc <- unlist(lapply(ER_results$ER_model, function(x) returnAICc(model = x)))
SYM_AICc <- unlist(lapply(SYM_results$SYM_model, function(x) returnAICc(model = x)))
ARD_AICc <- unlist(lapply(ARD_results$ARD_model, function(x) returnAICc(model = x)))
bridge_only_AICc  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnAICc(model = x)))

## combine and plot

df1 <- data.frame(model = "ER max_crep", AICc = ER_AICc)
df2 <- data.frame(model = "SYM max_crep", AICc = SYM_AICc)
df3 <- data.frame(model = "ARD max_crep", AICc = ARD_AICc)
df4 <- data.frame(model = "bridge_only max_crep", AICc = bridge_only_AICc)
df_crep <- rbind(df1, df2, df3, df4)

all_maxcrep_plot <- ggplot(df_crep, aes(x = fct_inorder(model), y = AICc)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AICc") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Cetaceans")

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_three_state_maxcrep_AICc_all_models_1k_trees_plot.png", units = "cm", height = 15, width = 15, res = 225)
all_maxcrep_plot
dev.off()

#extract and plot the AIC scores
ER_results <- readRDS(here("cetaceans_tree_max_crep_traits_ER_models.rds"))
SYM_results <- readRDS(here("cetaceans_tree_max_crep_traits_SYM_models.rds"))
ARD_results <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))


ER_AIC <- unlist(lapply(ER_results$ER_model, function(x) returnAIC(model = x)))
SYM_AIC <- unlist(lapply(SYM_results$SYM_model, function(x) returnAIC(model = x)))
ARD_AIC <- unlist(lapply(ARD_results$ARD_model, function(x) returnAIC(model = x)))
bridge_only_AIC  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnAIC(model = x)))

## combine and plot
df1 <- data.frame(model = "ER max_crep", AIC = ER_AIC)
df2 <- data.frame(model = "SYM max_crep", AIC = SYM_AIC)
df3 <- data.frame(model = "ARD max_crep", AIC = ARD_AIC)
df4 <- data.frame(model = "bridge_only max_crep", AIC = bridge_only_AIC)
df_crep <- rbind(df1, df2, df3, df4)

all_maxcrep_plot <- ggplot(df_crep, aes(x = fct_inorder(model), y = AIC)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AIC") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Cetaecans")

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_three_state_maxcrep_AIC_all_models_1k_trees_plot.png", units = "cm", height = 15, width = 15, res = 225)
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

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_crep_ARD/all_rates_ARD.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 2, nrow = 3)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_crep_ARD/all_rates_ARD_violinplot_median.png", width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Cetaceans") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
dev.off()


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

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_crep_bridge_only/all_rates_bridge_only.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 2, nrow = 2)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_crep_bridge_only/all_rates_bridge_only_violinplot_median.png", width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Cetaceans") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
dev.off()

#save out as png, remember to change file name
for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ARD/rates_plots", i, ".png"), width=30,height=5,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 5, nrow = 1, page = i))
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

#extract and plot AICc scores from 1k models
ER_SYM_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models.rds"))
ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models.rds"))

ER_AICc <- unlist(lapply(ER_SYM_results$ER_model, function(x) returnAICc(model = x)))
SYM_AICc  <- unlist(lapply(ER_SYM_results$SYM_model, function(x) returnAICc(model = x)))
ARD_AICc  <- unlist(lapply(ARD_results$ARD_model, function(x) returnAICc(model = x)))
bridge_only_AICc  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnAICc(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", AICc = ER_AICc)
df2 <- data.frame(model = "SYM", AICc = SYM_AICc)
df3 <- data.frame(model = "ARD", AICc = ARD_AICc)
df4 <- data.frame(model = "Bridge_only", AICc = bridge_only_AICc)
df_s6x <- rbind(df1, df2, df3, df4)

AICc_cet_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = AICc)) + geom_jitter(alpha = 0.6, color = "#00BFC4") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Cetaceans") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_three_state_max_dinoc_all_models_AICc.png",  units = "cm", height = 15, width = 15, res = 225)
AICc_cet_plot
dev.off()

#extract and plot AIC
ER_SYM_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models.rds"))
ARD_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models.rds"))
bridge_only_results <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models.rds"))

ER_AIC <- unlist(lapply(ER_SYM_results$ER_model, function(x) returnAIC(model = x)))
SYM_AIC  <- unlist(lapply(ER_SYM_results$SYM_model, function(x) returnAIC(model = x)))
ARD_AIC  <- unlist(lapply(ARD_results$ARD_model, function(x) returnAIC(model = x)))
bridge_only_AIC  <- unlist(lapply(bridge_only_results$bridge_only_model, function(x) returnAIC(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", AIC = ER_AIC)
df2 <- data.frame(model = "SYM", AIC = SYM_AIC)
df3 <- data.frame(model = "ARD", AIC = ARD_AIC)
df4 <- data.frame(model = "Bridge only", AIC = bridge_only_AIC)
df_s6x <- rbind(df1, df2, df3, df4)

AIC_cet_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = AIC)) + geom_jitter(alpha = 0.6, color = "#00BFC4") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "AIC score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Cetaceans") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/cetaceans_three_state_maxdinoc_AIC.png",  units = "cm", height = 15, width = 15, res = 225)
AIC_cet_plot
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

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_dinoc_SYM/all_rates_SYM.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 3, nrow = 2)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_dinoc_SYM/all_rates_SYM_violinplot_median.png", width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Cetaceans") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
dev.off()

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

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_dinoc_ARD/all_rates_ARD.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 3, nrow = 2)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/cetacean_max_dinoc_ARD/all_rates_ARD_violinplot_median.png", width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Cetaceans") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
dev.off()

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
likelihood <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 6, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_crep_ARD/likelihood_rates_plot.png", width=20,height=10,units="cm",res=200)
ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
dev.off()

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
likelihood<- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 4, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_crep_bridge_only/likelihood_rates_plot.png", width=20,height=10,units="cm",res=200)
ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
dev.off()


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
likelihood<- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 6, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_ARD/likelihood_rates_plot.png", width=20,height=12,units="cm",res=200)
ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
dev.off()

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
likelihood<- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
likelihood <- as.data.frame(likelihood)
index <- rep(1:nrow(likelihood), times = c(rep(x = 4, times = 1000)))
likelihood <- likelihood[index,]

#add this column to the rates dataframe
rates_df <- cbind(rates_df, likelihood)
rates_df$likelihood <- as.integer(rates_df$likelihood)
rates_df$likelihood <- as.factor(rates_df$likelihood)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_bridge_only/likelihood_rates_plot.png", width=20,height=14,units="cm",res=200)
ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
dev.off()


# Section 4: Artiodactyla, six_state--------------------------------------
#use below for Cox et al, data
# ER_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ER_SYM_models.rds"))
# SYM_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ER_SYM_models.rds"))
# ARD_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_ARD_models.rds"))
# bridge_only_results <- readRDS(here("artiodactyla_Cox_max_dinoc_traits_bridge_only_models.rds"))

#use below for new primary source data
model_results <- readRDS(here("artiodactyla_200_test_six_state_traits_ER_SYM_ARD_bridge_only_models.rds"))

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

all_six_states_artio_plot <- ggplot(df_dinoc, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "purple") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") #+ scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/new_artio_six_state_all_models_plot.png",  units = "cm", height = 15, width = 15, res = 225)
all_six_states_artio_plot
dev.off()



# ###Artio_six_state_ER_rates_NOT DONE### ----------------------------------------------------------------

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

# ###Artio_six_state_SYM_rates_NOT DONE### ----------------------------------------------------------------

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

# ###Artio_six_state_ARD_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_200_test_six_state_traits_ER_SYM_ARD_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]
row.names(rates_df) <- 1:length(rates_df$rates)

rates_df$solution <- rep(c("Cath/crep -> Cath", "Di-> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 20)

rates_colours <- rep(c("darkorchid1", "deeppink3", "maroon1", "dodgerblue2", "mediumpurple1", "darkorchid4","orchid1", "firebrick4", "slateblue1", "thistle3", "deeppink4", "orchid3", "darkorange","seagreen3","springgreen1","maroon3","firebrick1","darkorange3","olivedrab1","brown","dodgerblue4","slateblue4","seagreen4","olivedrab3","gold","mediumpurple4","thistle4","springgreen4","brown1","gold3"), 20)

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_six_state_ARD/all_rates.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)
dev.off()

#save out as png, remember to change file name
for(i in 1:30){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_six_state_ARD/rates_plots", i, ".png"), width=15,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ##To colour each of the rates by their model likelihood
# #there are six rates for each model (except bridge_only which has 4), so duplicate the likelihood of each model 6 times
# likelihood<- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
# likelihood <- as.data.frame(likelihood)
# index <- rep(1:nrow(likelihood), times = c(rep(x = 6, times = 1000)))
# likelihood <- likelihood[index,]
# 
# #add this column to the rates dataframe
# rates_df <- cbind(rates_df, likelihood)
# rates_df$likelihood <- as.integer(rates_df$likelihood)
# rates_df$likelihood <- as.factor(rates_df$likelihood)
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_ARD/likelihood_rates_plot.png", width=20,height=12,units="cm",res=200)
# ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
# dev.off()

# ###Artio_six_state_bridge_only_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("artiodactyla_200_test_six_state_traits_ER_SYM_ARD_bridge_only_models.rds"))

rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "max_dinoc_bridge_only"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]
rownames(rates_df) <- 1:length(rates_df$rates)

#exclude transitions from Di <-> Noc, Di <-> Noc/crep, Noc <-> Di/crep, Di/crep <-> Noc/crep
rates_df$solution <- rep(c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Noc -> Noc/crep"), 20)
rates_colours <- rep(c("darkorchid1", "deeppink3", "maroon1", "dodgerblue2", "mediumpurple1", "darkorchid4","orchid1", "slateblue1", "thistle3", "deeppink4", "orchid3", "darkorange","maroon3","darkorange3","brown","dodgerblue4","slateblue4","gold","mediumpurple4","thistle4","brown1","gold3"), 20)

#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution)

for(i in 1:22){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_six_state_bridge_only/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i))
  dev.off()
}

# ##To colour each of the rates by their model likelihood
# #there are four rates for each bridge_only model, so duplicate the likelihood of each model 4 times
# likelihood<- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
# likelihood <- as.data.frame(likelihood)
# index <- rep(1:nrow(likelihood), times = c(rep(x = 4, times = 1000)))
# likelihood <- likelihood[index,]
# 
# #add this column to the rates dataframe
# rates_df <- cbind(rates_df, likelihood)
# rates_df$likelihood <- as.integer(rates_df$likelihood)
# rates_df$likelihood <- as.factor(rates_df$likelihood)
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_bridge_only/likelihood_rates_plot.png", width=20,height=14,units="cm",res=200)
# ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
# dev.off()


# Section 4.5: Ruminants, six_state--------------------------------------

#plot likelihoods of all 1k model results
model_results <- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))

ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD", likelihoods = ARD_likelihoods)
df_s6x <- rbind(df1, df2, df3)

six_states_likelihood_artio_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "salmon") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "Log likelihood") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Ruminantia") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/ruminants_six_state_all_models_plot.png",  units = "cm", height = 15, width = 15, res = 225)
six_states_likelihood_artio_plot
dev.off()

#extract and plot AICc scores from 1k models
model_results <- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))

ER_AICc <- unlist(lapply(model_results$ER_model, function(x) returnAICc(model = x)))
SYM_AICc  <- unlist(lapply(model_results$SYM_model, function(x) returnAICc(model = x)))
ARD_AICc  <- unlist(lapply(model_results$ARD_model, function(x) returnAICc(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", AICc = ER_AICc)
df2 <- data.frame(model = "SYM", AICc = SYM_AICc)
df3 <- data.frame(model = "ARD", AICc = ARD_AICc)
df_s6x <- rbind(df1, df2, df3)

six_states_AICc_artio_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = AICc)) + geom_jitter(alpha = 0.6, color = "salmon") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Ruminantia") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/ruminants_six_state_all_models_AICc.png",  units = "cm", height = 15, width = 15, res = 225)
six_states_AICc_artio_plot
dev.off()

#extract and plot AIC
model_results <- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))

ER_AIC <- unlist(lapply(model_results$ER_model, function(x) returnAIC(model = x)))
SYM_AIC  <- unlist(lapply(model_results$SYM_model, function(x) returnAIC(model = x)))
ARD_AIC  <- unlist(lapply(model_results$ARD_model, function(x) returnAIC(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", AIC = ER_AIC)
df2 <- data.frame(model = "SYM", AIC = SYM_AIC)
df3 <- data.frame(model = "ARD", AIC = ARD_AIC)
df_s6x <- rbind(df1, df2, df3)

six_states_AIC_artio_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = AIC)) + geom_jitter(alpha = 0.6, color = "salmon") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "AIC score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only")) + ggtitle("Ruminantia") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/ruminants_six_state_all_models_AIC.png",  units = "cm", height = 15, width = 15, res = 225)
six_states_AIC_artio_plot
dev.off()

# ###Ruminants_six_state_ER_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Cath/crep -> Cath", "Di-> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 20)
rates_colours <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf", "#bd5c35", "#d37a57", "#e79979", "#fbb89d",  "#a63d13", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")

#alternative colour palette
#rates_colours <- rep(c("darkorchid1", "deeppink3", "maroon1", "dodgerblue2", "mediumpurple1", "darkorchid4","orchid1", "firebrick4", "slateblue1", "thistle3", "deeppink4", "orchid3", "darkorange","seagreen3","springgreen1","maroon3","firebrick1","darkorange3","olivedrab1","brown","dodgerblue4","slateblue4","seagreen4","olivedrab3","gold","mediumpurple4","thistle4","springgreen4","brown1","gold3"), 20)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ER/all_rates_ER.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 5, nrow = 6)
dev.off() 

for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 1, nrow = 1, page = i)) 
  dev.off()
}

# ###Ruminants_six_state_SYM_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- rep(c("Cath/crep -> Cath", "Di-> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 20)
rates_colours <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf", "#bd5c35", "#d37a57", "#e79979", "#fbb89d",  "#a63d13", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")

#alternative colour palette
#rates_colours <- rep(c("darkorchid1", "deeppink3", "maroon1", "dodgerblue2", "mediumpurple1", "darkorchid4","orchid1", "firebrick4", "slateblue1", "thistle3", "deeppink4", "orchid3", "darkorange","seagreen3","springgreen1","maroon3","firebrick1","darkorange3","olivedrab1","brown","dodgerblue4","slateblue4","seagreen4","olivedrab3","gold","mediumpurple4","thistle4","springgreen4","brown1","gold3"), 20)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_SYM/all_rates_SYM.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 5, nrow = 6)
dev.off()

for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 5, nrow = 1, page = i))
  dev.off()
}

# ###Ruminants_six_state_ARD_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]
row.names(rates_df) <- 1:length(rates_df$rates)

rates_df$solution <- rep(c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 20)
rates_colours <- rep(c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab"),1000)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 5, nrow = 6)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ARD/all_rates_ARD.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 5, nrow = 6)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ARD/all_rates_ARD_violinplot_median.png", width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Ruminantia") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
dev.off()

#save out as png, remember to change file name
for(i in 1:6){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ARD/rates_plots", i, ".png"), width=30,height=5,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 5, nrow = 1, page = i))
  dev.off()
}

#get the mean rate out of the 1000 model results for each of the 30 rates
#group by solution, calculate the mean
rates_colours_30 <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf", "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/ruminants_six_state_ARD/all_rates_mean_ARD.png", width = 30, height = 15, units = "cm", res = 600)
rates_df %>% group_by(solution) %>% summarise(mean_rates = mean(rates)) %>% ggplot(., aes(x = solution, y = mean_rates)) + geom_col(fill = rates_colours_30) + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95, size = 10)) + ggtitle("Ruminantia") + labs(x = "Transition", y = "Mean rate")
dev.off()

test <- rates_df %>% group_by(solution) %>% summarise(mean_rates = mean(rates)) 

# ##To colour each of the rates by their model likelihood
# #there are six rates for each model (except bridge_only which has 4), so duplicate the likelihood of each model 6 times
# likelihood<- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
# likelihood <- as.data.frame(likelihood)
# index <- rep(1:nrow(likelihood), times = c(rep(x = 6, times = 1000)))
# likelihood <- likelihood[index,]
# 
# #add this column to the rates dataframe
# rates_df <- cbind(rates_df, likelihood)
# rates_df$likelihood <- as.integer(rates_df$likelihood)
# rates_df$likelihood <- as.factor(rates_df$likelihood)
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/new_artiodactyla_max_dinoc_ARD/likelihood_rates_plot.png", width=20,height=12,units="cm",res=200)
# ggplot(rates_df, aes(x= rates, fill = likelihood)) + geom_histogram() + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank()) + facet_wrap(~solution)
# dev.off()

# Section 4.5: Whippomorpha, five_state--------------------------------------

#extract and plot likelihood scores
model_results <- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))

ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", likelihoods = ER_likelihoods)
df2 <- data.frame(model = "SYM", likelihoods = SYM_likelihoods)
df3 <- data.frame(model = "ARD", likelihoods = ARD_likelihoods)
df_s6x <- rbind(df1, df2, df3)

all_six_states_whippo_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "blue") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "Log likelihood") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different")) + ggtitle("Whippomorpha") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/whippomorpha_six_state_all_models_plot.png",  units = "cm", height = 15, width = 15, res = 225)
all_six_states_whippo_plot
dev.off()

#extract and plot AICc scores from 1k models
model_results <- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))

ER_AICc <- unlist(lapply(model_results$ER_model, function(x) returnAICc(model = x)))
SYM_AICc  <- unlist(lapply(model_results$SYM_model, function(x) returnAICc(model = x)))
ARD_AICc  <- unlist(lapply(model_results$ARD_model, function(x) returnAICc(model = x)))
bridge_only_AICc  <- unlist(lapply(model_results$bridge_only_model, function(x) returnAICc(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", AICc = ER_AICc)
df2 <- data.frame(model = "SYM", AICc = SYM_AICc)
df3 <- data.frame(model = "ARD", AICc = ARD_AICc)
df_s6x <- rbind(df1, df2, df3)

six_states_AICc_whippo_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = AICc)) + geom_jitter(alpha = 0.6, color = "blue") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different")) + ggtitle("Whippomorpha") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/whippomorpha_six_state_all_models_AICc.png",  units = "cm", height = 15, width = 15, res = 225)
six_states_AICc_whippo_plot
dev.off()

#extract and plot AIC
model_results <- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))

ER_AIC <- unlist(lapply(model_results$ER_model, function(x) returnAIC(model = x)))
SYM_AIC  <- unlist(lapply(model_results$SYM_model, function(x) returnAIC(model = x)))
ARD_AIC  <- unlist(lapply(model_results$ARD_model, function(x) returnAIC(model = x)))
bridge_only_AIC  <- unlist(lapply(model_results$bridge_only_model, function(x) returnAIC(model = x)))

## combine and plot

df1 <- data.frame(model = "ER", AIC = ER_AIC)
df2 <- data.frame(model = "SYM", AIC = SYM_AIC)
df3 <- data.frame(model = "ARD", AIC = ARD_AIC)
df_s6x <- rbind(df1, df2, df3)

six_states_AIC_whippo_plot <- ggplot(df_s6x, aes(x = fct_inorder(model), y = AIC)) + geom_jitter(alpha = 0.6, color = "blue") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =10), axis.text.y = element_text(size =10)) +
  labs(x = "Model", y = "AIC score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different")) + ggtitle("Whippomorpha") 

#save plot as a PNG, remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/whippo_six_state_all_models_AIC.png",  units = "cm", height = 15, width = 15, res = 225)
six_states_AIC_whippo_plot
dev.off()

# ###Whippomorpha_five_state_ER_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))

rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_ER"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
rates_colours_30 <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#a63d13", "#bd5c35", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ER/all_rates_ER.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 4, nrow = 5)
dev.off() 

for(i in 1:5){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ER/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 4, nrow = 1, page = i)) 
  dev.off()
}

# ###Whippomorpha_six_state_SYM_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))

rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_SYM"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]

rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
rates_colours_30 <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#a63d13", "#bd5c35", "#e79979", "#fbb89d",  "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_SYM/all_rates_SYM.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 4, nrow = 5)
dev.off()

for(i in 1:5){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_SYM/rates_plots", i, ".png"), width=20,height=10,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 4, nrow = 1, page = i))
  dev.off()
}

# ###Whippomorpha_five_state_ARD_rates### ----------------------------------------------------------------

#use below for Amelia primary source data
model_results <- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))

rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#we want to format this into a dataframe with three columns
#model name, rates, solution number
#this will allow us to plot the rates by which model they came from and by which 
rates_df <- as.data.frame(rates)
rates_df$model <- "six_state_ARD"

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$rates)),]
row.names(rates_df) <- 1:length(rates_df$rates)

rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
rates_colours <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#a63d13", "#bd5c35", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")

ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 4, nrow = 5)

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ARD/all_rates_ARD.png", width = 30, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution, ncol = 4, nrow = 5)
dev.off()

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ARD/all_rates_ARD_violinplot_median.png", width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Whippomorpha") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
dev.off()

rates_df <- rates_df %>% mutate(log_rates = log10(rates) + 10)
ggplot(rates_df, aes(x= solution, y = log(log_rates), group = solution, fill = solution, color = solution)) + geom_jitter(aes(color = solution, alpha = 0.2)) + geom_violin(colour = "black")

ggplot(rates_df, aes(x= solution, y = log(rates), group = solution)) + geom_boxplot(outlier.shape = NA) #+ theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Whippomorpha") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 
ggplot(rates_df, aes(x= solution, y = log_rates, group = solution, fill = solution, color = solution)) + geom_jitter(aes(color = solution, alpha = 0.2)) + geom_violin(colour = "black") + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_colours) + theme(legend.position = "none") + ggtitle("Whippomorpha") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=3, colour = "red") 

#save out as png, remember to change file name
for(i in 1:5){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ARD/rates_plots", i, ".png"), width=30,height=5,units="cm",res=200)
  print(ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap_paginate(~solution, ncol = 4, nrow = 1, page = i))
  dev.off()
}

#get the mean rate out of the 1000 model results for each of the 30 rates
#group by solution, calculate the mean
rates_colours <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#a63d13", "#bd5c35", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")

png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ARD/all_rates_mean_ARD.png", width = 30, height = 15, units = "cm", res = 600)
rates_df %>% group_by(solution) %>% summarise(mean_rates = mean(rates)) %>% ggplot(., aes(x = solution, y = mean_rates)) + geom_col(fill = rates_colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95, size = 10)) + ggtitle("Whippomorpha") + labs(x = "Transition", y = "Mean rate") 
dev.off()

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
bridge_only_results <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))
#use below for artiodactyla -primary source data 
#bridge_only_results <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

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

#coloured by likelihood
# scatter_dicath <- ggplot(rates_pivot, aes(x = CathDi, y = DiCath, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "deeppink4", high = "plum") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 
# scatter_noccath <- ggplot(rates_pivot, aes(x = NocCath, y = CathNoc, colour = likelihoods)) + geom_point() + scale_colour_gradient(low = "dodgerblue4", high = "cyan") +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

#Not coloured by likelihood
scatter_dicath <- ggplot(rates_pivot, aes(x = CathDi, y = DiCath)) + geom_point(color = "deeppink3", alpha = 0.5) +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 
scatter_noccath <- ggplot(rates_pivot, aes(x = NocCath, y = CathNoc)) + geom_point(color = "dodgerblue3", alpha = 0.5) +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) 

#to be continued
scatter_1 <- ggplot(rates_pivot, aes(x = NocCath, y = DiCath)) + geom_point(colour = "purple1", alpha = 0.5) +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))
scatter_2 <- ggplot(rates_pivot, aes(x = NocCath, y = CathDi)) + geom_point(colour = "purple2", alpha = 0.5) +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))
scatter_3 <- ggplot(rates_pivot, aes(x = CathNoc, y = DiCath)) + geom_point(colour = "purple3", alpha = 0.5) +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))
scatter_4 <- ggplot(rates_pivot, aes(x = CathNoc, y = CathDi)) + geom_point(colour = "purple4", alpha = 0.5) +  theme(plot.title = element_text(size = 24), axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"))


#save out as a png, remember to change file name

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/cet_max_crep_bridge_only_scatter_dicath_v_cathdi.png", units = "cm", height = 20, width = 20, res = 500)
scatter_dicath
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/cet_max_crep_bridge_only_scatter_noccath_v_cathnoc.png", units = "cm", height = 20, width = 20, res = 500)
scatter_noccath
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/cet_max_crep_bridge_only_scatter_noccath_v_dicath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_1
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/cet_max_crep_bridge_only_scatter_noccath_v_cathdi.png", units = "cm", height = 20, width = 20, res = 500)
scatter_2
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/cet_max_crep_bridge_only_scatter_cathnoc_v_dicath.png", units = "cm", height = 20, width = 20, res = 500)
scatter_3
dev.off()

png(filename = "C:/Users/ameli/OneDrive/Documents/R_projects/Rates_scatterplots/cet_max_crep_bridge_only_scatter_cathnoc_v_cathdi.png", units = "cm", height = 20, width = 20, res = 500)
scatter_4
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





# #Section 9: Plot max_clade_cred tree results------------------------------------
model_results <- readRDS(here("artiodactyla_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_models.rds"))

# ER_results <- model_results$ER_model
# SYM_results <- model_results$SYM_model
ARD_results <- model_results$ARD_model
plotMKmodel(ARD_results)

model_results <- readRDS(here("artiodactyla_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_models.rds"))

# create a new plotting window and set the plotting area into a 1*3 array
par(mfrow = c(2, 2))

# plot a bar chart for max.temp
plotMKmodel(model_results$ER_model)
plotMKmodel(model_results$SYM_model)
plotMKmodel(model_results$ARD_model)

rates_df <- as.data.frame(ARD_results$solution)

#use below to plot rates from 6 state model
rates_df <- pivot_longer(rates_df, cols = 1:ncol(rates_df), names_to = "rates")

#drop rows with no transition rates (ie Di -> Di, Noc -> Noc)
rates_df <- rates_df[!(is.na(rates_df$value)),]
row.names(rates_df) <- 1:length(rates_df$rates)

#add a column indicating what transition rate the value is for
#use below for 6 states
# rates_df$solution <- c("Cath -> Cath/crep", "Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep", "Cath/crep -> Cath", "Cath/crep -> Di", "Cath/crep -> Di/crep", "Cath/crep -> Noc", "Cath/crep -> Noc/crep", "Di -> Cath", "Di -> Cath/crep", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Cath/crep", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Cath/crep", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Cath/crep", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
# rates_colours_30 <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf", "#bd5c35", "#d37a57", "#e79979", "#fbb89d",  "#a63d13", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")

#use below for 5 states
rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
rates_colours_30 <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#bd5c35", "#e79979", "#fbb89d",  "#a63d13", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")

#remember to change file name
png("C:/Users/ameli/OneDrive/Documents/R_projects/rates_dump/whippomorpha_five_state_ARD/all_rates_max_clade_ARD.png", width = 30, height = 15, units = "cm", res = 600)
ggplot(rates_df, aes(x = solution, y = value)) + geom_col(fill = rates_colours_30) + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95, size = 8))
dev.off()

model_results <- readRDS(here("artiodactyla_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_models.rds"))

log_likelihoods <- unlist(lapply(model_results, function(x) returnLikelihoods(model = x)))
likelihoods <- as.data.frame(log_likelihoods)
likelihoods$model <- rownames(likelihoods)
likelihoods$AICc_scores <-AICc_scores 
likelihoods$AIC_scores <-AIC_scores 

likelihood_plot <- ggplot(likelihoods, aes(x = model, y = log_likelihoods)) + geom_point()
AICc_scores <- unlist(lapply(model_results, function(x) returnAICc(model = x)))
AICc_plot <- ggplot(likelihoods, aes(x = model, y = AICc_scores)) + geom_point()
AIC_scores <- unlist(lapply(model_results, function(x) returnAIC(model = x)))

likelihoods <- likelihoods %>% pivot_longer(!model, names_to = "model_metric", values_to = "model_value")

ggplot(likelihoods, aes(x = model, y = model_value, color = model)) + geom_point(size = 5) + scale_color_brewer(palette = "Accent")+ facet_wrap(~model_metric, scales = "free")

png("C:/Users/ameli/OneDrive/Documents/R_projects/max_clade_cred_likelihoods/likelihood_metrics.png", width = 30, height = 15, units = "cm", res = 600)
ggplot(likelihoods, aes(x = model, y = model_value, color = model)) + geom_point(size = 5) + scale_color_brewer(palette = "Accent")+ facet_wrap(~model_metric, scales = "free")
dev.off()

ARD_matrix <- as.matrix(ARD_results$data)
make.simmap(ARD_results$phy, as.matrix(ARD_results$data), model = ARD)

# #Section 10: Statistical significance -----------------------------------

#we want to see if the difference between the mean likelihood of one model is greater than another

#ie is the ARD model statistically more likely than the bridge_only model

# Compute the analysis of variance
res.aov <- aov(weight ~ group, data = my_data)
# Summary of the analysis
summary(res.aov)
