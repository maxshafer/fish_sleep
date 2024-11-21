#Import and run packages -------------------------------------------------
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
source("scripts/Amelia_plotting_functions.R")

# Set the working directory and source the functions (not used yet)
setwd("C:/Users/ameli/OneDrive/Documents/R_projects/fish_sleep/1k_model_results")


#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))


# #Section 1: max_clade_cred likelihood metrics ---------------------------


#set file name
filename <- "testC_whippomorpha_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#returns a dataframe of all three metrics for all models
likelihood_metrics <- max_clade_metrics(readRDS(here(filename)))

#generate and save out plots
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/max_clade_cred_likelihoods/likelihood_metrics_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)
ggplot(likelihood_metrics, aes(x = model, y = model_value, color = model)) + geom_point(size = 5) + scale_color_brewer(palette = "Accent") + theme(legend.position = "none")  + facet_wrap(~model_metric, scales = "free")
dev.off()

# #Section 2: max_clade_cred rates ---------------------------

#set file name
filename <- "artiodactyla_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_models.rds"

model_results <- readRDS(here(filename))

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/max_clade_cred_likelihoods/max_clade_cred_rates_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)

# create a new plotting window and set the plotting area into a 2*2 array
par(mfrow = c(2, 2))

# plot a bar chart for max.temp
plotMKmodel(model_results$ER_model)
plotMKmodel(model_results$SYM_model)
plotMKmodel(model_results$ARD_model)

dev.off()

# #Section 3: Plot likelihoods from 1k model results ----------------------
filename <- "test2_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#function returns a dataframe of the likelihoods for all 1k trees x number of Mk models
df_full <- plot1kLikelihoods(readRDS(here(filename)), 4)

#plot and save out
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/likelihoods_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)
ggplot(df_full, aes(x = fct_inorder(model), y = likelihoods)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))
dev.off()

# #Section 4: Plot AIC scores from 1k model results ----------------------
#filename <- "finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#returns a df of the AIC scores for all 1k trees x number of Mk models
df_full <- plot1kAIC(readRDS(here(filename)), 4)

#plot and save out
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/likelihoods_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)
ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AIC scores") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))

dev.off()

# #Section 5: Plot AICc scores from 1k model results ----------------------
#filename <- "finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#returns a df of the AICc scores for all 1k trees x number of Mk models
df_full <- plot1kAICc(readRDS(here(filename)), 4)

#plot and save out
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/1k_trees/likelihoods_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)
ggplot(df_full, aes(x = fct_inorder(model), y = AICc_score)) + geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Bridge only"))
dev.off()

# #Section 6: Plot transition rates from 1k model results ----------------------
#filename <- "finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename, number of Mk models and the number of states in the model
#returns a dataframe of the rates from each of the models, for each of the 
rates_df <- plot1kTransitionRates(readRDS(here(filename)), 4)

#plot as a violin plot

#plot as full histograms
ER_rates_plot 
ARD_rates_plot <- ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_colours)+ scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + facet_wrap(~solution)

#plot as individual histograms