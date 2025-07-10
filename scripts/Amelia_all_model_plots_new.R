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
#for faceting plots over multiple pages

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
install.packages("ggdist")
library(ggdist)
library(knitr)
#install.packages("kableExtra")
library(kableExtra)
#install.packages("webshot")
library(webshot)

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
filename <- "cetaceans_fixed_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

#returns a dataframe of all three metrics for all models
likelihood_metrics <- max_clade_metrics(readRDS(here(filename)))
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)

#generate and save out plots
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/max_clade_cred_likelihoods/likelihood_metrics_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)
ggplot(likelihood_metrics, aes(x = fct_inorder(model), y = model_value, color = model)) + geom_point(size = 5) + scale_color_brewer(palette = "Accent") + theme(legend.position = "none") + labs(x = "model", y = "metric value") + facet_wrap(~model_metric, scales = "free")
dev.off()

knitr::kable(likelihood_metrics, format = "html", digits = 2, caption = filename) %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("likelihood_table.html")
#webshot::install_phantomjs(force = TRUE)
webshot("likelihood_table.html", file = paste0("C:/Users/ameli/OneDrive/Documents/R_projects/max_clade_cred_likelihoods/likelihood_table_", filename, ".png"))

# #Section 2: max_clade_cred rates ---------------------------

#set file name
filename <- "fixed_whippomorpha_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_CONSYM_models.rds"

model_results <- readRDS(here(filename))

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/max_clade_cred_likelihoods/max_clade_cred_rates_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)

# create a new plotting window and set the plotting area into a 2*2 array
par(mfrow = c(2, 2))

# plot a bar chart for max.temp
plotMKmodel(model_results$ER_model)
plotMKmodel(model_results$SYM_model)
plotMKmodel(model_results$ARD_model)
plotMKmodel(model_results$bridge_only)
plotMKmodel(model_results$CONSYM_model)
dev.off()

# #Section 3: Plot likelihoods from 1k model results ----------------------
model_6 <- readRDS(here(paste0("finalized_1k_models/", "fixed_cetaceans_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models_1-5.rds")))
model_10 <- readRDS(here(paste0("finalized_1k_models/", "fixed_cetaceans_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models_6-10.rds")))
model_20 <- readRDS(here(paste0("finalized_1k_models/", "fixed_cetaceans_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models_10-20.rds")))
model_results <- c(model_6, model_10, model_20)
  
filename <- "fixed_cetaceans_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models_10-20.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#function returns a dataframe of the likelihoods for all 1k trees x number of Mk models
df_full <- plot1kLikelihoods(readRDS(here(paste0("finalized_1k_models/", filename))), 5)

df_full <- plot1kLikelihoods(model_results, 5)

#check if data fits the ANOVA assumptions
#anovaAssumptions(df_full, df_full$likelihoods)

#determine the number of comparisons
if(length(unique(df_full$model)) == 5){
  my_comparisons <- list( c("ER", "SYM"), c("ER", "ARD"), c("ER", "bridge_only"), c("SYM", "ARD"), c("SYM", "bridge_only"), c("ARD", "bridge_only"), c("ER", "CONSYM"), c("SYM", "CONSYM"), c("ARD", "CONSYM"), c("bridge_only", "CONSYM"))
} 

if(length(unique(df_full$model)) == 4){
  my_comparisons <- list( c("ER", "SYM"), c("ER", "ARD"), c("ER", "bridge_only"), c("SYM", "ARD"), c("SYM", "bridge_only"), c("ARD", "bridge_only"))
} 

if(length(unique(df_full$model)) == 3){
  my_comparisons <- list( c("ER", "SYM"), c("ER", "ARD"), c("SYM", "ARD"))
}

means <- aggregate(likelihoods ~  model, df_full, mean)
means$likelihoods <- round(means$likelihoods, digits = 2)

#plot and save out
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, " _likelihoods", ".png"), width = 20, height = 20, units = "cm", res = 400)
ggplot(df_full, aes(x = fct_inorder(model), y = likelihoods)) + 
  geom_jitter(alpha = 0.6, color = "#F8766D") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "Log likelihood") + 
  scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "Constrained ARD", "Constrained SYM")) + 
  geom_text(data = means, aes(label = likelihoods, y = likelihoods, vjust = -0.5), parse = TRUE) + 
  ggtitle(filename)  
dev.off()

# #Section 4: Plot AIC scores from 1k model results ----------------------
#filename <- "fixed_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_CONSYM_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#returns a df of the AIC scores for all 1k trees x number of Mk models
df_full <- plot1kAIC(readRDS(here(paste0("finalized_1k_models/", filename))), 5)

df_full <- plot1kAIC(model_results, 5)

means <- aggregate(AIC_score ~  model, df_full, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

#plot and save out - boxplot
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_AIC", ".png"), width = 20, height = 20, units = "cm", res = 400)
ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score)) + geom_jitter(alpha = 0.6, color = "royalblue1") + #other colour option 766df8 and 6daaf8
  geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size = 10), axis.title = element_text(size = 12)) +
  labs(x = "Model", y = "AIC score") +
  scale_x_discrete(labels = c("ER", "SYM", "ARD", "CON-ARD", "CON-SYM")) + 
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, vjust = -0.5), parse = TRUE) +
  ggtitle(filename) 
dev.off()

means <- aggregate(AIC_score ~  model, df_full, mean)
means$AIC_score <- as.integer(means$AIC_score)

#plot and save out - raincloud plot, ruminant format
filename1 <- "fixed_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_CONSYM_models.rds"
df_full1 <- plot1kAIC(readRDS(here(paste0("finalized_1k_models/", filename1))), 5)
means <- aggregate(AIC_score ~  model, df_full1, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", "1_test_", filename1, "_AIC", ".png"), width = 20, height = 15, units = "cm", res = 600)
ggplot(df_full1, aes(x = fct_inorder(model), y = AIC_score, fill = fct_inorder(model))) + 
  scale_fill_manual(values = c("blue", "royalblue","slateblue", "mediumpurple", "orchid")) +
  ggdist::stat_halfeye(alpha = 0.6, adjust = .5, width = .6, justification = -.3, .width = 0, point_colour = NA) +
  geom_boxplot(alpha = 0.2, width = .25, colour = "black",outlier.shape = NA) + 
  geom_point(aes(color = fct_inorder(model)), stroke = 1, size = 1, alpha = .2, position = position_jitter(seed = 1, width = .15)) +
  scale_color_manual(values = c("blue", "royalblue","slateblue", "mediumpurple", "orchid")) +
  labs(x = "Model", y = "AIC scores")  +
  scale_x_discrete(labels = c("ER", "SYM", "ARD", "CON-ARD")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18), legend.position = "none") +
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, hjust = -1), parse = TRUE) +
  ggtitle(filename1) #+ 
  #coord_cartesian(xlim = c(1.4, 4.3), ylim = c(342, 451))
dev.off()

#plot and save out - raincloud plot, whippo format
filename2 <- "fixed_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_CONSYM_models.rds"
df_full2 <- plot1kAIC(readRDS(here(paste0("finalized_1k_models/", filename2))), 4)
means <- aggregate(AIC_score ~  model, df_full2, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", "1_test_", filename2, "_AIC", ".png"), width = 20, height = 15, units = "cm", res = 600)
ggplot(df_full2, aes(x = fct_inorder(model), y = AIC_score, fill = fct_inorder(model))) + 
  scale_fill_manual(values = c("blue", "royalblue","slateblue", "mediumpurple")) +
  ggdist::stat_halfeye(alpha = 0.6, adjust = .5, width = 0.6, justification = -.3, .width = 0, point_colour = NA) +
  geom_boxplot(alpha = 0.2, width = .25, colour = "black",outlier.shape = NA) + 
  geom_point(aes(color = fct_inorder(model)), stroke = 1, size = 1, alpha = .2, position = position_jitter(seed = 1, width = .15)) +
  scale_color_manual(values = c("blue", "royalblue","slateblue", "mediumpurple")) +
  labs(x = "Model", y = "AIC scores")  +
  scale_x_discrete(labels = c("ER", "SYM", "ARD", "CON-ARD")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18), legend.position = "none") +
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, hjust = -1)) +
  ggtitle(filename2) + 
  coord_cartesian(xlim = c(1.4, 4.0), ylim = c(202, 233))
dev.off()

ggplot(df_full2, aes(y = AIC_score, x = model))+ geom_boxplot() + stat_summary(fun = "mean", geom="point", size=2, colour = "red")

# #Section 5: Plot AICc scores from 1k model results ----------------------
#filename <- "finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#returns a df of the AICc scores for all 1k trees x number of Mk models
df_full <- plot1kAICc(readRDS(here(paste0("finalized_1k_models/", filename))), 5)

means <- aggregate(AICc_score ~  model, df_full, mean)
means$AICc_score <- round(means$AICc_score, digits = 2)

#plot and save out
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_AICc", ".png"), width = 20, height = 20, units = "cm", res = 400)
ggplot(df_full, aes(x = fct_inorder(model), y = AICc_score)) + geom_jitter(alpha = 0.6, color = "#6daaf8") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "CON-ARD", "CON-SYM")) + 
  geom_text(data = means, aes(label = AICc_score, y = AICc_score, vjust = -0.5), parse = TRUE) +
  ggtitle(filename) 
dev.off()

# #Section 6: Plot transition rates from 1k model results ----------------------
#filename <- "artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename, the number of states in the model and the number of Mk models 
#returns a dataframe of the rates from each of the Mk models, for each of the 1k trees
#filename <- "ruminants_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"
#reminder: enter 5 states for the 6 state cetacean/whippomorpha
rates_df <- plot1kTransitionRates(readRDS(here(paste0("finalized_1k_models/", filename))), 4, 4)
rates_df <- plot1kTransitionRates(model_results, 4, 4)

#plot as a violin plot
#important note: the colours column doesn't line up with the correct solution in the df but if we plot the solutions in alphabetical order and then colouring them with the pallette in the colours column is in the correct order
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_violin_rate_plot_all", ".png"), width = 40, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df$colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_df$colours) + theme(legend.position = "none") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=2, colour = "red") + ggtitle(filename) + facet_wrap(~fct_inorder(model), ncol = 2, nrow = 2) 
dev.off()

#violin plots for each Mk model
for(i in 1:length(unique(rates_df$model))){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_violin_rate_plot_", unique(rates_df$model)[i], ".png"), width = 25, height = 12, units = "cm", res = 600)
  print(ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df$colours) + geom_violin(color = "black", scale = "width") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_df$colours) + theme(legend.position = "none") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=2, colour = "red") + ggtitle(filename) + facet_wrap_paginate(~fct_inorder(model), ncol = 1, nrow = 1, page = i)) 
  dev.off()
}

#plot as full rates histograms
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_histogram_rate_plot_all",".png"), width = 50, height = 20, units = "cm", res = 600)
ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_df$colours) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + labs(title = filename, y = "count" , x = "transition rate per unit of evolutionary time") + facet_wrap(~fct_inorder(model) + solution, nrow = length(unique(rates_df$model)), ncol = length(unique(rates_df$solution)))
dev.off()

#plot rates histograms separated by Mk model
#cannot find a way to facet_wrap_paginate by multiple variables, instead I will cycle through all the models and plot them with a for loop

for(i in 1:(length(unique(rates_df$model)))){
  png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_rate_plot_", unique(rates_df$model)[i], ".png"), width = 50, height = 20, units = "cm", res = 600)
  print(rates_df %>% filter(model == unique(rates_df$model)[i]) %>% ggplot(., aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_df$colours) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + labs(title = filename, y = "count" , x = "transition rate per unit of evolutionary time") + facet_wrap(~solution))
  dev.off()
}



