setwd(here())
source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")
source("scripts/Amelia_plotting_functions.R")

#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

# Section 1: max_clade_cred likelihood metrics ---------------------------

#set file name
filename <- "artiodactyla_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
#filename <- "whippomorpha_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
#filename <- "ruminants_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"

#filename <- "artiodactyla_finalized_max_clade_cred_six_state_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
#filename <- "whippomorpha_finalized_max_clade_cred_six_state_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
#filename <- "ruminants_finalized_max_clade_cred_six_state_traits_ER_SYM_CONSYM_ARD_bridge_only_models"

#returns a dataframe of all three metrics for all models
likelihood_metrics <- max_clade_metrics(readRDS(here(paste0(filename, ".rds"))))
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)
likelihood_metrics$most_likely <- ""  
likelihood_metrics[which(likelihood_metrics$AIC_scores == min(likelihood_metrics$AIC_scores)), "most_likely"] <- "**"

knitr::kable(likelihood_metrics, format = "html", digits = 2, caption = filename) %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("likelihood_table.html")
webshot("likelihood_table.html", file = paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/likelihood_table_", filename, ".png"), vwidth = 992, vheight = 300)

# Section 2: max_clade_cred rates ---------------------------

model_results <- readRDS(here(paste0(filename, ".rds")))

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_clade_cred_rates_", filename, ".png"), width = 30, height = 15, units = "cm", res = 600)

# create a new plotting window and set the plotting area into a 2*2 array
par(mfrow = c(2, 3))
plotMKmodel(model_results$ER_model)
plotMKmodel(model_results$SYM_model)
plotMKmodel(model_results$ARD_model)
plotMKmodel(model_results$bridge_only)
plotMKmodel(model_results$CONSYM_model)
dev.off()


# Section 3: Load in and format data --------------------------------------
model_1_500 <- readRDS(here("august_12__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
model_501_1000 <- readRDS(here("august_14__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
model_results <- c(model_1_500, model_501_1000)
saveRDS(model_results, here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

# Section 4: Plot AIC scores from 1k model results ----------------------
#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#returns a df of the AIC scores for all 1k trees x number of Mk models
df_full <- plot1kAIC(readRDS(here(filename)), 5)

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
df_full <- plot1kAIC(readRDS(here(filename)), 5)
means <- aggregate(AIC_score ~  model, df_full1, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", "1_test_", filename1, "_AIC", ".png"), width = 20, height = 15, units = "cm", res = 600)
ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score, fill = fct_inorder(model))) + 
  scale_fill_manual(values = c("blue", "royalblue","slateblue", "mediumpurple", "orchid")) +
  ggdist::stat_halfeye(alpha = 0.6, adjust = .5, width = .6, justification = -.3, .width = 0, point_colour = NA) +
  geom_boxplot(alpha = 0.2, width = .25, colour = "black",outlier.shape = NA) + 
  geom_point(aes(color = fct_inorder(model)), stroke = 1, size = 1, alpha = .2, position = position_jitter(seed = 1, width = .15)) +
  scale_color_manual(values = c("blue", "royalblue","slateblue", "mediumpurple", "orchid")) +
  labs(x = "Model", y = "AIC scores")  +
  scale_x_discrete(labels = c("ER", "SYM", "ARD", "CON-ARD", "CON-SYM")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18), legend.position = "none") +
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, hjust = -1), parse = TRUE) +
  ggtitle(filename) #+ 
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

# Section 5: Plot AICc scores from 1k model results ----------------------
#filename <- "finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, bridge_only)
#returns a df of the AICc scores for all 1k trees x number of Mk models
df_full <- plot1kAICc(readRDS(here(paste0("finalized_1k_models/", filename))), 5)
df_full <- plot1kAICc(model_results, 5)

means <- aggregate(AICc_score ~  model, df_full, mean)
means$AICc_score <- round(means$AICc_score, digits = 2)

#plot and save out
png(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/finalized_results_plots/", filename, "_AICc", ".png"), width = 20, height = 20, units = "cm", res = 400)
ggplot(df_full, aes(x = fct_inorder(model), y = AICc_score)) + geom_jitter(alpha = 0.6, color = "#6daaf8") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "CON-ARD", "CON-SYM")) + 
  geom_text(data = means, aes(label = AICc_score, y = AICc_score, vjust = -0.5), parse = TRUE) +
  ggtitle(filename) 
dev.off()

# Section 6: Plot transition rates from 1k model results ----------------------
#filename <- "artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"

#requires the filename, the number of states in the model and the number of Mk models 
#returns a dataframe of the rates from each of the Mk models, for each of the 1k trees
#filename <- "ruminants_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"
#reminder: enter 5 states for the 6 state cetacean/whippomorpha
rates_df <- plot1kTransitionRates(model_results, 4, 5)

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





# Section 7: Ancestral reconstruction -----------------------------------

#load in model file

#filename <- "artiodactyla_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
filename <- "whippomorpha_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
#filename <- "ruminants_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"

all_model_results <- readRDS(here(paste0(filename, ".rds")))


#separate the results by the model types we want to use (ER, SYM, ARD, bridge_only)
#uncomment the model you want to plot

# model_results <- all_model_results$ER_model
# model_name <- "ER"

# model_results <- all_model_results$SYM_model
# model_name <- "SYM"

model_results <- all_model_results$CONSYM_model
model_name <- "CONSYM"

# model_results <- all_model_results$ARD_model
# model_name <- "ARD"

# model_results <- all_model_results$bridge_only
# model_name <- "bridge_only"

#from the model results file, tip states describes the trait states at the tips, states describes the trait states at the nodes
lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
phylo_tree <- model_results$phy
#associate each of these species and their trait states with its node
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

#plot the ancestral reconstruction, displaying each of the three trait states (cathemeral, diurnal, nocturnal)
ancestral_plot_di <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_colour_gradientn(colours = c("white", "#BF491B"))  + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
ancestral_plot_di <- ancestral_plot_di + theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_noc <- ggtree(phylo_tree, layout = "circular", size =2) %<+% lik.anc + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_colour_gradientn(colours = c("white", "#297D63")) + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
ancestral_plot_noc <- ancestral_plot_noc + theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_cath <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_colour_gradientn(colours = c("white", "#A747B3")) + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
ancestral_plot_cath <- ancestral_plot_cath + theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_crep <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = crepuscular) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5) + scale_colour_gradientn(colours = c("white", "#856A54")) + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
ancestral_plot_crep <- ancestral_plot_crep + theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_", file_name, "_", model_name, ".pdf", sep = ""), width=18,height=17, bg = "transparent")
grid.arrange(ancestral_plot_di,ancestral_plot_noc, ancestral_plot_cath, ancestral_plot_crep, ncol = 2, nrow = 2)
dev.off()

ancestral_plot <- ggtree(phylo_tree, layout = "circular", size = 2) + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_text(aes(label=node, colour = "red"), hjust=-.2, size = 3)
ancestral_plot
#510 is the LCA of whippo, 511 is the LCA of cetaceans, 308 LC of ruminants in the artiodactyla tree
lik.anc %>% filter(node %in% c(510, 511, 308))

#in the whippomorpha tree the LCA node is 70, LCA of cetaceans is 79
lik.anc %>% filter(node %in% c(78, 79))

#in the ruminant tree the LCA node is 204
lik.anc %>% filter(node %in% c(204))

trait.data <- read.csv(here("ruminants_full.csv"))
#trait.data <- read.csv(here("whippomorpha.csv"))
trait.data <- trait.data[!is.na(trait.data$max_crep), c("tips", "max_crep")]
phylo_trees <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#subset trait data to only include species that are in the tree
trait.data <- trait.data[trait.data$tips %in% phylo_trees$tip.label,]
# this selects a tree that is only the subset with data (mutual exclusive)
phylo_trees <- keep.tip(phylo_trees, tip = trait.data$tips)
# bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), node.states = "marginal")
# bridge_only_crep_root <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), node.states = "marginal", root.p = c(0.19,0.48,0.18,0.15))

bridge_only_crep_root2 <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), node.states = "marginal", root.p = c(0,1,0,0))
bridge_only_di_root <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), node.states = "marginal", root.p = c(0,0,1,0))
bridge_only_cath_root <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), node.states = "marginal", root.p = c(1,0,0,0))
bridge_only_noc_root <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), node.states = "marginal", root.p = c(0,0,0,1))

#test_list <- list(bridge_only, bridge_only_crep_root, bridge_only_crep_root2)
#names(test_list) <- c("yang_root", "crep_root_50", "crep_root_100")
test_list <- list(bridge_only_crep_root2, bridge_only_di_root, bridge_only_cath_root, bridge_only_noc_root)
names(test_list) <- c("crep_root", "di_root", "cath_root", "noc_root")
likelihood_metrics <- max_clade_metrics(test_list)
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)
likelihood_metrics$most_likely <- ""  
likelihood_metrics[which(likelihood_metrics$AIC_scores == min(likelihood_metrics$AIC_scores)), "most_likely"] <- "**"

#create the name of the file by pasting together ancestral recon, the diel state and the file_name 
# pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_diurnal_", file_name, "_", model_name, ".pdf", sep = ""), width=17,height=16)
# ancestral_plot_di
# dev.off()
# 
# pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_nocturnal_", file_name, "_", model_name, ".pdf", sep = ""), width=17,height=16)
# ancestral_plot_noc
# dev.off()
# 
# pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_cathemeral_", file_name, "_", model_name,  ".pdf", sep = ""), width=17,height=16)
# ancestral_plot_cath
# dev.off()
# 
# pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_crepuscular_", file_name, "_", model_name,  ".pdf", sep = ""), width=17,height=16)
# ancestral_plot_crep
# dev.off()



# Section 10: Compare old and new results ---------------------------------

filename <- "whippomorpha_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"
model_results <- readRDS(here(paste0("finalized_1k_models/", filename)))
df_full <- plot1kAIC(readRDS(here(paste0("finalized_1k_models/", filename))), 4)

means <- aggregate(AIC_score ~  model, df_full, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

#plot 
ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score)) + geom_jitter(alpha = 0.6, color = "#6daaf8") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AIC score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "CON-ARD", "CON-SYM")) + 
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, vjust = -0.5), parse = TRUE) +
  ggtitle(filename) 

model_1_500 <- readRDS(here("august_12__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
model_501_1000 <- readRDS(here("august_14__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
model_results <- c(model_1_500, model_501_1000)

saveRDS(model_results, here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))

filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
model_results <- readRDS(here(filename))
df_full <- plot1kAIC(readRDS(here(filename)), 5)

df_full1 <- plot1kAIC(model_1_500, 5)
df_full2 <- plot1kAIC(model_501_1000, 5)
df_full <- plot1kAIC(model_results, 5)

df_full <- rbind(df_full1, df_full2)
means <- aggregate(AIC_score ~  model, df_full, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

#plot 
ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score)) + geom_jitter(alpha = 0.6, color = "#6daaf8") + geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95)) +
  labs(x = "Model", y = "AICc score") + scale_x_discrete(labels = c("Equal rates", "Symmetrical rates", "All rates different", "CON-ARD", "CON-SYM")) + 
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, vjust = -0.5), parse = TRUE) +
  ggtitle(filename) 

