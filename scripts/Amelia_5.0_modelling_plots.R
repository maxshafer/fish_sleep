
setwd(here())
source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")
source("scripts/Amelia_plotting_functions.R")

#load in mammal tree and cetacean dataframe
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

# Section 1: max_clade_cred likelihood metrics ---------------------------

#set file name
#filename <- "artiodactyla_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
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
model_results1 <- readRDS(here("august_12__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
model_results2 <- readRDS(here("august_14__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
model_results1$ER_model <- c(model_results1$ER_model, model_results2$ER_model)
model_results1$SYM_model <- c(model_results1$SYM_model, model_results2$SYM_model)
model_results1$ARD_model <- c(model_results1$ARD_model, model_results2$ARD_model)
model_results1$bridge_only_model <- c(model_results1$bridge_only_model, model_results2$bridge_only_model)
model_results1$CONSYM_model <- c(model_results1$CONSYM_model, model_results2$CONSYM_model)

saveRDS(model_results1, here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
# Section 4: Plot AIC scores from 1k model results ----------------------

filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_artiodactyla_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, CONARD, 5: ER, SYM, ARD, bridge_only, CONSYM)
#returns a df of the AIC scores for all 1k trees x number of Mk models
df_full <- plot1kAIC(readRDS(here(filename)), 5)

means <- aggregate(AIC_score ~  model, df_full, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

#plot boxplot
ggplot(df_full, aes(x = fct_inorder(model), y = AIC_score)) + geom_jitter(alpha = 0.6, color = "royalblue1") + #other colour option 766df8 and 6daaf8
  geom_boxplot(alpha = 0.5, outlier.shape = NA, colour = "black") + theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0, size = 10), axis.title = element_text(size = 12)) +
  labs(x = "Model", y = "AIC score") +
  scale_x_discrete(labels = c("ER", "SYM", "ARD", "CON-ARD", "CON-SYM")) + 
  geom_text(data = means, aes(label = AIC_score, y = AIC_score, vjust = -0.5), parse = TRUE) +
  ggtitle(filename) 

custom.colours <- c("#013873", "#04549F", "#056CCC", "#60B0FB", "#8DC6FC")

df_full$model <- factor(df_full$model, levels = c("ER", "SYM", "CONSYM", "ARD", "bridge_only"))
#plot and save out - raincloud plot
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "raincloud", filename, "_AIC", ".pdf"), width = 7, height = 5)
ggplot(df_full, aes(x = model, y = AIC_score, fill = model)) + 
  scale_fill_manual(values = custom.colours) +
  ggdist::stat_halfeye(alpha = 0.6, adjust = .5, width = .6, justification = -.3, .width = 0, point_colour = NA) +
  geom_boxplot(alpha = 0.2, width = .25, colour = "black",outlier.shape = NA) + 
  geom_point(aes(color = model), stroke = 1, size = 1, alpha = .2, position = position_jitter(seed = 1, width = .15)) +
  scale_color_manual(values = custom.colours) + geom_text(data = means, aes(label = AIC_score, y = AIC_score), hjust = -0.45, size = 5, parse = TRUE) +
  labs(x = "Model", y = "AIC score")  + theme_bw() +
  scale_x_discrete(labels = c("ER", "SYM", "CON-SYM", "ARD", "CON-ARD")) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 16), legend.position = "none") #+ coord_cartesian(ylim = c(398, 525)) #there is one extreme outlier in the ER models of 564
dev.off()

#what is the most likely model on the most likely tree?
min(df_full$AIC_score)

# Section 6: Plot transition rates from 1k model results ----------------------
#requires the filename, the number of states in the model and the number of Mk models 
#returns a dataframe of the rates from each of the Mk models, for each of the 1k trees
#rates_df <- plot1kTransitionRates(readRDS(here(filename)), 4, 5)

rates_df <- plot1kTransitionRates4state(readRDS(here(filename)), 5)

#plot as a violin plot
#important note: the colours column doesn't line up with the correct solution in the df but if we plot the solutions in alphabetical order and then colouring them with the palette in the colours column is in the correct order
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df$colours) + geom_violin(color = "black", scale = "width") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_df$colours) + theme(legend.position = "none") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=2, colour = "red") + ggtitle(filename) + facet_wrap(~fct_inorder(model), ncol = 3, nrow = 2) 

model_selection <- "SYM"
#pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", filename, "_violin_rate_plot_", model_selection, ".pdf"), width = 10, height = 7)
rates_df1 <- rates_df %>% filter(model == model_selection) 
ggplot(rates_df1, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + 
  geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df1$colours) +
  geom_violin(color = "black", scale = "width") +theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  +
  scale_fill_manual(values = rates_df1$colours) + theme(legend.position = "none") +
  labs(x = "Transition", y = "Log(transition rate)") + stat_summary(fun=median, geom="point", size=2, colour = "red") +
  ggtitle(filename) 
#dev.off()

#extract just the starting state
rates_df1$start_state <- word(rates_df1$solution, 1)
rates_df1$start_state <- paste(rates_df1$start_state, "to", sep = " ")
rates_df1$end_state <- word(rates_df1$solution, 3)

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", filename, "_separated_violin_rate_plot_", model_selection, ".pdf"), width = 8, height = 4)
ggplot(rates_df1, aes(x= end_state, y = log(rates), group = solution, fill = solution, colour = solution)) + 
  geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df1$colours) +
  geom_violin(color = "black", scale = "width") +theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =10), axis.text.y = element_text(size =10))  +
  scale_fill_manual(values = rates_df1$colours) + theme(legend.position = "none") + labs(x = "Transition", y = "Log(rates)") + 
  stat_summary(fun=median, geom="point", size=2, colour = "red") + facet_wrap(~start_state, scales = "free_x", nrow = 1, ncol = 4)
dev.off()

rates_df1$start_state <- factor(rates_df1$start_state, levels = c("Cathemeral to", "Diurnal to", "Crepuscular to", "Nocturnal to"))

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", filename, "_square_violin_rate_plot_", model_selection, ".pdf"), width = 7, height = 7)
ggplot(rates_df1, aes(x= end_state, y = log(rates), group = solution, fill = solution, colour = solution)) + 
  geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df1$colours) +
  geom_violin(color = "black", scale = "width", alpha = 0.8) + theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10), axis.text.y = element_text(size =10), axis.title = element_text(size = 16), strip.background = element_rect(fill = "grey95"))  +
  scale_fill_manual(values = rates_df1$colours) + theme(legend.position = "none") + labs(x = "Transition", y = "Log(transition rate)") + 
  stat_summary(fun=median, geom="point", size=2, colour = "firebrick3") + facet_wrap(~start_state, scales = "free_x", nrow = 2, ncol = 2)
dev.off()

#histogram of rates
#ggplot(rates_df, aes(x= rates, fill = solution)) + geom_histogram() + scale_fill_manual(values = rates_df$colours) + scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + theme_bw() + theme(plot.title = element_blank(), legend.position = "none") + labs(title = filename, y = "count" , x = "transition rate per unit of evolutionary time") + facet_wrap(~fct_inorder(model) + solution, nrow = length(unique(rates_df$model)), ncol = length(unique(rates_df$solution)))

#scatterplot of 
ggplot(rates_df1, aes(x = start_state, y = log(rates), colour = end_state)) + geom_boxplot(outlier.shape = NA) 
#does starting in one state give you a higher rate than starting in another
#plotting all the points doesn't show which model they're from, could calculate


# Section 7: Plot transition rates vs proportion of crepuscular sps -------
filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
rates_df <- plot1kTransitionRates4state(readRDS(here(filename)), 5)
rates_df$start_state <- word(rates_df$solution, 1)
rates_df$end_state <- word(rates_df$solution, 3)
rates_df  <- filter(rates_df, model == "ARD")
#rates_df  <- filter(rates_df, model == "Bridge_only")

#every model has 12 rates, constrained models have 10 rates. Every 12 rows is a new model
rates_df$tree_n <- rep(1:1000, each = 12)
ggplot(rates_df, aes(x = rates, y = solution, colour = tree_n)) + geom_jitter()

#for each model, what is the mean rate, how much do all the transition rates differ from that mean
#can allow us to see if cath and crep tend to have faster rates (positive numbers) vs di and noc as starting states (negative numbers)
mean_rates_df <- rates_df %>% group_by(tree_n) %>% summarize(mean_rates = mean(rates))
rates_df <- merge(mean_rates_df,rates_df, all = TRUE)
rates_df <- rates_df %>% mutate(rate_difference = rates - mean_rates)

#transitions out of diurnality tend to occur the fastest for whippo
ggplot(rates_df, aes(x = start_state, y = rate_difference)) + geom_boxplot()
ggplot(rates_df, aes(y = start_state, x = rate_difference)) + geom_violin() +  
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) 

#transitions into cathemerality tend to occur the fastest for whippo
ggplot(rates_df, aes(x = end_state, y = rate_difference)) + geom_boxplot() 

dat <- rates_df %>% group_by(start_state) %>% summarize(mean_rate_difference = mean(rate_difference))
dat2<- rates_df %>% group_by(end_state) %>% summarize(mean_rate_difference = mean(rate_difference))

#plot it like an odds ratio
ggplot(dat, aes(y = start_state, x = mean_rate_difference)) + geom_point(shape = 18, size = 5) +  
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) 

ggplot(dat2, aes(y = end_state, x = mean_rate_difference)) + geom_point(shape = 18, size = 5) +  
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) 

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

#option 2: use the most likely model + tree from 1k trees
#filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
filename <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

model_results <- readRDS(here(filename))
model_results <- model_results$bridge_only_model
#model_results <- model_results[474] #whippomorpha
model_results <- model_results[911] #ruminants
model_results <- model_results$UNTITLED
model_name <- "most_likely_bridge_only"

#from the model results file, tip states describes the trait states at the tips, states describes the trait states at the nodes
lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
phylo_tree <- model_results$phy

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



# Section 8: Compare old and new results ---------------------------------

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


# Section 9: Comparison of rate magnitude
filename1 <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
filename2 <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, CONARD, 5: ER, SYM, ARD, bridge_only, CONSYM)
#returns a df of the AIC scores for all 1k trees x number of Mk models
rates_df1 <- plot1kTransitionRates4state(readRDS(here(filename1)), 5)
rates_df1 <- filter(rates_df1, model == "ARD")
rates_df2 <- plot1kTransitionRates4state(readRDS(here(filename2)), 5)
rates_df2 <- filter(rates_df2, model == "ARD")

#need to compare rates from the same trees 
#each model has 12 rates, filter for one model (ARD), label each tree 
rates_df1$tree_n <- rep(1:1000, each = 12)
#subtract the cetacean rate from the ruminant rate, is it faster (negative number) or slower (positive number)
rates_df2$tree_n <- rep(1:1000, each = 12)

rates_df <- cbind(rates_df1, rates_df2)
colnames(rates_df) <- c("whippo_rates", "model", "solution", "colours", "tree_n", "rumi_rates", "model", "solution", "colours", "tree_n")
rates_df <- rates_df[, c("whippo_rates", "solution", "tree_n", "rumi_rates")]
rates_df$difference <- rates_df$whippo_rates - rates_df$rumi_rates
#difference is negative or small -whippo is much faster, difference is positive or large ruminants have similar rates

ggplot(rates_df, aes(x = solution, y = difference)) + geom_jitter()
ggplot(rates_df, aes(x = whippo_rates, y = rumi_rates, colour = solution)) + geom_point() + facet_wrap(~solution)
ggplot(rates_df, aes(x = whippo_rates, fill = solution)) + geom_histogram() + facet_wrap(~solution)

rates_df %>% group_by(solution) %>% summarize(mean_rates = mean(whippo_rates)) %>% ggplot(., aes(x = mean_rates, y = solution)) + geom_point()
rates_df %>% group_by(solution) %>% summarize(mean_rates = mean(rumi_rates)) %>% ggplot(., aes(x = mean_rates, y = solution)) + geom_point()
