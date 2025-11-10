
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

filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_artiodactyla_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

rates_df <- plot1kTransitionRates4state(readRDS(here(filename)), 5)

#plot as a violin plot
#important note: the colours column doesn't line up with the correct solution in the df but if we plot the solutions in alphabetical order and then colouring them with the palette in the colours column is in the correct order
ggplot(rates_df, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df$colours) + geom_violin(color = "black", scale = "width") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  + scale_fill_manual(values = rates_df$colours) + theme(legend.position = "none") + labs(x = "Transition", y = "Log(rates)") + stat_summary(fun=median, geom="point", size=2, colour = "red") + ggtitle(filename) + facet_wrap(~fct_inorder(model), ncol = 3, nrow = 2) 

model_selection <- "Bridge_only"
#pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", filename, "_violin_rate_plot_", model_selection, ".pdf"), width = 10, height = 7)
rates_df1 <- rates_df %>% filter(model == model_selection) 
ggplot(rates_df1, aes(x= solution, y = log(rates), group = solution, fill = solution, colour = solution)) + 
  geom_jitter(aes(alpha = 0.1)) + scale_color_manual(values = rates_df1$colours) +
  geom_violin(color = "black", scale = "width") +theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size =10), axis.text.y = element_text(size =10))  +
  scale_fill_manual(values = rates_df1$colours) + theme(legend.position = "none") +
  labs(x = "Transition", y = "Log(transition rate)") + stat_summary(fun=median, geom="point", size=2, colour = "red") +
  stat_summary(geom = "point", fun.y = "mean") + stat_summary(aes(label=..y..), fun.y=mean, geom="text", size=2, colour = "black") +
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

# Section 8: Comparison of rate magnitude -------------------------------

#get rates
rates_df1 <- plot1kTransitionRates4state(readRDS(here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df1 <- filter(rates_df1, model == "Bridge_only")
rates_df2 <- plot1kTransitionRates4state(readRDS(here("august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df2 <- filter(rates_df2, model == "Bridge_only")

#need to compare rates from the same trees 
#each model has 12 rates, filter for one model (ARD), label each tree 
rates_df1$tree_n <- rep(1:1000, each = 10)
#subtract the cetacean rate from the ruminant rate, is it faster (negative number) or slower (positive number)
rates_df2$tree_n <- rep(1:1000, each = 10)

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

rates_df %>% group_by(solution) %>% 
  summarize(mean_whippo_rates = mean(whippo_rates), mean_rumi_rates = mean(rumi_rates)) %>% 
  ggplot(., aes(y = solution)) + geom_point(aes(x=mean_whippo_rates), colour = "blue") +
  geom_point(aes(x = mean_rumi_rates), colour = "red")

#make a forest-ish plot
ggplot(rates_df, aes(y = solution)) + geom_point(aes(x=whippo_rates), colour = "blue") +
  geom_point(aes(x = rumi_rates), colour = "red")


rates_df1 <- plot1kTransitionRates4state(readRDS(here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df1 <- filter(rates_df1, model == "Bridge_only")
rates_df2 <- plot1kTransitionRates4state(readRDS(here("august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df2 <- filter(rates_df2, model == "Bridge_only")

rates_df1$clade <- "whippomorpha"
rates_df2$clade <- "ruminants"

rates_df <- rbind(rates_df1, rates_df2)

df <- rates_df %>% group_by(clade, solution) %>% summarize(mean_rates = mean(rates), 
                                                           lci = t.test(rates, conf.level = 0.95)$conf.int[1],
                                                           uci = t.test(rates, conf.level = 0.95)$conf.int[2])

ggplot(df, aes(x = mean_rates, y = solution, colour = clade)) + geom_point() + geom_errorbar(aes(y = solution, xmin = lci, xmax =uci),width = 0.4)

# Section 9: Comparing rates between different clades ---------------------

rates_df1 <- plot1kTransitionRates4state(readRDS(here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)

#find the mean scores
means_whippo <- aggregate(rates ~  model + solution, rates_df, mean)
means_whippo$clade <- "whippomorpha"
#means$AIC_score <- round(means$AIC_score, digits = 2)

#repeat for ruminants
rates_df2 <- plot1kTransitionRates4state(readRDS(here("august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
means_rumi <- aggregate(rates ~  model + solution, rates_df, mean)
means_rumi$clade <- "ruminantia"

means <- rbind(means_whippo, means_rumi)

ggplot(means, aes(y = rates, x = solution, colour = clade)) + geom_point() + facet_wrap(~model)


# #Section 10: Total garbage test ------------------------------------------

filename <- "whippomorpha_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
model_results <- readRDS(here(paste0(filename, ".rds")))

trait.data <- model_results$ER_model$data
table(trait.data$max_crep)

#since there are four trait states
#likelihood <- ((factorial(n))/(factorial(n1) * factorial(n2) * factorial(n3) *factorial(n4))) * (p1^n1)  * (p2^n2)  * (p3^n3)  * (p4^n4) 

#using the natural log
n1 = 27
n2 = 22
n3 = 7
n4 = 21
n = 77

lnL_garb = n1 * log(n1 / n) + n2 * log(n2 / n) + n3 * log(n3 / n) + n4 * log(n4 / n)
#ln likelihood is -99.926

#compared to the actual likelihood
model_results$bridge_only_model$loglik #-91.0821
model_results$ER_model$loglik #-105.561
model_results$SYM_model$loglik #-97.924
model_results$ARD_model$loglik #-91.025
model_results$CONSYM_model$loglik # -97.46

likelihood_metrics <- max_clade_metrics(readRDS(here(paste0(filename, ".rds"))))
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)
likelihood_metrics <- rbind(likelihood_metrics, data.frame(model = "Total garbage", log_likelihoods = lnL_garb, AICc_scores = NA,  AIC_scores = NA))
ggplot(likelihood_metrics, aes(y = log_likelihoods, x = model, fill = model))  + geom_bar(stat = "identity")

#so the log likelihood is similar, but the actual model is more likely (higher log lik)

#for ruminants
filename <- "ruminants_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
model_results <- readRDS(here(paste0(filename, ".rds")))

trait.data <- model_results$ER_model$data
table(trait.data$max_crep)

n1 = 21
n2 = 125
n3 = 34
n4 = 23
n = 203

lnL_garb = n1 * log(n1 / n) + n2 * log(n2 / n) + n3 * log(n3 / n) + n4 * log(n4 / n)

#the garbage ln likelihood is -219 which is different than our most likely model (-199.78)

model_results$bridge_only_model$loglik #-199.7884
model_results$ER_model$loglik # -230.505
model_results$SYM_model$loglik #-219.0621
model_results$ARD_model$loglik #-199.8458
model_results$CONSYM_model$loglik # -219.0807

likelihood_metrics <- max_clade_metrics(readRDS(here(paste0(filename, ".rds"))))
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)
likelihood_metrics <- rbind(likelihood_metrics, data.frame(model = "Total garbage", log_likelihoods = lnL_garb, AICc_scores = NA,  AIC_scores = NA))
ggplot(likelihood_metrics, aes(y = log_likelihoods, x = model, fill = model))  + geom_bar(stat = "identity") + ggtitle(filename)

likelihood_metrics$length <- 100
ggplot(likelihood_metrics, aes(y = log_likelihoods, x = length, fill = model)) +
  geom_line() + ggtitle(filename)
