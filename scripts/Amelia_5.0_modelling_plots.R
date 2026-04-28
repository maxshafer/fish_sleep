
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
# model_results1 <- readRDS(here("august_12__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
# model_results2 <- readRDS(here("august_14__whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))
# model_results1$ER_model <- c(model_results1$ER_model, model_results2$ER_model)
# model_results1$SYM_model <- c(model_results1$SYM_model, model_results2$SYM_model)
# model_results1$ARD_model <- c(model_results1$ARD_model, model_results2$ARD_model)
# model_results1$bridge_only_model <- c(model_results1$bridge_only_model, model_results2$bridge_only_model)
# model_results1$CONSYM_model <- c(model_results1$CONSYM_model, model_results2$CONSYM_model)
# 
# saveRDS(model_results1, here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"))

# Section 4: Plot AIC scores from 1k model results ----------------------

filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_artiodactyla_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

#requires the filename and the number of Mk models (3: ER, SYM, ARD or 4: ER, SYM, ARD, CONARD, 5: ER, SYM, ARD, bridge_only, CONSYM)
#returns a df of the AIC scores for all 1k trees x number of Mk models
df_full <- plot1kAIC(readRDS(here(filename)), 5)
df_full$model <- factor(df_full$model, levels = c("ER", "SYM", "CONSYM", "ARD", "bridge_only"))

means <- aggregate(AIC_score ~  model, df_full, mean)
means$AIC_score <- round(means$AIC_score, digits = 2)

#calculate the delta AIC (difference between the AIC scores)
means <- means %>% mutate(delta_AIC = round(AIC_score - min(AIC_score), digits = 2))
#allows for the trailing zeros
means$delta_AIC <- formatC(means$delta_AIC, format = "f", digits = 2)

custom.colours <- c("#013873", "#04549F", "#056CCC", "#60B0FB", "#8DC6FC")

#plot and save out - raincloud plot
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "raincloud", filename, "_AIC", ".pdf"), width = 7, height = 5)
ggplot(df_full, aes(x = model, y = AIC_score, fill = model)) + 
  scale_fill_manual(values = custom.colours) +
  ggdist::stat_halfeye(alpha = 0.6, adjust = .5, width = .6, justification = -.3, .width = 0, point_colour = NA) +
  geom_boxplot(alpha = 0.2, width = .25, colour = "black",outlier.shape = NA) + 
  geom_point(aes(color = model), stroke = 1, size = 1, alpha = .2, position = position_jitter(seed = 1, width = .15)) +
  scale_color_manual(values = custom.colours) + 
  geom_text(data = means, aes(label = delta_AIC, y = AIC_score), hjust = -0.55, size = 5) +
  labs(x = "Model", y = "AIC score")  + theme_bw() +
  scale_x_discrete(labels = c("ER", "SYM", "CON-SYM", "ARD", "CON-ARD")) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 16), legend.position = "none") + coord_cartesian(ylim = c(398, 525)) #there is one extreme outlier in the ER models of 564
dev.off()

# Section 5: Plot transition rates from 1k model results ----------------------

filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_artiodactyla_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

#requires the filename, the number of states in the model and the number of Mk models 
#returns a dataframe of the rates from each of the Mk models, for each of the 1k trees
rates_df <- plot1kTransitionRates4state(readRDS(here(filename)), 5)

model_selection <- "SYM"
#colour by unique transition (solution)
rate_colours_sol <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#AD9680", "#D1B49B","#EECBAD",  "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
#colour by ending state
rate_colours_end = c("#AD9680","#FA4A05","#3C967E","#A024AE", "#FA4A05","#3C967E", "#A024AE","#AD9680", "#3C967E","#A024AE","#AD9680", "#FA4A05")

#if model selection is bridge only, adjust colours
model_selection <- "Bridge_only"
#colour by ending state
rate_colours_end = c("#AD9680","#FA4A05","#3C967E","#A024AE", "#FA4A05","#3C967E", "#A024AE","#AD9680", "#A024AE","#AD9680")

#filter by the model you're plotting
rates_df1 <- rates_df %>% filter(model == model_selection) 

#extract just the starting state
rates_df1$start_state <- word(rates_df1$solution, 1)
rates_df1$start_state <- paste(rates_df1$start_state, "to", sep = " ")
rates_df1$end_state <- word(rates_df1$solution, 3)

rates_df1$start_state <- factor(rates_df1$start_state, levels = c("Cathemeral to", "Diurnal to", "Crepuscular to", "Nocturnal to"))

rates_plot <- 
  ggplot(rates_df1, aes(x= end_state, y = log(rates), group = solution, fill = solution, colour = solution)) + 
  geom_quasirandom(alpha = 0.3, width = 0.5, method = "quasirandom") + 
  scale_color_manual(values = rate_colours_end) +
  geom_violin(color = "black", scale = "width", alpha = 0.5) + theme_bw() +
  scale_fill_manual(values = rate_colours_end) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10), axis.text.y = element_text(size =10), axis.title = element_text(size = 12), strip.background = element_rect(fill = "grey90"), legend.position = "none")  +
  labs(x = "\n Transition", y = "Log(transition rate)") + 
  stat_summary(fun=median, geom="point", size=3, colour = "black", alpha = 0.2) 

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", substr(filename, start = 1, stop = 15), "_square_violin_rate_plot_", model_selection, ".pdf"), width = 7, height = 7)
rates_plot +
  facet_wrap(~start_state, scales = "free_x", nrow = 2, ncol = 2)
dev.off()

density_plot <-
  rates_df1 %>%
  ggplot(., aes(y = log(rates), group = solution, fill = end_state)) +
  geom_density(alpha = 0.4) + theme_bw() +
  scale_fill_manual(values = c( "#A024AE","#AD9680",  "#FA4A05","#3C967E")) 

#kernel density estimates of rates
# pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", filename, "rate_kernel_density_plot", model_selection, ".pdf"), width = 7, height = 7)
density_plot  +
  facet_wrap(~start_state, nrow = 2, ncol = 2)
# dev.off()

#combined ggplot and violin plot in column format
rates_column_plot <-
  rates_plot +
  facet_wrap(~start_state, scales = "free_x", nrow = 4, ncol = 1)

density_column_plot <- 
  density_plot + 
  facet_wrap(~start_state, nrow = 4, ncol = 1) +
  labs(x = "\n Density", y = "") +
  #for the ARD ruminant plot only either set scales to be free or limit x axis to 0.3, removes the large number of near zeros skewing the plot
  #xlim(0, 0.4) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10), axis.text.y = element_text(size =10), axis.title = element_text(size = 12), strip.background = element_rect(fill = "grey90"), legend.position = "none") 
  
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", substr(filename, start = 1, stop = 15), "_combined_rate_plot_column_trimmed_", model_selection, ".pdf"), width = 7, height = 10)
ggarrange(rates_column_plot, density_column_plot)
dev.off()

#combined ggplot and violin plot in row format
# rates_column_plot <-
#   rates_plot +
#   facet_wrap(~start_state, scales = "free_x", nrow = 1, ncol = 4)
# 
# density_column_plot <- 
#   density_plot + 
#   facet_wrap(~start_state, nrow = 1, ncol = 4) +
#   labs(x = "\n Density", y = "") +
#   theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10), axis.text.y = element_text(size =10), axis.title = element_text(size = 12), strip.background = element_rect(fill = "grey90"), legend.position = "none") 
# 
# pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", filename, "_combined_rate_plot_row_", model_selection, ".pdf"), width = 10, height = 7)
# ggarrange(rates_column_plot, density_column_plot,
#           ncol = 1, nrow = 2)
# dev.off()


# Plot rates from 100 most likely models ----------------------------------

filename <- "august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"
#filename <- "august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds"

rates_df <- plot1kTransitionRates4state(readRDS(here(filename)), 5)

model_selection <- "SYM"

#filter by the model you're plotting
rates_df1 <- rates_df %>% filter(model == model_selection) 

df_full <- plot1kAIC(readRDS(here(filename)), 5)
df_full$model <- factor(df_full$model, levels = c("ER", "SYM", "CONSYM", "ARD", "bridge_only"))

df_full <- df_full %>% filter(model == model_selection)

df_full$model_number <- 1:nrow(df_full)

rates_df1$model_number <- rep(1:1000, each = (nrow(rates_df1)/1000))

#merge by model number

rates_df1 <- merge(rates_df1, df_full[, c("AIC_score", "model_number")], by = "model_number", all = TRUE)

#create list of 100 best models (lowest AIC score)
lowest_100 <- rates_df1 %>% arrange(AIC_score) %>% select(model_number) %>% slice(1:(nrow(rates_df1)/10))

#filter by list
rates_df1 <- rates_df1 %>% filter(model_number %in% unique(lowest_100$model_number))

#plot as usual

#extract just the starting state
rates_df1$start_state <- word(rates_df1$solution, 1)
rates_df1$start_state <- paste(rates_df1$start_state, "to", sep = " ")
rates_df1$end_state <- word(rates_df1$solution, 3)

rates_df1$start_state <- factor(rates_df1$start_state, levels = c("Cathemeral to", "Diurnal to", "Crepuscular to", "Nocturnal to"))

#ARD colours
rate_colours_end = c("#AD9680","#FA4A05","#3C967E","#A024AE", "#FA4A05","#3C967E", "#A024AE","#AD9680", "#3C967E","#A024AE","#AD9680", "#FA4A05")
#bridge colours
#rate_colours_end = c("#AD9680","#FA4A05","#3C967E","#A024AE", "#FA4A05","#3C967E", "#A024AE","#AD9680", "#A024AE","#AD9680")

rates_plot <- 
  ggplot(rates_df1, aes(x= end_state, y = log(rates), group = solution, fill = solution, colour = solution)) + 
  geom_quasirandom(alpha = 0.8, width = 0.5, method = "quasirandom") + 
  scale_color_manual(values = rate_colours_end) +
  geom_violin(color = "black", scale = "width", alpha = 0.5) + theme_bw() +
  scale_fill_manual(values = rate_colours_end) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10), axis.text.y = element_text(size =10), axis.title = element_text(size = 12), strip.background = element_rect(fill = "grey90"), legend.position = "none")  +
  labs(x = "\n Transition", y = "Log(transition rate)") + 
  stat_summary(fun=median, geom="point", size=3, colour = "black", alpha = 0.2) +
  facet_wrap(~start_state, scales = "free_x", nrow = 2, ncol = 2)

rates_plot

