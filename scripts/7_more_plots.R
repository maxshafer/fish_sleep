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
library(viridis)

#### This script to make some extra figures, based on extinction, fossil, and other data

setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

index_list <- list()
index_list[[1]] <- c("only_ingroup")
index_list[[2]] <- c("only_cartilaginous")
index_list[[3]] <- c("all")
index_list[[4]] <- c("sauropsids")
index_list[[5]] <- c("amphibians")
names(index_list) <- c("fish", "fish", "mammals", "tetrapods", "tetrapods")

# Set simulation parameters (type and #)
model_types <- "HR"
sim_numb <- 500

recon <- "marg"

smooth <- 5

anc_states <- list()
anc_rates <- list()
switch.ratio <- list()
switch.ratio.rates <- list()

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
      anc_rates[[i]] <- readRDS(file = paste("diel_ancestral_rates", dataset_variable, name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      
      switch.ratio[[i]] <- switchRatio(ancestral_states = anc_states[[i]], phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = smooth)
      switch.ratio.rates[[i]] <- switchRatio(ancestral_states = anc_rates[[i]], phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = smooth)
      
}


# ## Load data from Rohde et al 2005
# url <- 'https://docs.google.com/spreadsheets/d/1mdd1tQ46U6MFIoFkqilHgpQXf24ue93RFoD6M6iZODo/edit#gid=0'
# rohde_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
# rohde_data <- rohde_data[,c("Time_Ma", "Extinction_Intensity", "Origination_Intensity")]
# colnames(rohde_data) <- c("node.age", "extinction", "origination")

## Smooth the data by 5 my (same as mine)
# rohde_data$node.age <- cut_width(rohde_data$node.age, smooth, label = FALSE)
# rohde_data <- rohde_data %>% group_by(node.age) %>% summarise(extinction = mean(extinction), origination = mean(origination))
# rohde_data$node.age <- (rohde_data$node.age*smooth)#-smooth

# ## Make extinction plots
# extinction_rohde <- ggplot(rohde_data, aes(x = node.age, y = extinction)) + geom_line() + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + theme_classic()
# origination <- ggplot(rohde_data, aes(x = node.age, y = origination)) + geom_line(colour = "red") + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + ylim(0,30) + theme_classic()


## Load data from Song et al 2021
url <- 'https://docs.google.com/spreadsheets/d/14B_TtvUbl0NaXV2kw50qDkKBz38iG7Vx6MSxHUC78O8/edit#gid=0'
song_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
song_data <- song_data[,c("Age..base.", "Time.span", "ΔT...C.", "GF", "X3T")]
colnames(song_data) <- c("node.age", "time.span", "temp", "extinction_gf", "extinction_3t")

## Make a combined plot
names(switch.ratio) <- unlist(lapply(seq_along(index_list), function(x) paste(names(index_list)[x], index_list[x], sep = "_")))
combined.data <- Reduce(rbind, lapply(seq_along(switch.ratio), function(x) {
  df <- switch.ratio[[x]]$data
  df$group <- names(switch.ratio)[[x]]
  return(df)
}))

names(switch.ratio.rates) <- unlist(lapply(seq_along(index_list), function(x) paste(names(index_list)[x], index_list[x], sep = "_")))
combined.data.rates <- Reduce(rbind, lapply(seq_along(switch.ratio.rates), function(x) {
  df <- switch.ratio.rates[[x]]$data
  df$group <- names(switch.ratio.rates)[[x]]
  return(df)
}))

## Make plots for diel transitions and extinction/temp

theme_set(theme_classic(base_size = 8))

## Make test figure 2a
combined.plot <- ggplot(combined.data[combined.data$group %in% c("fish_only_ingroup", "mammals_all", "tetrapods_amphibians", "tetrapods_sauropsids"),], aes(x = node.age, y = ratio, group = group, colour = group)) + scale_color_viridis(discrete = TRUE, option = "D") + theme_classic() + scale_x_reverse()
geo_scale <- gggeo_scale(combined.plot, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", axis.text.x = element_text(colour = "black"), axis.title.x = element_text(colour = "black"))  + coord_cartesian(xlim = abs(layer_scales(combined.plot)$x$range$range)) 
geo_scale <- geo_scale + xlab("Millions of years ago (mya)")

combined.plot <- combined.plot + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + geom_line(size = 2, alpha = 0.75)
combined.plot <- combined.plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(colour = "black"))
combined.plot <- combined.plot + ylab("Cummulative\ntransitions / lineages")  + theme(legend.position = "none")

combined.plot.rates <- ggplot(combined.data.rates[combined.data.rates$group %in% c("fish_only_ingroup", "mammals_all", "tetrapods_amphibians", "tetrapods_sauropsids"),], aes(x = node.age, y = ratio, group = group, colour = group)) + scale_color_viridis(discrete = TRUE, option = "D") + theme_classic() + scale_x_reverse()
combined.plot.rates <- combined.plot.rates + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + geom_line(size = 2, alpha = 0.75)
combined.plot.rates <- combined.plot.rates + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(colour = "black"))
combined.plot.rates <- combined.plot.rates + ylab("Cummulative\ntransitions / lineages")  + theme(legend.position = "none")

extinction_song <- ggplot(song_data[complete.cases(song_data$extinction_gf),], aes(x = node.age-time.span, y = extinction_gf)) + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + theme_classic() + coord_cartesian(xlim = abs(layer_scales(combined.plot)$x$range$range))
extinction_song <- extinction_song + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + geom_line(size = 1) 
extinction_song <- extinction_song + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(colour = "black")) + ylab("Extinction\nrate")

temp_song <- ggplot(song_data[complete.cases(song_data$extinction_gf),], aes(x = node.age-time.span, y = abs(temp))) + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + theme_classic() + coord_cartesian(xlim = abs(layer_scales(combined.plot)$x$range$range))
temp_song <- temp_song + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + geom_line(size = 1, colour = "red")
temp_song <- temp_song + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(colour = "black")) + ylab("Temperature\nchange")

pdf(file = here("outs/Figures/plot_XX_transitions_vs_song_data.pdf"), width = 4, height = 6)
combined.plot / extinction_song / temp_song / geo_scale + plot_layout(nrow = 4, heights = c(10,5,5,1.5), guides = "collect")
dev.off()

pdf(file = here("outs/Figures/plot_XX_transitions_rates.pdf"), width = 4, height = 4)
combined.plot.rates / geo_scale + plot_layout(nrow = 2, heights = c(10,1.5), guides = "collect")
dev.off()

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
db_data$group <- factor(db_data$group, levels = c("actinopterygii", "sauropsida", "amphibia", "mammalia"))

fossil_plot <- ggplot(db_data[db_data$group %in% c("actinopterygii", "amphibia", "mammalia", "sauropsida"),], aes(x = node.age, y = fossil, group = group, colour = group)) + geom_line() + scale_x_reverse() + xlim(c(max(combined.plot$data$node.age),0)) + theme_classic()
fossil_plot <- fossil_plot + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
fossil_plot <- fossil_plot + facet_wrap(~group, scales = "free_y", nrow = 4, strip.position = "bottom") + ylab("Sampled fossil diversity") + theme(legend.position = "none")

geo_scale <- gggeo_scale(fossil_plot, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", axis.text.x = element_text(colour = "black"), axis.title.x = element_text(colour = "black"))  + coord_cartesian(xlim = abs(layer_scales(fossil_plot)$x$range$range)) 
geo_scale <- geo_scale + xlab("Millions of years ago (mya)")

pdf(file = here("outs/Figures/plot_XX_fossils.pdf"), width = 4, height = 6)
fossil_plot / geo_scale + plot_layout(nrow = 2, heights = c(10,1))
dev.off()

