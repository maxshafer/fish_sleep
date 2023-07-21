library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)
library(geiger)
library(phytools)
library(rfishbase)
library(xlsx)
library(RColorBrewer)
library(here)

setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

index_list <- list()
index_list[[1]] <- c("all", "only_highqual", "only_cartilaginous", "only_ingroup")
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

for (i in 1:length(index_list)) {
  dataset_variable <- names(index_list)[[i]]
  
  for (j in 1:length(index_list[[i]])) {
    name_variable <- index_list[[i]][[j]]
    
    ## Load in the tree
    trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c, only_model_data = FALSE)
    trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable, only_model_data = FALSE)
    colnames(trait.data_n) <- tolower(colnames(trait.data_n))
    
    ## Make Figures
    ## #1 Make a plot showing species names and diel activity
    diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data_n[,c("species", "diel2", "order")]
    diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(y=y, x=x+15, fill = diel2), width = 30, inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = c("#abdda4", "#d7191c", "#2c7bb6", "#fdae61"))
    diel.plot.labelled <- diel.plot + geom_tiplab(color = "black", size = 1.5, offset = 30)
    
    pdf(file = here(paste("outs/Figures/plot_1_AllData", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 40, height = 40)
    print(diel.plot.labelled)
    dev.off()
    
    
    
    ## #2 Make a plot showing diel activity and highlighting orders
    orders <- unique(trait.data_n$order)[!(is.na(unique(trait.data_n$order)))]
    my_colors <- hue_pal()(length(orders))
    diel.plot.orders <- addOrderLabels(diel.plot.orders = diel.plot, colours = my_colors, orders = orders, resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 45, alpha = 0.25)
    
    pdf(file = here(paste("outs/Figures/plot_2_AllData_Orders", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 40, height = 40)
    print(diel.plot.orders + xlim(0,600))
    dev.off()
    
    ## #3 Make diurnal/nocturnal only plots
    ## For loading the models, need to have subsets
    trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c, only_model_data = TRUE)
    trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable, only_model_data = TRUE)
    
    models <- readRDS(file = paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
    model <- models$HMM_2state_2rate_marg
    
    lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
    colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2", "diurnal_R3", "nocturnal_R3", "diurnal_R4", "nocturnal_R4")[1:ncol(lik.anc)]
    
    lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
    lik.anc$noc.sum <- rowSums(lik.anc[,grep("nocturnal", colnames(lik.anc))])
    lik.anc$di.sum <- rowSums(lik.anc[,grep("diurnal", colnames(lik.anc))])
    lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
    lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])
    
    if (any(grepl("R3", colnames(lik.anc)))) {
      lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
      lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])
      lik.anc$R3.sum <- rowSums(lik.anc[,c(5,6)])
    }
    
    rate_plots <- lapply(grep("urnal_R",colnames(lik.anc)), function(x) {
      palette <- c("Reds", "Blues", "OrRd", "Purples", "YlOrBr", "BuPu")[x]
      variable <- colnames(lik.anc)[x]
      plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = .data[[variable]]) + scale_color_distiller(palette = palette, direction = 1) + geom_tippoint(aes(color = .data[[variable]]), shape = 16, size = 1.5) + scale_color_distiller(palette = palette, direction = 1)
      return(plot)
    })
    
    rate_plots <- wrap_plots(rate_plots, ncol = 2)
    ancestral_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
    rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")
    
    # Add markers for those nodes where the reconstruction is confident (80% or higher)
    lik.anc$cutoff <- ifelse(lik.anc$noc.sum > 0.8 | lik.anc$di.sum > 0.8, TRUE, FALSE)
    ancestral_plot$data$cutoff <- ifelse(!(ancestral_plot$data$isTip) & lik.anc$cutoff, TRUE, FALSE)
    
    #### Save the plots
    # pdf(file = here(paste("outs/Figures/plot_3_AncRecon_CutoffNodesLabel", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 40, height = 40)
    # print(ancestral_plot + geom_text2(aes(subset=cutoff, label = node), colour = "black"))
    # dev.off()
    
    pdf(file = here(paste("outs/Figures/plot_3_AncRecon_Rate", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 40, height = 40)
    print(rate_plot)
    dev.off()
    
    pdf(file = here(paste("outs/Figures/plot_4_AncRecon_StateRate", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 40, height = 40)
    print(rate_plots)
    dev.off()
    
    pdf(file = here(paste("outs/Figures/plot_5_AncRecon_CutoffNodes", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 20, height = 20)
    print(ancestral_plot + geom_point(data = subset(ancestral_plot$data, subset = cutoff), aes(x = x, y = y), size = 1.5, colour = "black"))
    dev.off()
    
    pdf(file = here(paste("outs/Figures/plot_6_AncRecon", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 20, height = 20)
    print(ancestral_plot)
    dev.off()
    
    pdf(file = here(paste("outs/Figures/plot_7_AncRecon_labelled", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 40, height = 40)
    print(ancestral_plot + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5))
    dev.off()
    
    pdf(file = here(paste("outs/Figures/plot_8_AncRecon_small", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 10, height = 10)
    print(ancestral_plot)
    dev.off()
    
    ## #4 Add order highlights to ancestral plot
    
    orders <- c("Anguilliformes", "Cypriniformes", "Characiformes", "Siluriformes", "Kurtiformes", "Gobiiformes", "Syngnathiformes", "Pleuronectiformes", "Cichliformes", "Perciformes/Scorpaenoidei", "Tetraodontiformes", "Carcharhiniformes", "Holocentriformes", "Blenniiformes")
    my_colors <- hue_pal()(length(orders))
    ancestral.plot <- addOrderLabels(diel.plot.orders = ancestral_plot, colours = my_colors, orders = orders, resolved_names = trait.data_n, tr.calibrated = trpy_n, highlight = FALSE, offset = 15, alpha = 0.25)
    
    pdf(file = here(paste("outs/Figures/plot_9_AncRecon_orders", dataset_variable, name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_")), width = 20, height = 20)
    print(ancestral.plot + xlim(0,500))
    dev.off()
    
  }
}
    
    
   





