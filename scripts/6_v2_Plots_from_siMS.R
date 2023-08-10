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

#### This script loads in the results from simulations and SIMMAPs to make figures and combined figures
#### For each dataset/comparison there should be a 'simulations' and a "SIMMAPs' .rds file and 'cumsums' rds files


setwd(here())

source(here("scripts/Fish_sleep_functions.R"))

index_list <- list()
index_list[[1]] <- c("all", "only_cartilaginous", "only_ingroup") # , "only_highqual"
index_list[[2]] <- c("all")
index_list[[3]] <- c("all", "not_mammals", "amphibians", "sauropsids", "lepidosauria", "testudines", "aves")
index_list[[4]] <- c("all")
names(index_list) <- c("fish", "mammals", "tetrapods", "AllGroups")

# Set simulation parameters (type and #)
model_types <- c("ARD", "HR")
sim_numb <- 500

recon <- "marg"

theme_set(theme_classic(base_size = 8))

sIMs_states_plots <- list()
sIMs_rates_plots <- list()

for (y in 1:length(model_types)) {
  model_type <- model_types[[y]]
  
  for (i in 1:length(index_list)) {
    dataset_variable <- names(index_list)[[i]]
    
    sIMs_states_plots[[i]] <- list()
    sIMs_rates_plots[[i]] <- list()
    
    for (j in 1:length(index_list[[i]])) {
      name_variable <- index_list[[i]][[j]]
      
      setwd(here())
      
      ###############################################################################################################################################
      ### Load files ### 
      ###############################################################################################################################################
      
      ## Load in the tree
      trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
      
      ## Load the best model, which is the HMM 2 state 2 rate model
      models <- readRDS(file = here(paste("marginal_and_joint_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_")))
      
      ## I think I will always take the 2 rate model, 3 is too hard to comprehend
      
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
      anc_states <- readRDS(file = paste("diel_ancestral_states", dataset_variable, name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      anc_rates <- readRDS(file = paste("diel_ancestral_rates", dataset_variable, name_variable, Ntip(trpy_n), "species.rds", sep = "_"))
      
      
      ## Load cumsum files for states and rates, and for sims and SIMs
      simulated_cumsums <- readRDS(file = here(paste("diel_switch_simulations_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
      
      SIMMAP_cumsums <- readRDS(file = here(paste("diel_switch_SIMMAP_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
      
      if (model_type == "HR") {
        simulated_cumsums_rates <- readRDS(file = here(paste("diel_switch_simulations_rates_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
        SIMMAP_cumsums_rates <- readRDS(file = here(paste("diel_switch_SIMMAP_rates_cumsums", dataset_variable, name_variable, Ntip(trpy_n), "species", model_type, sim_numb, "rTraitDisc.rds", sep = "_")))
      }
      
      ###############################################################################################################################################
      ### Plot the results of the simulation, along with the actual data ### 
      ###############################################################################################################################################
      
      ## Plot the results of the simulation, either the summary (mean +/- stdev), the simulations, or both ('overlay')
      plot_simulation <- simulatedSwitchRatio(simulated_cumsums = simulated_cumsums, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "red", error = "sd", smooth = 5) #+ ylim(c(0,0.15))
      
      ## Can I add the computed value? Actual transitions
      ## Easy way is to make the plot, then use the data from it to add to the above
      switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = 5)
      
      plot <- plot_simulation + geom_line(data = switch.ratio$data, aes(x=node.age,y=ratio), colour = "black")
      plot <- plot + ggtitle(paste(model_type, "model simulation,", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", sim_numb, ", Avg. transitions:", round(mean(colMax(simulated_cumsums)),1), "+/-", round(sd(colMax(simulated_cumsums)),1), sep = " "))
      
      # This 'zooms-in', which preserves the ribbon (setting limits removes part of the ribbon stdev)
      recon_plot <- plot + coord_cartesian(ylim = c(0, layer_scales(plot)$y$range$range[[2]]), xlim = abs(layer_scales(plot)$x$range$range)) 
      
      # Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
      geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot)$x$range$range)) 
      
      simmulation_recon_plot <- recon_plot / geo_scale + plot_layout(nrow = 2, heights = c(5,0.5))
      
      if (model_type == "HR") {
        plot_simulation_rates <- simulatedSwitchRatio(simulated_cumsums = simulated_cumsums_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "red", error = "sd", smooth = 5) #+ ylim(c(0,0.15))
        switch.ratio.rates <- switchRatio(ancestral_states = anc_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = 5)
        plot.rates <- plot_simulation_rates + geom_line(data = switch.ratio.rates$data, aes(x=node.age,y=ratio), colour = "black")
        plot.rates <- plot.rates + ggtitle(paste(model_type, "model simulation (rates),", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", sim_numb, ", Avg. transitions:", round(mean(colMax(simulated_cumsums_rates)),1), "+/-", round(sd(colMax(simulated_cumsums_rates)),1), sep = " "))
        recon_plot_rates <- plot.rates + coord_cartesian(ylim = c(0, layer_scales(plot.rates)$y$range$range[[2]]), xlim = abs(layer_scales(plot.rates)$x$range$range)) 
        geo_scale_rates <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot.rates)$x$range$range))
        
        simmulation_recon_plot_rates <- recon_plot_rates / geo_scale_rates + plot_layout(nrow = 2, heights = c(5,0.5))
      }
      
      
      ###############################################################################################################################################
      ### Plot the results of the simulation, along with the actual data ### 
      ###############################################################################################################################################
      
      ## Plot the results of the simulation, either the summary (mean +/- stdev), the simulations, or both ('overlay')
      plot_simmap <- simulatedSwitchRatio(simulated_cumsums = SIMMAP_cumsums, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "blue", error = "sd", smooth = 5) #+ ylim(c(0,0.15))
      
      ## Can I add the computed value? Actual transitions
      ## Easy way is to make the plot, then use the data from it to add to the above
      switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = 5)
      
      plot <- plot_simmap + geom_line(data = switch.ratio$data, aes(x=node.age,y=ratio), colour = "black")
      plot <- plot + ggtitle(paste(model_type, "model SIMMAP,", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", sim_numb, ", Avg. transitions:", round(mean(colMax(SIMMAP_cumsums)),1), "+/-", round(sd(colMax(SIMMAP_cumsums)),1), sep = " "))
      
      # This 'zooms-in', which preserves the ribbon (setting limits removes part of the ribbon stdev)
      recon_plot <- plot + coord_cartesian(ylim = c(0, layer_scales(plot)$y$range$range[[2]]), xlim = abs(layer_scales(plot)$x$range$range)) 
      
      # Make an ggplot that is just the scale (plot on top of switch.ratio and cut the scales with blank.gg = TRUE)
      geo_scale <- gggeo_scale(switch.ratio, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = c(415,0)) 
      
      simmap_recon_plot <- recon_plot / geo_scale + plot_layout(nrow = 2, heights = c(5,0.5))
      
      if (model_type == "HR") {
        plot_simmap_rates <- simulatedSwitchRatio(simulated_cumsums = SIMMAP_cumsums_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, plot_type = "summary", highlight_colour = "blue", error = "sd", smooth = 5) #+ ylim(c(0,0.15))
        switch.ratio.rates <- switchRatio(ancestral_states = anc_rates, phylo_tree = trpy_n, node.age.cutoff = 0.02, smooth = 5)
        plot.rates <- plot_simmap_rates + geom_line(data = switch.ratio.rates$data, aes(x=node.age,y=ratio), colour = "black")
        plot.rates <- plot.rates + ggtitle(paste(model_type, "model SIMMAP (rates),", name_variable, Ntip(trpy_n), "species", sep = " "), subtitle = paste("# Simulations:", sim_numb, ", Avg. transitions:", round(mean(colMax(SIMMAP_cumsums_rates)),1), "+/-", round(sd(colMax(SIMMAP_cumsums_rates)),1), sep = " "))
        recon_plot_rates <- plot.rates + coord_cartesian(ylim = c(0, layer_scales(plot.rates)$y$range$range[[2]]), xlim = abs(layer_scales(plot.rates)$x$range$range)) 
        geo_scale_rates <- gggeo_scale(switch.ratio.rates, pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme_void() + coord_cartesian(xlim = abs(layer_scales(plot.rates)$x$range$range))
        
        simmap_recon_plot_rates <- recon_plot_rates / geo_scale_rates + plot_layout(nrow = 2, heights = c(5,0.5))
      }
      
      ## Save out the plots
      pdf(file = here(paste("outs/Figures/plot_16_simulations_states", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, sim_numb, "zoom.pdf", sep = "_")), width = 10, height = 10)
      print(simmulation_recon_plot)
      dev.off()

      pdf(file = here(paste("outs/Figures/plot_17_SIMMAPs_states", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, sim_numb, "zoom.pdf", sep = "_")), width = 10, height = 10)
      print(simmap_recon_plot)
      dev.off()
      
      if (model_type == "HR") {
        pdf(file = here(paste("outs/Figures/plot_18_simulations_rates", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, sim_numb, "zoom.pdf", sep = "_")), width = 10, height = 10)
        print(simmulation_recon_plot_rates)
        dev.off()
        pdf(file = here(paste("outs/Figures/plot_19_SIMMAPs_rates", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, sim_numb, "zoom.pdf", sep = "_")), width = 10, height = 10)
        print(simmap_recon_plot_rates)
        dev.off()
      }
      
      ###############################################################################################################################################
      ### Plot simulations vs SIMMAPs ### 
      ###############################################################################################################################################
      
      merge_data <- merge(plot_simulation$data, plot_simmap$data, by = "node.age", suffixes = c("_simulation", "_simmap"))
      merge_plot <- ggplot(merge_data, aes(x = node.age)) + geom_line(aes(x=node.age,y=mean_simulation), colour = "red") + geom_line(aes(x=node.age, y=mean_simmap), colour = "blue") + geom_ribbon(aes(ymin = ifelse(mean_simulation-sd_simulation < 0, 0, mean_simulation-sd_simulation), ymax = mean_simulation+sd_simulation), fill = "red", alpha = 0.1) + geom_ribbon(aes(ymin = ifelse(mean_simmap-sd_simmap < 0, 0, mean_simmap-sd_simmap), ymax = mean_simmap+sd_simmap), fill = "blue", alpha = 0.1)
      merge_plot <- merge_plot + theme_classic() + scale_x_reverse() + ylab("Fraction of lineages transitioning") + xlab("Millions of years ago") + xlim(c(415,0))
      merge_plot <- merge_plot + ggtitle(paste(model_type, "model STATES: Simulations vs SIMMAP", sep = " "), subtitle = paste("# Simulations:", sim_numb, "-", dataset_variable, name_variable, Ntip(trpy_n), "species", sep = " "))
      sIMs_states_plots[[i]][[j]] <- merge_plot
      
      merge_plot <- merge_plot / geo_scale + plot_layout(nrow = 2, heights = c(5,0.5))
      
      pdf(file = here(paste("outs/Figures/plot_20_simulations_vs_SIMMAPs_states", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, sim_numb, "zoom.pdf", sep = "_")), width = 10, height = 10)
      print(merge_plot)
      dev.off()
      
      
      
      if (model_type == "HR") {
        merge_data_rates <- merge(plot_simulation_rates$data, plot_simmap_rates$data, by = "node.age", suffixes = c("_simulation", "_simmap"))
        merge_plot_rates <- ggplot(merge_data_rates, aes(x = node.age)) + geom_line(aes(x=node.age,y=mean_simulation), colour = "red") + geom_line(aes(x=node.age, y=mean_simmap), colour = "blue") + geom_ribbon(aes(ymin = ifelse(mean_simulation-sd_simulation < 0, 0, mean_simulation-sd_simulation), ymax = mean_simulation+sd_simulation), fill = "red", alpha = 0.1) + geom_ribbon(aes(ymin = ifelse(mean_simmap-sd_simmap < 0, 0, mean_simmap-sd_simmap), ymax = mean_simmap+sd_simmap), fill = "blue", alpha = 0.1)
        merge_plot_rates <- merge_plot_rates + theme_classic() + scale_x_reverse() + ylab("Fraction of lineages transitioning") + xlab("Millions of years ago") + xlim(c(415,0))
        merge_plot_rates <- merge_plot_rates + ggtitle(paste(model_type, "model RATES: Simulations vs SIMMAP", sep = " "), subtitle = paste("# Simulations:", sim_numb, "-", dataset_variable, name_variable, Ntip(trpy_n), "species", sep = " "))
        sIMs_rates_plots[[i]][[j]] <- merge_plot_rates
        
        merge_plot_rates <- merge_plot_rates / geo_scale_rates + plot_layout(nrow = 2, heights = c(5,0.5))
        
        pdf(file = here(paste("outs/Figures/plot_21_simulations_vs_SIMMAPs_rates", dataset_variable, name_variable, length(trpy_n$tip.label), "species", model_type, sim_numb, "zoom.pdf", sep = "_")), width = 10, height = 10)
        print(merge_plot_rates)
        dev.off()
      }
      

    }
  }
  
}


## Make a combined plot for figures
actinos <- sIMs_states_plots[[1]][[3]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
actinos <- actinos + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
mammals <- sIMs_states_plots[[2]][[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
mammals <- mammals + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
amphibians <- sIMs_states_plots[[3]][[3]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
amphibians <- amphibians + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
sauropsids <- sIMs_states_plots[[3]][[4]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
sauropsids <- sauropsids + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)

geo_scale <- gggeo_scale(sIMs_states_plots[[1]][[3]], pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", axis.text.x = element_text(colour = "black"), axis.title.x = element_text(colour = "black"))  + coord_cartesian(xlim = abs(layer_scales(actinos)$x$range$range)) 
geo_scale <- geo_scale + xlab("Millions of years ago (mya)") + ggtitle(element_blank(), subtitle = element_blank())


pdf(file = here("outs/Figures/plot_XX_models_combined.pdf"), width = 4, height = 6)
actinos + amphibians + mammals + sauropsids + geo_scale + plot_layout(nrow = 5, heights = c(5,5,5,5,1.5))
dev.off()


## Make a combined plot for figures
actinos <- sIMs_rates_plots[[1]][[3]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
actinos <- actinos + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
mammals <- sIMs_rates_plots[[2]][[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
mammals <- mammals + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
amphibians <- sIMs_rates_plots[[3]][[3]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
amphibians <- amphibians + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)
sauropsids <- sIMs_rates_plots[[3]][[4]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black")) + ggtitle(element_blank(), subtitle = element_blank())
sauropsids <- sauropsids + annotate("rect", xmin = 145-15, xmax = 145+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + annotate("rect", xmin = 66-15, xmax = 66+15, ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1)

geo_scale <- gggeo_scale(sIMs_rates_plots[[1]][[3]], pos = "top", blank.gg = TRUE) + scale_x_reverse() + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", axis.text.x = element_text(colour = "black"), axis.title.x = element_text(colour = "black"))  + coord_cartesian(xlim = abs(layer_scales(actinos)$x$range$range)) 
geo_scale <- geo_scale + xlab("Millions of years ago (mya)") + ggtitle(element_blank(), subtitle = element_blank())


pdf(file = here("outs/Figures/plot_XX_models_combined_rates.pdf"), width = 4, height = 6)
actinos + amphibians + mammals + sauropsids + geo_scale + plot_layout(nrow = 5, heights = c(5,5,5,5,1.5))
dev.off()
