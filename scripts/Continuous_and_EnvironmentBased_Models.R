library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(rfishbase)
library(stringr)
library(ggtree)
library(ggplot2)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/scripts/Trait_rate_matrix_figure_script.R")

## Load files

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
tr.calibrated <- readRDS("calibrated_phylo.rds")
trait.data <- readRDS("trait_data.rds")

name_variable <- "all"

# Remove low quality data
trait.data <- trait.data[trait.data$confidence > 1,]
name_variable <- "only_highqual"

# Remove cartilaginous fish (just actinopterygii)
# Common ancestor of Lepisosteus_osseus (Order: Lepisosteiformes) Lepidosiren_paradoxa (lepidosiren) and Lutjanus_fulvus (Order: Perciformes)
# node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Equulites_elongatus", "Lutjanus_fulvus"))
# tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
# trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
# name_variable <- "only_ingroup"

###############################################################################################################################################
### Now subset the tree based on which traits you want to test ### 
###############################################################################################################################################

# For continous
trait.vector_c <- trait.data$diel_continuous
names(trait.vector_c) <- trait.data$species
trait.vector_c <- trait.vector_c[trait.vector_c %in% c(-5:5)]
trpy_c <- keep.tip(tr.calibrated, tip = names(trait.vector_c))
trpy_c$edge.length[trpy_c$edge.length == 0] <- 0.001


###############################################################################################################################################
### Run standard tests from ape ### 
###############################################################################################################################################

BM_cont <- ace(trait.vector_c, trpy_c, type = "continuous", model = "BM")

saveRDS(BM_cont, file = paste("best_fit_model_cont", name_variable, length(trpy_c$tip.label), "species.rds", sep = "_"))

BM_cont <- readRDS(file = paste("best_fit_model_cont", name_variable, length(trpy_c$tip.label), "species.rds", sep = "_"))

###############################################################################################################################################
### Plot the Hidden rate ancestral reconstructions (of rate and state) ### 
###############################################################################################################################################

lik.anc.cont <- data.frame(diel = c(trait.vector_c, BM_cont$ace), tips = c(names(trait.vector_c), names(BM_cont$ace)))
lik.anc.cont$node <- c(1:length(trpy_c$tip.label), (length(trpy_c$tip.label) + 1):(trpy_c$Nnode + length(trpy_c$tip.label)))


# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
BM_cont_plot <- ggtree(trpy_c, layout = "circular") %<+% lik.anc.cont + aes(color = diel) + geom_tippoint(aes(color = diel), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")


BM_cont_plot_rect <- ggtree(trpy_n) %<+% lik.anc.cont + aes(color = diel) + geom_tippoint(aes(color = diel), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")





###############################################################################################################################################
### Try the climate model OU ### 
###############################################################################################################################################


extinction <- read.csv("Speciation_paleo.csv", header = F)
colnames(extinction) <- c("Time", "speciation")

setwd("~/Downloads/R_code_data/")

# FIT OU CLIMATE MODEL FOR SCOTESE GLOBAL AVERAGE TEMP
source("fit_ou_process_simul_3.R")
source("fit_ou_env_v1.1.R")

## Fit a spline for the extinction rates
spline_result <- sm.spline(x = extinction[,1], y = extinction[,2], df=50)
env_func <- function(t){predict(spline_result,t)}
t <- unique(extinction[,1])
Extinction <- splinefun(t,env_func(t))
par(mar=c(5,5,3,5))

root_age <- nodeheight(trpy_c, node=1)
curve(Extinction, 0,root_age,col="blue", xlab="Age", ylab="Temperature C")


fun_temp <- function(x, par, theta0) theta0 + par*Extinction(root_age-x)

output <- fit_OU_trend(trpy_c, data = trait.vector_c, fun=fun_temp, startvalues = 0.1, method = "Nelder-Mead", control=list(maxit=2000), nuisance = TRUE)


output_2 <- fit_OU_trend(trpy_c, data = trait.vector_c, fun=fun_temp, startvalues = -1, method = "Nelder-Mead", control=list(maxit=2000), nuisance = TRUE)

options(mc.cores = 4)
fitCont <- fitContinuous(phy = trpy_c, dat = trait.vector_c, model = "OU", )

























## Compare to ER model for discrete

ER <- standard_tests[[1]]


lik.anc.dis <- trait.data_n[,c(3,3)]
colnames(lik.anc.dis) <- c("diurnal", "nocturnal")
lik.anc.dis$diurnal <- ifelse(lik.anc.dis$diurnal == "diurnal", 1, 0)
lik.anc.dis$nocturnal <- ifelse(lik.anc.dis$diurnal == "diurnal", 0, 1)
lik.anc.dis <- as.data.frame(rbind(lik.anc.dis, ER$lik.anc))

lik.anc.dis$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
ER_dis_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc.dis + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
ER_dis_plot_rect <- ggtree(trpy_n) %<+% lik.anc.dis + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")



comb <- cbind(lik.anc.cont[,c(1,2)], lik.anc.dis[,c(1,2)])

# Plot correlation between discrete and continuous
ggplot(comb, aes(x = diel, y = diurnal)) + geom_point() + theme_bw()

# Just reconstructed nodes
ggplot(comb[(length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)),], aes(x = diel, y = diurnal)) + geom_point() + theme_bw()


