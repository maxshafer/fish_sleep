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

# # Remove low quality data
# trait.data <- trait.data[trait.data$confidence > 1,]
# name_variable <- "only_highqual"

# # Remove cartilaginous fish (just actinopterygii)
# # Common ancestor of Lepisosteus_osseus (Order: Lepisosteiformes) Lepidosiren_paradoxa (lepidosiren) and Lutjanus_fulvus (Order: Perciformes)
# node_of_interest <- getMRCA(phy = tr.calibrated, tip = c("Lepidosiren_paradoxa", "Lutjanus_fulvus"))
# tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
# trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
# name_variable <- "only_ingroup"

# Keep only cartilaginous fish (no actinopterygii)
node_of_interest <- getMRCA(tr.calibrated, tip = c("Rhizoprionodon_terraenovae", "Rhynchobatus_djiddensis"))
tr.calibrated <- extract.clade(phy = tr.calibrated, node = node_of_interest)
trait.data <- trait.data[trait.data$species %in% tr.calibrated$tip.label,]
name_variable <- "only_cartilaginous"

###############################################################################################################################################
### Now subset the tree based on which traits you want to test ### 
###############################################################################################################################################

# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]

# Diurnal/nocturnal plus marine/fresh
trait.vector_m <- paste(trait.data$diel1, trait.data$marine, sep = "_")
names(trait.vector_m) <- trait.data$species
trait.vector_m <- trait.vector_m[trait.vector_m %in% c("diurnal_freshwater", "diurnal_saltwater", "nocturnal_freshwater", "nocturnal_saltwater")]
trpy_m <- keep.tip(tr.calibrated, tip = names(trait.vector_m))
trpy_m$edge.length[trpy_m$edge.length == 0] <- 0.001 
trait.data_m <- trait.data[trait.data$species %in% trpy_m$tip.label,]


###############################################################################################################################################
### Run standard tests from ape ### 
###############################################################################################################################################

standard_tests <- list()
# Equal rates, symmetric rates (same as ER), and All rates different
standard_tests[[1]] <- ace(trait.vector_n, trpy_n, model = "ER", type = "discrete")
standard_tests[[2]] <- ace(trait.vector_n, trpy_n, model = "SYM", type = "discrete")
standard_tests[[3]] <- ace(trait.vector_n, trpy_n, model = "ARD", type = "discrete")

# Equal rates, symmetric rates (same as ER), and All rates different WITH marine
standard_tests[[4]] <- ace(trait.vector_m, trpy_m, model = "ER", type = "discrete")
standard_tests[[5]] <- ace(trait.vector_m, trpy_m, model = "SYM", type = "discrete")
standard_tests[[6]] <- ace(trait.vector_m, trpy_m, model = "ARD", type = "discrete")


names(standard_tests) <- c("fitER", "fitSYM", "fitARD", "fitERmarine", "fitSYMmarine", "fitARDmarine")


# Symmetric rates, without a rate for switching both at once
standard_tests[[7]] <- ace(trait.vector_m, trpy_m, model = matrix(c(0,1,2,0, 1,0,0,3, 2,0,0,4, 0,3,4,0), 4), type = "discrete")
# ARD rates, without a rate for switching both at once
standard_tests[[8]] <- ace(trait.vector_m, trpy_m, model = matrix(c(0,1,3,0, 2,0,0,5, 4,0,0,7, 0,6,8,0), 4), type = "discrete")

# Only rates for M <-> F or D <-> N
trait.vector_m2 <- str_replace(trait.vector_m, "saltwater", "freshwater")
standard_tests[[9]] <- ace(trait.vector_m2, trpy_m, model = "ER", type = "discrete")
standard_tests[[10]] <- ace(trait.vector_m2, trpy_m, model = "ARD", type = "discrete")

trait.vector_m2 <- str_replace(trait.vector_m, "nocturnal", "diurnal")
standard_tests[[11]] <- ace(trait.vector_m2, trpy_m, model = "ER", type = "discrete")
standard_tests[[12]] <- ace(trait.vector_m2, trpy_m, model = "ARD", type = "discrete")

names(standard_tests) <- c("fitER", "fitSYM", "fitARD", "fitERmarine", "fitSYMmarine", "fitARDmarine", "fitSINGLE.SYM", "fitSINGLE.ARD", "fitDNonlyER", "fitDNonlyARD", "fitMFonlyER", "fitMFonlyARD")

unlist(lapply(standard_tests, function(x) x[names(x[grep("loglik", names(x))])]))

View(unlist(lapply(standard_tests, function(x) x[names(x[grep("loglik", names(x))])])))

# Very clear from the above results, that adding in Marine doesn't help the ancestral reconstruction
# lowest loglik values are when you only consider Diurnal/Nocturnal
# Lower also if you remove low quality data

saveRDS(standard_tests, file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

###############################################################################################################################################
### Run Hidden Rates Models ### 
###############################################################################################################################################

# Run corHMM for just diel or diel and marine, with 1 rat category (Markov model)
MK_2state <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = "ARD")
MK_4state_m <- corHMM(phy = trpy_m, data = trait.data_m[trpy_m$tip.label, c("species", "diel1", "marine")], rate.cat = 1, model = "ARD")

# Run corHMM for just diel, with 2 or 3 rate categories (large improvement for 2, diminishing returns for 3)
HMM_2state_2rate <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, model = "ARD", get.tip.states = TRUE)
HMM_2state_3rate <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 3, model = "ARD", get.tip.states = TRUE)
HMM_2state_4rate <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 4, model = "ARD", get.tip.states = TRUE)

# # Run corHMM for diel x marine or diel x acanthomorpha, with 2 or 3 rate categories
# HMM_4state_2rate_m <- corHMM(phy = trpy_m, data = trait.data_m[trpy_m$tip.label, ], rate.cat = 2, model = "ARD", get.tip.states = TRUE)
# HMM_4state_3rate_m <- corHMM(phy = trpy_m, data = trait.data_m[trpy_m$tip.label, ], rate.cat = 3, model = "ARD", get.tip.states = TRUE)

###############################################################################################################################################

# Let's do a 2state_2rate model, where once you transition into a higher rate class, you can't switch back (testing the 'origin' of fast-switching)
# Let's say R2 cannot transition to R1, but can go D->N and N->D, R1 can go to R2, and D->N and N->D
# NoGoingBack (NGB) model

DN_LegendAndRate <- getStateMat4Dat(trait.data_n[trpy_n$tip.label, c("species", "diel1")])
DN_R1 <- DN_LegendAndRate$rate.mat
DN_R2 <- DN_LegendAndRate$rate.mat

DN_ObsStateClasses <- list(DN_R1, DN_R2)
DN_RateClassMat <- getRateCatMat(2)
DN_RateClassMat <- dropStateMatPars(DN_RateClassMat, c(1))

DN_FullMat <- getFullMat(DN_ObsStateClasses, DN_RateClassMat)
# plotMKmodel(DN_FullMat, rate.cat = 2, display = "square", text.scale = 0.9)

HMM_2state_2rate_NGB <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, rate.mat = DN_FullMat)

###############################################################################################################################################

# 2state_2rate model, where R1 cannot switch D<->N, but R2 can
# NoSwitchRate (NSR) model

DN_LegendAndRate <- getStateMat4Dat(trait.data_n[trpy_n$tip.label, c("species", "diel1")])
DN_R1 <- dropStateMatPars(DN_LegendAndRate$rate.mat, c(1,2))
DN_R2 <- DN_LegendAndRate$rate.mat

DN_ObsStateClasses <- list(DN_R1, DN_R2)
DN_RateClassMat <- getRateCatMat(2)

DN_FullMat <- getFullMat(DN_ObsStateClasses, DN_RateClassMat)
# plotMKmodel(DN_FullMat, rate.cat = 2, display = "square", text.scale = 0.9)

HMM_2state_2rate_NSR <- corHMM(phy = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 2, rate.mat = DN_FullMat, node.states)

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# 2 state 2 rate is best fit that still makes sense (3 rates starts to maybe overfit)
standard_tests <- append(standard_tests, list(MK_2state, MK_4state_m, HMM_2state_2rate, HMM_2state_3rate, HMM_2state_4rate, HMM_2state_2rate_NGB, HMM_2state_2rate_NSR))
names(standard_tests) <- c("fitER", "fitSYM", "fitARD", "fitERmarine", "fitSYMmarine", "fitARDmarine", "fitSINGLE.SYM", "fitSINGLE.ARD", "fitDNonlyER", "fitDNonlyARD", "fitMFonlyER", "fitMFonlyARD", "MK_2state", "MK_4state_m", "HMM_2state_2rate", "HMM_2state_3rate", "HMM_2state_4rate", "HMM_2state_2rate_NGB", "HMM_2state_2rate_NSR")
View(unlist(lapply(standard_tests, function(x) x[names(x[grep("loglik", names(x))])])))
View(unlist(lapply(standard_tests, function(x) x[names(x[grep("AIC", names(x))])])))

saveRDS(standard_tests, file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

saveRDS(HMM_2state_2rate, file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

# Look at models
plotMKmodel(MK_2state) # High Diurnal -> Nocturnal
plotMKmodel(MK_4state_m)
plotMKmodel(HMM_2state_2rate) # 1 rate is high transition rate (both ways), other is low both ways, with slightly higher rate to go from high to low transition rate
plotMKmodel(HMM_2state_2rate_NGB)
plotMKmodel(HMM_2state_2rate_NSR)
plotMKmodel(HMM_2state_3rate) # 1) 3x higher Noc->Di (high bow), 2) 3x higher Di-> Noc (medium both), 3) Low both (less Di->Noc) - low transition rates all around
plotMKmodel(HMM_2state_4rate)

###############################################################################################################################################

# Re-load the best fit model
# HMM_2state_2rate <- readRDS(file = paste("best_fit_model_", length(trpy_n$tip.label), "_species.rds", sep = ""))
HMM_2state_2rate <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))

###############################################################################################################################################
### Plot the Hidden rate ancestral reconstructions (of rate and state) ### 
###############################################################################################################################################

model <- HMM_2state_2rate

lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2")
# colnames(lik.anc) <- c("diurnal_extinction", "diurnal_radiation", "nocturnal_extinction", "nocturnal_radiation")
lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

lik.anc$noc.sum <- rowSums(lik.anc[,c(2,4)])
lik.anc$di.sum <- rowSums(lik.anc[,c(1,3)])

lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])

lik.anc$noc.sum <- rowSums(lik.anc[,c(3,4)])
lik.anc$di.sum <- rowSums(lik.anc[,c(1,2)])

lik.anc$Ext.sum <- rowSums(lik.anc[,c(1,3)])
lik.anc$Rad.sum <- rowSums(lik.anc[,c(2,4)])

# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
diel_2state_2rate <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
diel_2state_2rate_rect <- ggtree(trpy_n) %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")

rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")


dir1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R1) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1)
dir2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R2) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
nocr1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R1) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)
nocr2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R2) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)
dir3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R3) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "YlOrBr", direction = 1)
nocr3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R3) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "BuPu", direction = 1)

two_rate_plots <- dir1 + dir2 + nocr1 + nocr2 + plot_layout(nrow = 2)

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_diel_rates", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width =20, height = 20)
two_rate_plots
dev.off()

pdf(file = paste("outs/Figures/fish_phylogeny_diel_ancestral_rates", name_variable, length(trpy_n$tip.label), "species.pdf", sep = "_"), width = 20, height = 20)
rate_plot
dev.off()


###############################################################################################################################################
###############################################################################################################################################

# lik.anc <- as.data.frame(rbind(HMM_2state_4rate$tip.states, HMM_2state_4rate$states))
# colnames(lik.anc) <- c("diurnal_R1", "nocturnal_R1", "diurnal_R2", "nocturnal_R2", "diurnal_R3", "nocturnal_R3", "diurnal_R4", "nocturnal_R4")
# lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))
# 
# lik.anc$noc.sum <- rowSums(lik.anc[,c(2,4,6,8)])
# lik.anc$di.sum <- rowSums(lik.anc[,c(1,3,5,7)])
# 
# lik.anc$R1.sum <- rowSums(lik.anc[,c(1,2)])
# lik.anc$R2.sum <- rowSums(lik.anc[,c(3,4)])
# lik.anc$R3.sum <- rowSums(lik.anc[,c(5,6)])
# lik.anc$R4.sum <- rowSums(lik.anc[,c(7,8)])
# 
# # Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
# diel_2state_3rate <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
# 
# r1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R1.sum) + geom_tippoint(aes(color = R1.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1) + scale_color_distiller(palette = "Reds", direction = 1)
# r2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R2.sum) + geom_tippoint(aes(color = R2.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1) + scale_color_distiller(palette = "OrRd", direction = 1)
# r3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R3.sum) + geom_tippoint(aes(color = R3.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1) + scale_color_distiller(palette = "OrRd", direction = 1)
# r4 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = R4.sum) + geom_tippoint(aes(color = R3.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1) + scale_color_distiller(palette = "OrRd", direction = 1)
# 
# r1 + r2 + r3 + r4 + plot_layout(nrow = 1)
# 
# dir1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R1) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1)
# dir2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R2) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
# dir3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = diurnal_R3) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1)
# nocr1 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R1) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_R1), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)
# nocr2 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R2) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R2), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)
# nocr3 <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = nocturnal_R3) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_R3), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1)
# 
# three_rate_plots <- dir1 + dir2 + dir3 + nocr1+ nocr2 + nocr3 + plot_layout(nrow = 2)

###############################################################################################################################################
###############################################################################################################################################














# get simmap inputs from corhmm outputs
phy = HMM_2state_3rate$phy
data = HMM_2state_3rate$data
model = HMM_2state_3rate$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
# run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree = trpy, data = trait.data[,c(1,2)], model = model, rate.cat = 3, nSim = 1, nCores = 1)
# we import phytools plotSimmap for plotting
cols<-setNames(c("green","#E4D96F","darkgreen",
                 "brown","black","darkgrey"),
               c("CG","GB","TC","TG","Tr","Tw"))
phytools::plotSimmap(simmap[[1]], fsize = 0.5, colors = setNames(c("#8B0000", "#00008B", "#FF0000", "#0000FF", "#FFC0CB", "#ADD8E6"), c("1","2", "3", "4", "5", "6")))


test <- ggtree(simmap[[1]], layout = "circular")




















library(ape)
library(geiger)
library(phytools)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

## Load files for making figures

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)

trpy <- readRDS("calibrated_phylo_tree_ancestral.rds") 


trait.vector <- resolved_names$diel[match(trpy$tip.label, resolved_names$tips)]
trait.vector <- str_replace(trait.vector, "crepuscular/", "")
names(trait.vector) <- trpy$tip.label

trait.vector <- str_replace(trait.vector_n, "_saltwater", "")
trait.vector <- str_replace(trait.vector, "_freshwater", "")
names(trait.vector) <- names(trait.vector_n)
fish_diel <- treedata(trpy_n, trait.vector)
fish_diel_phy <- fish_diel$phy
fish_diel_dat <- as.character(fish_diel$data)
names(fish_diel_dat) <- rownames(fish_diel$data)



ngen <- 5e+05 # number of generations to run
burnin <- 0.2 * ngen # approximate number of initial generations to eliminate
sample <- 500 # sample posterior distribution every 1000 samples


mcmc_diel <- ancThresh(fish_diel_phy, fish_diel_dat[1:length(fish_diel_phy$tip.label)], ngen = ngen, sequence = c("nocturnal", "diurnal"), model = "lambda", control = list(sample = sample, plot = TRUE))

plot(mcmc_diel$par[, "logLik"], type = "l", xlab = "generation", ylab = "logLik")

colMeans(mcmc_diel$par[(0.2 * ngen/sample):(ngen/sample) + 1, c("nocturnal", "diurnal")])

matrix_diel_colors <-to.matrix(fish_diel_dat, c("nocturnal", "diurnal"))


phylogeny_name_order <- as.data.frame(fish_diel_phy$tip.label) # creates a data frame for species names in correct order


names(phylogeny_name_order) <- c("species_names") # changes the name of the column
rownames(phylogeny_name_order) <- phylogeny_name_order$species_names # assigns the namesto rows

matrix_diel_colors_ordered <- matrix_diel_colors[rownames(phylogeny_name_order),,drop=FALSE] # reorders matrix names


plotTree(fish_diel_phy, type = "fan", setEnv = TRUE, ftype = "off", direction = "downwards")
tiplabels(pie = matrix_diel_colors_ordered, piecol = palette()[1:3], cex = 0.2)
nodelabels(pie = mcmc_diel$ace, piecol = palette()[1:3], cex = 0.2)


saveRDS(mcmc_diel, file = "mcmc_diel_500kgen.rds")
