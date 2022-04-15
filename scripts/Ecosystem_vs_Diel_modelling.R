library(rfishbase)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/scripts/Trait_rate_matrix_figure_script.R")

## Load files

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
tr.calibrated <- readRDS("calibrated_phylo.rds")
trait.data <- readRDS("trait_data.rds")

## Load fishbase data

fishbase_species <- rfishbase::load_taxa()

fishbase_diet <- diet()
fishbase_ecology <- ecology() # Has troph and 
fishbase_morph <- morphometrics()
fishbase_ecosystem <- ecosystem() # Salinity is here duh



#resolved_names$marine <- fishbase_ecosystem$Salinity[match(resolved_names$unique_name, fishbase_ecosystem$Species)]

### Compare marine and diel together as discrete traits (4 of them)


# First create a vector for the trait (diurnal or nocturnal + marine)
trait.vector <- paste(resolved_names$diel[match(trpy$tip.label, resolved_names$tips)], resolved_names$marine[match(trpy$tip.label, resolved_names$tips)], sep = "_")
trait.vector <- str_replace(trait.vector, "crepuscular/", "")
names(trait.vector) <- trpy$tip.label
# Only those with diurnal/nocturnal or fresh/salt water
trait.vector.marine <- trait.vector[unique(c(grep("diurnal_freshwater", trait.vector), grep("diurnal_saltwater", trait.vector), grep("nocturnal_freshwater", trait.vector), grep("nocturnal_saltwater", trait.vector)))]

# diel.vector2 <- vector(length = length(trait.vector))
# diel.vector2[grep("diurnal", trait.vector)] <- "diurnal"
# diel.vector2[grep("nocturnal", trait.vector)] <- "nocturnal"
# 
# marine.vector2 <- vector(length = length(trait.vector))
# marine.vector2[grep("saltwater", trait.vector)] <- "saltwater"
# marine.vector2[grep("freshwater", trait.vector)] <- "freshwater"


# # First create a vector for the trait (diurnal or nocturnal + climate)
# trait.vector <- paste(resolved_names$diel[match(trpy$tip.label, resolved_names$tips)], resolved_names$Climate[match(trpy$tip.label, resolved_names$tips)], sep = "_")
# trait.vector <- str_replace(trait.vector, "crepuscular/", "")
# names(trait.vector) <- trpy$tip.label
# trait.vector.climate <- trait.vector[unique(c(grep("diurnal_subtropical", trait.vector), grep("diurnal_temperate", trait.vector), grep("diurnal_tropical", trait.vector), grep("diurnal_boreal", trait.vector), grep("nocturnal_subtropical", trait.vector), grep("nocturnal_temperate", trait.vector), grep("nocturnal_tropical", trait.vector), grep("nocturnal_boreal", trait.vector)))]

# fit an ARD based Mk model for discrete character evolution (allows different forward and backward rates), which returns ancestral liklihoods for each node
# ARD model has the lowest log-likelihood compared with other models for discrete traits, so it fits best

trait.vector_n <- trait.vector.marine # Specify which trait vector to use below
trpy_n <- keep.tip(trpy, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 


fitER <- ace(trait.vector_n, trpy_n, model = "ER", type = "discrete")
fitSYM <- ace(trait.vector_n, trpy_n, model = "SYM", type = "discrete")
fitARD <- ace(trait.vector_n, trpy_n, model = "ARD", type = "discrete")

plotRateIndex(model = fitSYM)

# ace(diel.vector2, trpy_n, model = "ER", type = "discrete")
# ace(diel.vector2, trpy_n, model = "SYM", type = "discrete")
# ace(diel.vector2, trpy_n, model = "ARD", type = "discrete")


# Symmetric rates, without a rate for switching both at once
fitSINGLE.SYM <- ace(trait.vector_n, trpy_n, model = matrix(c(0,1,2,0, 1,0,0,3, 2,0,0,4, 0,3,4,0), 4), type = "discrete")
# ARD rates, without a rate for switching both at once
fitSINGLE.ARD <- ace(trait.vector_n, trpy_n, model = matrix(c(0,1,3,0, 2,0,0,5, 4,0,0,7, 0,6,8,0), 4), type = "discrete")

# Only rates for M <-> F or D <-> N
fitMFonly <- ace(trait.vector_n, trpy_n, model = matrix(c(0,2,0,2, 1,0,1,0, 0,2,0,2, 1,0,1,0), 4), type = "discrete")
fitDNonly <- ace(trait.vector_n, trpy_n, model = matrix(c(0,0,2,2, 0,0,2,2, 1,1,0,0, 1,1,0,0), 4), type = "discrete")


plotRateIndex(model = fitMFonly)


tv <- data.frame(tips = names(trait.vector_n), diel_marine = trait.vector_n)
diel.plot <- ggtree(trpy2, layout = "circular") %<+% tv + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diel_marine), shape = 16, size = 1) + scale_color_manual(values = c("red1", "red4", "blue1", "blue4"))

lik.anc2 <- as.data.frame(fitSYM$lik.anc)
lik.anc2$node <- (1:nrow(lik.anc2)) + length(trpy_n$tip.label)

node.data2 <- data.frame(diurnal_freshwater = ifelse(trait.vector_n == "diurnal_freshwater", 1, 0), diurnal_saltwater = ifelse(trait.vector_n == "diurnal_saltwater", 1, 0), nocturnal_freshwater = ifelse(trait.vector_n == "nocturnal_freshwater", 1, 0), nocturnal_saltwater = ifelse(trait.vector_n == "nocturnal_saltwater", 1, 0), node = 1:length(trpy_n$tip.label))
# node.data <- data.frame(diurnal = ifelse(diel.vector == "diurnal", 1, 0), nocturnal = ifelse(diel.vector == "diurnal", 0, 1), unclear = ifelse(diel.vector == "unclear", 0, 1), crepuscular = ifelse(diel.vector == "crepuscular", 0, 1), node = 1:length(trpy$tip.label))
node.data2 <- rbind(node.data2, lik.anc2)


rate_index <- plotRateIndex(model = fitSYM)

di_fresh <- ggtree(trpy2, layout = "circular") %<+% node.data2 + aes(color = diurnal_freshwater) + scale_color_distiller(palette = "Reds", direction = 1) + geom_tippoint(aes(color = diurnal_freshwater), shape = 16, size = 1.5) + scale_color_distiller(palette = "Reds", direction = 1) + ggtitle("Diurnal Freshwater")
di_salt <- ggtree(trpy2, layout = "circular") %<+% node.data2 + aes(color = diurnal_saltwater) + scale_color_distiller(palette = "OrRd", direction = 1) + geom_tippoint(aes(color = diurnal_saltwater), shape = 16, size = 1.5) + scale_color_distiller(palette = "OrRd", direction = 1) + ggtitle("Diurnal Saltwater")

noc_fresh <- ggtree(trpy2, layout = "circular") %<+% node.data2 + aes(color = nocturnal_freshwater) + scale_color_distiller(palette = "Blues", direction = 1) + geom_tippoint(aes(color = nocturnal_freshwater), shape = 16, size = 1.5) + scale_color_distiller(palette = "Blues", direction = 1)  + ggtitle("Nocturnal Freshwater")
noc_salt <- ggtree(trpy2, layout = "circular") %<+% node.data2 + aes(color = nocturnal_saltwater) + scale_color_distiller(palette = "Purples", direction = 1) + geom_tippoint(aes(color = nocturnal_saltwater), shape = 16, size = 1.5) + scale_color_distiller(palette = "Purples", direction = 1) + ggtitle("Nocturnal Saltwater")

di_fresh + di_salt + rate_index + noc_fresh + noc_salt + plot_layout(nrow =2, ncol = 3)














