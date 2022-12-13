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


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/Fish_sleep_functions.R")

## Load files

resolved_names <- readRDS("resolved_names_AllGroups.rds")
tr.calibrated <- readRDS("tr_tree_calibrated_AllGroups.rds")
trait.data <- readRDS("trait_data_AllGroups.rds")

name_variable <- "AllGroups"


# Just diurnal/nocturnal
trait.vector_n <- trait.data$diel1
names(trait.vector_n) <- trait.data$species
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("diurnal", "nocturnal")]
trpy_n <- keep.tip(tr.calibrated, tip = names(trait.vector_n))
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]


# Load the best model, which is the HMM 2 state 2 rate model
models <- readRDS(file = paste("standard_tests", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
model <- models[[3]]

###############################################################################################################################################
### Create the object with lineages through time and diel switches data ### 
###############################################################################################################################################

# First, extract the ancestral states from the best fit model
anc_states <- returnAncestralStates(phylo_model = model, phylo_tree = trpy_n)

# Then, calculate transitions between states
anc_states <- calculateStateTranstitions(ancestral_states = anc_states, phylo_tree = trpy_n)

# Determine transition histories (types of lineages)
anc_states <- calculateLinTransHist(ancestral_states = anc_states, phylo_tree = trpy_n)

# Calculate cumsums through time (for ltt plots)
anc_states <- returnCumSums(ancestral_states = anc_states, phylo_tree = trpy_n)


###############################################################################################################################################
### Make Plots! ### 
###############################################################################################################################################

switch.histo <- switchHisto(ancestral_states = anc_states, replace_variable_names = TRUE)

switch.ratio <- switchRatio(ancestral_states = anc_states, phylo_tree = trpy_n, use_ltt = F)







node.data3 <- gather(anc_states$cumsums, "variable", "Cummulative_ratio", -node.age, -Lineage_Cumsum)

ltt.diel.switch.dir.plot <- ggplot(node.data3, aes(x = node.age, y = Cummulative_ratio/Lineage_Cumsum, colour = variable)) + geom_line(size = 1.5) + theme_classic() #+ scale_colour_manual(values = c("switch.ND.Cum" = "#d6604d", "switch.DN.Cum" = "#4393c3"))
ltt.diel.switch.dir.plot <- ltt.diel.switch.dir.plot + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10)) + scale_x_reverse() + xlab("Millions of years ago") # + ylab("log(lineages-through-time)") 

node.data3$CRLC <- node.data3$Cummulative_ratio#/node.data3$Lineage_Cumsum
node.data3 <- node.data3 %>% group_by(node.age, variable) %>% summarise(n = sum(CRLC)) %>% mutate(percentage = n / sum(n))

## Should come up with function to make this? Yeah, and auto-find colors maybe from the red and blue scales
switch_histo <- ggplot(node.data3, aes(x = node.age, y = percentage, fill = variable)) + theme_classic() + geom_area() + scale_x_reverse() + xlab("Millions of years ago") #+ scale_fill_manual(values = c("blue4", "red4", "blue3", "red3", "blue1", "red1"))


# Plots
# switch_histo <- ggplot(lik.anc, aes(x = node.age, fill = as.character(switch.type))) + geom_histogram(position = "fill", bins = 100) + scale_x_reverse()
switch_ratio <- ggplot(lik.anc, aes(x = node.age, y = switch2/Lineage_Cumsum2)) + geom_line() + theme_classic() + scale_x_reverse() # Lineage_Cumsum2 would be the number of times it was possible to switch? Or something like that?
ltt.switch.plot <- ggplot(lik.anc, aes(x = node.age, y = log(switch2))) + geom_line(size = 1.5) + theme_classic() + scale_x_reverse() + ylab("log(switches-through-time)") + xlab("Millions of years ago")




