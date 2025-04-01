##Packages we will use ---------------------------------------------------
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)
#open tree of life
library(rotl)
#adds timescale
library(deeptime)
#colours
library(RColorBrewer)
#apply two separate colour palettes
library(ggnewscale)
#more colours
library(pals)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")



# Section 1: Plot data on tree --------------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

diel_full <- diel_full %>% filter(Parvorder == "Mysticeti")

diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label, -c(2, 15:21)]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_3", "New_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_3), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") 
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "New_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/new_cetacean_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()