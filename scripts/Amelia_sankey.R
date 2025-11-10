# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
#install.packages("deeptime")
library(deeptime)
#update.packages("ggplot2")
library(ggplot2)
setwd(here())
library(networkD3)
library(tibble)

source("scripts/Amelia_functions.R")

# Section 1: Sankey tutorial ----------------------------------------------
#tutorial from https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html

#example data, includes both a source and target column
links <- data.frame(
  source=c("group_A","group_A", "group_B", "group_C", "group_C", "group_E"), 
  target=c("group_C","group_D", "group_E", "group_F", "group_G", "group_H"), 
  value=c(2,3, 2, 3, 1, 3)
)

#From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p

#we will use the mammals_df dataframe
links_mam <- mammals_df[, c("Maor_diel", "Bennie_diel")]
#we need to change this dataframe to include the number of species that start in one category and end in another
#ie 20 species start in diurnal and end in nocturnal

#

#fake data
links_mam <- data.frame(source = c("diurnal-A", "diurnal-A", "diurnal-A", "diurnal-A", "nocturnal-A", "nocturnal-A", "nocturnal-A", "nocturnal-A", "crepuscular-A", "crepuscular-A","crepuscular-A","crepuscular-A", "cathemeral-A", "cathemeral-A", "cathemeral-A", "cathemeral-A"), target = c("diurnal-B", "nocturnal-B", "crepuscular-B", "cathemeral-B"))
links_mam$value <- c(10, 5, 4, 17, 27, 37, 44, 25, 4, 13, 22, 25, 4, 6, 10, 15)

nodes_mam <- data.frame(name = c(as.character(links_mam$source), as.character(links_mam$target)) %>% unique())

links_mam$IDsource <- match(links_mam$source, nodes_mam$name)-1 
links_mam$IDtarget <- match(links_mam$target, nodes_mam$name)-1

p <- sankeyNetwork(Links = links_mam, Nodes = nodes_mam,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p

# Section 2: load in datasets ---------------------------------------------

#this is a merged dataframe of the Maor et al and Bennie et al diel activity dataframes
#only includes species that are present in both
mammals_df <- read.csv(here("sleepy_mammals.csv"))

#interested only in the species where there is also primary source data, meaning artiodactyls
artio_df <- read.csv(here("sleepy_artiodactyla_full.csv"))
#filter, leaves 318 species
artio_df <- artio_df[artio_df$Diel_Pattern_4 %in% c("diurnal", "nocturnal", "crepuscular", "cathemeral"), c("Diel_Pattern_2", "tips")]
colnames(artio_df) <- c("Amelia_diel", "tips")

#merge my artiodactyla data with the mammal data
mammals_df <- merge(mammals_df, artio_df, by = "tips") #leaves 149 species


# Section 3:Better sankey tutorial -----------------------------
#from https://r-charts.com/flow/sankey-diagram-ggplot2/ 
# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")

library(ggsankey)

#ggsankey comes with a function make_long that converts your data into the correct format
df <- mtcars %>% make_long(cyl, vs, am, gear, carb)

#creates plot based on the starting node and end node
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) + geom_sankey() + theme_sankey(base_size = 16)

#try with our data 
df <- mammals_df %>% make_long(Bennie_diel, Maor_diel)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) + geom_sankey() + theme_sankey(base_size = 16)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha= 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option = "D", alpha = 0.95) + theme_sankey(base_size = 16) + guides(fill = guide_legend(title = "Temporal activity pattern")) +
  theme(legend.position = "none") + labs(x = NULL)

#move my data to centre so its easier to compare my data to both existing datasets
mammals_df <- mammals_df %>% relocate(Maor_diel, .after = last_col())

#try with three data sources  
mammals_df$Amelia_diel <- str_replace(mammals_df$Amelia_diel, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
df <- mammals_df %>% make_long(Bennie_diel, Amelia_diel, Maor_diel,)
test <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha= 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option = "D", alpha = 0.95) + theme_sankey(base_size = 16) + guides(fill = guide_legend(title = "Temporal activity pattern")) +
  theme(legend.position = "none") + labs(x = NULL)

#save out to figure folder
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "sankey_all.pdf"))
test
dev.off()

#max crep
mammals_df$max_crep <- mammals_df$Amelia_diel
mammals_df$max_crep <- str_replace(mammals_df$max_crep, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
mammals_df$max_crep <- str_replace(mammals_df$max_crep, pattern = "diurnal/crepuscular", replacement = "crepuscular")
mammals_df$max_crep <- str_replace(mammals_df$max_crep, pattern = "nocturnal/crepuscular", replacement = "crepuscular")


df <- mammals_df %>% make_long(Bennie_diel, max_crep, Maor_diel)
max_crep <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha= 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option = "D", alpha = 0.95) + theme_sankey(base_size = 16) + guides(fill = guide_legend(title = "Temporal activity pattern")) +
  theme(legend.position = "none") + labs(x = NULL)

#save out to figure folder
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "sankey_max_crep.pdf"))
max_crep
dev.off()

#max dinoc
mammals_df$max_dinoc <- mammals_df$Amelia_diel
mammals_df$max_dinoc <- str_replace(mammals_df$max_dinoc, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
mammals_df$max_dinoc <- str_replace(mammals_df$max_dinoc, pattern = "diurnal/crepuscular", replacement = "diurnal")
mammals_df$max_dinoc <- str_replace(mammals_df$max_dinoc, pattern = "nocturnal/crepuscular", replacement = "nocturnal")

df <- mammals_df %>% make_long(Bennie_diel, max_dinoc, Maor_diel)
max_dinoc <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha= 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option = "D", alpha = 0.95) + theme_sankey(base_size = 16) + guides(fill = guide_legend(title = "Temporal activity pattern")) +
  theme(legend.position = "none") + labs(x = NULL)

#save out to figure folder
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "sankey_max_dinoc.pdf"))
max_dinoc
dev.off()