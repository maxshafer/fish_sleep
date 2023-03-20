library(phytools)
library(corHMM)
library(ggtree)
library(tidyr)
library(dplyr)
library(tibble)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/")

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/fish_sleep/scripts/fish_sleep_functions.R")

# Which tree?
name_variable <- "only_cartilaginous" # all, only_highqual, only_cartilaginous, or only_ingroup
dataset_variable <- "fish" # fish or AllGroups

## Load in the tree
trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c)
trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable)

## Load the best model, which is the HMM 2 state 2 rate model
models <- readRDS(file = paste("standard_tests", dataset_variable, name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
if (dataset_variable == "AllGroups") {
  model <- models[[2]]
} else {
  model <- readRDS(file = paste("best_fit_model", name_variable, length(trpy_n$tip.label), "species.rds", sep = "_"))
}


### Run simmap

simmaps <- makeSimmap(tree = trpy_n, data = trait.data_n[trpy_n$tip.label, c("species", "diel1")], rate.cat = 1, model = model$solution, nSim = 10, nCores = 4)


plotSimmap(simmaps[[5]], node.numbers = TRUE)



### This is a solution from phytools guy
### It works for simple cases, but not the back-and-forths I see
### I will have to first transform things into each other (both fast rate states as 1)
## Then run something like the following

tree <- simmaps[[5]]

states <- sort(unique(getStates(tree)))

colors <- setNames(palette()[1:length(states)],states)

obj <- get("last_plot.phylo", envir=.PlotPhyloEnv)

nc <- sapply(tree$maps,length)-1 # This is whenever there is more than 1 state per edge

ii <- which(nc>0) # which edges have transitions

nc <- nc[ii] # the number of transitions (if one, then 2 states, j and j+1)

h <- vector()

for(i in 1:length(ii)){
  for(j in 1:nc[i]){
    ss <- names(tree$maps[[ii[i]]])[j+1] # this returns the second state
    mm <- tree$edge[ii[i],1] # This returns the parental node for the edge
    dd <- tree$edge[ii[i],2] # descendant node
    
    x <- rep(obj$xx[mm]+cumsum(tree$maps[[ii[i]]])[j],2) # This is the x position of the transition (the age)
    # Instead of querying the plot ('obj$xx[mm]'), I can just query the node ages
    
    y<-c(obj$yy[dd]-0.5*mean(strheight(LETTERS)*cex),
         obj$yy[dd]+0.5*mean(strheight(LETTERS)*cex)) # I suppose this is the y position
    
    lines(x,y,lwd=lwd,col=colors[ss],lend=2)
    
    h<-c(h,x[1])
    
  }
}
invisible(h)



#### This is stuff from the paper https://www.nature.com/articles/s41467-023-36501-4.epdf?sharing_token=QfMaSEdzTqF5T6YrXpW1pNRgN0jAjWel9jnR3ZoTv0OUoPP6nlC3aeT4B3x29l9M9UoeMzYPnHHHYj8ZAMisGmVLEM_TKDrD1V848azrWE5wnd64ljbzUMAuqmWIDS5MneI61VcU-gi5yqmh72WY59Y-WFX82AP_OoeNsEIEQyY%3D

## https://github.com/stfriedman/Depth-transitions-paper/blob/main/Rcode/02_transition_analysis.R

### What this does is take the 'mapped_edge' slot of the simmap, which is the sum of the time spent in each state for each edge, and ask how many times
### for each descendant node of an mrca, there is more than 1 state on that edge. If more than 1 state, that means there was a transition on that edge
### This should work, but could fail if you expect transitions from A -> B -> C on a branch, as this would only be counted as 1 transition.
### At this point, it avoids double counting nodes, by only searching for each node in the "node_1" column

### So this works for my data, and I could adapt it so that it combines X1 & X3, and X2 & X4, to be clear it is correct  for example, only X1 and X3 might appear 
### on an edge, and be counted as a transition when it actually isn't). Easy to modify this to count transitions between rates as well
### However, it does not know when the transition took place, just that it took place on that branch (nothing in 'mapped_edge' lets you know which state was first, etc)
### For this, I would need to access the maps slot, which has the ordered list of states
### I could use the same framework, but then instead of looking at the 'mapped_edges', I use that to index 'maps', then the strings of fast rate transitions can be 'averaged'
### and I can take the midpoint + the age of node_1 as the transition time.


### grep the states in the 'maps' vector, take the mean(), and use that to order the df?

## These are some functions, they use purrr, not lapply

num_trans <- function(mrca, simtree) {
  trans_count <- function(desc) {
    edges <- get_mapped_edge(simtree) %>%
      mutate_at(vars(node_1:node_2), as.numeric)
    edges <- edges %>% mutate(index = row_number(), Di = X1, Noc = X2)

    
    df <- edges %>%
            filter(node_1 %in% desc) %>%
            select(-c(X1, X2)) %>%
            pivot_longer(-c(node_1, node_2, index), names_to = "state", values_to = "value") %>%
            filter(value != 0) %>%
            group_by(node_1, node_2) %>%
            add_count() %>%
            select(node_1, node_2, value, state, index, n) %>%
            filter(n != 1)
    
    ## At this point, I know that there was a transition, but not when (9.14 after 153, or 16.1 after 153?)
    ## If there is a transition, I can just ask if the first and last states differ,
    ## if so, then the time is node + the time spent in the original state. If they are the same, but this says there is a transition
    ## then there were two transitions?? How to deal with this?
    ## in this case, it starts at 4, ends as 4, but spends 9 million years as diurnal 3?
    
    grepl("1|3", names(simtree$maps[unique(df$index)][[1]])[1]) # not diurnal, therefore starts as noc
    
    
  }
  tibble(mrca) %>%
    mutate(
      desc = map(mrca, ~ c(.x, getDescendants(simtree, .x))),
      trans_num = map_dbl(desc, trans_count)
    )
}

get_mapped_edge <- function(x) {
  data.frame(x$mapped.edge) %>%
    rownames_to_column("edges") %>%
    as_tibble() %>%
    separate(edges, c("node_1", "node_2"), sep = ",")
}


## This is the actual command that returns the averaged transitions

clade_transitions <- tibble(tree_id = 1:length(simmaps)) %>%
  mutate(simmap = simmaps,
         trans = map(simmaps, ~num_trans(mrcas$mrca, .x))) %>%
  unnest(trans) %>%
  select(-simmap, -desc) %>%
  group_by(mrca) %>%
  summarize(trans_num = median(trans_num)) 





## OK, so 1 and 2 are the fast rate transitions

test <- simmaps[[1]]$maps[6][[1]]


test2 <- c(test)






## keep
num_trans <- function(mrca, simtree) {
  trans_count <- function(desc) {
    edges <- get_mapped_edge(simtree) %>%
      mutate_at(vars(node_1:node_2), as.numeric)
    
    edges %>%
      filter(node_1 %in% desc) %>%
      pivot_longer(-c(node_1, node_2), names_to = "state", values_to = "value") %>%
      filter(value != 0) %>%
      group_by(node_1, node_2) %>%
      add_count() %>%
      filter(n > 1) %>%
      select(node_1, node_2) %>%
      unique() %>%
      nrow(.)
  }
  tibble(mrca) %>%
    mutate(
      desc = map(mrca, ~ c(.x, getDescendants(simtree, .x))),
      trans_num = map_dbl(desc, trans_count)
    )
}