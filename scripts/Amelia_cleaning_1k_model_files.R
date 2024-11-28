#Import and run packages -------------------------------------------------
# For retreiving data from the open tree of life
library(rotl)
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

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

#extra
library(ggplot2)
library(scales)
library(RColorBrewer)
#install.packages("ggforce")
library(ggforce)
library(forcats)
library(tidyr)

setwd(here())
source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")

# Set the working directory and source the functions (not used yet)
setwd("C:/Users/ameli/OneDrive/Documents/R_projects/fish_sleep/1k_model_results")

# Section 1: Three-state, cleaning and standardizing the 1k model results --------------

#we want to combine all the model results run on the cluster
#so that for each set of variables there is one file containing the ER, SYM, ARD and bridge only results

#new artiodactyla max crep three traits
ARD <- readRDS(here("artiodactyla_new_max_crep_traits_ARD_models.rds"))
ER_SYM_bridge_only <- readRDS(here("artiodactyla_new_max_crep_traits_ER_SYM_bridge_only_models.rds"))

all_models <- append(ARD, ER_SYM_bridge_only)

saveRDS(all_models, here("finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))

#new artiodactyla max dinoc three traits
ARD <- readRDS(here("artiodactyla_new_max_dinoc_traits_ARD_models.rds"))
ER_SYM_bridge_only <- readRDS(here("artiodactyla_new_max_dinoc_traits_ER_SYM_bridge_only_models.rds"))

all_models <- append(ARD, ER_SYM_bridge_only)

saveRDS(all_models, here("finalized_1k_models/artiodactyla_three_state_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))

#cetaceans max crep three traits
ARD <- readRDS(here("cetaceans_tree_max_crep_traits_ARD_models.rds"))
ER <- readRDS(here("cetaceans_tree_max_crep_traits_ER_models.rds"))
SYM <- readRDS(here("cetaceans_tree_max_crep_traits_SYM_models.rds"))
bridge_only <- readRDS(here("cetaceans_tree_max_crep_traits_bridge_only_models.rds"))

all_models <- append(ER, SYM)
all_models <- append(all_models, ARD)
all_models <- append(all_models, bridge_only)

saveRDS(all_models, here("finalized_1k_models/cetaceans_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))

#cetaceans max dinoc three traits
ARD <- readRDS(here("cetaceans_tree_max_dinoc_traits_ARD_models.rds"))
ER_SYM <- readRDS(here("cetaceans_tree_max_dinoc_traits_ER_SYM_models.rds"))
bridge_only <- readRDS(here("cetaceans_tree_max_dinoc_traits_bridge_only_models.rds"))

all_models <- append(ARD, ER_SYM)
all_models <- append(all_models, bridge_only)

saveRDS(all_models, here("finalized_1k_models/cetaceans_three_state_max_dinoc_traits_ER_SYM_ARD_bridge_only_models.rds"))

# Section 2: Six-state, cleaning and standardizing the 1k model results --------------

#artiodactyla six state, NOT DONE YET, NO CONSTRAINED MODEL

#ER_SYM_ARD<- readRDS(here("artiodactyla_six_state_traits_ER_SYM_ARD_models.rds"))
#bridge_only <- readRDS(here("artiodactyla_six_state_traits_bridge_only_models.rds"))
#all_models <- append(ER_SYM_ARD, bridge_only)

saveRDS(ER_SYM_ARD, here("finalized_1k_models/artiodactyla_six_state_ER_SYM_ARD_models.rds"))

#ruminants six state, NO CONSTRAINED MODEL

ER_SYM_ARD<- readRDS(here("ruminants_six_state_traits_ER_SYM_ARD_models.rds"))
#bridge_only <- readRDS(here("ruminants_six_state_traits_bridge_only_models.rds"))
#all_models <- append(ER_SYM_ARD, bridge_only)

saveRDS(ER_SYM_ARD, here("finalized_1k_models/ruminants_six_state_ER_SYM_ARD_models.rds"))


#whippomorpha six state, NO CONSTRAINED MODEL

ER_SYM_ARD<- readRDS(here("whippomorpha_six_state_traits_ER_SYM_ARD_models.rds"))
#bridge_only <- readRDS(here("whippomorpha_six_state_traits_bridge_only_models.rds"))
#all_models <- append(ER_SYM_ARD, bridge_only)

saveRDS(ER_SYM_ARD, here("finalized_1k_models/whippomorpha_six_state_ER_SYM_ARD_models.rds"))

# Section 3: Four state, cleaning and standardizing the 1k model --------
#whippomorpha
ER_SYM <- readRDS(here("whippomorpha_four_state_max_crep_traits_ER_SYM_models.rds"))
ARD_bridge_only <- readRDS(here("whippomorpha_four_state_max_crep_traits_ARD_bridge_only_models.rds"))
all_models <- append(ER_SYM, ARD_bridge_only)
saveRDS(all_models, here("finalized_1k_models/whippomorpha_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"))

#ruminants
ER_SYM <- readRDS(here("ruminants_four_state_max_crep_traits_ER_SYM_models.rds"))
ARD_bridge_only <- readRDS(here("ruminants_four_state_max_crep_traits_ARD_bridge_only_models.rds"))
all_models <- append(ER_SYM, ARD_bridge_only)
saveRDS(all_models, here("finalized_1k_models/ruminants_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"))

# Section 4: Four state, max clade cred model ---------------------------


