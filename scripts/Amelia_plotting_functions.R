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


#load in mammal tree and cetacean dataframe
mammal_trees <- read.nexus(here("Cox_mammal_data/Complete_phylogeny.nex"))
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))




# #Function 1: max_clade_cred likelihood metrics --------------------------


max_clade_metrics <- function(model_results = readRDS(here("test_whippomorpha_max_clade_cred_four_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds"))) {
  # will take take the max clade cred tree result and plot the likelihood, AIC, and AICc score
  log_likelihoods <- unlist(lapply(model_results, function(x) returnLikelihoods(model = x)))
  likelihoods <- as.data.frame(log_likelihoods)
  likelihoods$modelname <- rownames(likelihoods)
  likelihoods <- separate(likelihoods, col = "modelname", into = c("model", "to_drop"), sep = "[.]", remove = TRUE)
  AICc_scores <- unlist(lapply(model_results, function(x) returnAICc(model = x)))
  likelihoods$AICc_scores <-AICc_scores 
  AIC_scores <- unlist(lapply(model_results, function(x) returnAIC(model = x)))
  likelihoods$AIC_scores <-AIC_scores 
  likelihoods <- likelihoods[,-3]
  likelihoods <- likelihoods %>% pivot_longer(!model, names_to = "model_metric", values_to = "model_value")
  
  return(likelihoods)
}


# # Function 2: Likelihoods from 1k model results (ER, SYM, ARD, bridge_only) -------------------------

plot1kLikelihoods <- function(model_results = readRDS(here("finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds")), number_of_models = 4) {
  # extract likelihoods
  if(number_of_models == 3){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
    
    df1 <- data.frame(model = "ER", likelihoods = ER_likelihoods)
    df2 <- data.frame(model = "SYM", likelihoods = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", likelihoods = ARD_likelihoods)
    df_full <- rbind(df1, df2, df3)
  }
  
  if(number_of_models == 4){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
    bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
    
    df1 <- data.frame(model = "ER", likelihoods = ER_likelihoods)
    df2 <- data.frame(model = "SYM", likelihoods = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", likelihoods = ARD_likelihoods)
    df4 <- data.frame(model = "bridge_only", likelihoods = bridge_only_likelihoods)
    df_full <- rbind(df1, df2, df3, df4)
  }
  
  return(df_full)
}

# # Function 3: AIC scores from 1k model results -------------------------
plot1kAIC <- function(model_results = readRDS(here("finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds")), number_of_models = 4) {
  # extract likelihoods
  if(number_of_models == 3){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnAIC(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnAIC(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnAIC(model = x)))
    
    df1 <- data.frame(model = "ER", AIC_score = ER_likelihoods)
    df2 <- data.frame(model = "SYM", AIC_score = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", AIC_score = ARD_likelihoods)
    df_full <- rbind(df1, df2, df3)
  }
  
  if(number_of_models == 4){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnAIC(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnAIC(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnAIC(model = x)))
    bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnAIC(model = x)))
    
    df1 <- data.frame(model = "ER", AIC_score = ER_likelihoods)
    df2 <- data.frame(model = "SYM", AIC_score = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", AIC_score = ARD_likelihoods)
    df4 <- data.frame(model = "bridge_only", AIC_score = bridge_only_likelihoods)
    df_full <- rbind(df1, df2, df3, df4)
  }
  
  return(df_full)
}

# # Function 4: AICc scores from 1k model results -------------------------
plot1kAICc <- function(model_results = readRDS(here("finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds")), number_of_models = 4) {
  # extract likelihoods
  if(number_of_models == 3){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnAICc(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnAICc(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnAICc(model = x)))
    
    df1 <- data.frame(model = "ER", AICc_score = ER_likelihoods)
    df2 <- data.frame(model = "SYM", AICc_score = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", AICc_score = ARD_likelihoods)
    df_full <- rbind(df1, df2, df3)
  }
  
  if(number_of_models == 4){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnAICc(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnAICc(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnAICc(model = x)))
    bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnAICc(model = x)))
    
    df1 <- data.frame(model = "ER", AICc_score = ER_likelihoods)
    df2 <- data.frame(model = "SYM", AICc_score = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", AICc_score = ARD_likelihoods)
    df4 <- data.frame(model = "bridge_only", AICc_score = bridge_only_likelihoods)
    df_full <- rbind(df1, df2, df3, df4)
  }
  
  return(df_full)
}
  

# # Function 6: Transition rates from 1k model results -------------------------

# plot1kTransitionRates(model_results <- readRDS(here("finalized_1k_models/artiodactyla_three_state_max_crep_traits_ER_SYM_ARD_bridge_only_models.rds")), states_in_model = 3, number_of_models = 4){
#   if("ER" %in% number_of_models){
#     rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
#     ER_rates_df <- as.data.frame(rates)
#     ER_rates_df$model <- "ER"
#     ER_rates_df <- ER_rates_df[!(is.na(ER_rates_df$rates)),]
#     
#     if(number_of_states = 3){
#       ER_rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
#       ER_rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
#     }
#     
#     if(number_of_states = 4){
#       #TBD
#     }
#     
#     if(number_of_states = 6){
#       ER_rates_df$solution <- rep(c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 1000)
#       ER_rates_colours <- rep(c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab"),1000)
#     }
#   }
#   
#   if("SYM" %in% number_of_models){
#     rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
#     SYM_rates_df <- as.data.frame(rates)
#     SYM_rates_df$model <- "SYM"
#     SYM_rates_df <- SYM_rates_df[!(is.na(SYM_rates_df$rates)),]
#     
#     if(number_of_states = 3){
#       SYM_rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
#       SYM_rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
#     }
#     
#     if(number_of_states = 4){
#       #TBD
#     }
#     
#     if(number_of_states = 6){
#       SYM_rates_df$solution <- rep(c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 1000)
#       SYM_rates_colours <- rep(c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab"),1000)
#     }
#   }
#   
#   if("ARD" %in% number_of_models){
#     rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
#     ARD_rates_df <- as.data.frame(rates)
#     ARD_rates_df$model <- "ARD"
#     ARD_rates_df <- ARD_rates_df[!(is.na(ARD_rates_df$rates)),]
#     
#     if(number_of_states = 3){
#       ARD_rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc"), 1000)
#       ARD_rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
#     }
#     
#     if(number_of_states = 4){
#     #TBD
#     }
#     
#     if(number_of_states = 6){
#       rates_df$solution <- rep(c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep"), 1000)
#       rates_colours <- rep(c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab"),1000)
#     }
#   }
#     
#   if("bridge_only" %in% number_of_models){
#     rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
#     bridge_only_rates_df <- as.data.frame(rates)
#     bridge_only_rates_df$model <- "Bridge_only"
#     bridge_only_rates_df <- bridge_only_rates_df[!(is.na(bridge_only_rates_df$rates)),]
#     
#     if(number_of_states = 3){
#       bridge_only_rates_df$solution <- rep(c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc"), 1000)
#       bridge_only_rates_colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")
#       }
#     
#     if(number_of_states = 4){
#       #TBD
#     }
#     
#     if(number_of_states = 6){
#       #TBD
#       }
#     }
#  
#   if(number_of_models = 3){
#     rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df)
#     }
#   
#   if(number_of_models = 4){
#     rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df, bridge_only_rates_df)
#     }
#   
#   return(rates_full)
# }

# # Function 7: Transition rate scatterplots from 1k model results -------------------------

