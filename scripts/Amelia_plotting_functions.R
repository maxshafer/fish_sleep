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

max_clade_rates1 <- function(data = mam_model_result_list$Carnivora){
  if(length(unique(data$data$max_crep)) == 4){
    rates <- as.data.frame(data$solution)
    colnames(rates) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
    row.names(rates) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
    rates$start_state <- row.names(rates)
    rates <- pivot_longer(rates, cols = !start_state, names_to = "end_state", values_to = "rate")
    rates <- as.data.frame(rates)
    rates$solution <- paste(rates$start_state, "to", rates$end_state, sep = "-")
    rates <- rates[!is.na(rates$rate), ]
  }
  if(length(unique(data$data$max_crep)) < 4){
    stop("Less than 4 states in the model")
  }
  
  return(rates)
}

max_clade_rates2 <- function(data = mam_model_result_list$Carnivora){
  if(length(unique(data$data$Diel_Pattern)) == 4){
    rates <- as.data.frame(data$solution)
    colnames(rates) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
    row.names(rates) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
    rates$start_state <- row.names(rates)
    rates <- pivot_longer(rates, cols = !start_state, names_to = "end_state", values_to = "rate")
    rates <- as.data.frame(rates)
    rates$solution <- paste(rates$start_state, "to", rates$end_state, sep = "-")
    rates <- rates[!is.na(rates$rate), ]
  }
  if(length(unique(data$data$Diel_Pattern)) < 4){
    stop("Less than 4 states in the model")
  }
  
  return(rates)
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
  
  if(number_of_models == 5){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnLikelihoods(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnLikelihoods(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnLikelihoods(model = x)))
    bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnLikelihoods(model = x)))
    CONSYM_likelihoods  <- unlist(lapply(model_results$CONSYM_model, function(x) returnLikelihoods(model = x)))
    
    df1 <- data.frame(model = "ER", likelihoods = ER_likelihoods)
    df2 <- data.frame(model = "SYM", likelihoods = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", likelihoods = ARD_likelihoods)
    df4 <- data.frame(model = "bridge_only", likelihoods = bridge_only_likelihoods)
    df5 <- data.frame(model = "CONSYM", likelihoods = CONSYM_likelihoods)
    df_full <- rbind(df1, df2, df3, df4, df5)
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
  
  if(number_of_models == 5){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnAIC(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnAIC(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnAIC(model = x)))
    bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnAIC(model = x)))
    CONSYM_likelihoods  <- unlist(lapply(model_results$CONSYM_model, function(x) returnAIC(model = x)))
    
    df1 <- data.frame(model = "ER", AIC_score = ER_likelihoods)
    df2 <- data.frame(model = "SYM", AIC_score = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", AIC_score = ARD_likelihoods)
    df4 <- data.frame(model = "bridge_only", AIC_score = bridge_only_likelihoods)
    df5 <- data.frame(model = "CONSYM", AIC_score = CONSYM_likelihoods)
    df_full <- rbind(df1, df2, df3, df4, df5)
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
  
  if(number_of_models == 5){
    ER_likelihoods <- unlist(lapply(model_results$ER_model, function(x) returnAICc(model = x)))
    SYM_likelihoods  <- unlist(lapply(model_results$SYM_model, function(x) returnAICc(model = x)))
    ARD_likelihoods  <- unlist(lapply(model_results$ARD_model, function(x) returnAICc(model = x)))
    bridge_only_likelihoods  <- unlist(lapply(model_results$bridge_only_model, function(x) returnAICc(model = x)))
    CONSYM_likelihoods  <- unlist(lapply(model_results$CONSYM_model, function(x) returnAICc(model = x)))
    
    df1 <- data.frame(model = "ER", AICc_score = ER_likelihoods)
    df2 <- data.frame(model = "SYM", AICc_score = SYM_likelihoods)
    df3 <- data.frame(model = "ARD", AICc_score = ARD_likelihoods)
    df4 <- data.frame(model = "bridge_only", AICc_score = bridge_only_likelihoods)
    df5 <- data.frame(model = "CONSYM", AICc_score = CONSYM_likelihoods)
    df_full <- rbind(df1, df2, df3, df4, df5)
  }
  
  return(df_full)
}
  


# # Function 5: ANOVA assumptions -----------------------------------------
anovaAssumptions <- function(results_df = df_full, metric = results_df$likelihoods) {
  # takes the dataframe of all model results, determines if it means the assumptions of ANOVA
  # verify equality of variances
  lav_test <- leveneTest(metric ~ factor(model),
                         data = results_df,
  )
  #if the p value is less than 0.05 the residuals fit the assumption of equal variance
  if(lav_test$`Pr(>F)`[1] < 0.05){
    residual_variance <- "Equal-variance"  
  } else{
    #error variance not equal
    residual_variance <- "Not-equal-variance"
  }
  
  #run an ANOVA to check if the residuals are normally distributed
  res_aov <- aov(metric ~ model,
                 data = results_df
  )
  #check for normal distribution with the shapiro test
  norm_test <- shapiro.test(res_aov$residuals)
  if(norm_test$p.value < 0.05){
    normality <- "normal"
  } else{
    normality <- "not-normal"
  }
  
  model=lm(metric ~ results_df$model )
  ANOVA=aov(model)
  
  # Tukey test to compare means between test groups
  TUKEY <- TukeyHSD(x=ANOVA, 'results_df$model', conf.level=0.95)
  
  results_list <- c(residual_variance, normality, TUKEY)
  if("Not-equal-variance" %in% results_list) warning("Variance of residuals not equal")
  if("not-normal" %in% results_list) warning("Data is not normally distributed")
  
  return(results_list)
}

# # Function 6: Transition rates from 1k model results -------------------------
plot1kTransitionRates <- function(model_results = readRDS(here("finalized_1k_models/ruminants_six_state_ER_SYM_ARD_models.rds")), states_in_model = 6, number_of_models = 3){

  if(number_of_models == 3){
    models_in_file = c("ER","SYM","ARD")
  }

  if(number_of_models == 4){
    models_in_file = c("ER","SYM","ARD","bridge_only")
  }
  
  if(number_of_models == 5){
    models_in_file = c("ER","SYM","ARD","bridge_only", "CONSYM")
  }

  if("ER" %in% models_in_file){
    rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
    ER_rates_df <- as.data.frame(rates)
    ER_rates_df$model <- "ER"
    ER_rates_df <- ER_rates_df[!(is.na(ER_rates_df$rates)),]

    if(states_in_model == 3){
      ER_rates_df$solution <- c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc")
      ER_rates_df$colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
    }

    if(states_in_model == 4){
      ER_rates_df$solution <- c("Crep -> Cath", "Di -> Cath", "Noc -> Cath", "Cath -> Crep", "Di -> Crep", "Noc -> Crep",  "Cath -> Di", "Crep -> Di", "Noc -> Di", "Cath -> Noc", "Crep -> Noc", "Di -> Noc")
      ER_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#9F7C60", "#BFA895","#D3C3B6",  "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
    }
    
    if(states_in_model == 5){
      ER_rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
      ER_rates_df$colours <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#a63d13", "#bd5c35", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")
    }

    if(states_in_model == 6){
      ER_rates_df$solution <- c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep")
      ER_rates_df$colours <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")
    }
  }

  if("SYM" %in% models_in_file){
    rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
    SYM_rates_df <- as.data.frame(rates)
    SYM_rates_df$model <- "SYM"
    SYM_rates_df <- SYM_rates_df[!(is.na(SYM_rates_df$rates)),]

    if(states_in_model == 3){
      SYM_rates_df$solution <- c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc")
      SYM_rates_df$colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
    }

    if(states_in_model == 4){
      SYM_rates_df$solution <- c("Crep -> Cath", "Di -> Cath", "Noc -> Cath", "Cath -> Crep", "Di -> Crep", "Noc -> Crep",  "Cath -> Di", "Crep -> Di", "Noc -> Di", "Cath -> Noc", "Crep -> Noc", "Di -> Noc")
      SYM_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#9F7C60", "#BFA895","#D3C3B6",  "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
    }

    if(states_in_model == 5){
      SYM_rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
      SYM_rates_df$colours <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#a63d13", "#bd5c35", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")
    }
    
    if(states_in_model == 6){
      SYM_rates_df$solution <- c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep")
      SYM_rates_df$colours <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")
    }
  }

  if("ARD" %in% models_in_file){
    rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
    ARD_rates_df <- as.data.frame(rates)
    ARD_rates_df$model <- "ARD"
    ARD_rates_df <- ARD_rates_df[!(is.na(ARD_rates_df$rates)),]

    if(states_in_model == 3){
      ARD_rates_df$solution <- c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Noc -> Di", "Cath/crep -> Noc", "Di -> Noc")
      ARD_rates_df$colours <- c("deeppink4", "dodgerblue4", "deeppink3", "seagreen4", "dodgerblue2", "seagreen3")
    }

    if(states_in_model == 4){
      ARD_rates_df$solution <- c("Crep -> Cath", "Di -> Cath", "Noc -> Cath", "Cath -> Crep", "Di -> Crep", "Noc -> Crep",  "Cath -> Di", "Crep -> Di", "Noc -> Di", "Cath -> Noc", "Crep -> Noc", "Di -> Noc")
      ARD_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#9F7C60", "#BFA895","#D3C3B6",  "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
    }
    
    if(states_in_model == 5){
      ARD_rates_df$solution <- c("Cath -> Di", "Cath -> Di/crep", "Cath -> Noc", "Cath -> Noc/crep",  "Di -> Cath", "Di -> Di/crep", "Di -> Noc", "Di -> Noc/crep", "Di/crep -> Cath", "Di/crep -> Di", "Di/crep -> Noc", "Di/crep -> Noc/crep", "Noc -> Cath", "Noc -> Di", "Noc -> Di/crep", "Noc -> Noc/crep", "Noc/crep -> Cath", "Noc/crep -> Di", "Noc/crep -> Di/crep", "Noc/crep -> Noc" )
      ARD_rates_df$colours <- c("#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3", "#a63d13", "#bd5c35", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#e9bb65","#facf80","#176d56","#629884","#82ae9d","#a2c4b6", "#5c8816", "#92b264", "#adc887", "#c9deab")
    }

    if(states_in_model == 6 ){
      ARD_rates_df$solution <- c("Cath/crep -> Cath", "Di -> Cath", "Di/crep -> Cath", "Noc -> Cath", "Noc/crep -> Cath", "Cath -> Cath/crep", "Di -> Cath/crep", "Di/crep -> Cath/crep", "Noc -> Cath/crep", "Noc/crep -> Cath/crep", "Cath -> Di", "Cath/crep -> Di", "Di/crep -> Di", "Noc -> Di", "Noc/crep -> Di", "Cath -> Di/crep", "Cath/crep -> Di/crep", "Di -> Di/crep", "Noc -> Di/crep", "Noc/crep -> Di/crep", "Cath -> Noc", "Cath/crep -> Noc", "Di -> Noc", "Di/crep -> Noc", "Noc/crep -> Noc", "Cath -> Noc/crep", "Cath/crep -> Noc/crep", "Di -> Noc/crep", "Di/crep -> Noc/crep", "Noc -> Noc/crep")
      ARD_rates_df$colours <- c("#ac00b6", "#bb46c2", "#ca6ccd", "#d78fd8", "#e3b0e3","#2a2956", "#383e6f", "#47558a", "#556ca4", "#6385bf",  "#a63d13", "#bd5c35", "#d37a57", "#e79979", "#fbb89d", "#b48204", "#c6952e", "#d7a84a", "#e9bb65","#facf80","#176d56", "#40826d","#629884","#82ae9d","#a2c4b6", "#5c8816", "#779d40", "#92b264", "#adc887", "#c9deab")
    }
  }
  

  if("bridge_only" %in% models_in_file){
    rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
    bridge_only_rates_df <- as.data.frame(rates)
    bridge_only_rates_df$model <- "Bridge_only"
    bridge_only_rates_df <- bridge_only_rates_df[!(is.na(bridge_only_rates_df$rates)),]

    if(states_in_model == 3){
      bridge_only_rates_df$solution <- c("Di -> Cath/crep", "Noc -> Cath/crep", "Cath/crep -> Di", "Cath/crep -> Noc")
      bridge_only_rates_df$colours <- c("deeppink4", "dodgerblue4", "deeppink3", "dodgerblue2")
    }

    if(states_in_model == 4){
      bridge_only_rates_df$solution <- c("Crep -> Cath", "Di -> Cath", "Noc -> Cath", "Cath -> Crep", "Di -> Crep", "Noc -> Crep",  "Cath -> Di", "Crep -> Di", "Cath -> Noc", "Crep -> Noc")
      bridge_only_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#9F7C60", "#BFA895","#D3C3B6",  "#FA4A05", "#FC8D62","#3C967E", "#66C2A5")
    }

    # if(number_of_states == 6){
    #   print("six_not_done")
    #   }
  }
  
  if("CONSYM" %in% models_in_file){
    rates <- unlist(lapply(model_results$CONSYM_model, function(x) returnRates(model = x)))
    CONSYM_rates_df <- as.data.frame(rates)
    CONSYM_rates_df$model <- "CONSYM"
    CONSYM_rates_df <- CONSYM_rates_df[!(is.na(CONSYM_rates_df$rates)),]
    
    if(states_in_model == 4){
      CONSYM_rates_df$solution <- c("Crep -> Cath", "Di -> Cath", "Noc -> Cath", "Cath -> Crep", "Di -> Crep", "Noc -> Crep",  "Cath -> Di", "Crep -> Di", "Cath -> Noc", "Crep -> Noc")
      CONSYM_rates_df$colours <-  c( "#A024AE", "#DD8AE7","#EEC4F3", "#9F7C60", "#BFA895","#D3C3B6",  "#FA4A05", "#FC8D62","#3C967E", "#66C2A5")
    }
  }

  if(number_of_models == 3){
    rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df)
  }

  if(number_of_models == 4){
    rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df, bridge_only_rates_df)
    
  }
  
  if(number_of_models == 5){
    rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df, bridge_only_rates_df, CONSYM_rates_df)
  }

  return(rates_full)
}


# #Function 6.5: Transition rates from 1k model only for 4 state  --------
plot1kTransitionRates4state <- function(model_results = readRDS(here("finalized_1k_models/ruminants_six_state_ER_SYM_ARD_models.rds")), number_of_models = 5){
  
  if(number_of_models == 3){
    models_in_file = c("ER","SYM","ARD")
  }
  
  if(number_of_models == 4){
    models_in_file = c("ER","SYM","ARD","bridge_only")
  }
  
  if(number_of_models == 5){
    models_in_file = c("ER","SYM","ARD","bridge_only", "CONSYM")
  }
  
  if("ER" %in% models_in_file){
    rates <- unlist(lapply(model_results$ER_model, function(x) returnRates(model = x)))
    ER_rates_df <- as.data.frame(rates)
    ER_rates_df$model <- "ER"
    ER_rates_df <- ER_rates_df[!(is.na(ER_rates_df$rates)),]
    ER_rates_df$solution <- c("Crepuscular -> Cathemeral", "Diurnal -> Cathemeral", "Nocturnal -> Cathemeral", "Cathemeral -> Crepuscular", "Diurnal -> Crepuscular", "Nocturnal -> Crepuscular",  "Cathemeral -> Diurnal", "Crepuscular -> Diurnal", "Nocturnal -> Diurnal", "Cathemeral -> Nocturnal", "Crepuscular -> Nocturnal", "Diurnal -> Nocturnal")
    ER_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#AD9680", "#D1B49B","#EECBAD",  "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
  }
  
  if("SYM" %in% models_in_file){
    rates <- unlist(lapply(model_results$SYM_model, function(x) returnRates(model = x)))
    SYM_rates_df <- as.data.frame(rates)
    SYM_rates_df$model <- "SYM"
    SYM_rates_df <- SYM_rates_df[!(is.na(SYM_rates_df$rates)),]
    SYM_rates_df$solution <- c("Crepuscular -> Cathemeral", "Diurnal -> Cathemeral", "Nocturnal -> Cathemeral", "Cathemeral -> Crepuscular", "Diurnal -> Crepuscular", "Nocturnal -> Crepuscular",  "Cathemeral -> Diurnal", "Crepuscular -> Diurnal", "Nocturnal -> Diurnal", "Cathemeral -> Nocturnal", "Crepuscular -> Nocturnal", "Diurnal -> Nocturnal")
    SYM_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3","#AD9680", "#D1B49B","#EECBAD", "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
  }
  
  if("ARD" %in% models_in_file){
    rates <- unlist(lapply(model_results$ARD_model, function(x) returnRates(model = x)))
    ARD_rates_df <- as.data.frame(rates)
    ARD_rates_df$model <- "ARD"
    ARD_rates_df <- ARD_rates_df[!(is.na(ARD_rates_df$rates)),]
    ARD_rates_df$solution <- c("Crepuscular -> Cathemeral", "Diurnal -> Cathemeral", "Nocturnal -> Cathemeral", "Cathemeral -> Crepuscular", "Diurnal -> Crepuscular", "Nocturnal -> Crepuscular",  "Cathemeral -> Diurnal", "Crepuscular -> Diurnal", "Nocturnal -> Diurnal", "Cathemeral -> Nocturnal", "Crepuscular -> Nocturnal", "Diurnal -> Nocturnal")
    ARD_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#AD9680", "#D1B49B","#EECBAD",  "#FA4A05", "#FC8D62", "#FECCB9","#3C967E", "#66C2A5","#ABDECE")
  }
  
  if("bridge_only" %in% models_in_file){
    rates <- unlist(lapply(model_results$bridge_only_model, function(x) returnRates(model = x)))
    bridge_only_rates_df <- as.data.frame(rates)
    bridge_only_rates_df$model <- "Bridge_only"
    bridge_only_rates_df <- bridge_only_rates_df[!(is.na(bridge_only_rates_df$rates)),]
    bridge_only_rates_df$solution <- c("Crepuscular -> Cathemeral", "Diurnal -> Cathemeral", "Nocturnal -> Cathemeral", "Cathemeral -> Crepuscular", "Diurnal -> Crepuscular", "Nocturnal -> Crepuscular",  "Cathemeral -> Diurnal", "Crepuscular -> Diurnal", "Cathemeral -> Nocturnal", "Crepuscular -> Nocturnal")
    bridge_only_rates_df$colours <- c( "#A024AE", "#DD8AE7","#EEC4F3", "#AD9680", "#D1B49B","#EECBAD",  "#FA4A05", "#FC8D62","#3C967E", "#66C2A5")
  }
  
  if("CONSYM" %in% models_in_file){
    rates <- unlist(lapply(model_results$CONSYM_model, function(x) returnRates(model = x)))
    CONSYM_rates_df <- as.data.frame(rates)
    CONSYM_rates_df$model <- "CONSYM"
    CONSYM_rates_df <- CONSYM_rates_df[!(is.na(CONSYM_rates_df$rates)),]
    CONSYM_rates_df$solution <- c("Crepuscular -> Cathemeral", "Diurnal -> Cathemeral", "Nocturnal -> Cathemeral", "Cathemeral -> Crepuscular", "Diurnal -> Crepuscular", "Nocturnal -> Crepuscular",  "Cathemeral -> Diurnal", "Crepuscular -> Diurnal", "Cathemeral -> Nocturnal", "Crepuscular -> Nocturnal")
    CONSYM_rates_df$colours <-  c( "#A024AE", "#DD8AE7","#EEC4F3", "#9F7C60", "#BFA895","#D3C3B6",  "#FA4A05", "#FC8D62","#3C967E", "#66C2A5")
    
  if(number_of_models == 3){
    rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df)
  }
  
  if(number_of_models == 4){
    rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df, bridge_only_rates_df)
    
  }
  
  if(number_of_models == 5){
    rates_full <- rbind(ER_rates_df, SYM_rates_df, ARD_rates_df, bridge_only_rates_df, CONSYM_rates_df)
  }
  
  return(rates_full)
  }
  }

# # Function 7: Transition rate scatterplots from 1k model results -------------------------

model_results <- readRDS(here("finalized_1k_models/whippomorpha_four_state_max_crep_ER_SYM_ARD_bridge_only_models.rds"))


# Function 8: Making split violin splots ----------------------------------

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

    

# Max's functions (customized) --------------------------------------------

### Write a function to extract ancestral likelihoods from a model
returnAncestralStates <- function(phylo_model = model, phylo_tree = trpy_n, rate.cat = FALSE, recon = c("joint", "marg")) {
  
  # Create a data frame with trait values from reconstruction
  
  # Joint reconstruction gives a state, rather than the likelihood of each state (as in marginal)
  # Easier to probably make joint look like marginal, to fit the rest of my functions
  if (recon == "joint") {
    tip_states <- data.frame(one_R1 = rep(NA, length(phylo_model$tip.states)), two_R1 = rep(NA, length(phylo_model$tip.states)), one_R2 = rep(NA, length(phylo_model$tip.states)), two_R2 = rep(NA, length(phylo_model$tip.states)))
    tip_states$one_R1 <- ifelse(phylo_model$tip.states == 1, 1, 0)
    tip_states$two_R1 <- ifelse(phylo_model$tip.states == 2, 1, 0)
    tip_states$one_R2 <- ifelse(phylo_model$tip.states == 3, 1, 0)
    tip_states$two_R2 <- ifelse(phylo_model$tip.states == 4, 1, 0)
    states <- data.frame(one_R1 = rep(NA, length(phylo_model$states)), two_R1 = rep(NA, length(phylo_model$states)), one_R2 = rep(NA, length(phylo_model$states)), two_R2 = rep(NA, length(phylo_model$states)))
    states$one_R1 <- ifelse(phylo_model$states == 1, 1, 0)
    states$two_R1 <- ifelse(phylo_model$states == 2, 1, 0)
    states$one_R2 <- ifelse(phylo_model$states == 3, 1, 0)
    states$two_R2 <- ifelse(phylo_model$states == 4, 1, 0)
    lik.anc <- rbind(tip_states, states)
    row.names(lik.anc) <- c(phylo_tree$tip.label, c((Ntip(phylo_tree) + 1):(Ntip(phylo_tree) + Nnode(phylo_tree))))
  }
  
  if (recon == "marg") {
    # state 2 is nocturnal, state 1 is diurnal
    lik.anc <- as.data.frame(rbind(phylo_model$tip.states, phylo_model$states))
    row.names(lik.anc) <- c(row.names(lik.anc)[1:length(phylo_tree$tip.label)], (Ntip(phylo_tree) + 1):(Ntip(phylo_tree) + Nnode(phylo_tree)))
  }
  
  if (phylo_model$rate.cat == 1) {
    print("model has single rate category")
    states <- unique(phylo_model$data[,2])
    names(states) <- as.numeric(unique(phylo_model$data.legend[,2]))
    rate_states <- sort(states)
    states <- sort(states)
  }
  if (phylo_model$rate.cat == 2) {
    print("model has 2 rate categories")
    states <- unique(phylo_model$data[,2])
    names(states) <- as.numeric(unique(phylo_model$data.legend[,2]))
    states <- sort(states)
    rate_states <- c(paste(states, "R1", sep = "_"), paste(states, "R2", sep = "_"))
  }
  if (phylo_model$rate.cat == 3) {
    print("model has 3 rate categories")
    states <- unique(phylo_model$data[,2])
    names(states) <- as.numeric(unique(phylo_model$data.legend[,2]))
    states <- sort(states)
    rate_states <- c(paste(states, "R1", sep = "_"), paste(states, "R2", sep = "_"), paste(states, "R3", sep = "_"))
  }
  
  colnames(lik.anc) <- rate_states
  
  ancestral_states <- list()
  ancestral_states$lik.anc <- lik.anc
  ancestral_states$node <- 1:(length(phylo_tree$tip.label) + phylo_tree$Nnode) # The only thing that isn't correct is the row.names on lik.anc, but I don't use those!
  ancestral_states$states <- states
  ancestral_states$rate_states <- rate_states
  print(paste("returning ancestral states for", phylo_tree$Nnode, "internal nodes corresponding to", length(phylo_tree$tip.label), "tips in tree provided", sep = " "))
  return(ancestral_states)
}

### Function to determine the number of transitions between states
calculateStateTransitions <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, rate.cat = FALSE) {
  
  states <- ancestral_states$states
  rate_states <- ancestral_states$rate_states
  
  # Determine the number of transitions or diurnal/nocturnal taxa by age
  
  # Return the probability of diurnal (or state 1)
  if (length(ancestral_states$states) == length(ancestral_states$rate_states)){
    ancestral_states$recon_states <- ancestral_states$lik.anc[,states[[1]]]
  } else {
    ancestral_states$recon_states <- as.numeric(rowSums(ancestral_states$lik.anc[,rate_states[grep(states[[1]], rate_states)]]))
  }
  
  ## Added functionality to do this for rate categories
  if (rate.cat) {
    print("returning rate categories instead of states")
    ancestral_states$recon_states <- as.numeric(rowSums(ancestral_states$lik.anc[,rate_states[grep("R1", rate_states)]]))
  }
  
  # Determine ages of each node
  node_heights <- nodeHeights(phylo_tree)
  
  ancestor <- round(ancestral_states$recon_states[Ntip(phylo_tree)+1])
  
  ancestral_states$node.age <- node_heights[match(ancestral_states$node, phylo_tree$edge[,2]),2] # This extracts all but one that is missing from phylo_tree$edge[,2], and generates an NA
  ancestral_states$node.age[ancestral_states$node[is.na(ancestral_states$node.age)]] <- 0 # This is the root
  # ancestral_states$node.age[is.na(ancestral_states$node.age)] <- 0
  ancestral_states$node.age <- (ancestral_states$node.age - max(ancestral_states$node.age))*-1
  
  # ancestral_states <- ancestral_states[order(ancestral_states$Time),]
  
  # ID parental nodes and parental states
  ancestral_states$parental.node <- unlist(lapply(ancestral_states$node, function(x) Ancestors(phylo_tree, x, type = "parent")))
  ancestral_states$parent.diel <- unlist(lapply(ancestral_states$parental.node, function(x) ancestral_states$recon_states[match(x, ancestral_states$node)]))
  # parent of the root is NA
  ancestral_states$parent.diel[is.na(ancestral_states$parent.diel)] <- as.numeric(ancestor)
  
  ancestral_states$transition <- ifelse(ifelse(ancestral_states$parent.diel > 0.5, 1, 0) != ifelse(ancestral_states$recon_states > 0.5, 1, 0), 1, 0)
  # For those with transitions, I can just ask what they are, and that's the switch type!
  ancestral_states$trans.ND <- ifelse(ancestral_states$transition == 1, ifelse(ancestral_states$recon_states > 0.5, 1, 0),0)
  ancestral_states$trans.DN <- ifelse(ancestral_states$transition == 1, ifelse(ancestral_states$recon_states < 0.5, 1, 0),0)
  # ancestral_states$transition[is.na(ancestral_states$transition)] <- 0
  
  print(paste("Identified", table(ancestral_states$transition)[2], "transitions between", ancestral_states$states[1], "and", ancestral_states$states[2], sep = " "))
  return(ancestral_states)
}

### Function to calculate transition history on tree for the state
calculateLinTransHist2 <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n) {
  
  # ancestors <- lapply(c(1:length(phylo_tree$tip.label)), function(x) Ancestors(phylo_tree, x, type = "all"))
  # ancestors <- lapply(seq_along(ancestors), function(x) append(c(1:length(phylo_tree$tip.label))[[x]], ancestors[[x]]))
  ## The above only did it for tips, but below does it for all internal nodes as well
  ## 
  ancestors <- lapply(ancestral_states$node, function(x) Ancestors(phylo_tree, x, type = "all"))
  ancestors <- lapply(seq_along(ancestors), function(x) append(ancestral_states$node[[x]], ancestors[[x]]))
  print("calculating ancestral states for all nodes")
  recon_states <- ifelse(ancestral_states$recon_states > 0.5, 1, 0)
  ancestors.diel <- lapply(ancestors, function(x) lapply(x, function(y) recon_states[match(y, ancestral_states$node)]))
  
  print("identifying switch types")
  # This works, can I simplify it so it works on a vector?
  switch.type <- unlist(lapply(ancestors.diel, function(x) {
    df <- data.frame(test = unlist(x))
    df <- df[with(df, c(test[-1]!= test[-nrow(df)], TRUE)),]
    return(paste(df, collapse = ""))
  }))
  
  ancestral_states$switch.type <- as.character(switch.type)
  
  print("Identified the following switch types")
  print(table(switch.type))
  
  return(ancestral_states)
}



### Function to calculate the cummulative sums (of switches and switch types)
returnCumSums <- function(ancestral_states = ancestral_states, phylo_tree = trpy_n, use.height = TRUE) {
  # This is all nodes (for transitions only)
  node_order <- order(ancestral_states$node.age, decreasing = T) # Oldest node first (root)
  
  ancestral_states$transition_cumsum <- cumsum(ancestral_states$transition[node_order])
  ancestral_states$trans.DN_cumsum <- cumsum(ancestral_states$trans.DN[node_order])
  ancestral_states$trans.ND_cumsum <- cumsum(ancestral_states$trans.ND[node_order])
  
  # ## This uses switch types of both internal nodes and tips
  # df <- list()
  # for (i in unique(ancestral_states$switch.type)) {
  #   df[[i]] <- (ifelse(ancestral_states$switch.type == i, 1, 0))
  # }
  # df <- Reduce(cbind, df)
  # colnames(df) <- unique(ancestral_states$switch.type)
  # df <- as.data.frame(df)
  # df$node.age <- ancestral_states$node.age
  
  ## This uses only internal nodes
  df <- list()
  for (i in unique(ancestral_states$switch.type)) {
    df[[i]] <- (ifelse(ancestral_states$switch.type == i, 1, 0))[(Ntip(phylo_tree)+1):length(ancestral_states$switch.type)]
  }
  df <- Reduce(cbind, df)
  colnames(df) <- unique(ancestral_states$switch.type)
  df <- as.data.frame(df)
  df$node.age <- ancestral_states$node.age[(Ntip(phylo_tree)+1):length(ancestral_states$switch.type)]
  
  cumsums <- df
  cumsums <- cumsums[order(cumsums$node.age, decreasing = T),]
  
  for (i in 1:(ncol(cumsums)-1)) {
    cumsums[,i] <- cumsum(cumsums[,i])
  }
  

  ## OK, for this graph, I want to index the first half of node_heights!
  
  # First column from 'edge' is the higher node, second is lower
  # both columns from 'node_heights' should match by the index ('edge')
  # If I find the tip in edge[,2], but pull node_heights[,1], I get the age of the node connected to the tip (the branch point of the tip)

  # node_heights <- nodeHeights(phylo_tree)
  # node.age2 <- node_heights[match(c(1:length(phylo_tree$tip.label)), phylo_tree$edge[,2]),1] # This finds the position of each node # in the 2nd column of 'edge', and returns the first column of 'edge', which is the age of the split
  # 
  # node.age2[is.na(node.age2)] <- 0
  # node.age2 <- (node.age2 - max(node.age2))*-1
  # 
  # cumsums <- data.frame(row.names = phylo_tree$tip.label[order(node.age2, decreasing = T)])
  # cumsums$node.age <- node.age2[order(node.age2, decreasing = T)]
  # # cumsums$node.age <- ancestral_states$node.age[node_order_2]
  # cumsums$Lineage_Cumsum <- cumsum(rep(1, length(phylo_tree$tip.label)))
  # 
  # for (i in unique(ancestral_states$switch.type)) {
  #   cumsums[,i] <- cumsum(ifelse(ancestral_states$switch.type[order(node.age2, decreasing = T)] == i, 1, 0))
  # }
  
  ancestral_states$cumsums <- cumsums
  return(ancestral_states)
}

