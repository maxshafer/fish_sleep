#run the bridge only model for all mammal suborders in bennie et al data

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data <- trait.data[!duplicated(trait.data$tips), c("max_crep", "tips", "Order")]
trait.data <- trait.data[trait.data$max_crep %in% c("cathemeral", "nocturnal", "diurnal", "crepuscular"), ]
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #4228 mammals in final tree
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips) #we will subset again by order later
table(trait.data$Order) 

table(trait.data$Order, trait.data$max_crep) 

ggplot(trait.data, aes(x = Order, fill = max_crep)) + geom_bar()

#optional: filter for orders with 4 states
trait.data <- trait.data %>% group_by(Order) %>% filter("diurnal" %in% max_crep & "nocturnal" %in% max_crep & "cathemeral" %in% max_crep & "crepuscular" %in% max_crep)

#filter for orders that have over 100 species
#trait.data <- trait.data %>% group_by(Order) %>% filter(n() > 1000) %>% ungroup()
#leaves 7 orders = artio, carnivora, dipto, eulipo, primates, rodentia
#remove Bennie et al artiodactyla data

#trait.data <- filter(trait.data, Order != "Artiodactyla")

trait.data <- as.data.frame(trait.data)

#sub in our artiodactyla data
artio <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio <- artio[!is.na(artio$max_crep), c("tips", "max_crep", "Order")]
artio$Order <- "Amelia_artiodactyla"

trait.data <- rbind(trait.data, artio)

calcMaxCladeCred <- function(trait.data = trait.data, Order_name = "Primates", phy = mam.tree){
  trait.data.1 <- trait.data[trait.data$Order == Order_name,]
  trait.data.1 <- trait.data.1[trait.data.1$tips %in% phy$tip.label,]
  mam.tree.1 <- keep.tip(phy, tip = trait.data.1$tips)
  trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep), c("tips", "max_crep")]
  model_result <- corHMM(phy = mam.tree.1, data = trait.data.1, node.states = "marginal", rate.cat = 1, rate.mat = NULL, model = "ARD")
  return(model_result)
}

#mam_model_result_list <- lapply(unique(trait.data$Order), function(x) calcMaxCladeCred(trait.data = trait.data, Order_name = x, phy = mam.tree))
#names(mam_model_result_list) <- unique(trait.data$Order)

#saveRDS(mam_model_result_list, here("mam_model_result_list.rds")) #add rodents
mam_model_result_list <- readRDS(here("mam_model_result_list.rds"))
#save out the result
#saveRDS(mam_model_result_list, here("mam_model_result_list.rds"))
mam_model_result_list2 <- readRDS(here("mam_model_result_list2.rds"))

mam_model_result_list <- append(mam_model_result_list, mam_model_result_list2)

setwd(here())
source("scripts/Amelia_plotting_functions.R")

likelihood_metrics <- max_clade_metrics(mam_model_result_list)
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)

#can't compare likelihood for differnt groups
ggplot(likelihood_metrics, aes(x = model, y = AIC_scores)) + geom_bar(stat = "identity") + geom_text(aes(label = round(AIC_scores, digits = 3)), vjust = -0.2)

#but we can compare rates
rates_list <- lapply(mam_model_result_list, function(x) max_clade_rates(data = x))

rates_df <- bind_rows(rates_list, .id = "Order")

ggplot(rates_df, aes(x = solution, y = rate, fill = Order)) + geom_bar(stat = "identity", position = "dodge") + facet_wrap(~start_state, scales = "free") #+ geom_text(aes(label = round(rate, digits = 3)), vjust = -0.2)

rates_df1 <- rates_df %>% group_by(Order) %>% summarize(mean_rates = mean(rate))

ggplot(rates_df1, aes(x = Order, y = mean_rates, fill = Order)) + geom_bar(stat = "identity", position = "dodge") + geom_text(aes(label = round(mean_rates, digits = 3)), vjust = -0.2)

#how do the rates compare to the percentage of cathemeral or crepuscular species

trait.data.1 <- trait.data %>% group_by(Order, max_crep) %>% count(max_crep)
trait.data.2 <- trait.data %>% group_by(Order) %>% count(Order)
colnames(trait.data.2) <- c("Order", "totals")

trait.data.1 <- merge(trait.data.1, trait.data.2, all = TRUE)
trait.data.1$proportion <- trait.data.1$n/trait.data.1$totals

trait.data.1 <- merge(trait.data.1, rates_df1, all = TRUE)

trait.data.1 %>% filter(max_crep == "crepuscular") %>% ggplot(., aes(x = mean_rates, y = proportion)) + geom_point(size = 3, aes(colour = Order)) + geom_smooth(method = "lm") + stat_poly_eq()
trait.data.1 %>% filter(max_crep == "cathemeral") %>% ggplot(., aes(x = mean_rates, y = proportion)) + geom_point(size = 3, aes(colour = Order)) + geom_smooth(method = "lm")+ stat_poly_eq()
trait.data.1 %>% filter(max_crep == "nocturnal") %>% ggplot(., aes(x = mean_rates, y = proportion)) + geom_point(size = 3, aes(colour = Order)) + geom_smooth(method = "lm")+ stat_poly_eq()
trait.data.1 %>% filter(max_crep == "diurnal") %>% ggplot(., aes(x = mean_rates, y = proportion)) + geom_point(size = 3, aes(colour = Order)) + geom_smooth(method = "lm")+ stat_poly_eq()

ggplot(trait.data.1, aes(x = mean_rates, y = proportion)) + 
  geom_point(size = 3, aes(colour = Order)) + geom_smooth(method = "lm") +
  stat_poly_eq() + facet_wrap(~max_crep)

trait.data.2 <- merge(trait.data.1, rates_df, all = TRUE)

#filtering for orders with more than 20 species and replacing artio with my artio data
trait.data.2 %>% filter(Order %in% c("Amelia_artiodactyla", "Carnivora", "Didelphimorpha", "Lagomorpha", "Primates", "Rodentia")) %>% 
  ggplot(., aes(x = mean_rates, y = proportion)) + 
  geom_jitter(size = 3, aes(colour = Order)) + geom_smooth(method = "lm") +
  stat_poly_eq() + facet_wrap(~max_crep)

trait.data.2 %>% filter(Order %in% c("Amelia_artiodactyla", "Carnivora", "Didelphimorpha", "Lagomorpha", "Primates", "Rodentia")) %>% 
  ggplot(., aes(x = rate, y = proportion)) + 
  geom_jitter(size = 2, aes(colour = Order)) + geom_smooth(method = "lm") +
  stat_poly_eq() + facet_wrap(~max_crep)
