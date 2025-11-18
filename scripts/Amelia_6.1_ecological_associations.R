source("scripts/fish_sleep_functions.R")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
# Continuous traits -------------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#chose the variable we will look at
# continuous_trait <- "Orbit_ratio"
# trait.data <- trait.data[!is.na(trait.data$Orbit_ratio),]

continuous_trait <- "Dive_depth_m"
trait.data <- trait.data[!is.na(trait.data$Dive_depth_m),]

# continuous_trait <- "Body_mass_kg"
# trait.data <- trait.data[!is.na(trait.data$Body_mass_kg),]

#use below instead for artiodactyla orbit size
# trait.data <- read.csv(here("artio_orbit_ratio.csv"))
# continuous_trait <- "Orbit_ratio"
# trait.data <- trait.data[!is.na(trait.data$Orbit_ratio),]
# trait.data <- filter(trait.data, Order == "Artiodactyla")

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]

#use below for just one Parvorder
trait.data <- trait.data %>% filter(Parvorder == "Odontoceti")
# trait.data <- trait.data %>% filter(Parvorder == "Mysticeti")

#filter for species in the final tree
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#calculate the phyogenetic signal
#requires a vector of trait data in the same species order as in the tree
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
trait.vector <- test[, c(continuous_trait)]
names(trait.vector) <- test$tips

#y trait is the continuous (response) variable 
trait.y <- test[, c(continuous_trait)]
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
lambda <- round(lambda$lambda, digits = 2)

#custom.colours <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", continuous_trait)]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = trait.data[,c(continuous_trait)]), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = paste(continuous_trait, "\n", "λ = ", lambda))
diel.plot <- diel.plot + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

comparison_list <- list(c("cathemeral", "crepuscular"), c("crepuscular", "diurnal"), c("diurnal", "nocturnal"), c("cathemeral", "diurnal"), c("cathemeral", "nocturnal"),  c("crepuscular", "nocturnal"))

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
# boxplot_KW <- ggplot(trait.data, aes(x = Diel_Pattern, y = trait.data[,c(continuous_trait)])) + geom_boxplot(outlier.shape = NA) +  
#   geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
#   labs(x = "Temporal activity pattern", y = continuous_trait) + scale_fill_manual(values=unique(trait.data$fam_colours)) +
#   stat_compare_means(comparisons = comparison_list) + stat_compare_means(label.y = 0) #+ facet_wrap(~Parvorder)
# boxplot_KW

#phylogenetically corrected ANOVA
phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant for cetaceans, p value of F statistic = 0.441
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.002

boxplot_anova <- ggplot(trait.data, aes(x = Diel_Pattern, y = trait.data[,c(continuous_trait)])) +
  geom_boxplot(outlier.shape = NA, aes(fill = Diel_Pattern)) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = continuous_trait) + scale_fill_manual(values=unique(trait.data$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 3200, method = "anova") +annotate("text", x = 1.15, y = 3500, label = paste("phylANOVA, p =", phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_anova

#cet orbit ratio: label.y = 
#save out the plots
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_vs_diel_cetaceans.pdf"))
diel.plot
dev.off()

# pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_boxplots_KW_cetaceans.pdf"))
# boxplot_KW
# dev.off()

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_anova
dev.off()

# pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_boxplots_anova_artio.pdf"))
# boxplot_anova
# dev.off()

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
orbit_ANOVA <- aov(Orbit_ratio ~ Diel_Pattern, data = trait.data)
summary(orbit_ANOVA) #statistically significant for cetaceans, p value of F statistic = 0.014
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.00213

#use the tukey HSD test 
#library(multcomp)
post_test <- glht(orbit_ANOVA, linfct = mcp(Diel_Pattern = "Tukey"))

summary(post_test)


#pairwise t test
tapply(X = trait.data$Orbit_ratio, INDEX = list(trait.data$Diel_Pattern), FUN = mean)
pairwise.t.test(trait.data$Orbit_ratio, trait.data$Diel_Pattern,  p.adj = "none")


# Discrete traits ---------------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]

#chose the variable we will look at
# discrete_trait <- "Habitat"
# trait.data <- trait.data[!is.na(trait.data$Habitat),]

# discrete_trait <- "Feeding_method"
# trait.data <- trait.data[!is.na(trait.data$Feeding_method),]
# 
# discrete_trait <- "Prey_capture" #less species than feeding method (60 vs 76) but come from one source (Churchill)
# trait.data <- trait.data[!is.na(trait.data$Prey_capture),]

discrete_trait <- "Diet"
trait.data <- trait.data[!is.na(trait.data$Diet),]

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#filter for species in the final tree
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
#custom.colours.2 <- c( "grey90", "grey40", "grey66", "black", "red")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", discrete_trait)]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = trait.data[, c(discrete_trait)]), inherit.aes = FALSE, colour = "transparent") + scale_fill_grey(name = discrete_trait)
diel.plot <- diel.plot
diel.plot

test_df <- trait.data[, c("Diel_Pattern", discrete_trait)]
test <- fisher.test(table(test_df))

mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 

ggplot(trait.data, aes(x = Diel_Pattern, fill = trait.data[, c(discrete_trait)])) + geom_histogram(stat = "count", position = "dodge") + scale_fill_grey(name = discrete_trait) +theme_minimal()

#test to see if aquatic lifestyle is associated with cathemerality
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern),]
artio_full$enviro <- "aquatic"
for(i in 1:nrow(artio_full)){
  if(artio_full[i, "Parvorder"] == "non-cetacean"){
    artio_full[i, "enviro"] <- "terrestrial"
  }
}

discrete_trait <- "enviro"

test_df <- artio_full[, c("max_crep", discrete_trait)]

test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 9))) 
ggplot(artio_full, aes(x = max_crep, fill = artio_full[, c(discrete_trait)])) + geom_histogram(stat = "count", position = "dodge") + scale_fill_grey(name = discrete_trait) +theme_minimal()


# Multivariate analysis ---------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

trait.data <- trait.data[!is.na(trait.data$Habitat),]

#we do not have dive depth data for any riverine species
#given them a small dive depth (50m)
#trait.data[trait.data$Habitat == "riverine", c("Dive_depth_m")] <- 50

trait.data <- trait.data[!is.na(trait.data$Dive_depth_m),]


ggplot(trait.data, aes(x = Habitat, y = Dive_depth_m)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Diel_Pattern), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + scale_fill_manual(values=custom.colours) +
  stat_compare_means(label = "p.format", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 2000, method = "anova") + facet_wrap(~ Diel_Pattern)


ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth_m)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Habitat), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + scale_fill_manual(values=custom.colours) +
  stat_compare_means(label = "p.format", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 2000, method = "anova") + facet_wrap(~ Habitat)

#two way ANOVA of habitat and diel pattern on dive depth

aggregate(Dive_depth_m ~ Habitat + Diel_Pattern, data = trait.data, FUN = mean)

ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth_m, fill = Habitat)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=custom.colours)
 
model <- aov(Dive_depth_m ~ Habitat + Diel_Pattern, data = trait.data)
summary(model)


# # Function to create trait vector for phylANOVA -------------------------

calculatePhylANOVA <- function(trait.data = trait.data, continuous_trait = "Orbit_ratio"){
  trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
  trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
  sps_order <- as.data.frame(trpy_n$tip.label)
  colnames(sps_order) <- "tips"
  sps_order$id <- 1:nrow(sps_order)
  test <- merge(trait.data, sps_order, by = "tips")
  trait.vector <- test[, c(continuous_trait)]
  names(trait.vector) <- test$tips
  #y trait is the continuous (response) variable 
  trait.y <- test[, c(continuous_trait)]
  names(trait.y) <- test$tips
  #x trait is the categorical variable 
  trait.x <- test$Diel_Pattern
  names(trait.x) <- test$tips
  phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
  return(phylANOVA)
}


# Save out all continuous trait plots -------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))
#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]

#chose the variable we will look at
trait.data.1 <- trait.data[!is.na(trait.data$Orbit_ratio),]

ggplot(trait.data.1, aes(x = log(Bizygomatic_width), y = log(Average_orbit_length), colour = max_crep)) +
  geom_point() + geom_smooth(method = "lm", na.rm = T, se = F, formula = y~x, aes(colour =max_crep))

#perform the one-way anova
orbit_model <- aov(Orbit_ratio ~ Diel_Pattern, data = trait.data.1)
summary(orbit_model)

#perform the post-hoc tukey test
TukeyHSD(orbit_model, conf.level = .95)
#plot(TukeyHSD(orbit_model, conf.level=.95), las = 2)

#perform the phylogenetically corrected one-way anova
orbit_phylANOVA <- calculatePhylANOVA(trait.data.1, "Orbit_ratio")

#plot out the group means
boxplot_orbit <- ggplot(trait.data.1, aes(x = Diel_Pattern, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = Diel_Pattern), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit size") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 38, method = "anova") +annotate("text", x = 1.15, y = 40, label = paste("phylANOVA, p =", orbit_phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_orbit

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Orbit_ratio", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_orbit
dev.off()

#chose the variable we will look at
trait.data.1 <- trait.data[!is.na(trait.data$Body_mass_kg),]

#perform the one-way anova
mass_model <- aov(Body_mass_kg ~ Diel_Pattern, data = trait.data.1)
summary(mass_model)

#perform the post-hoc tukey test
TukeyHSD(mass_model, conf.level = .95)
#plot(TukeyHSD(orbit_model, conf.level=.95), las = 2)

#perform the phylogenetically corrected one-way anova
mass_phylANOVA <- calculatePhylANOVA(trait.data.1, "Body_mass_kg")

#plot out the group means
boxplot_mass <- ggplot(trait.data.1, aes(x = Diel_Pattern, y = log(Body_mass_kg))) +
  geom_boxplot(aes(fill = Diel_Pattern), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Log (Body mass kg)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 12.5, method = "anova") +annotate("text", x = 1.15, y = 12, label = paste("phylANOVA, p =", mass_phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_mass

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Body_mass", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_mass
dev.off()

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$Dive_depth_m),]

#perform the one-way anova
dive_model <- aov(Dive_depth_m ~ Diel_Pattern, data = trait.data.1)
summary(dive_model)

#perform the post-hoc tukey test
TukeyHSD(dive_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
dive_phylANOVA <- calculatePhylANOVA(trait.data.1, "Dive_depth_m")

boxplot_dive <- ggplot(trait.data.1, aes(x = Diel_Pattern, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = Diel_Pattern), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Maximum dive depth (m)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
   stat_compare_means(label.y = 3400, method = "anova") +annotate("text", x = 1.15, y = 3300, label = paste("phylANOVA, p =", dive_phylANOVA$Pf)) + facet_wrap(~Parvorder)
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Dive_depth", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_dive
dev.off()

trait.data.od <- filter(trait.data, Parvorder == "Odontoceti")
trait.data.od <- trait.data.od[!is.na(trait.data.od$Dive_depth_m),]

trait.data.my <- filter(trait.data, Parvorder == "Mysticeti")
trait.data.my <- trait.data.my[!is.na(trait.data.my$Dive_depth_m),]

#perform the phylogenetically corrected one-way anova
odonto_phylANOVA <- calculatePhylANOVA(trait.data.od, "Dive_depth_m")
mystic_phylANOVA <- calculatePhylANOVA(trait.data.my, "Dive_depth_m")

dat_text <- data.frame(label = c(paste("phylANOVA =", mystic_phylANOVA$Pf, sep = " "), paste("phylANOVA =", odonto_phylANOVA$Pf, sep = " ")), Parvorder = c("Mysticeti", "Odontoceti"))

boxplot_dive <- ggplot(trait.data.1, aes(x = Diel_Pattern, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = Diel_Pattern), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Maximum dive depth (m)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.x = 0.8, label.y = 3400, method = "anova") + 
  geom_text(data= dat_text,mapping = aes(x = 1, y = 3300, label = label)) +
  facet_wrap(~Parvorder, scales = "free_x")
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Dive_depth", "_boxplots_anova_parvorder.pdf"), width = 12, height = 7)
boxplot_dive
dev.off()

trait.data.art <- read.csv(here("artio_orbit_ratio.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$Diel_Pattern),]
trait.data.art$Diel_Pattern <- str_replace(trait.data.art$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data.art$Diel_Pattern <- str_replace(trait.data.art$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data.art$Diel_Pattern <- str_replace(trait.data.art$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")
#use below to remove 6 perissodactyla sps
trait.data.art <- filter(trait.data.art, Order == "Artiodactyla")
trait.data.art <- trait.data.art[!is.na(trait.data.art$Orbit_ratio), c("tips", "Orbit_ratio", "Diel_Pattern", "Family", "fam_colours")]

#perform the one-way anova
art_model <- aov(Orbit_ratio ~ Diel_Pattern, data = trait.data.art)
summary(art_model)

#perform the post-hoc tukey test
TukeyHSD(art_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
art_phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#cannot calculate the pairwise corrected p-values likley because the cathemeral sample size is too small
#remove cathemeral if artio only and rerun
#trait.data.art <- filter(trait.data.art, Diel_Pattern != "cathemeral")
# art_phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#add p values manually
library(rstatix)
#for artio
stat.test <- data.frame(group1 = c("crepuscular", "crepuscular", "diurnal"),
                        group2 = c("diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(art_phylANOVA$Pt[2], art_phylANOVA$Pt[3], art_phylANOVA$Pt[6]),
                        y.position = c(0.95, 1.0, 0.994))

stat.test <- stat.test %>% add_x_position(x = "Diel_Pattern")

boxplot_art <- ggplot(trait.data.art, aes(x = Diel_Pattern, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = Diel_Pattern), alpha=0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Corneal diameter: axial length") + scale_fill_manual(values=unique(trait.data.art$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 1.01, method = "anova") + annotate("text", x = 1.15, y = 1.007, label = paste("phylANOVA, p =", art_phylANOVA$Pf)) +
  stat_pvalue_manual(stat.test, label = "p.adj")
boxplot_art

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Orbit_ratio", "_boxplots_anova_artio.pdf"), width = 8, height = 7)
boxplot_art
dev.off()

# continuous_trait <- "Body_mass_kg"
# trait.data <- trait.data[!is.na(trait.data$Body_mass_kg),]

# continuous_trait <- "Body_length_m"
# trait.data <- trait.data[!is.na(trait.data$Body_length_m),]

#use below instead for artiodactyla orbit size
# trait.data <- read.csv(here("artio_orbit_ratio.csv"))
# continuous_trait <- "Orbit_ratio"
# trait.data <- trait.data[!is.na(trait.data$Orbit_ratio),]
#trait.data <- filter(trait.data, Order == "Artiodactyla")


boxplot_anova <- ggplot(trait.data, aes(x = Diel_Pattern, y = trait.data[,c(continuous_trait)])) +
  geom_boxplot(outlier.shape = NA, aes(fill = Diel_Pattern)) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = continuous_trait) + scale_fill_manual(values=unique(trait.data$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label = "p.format", label.y = 37, method = "t.test",ref.group = ".all.") + 
  stat_compare_means(label.y = 39, method = "anova") +annotate("text", x = 1.15, y = 41, label = paste("phylANOVA, p =", phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_anova
