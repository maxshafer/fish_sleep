source("scripts/fish_sleep_functions.R")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
# Continuous traits -------------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#chose the variable we will look at
# continuous_trait <- "Orbit_ratio"
# trait.data <- trait.data[!is.na(trait.data$Orbit_ratio),]

# continuous_trait <- "Dive_depth_m"
# trait.data <- trait.data[!is.na(trait.data$Dive_depth_m),]

# continuous_trait <- "Body_mass_kg"
# trait.data <- trait.data[!is.na(trait.data$Body_mass_kg),]

# continuous_trait <- "Body_length_m"
# trait.data <- trait.data[!is.na(trait.data$Body_length_m),]

#use below instead for artiodactyla orbit size
# trait.data <- read.csv(here("artio_orbit_ratio.csv"))
# continuous_trait <- "Orbit_ratio"
# trait.data <- trait.data[!is.na(trait.data$Orbit_ratio),]
#trait.data <- filter(trait.data, Order == "Artiodactyla")

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]

#use below for just one Parvorder
#trait.data <- trait.data %>% filter(Parvorder == "Odontoceti")
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
#k = 1.97676 for orbit ratio, p-value = 0.001

lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
lambda <- round(lambda$lambda, digits = 2)
#lambda = 1.01589 for orbit ratio, p-value = 1.297 e-12

#custom.colours <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", continuous_trait)]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = trait.data[,c(continuous_trait)]), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = paste(continuous_trait, "\n", "λ = ", lambda))
diel.plot <- diel.plot 
diel.plot

comparison_list <- list(c("cathemeral", "crepuscular"), c("crepuscular", "diurnal"), c("diurnal", "nocturnal"), c("cathemeral", "diurnal"), c("cathemeral", "nocturnal"),  c("crepuscular", "nocturnal"))

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
boxplot_KW <- ggplot(trait.data, aes(x = Diel_Pattern, y = trait.data[,c(continuous_trait)])) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = continuous_trait) + scale_fill_manual(values=unique(trait.data$fam_colours)) +
  stat_compare_means(comparisons = comparison_list) + stat_compare_means(label.y = 0) #+ facet_wrap(~Parvorder)
boxplot_KW
#phylogenetically corrected ANOVA
phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant for cetaceans, p value of F statistic = 0.441
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.002

boxplot_anova <- ggplot(trait.data, aes(x = Diel_Pattern, y = trait.data[,c(continuous_trait)])) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = continuous_trait) + scale_fill_manual(values=unique(trait.data$fam_colours))  + theme_minimal() +
  theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent')) + 
  stat_compare_means(label = "p.format", label.y = 3.5, method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 4.15, method = "anova") + 
  annotate("text", x = 1.15, y = 4, label = paste("phylANOVA, p =", phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_anova

#save out the plots
pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_vs_diel_cetaceans.pdf"))
diel.plot
dev.off()

# pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_boxplots_KW_cetaceans.pdf"))
# boxplot_KW
# dev.off()

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", continuous_trait, "_boxplots_anova_cetaceans.pdf"))
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

#chose the variable we will look at
# discrete_trait <- "Habitat"
# trait.data <- trait.data[!is.na(trait.data$Habitat),]

# discrete_trait <- "Feeding_method"
# trait.data <- trait.data[!is.na(trait.data$Feeding_method),]

# discrete_trait <- "Prey_capture"
# trait.data <- trait.data[!is.na(trait.data$Prey_capture),]

# discrete_trait <- "Diet"
# trait.data <- trait.data[!is.na(trait.data$Diet),]

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]

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



# Multivariate analysis ---------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

trait.data <- trait.data[!is.na(trait.data$Habitat),]

#we do not have dive depth data for any riverine species
#given them a small dive depth (50m)
#trait.data[trait.data$Habitat == "riverine", c("Dive_depth")] <- 50

trait.data <- trait.data[!is.na(trait.data$Dive_depth),]


ggplot(trait.data, aes(x = Habitat, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Diel_Pattern), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + scale_fill_manual(values=custom.colours) +
  stat_compare_means(label = "p.format", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 2000, method = "anova") + facet_wrap(~ Diel_Pattern)


ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Habitat), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + scale_fill_manual(values=custom.colours) +
  stat_compare_means(label = "p.format", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 2000, method = "anova") + facet_wrap(~ Habitat)

#two way ANOVA of habitat and diel pattern on dive depth

aggregate(Dive_depth ~ Habitat + Diel_Pattern, data = trait.data, FUN = mean)

ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth, fill = Habitat)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=custom.colours)
 
model <- aov(Dive_depth ~ Habitat + Diel_Pattern, data = trait.data)
summary(model)