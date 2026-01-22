source("scripts/fish_sleep_functions.R")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

# Multivariate analysis ---------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

trait.data <- trait.data[!is.na(trait.data$max_crep),]

trait.data <- trait.data[!is.na(trait.data$Habitat),]

#we do not have dive depth data for any riverine species
#given them a small dive depth (50m)
#trait.data[trait.data$Habitat == "riverine", c("Dive_depth_m")] <- 50

trait.data <- trait.data[!is.na(trait.data$Dive_depth_m),]


ggplot(trait.data, aes(x = Habitat, y = Dive_depth_m)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = max_crep), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + scale_fill_manual(values=custom.colours) +
  stat_compare_means(label = "p.format", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 2000, method = "anova") + facet_wrap(~ max_crep)


ggplot(trait.data, aes(x = max_crep, y = Dive_depth_m)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Habitat), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + scale_fill_manual(values=custom.colours) +
  stat_compare_means(label = "p.format", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 2000, method = "anova") + facet_wrap(~ Habitat)

#two way ANOVA of habitat and diel pattern on dive depth

aggregate(Dive_depth_m ~ Habitat + max_crep, data = trait.data, FUN = mean)

ggplot(trait.data, aes(x = max_crep, y = Dive_depth_m, fill = Habitat)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=custom.colours)
 
model <- aov(Dive_depth_m ~ Habitat + max_crep, data = trait.data)
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
  trait.x <- test$max_crep
  names(trait.x) <- test$tips
  phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
  return(phylANOVA)
}


# Save out all continuous trait plots -------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$max_crep),]

### Cetacean orbit ratio

trait.data.1 <- trait.data[!is.na(trait.data$Orbit_ratio),]

#perform the one-way anova
orbit_model <- aov(Orbit_ratio ~ max_crep, data = trait.data.1)
summary(orbit_model)

#perform the post-hoc tukey test
TukeyHSD(orbit_model, conf.level = .95)
#plot(TukeyHSD(orbit_model, conf.level=.95), las = 2)

#perform the phylogenetically corrected one-way anova
orbit_phylANOVA <- calculatePhylANOVA(trait.data.1, "Orbit_ratio")

#plot out the group means
boxplot_orbit <- ggplot(trait.data.1, aes(x = max_crep, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit size") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 38, method = "anova") +annotate("text", x = 1.15, y = 40, label = paste("phylANOVA, p =", orbit_phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_orbit

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Orbit_ratio", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_orbit
dev.off()

### Cetacean body mass

trait.data.1 <- trait.data[!is.na(trait.data$Body_mass_kg),]

#perform the one-way anova
mass_model <- aov(Body_mass_kg ~ max_crep, data = trait.data.1)
summary(mass_model)

#perform the post-hoc tukey test
TukeyHSD(mass_model, conf.level = .95)
#plot(TukeyHSD(orbit_model, conf.level=.95), las = 2)

#perform the phylogenetically corrected one-way anova
mass_phylANOVA <- calculatePhylANOVA(trait.data.1, "Body_mass_kg")

#plot out the group means
boxplot_mass <- ggplot(trait.data.1, aes(x = max_crep, y = log(Body_mass_kg))) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Log (Body mass kg)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 12.5, method = "anova") +annotate("text", x = 1.15, y = 12, label = paste("phylANOVA, p =", mass_phylANOVA$Pf)) #+ facet_wrap(~Parvorder)
boxplot_mass

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Body_mass", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_mass
dev.off()

#body mass separated by parvorder
trait.data %>% filter(!is.na(Body_mass_kg)) %>% ggplot(., aes(x = max_crep, y = log(Body_mass_kg))) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + facet_wrap(~ Parvorder)

### Cetacean dive depth
trait.data.1 <- trait.data[!is.na(trait.data$Dive_depth_m),]

#perform the one-way anova
dive_model <- aov(Dive_depth_m ~ max_crep, data = trait.data.1)
summary(dive_model)

#perform the post-hoc tukey test
TukeyHSD(dive_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
dive_phylANOVA <- calculatePhylANOVA(trait.data.1, "Dive_depth_m")

boxplot_dive <- ggplot(trait.data.1, aes(x = max_crep, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Maximum dive depth (m)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
   stat_compare_means(label.y = 3400, method = "anova") +annotate("text", x = 1.15, y = 3300, label = paste("phylANOVA, p =", dive_phylANOVA$Pf))
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Dive_depth", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot_dive
dev.off()

### Cetacean dive depth by parvorder
trait.data.od <- filter(trait.data, Parvorder == "Odontoceti")
trait.data.od <- trait.data.od[!is.na(trait.data.od$Dive_depth_m),]

trait.data.my <- filter(trait.data, Parvorder == "Mysticeti")
trait.data.my <- trait.data.my[!is.na(trait.data.my$Dive_depth_m),]

#perform the phylogenetically corrected one-way anova
odonto_phylANOVA <- calculatePhylANOVA(trait.data.od, "Dive_depth_m")
mystic_phylANOVA <- calculatePhylANOVA(trait.data.my, "Dive_depth_m")

dat_text <- data.frame(label = c(paste("phylANOVA =", mystic_phylANOVA$Pf, sep = " "), paste("phylANOVA =", odonto_phylANOVA$Pf, sep = " ")), Parvorder = c("Mysticeti", "Odontoceti"))

boxplot_dive <- ggplot(trait.data.1, aes(x = max_crep, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Maximum dive depth (m)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.x = 0.8, label.y = 3400, method = "anova") + 
  geom_text(data= dat_text,mapping = aes(x = 1, y = 3300, label = label)) +
  facet_wrap(~Parvorder, scales = "free_x")
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Dive_depth_by_Parvorder", "_boxplots_anova_parvorder.pdf"), width = 12, height = 7)
boxplot_dive
dev.off()

### Cetacean mean latitude

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$mean_lat),]
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep),]

#perform the one-way anova
model <- aov(mean_lat ~ max_crep, data = trait.data.1)
summary(model)

#perform the post-hoc tukey test
TukeyHSD(model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "mean_lat")

boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = mean_lat)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Mean latitude") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 6100000, method = "anova") +annotate("text", x = 1.15, y = 6500000, label = paste("phylANOVA, p =", phylANOVA$Pf))
boxplot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "mean_latitude", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot
dev.off()

### Cetacean max latitude

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$max_lat),]
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep),]

#perform the one-way anova
model <- aov(mean_lat ~ max_crep, data = trait.data.1)
summary(model)

#perform the post-hoc tukey test
TukeyHSD(model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "mean_lat")

boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = max_lat)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Mean latitude") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 7100000, method = "anova") +annotate("text", x = 1.15, y = 7500000, label = paste("phylANOVA, p =", phylANOVA$Pf))
boxplot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "max_latitude", "_boxplots_anova_cetaceans.pdf"), width = 8, height = 7)
boxplot
dev.off()

### Artiodactyla orbit size
trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$Orbit_ratio), c("tips", "Orbit_ratio", "max_crep", "Family", "fam_colours")]

#perform the one-way anova
art_model <- aov(Orbit_ratio ~ max_crep, data = trait.data.art)
summary(art_model)

#perform the post-hoc tukey test
TukeyHSD(art_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
art_phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#cannot calculate the pairwise corrected p-values likley because the cathemeral sample size is too small
#remove cathemeral if artio only and rerun
#trait.data.art <- filter(trait.data.art, max_crep != "cathemeral")
#art_phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#add p values manually
library(rstatix)
#for artio
stat.test <- data.frame(group1 = c("crepuscular", "crepuscular", "diurnal"),
                        group2 = c("diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(art_phylANOVA$Pt[2], art_phylANOVA$Pt[3], art_phylANOVA$Pt[6]),
                        y.position = c(0.95, 1.0, 0.994))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

boxplot_art <- ggplot(trait.data.art, aes(x = max_crep, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = max_crep), alpha=0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Corneal diameter: axial length") + scale_fill_manual(values=unique(trait.data.art$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 1.01, method = "anova") + annotate("text", x = 1.15, y = 1.007, label = paste("phylANOVA, p =", art_phylANOVA$Pf)) +
  stat_pvalue_manual(stat.test, label = "p.adj")
boxplot_art

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Orbit_ratio", "_boxplots_anova_artio.pdf"), width = 8, height = 7)
boxplot_art
dev.off()

#only ruminants
trait.data.art %>% filter(Family %in% c("Bovidae", "Cervidae", "Giraffidae", "Tragulidae")) %>% 
  ggplot(., aes(x = max_crep, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = max_crep), alpha=0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Corneal diameter: axial length") + scale_fill_manual(values=unique(trait.data.art$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 1.01, method = "anova") + annotate("text", x = 1.15, y = 1.007, label = paste("phylANOVA, p =", art_phylANOVA$Pf)) +
  stat_pvalue_manual(stat.test, label = "p.adj")

### Artiodactyla body mass 
trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$AdultBodyMass_g), c("tips", "AdultBodyMass_g", "max_crep", "Family", "fam_colours")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ]

#perform the one-way anova
art_model <- aov(AdultBodyMass_g ~ max_crep, data = trait.data.art)
summary(art_model)

#perform the post-hoc tukey test
TukeyHSD(art_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
art_phylANOVA <- calculatePhylANOVA(trait.data.art, "AdultBodyMass_g")

boxplot_art <- ggplot(trait.data.art, aes(x = max_crep, y = log(AdultBodyMass_g))) +
  geom_boxplot(aes(fill = max_crep), alpha=0.8, outlier.shape = NA) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Body mass (g)") + scale_fill_manual(values=unique(trait.data.art$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 14, method = "anova") + annotate("text", x = 1.15, y = 14.5, label = paste("phylANOVA, p =", art_phylANOVA$Pf))
boxplot_art   

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Body_mass", "_boxplots_anova_artio.pdf"), width = 8, height = 7)
boxplot_art
dev.off()

### Artiodactyla latitude
trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$GR_MidRangeLat_dd), c("tips", "GR_MidRangeLat_dd", "max_crep", "Family", "fam_colours")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ]

#perform the one-way anova
art_model <- aov(GR_MidRangeLat_dd ~ max_crep, data = trait.data.art)
summary(art_model)

#perform the post-hoc tukey test
TukeyHSD(art_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
art_phylANOVA <- calculatePhylANOVA(trait.data.art, "GR_MidRangeLat_dd")

boxplot_art <- ggplot(trait.data.art, aes(x = max_crep, y = log(GR_MidRangeLat_dd))) +
  geom_boxplot(aes(fill = max_crep), alpha=0.8, outlier.shape = NA) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Mid latitude range") + scale_fill_manual(values=unique(trait.data.art$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 5, method = "anova") + annotate("text", x = 1.15, y = 5.5, label = paste("phylANOVA, p =", art_phylANOVA$Pf))
boxplot_art   

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "mean_latitude", "_boxplots_anova_artio.pdf"), width = 8, height = 7)
boxplot_art
dev.off()

### Artiodactyla max latitude
trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$GR_MaxLat_dd), c("tips", "GR_MaxLat_dd", "max_crep", "Family", "fam_colours")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ]

#perform the one-way anova
art_model <- aov(GR_MaxLat_dd ~ max_crep, data = trait.data.art)
summary(art_model)

#perform the post-hoc tukey test
TukeyHSD(art_model, conf.level = .95)

#perform the phylogenetically corrected one-way anova
art_phylANOVA <- calculatePhylANOVA(trait.data.art, "GR_MaxLat_dd")

boxplot_art <- ggplot(trait.data.art, aes(x = max_crep, y = log(GR_MaxLat_dd))) +
  geom_boxplot(aes(fill = max_crep), alpha=0.8, outlier.shape = NA) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "GR_MaxLat_dd") + scale_fill_manual(values=unique(trait.data.art$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.y = 5, method = "anova") + annotate("text", x = 1.15, y = 5.5, label = paste("phylANOVA, p =", art_phylANOVA$Pf))
boxplot_art   

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "max_latitude", "_boxplots_anova_artio.pdf"), width = 8, height = 7)
boxplot_art
dev.off()

# Save out all discrete traits plots ---------------------------------------------------------

#### Cetacean diet
trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$max_crep),]
trait.data <- trait.data[!is.na(trait.data$Diet),]

# #filter for species in the final tree
# trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
# #custom.colours.2 <- c( "grey90", "grey40", "grey66", "black", "red")
# diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "Diet")]
# diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
# diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Diet), inherit.aes = FALSE, colour = "transparent") + scale_fill_grey(name = Diet)
# diel.plot <- diel.plot
# diel.plot

test_df <- trait.data[, c("max_crep", "Diet")]
test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 


####Cetacean habitat
trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$max_crep),]
trait.data <- trait.data[!is.na(trait.data$Habitat),]

# #filter for species in the final tree
# trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
# #custom.colours.2 <- c( "grey90", "grey40", "grey66", "black", "red")
# diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "Habitat")]
# diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
# diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Habitat), inherit.aes = FALSE, colour = "transparent") + scale_fill_grey(name = "Habitat")
# diel.plot <- diel.plot
# diel.plot

test_df <- trait.data[, c("max_crep", "Habitat")]
test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 

####Cetacean feeding mechanism 

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

# "Feeding_method" takes data from 4 databases
# "Prey_capture" #less species than feeding method (60 vs 76) but come from one source (Churchill)

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$max_crep),]
trait.data <- trait.data[!is.na(trait.data$Feeding_method),]

# #filter for species in the final tree
# trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
# #custom.colours.2 <- c( "grey90", "grey40", "grey66", "black", "red")
# diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "Feeding_method")]
# diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
# diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Feeding_method), inherit.aes = FALSE, colour = "transparent") + scale_fill_grey(name = "Feeding method")
# diel.plot <- diel.plot
# diel.plot

test_df <- trait.data[, c("max_crep", "Feeding_method")]
test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 


#### Terresterial vs aquatic 

#test to see if aquatic lifestyle is associated with cathemerality
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$max_crep),]
artio_full$enviro <- "aquatic"
for(i in 1:nrow(artio_full)){
  if(artio_full[i, "Parvorder"] == "non-cetacean"){
    artio_full[i, "enviro"] <- "terrestrial"
  }
}

test_df <- artio_full[, c("max_crep", "enviro")]
test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 9)))


#### Habitat and diet terrestrial artiodactyla

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep),]

trait.data.1 <- trait.data.art[!is.na(trait.data.art$TrophicLevel),]
test_df <- trait.data.1[, c("max_crep", "TrophicLevel")]
test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 9)))

trait.data.1 <- trait.data.art[!is.na(trait.data.art$DietBreadth),]
test_df <- trait.data.1[, c("max_crep", "DietBreadth")]
test <- fisher.test(table(test_df))
mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 9)))

### Possible multivariate associations

ggplot(trait.data.art, aes(x = max_crep, y = log(AdultBodyMass_g))) + geom_boxplot() + facet_wrap(~TrophicLevel)
