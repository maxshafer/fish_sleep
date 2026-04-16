source("scripts/fish_sleep_functions.R")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))


# Function to create trait vector for phylANOVA -------------------------

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
  phylANOVA <- phytools::phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
  return(phylANOVA)
}


# Plot formatting -------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

#filter for species with activity pattern data
trait.data <- trait.data[!is.na(trait.data$max_crep),]

#save formatting
custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
boxplot_theme <-  theme_minimal() +
  theme(panel.background = element_rect(fill='transparent', colour = "transparent"), 
                        plot.background = element_rect(fill='transparent', color=NA), 
                        legend.background = element_rect(fill='transparent'),
                        legend.position = "bottom",
                        panel.border = element_rect(colour = "black", fill = "transparent"),
                        axis.text = element_text(size = 11), axis.title = element_text(size = 12))  

# Cetacean orbit ratio ----------------------------------------------------

trait.data.1 <- trait.data[!is.na(trait.data$Orbit_ratio),]

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "Orbit_ratio")

#set the family colours for consistency
family.colours <- c("grey", "black","dodgerblue", "royalblue3",  "cyan","darkgreen", "springgreen", "darkolivegreen1", "gold", "orange", "firebrick2", "brown", "orchid1", "darkorchid")

#plot out the group means

boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Proportional orbit size\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 40, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme + 
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 
  
boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Orbit_ratio_boxplot_cetaceans.pdf", width = 7.5, height = 7)
boxplot
dev.off()

# Cetacean mean latitude --------------------------------------------------

family.colours <- c("dodgerblue", "darkolivegreen1", "orange")

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$mean_lat),]
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep),]

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "mean_lat")

boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = mean_lat)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Mean latitude range\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 7500000, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme + 
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Mean_latitude_boxplot_cetaceans.pdf", width = 7.5, height = 7)
boxplot
dev.off()

# Cetacean max latitude ---------------------------------------------------
family.colours <- c("dodgerblue", "darkolivegreen1", "orange")

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$max_lat),]
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep),]

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "max_lat")

boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = max_lat)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Maximum latitude range\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 7500000, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme + 
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Max_latitude_boxplot_cetaceans.pdf", width = 7.5, height = 7)
boxplot
dev.off()

# Cetacean dive depth -----------------------------------------------------
family.colours <- c("grey", "black","dodgerblue", "cyan", "darkgreen","darkolivegreen1", "orange", "firebrick2", "orchid1", "darkorchid")

trait.data.1 <- trait.data[!is.na(trait.data$Dive_depth_m),]

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "Dive_depth_m")

boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Maximum dive depth (m)\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 3300, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme + 
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Dive_depth_boxplot_cetaceans.pdf", width = 7.5, height = 7)
boxplot
dev.off()

# Delphinidae dive depth --------------------------------------------------
trait.data.1 <- trait.data[!is.na(trait.data$Dive_depth_m),]
trait.data.delph <- trait.data.1 %>% filter(Family == "Delphinidae")

phylANOVA <- calculatePhylANOVA(trait.data.delph, "Dive_depth_m")

stat.test <- data.frame(group1 = c("cathemeral", "cathemeral", "cathemeral", "crepuscular", "crepuscular", "diurnal"),
                        group2 = c("crepuscular", "diurnal", "nocturnal", "diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[4], phylANOVA$Pt[7], phylANOVA$Pt[8], phylANOVA$Pt[12]),
                        y.position = c(7.2, 7.6, 8, 8.4, 8.8, 9.3))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

boxplot <- ggplot(trait.data.delph, aes(x = max_crep, y = log(Dive_depth_m))) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Log (maximum dive depth (m))\n") + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = 'black', fill = "dodgerblue", pch = 21) +
  annotate("text", x = 1.15, y = 9, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme +
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Dive_depth_boxplot_delphinidae.pdf", width = 7, height = 7.5)
boxplot
dev.off()



# Ruminant orbit size -----------------------------------------------------

family.colours <- c("dodgerblue", "springgreen",  "gold", "darkorchid")

#load in the data
trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$Orbit_ratio), c("tips", "Orbit_ratio", "max_crep", "Family", "fam_colours")]

#filter for just ruminants
trait.data.art <-  filter(trait.data.art, Family %in% c("Bovidae", "Cervidae", "Giraffidae", "Tragulidae"))

#calculate phylogenetic anova
phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#cannot calculate the pairwise corrected p-values likely because the cathemeral sample size is too small
#remove cathemeral and rerun
trait.data.art <- filter(trait.data.art, max_crep != "cathemeral")
phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#add p values manually
stat.test <- data.frame(group1 = c("crepuscular", "crepuscular", "diurnal"),
                        group2 = c("diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[6]),
                        y.position = c(0.99, 1, 1.01))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = Orbit_ratio)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Corneal diameter: axial length\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 1.01, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme + 
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_x_discrete(labels = c("crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Orbit_ratio_boxplot_ruminants.pdf", width = 7, height = 7.5)
boxplot
dev.off()

# Ruminant mean latitude -------------------------------------------------------
family.colours <- c("black","dodgerblue", "springgreen",  "gold", "firebrick2", "darkorchid")

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$GR_MidRangeLat_dd), c("tips", "GR_MidRangeLat_dd", "max_crep", "Family", "fam_colours")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ] %>% filter(Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.art, "GR_MidRangeLat_dd")

#add p values manually
stat.test <- data.frame(group1 = c("cathemeral", "cathemeral", "cathemeral", "crepuscular", "crepuscular", "diurnal"),
                        group2 = c("crepuscular", "diurnal", "nocturnal", "diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[4], phylANOVA$Pt[7], phylANOVA$Pt[8], phylANOVA$Pt[12]),
                        y.position = c(80, 90, 100, 110, 120, 130))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = GR_MidRangeLat_dd)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Mean latitude range\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 115, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme +
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal"))

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/mean_latitude_boxplots_ruminant.pdf", width = 7, height = 7.5)
boxplot
dev.off()


# Ruminant max latitude ---------------------------------------------------
family.colours <- c("black","dodgerblue", "springgreen",  "gold", "firebrick2", "darkorchid")

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$GR_MaxLat_dd), c("tips", "GR_MaxLat_dd", "max_crep", "Family", "fam_colours")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ] %>% filter(Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.art, "GR_MaxLat_dd")

#add p values manually
stat.test <- data.frame(group1 = c("cathemeral", "cathemeral", "cathemeral", "crepuscular", "crepuscular", "diurnal"),
                        group2 = c("crepuscular", "diurnal", "nocturnal", "diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[4], phylANOVA$Pt[7], phylANOVA$Pt[8], phylANOVA$Pt[12]),
                        y.position = c(80, 90, 100, 110, 120, 130))

stat.test <- stat.test %>% add_x_position(x = "max_crep")


boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = GR_MaxLat_dd)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "\nTemporal activity pattern", y = "Maximum latitude range\n") + 
  scale_fill_manual(values=family.colours)  + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) +
  annotate("text", x = 1.15, y = 115, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme +
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal"))

boxplot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_latitude_boxplots_ruminant.pdf", width = 7.5, height = 7)
boxplot
dev.off()
