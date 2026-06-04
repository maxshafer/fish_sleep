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
                        legend.position = "none",
                        panel.border = element_rect(colour = "black", fill = "transparent"),
                        axis.text = element_text(size = 9), axis.title = element_text(size = 11))  

#save out the standard plot elements
boxplot_format <- list(geom_boxplot(aes(fill = max_crep), alpha = 0.5, outlier.shape = NA),
                       scale_fill_manual(values = custom.colours, guide = "none"),
                       geom_jitter(aes(fill = max_crep), size = 4, width = 0.1, height = 0, colour = "black", pch = 21, alpha = 0.8), boxplot_theme,
                       scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")))


# Cetacean orbit ratio ----------------------------------------------------

trait.data.1 <- trait.data[!is.na(trait.data$Orbit_ratio),]

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "Orbit_ratio")

#set the family colours for consistency
#family.colours <- c("grey", "black","dodgerblue", "royalblue3",  "cyan","darkgreen", "springgreen", "darkolivegreen1", "gold", "orange", "firebrick2", "brown", "orchid1", "darkorchid")

#plot out the group means
cet_orbit_boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = Orbit_ratio)) + 
  boxplot_format +
  labs(x = "Temporal activity pattern", y = "Relative eye size") + 
  annotate("text", x = 1.4, y = 40, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

cet_orbit_boxplot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Orbit_ratio_boxplot_cetaceans.pdf", width = 7.5, height = 7)
# boxplot
# dev.off()

# Cetacean mean latitude --------------------------------------------------

#family.colours <- c("dodgerblue",  "cyan","darkgreen", "springgreen", "darkolivegreen1", "gold", "orange", "firebrick2", "brown", "orchid1", "darkorchid")

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$mean_lat),]
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep),]

trait.data.1 <- filter(trait.data.1, Parvorder == "Odontoceti")

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "mean_lat")

cet_lat_boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = mean_lat)) + 
  boxplot_format +
  labs(x = "Temporal activity pattern", y = "Mean latitude range") + 
  annotate("text", x = 1.4, y = 100, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Mean_latitude_boxplot_cetaceans.pdf", width = 7, height = 7.5)
# boxplot
# dev.off()

# Cetacean max latitude ---------------------------------------------------
#family.colours <- c("dodgerblue",  "cyan","darkgreen", "springgreen", "darkolivegreen1", "gold", "orange", "firebrick2", "brown", "orchid1", "darkorchid")

#select the variable we're looking at
trait.data.1 <- trait.data[!is.na(trait.data$max_lat),]
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep),]

trait.data.1 <- filter(trait.data.1, Parvorder == "Odontoceti")

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "max_lat")

cet_maxlat_boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = max_lat)) + 
  boxplot_format +
  labs(x = "Temporal activity pattern", y = "Maximum latitude range") + 
  annotate("text", x = 1.4, y = 100, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Max_latitude_boxplot_cetaceans.pdf", width = 7, height = 7.5)
# cet_max_latboxplot
# dev.off()

# Cetacean body mass ------------------------------------------------------

trait.data.1 <- trait.data[!is.na(trait.data$Body_mass_kg),]

#set the family colours for consistency
#family.colours <- c("grey", "black","dodgerblue", "royalblue3",  "cyan","darkgreen", "springgreen", "darkolivegreen1", "gold", "orange", "firebrick2", "brown", "orchid1", "darkorchid")

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "Body_mass_kg")

cet_mass_boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = log(Body_mass_kg))) + 
  boxplot_format +
  labs(x = "Temporal activity pattern", y = "Log (body mass (kg))") + 
  annotate("text", x = 1.4, y = 13, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Body_mass_boxplots_anova_cetaceans.pdf", width = 8, height = 7)
# cet_mass_boxplot
# dev.off()

# Cetacean dive depth -----------------------------------------------------
#family.colours <- c("grey", "black","dodgerblue", "cyan", "darkgreen","darkolivegreen1", "orange", "firebrick2", "orchid1", "darkorchid")

trait.data.1 <- trait.data[!is.na(trait.data$Dive_depth_m),]

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.1, "Dive_depth_m")

cet_dive_boxplot <- ggplot(trait.data.1, aes(x = max_crep, y = log(Dive_depth_m))) + 
  boxplot_format +
  labs(x = "Temporal activity pattern", y = "Log (maximum dive depth (m))") + 
  annotate("text", x = 1.4, y = 9, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

cet_dive_boxplot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Dive_depth_boxplot_cetaceans.pdf", width = 7.5, height = 7)
# boxplot
# dev.off()

# Ruminant orbit size -----------------------------------------------------

#family.colours <- c("dodgerblue", "springgreen",  "gold", "darkorchid")

#load in the data
trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$Orbit_ratio), c("tips", "Orbit_ratio", "max_crep", "Family")]

#filter for just ruminants
trait.data.art <-  filter(trait.data.art, Family %in% c("Bovidae", "Cervidae", "Giraffidae", "Tragulidae"))

#calculate phylogenetic anova
phylANOVA <- calculatePhylANOVA(trait.data.art, "Orbit_ratio")

#cannot calculate the pairwise corrected p-values likely because the cathemeral sample size is too small
#remove cathemeral and rerun
trait.data.art.subset <- filter(trait.data.art, max_crep != "cathemeral")
phylANOVA <- calculatePhylANOVA(trait.data.art.subset, "Orbit_ratio")

#add p values manually
stat.test <- data.frame(group1 = c("crepuscular", "crepuscular", "diurnal"),
                        group2 = c("diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[6]),
                        y.position = c(0.99, 1, 1.01))

stat.test <- stat.test %>% add_x_position(x = "max_crep")
stat.test <- stat.test %>% mutate(xmin = xmin + 1, xmax = xmax + 1)

rum_orbit_boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = Orbit_ratio)) + 
  boxplot_format +
  stat_pvalue_manual(stat.test, label = "p.adj", hjust = 0.7) +
  labs(x = "Temporal activity pattern", y = "Relative eye size") + 
  annotate("text", x = 1.4, y = 1.01, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

rum_orbit_boxplot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Orbit_ratio_boxplot_ruminants.pdf", width = 7, height = 7.5)
# rum_orbit_boxplot
# dev.off()

# Ruminant mean latitude -------------------------------------------------------
#family.colours <- c("black","dodgerblue", "springgreen",  "gold", "firebrick2", "darkorchid")

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
#use below for IUCN data
trait.data.art <- trait.data.art[!is.na(trait.data.art$mean_lat), c("tips", "mean_lat", "max_crep", "Family")]
#use below for pantheria data
#trait.data.art <- trait.data.art[!is.na(trait.data.art$GR_MidRangeLat_dd), c("tips", "GR_MidRangeLat_dd", "max_crep", "Family")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ] %>% filter(Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))

#perform the phylogenetically corrected one-way anova
#phylANOVA <- calculatePhylANOVA(trait.data.art, "GR_MidRangeLat_dd")
phylANOVA <- calculatePhylANOVA(trait.data.art, "mean_lat")

#add p values manually
stat.test <- data.frame(group1 = c("cathemeral", "cathemeral", "cathemeral", "crepuscular", "crepuscular", "diurnal"),
                        group2 = c("crepuscular", "diurnal", "nocturnal", "diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[4], phylANOVA$Pt[7], phylANOVA$Pt[8], phylANOVA$Pt[12]),
                        y.position = c(80, 90, 100, 110, 120, 135))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

rum_lat_boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = mean_lat)) + 
  boxplot_format +
  stat_pvalue_manual(stat.test, label = "p.adj", hjust = 0.7) +
  labs(x = "Temporal activity pattern", y = "Mean latitude range") + 
  annotate("text", x = 1.4, y = 140, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

rum_lat_boxplot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/mean_latitude_boxplots_ruminant.pdf", width = 7, height = 7.5)
# boxplot
# dev.off()


# Ruminant max latitude ---------------------------------------------------
#family.colours <- c("black","dodgerblue", "springgreen",  "gold", "firebrick2", "darkorchid")

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
#use below for IUCN data
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_lat), c("tips", "max_lat", "max_crep", "Family")]
#use below for pantheria data
#trait.data.art <- trait.data.art[!is.na(trait.data.art$GR_MaxLat_dd), c("tips", "GR_MaxLat_dd", "max_crep", "Family")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ] %>% filter(Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))

#perform the phylogenetically corrected one-way anova
#phylANOVA <- calculatePhylANOVA(trait.data.art, "GR_MaxLat_dd")
phylANOVA <- calculatePhylANOVA(trait.data.art, "max_lat")

#add p values manually
stat.test <- data.frame(group1 = c("cathemeral", "cathemeral", "cathemeral", "crepuscular", "crepuscular", "diurnal"),
                        group2 = c("crepuscular", "diurnal", "nocturnal", "diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[4], phylANOVA$Pt[7], phylANOVA$Pt[8], phylANOVA$Pt[12]),
                        y.position = c(80, 90, 100, 110, 120, 135))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

rum_max_lat_boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = max_lat)) + 
  boxplot_format +
  stat_pvalue_manual(stat.test, label = "p.adj", hjust = 0.7) +
  labs(x = "Temporal activity pattern", y = "Maximum latitude range") + 
  annotate("text", x = 1.4, y = 135, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

rum_max_lat_boxplot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_latitude_boxplots_ruminant.pdf", width = 7.5, height = 7)
# rum_max_lat_boxplot
# dev.off()


# Ruminant body mass ------------------------------------------------------
#family.colours <- c("black","dodgerblue", "springgreen",  "gold", "firebrick2", "darkorchid")

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$AdultBodyMass_g), c("tips", "AdultBodyMass_g", "max_crep", "Family", "Diel_Pattern")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ] %>% filter(Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))

trait.data.art <- trait.data.art %>% mutate(AdultBodyMass_kg = AdultBodyMass_g/1000)

#perform the phylogenetically corrected one-way anova
phylANOVA <- calculatePhylANOVA(trait.data.art, "AdultBodyMass_kg")

rum_mass_boxplot <- ggplot(trait.data.art, aes(x = max_crep, y = log(AdultBodyMass_kg))) + 
  boxplot_format +
  labs(x = "Temporal activity pattern", y = "Log (body mass (kg))") + 
  annotate("text", x = 1.4, y = 7.5, label = paste("phylANOVA, p =", phylANOVA$Pf)) 

# pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Body_mass_boxplots_anova_ruminants.pdf"), width = 8, height = 7)
# rum_mass_boxplot
# dev.off()



# Arranging and saving out the plots --------------------------------------

boxplots <- (cet_orbit_boxplot + rum_orbit_boxplot) /
  (cet_lat_boxplot + rum_lat_boxplot) /
  (cet_dive_boxplot + plot_spacer()) +
  plot_annotation(tag_levels = 'a') + theme(plot.tag = element_text(size = 14))

boxplots

#alternative plot arrangement
boxplots <- (cet_orbit_boxplot + cet_lat_boxplot + cet_mass_boxplot) /
  (rum_orbit_boxplot + rum_lat_boxplot + rum_mass_boxplot) /
  plot_annotation(tag_levels = 'a') + theme(plot.tag = element_text(size = 14))

boxplots


pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/ecomorphological_boxplots.pdf"), width = 8.5, height = 8)
boxplots
dev.off()


