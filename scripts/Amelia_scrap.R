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


# Multivariate analysis ---------------------------------------------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

trait.data <- trait.data %>% filter(!is.na(max_crep) & !is.na(Habitat) & !is.na(Dive_depth_m))

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



#From section 6.1
# Artiodactyla latitude ---------------------------------------------------

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


# Artiodactyla max latitude -----------------------------------------------

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


# Artiodactyla orbit size -------------------------------------------------

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


# Artiodactyla body mass --------------------------------------------------

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




# Cetacean body mass ------------------------------------------------------

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


# Ruminant body mass ------------------------------------------------------

trait.data.art <- read.csv(here("artiodactyla_ecomorphology_dataset.csv"))
trait.data.art <- trait.data.art[!is.na(trait.data.art$AdultBodyMass_g), c("tips", "AdultBodyMass_g", "max_crep", "Family", "Diel_Pattern")]
trait.data.art <- trait.data.art[!is.na(trait.data.art$max_crep), ] %>% filter(Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))

#perform the phylogenetically corrected one-way anova
art_phylANOVA <- calculatePhylANOVA(trait.data.art, "AdultBodyMass_g")

boxplot_art <- ggplot(trait.data.art, aes(x = Diel_Pattern, y = log(AdultBodyMass_g))) +
  geom_boxplot(aes(fill = max_crep), alpha=0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours) +
  new_scale_fill() + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Body mass (g)") +
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  annotate("text", x = 1.15, y = 14.5, label = paste("phylANOVA, p =", art_phylANOVA$Pf)) 
boxplot_art   

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Body_mass", "_boxplots_anova_ruminants.pdf"), width = 8, height = 7)
boxplot_art
dev.off()

#filter for species in the final tree
# trait.data <- trait.data.art[trait.data.art$tips %in% mam.tree$tip.label,]
# trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
# 
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
# diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep", "AdultBodyMass_g")]
# diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
# diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = log(AdultBodyMass_g)), inherit.aes = FALSE, colour = "transparent") + scale_fill_viridis_c(option = "magma")
# diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 4) 
# diel.plot


# Cetacean dive data in detail ---------------------------------------------------

dive.data <- read.csv(here("cetacean_dive_depth_all_sources.csv"))
trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

dive.data <- merge(dive.data, trait.data[, -c(14)], by = "tips", all = TRUE)

dive.data %>% ggplot(., aes(x = Diel_Pattern, y = Final_dive_depth)) + geom_boxplot() + geom_point() + facet_wrap(~Parvorder)

#how deep is the deep scattering layer? 
#Wikipedia says they rise to 100m at night and descend to 800-1000m during day
#coastal species also wouldn't be involved in this since they aren't diving that deep

dive.data %>% filter(Habitat %in% c("coastal/pelagic", "pelagic")) %>% filter(!is.na(max_crep)) %>%
  ggplot(., aes(x = max_crep, y = Mean_dive_depth)) + geom_boxplot() +
  geom_point() + facet_wrap(~Parvorder) + stat_compare_means(method = "anova")

dive.data %>% filter(Habitat %in% c("coastal/pelagic", "pelagic")) %>% filter(!is.na(max_crep)) %>%
  ggplot(., aes(x = max_crep, y = Final_dive_depth)) + geom_boxplot() +
  geom_point() + facet_wrap(~Parvorder) + stat_compare_means(method = "anova")


#check for phylogenetic significance

#maximum dive depth
trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))
trait.data <- trait.data[!is.na(trait.data$max_crep),]

trait.data <- filter(trait.data, Habitat %in% c("coastal/pelagic", "pelagic"))

trait.data.od <- filter(trait.data, Parvorder == "Odontoceti")
trait.data.od <- filter(trait.data, Family == "Delphinidae")

trait.data.od <- trait.data.od[!is.na(trait.data.od$Dive_depth_m),]

#perform the phylogenetically corrected one-way anova
odonto_phylANOVA <- calculatePhylANOVA(trait.data.od, "Dive_depth_m")

boxplot_dive <- ggplot(trait.data.od, aes(x = max_crep, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Dive_depth_m") + scale_fill_manual(values=unique(trait.data$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.x = 0.8, label.y = 3400, method = "anova") + 
  annotate("text", x = 1.15, y = 3300, label = paste("phylANOVA, p =", odonto_phylANOVA$Pf))
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Maximum_dive_depth", "_boxplots_anova_pelagic_odontocetes.pdf"), width = 8, height = 7)
boxplot_dive
dev.off()

#with mean dive depth 
trait.data.od <- filter(trait.data, Parvorder == "Odontoceti")
trait.data.od <- trait.data.od[!is.na(trait.data.od$Mean_dive_depth_m),]

#perform the phylogenetically corrected one-way anova
odonto_phylANOVA <- calculatePhylANOVA(trait.data.od, "Mean_dive_depth_m")

boxplot_dive <- ggplot(trait.data.od, aes(x = max_crep, y = Mean_dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Mean_dive_depth_m") + scale_fill_manual(values=unique(trait.data$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.x = 0.8, label.y = 3400, method = "anova") + 
  annotate("text", x = 1.15, y = 3300, label = paste("phylANOVA, p =", odonto_phylANOVA$Pf))
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Mean_dive_depth", "_boxplots_anova_pelagic_odontocetes.pdf"), width = 8, height = 7)
boxplot_dive
dev.off()

#both parvorders on same plot
trait.data.od <- filter(trait.data, Parvorder == "Odontoceti")
trait.data.od <- trait.data.od[!is.na(trait.data.od$Dive_depth_m),]

trait.data.my <- filter(trait.data, Parvorder == "Mysticeti")
trait.data.my <- trait.data.my[!is.na(trait.data.my$Dive_depth_m),]

#perform the phylogenetically corrected one-way anova
odonto_phylANOVA <- calculatePhylANOVA(trait.data.od, "Dive_depth_m")
mystic_phylANOVA <- calculatePhylANOVA(trait.data.my, "Dive_depth_m")

dat_text <- data.frame(label = c(paste("phylANOVA =", mystic_phylANOVA$Pf, sep = " "), paste("phylANOVA =", odonto_phylANOVA$Pf, sep = " ")), Parvorder = c("Mysticeti", "Odontoceti"))

boxplot_dive <- ggplot(trait.data, aes(x = max_crep, y = Dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Maximum dive depth (m)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.x = 0.8, label.y = 3400, method = "anova") + 
  geom_text(data= dat_text,mapping = aes(x = 1, y = 3300, label = label)) +
  facet_wrap(~Parvorder, scales = "free_x")
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Maximum_dive_depth", "_boxplots_anova_pelagic_odontocetes_mysticeties.pdf"), width = 8, height = 7)
boxplot_dive
dev.off()

#both parvorders on same plot for mean
trait.data.od <- filter(trait.data, Parvorder == "Odontoceti")
trait.data.od <- trait.data.od[!is.na(trait.data.od$Mean_dive_depth_m),]

trait.data.my <- filter(trait.data, Parvorder == "Mysticeti")
trait.data.my <- trait.data.my[!is.na(trait.data.my$Mean_dive_depth_m),]

#perform the phylogenetically corrected one-way anova
odonto_phylANOVA <- calculatePhylANOVA(trait.data.od, "Mean_dive_depth_m")
mystic_phylANOVA <- calculatePhylANOVA(trait.data.my, "Mean_dive_depth_m")

dat_text <- data.frame(label = c(paste("phylANOVA =", mystic_phylANOVA$Pf, sep = " "), paste("phylANOVA =", odonto_phylANOVA$Pf, sep = " ")), Parvorder = c("Mysticeti", "Odontoceti"))

boxplot_dive <- ggplot(trait.data, aes(x = max_crep, y = Mean_dive_depth_m)) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8) + scale_fill_manual(values = custom.colours) +
  new_scale_fill() + geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Mean dive depth (m)") + scale_fill_manual(values=unique(trait.data.1$fam_colours))  + 
  theme_minimal() + theme(panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'), panel.border = element_rect(colour = "black", fill = "transparent")) + 
  stat_compare_means(label.x = 0.8, label.y = 3400, method = "anova") + 
  geom_text(data= dat_text,mapping = aes(x = 1, y = 3300, label = label)) +
  facet_wrap(~Parvorder, scales = "free_x")
boxplot_dive

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Mean_dive_depth", "_boxplots_anova_pelagic_odontocetes_mysticeties.pdf"), width = 8, height = 7)
boxplot_dive
dev.off()

# Cetacean dive depth parvorder -------------------------------------------

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
#From section 5.0
# Section 7: Comparison of rate magnitude -------------------------------

#get rates
rates_df1 <- plot1kTransitionRates4state(readRDS(here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df1 <- filter(rates_df1, model == "Bridge_only")
rates_df2 <- plot1kTransitionRates4state(readRDS(here("august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df2 <- filter(rates_df2, model == "Bridge_only")

#need to compare rates from the same trees 
#each model has 12 rates, filter for one model (ARD), label each tree 
rates_df1$tree_n <- rep(1:1000, each = 10)
#subtract the cetacean rate from the ruminant rate, is it faster (negative number) or slower (positive number)
rates_df2$tree_n <- rep(1:1000, each = 10)

rates_df <- cbind(rates_df1, rates_df2)
colnames(rates_df) <- c("whippo_rates", "model", "solution", "colours", "tree_n", "rumi_rates", "model", "solution", "colours", "tree_n")
rates_df <- rates_df[, c("whippo_rates", "solution", "tree_n", "rumi_rates")]
rates_df$difference <- rates_df$whippo_rates - rates_df$rumi_rates
#difference is negative or small -whippo is much faster, difference is positive or large ruminants have similar rates

ggplot(rates_df, aes(x = solution, y = difference)) + geom_jitter()
ggplot(rates_df, aes(x = whippo_rates, y = rumi_rates, colour = solution)) + geom_point() + facet_wrap(~solution)
ggplot(rates_df, aes(x = whippo_rates, fill = solution)) + geom_histogram() + facet_wrap(~solution)

rates_df %>% group_by(solution) %>% summarize(mean_rates = mean(whippo_rates)) %>% ggplot(., aes(x = mean_rates, y = solution)) + geom_point()
rates_df %>% group_by(solution) %>% summarize(mean_rates = mean(rumi_rates)) %>% ggplot(., aes(x = mean_rates, y = solution)) + geom_point()

rates_df %>% group_by(solution) %>% 
  summarize(mean_whippo_rates = mean(whippo_rates), mean_rumi_rates = mean(rumi_rates)) %>% 
  ggplot(., aes(y = solution)) + geom_point(aes(x=mean_whippo_rates), colour = "blue") +
  geom_point(aes(x = mean_rumi_rates), colour = "red")

#make a forest-ish plot
ggplot(rates_df, aes(y = solution)) + geom_point(aes(x=whippo_rates), colour = "blue") +
  geom_point(aes(x = rumi_rates), colour = "red")


rates_df1 <- plot1kTransitionRates4state(readRDS(here("august_whippomorpha_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df1 <- filter(rates_df1, model == "Bridge_only")
rates_df2 <- plot1kTransitionRates4state(readRDS(here("august_ruminants_four_state_max_crep_traits_ER_SYM_ARD_CONSYM_bridge_only_models.rds")), 5)
rates_df2 <- filter(rates_df2, model == "Bridge_only")

rates_df1$clade <- "whippomorpha"
rates_df2$clade <- "ruminants"

rates_df <- rbind(rates_df1, rates_df2)

df <- rates_df %>% group_by(clade, solution) %>% summarize(mean_rates = mean(rates), 
                                                           lci = t.test(rates, conf.level = 0.95)$conf.int[1],
                                                           uci = t.test(rates, conf.level = 0.95)$conf.int[2])

ggplot(df, aes(x = mean_rates, y = solution, colour = clade)) + geom_point() + geom_errorbar(aes(y = solution, xmin = lci, xmax =uci),width = 0.4)


# Section 8: Total garbage test ------------------------------------------

filename <- "whippomorpha_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
model_results <- readRDS(here(paste0(filename, ".rds")))

trait.data <- model_results$ER_model$data
table(trait.data$max_crep)

#since there are four trait states
#likelihood <- ((factorial(n))/(factorial(n1) * factorial(n2) * factorial(n3) *factorial(n4))) * (p1^n1)  * (p2^n2)  * (p3^n3)  * (p4^n4) 

#using the natural log
n1 = 27
n2 = 22
n3 = 7
n4 = 21
n = 77

lnL_garb = n1 * log(n1 / n) + n2 * log(n2 / n) + n3 * log(n3 / n) + n4 * log(n4 / n)
#ln likelihood is -99.926

#compared to the actual likelihood
model_results$bridge_only_model$loglik #-91.0821
model_results$ER_model$loglik #-105.561
model_results$SYM_model$loglik #-97.924
model_results$ARD_model$loglik #-91.025
model_results$CONSYM_model$loglik # -97.46

likelihood_metrics <- max_clade_metrics(readRDS(here(paste0(filename, ".rds"))))
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)
likelihood_metrics <- rbind(likelihood_metrics, data.frame(model = "Total garbage", log_likelihoods = lnL_garb, AICc_scores = NA,  AIC_scores = NA))
ggplot(likelihood_metrics, aes(y = log_likelihoods, x = model, fill = model))  + geom_bar(stat = "identity")

#so the log likelihood is similar, but the actual model is more likely (higher log lik)

#for ruminants
filename <- "ruminants_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models"
model_results <- readRDS(here(paste0(filename, ".rds")))

trait.data <- model_results$ER_model$data
table(trait.data$max_crep)

n1 = 21
n2 = 125
n3 = 34
n4 = 23
n = 203

lnL_garb = n1 * log(n1 / n) + n2 * log(n2 / n) + n3 * log(n3 / n) + n4 * log(n4 / n)

#the garbage ln likelihood is -219 which is different than our most likely model (-199.78)

model_results$bridge_only_model$loglik #-199.7884
model_results$ER_model$loglik # -230.505
model_results$SYM_model$loglik #-219.0621
model_results$ARD_model$loglik #-199.8458
model_results$CONSYM_model$loglik # -219.0807

likelihood_metrics <- max_clade_metrics(readRDS(here(paste0(filename, ".rds"))))
likelihood_metrics <- pivot_wider(likelihood_metrics, names_from = model_metric, values_from = model_value)
likelihood_metrics <- rbind(likelihood_metrics, data.frame(model = "Total garbage", log_likelihoods = lnL_garb, AICc_scores = NA,  AIC_scores = NA))
ggplot(likelihood_metrics, aes(y = log_likelihoods, x = model, fill = model))  + geom_bar(stat = "identity") + ggtitle(filename)

likelihood_metrics$length <- 100
ggplot(likelihood_metrics, aes(y = log_likelihoods, x = length, fill = model)) +
  geom_line() + ggtitle(filename)

# Section 7: McCurry et al latitude ---------------------------------------
#https://doi.org/10.1093/biolinnean/blac128

McCurry <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\McCurry_2023.xlsx")

McCurry <- as.data.frame(McCurry[, c(6, 8:51)])

test <- McCurry %>% pivot_longer(cols = !`Absolute latitude`, names_to = "Species", values_to = "count")

test <- test %>% filter(count > 0)
test$Species <- str_replace_all(test$Species, pattern = "'", replacement = "")
test$Species <- str_replace_all(test$Species, pattern = " ", replacement = "_")
colnames(test) <- c("Absolute_latitude", "Species", "count")

test <- test %>% group_by(Species) %>% summarize(max_lat = max(Absolute_latitude), mean_lat = mean(Absolute_latitude), min_lat = min(Absolute_latitude))
names <- read.csv(here("cetaceans_full.csv"))

names <- names %>% separate(col = tips, into = c("Genus", "Species"), sep = "_")
test <- test %>% separate(col = Species, into = c("Genus", "Species"), sep = "_")

test <- merge(names, test, by = "Species", all.y =TRUE)
#two species have the species name attentuata, glacialis, australis and hectori. Drop the duplicates
latitude_df <- test[-c(5,8,9,11,31), ]

latitude_df %>% filter(!is.na(max_crep)) %>% ggplot(., aes(x = max_crep, y = max_lat)) + geom_boxplot() + stat_compare_means(method = "anova")

latitude_df$tips <- str_replace(latitude_df$Species_name, pattern = " ", replacement = "_")
latitude_df <- latitude_df[, c("tips", "max_lat", "mean_lat", "min_lat")]

#check if species names are spelled correctly
latitude_df[!latitude_df$tips %in% mam.tree$tip.label,]

#save out 
write.csv(latitude_df, here("cetacean_latitude_df.csv"), row.names = FALSE)

# Section 6: Groot et al ignore for now --------------------------------------------------

groot <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Groot_et_al_2023.xlsx")
#this data is coded so need to decode it with the original paper
#contains trait data that is in the other dataframes
#lifespan, length, mass, brain mass, EQ, age to reproduction, group size, gestation, sociality, group foraging, learned foraging, communication

