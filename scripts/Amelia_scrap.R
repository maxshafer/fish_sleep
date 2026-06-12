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



# Section: Random diel plots ----------------------------------------------

#all diel patterns
custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "gold", "#66C2A5", "#A6D854","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+1.5, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent", width = 3) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot <- diel.plot + geom_tiplab(size = 3, offset = 3.2) 
diel.plot

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", clade_name, "_six_state_plot_labelled.pdf"), width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 


cetaceans_full <- read.csv(here("whippomorpha.csv"))
cetaceans_full <- cetaceans_full[!is.na(cetaceans_full$max_crep), ]
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- cetaceans_full[cetaceans_full$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#add clade labels
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Mysticeti")
findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 4, taxonomic_level_name = "Odontoceti")

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular", fill = "transparent") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + 
  geom_cladelab(node = 137, label = "Mysticeti", align = TRUE, geom = "label", offset=1, align=TRUE, offset.text=1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot <- diel.plot + 
  geom_cladelab(node = 77, label = "Odontoceti", align = FALSE, geom = "label", offset=1, align=FALSE, offset.text=1, hjust = 1, barsize=2, fontsize=3, fill = "grey", barcolour = "grey", textcolour = "black")
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/max_crep_plot_no_Na_unlabelled_cladelabels.pdf", bg = "transparent")
diel.plot
dev.off()

# Section 7: Mammal tree ----------------------------------------------------

#make plot of bennie et al data with the artiodactyla data replaced with my own
new_mammals <- read.csv(here("Bennie_mam_data.csv")) #data from Bennie et al, 2014, 4732 species
new_mammals <- new_mammals %>% filter(Order != "Artiodactyla") #4492 species (removes 240 artios)
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv")) 
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ] #317 species (82 cetaceans, 235 non cetaceans)
artio_full <- artio_full %>% select(Species_name, Order, Family, max_crep)
new_mammals <- new_mammals %>% select(Species_name, Order, Family, max_crep)
diel_full <- rbind(new_mammals, artio_full) #4809 species
diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label,] #should be 4,400 species in tree (other 400 misnamed)
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5","grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x+3, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 6) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/clade_name/sleepy_mammals_max_crep_plot_unlabelled.pdf", width = 9, height = 8, bg = "transparent")
diel.plot
dev.off() 

order_list <- trait.data %>% group_by(Order) %>% filter(n()>10)
order_list <- unique(order_list$Order)
node_labels <- lapply(order_list, function(x){findMRCANode2(phylo = trpy_n, trait.data = trait.data, taxonomic_level_col = 2, taxonomic_level_name = x)})
node_labels <- do.call(rbind.data.frame, node_labels)
node_labels$barsize <- 2
node_labels$vjust <- 0.5
node_labels_left <- node_labels[node_labels$clade_name %in% c("Lagomorpha", "Scandentia", "Primates", "Artiodactyla", "Carnivora", "Perissodactyla", "Dermoptera", "Pholidota"),]
#node_labels_left[node_labels_left$clade_name %in% c("Dermoptera", "Pholidota"), "vjust"] <- -2
#node_labels_left[node_labels_left$clade_name %in% c("Perissodactyla"), "vjust"] <- -2
node_labels_right <- node_labels[!node_labels$clade_name %in% c("Lagomorpha", "Scandentia", "Primates", "Artiodactyla", "Carnivora", "Perissodactyla", "Dermoptera", "Pholidota"),]
#node_labels_right[node_labels_right$clade_name %in% c("Proboscidea", "Monotremata"), "vjust"] <- -2
node_labels_right[node_labels_right$clade_name %in% c("Cingulata"), "vjust"] <- 0

diel.plot <- ggtree(trpy_n, layout = "circular", fill = "transparent") %<+% trait.data[,c("tips", "max_crep")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = max_crep), inherit.aes = FALSE, colour = "transparent", width = 6) + scale_fill_manual(values = custom.colours)
diel.plot <- diel.plot + geom_cladelab(barsize = 1.5, barcolor = "grey50", node = node_labels_left$node_number, label = node_labels_left$clade_name, hjust = 1, offset = 3, vjust = node_labels_left$vjust, offset.text = 2)
diel.plot <- diel.plot + geom_cladelab(barsize = 1.5, barcolor = "grey50", node = node_labels_right$node_number, label = node_labels_right$clade_name, offset = 3, vjust = node_labels_right$vjust, offset.text = 2)
diel.plot <- diel.plot + theme(legend.position = "none", panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
diel.plot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/mammals_max_crep_plot_cladelabels.pdf", bg = "transparent", width = 10, height = 10)
# diel.plot
# dev.off()




# Comparison to Maor and Bennie artio datasets ----------------------------

#compare my ruminant data to Maor and Bennie datasets
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- filter(artio_full, Family %in% c("Bovidae", "Cervidae", "Antilocapridae", "Giraffidae", "Tragulidae", "Moschidae"))
table(artio_full$max_crep)
table(artio_full$Diel_Pattern)

Bennie_mam_data <- read.csv(here("Bennie_mam_data.csv")) #data from Bennie et al, 2014
Bennie_mam_data <- Bennie_mam_data[Bennie_mam_data$Species_name %in% artio_full$Species_name,]
table(Bennie_mam_data$max_crep)

maor_mam_data <- read.csv(here("Maor_artio_full.csv")) #data from Maor et al, 2017
maor_mam_data <- maor_mam_data[maor_mam_data$tips %in% artio_full$tips,]
maor_mam_data$Diel_pattern <- str_replace(maor_mam_data$Diel_pattern, pattern = c("Cathemeral/Crepuscular"), replacement = c("Crepuscular"))
maor_mam_data$Diel_pattern <- str_replace(maor_mam_data$Diel_pattern, pattern = c("Diurnal/Crepuscular"), replacement = c("Crepuscular"))
maor_mam_data$Diel_pattern <- str_replace(maor_mam_data$Diel_pattern, pattern = c("Nocturnal/Crepuscular"), replacement = c("Crepuscular"))
maor_mam_data$Diel_pattern <- str_replace(maor_mam_data$Diel_pattern, pattern = c("Diurnal/Cathemeral"), replacement = c("Cathemeral"))
maor_mam_data$Diel_pattern <- str_replace(maor_mam_data$Diel_pattern, pattern = c("Nocturnal/Cathemeral"), replacement = c("Cathemeral"))
table(maor_mam_data$Diel_pattern)

# Section 9: Phylogenetic signal lambda -wrong  ----------

#requires a vector of the trait data in the same order as phy$tip.label
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
mam.tree <- keep.tip(mam.tree, tip = trait.data$tips)

#all branches need to have a positive length
#replace branches with length 0 with 1% of the 1% quantile (replace it with a number very close to zero)
mam.tree$edge.length[mam.tree$edge.length == 0] <- quantile(mam.tree$edge.length, 0.1)*0.1

#function phylo.signal isn't working so load in the base code
#rep doesn't do anything?
phylo.signal <- function(trait, phy, rep = 999) {
  if (length(attributes(factor(trait))$levels) == length(trait)) 
    stop("Are you sure this variable is categorical?")
  
  phy <- keep.tip(phy, tip = names(trait))
  
  # calculate likelihood corresponding to maximum likelihood value of lambda
  obs <- fitDiscrete(phy, trait, transform="lambda")
  
  # calculate likelihood of model with no phylogenetic signal
  #null <- fitDiscrete(transform(phylo, "lambda", 0), trait)
  # this wasn't working so I tested the likelihood of the lambda against a tree with traits randomly distributed
  trait.random <- sample(trait)
  names(trait.random) <- names(trait)
  null <- fitDiscrete(phy, trait.random, transform = "lambda")
  
  # calculate the likelihood ratio between the two models
  LLR <- -2*(null$opt$lnL - obs$opt$lnL)
  
  # what is the p value of this likelihood ratio?
  p <- pchisq(LLR, df=1, lower.tail=FALSE)
  
  return(data.frame(row.names=NULL, lambda=obs$opt$lambda, obs=obs$opt$lnL, null=null$opt$lnL, LLR=LLR, p=p))
}

#function to create a vector of trait data, with species in same order as in tree (mam.tree$tip.label)
makeTraitVector <- function(trait.data = trait.data, taxonomic_level = "Order", clade_name = "Primates"){
  trait.data <- trait.data[trait.data[,taxonomic_level] == clade_name, ]
  sps_order <- as.data.frame(mam.tree$tip.label)
  colnames(sps_order) <- "tips"
  sps_order$id <- 1:nrow(sps_order)
  trait.data <- merge(trait.data, sps_order, by = "tips")
  trait.data <- trait.data[order(trait.data$id), ]
  trait <- trait.data$max_crep
  names(trait) <- trait.data$tips
  return(trait)
}

#calculate signal for all mammals
trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data$Kingdom <- "Mammals"
trait <- makeTraitVector(trait.data = trait.data, taxonomic_level = "Kingdom", clade_name = "Mammals")
mam_signal <- phylo.signal(trait = trait, phy = mam.tree, rep = 999)
mam_signal$clade <- "Mammals"
row.names(mam_signal) <- "Mammals"

#calculate signal for artiodactyla suborders
trait.data <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.vector.list <- lapply(unique(trait.data$Suborder), function(x) makeTraitVector(trait = trait.data, taxonomic_level = "Suborder", clade_name = x))
names(trait.vector.list) <- unique(trait.data$Suborder)
phylo.sig.list <- lapply(trait.vector.list, function(x) phylo.signal(trait = x, phy = mam.tree, rep = 999))
suborder_df <- do.call(rbind.data.frame, phylo.sig.list)
suborder_df$clade <- row.names(suborder_df)

#write.csv(suborder_df, here("phylogenetic_signal_artio_suborder.csv"), row.names = FALSE)

#calculate phylogenetic signal for all families with more than 100 species
trait.data <- read.csv(here("Bennie_mam_data.csv"))
trait.data <- trait.data[, c("max_crep", "tips", "Order")]

#to use my artio data instead (or in addition to)
#trait.data <- filter(trait.data, Order != "Artiodactyla")
trait.data.1 <- read.csv(here("sleepy_artiodactyla_full.csv"))
trait.data.1 <- trait.data.1[!is.na(trait.data.1$max_crep), c("max_crep", "tips", "Order")]
trait.data.1$Order <- str_replace(trait.data.1$Order, pattern = "Artiodactyla", replacement = "Amelia_artiodactyla")

trait.data <- rbind(trait.data, trait.data.1) #4459 species

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,] #4303 mammals in final tree
table(trait.data$Order, trait.data$max_crep)

#filter for orders that have over x number of species
#trait.data <- trait.data %>% group_by(Order) %>% filter(n() > 70)

#filter for species with all four diel categories
trait.data <- trait.data %>% group_by(Order) %>% filter(length(unique(max_crep)) ==4)

trait.vector.list <- lapply(unique(trait.data$Order), function(x) makeTraitVector(trait = trait.data, taxonomic_level = , clade_name = x))
names(trait.vector.list) <- unique(trait.data$Order)

phylo.sig.list <- lapply(trait.vector.list, function(x) phylo.signal(trait = x, phy = mam.tree))
phylo.sig.df <- do.call(rbind.data.frame, phylo.sig.list)
phylo.sig.df$clade <- row.names(phylo.sig.df)

#save out 
#write.csv(phylo.sig.df, here("phylogenetic_signal_mammals.csv"), row.names = FALSE)

#final comparison
phylo_signal_final <- rbind(phylo.sig.df, suborder_df, mam_signal)
write.csv(phylo_signal_final, here("phylogenetic_signal.csv"), row.names = FALSE)

#i think chiroptera only has nocturnal species so the random reordering isn't any different than the actual data on the tree
ggplot(phylo_signal_final, aes(x = clade, y = lambda, fill = log(p))) + geom_bar(stat = "identity") + geom_text(aes(label = round(lambda, digits = 3)), vjust = -0.2)

knitr::kable(phylo_signal_final, format = "html", digits = 3, caption = "Table X") %>%  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% save_kable("phylosig_table_longer.html")
webshot("phylosig_table_longer.html", file = "C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/phylogenetic_signal.pdf")





# Latitude maps -----------------------------------------------------------
lat.df <- read.csv(here("cetacean_latitude_df.csv"))

cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- cetaceans_full[, c("Parvorder", "Family", "Diel_Pattern", "max_crep", "Confidence", "tips")]

lat.df <- merge(cetaceans_full, lat.df, by = "tips", all = TRUE)
lat.df[lat.df == ""] <- NA

lat.df <- lat.df[!is.na(lat.df$max_lat),]
lat.df <- lat.df[!is.na(lat.df$max_crep),]

#what species have the smallest ranges? Largest ranges?
ggplot(lat.df, aes(x = (max_lat - min_lat), y = reorder(tips, (max_lat - min_lat)), fill = max_crep)) +
  geom_col() + facet_wrap(~Parvorder, scales = "free")

ggplot(lat.df, aes(y = (max_lat - min_lat), x = Family)) +
  geom_boxplot() 


range <- read_sf("C:/Users/ameli/Downloads/redlist_species_data_b8eeb8cf-3383-4314-bfb2-55dad2b8fec3/data_0.shp")

unique(range$SCI_NAME)

#filter for species with smallest ranges 
small_range_list <- lat.df %>% filter((max_lat - min_lat) < 30) %>% pull(Species_name)

sps_range <- range %>% filter(SCI_NAME %in% small_range_list)

ggplot(sps_range) +
  geom_sf(aes(fill = SCI_NAME), color = "black")


#filter for all the delphinidae ranges
delphinid_list <- lat.df %>% filter(Family == "Delphinidae") %>% pull(Species_name)
sps_range <- range %>% filter(SCI_NAME %in% delphinid_list)

ggplot(sps_range) +
  geom_sf(aes(fill = SCI_NAME), color = "black")


#filter for species by activity pattern
diel_list <- lat.df %>% filter(max_crep == "diurnal") %>% pull(Species_name)

sps_range <- range %>% filter(SCI_NAME %in% diel_list)

ggplot(sps_range) +
  geom_sf(aes(fill = SCI_NAME), color = "black")


#Do sympatric species show temporal niche partitioning?

lat.df %>% filter(Species_name %in% small_range_list)

##function to take max and min longitude

extractLongitude <-function(species_name){
  
  sps_range <- range %>% filter(SCI_NAME == species_name)
  
  xmin <- extent(sps_range)@xmin
  xmax <- extent(sps_range)@xmax
  
  latitude_list <- c(xmin, xmax)
  return(latitude_list)
}

longitude_list <- lapply(unique(range$SCI_NAME), function(x) extractLongitude(species_name = x))

names(longitude_list) <- unique(range$SCI_NAME)

lon.df <- data.frame(coords = unlist(longitude_list))
lon.df$Species_name <- rownames(lon.df)
lon.df <- lon.df %>% separate(Species_name, into = c("Species_name", "minmax"), sep = "\\.")

lon.df <- pivot_wider(lon.df, names_from = minmax, values_from = coords)
colnames(lon.df) <- c("Species_name", "min_lon", "max_lon")
lon.df$tips <- str_replace(lon.df$Species_name, pattern = " ", replacement = "_")
lon.df$mean_lon <- (lon.df$min_lon + lon.df$max_lon)/2

lon.df$tips <- str_replace(lon.df$Species_name, pattern = " ", replacement = "_")

trait.data <- merge(lat.df, lon.df, by = "tips", all = TRUE)

trait.data <- trait.data[!is.na(trait.data$max_crep),]

#install.packages("amt")
library(amt)







# Section X: Cetacean sleep duration --------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1F_m52NE1IWQRjfwgwF1yzVV4TaR4MEkKJUr2aDJt6LE/edit?usp=sharing'
duration <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

duration <- duration %>% filter(!is.na(Sleep_duration)) %>% select(tips, Sleep_duration)

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

trait.data <- merge(trait.data, duration, by = "tips", all = TRUE)

trait.data %>% filter(!is.na(Sleep_duration)) %>%
  ggplot(., aes(x = max_crep, y = Sleep_duration)) + 
  geom_boxplot() + geom_jitter()

trait.data %>% filter(!is.na(Sleep_duration)) %>%
  ggplot(., aes(x = log(Body_mass_kg), y = Sleep_duration)) + 
  geom_point() + geom_smooth(method = "lm")






# Section X: The testing ground -------------------------------------------

#Investigating the Bayestrait packages
library("devtools")
install_github("rgriff23/btw")
library("btw")




#Is there a hidden rate in the artiodactyla transition rates: are cetaceans transitioning faster than ruminants

trait.data <- read.csv(here("Sleepy_artiodactyla_full.csv"))
trait.data <- trait.data[!is.na(trait.data$max_crep), c("tips", "max_crep")]
phylo_trees <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#subset trait data to only include species that are in the tree
trait.data <- trait.data[trait.data$tips %in% phylo_trees$tip.label,]
# this selects a tree that is only the subset with data (mutual exclusive)
phylo_trees <- keep.tip(phylo_trees, tip = trait.data$tips)

hidden_rate_ARD <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 2, model = "ARD", node.states = "marginal")
plotMKmodel(hidden_rate_ARD)

ARD <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, model = "ARD", node.states = "marginal")
plotMKmodel(ARD)

#whippo bridge model
trait.data <- read.csv(here("whippomorpha.csv"))
#ruminant bridge model
trait.data <- read.csv(here("ruminants_full.csv"))
trait.data <- trait.data[!is.na(trait.data$max_crep), c("tips", "max_crep")]
phylo_trees <- readRDS(here("maxCladeCred_mammal_tree.rds"))

bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), model = "ARD", node.states = "marginal")
bridge_only_HR <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 2, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), model = "ARD", node.states = "marginal")
plotMKmodel(bridge_only)
plotMKmodel(brdige_only_HR)

#subset trait data to only include species that are in the tree
trait.data <- trait.data[trait.data$tips %in% phylo_trees$tip.label,]
# this selects a tree that is only the subset with data (mutual exclusive)
phylo_trees <- keep.tip(phylo_trees, tip = trait.data$tips)

bridge_only <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, rate.mat = matrix(c(0,1,2,3,4,0,5,6,7,8,0,0,10,11,0,0), ncol = 4, nrow = 4), model = "ARD", node.states = "marginal")

#simmap: plot lineages through time based on max clade cred tree: compare cetaceans and ruminants
trait.data <- read.csv(here("whippomorpha.csv"))
trait.data <- trait.data[!is.na(trait.data$max_crep), c("tips", "max_crep")]
phylo_trees <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#subset trait data to only include species that are in the tree
trait.data <- trait.data[trait.data$tips %in% phylo_trees$tip.label,]
# this selects a tree that is only the subset with data (mutual exclusive)
phylo_trees <- keep.tip(phylo_trees, tip = trait.data$tips)

model_ARD <- corHMM(phy = phylo_trees, data = trait.data, rate.cat = 1, model = "ARD", node.states = "marginal")
  
simmaps <- corHMM::makeSimmap(tree = phylo_trees, data = trait.data, rate.cat = 1, model = model_ARD$solution, nSim = 2, nCores = 1)

plotSimmap(simmaps)


# Concordance for each confidence level -----------------------------------

function to plot the concordance for each of the confidence levels
plotConcordance = function(set_column = "Conf2"){
  diel_full_filtered <- diel_full %>% filter(column == set_column)
  #need to filter for species with more than one entry or else concordance will always be 100%
  mulitple_sources <- diel_full_filtered %>% count(Species_name) %>% filter(n>1)
  diel_full_filtered <- diel_full_filtered[diel_full_filtered$Species_name %in% mulitple_sources$Species_name,]
  concordance <- as.data.frame(table(diel_full_filtered$max_crep, diel_full_filtered$value))
  colnames(concordance) <- c("actual", "predicted", "freq")
  totals_df <- aggregate(concordance$freq, by=list(Category=concordance$actual), FUN=sum)
  colnames(totals_df) <- c("actual", "total")
  concordance <- merge(concordance, totals_df, by = "actual")
  concordance$percent <- round(concordance$freq / concordance$total * 100, 1)
  return(concordance)
}

concordance_list <- lapply(sort(unique(diel_full$column)), function(x){plotConcordance(x)})

use below if
for(i in seq_along(concordance_list)){
  pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "artio", "confidence", i, "_confusion_matrix.pdf"))
  print(ggplot(as.data.frame(concordance_list[i]), aes(actual, predicted, fill = percent)) + geom_tile() + geom_text(aes(label = percent)) +
          scale_fill_gradient(low = "white", high = "dodgerblue") + labs(x = "Actual", y = "Predicted") +
          ggtitle(paste("Confidence level ", i, " concordance")))
  dev.off()
}





#Stochastic mapping with bayesian framework
#tutorial from Liam Revell, 2017 https://blog.phytools.org/2017/11/visualizing-rate-of-change-in-discrete.html
library(phytools)
tree

x






# Section: Transition rate scrap ------------------------------------------

#how do the principal components relate to the data 

ggplot(rates_df1, aes(y = log(rates), x = model_number, colour = solution)) + 
  theme_minimal() +
  geom_point()

#correlations between rates (I've done this before)
rates_wider <- rates_df1 %>% select(rates, solution, model_number) %>%
  pivot_wider(., names_from = solution, values_from = rates)

ggplot(rates_wider, aes(x = log(`Crepuscular -> Cathemeral`), y = log(`Diurnal -> Cathemeral`))) +
  geom_point() + geom_smooth(method = "lm")

#round values to one digit
ggplot(rates_df1, aes(x = log(round(rates)))) + 
  #theme_minimal() +
  geom_density() + facet_wrap(~solution)

#can take the mean because there are no zeros (only values very close to zero)
rates_df1 %>% group_by(solution) %>% summarize(mean_rates = mean(rates), SD_rates = sd(rates)) %>%
  ggplot(., aes(x = solution, y = mean_rates, fill = solution)) + geom_col() + 
  geom_errorbar(aes(ymin = mean_rates, ymax = mean_rates + SD_rates))




# Section: Does cetacean orbit size associate with eye size ---------------

trait.data <- read.csv(here("cetacean_ecomorphology_dataset.csv"))

trait.data %>% filter(!is.na(Orbit_ratio) & !is.na(Mean_dive_depth_m)) %>%
  ggplot(., aes(x = Orbit_ratio, y = Dive_depth_m)) +
  geom_point() + geom_smooth(method = "lm") + stat_poly_eq() +
  facet_wrap(~Family)


# Delphinidae dive depth --------------------------------------------------
trait.data.1 <- trait.data[!is.na(trait.data$Dive_depth_m),]
trait.data.delph <- trait.data.1 %>% filter(Family == "Delphinidae")

phylANOVA <- calculatePhylANOVA(trait.data.delph, "Dive_depth_m")

stat.test <- data.frame(group1 = c("cathemeral", "cathemeral", "cathemeral", "crepuscular", "crepuscular", "diurnal"),
                        group2 = c("crepuscular", "diurnal", "nocturnal", "diurnal", "nocturnal", "nocturnal"),
                        p.adj = c(phylANOVA$Pt[2], phylANOVA$Pt[3], phylANOVA$Pt[4], phylANOVA$Pt[7], phylANOVA$Pt[8], phylANOVA$Pt[12]),
                        y.position = c(7.2, 7.6, 8, 8.4, 8.8, 9.3))

stat.test <- stat.test %>% add_x_position(x = "max_crep")

delph_dive_boxplot <- ggplot(trait.data.delph, aes(x = max_crep, y = log(Dive_depth_m))) +
  geom_boxplot(aes(fill = max_crep), alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = custom.colours, guide = "none") +
  new_scale_fill() + 
  labs(x = "Temporal activity pattern", y = "Log (maximum dive depth (m))") + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = 'black', fill = "dodgerblue", pch = 21) +
  annotate("text", x = 1.3, y = 9, label = paste("phylANOVA, p =", phylANOVA$Pf)) +
  boxplot_theme +
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_x_discrete(labels = c("cathemeral" = "Cathemeral", "crepuscular" = "Crepuscular", "diurnal" = "Diurnal", "nocturnal" = "Nocturnal")) 

delph_dive_boxplot

# pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/Dive_depth_boxplot_delphinidae.pdf", width = 7, height = 7.5)
# delph_dive_boxplot
# dev.off()



# Section 6: Proportion plots OLD----------------------------

new_mammals <- read.csv(here("Bennie_mam_data.csv")) #data from Bennie et al, 2014, 4477 sps
new_mammals <- new_mammals[!is.na(new_mammals$max_crep), ] 

#add in my primary source data 
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_full <- artio_full[!is.na(artio_full$Diel_Pattern), ]

#we want to compare all mammals, vs all artiodactyla vs cetaceans/ruminants
new_mammals$mammals <- "Mammals"

custom.colours <- c("#dd8ae7","#EECBAD", "#FC8D62", "#66C2A5")
mammals_plot <- ggplot(new_mammals, aes(x = mammals, fill = max_crep)) + geom_bar(position = "fill", width = 0.75) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank())
ruminantia_plot <- artio_full %>% filter(Suborder == "Ruminantia") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())
whippomorpha_plot <- artio_full %>% filter(Suborder == "Whippomorpha") %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())

artiodactyla_plot <- artio_full %>% ggplot(., aes(x = Order, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())
#suborders_plot <- 

#plot the full order proportions and the suborders in the same barplot
artio_full %>% mutate(Suborder = "Artiodactyla") %>% rbind(., artio_full) %>%
  ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) +
  theme_bw() +
  #facet_wrap(~Suborder, nrow = 1, scales = "free") +
  labs(y = "Proportion", x = "Clade")
#theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())

#suborders_plot <- artio_full %>% filter(Suborder %in% c("Whippomorpha", "Ruminantia")) %>% ggplot(., aes(x = Suborder, fill = max_crep)) + geom_bar(position = "fill", width = 0.6) + scale_fill_manual(values = custom.colours) + theme_minimal() #+ theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/barplot_percentages.pdf", width = 5, height = 2, bg = "transparent")
(artiodactyla_plot + plot_spacer() + plot_spacer() + plot_spacer()) / 
  suborders_plot
dev.off()

mammals_plot <- 
  new_mammals %>% group_by(max_crep) %>% summarize(count = n()) %>%
  ggplot(., aes(x = "", y = count, fill = max_crep)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank()) + 
  scale_fill_manual(values = custom.colours) + geom_text(aes(label = round((count/4477)*100, digits = 1)), position = position_stack(vjust = 0.5))

ruminantia_plot <- 
  artio_full %>% filter(Suborder == "Ruminantia") %>% 
  group_by(max_crep) %>% summarize(count = n()) %>%
  ggplot(., aes(x = "", y = count, fill = max_crep)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank()) + 
  scale_fill_manual(values = custom.colours) + geom_text(aes(label = round((count/206)*100, digits = 1)), position = position_stack(vjust = 0.5))

whippomorpha_plot <- 
  artio_full %>% filter(Suborder == "Whippomorpha") %>% 
  group_by(max_crep) %>% summarize(count = n()) %>%
  ggplot(., aes(x = "", y = count, fill = max_crep)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") + 
  theme_void() + theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank()) + 
  scale_fill_manual(values = custom.colours) + geom_text(aes(label = round((count/84)*100, digits = 1)), position = position_stack(vjust = 0.5))

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/piechart_percentages.pdf", width = 10, height = 5, bg = "transparent")
grid.arrange(mammals_plot, ruminantia_plot, whippomorpha_plot, nrow = 1)
dev.off()


# Section : Max crep sankey with labels ----------------------------------
#add the concordance values to the sankey diagram
df <- mammals_df1 %>% make_long(Bennie_diel, Amelia_diel, Maor_diel,)

#with labels
concordance <- as.data.frame(table(mammals_df1$Bennie_diel, mammals_df1$Amelia_diel))
colnames(concordance) <- c("Bennie_data", "Amelia_data", "freq")
totals_df <- aggregate(concordance$freq, by=list(Category=concordance$Bennie_data), FUN=sum)
colnames(totals_df) <- c("Bennie_data", "total")
concordance <- merge(concordance, totals_df, by = "Bennie_data")
concordance$percent <- round(concordance$freq / concordance$total * 100, 1)
concordance$freq_count <- paste0(round(concordance$percent, 2), "%", " ", "(n=", concordance$freq, ")")
Bennie_freq_count <- concordance$freq_count
Bennie_freq_count <- Bennie_freq_count[-5] #remove zero for crep to cath

concordance <- as.data.frame(table(mammals_df1$Maor_diel, mammals_df1$Amelia_diel))
colnames(concordance) <- c("Maor_data", "Amelia_data", "freq")
totals_df <- aggregate(concordance$freq, by=list(Category=concordance$Maor_data), FUN=sum)
colnames(totals_df) <- c("Maor_data", "total")
concordance <- merge(concordance, totals_df, by = "Maor_data")
concordance$percent <- round(concordance$freq / concordance$total * 100, 1)
concordance$freq_count <- paste0(round(concordance$percent, 2), "%", " ", "(n=", concordance$freq, ")")
Maor_freq_count <- concordance$freq_count

sankey <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha= 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") + scale_fill_manual(values = c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5")) +
  theme_sankey(base_size = 16) + scale_x_discrete(labels = c("Bennie_diel" = "Existing database \n (Bennie et al)", "Amelia_diel" = "Current database \n (Mesich et al)", "Maor_diel" = "Existing database \n (Maor et al)")) +
  theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent')) + labs(y = NULL, x = NULL) + 
  annotate("text", x = 1.27, y = c(-108,-98,-88,-84,-60,-49,-45,-31.5,-4,38,50,64,80,93.5,109),
           label = Bennie_freq_count, size = 3.5, colour = "grey25") + 
  annotate("text", x = 2.73, y = c(-108,-90,-73,-67,-54,-30,-4,4,18,35,63,73,86,95,102,109),
           label = Maor_freq_count, size =3.5, colour = "grey25")
sankey

pdf(paste0("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "Maor_Bennie_sankey_labelled.pdf"), width = 10.5, height = 8)
test
dev.off()

