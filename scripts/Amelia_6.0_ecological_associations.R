source("scripts/fish_sleep_functions.R")

# Churchill et al, traits: orbit ratio ---------------------------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
#orbit information
church_pt1 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt1.xlsx")
#we are interested in orbit ratio, length of orbit relative to the bigozymatic width (proxy for body size)
#this is a continuous trait

trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern)), c("tips", "Diel_Pattern", "Parvorder")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular", "cathemeral/crepuscular"), ]

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

colnames(church_pt1) <- c("Species", "Family", "Specimen_number", "Left_orbit_length", "Right_orbit_length", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio")

church_pt1 <- church_pt1[!(is.na(church_pt1$Orbit_ratio)), c("Species", "Family", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio")]
church_pt1$tips <- str_replace(church_pt1$Species, " ", "_")

#filter the orbit size data by the species in our trait data, leaves 54 species
church_pt1 <- filter(church_pt1, tips %in% trait.data$tips)

#take the mean orbit size since there are some species with duplicates
church_pt1 <- church_pt1 %>% group_by(tips) %>% mutate(Orbit_ratio = mean(Orbit_ratio),  Bizygomatic_width = mean( Bizygomatic_width), Average_orbit_length = mean(Average_orbit_length))

#merge diel pattern data with orbit ratio data, leaves 54 species
trait.data <- merge(trait.data, church_pt1)

#52 species in the final tree,we lose M hotaula and S plumbea 
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#calculate the phyogenetic signal of orbit size
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
trait.vector <- test$Orbit_ratio
names(trait.vector) <- test$tips

orbit_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
#k = 1.97676 for cetaceans, p-value = 0.001 for cetaceans

orbit_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
lambda <- round(orbit_lambda$lambda, digits = 2)
#lambda = 1.01589 for cetaceans, p-value = 1.297 e-12 for cetaceans

#custom.colours <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Orbit_ratio")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Orbit_ratio), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = paste("Orbit Ratio", "\n", "λ = ", lambda))
diel.plot <- diel.plot #+ geom_tiplab(size = 2, offset = 3) 
diel.plot

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_vs_diel_cetaceans.pdf")
diel.plot
dev.off()

comparison_list <- list(c("cathemeral", "crepuscular"), c("crepuscular", "diurnal"), c("diurnal", "nocturnal"), c("cathemeral", "diurnal"), c("cathemeral", "nocturnal"),  c("crepuscular", "nocturnal"))

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
boxplot_KW <- ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(comparisons = comparison_list) + stat_compare_means(label.y = 60)

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_boxplots_KW_cetaceans.pdf")
boxplot_KW
dev.off()

boxplot_anova <- ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = 37) + stat_compare_means(label.y = 40, method = "anova")

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_boxplots_anova_cetaceans.pdf")
boxplot_anova
dev.off()

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
orbit_ANOVA <- aov(Orbit_ratio ~ Diel_Pattern, data = trait.data)
summary(orbit_ANOVA) #statistically significant for cetaceans, p value of F statistic = 0.014
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.00213

sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
#y trait is the continuous (response) variable 
trait.y <- test$Orbit_ratio
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

#phylogenetically corrected ANOVA
orbit_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant for cetaceans, p value of F statistic = 0.441
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.002

boxplot_phyloanova <- ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) +
  annotate("text", x = 1:4, y = 38, label = list("ns", "ns", "ns", "ns")) + annotate("text", x = 1.5, y = 41, label = paste("phylANOVA, p =", orbit_phylANOVA$Pf))

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_boxplots_phyloanova_cetaceans.pdf")
boxplot_phyloanova
dev.off()

# Baker et al artiodactyla axial length -----------------------------------
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
artio_eyes <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Baker_2019.xlsx")
artio_eyes <- artio_eyes[2: nrow(artio_eyes),]
colnames(artio_eyes) <- c("tips", "Order", "Corneal_diameter", "Axial_length", "Activity_pattern", "Source")
artio_eyes <- filter(artio_eyes, Order == "Artiodactyla")
artio_eyes <- artio_eyes[artio_eyes$Corneal_diameter != "n/a", c("tips", "Corneal_diameter", "Axial_length", "Activity_pattern") ]
sleepy_artio <- read.csv(here("sleepy_artiodactyla_full.csv"))
artio_eyes <- artio_eyes %>% filter(artio_eyes$tips %in% sleepy_artio$tips)
artio_eyes$Transverse_diameter <- "NA"
artio_eyes <- artio_eyes %>% relocate(Transverse_diameter, .before= Activity_pattern)

trait.data <- merge(sleepy_artio, artio_eyes, by = "tips")
trait.data$Axial_length <- as.numeric(trait.data$Axial_length)
trait.data$Corneal_diameter <- as.numeric(trait.data$Corneal_diameter)
trait.data <- trait.data[, c("tips", "Species_name", "Family", "Order", "Diel_Pattern", "Corneal_diameter", "Axial_length", "Activity_pattern")]

#we can look at the ratio of corneal diameter to axial length
#From Hall et al: https://doi.org/10.1098/rspb.2012.2258
#The ratio of corneal diameter to axial length of the eye is a useful measure of relative sensitivity 
#and relative visual acuity that has been used in previous studies as a way to compare animals of disparate size
trait.data$Orbit_ratio <- trait.data$Corneal_diameter/trait.data$Axial_length

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
#two artio species dropped from the tree Camelus bactrianus and Lama glama
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#calculate the phyogenetic signal of orbit size
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
trait.vector <- test$Orbit_ratio
names(trait.vector) <- test$tips

orbit_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
#k = 1.97676 for cetaceans, p-value = 0.001 for cetaceans

orbit_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
lambda <- round(orbit_lambda$lambda, digits = 2)
#lambda = 1.01589 for cetaceans, p-value = 1.297 e-12 for cetaceans

#custom.colours <- c("#dd8ae7", "peachpuff2", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Orbit_ratio")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +1.5, y=y, fill = Diel_Pattern), width = 3, inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +6, y=y, fill = Orbit_ratio), width = 3, inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = paste("Orbit Ratio", "\n", "λ = ", lambda))
diel.plot <- diel.plot #+ geom_tiplab(size = 2, offset = 3) 
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/artio_orbit_ratio_vs_diel_.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

comparison_list <- list(c("crepuscular", "diurnal"), c("crepuscular", "nocturnal"), c("diurnal", "nocturnal"))

#plot out  orbit ratio vs activity pattern (NOT phylogenetically collected)
boxplot_KW <- ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(comparisons = comparison_list) + stat_compare_means(label.y = 1.01)

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_boxplots_KW_artio.pdf")
boxplot_KW
dev.off()

boxplot_anova <- ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = 1) + stat_compare_means(label.y = 1.01, method = "anova")

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_boxplots_anova_artio.pdf")
boxplot_anova
dev.off()

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
orbit_ANOVA <- aov(Orbit_ratio ~ Diel_Pattern, data = trait.data)
summary(orbit_ANOVA) #statistically significant for cetaceans, p value of F statistic = 0.014
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.00213

sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
#y trait is the continuous (response) variable 
trait.y <- test$Orbit_ratio
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

#phylogenetically corrected ANOVA
orbit_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant for cetaceans, p value of F statistic = 0.441
#statistically significant for non-cetacean artiodactyls, p value of F statistic = 0.002

boxplot_phyloanova <- ggplot(trait.data, aes(x = Diel_Pattern, y = Orbit_ratio)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Orbit ratio") + scale_fill_manual(values=as.vector(polychrome(26))) +
  annotate("text", x = 1.5, y = 1, label = paste("phylANOVA, p =", orbit_phylANOVA$Pf))

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/orbit_ratio_boxplots_phyloanova_artio.pdf")
boxplot_phyloanova
dev.off()


# # Churchill et al: dive depth -------------------------------------------

#read in trait data, dive depth on 56 species
trait.data <- read.csv(here("cetacean_dive_depth.csv"))
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#use below to look just at odontocetes or mysticetes
#trait.data <- trait.data %>% filter(Parvorder == "Odontoceti")
trait.data <- trait.data %>% filter(Parvorder == "Mysticeti")

#53 species with dive data are in the tree, B brydei, B ricei and S plumbea are not 
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

trait.vector <- trait.data$Dive_depth
names(trait.vector) <- trait.data$tips

dive_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
#K: 0.714, p-value = 0.001

dive_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list(), niter=10)
lambda <- round(dive_lambda$lambda, digits = 2)
#lambda: 1.0027, p-value = 1.287 e-8

#custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5", "grey")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Dive_depth")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Dive_depth), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = paste("Maximum dive depth", "\n", "λ = ", lambda))
diel.plot <- diel.plot #+ geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/improved_dive_depth_vs_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#requires the trait data to be a named vector in the same species order as appears in trpy_n$tip.labels
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]
#y trait is the continuous (response) variable 
trait.y <- test$Dive_depth
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

comparison_list <- list(c("cathemeral", "crepuscular"), c("crepuscular", "diurnal"), c("diurnal", "nocturnal"), c("cathemeral", "diurnal"), c("cathemeral", "nocturnal"),  c("crepuscular", "nocturnal"))

#plot out  dive depth vs activity pattern (NOT phylogenetically collected)
boxplot_KW <- ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(comparisons = comparison_list) + stat_compare_means(label.y = 4000)

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/dive_depth_boxplots_KW.pdf")
boxplot_KW
dev.off()

boxplot_anova <- ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) +
  facet_wrap(~Parvorder) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 3500, method = "anova")

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/dive_depth_boxplots_anova.pdf")
boxplot_anova
dev.off()

#check for association between response variable y (dive depth) and categorical variable x (diel pattern)
#regular one way ANOVA
dive_ANOVA <- aov(Dive_depth ~ Diel_Pattern, data = trait.data)
summary(dive_ANOVA) #statistically significant, p value of F statistic = 0.0003
#statistically significant, p value of F statistic = 0.02

#pairwise t test
tapply(X = trait.data$Dive_depth, INDEX = list(trait.data$Diel_Pattern), FUN = mean)
pairwise.t.test(trait.data$Dive_depth, trait.data$Diel_Pattern,  p.adj = "none")

#phylogenetically corrected ANOVA
dive_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#statistically significant, p value of F statistic = 0.01
#not statistically significant, p value of F statistic = 0.227

boxplot_phyloanova <- ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) +
  annotate("text", x = 1.5, y = 700, label = paste("phylANOVA, p =", dive_phylANOVA$Pf))


# Section: Habitat vs activity pattern ------------------------------------

#we should see a similar pattern when we look at habitats (shallow habitats like rivers vs deep habitats like pelagic species)
trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern)), c("tips", "Diel_Pattern", "Family", "Parvorder")]

#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"),]
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace_all(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

echo <- read_xlsx("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_discrete_traits/Coombs_et_al_2021.xlsx")
echo <- echo[!(is.na(echo$Echo)), c("Museum ID", "Age", "Echo", "Diet", "Dentition", "FM", "Habitat")]
echo$tips <- str_replace(echo$`Museum ID`, " ", "_")
echo <- echo %>% filter(echo$tips %in% trait.data$tips)

trait.data <- merge(trait.data, echo, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
custom.colours.2 <- c( "grey90", "grey40", "grey66", "black")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Habitat")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() +  geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Habitat), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours.2, name = "Habitat")
diel.plot <- diel.plot #+ geom_tiplab(size = 2, offset = 3)
diel.plot

test_df <- trait.data[trait.data$Diel_Pattern %in% c("cathemeral", "diurnal","nocturnal", "crepuscular"), c("Diel_Pattern", "Habitat")]
test <- fisher.test(table(test_df))

mosaicplot(table(test_df), color = TRUE, main = paste0("Fischer's exact test: ", "p-value = ", round(test$p.value, 3))) 

dive_data <- read.csv(here("cetacean_dive_depth.csv"))
trait.data <- merge(trait.data, dive_data, by = "tips")

custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Habitat", "Dive_depth")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Habitat), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Habitat")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Dive_depth), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = paste("Maximum dive depth", "\n", "λ = ", lambda))
diel.plot <- diel.plot #+ geom_tiplab(size = 2, offset = 3)
diel.plot


boxplot_KW <- ggplot(trait.data, aes(x = Habitat, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(fill = Family.x), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Habitat", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(comparisons = comparison_list) + stat_compare_means(label.y = 4000)

pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/habitat_dive_depth_boxplots_KW.pdf")
boxplot_KW
dev.off()

boxplot_anova <- ggplot(trait.data, aes(x = Habitat, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family.x), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Habitat", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 3500, method = "anova")

ggplot(trait.data, aes(x = Habitat, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Diel_Pattern.x), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Habitat", y = "Dive depth") + scale_fill_manual(values=as.vector(polychrome(26))) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.") + stat_compare_means(label.y = 3500, method = "anova")


pdf("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/habitat_dive_depth_boxplots_anova.pdf")
boxplot_anova
dev.off()


# #Body mass from Churchill et all ----------------------------------------
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")

#use below to rerun with max-crep four state model
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

#filter for species with body mass data, leaves 62 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Mass)), c("Species", "Dive_depth", "Dive_duration", "Mass", "Total_body_length", "Assymetry_index", "Peak_frequency")]

#correct species spelling to match
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus commersoni", replacement = "Cephalorhynchus commersonii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Berardius bairdi", replacement = "Berardius bairdii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias australis", replacement = "Lagenorhynchus australis")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias cruciger", replacement = "Lagenorhynchus cruciger")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias obliquidens", replacement = "Lagenorhynchus obliquidens")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus heavisidi", replacement = "Cephalorhynchus heavisidii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Leucopleurus acutus", replacement = "Lagenorhynchus_acutus")

church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")

#filter the body mass by the species in our trait data 
#leaves 56 species, removes phocoena dioptrica, lagenorhynchus cruciger and some mesoplodon sps
church_pt2 <- church_pt2 %>% filter(church_pt2$tips %in% trait.data$tips)

trait.data <- merge(trait.data, church_pt2, by = "tips")
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]

#55 species with body mass are in the tree, Sousa_plumbea is not in the tree
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

#custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern", "Mass")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = Mass), inherit.aes = FALSE, colour = "transparent") + scale_fill_gradient(low = "#ceecff", high = "#07507e", name = "Body mass")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/body_size_vs_diel_.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

#plot out body size vs activity pattern (NOT phylogenetically corrected)
ggplot(trait.data, aes(x = Diel_Pattern, y = log(Mass))) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Parvorder), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  labs(x = "Temporal activity pattern", y = "Log body mass") + scale_fill_manual(values=as.vector(polychrome(26))) + facet_wrap(~Parvorder)

#identify the phylogenetic signal of dive depth
trait.vector <- trait.data$Mass
names(trait.vector) <- trait.data$tips

mass_K <- phylosig(trpy_n, trait.vector, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                   control=list(), niter=10)
#K: 1.96692
#p-value = 0.001

mass_lambda <- phylosig(trpy_n, trait.vector, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                        control=list(), niter=10)
#lambda: 1.01668
#p-value = 3.62866e-18 

#we can make a phenogram of body mass evolution
#requires the trait data to be a named vector in the same species order as appears in trpy_n$tip.labels
sps_order <- as.data.frame(trpy_n$tip.label)
colnames(sps_order) <- "tips"
sps_order$id <- 1:nrow(sps_order)
test <- merge(trait.data, sps_order, by = "tips")
test <- test[order(test$id), ]

#filter to examine mysticetes and odontocetes separately
#test <- filter(test, Parvorder == "Mysticeti")  #only 8 species
test <- filter(test, Parvorder == "Odontoceti")
trpy_n <- keep.tip(mam.tree, tip = test$tips)

#y trait is the continuous (response) variable 
trait.y <- test$Mass
names(trait.y) <- test$tips
#x trait is the categorical variable 
trait.x <- test$Diel_Pattern
names(trait.x) <- test$tips

#use function phenogram from phytools
phenogram(trpy_n, log(trait.y), ftype="reg", spread.labels = TRUE, spread.cost=c(1,0))

#colour the tips by the activity pattern
tiplabels(pie = to.matrix(trait.x, unique(trait.data$Diel_Pattern)), piecol = c("#dd8ae7", "#FC8D62", "#66C2A5","peachpuff2"), cex = 0.25)

#check for association between response variable y (orbit size) and categorical variable x (diel pattern)
#regular one way ANOVA
mass_ANOVA <- aov(Mass ~ Diel_Pattern, data = test)
summary(mass_ANOVA) #not statistically significant, p value of F statistic = 0.142

#phylogenetically corrected ANOVA
mass_phylANOVA <- phylANOVA(trpy_n, trait.x, trait.y, nsim=1000, posthoc=TRUE, p.adj="holm")
#not statistically significant, p value of F statistic = 0.393


# Section: Organize dive data ---------------------------------------------

#read in the Churchill et al data
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")
colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

#filter for species with dive data, leaves 62 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Dive_depth)), c("Species", "Dive_depth", "Dive_duration")]

#correct species spelling to match
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus commersoni", replacement = "Cephalorhynchus commersonii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Berardius bairdi", replacement = "Berardius bairdii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias australis", replacement = "Lagenorhynchus australis")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias cruciger", replacement = "Lagenorhynchus cruciger")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Sagmatias obliquidens", replacement = "Lagenorhynchus obliquidens")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Cephalorhynchus heavisidi", replacement = "Cephalorhynchus heavisidii")
church_pt2$Species <- str_replace(church_pt2$Species, pattern = "Leucopleurus acutus", replacement = "Lagenorhynchus_acutus")

church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")

#read in the additional dive data I collected
url <- 'https://docs.google.com/spreadsheets/d/1_0ZS_tbddOCckkcKn9H5HpVRDZty4jhkUU20Nc0YYQY/edit?usp=sharing'
dive_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#average dive depth is very inconsistent (deep vs shallow dives, day vs night dives)
#and measured inconsistently (average of means across individuals, average of medians etc) so exclude for now
dive_full <- dive_full[!is.na(dive_full$Diel_Pattern),1:9]

#filter to rows with an alternative value
dive_full_alt <- dive_full[!is.na(dive_full$Alt_Max_1), ]
dive_full <- dive_full[is.na(dive_full$Alt_Max_1), ]

#pick the bigger value between max dive depth and alt dive depth
dive_full_alt[is.na(dive_full_alt)] <- 0
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_1)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_2)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_3)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_4)

#combine with original values
dive_full <- rbind(dive_full, dive_full_alt)
dive_full <- dive_full[, c("Species_name", "Max_dive_depth_m")]
colnames(dive_full) <- c("Species", "Max_dive_depth_m")
church_pt2 <- merge(church_pt2, dive_full, all.x = TRUE, all.y = TRUE)
church_pt2[is.na(church_pt2)] <- 0
church_pt2$Final_dive_depth <- pmax(church_pt2$Max_dive_depth_m, church_pt2$Dive_depth)
church_pt2$tips <- str_replace(church_pt2$Species, pattern = " ", replacement = "_")

#more data
Laeta_data <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Laeta_2020_dive_data.xlsx")
Laeta_data$tips <- str_replace(Laeta_data$Species, " ", "_")
Laeta_data <- Laeta_data[!is.na(Laeta_data$Species), c("Depth (m)", "tips")]
colnames(Laeta_data) <- c("Dive_depth_Laeta", "tips")
Laeta_data <- filter(Laeta_data, Dive_depth_Laeta != "-" )
Laeta_data$Dive_depth_Laeta <- as.integer(Laeta_data$Dive_depth_Laeta)

Dive_depth <- merge(church_pt2, Laeta_data, all.x = TRUE, all.y = TRUE)
Dive_depth[is.na(Dive_depth)] <- 0
Dive_depth$Dive_depth_final <- pmax(Dive_depth$Final_dive_depth, Dive_depth$Dive_depth_Laeta)

#save out only the relevant columns, leaves 64 species
Dive_depth <- Dive_depth %>% filter(Dive_depth_final != 0)
Dive_depth <- Dive_depth[, c(1,2,8)]
colnames(Dive_depth) <- c("tips", "Species", 'Dive_depth')
#four entries only have genus (Delphinus_sp, Lagenorhynchus_cruciger, Lissodelphis_sp, Neophocaena_sp,)

trait.data <- read.csv(here("cetaceans_full.csv"))
trait.data <- trait.data[!(is.na(trait.data$Diel_Pattern)), c("tips", "Diel_Pattern", "Parvorder", "Family")]
#filer out species with an unknown activity pattern
trait.data <- trait.data[trait.data$Diel_Pattern %in% c("diurnal", "cathemeral", "diurnal/crepuscular", "nocturnal", "nocturnal/crepuscular"), ]

#use below to rerun with max-crep four state model
# trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "diurnal/crepuscular", replacement = "crepuscular")
# trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "nocturnal/crepuscular", replacement = "crepuscular")
# trait.data$Diel_Pattern <- str_replace(trait.data$Diel_Pattern, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#final dive data has 56 species
trait.data <- merge(trait.data, Dive_depth, by = "tips")

write.csv(trait.data, here("cetacean_dive_depth.csv"), row.names = FALSE)
