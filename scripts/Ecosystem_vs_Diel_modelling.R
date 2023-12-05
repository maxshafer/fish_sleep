library(rfishbase)
library(ggalt)
library(ggbiplot)
library(dplyr)
library(patchwork)
library(ggridges)
library(here)

setwd(here::here())

source(here::here("scripts/fish_sleep_functions.R"))

## Load files

resolved_names <- read.csv("resolved_names_local.csv", row.names = "X", header = TRUE)
# tr.calibrated <- readRDS("calibrated_phylo.rds")
trait.data <- readRDS("trait_data_fish.rds")
# trpy <- readRDS(file = "calibrated_phylo.rds")

#############  #############  #############  #############  #############  #############  #############  
#############  #############  #############  #############  #############  #############  #############  

## Load fishbase data

# fishbase_df <- rfishbase::load_taxa()
# fishbase_diet <- diet()
# fishbase_morph <- morphometrics()
# fishbase_maturity <- maturity()
# fishbase_spawning <- spawning()
# fishbase_swimming <- swimming()
# fishbase_larvae <- larvae()
# fishbase_dietitems <- diet_items()
# fishbase_fecundity <- fecundity()

fishbase_reproduction <- reproduction(version = "21.06")
fishbase_ecology <- ecology(version = "21.06") # Has troph and 
fishbase_length <- length_weight(version = "21.06")
fishbase_ecosystem <- ecosystem(version = "21.06") # Salinity is here duh
fishbase_species <- species(version = "21.06")

#############  #############  #############  #############  #############  #############  #############  
#############  #############  #############  #############  #############  #############  #############  

## OK need to collect my traits, then transform them, then do PCA

############# TL/SL for 11k species or actual length measurements ############# 
# standard_size <- fishbase_morph %>% group_by(Species) %>% dplyr::summarise(mean_TL = mean(TL), mean_SL = mean(SL))
# fishbase_length2 <- fishbase_length[complete.cases(fishbase_length[,c("LengthMax", "a", "b")])]
standard_size <- fishbase_length %>% group_by(Species) %>% dplyr::summarise(mean_LengthMax = mean(LengthMax), mean_Weight = mean(a*LengthMax^b))


# #############  Next is reproduction, rate, size, age ############# 
# 
# # Max fecundity for 456 species (sad), SpawningCycles for 589, gestation values for 30-40 (boo)
# fecundity <- fishbase_spawning %>% group_by(Species) %>% dplyr::summarise(mean_FecundityMax = mean(FecundityMax), mean_SpawningCycles = mean(SpawningCycles))
# # Can also get 'spawning ground' from fishbase_spawning, which is a character vector for ~1500 species
# spawning_ground <- fishbase_spawning %>% group_by(Species) %>% dplyr::summarise(spawning_ground = paste(unique(SpawningGround), sep = "_"))
# # This gives me the age at maturity for 1015 species
# age_maturity <- fishbase_maturity %>% group_by(Species) %>% dplyr::summarize(mean_AgeMatMin = mean(c(AgeMatMin, AgeMatMin2), na.rm = TRUE))


############# Next is trophic level, benthopelagic, and ecosystem breadth) #############  

trophic <- fishbase_ecology %>% group_by(Species) %>% dplyr::summarise(mean_DietTroph = mean(DietTroph), mean_FoodTroph = mean(FoodTroph), FeedingType = paste(unique(FeedingType), sep = "_"))
ecosystem <- fishbase_ecosystem %>% group_by(Species) %>% dplyr::summarise(ecosystem_type = paste(unique(EcosystemType), sep = "_"))
realm <- fishbase_ecosystem[fishbase_ecosystem$EcosystemType == "Zoogeographic realm",] %>% group_by(Species) %>% dplyr::summarise(realm = paste(unique(EcosystemName), sep = "_"))
benthopelagic <- fishbase_species %>% group_by(Species) %>% dplyr::summarise(benthopelagic = paste(unique(DemersPelag), sep = "_"), Fresh_Brack_Saltwater = paste(mean(Fresh), mean(Brack), mean(Saltwater), sep = "_"))


#############  Reproduction! Can use the Guild, which is a reflection of parental investment, and is super complete ############# 

reproduction <- fishbase_reproduction %>% group_by(Species) %>% dplyr::summarise(RepGuild1 = paste(unique(RepGuild1), sep = "_"), RepGuild2 = paste(unique(RepGuild2), sep = "_"), ParentalCare = paste(unique(ParentalCare), sep = "_"))

#############  #############  #############  #############  #############  #############  #############  
#############  #############  #############  #############  #############  #############  #############  

## Make a combined list
combined_ecology_metrics <- Reduce(merge, list(standard_size, trophic, ecosystem, realm, benthopelagic, reproduction))

## Find those species with 'complete.cases'
combined_ecology_metrics <- data.frame(lapply(combined_ecology_metrics, function(x) gsub("NA", NA, x)))
combined_ecology_metrics <- combined_ecology_metrics[!(is.na(combined_ecology_metrics$Species)),]

zoo <- combined_ecology_metrics[combined_ecology_metrics$ecosystem_type == "Zoogeographic realm",]
combined_ecology_metrics <- combined_ecology_metrics[combined_ecology_metrics$ecosystem_type != "Zoogeographic realm",]

#############  Calculate the ecosystem breadth (the number of ecosystem types x the number of zoogeographic realms)#############  

## This isn't working correctly atm

realm_vector <- lapply(unique(zoo$Species), function(x) {
  vector <- unique(tolower(zoo[zoo$Species %in% x, "realm"]))
  vector2 <- paste(unique(sort(vector)), collapse = "_")
  df <- data.frame(Species = x, realm_names = vector2, realm_count = length(vector))
  return(df)
})

realm_df <- Reduce(rbind, realm_vector)

ecosystem_vector <- lapply(unique(combined_ecology_metrics$Species), function(x) {
  vector <- unique(tolower(combined_ecology_metrics[combined_ecology_metrics$Species %in% x, "ecosystem_type"]))
  vector2 <- paste(unique(sort(vector)), collapse = "_")
  df <- data.frame(Species = x, ecosystem_names = vector2, ecosystem_count = length(vector))
  return(df)
})

ecosystem_df <- Reduce(rbind, ecosystem_vector)

#############  Combine back with all data ############# 

combined_ecology_metrics <- combined_ecology_metrics[,c("Species", "mean_LengthMax", "mean_Weight", "mean_DietTroph", "mean_FoodTroph", "benthopelagic", "RepGuild1", "RepGuild2")]
combined_ecology_metrics <- combined_ecology_metrics[!duplicated(combined_ecology_metrics),]
combined_ecology_metrics <- Reduce(merge, list(combined_ecology_metrics, ecosystem_df, realm_df))


#############  Calculate simplified Trophic and reproduction guild metrics #############  
combined_ecology_metrics$RepGuild <- ifelse(!(is.na(combined_ecology_metrics$RepGuild2)), combined_ecology_metrics$RepGuild2, combined_ecology_metrics$RepGuild1)
combined_ecology_metrics$mean_DietTroph <- as.numeric(combined_ecology_metrics$mean_DietTroph)
combined_ecology_metrics$mean_FoodTroph <- as.numeric(combined_ecology_metrics$mean_FoodTroph)
combined_ecology_metrics$Trophic <- ifelse(is.na(combined_ecology_metrics$mean_FoodTroph), combined_ecology_metrics$mean_DietTroph, ifelse(!(is.na(combined_ecology_metrics$mean_DietTroph)) & !(is.na(combined_ecology_metrics$mean_FoodTroph)), rowMeans(combined_ecology_metrics[,c("mean_DietTroph", "mean_FoodTroph")]), combined_ecology_metrics$mean_FoodTroph))


#############  Subset for only those with diel information #############
# combined_ecology_metrics_diel <- combined_ecology_metrics[combined_ecology_metrics$Species %in% gsub("_", " ", trait.data$species),]

combined_ecology_metrics_diel <- combined_ecology_metrics
combined_ecology_metrics_diel$diel <- trait.data$diel2[match(combined_ecology_metrics_diel$Species, gsub("_", " ", trait.data$species))]

combined_ecology_metrics_diel <- combined_ecology_metrics_diel[!(is.na(combined_ecology_metrics_diel$diel)),]

combined_ecology_metrics_diel$mean_Weight <- log(as.numeric(combined_ecology_metrics_diel$mean_Weight))
combined_ecology_metrics_diel$benthopelagic <- factor(combined_ecology_metrics_diel$benthopelagic, levels = c("demersal", "reef-associated", "benthopelagic", "pelagic-neritic", "pelagic", "pelagic-oceanic"))
combined_ecology_metrics_diel$RepGuild <- tolower(combined_ecology_metrics_diel$RepGuild)
combined_ecology_metrics_diel$RepGuild <- factor(combined_ecology_metrics_diel$RepGuild, levels = c("open water/substratum egg scatterers", "nonguarders", "brood hiders", "guarders", "clutch tenders", "nesters", "external brooders", "internal live bearers"))


# Make histograms

combined_ecology_metrics_diel$diel <- factor(combined_ecology_metrics_diel$diel, levels = c("nocturnal","diurnal","crepuscular","unclear"))

p_mean_LengthMax <- ggplot(combined_ecology_metrics_diel, aes(x = log(as.numeric(mean_LengthMax)), y = diel, colour = diel, fill = diel)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) + theme_classic() + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H")
p_mean_Weight <- ggplot(combined_ecology_metrics_diel, aes(x = mean_Weight, y = diel, colour = diel, fill = diel)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H") + theme_classic()
p_mean_FoodTroph <- ggplot(combined_ecology_metrics_diel, aes(x = as.numeric(Trophic), y = diel, colour = diel, fill = diel)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H") + theme_classic()
p_benthopelagic <- ggplot(combined_ecology_metrics_diel, aes(x = as.numeric(benthopelagic), y = diel, colour = diel, fill = diel)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H") + theme_classic()
p_ecosystem_count <- ggplot(combined_ecology_metrics_diel, aes(x = as.numeric(ecosystem_count), y = diel, colour = diel, fill = diel)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H") + theme_classic()
p_RepGuild <- ggplot(combined_ecology_metrics_diel, aes(x = as.numeric(RepGuild), y = diel, colour = diel, fill = diel)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H") + theme_classic()

## This shows that the distribution among the individual ecological traits is similar across diel niches
## Individually they have more species info

histograms <- ((p_mean_LengthMax + p_mean_Weight + p_mean_FoodTroph + p_benthopelagic + p_RepGuild + p_ecosystem_count) & theme(axis.title.y = element_blank())) + plot_layout(ncol = 2, guides = "collect")




#############  Subset to only complete cases #############  
complete_cases <- tibble(combined_ecology_metrics[,c("Species", "mean_LengthMax", "mean_Weight", "Trophic", "benthopelagic", "RepGuild", "ecosystem_count", "realm_count")])
complete_cases <- complete_cases[complete.cases(complete_cases),]
# complete_cases <- combined_ecology_metrics_diel[complete.cases(combined_ecology_metrics_diel[,c("Species","Trophic","benthopelagic")]), c("Species", "mean_LengthMax", "mean_Weight", "Trophic", "benthopelagic", "RepGuild", "ecosystem_count")]

complete_cases$diel <- trait.data$diel2[match(complete_cases$Species, gsub("_", " ", trait.data$species))]
complete_cases$mean_Weight <- log(as.numeric(complete_cases$mean_Weight))
# View(lapply(complete_cases, function(x) table(is.na(x))))

## Assign factor levels to character vectors

# benthopelagic, RepGuild (and combine), FeedingType
complete_cases$benthopelagic <- factor(complete_cases$benthopelagic, levels = c("demersal", "reef-associated", "benthopelagic", "pelagic-neritic", "pelagic", "pelagic-oceanic"))
complete_cases$RepGuild <- tolower(complete_cases$RepGuild)
complete_cases$RepGuild <- factor(complete_cases$RepGuild, levels = c("open water/substratum egg scatterers", "nonguarders", "brood hiders", "guarders", "clutch tenders", "nesters", "external brooders", "internal live bearers"))

# complete_cases$FeedingType <- factor(complete_cases$FeedingType, levels = c("variable", "grazing on aquatic plants", "browsing on substrate", "filtering plankton", "other", "selective plankton feeding", "hunting macrofauna (predator)", "feeding on a host (parasite)"))



# Make histograms

p_mean_LengthMax <- ggplot(complete_cases, aes(x = log(as.numeric(mean_LengthMax)))) + geom_histogram(color = "black", fill = "grey", binwidth = 1) + theme_classic()
p_mean_Weight <- ggplot(complete_cases, aes(x = as.numeric(mean_Weight))) + geom_histogram(color = "black", fill = "grey", binwidth = 1) + theme_classic()
p_mean_FoodTroph <- ggplot(complete_cases, aes(x = as.numeric(Trophic))) + geom_histogram(color = "red4", fill = "red1", binwidth = 0.25) + theme_classic()
# p_FeedingType <- ggplot(complete_cases, aes(x = as.numeric(FeedingType))) + geom_histogram(color = "gold4", fill = "gold1") + theme_classic()
p_benthopelagic <- ggplot(complete_cases, aes(x = as.numeric(benthopelagic))) + geom_histogram(color = "royalblue4", fill = "royalblue1", binwidth = 1) + theme_classic()
p_ecosystem_count <- ggplot(complete_cases, aes(x = as.numeric(ecosystem_count))) + geom_histogram(color = "seagreen", fill = "seagreen2", binwidth = 1) + theme_classic()
p_RepGuild <- ggplot(complete_cases, aes(x = as.numeric(RepGuild))) + geom_histogram(color = "mediumpurple4", fill = "mediumpurple1", binwidth = 1) + theme_classic()

histograms_complete <- p_mean_LengthMax + p_mean_Weight + p_mean_FoodTroph + p_benthopelagic + p_RepGuild + p_ecosystem_count + plot_layout(ncol = 1)

complete_cases_2 <- complete_cases[complete_cases$diel %in% c("crepuscular", "diurnal", "nocturnal", "unclear"),]
# can include "realm_count" here as well if wanted
pca.data <- as.data.frame(complete_cases_2[,c("mean_Weight", "Trophic", "benthopelagic", "RepGuild", "ecosystem_count")])
row.names(pca.data) <- complete_cases_2$Species
for (i in 1:ncol(pca.data)) {
  pca.data[,i] <- as.numeric(pca.data[,i])
}
diel.pca <- prcomp(pca.data, center = TRUE, scale. = TRUE)


complete_cases_2$diel <- factor(complete_cases_2$diel, levels = c("nocturnal","diurnal","crepuscular","unclear"))
biplot <- ggbiplot(diel.pca, groups = complete_cases_2$diel, choices = c(1,2)) + theme_classic() 
# Find the convex hull of the points plotted
hull <- biplot$data %>% group_by(groups) %>% slice(chull(xvar, yvar))

# biplot <- biplot + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35")) + geom_polygon(data = hull, aes(fill = groups, colour = groups), alpha = 0, show.legend = FALSE) + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35")) + theme(legend.position = c(0.15,0.15))

biplot2 <- biplot + geom_polygon(data = hull, aes(fill = groups, colour = groups), alpha = 0.05, show.legend = FALSE) + scale_color_viridis_d(option = "H") #+ theme(legend.position = c(0.15,0.15))


pdf(file = "outs/Figures/Ecological_economics_histograms_all.pdf", width = 8, height = 8)
histograms
dev.off()

pdf(file = "outs/Figures/Ecological_economics_biplot.pdf", width = 4, height = 4)
biplot2
dev.off()


pca.data <- as.data.frame(diel.pca$x)


nesters <- c("Cottus gobio", "Heteropneustes fossilis", "Butis koilomatodon")
scatterers <- c("Vimba vimba", "Oligoplites palometa", "Pangasianodon hypophthalmus")
big_nesters <- c("Protopterus aethiopicus", "Amia calva", "Wallago attu", "Centrarchus macropterus")
loweco_tenders <- c("Andinoacara pulcher", "Synodontis nigriventris", "Pterophyllum scalare")

combined_ecology_metrics_diel[combined_ecology_metrics_diel$Species %in% c(scatterers, big_nesters, loweco_tenders),]

resolved_names[resolved_names$unique_name %in% c(scatterers, big_nesters, loweco_tenders),]

biplot$data$sof <- ifelse(rownames(biplot$data) %in% c(scatterers, big_nesters, loweco_tenders), rownames(biplot$data), NA)

pdf(file = "outs/Figures/Ecological_economics_summary_figure_labels.pdf", width = 15, height = 15)
biplot + geom_label_repel(data = biplot$data, aes(x = xvar, y = yvar, label = sof))
dev.off()

biplot + geom_label_repel(data= biplot$data, aes(x = xvar, y = yvar, label = rownames(biplot$data))) + xlim(0.5,2) + ylim(0.5,2)

pca.data$diel <- complete_cases$diel

ggplot(pca.data, aes(x = PC1, group = diel, color = diel, fill = diel)) + geom_density(alpha = 0.25) + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))
ggplot(pca.data, aes(x = PC2, group = diel, color = diel, fill = diel)) + geom_density(alpha = 0.25) + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))
ggplot(pca.data, aes(x = PC3, group = diel, color = diel, fill = diel)) + geom_density(alpha = 0.25) + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))
ggplot(pca.data, aes(x = PC4, group = diel, color = diel, fill = diel)) + geom_density(alpha = 0.25) + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))
ggplot(pca.data, aes(x = PC5, group = diel, color = diel, fill = diel)) + geom_density(alpha = 0.25) + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))

# Make a ridgeline plot
pca.data.long <- pca.data %>% gather("PC_numb", "value", PC1, PC2, PC3, PC4, PC5)

ggplot(pca.data.long, aes(y = as.factor(PC_numb), x = value, fill = diel)) + geom_density_ridges(alpha = 0.25) + theme_classic()

# Make a plot of loadings
pca.loadings <- as.data.frame(diel.pca$rotation)
pca.loadings$category <- row.names(pca.loadings)
pca.loadings <- pca.loadings %>% gather("PC_numb", "value", PC1, PC2, PC3, PC4, PC5)

ggplot(pca.loadings, aes(x = PC_numb, y = value, color = category, group = category)) + geom_point() + geom_line() + theme_classic()

ggplot(complete_cases, aes(x = benthopelagic, group = diel, color = diel, fill = diel)) + geom_density(alpha = 0.25) + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))

tb.plot <- ggplot(complete_cases, aes(x = Trophic, y = as.numeric(benthopelagic), group = diel, color = diel, fill = diel)) + geom_point() + theme_classic() + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green")) + scale_fill_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35", "green"))

hull <- tb.plot$data %>% group_by(diel) %>% slice(chull(Trophic, as.numeric(benthopelagic)))
tb.plot + geom_polygon(data = hull, aes(fill = diel, colour = diel), alpha = 0, show.legend = FALSE) + scale_color_manual(values = c("royalblue4", "goldenrod1", "mediumpurple1", "grey35")) #+ theme(legend.position = c(0.15,0.15))

