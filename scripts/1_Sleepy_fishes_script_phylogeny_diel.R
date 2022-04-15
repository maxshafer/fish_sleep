library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)
library(geiger)
library(phytools)
library(rfishbase)
library(xlsx)


setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")


### LOAD IN GOOGLE SHEET DATA ### 

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit#gid=0'
sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
sleepy_fish$Diel.Pattern <- tolower(sleepy_fish$Diel.Pattern)

write.csv(sleepy_fish, file = "sleepy_fish_database_local.csv")

sleepy_fish <- sleepy_fish[sleepy_fish$Diel.Pattern != "",]
# sleepy_fish <- sleepy_fish[as.numeric(sleepy_fish[,8]) > 1,]
# sleepy_fish <- sleepy_fish[!(grepl("sp.", sleepy_fish$Species)),]
# sleepy_fish <- sleepy_fish[!(grepl("spp.", sleepy_fish$Species)),]

################################################################################################################################################
### FETCH DATA FROM TREE OF LIFE ### 
################################################################################################################################################

# Fetch species from tree of life using rotl package
# Seems that there are at least 77 species that are in sleepy_fish, in fishbase_df, and on OTL, but not found by tnrs_match_names??

resolved_names <- tnrs_match_names(sleepy_fish$Species.name)
resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]


# Print the number of matches
print(paste("rotl found matches for", nrow(resolved_names), "out of", nrow(sleepy_fish), "from the Sleepy fish database", sep = " "))

# Remove excess information, clean up, and add tip label ids that will match the tree
resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
resolved_names <- resolved_names[!duplicated(resolved_names$tips),]

# Fetch the taxonomic levels for all of the species and make it into a dataframe
# Get it from fishbase, much faster

fishbase_df <- load_taxa(collect = T)
fishbase_df <- as.data.frame(fishbase_df)

# combine with resolved_names

resolved_names$genus <- fishbase_df$Genus[match(resolved_names$unique_name, fishbase_df$Species)]
resolved_names$family <- fishbase_df$Family[match(resolved_names$unique_name, fishbase_df$Species)]
resolved_names$order <- fishbase_df$Order[match(resolved_names$unique_name, fishbase_df$Species)]

# Add data on traits (from sleepy_fish, and/or from fishbase)

resolved_names$diel <- sleepy_fish$Diel.Pattern[match(resolved_names$search_string, tolower(sleepy_fish$Species.name))]
resolved_names$diel <- factor(resolved_names$diel, levels = c("diurnal", "nocturnal", "unclear", "crepuscular", "crepuscular/diurnal", "crepuscular/nocturnal", "crepuscular/unclear", ""))
resolved_names$diel2 <- ifelse(resolved_names$diel == "diurnal", "diurnal", ifelse(resolved_names$diel == "nocturnal", "nocturnal", ifelse(resolved_names$diel == "crepuscular", "crepuscular", ifelse(resolved_names$diel == "crepuscular/diurnal", "crepuscular", ifelse(resolved_names$diel == "crepuscular/nocturnal", "crepuscular", ifelse(resolved_names$diel == "crepuscular/unclear", "crepuscular", "unknown"))))))
resolved_names$diel_confidence <- sleepy_fish$Confidence..5...good.behavioural.tracking.data..wild.or.lab...and.or.multiple.sources..4...some.behavioural.data..3...observational.recordings.in.wild..2...annectodatal..low.n..or.otherwise.poor.data..1..not.supported..but.purported.[match(resolved_names$search_string, tolower(sleepy_fish$Species.name))]
resolved_names$genome <- sleepy_fish$Genome..Ensembl..SRA.[match(resolved_names$search_string, tolower(sleepy_fish$Species.name))]

#### Add OR and lamellia data

lam_data <- read.xlsx("~/Downloads/msab145_supplementary_data/Supplementary Data 1.xlsx", 5, header = T)
lam_data <- lam_data[,c(1:8)]
colnames(lam_data) <- lam_data[1,]
lam_data <- lam_data[c(2:nrow(lam_data)),]

OR_data <- read.xlsx("~/Downloads/msab145_supplementary_data/Supplementary Data 1.xlsx", 2, header = T)
colnames(OR_data) <- OR_data[1,]
OR_data <- OR_data[c(2:nrow(OR_data)),]

resolved_names$lam_number <- lam_data$Lamellae_number[match(resolved_names$tips, lam_data[,1], nomatch = "")]
resolved_names$lam_discrete <- lam_data$`Binary_coding (0, 1, 2 lamellae = non-ML, 3 or more lamellae = ML)`[match(resolved_names$tips, lam_data[,1], nomatch = "")]

resolved_names$functional_OR_genes <- OR_data$Functionnal[match(resolved_names$tips, OR_data$Species, nomatch = "")]
resolved_names$pseudo_OR_genes <- OR_data$Pseudogene[match(resolved_names$tips, OR_data$Species, nomatch = "")]

# These break tol_induced_subtree
# c("incertae_sedis_inherited", "unplaced_inherited","incertae_sedis", "not_otu, incertae_sedis")
# These do not break tol_induced_subtree, but do break congruify.phylo
# c("infraspecific", "sibling_higher")
# resolved_names <- resolved_names[resolved_names$flags %out% c("incertae_sedis_inherited", "unplaced_inherited", "incertae_sedis", "not_otu, incertae_sedis", "extinct_inherited, incertae_sedis"),]

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

write.csv(resolved_names, file = "resolved_names_local.csv")

print(paste("Sleepy fish database covers ", round((length(unique(resolved_names$tips))/length(unique(fishbase_df$Species)))*100), "% of Species, ", round((length(unique(resolved_names$genus))/length(unique(fishbase_df$Genus)))*100), "% of Genuses, ", round((length(unique(resolved_names$family))/length(unique(fishbase_df$Family)))*100), "% of Families, and ", round((length(unique(resolved_names$order))/length(unique(fishbase_df$Order)))*100), "% of Orders", sep = ""))

# ## Figure out which Species, Genus, Family, Order, are missing from the datas
# '%out%' <- Negate('%in%')
# 
# species <- unique(fishbase_df$Species)[unique(fishbase_df$Species) %out% unique(resolved_names$unique_name)]
# genus <- unique(fishbase_df$Genus)[unique(fishbase_df$Genus) %out% unique(resolved_names$genus)]
# family <- unique(fishbase_df$Family)[unique(fishbase_df$Family) %out% unique(resolved_names$family)]
# order <- unique(fishbase_df$Order)[unique(fishbase_df$Order) %out% unique(resolved_names$order)]
# 
# # Find species which are from missing genuses, families, and orders
# # match the genus, then ask if it is in genus
# 
# species2 <- unique(species[fishbase_df$Genus[match(species, fishbase_df$Species)] %in% genus])
# species3 <- unique(species2[fishbase_df$Family[match(species2, fishbase_df$Species)] %in% family])
# species4 <- unique(species3[fishbase_df$Order[match(species3, fishbase_df$Species)] %in% order])
# 
# new_df <- data.frame(species = species4, genus = fishbase_df$Genus[match(species4, fishbase_df$Species)], family = fishbase_df$Family[match(species4, fishbase_df$Species)], order = fishbase_df$Order[match(species4, fishbase_df$Species)], data = species4 %in% sleepy_fish$Species.name)
# new_df_family <- unique(data.frame(genus = fishbase_df$Genus[match(species3, fishbase_df$Species)], family = fishbase_df$Family[match(species3, fishbase_df$Species)], order = fishbase_df$Order[match(species3, fishbase_df$Species)], data = species3 %in% sleepy_fish$Species.name))
# 
# new_df_family$scholar <- paste("https://scholar.google.com/scholar?start=0&q=%22", new_df_family$genus, "%22+AND+(%22circadian%22+OR+%22diel%22+OR+%22diurnal%22+OR+%22nocturnal%22+OR+%22crepuscular%22)&hl=en&as_sdt=0,5", sep = "")
# new_df_family$google <- paste("https://www.google.com/search?q=%22", new_df_family$genus, "%22+AND+%28circadian+OR+diel+OR+diurnal+OR+nocturnal+OR+crepuscular+%29&rlz=1C5CHFA_enCH822CH822&biw=1680&bih=831&sxsrf=AOaemvIvLVifugmAKccz0kD3lCJFQ9-laQ%3A1637151255625&ei=F_KUYYXQJYjjkgXZgb74BQ&oq=%22", new_df_family$genus, "%22+AND+%28circadian+OR+diel+OR+diurnal+OR+nocturnal+OR+crepuscular+%29&gs_lcp=Cgdnd3Mtd2l6EANKBAhBGAFQ5QJY5QJgpARoAnAAeACAAWqIAWqSAQMwLjGYAQCgAQKgAQHAAQE&sclient=gws-wiz&ved=0ahUKEwjFt6PYr5_0AhWIsaQKHdmAD18Q4dUDCA4&uact=5", sep = "")
# 
# library(clipr)
# write_clip(new_df_family)

# Checked all Gobiesociformes, Amiiformes

################################################################################################################################################
### FETCH AND TIME-CALIBRATE THE TREE ### 
################################################################################################################################################

# Fetch the combined tree from tree of life for the species ids found in resolved_names
# "sibling_higher" is the only flag that can be included where I can both fetch the tree and time calibrate it
'%out%' <- Negate('%in%')

# tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %out% c("incertae_sedis_inherited", "unplaced_inherited", "incertae_sedis", "not_otu", "not_otu, incertae_sedis", "extinct_inherited, incertae_sedis", "infraspecific")], label_format = "name")

tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "name")

# Time calibrate it using geiger and timetree.org
# First resolve polytomies ~randomly using multi2dr

tr <- multi2di(tr)

# Make the reference file
# Ensure that the rownames and tip.labels in the target match the species names in the reference 

reference.df <- resolved_names[resolved_names$tips %in% tr$tip.label,c("order", "family", "genus", "unique_name", "tips")]
rownames(reference.df) <- reference.df$tips

# some are duplicated, or have missing data, remove them
reference.df <- reference.df[!duplicated(reference.df$unique_name),]
reference.df <- reference.df[!is.na(reference.df$unique_name),]

# Load the timetree tree (genus level data works, but not species)
# Have download timetree data for species, genus, family, and order
# Genus level data has the most calibration points
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

timetree_order <- ape::read.tree("timetree_data/actinopterygii_order.nwk")
timetree_family <- ape::read.tree("timetree_data/actinopterygii_family.nwk")
timetree_genus <- ape::read.tree("timetree_data/actinopterygii_genus.nwk")

# Use geiger to congruify the tree, works with treePL
# This seems to work up to genus, but not species (by replacing tip.labels with the same names)
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")
geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")
geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")
geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

tr.calibrated <- geiger.genus$phy

## Save out files

saveRDS(tr.calibrated, file = "calibrated_phylo.rds")

## load back in to do by species manually

tr.calibrated <- readRDS("calibrated_phylo.rds")

# Below works if you modify the heights.phylo function
# trace(geiger:::heights.phylo, edit = TRUE)
# depth = max(xx[!(is.na(xx))])
# Also have to do it manually - which is ugh!
timetree_species <- ape::read.tree("timetree_data/actinopterygii_species.nwk")
timetree_species <- multi2di(timetree_species)
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/")

geiger.species <- congruify.phylo(reference = timetree_species, target = tr.calibrated, taxonomy = reference.df, tol = 0, scale = "treePL")
tr.calibrated <- geiger.species$phy

saveRDS(tr.calibrated, file = "calibrated_phylo.rds")


## Generate subset with Nocturnal vs Diurnal

trait.data <- data.frame(species = tr.calibrated$tip.label, diel = resolved_names$diel[match(tr.calibrated$tip.label, resolved_names$tips)]) # , marine = resolved_names$marine[match(tr.calibrated$tip.label, resolved_names$tips)]
# Create vectors including crepuscular/unclear, or not
trait.data$diel1 <- ifelse(trait.data$diel %in% c("diurnal", "crepuscular/diurnal"), "diurnal", ifelse(trait.data$diel %in% c("nocturnal", "crepuscular/nocturnal"), "nocturnal", "unclear/crepuscular"))
levels(trait.data$diel1) <- c("diurnal", "nocturnal", "unclear/crepuscular")
trait.data$diel2 <- ifelse(trait.data$diel %in% c("diurnal"), "diurnal", ifelse(trait.data$diel %in% c("nocturnal"), "nocturnal", ifelse(trait.data$diel %in% c("crepuscular", "crepuscular/diurnal", "crepuscular/nocturnal"), "crepuscular", "unclear")))
levels(trait.data$diel2) <- c("diurnal", "nocturnal", "crepuscular", "unclear")

trait.data <- trait.data[!(is.na(trait.data$diel)),]
rownames(trait.data) <- trait.data$species

trait.data$tips <- resolved_names$tips[match(trait.data$species, resolved_names$tips)]
trait.data$order <- resolved_names$order[match(trait.data$species, resolved_names$tips)]
trait.data$confidence <- resolved_names$diel_confidence[match(trait.data$species, resolved_names$tips)]

trait.data$crepuscular <- ifelse(trait.data$diel2 == "crepuscular", "crepuscular", ifelse(trait.data$confidence > 4, "non_crepuscular", NA))
trait.data$diel_continuous <- ifelse(trait.data$diel1 == "diurnal", trait.data$confidence, ifelse(trait.data$diel1 == "nocturnal", trait.data$confidence*-1, NA))

saveRDS(trait.data, file = "trait_data.rds")

 















# ######
# library(caper)
# 
# # fishbase_ecology$DietTroph
# 
# data_as_factors_df <- data.frame(diel = factor(diel.vector, levels = c("diurnal", "nocturnal")))
# rownames(data_as_factors_df) <- names(diel.vector)
# data_as_factors_df$Species <- rownames(data_as_factors_df)
# #data_as_factors_df$DietTroph <- fishbase_ecology$DietTroph[match(rownames(data_as_factors_df), str_replace(fishbase_ecology$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$SL <- fishbase_morph$SL[match(rownames(data_as_factors_df), str_replace(fishbase_morph$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$TL <- fishbase_morph$TL[match(rownames(data_as_factors_df), str_replace(fishbase_morph$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$EyeFrontX <- fishbase_morph$EyeFrontX[match(rownames(data_as_factors_df), str_replace(fishbase_morph$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$EyeEndX <- fishbase_morph$EyeEndX[match(rownames(data_as_factors_df), str_replace(fishbase_morph$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$EyeFrontY <- fishbase_morph$EyeFrontY[match(rownames(data_as_factors_df), str_replace(fishbase_morph$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$EyeEndY <- fishbase_morph$EyeEndY[match(rownames(data_as_factors_df), str_replace(fishbase_morph$Species, pattern = " ", replacement = "_"))]
# data_as_factors_df$EyeDiaX <- abs(data_as_factors_df$EyeFrontX - data_as_factors_df$EyeEndX)/data_as_factors_df$TL
# data_as_factors_df$EyeDiaY <- abs(data_as_factors_df$EyeFrontY - data_as_factors_df$EyeEndY)/data_as_factors_df$TL
# data_as_factors_df$EyeDia <- ifelse(data_as_factors_df$EyeDiaX == 0, data_as_factors_df$EyeDiaY, data_as_factors_df$EyeDiaX)
# #data_as_factors_df$lam_number <- resolved_names$lam_number[match(rownames(data_as_factors_df), resolved_names$tips)]
# 
# 
# fish_continuous_discrete_as_factors_caper <- comparative.data(trpy, data_as_factors_df, names.col = Species)
# fish_continuous_discrete_as_factors_caper
# 
# fish_troph_diel_brunch <- brunch(diel ~ DietTroph, data = fish_continuous_discrete_as_factors_caper)
# summary(fish_troph_diel_brunch)
# 
# fish_SL_diel_brunch <- brunch(diel ~ SL, data = fish_continuous_discrete_as_factors_caper)
# summary(fish_SL_diel_brunch)
# 
# fish_EyeDia_diel_brunch <- brunch(diel ~ EyeDia, data = fish_continuous_discrete_as_factors_caper)
# summary(fish_EyeDia_diel_brunch)
# 
# fish_lamnumb_diel_brunch <- brunch(diel ~ lam_number, data = fish_continuous_discrete_as_factors_caper)
# summary(fish_lamnumb_diel_brunch)
# 
# 
# ggplot(data_as_factors_df, aes(x = diel, y = EyeDia)) + geom_jitter() + geom_boxplot() + theme_classic()
# 
# 
# 
# 
# 
# ### Plot only those fish with genomes
# genome.tips <- resolved_names$tips[resolved_names$genome != ""]
# genome.tips <- genome.tips[genome.tips %in% tr.calibrated$tip.label]
# 
# trpy2 <- keep.tip(tr.calibrated, tip = genome.tips)
# 
# genome.tree <- ggtree(trpy2, layout = "circular") %<+% resolved_names[,c("tips", "diel", "order")] + geom_tiplab(color = "black", size = 3, offset = 10) + geom_tippoint(aes(color = diel), shape = 16, size = 3) + scale_color_manual(values = c("red", "blue", "yellow", "green3", "green1", "green4", "green2", "white"))
# 
# 
# 
# pdf(file = paste("outs/Figures/fish_phylogeny_diel_orders_ancestral_", length(trpy2$tip.label), "_genomes_only.pdf", sep = ""), width = 13, height = 10)
# genome.tree
# dev.off()
# 
# 
# 
# ## Check whether acanthomorpha have higher transition rates (unclear how to consistently find the LCA node #)
# 
# acanth <- extract.clade(phy = trpy, node = 1427)
# 
# ggtree(acanth, layout = "circular") %<+% resolved_names[,c("tips", "diel", "order")] + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diel), shape = 16, size = 1.5) + scale_color_manual(values = c("red", "blue", "yellow", "green3", "green1", "green4", "green2", "white"))
# 
# diel.vector <- resolved_names$diel[match(acanth$tip.label, resolved_names$tips)]
# names(diel.vector) <- acanth$tip.label
# # diel.vector[is.na(diel.vector)] <- ""
# # diel.vector[diel.vector == ""] <- "unknown"
# diel.vector[diel.vector == "crepuscular/diurnal"] <- "diurnal"
# diel.vector[diel.vector == "crepuscular/nocturnal"] <- "nocturnal"
# diel.vector[diel.vector == "crepuscular/unclear"] <- "unclear"
# diel.vector[is.na(diel.vector)] <- ""
# # diel.vector <- diel.vector[diel.vector == "diurnal" | diel.vector == "nocturnal" | diel.vector == "crepuscular"]
# diel.vector <- diel.vector[diel.vector == "diurnal" | diel.vector == "nocturnal"]
# diel.vector <- droplevels(diel.vector)
# # Subset the tree for only those tips where there is diel info
# 
# acanth <- keep.tip(acanth, tip = names(diel.vector))
# 
# # Some edge lengths are equal to 0, which causes errors with ace
# 
# acanth$edge.length[acanth$edge.length == 0] <- 0.001 
# 
# # fit an ARD based Mk model for discrete character evolution (allows different forward and backward rates), which returns ancestral liklihoods for each node
# # ARD model has the lowest log-likelihood compared with other models for discrete traits, so it fits best
# 
# fitER <- ace(diel.vector, acanth, model = "ER", type = "discrete")
# fitARD <- ace(diel.vector, acanth, model = "ARD", type = "discrete")
# 
# 
# 
# ### Other stuff, no longer used, including alternative pie chart plots, and other methods for determining ancestral states
# 
# # pies <- nodepie(lik.anc, cols=1:2, color=c("red","blue"), alpha=0.8)
# # test <- inset(phylo.plot.age.2, pies, width = 0.5, height = 0.5)
# # 
# # # plot discrete ASRs with pie charts
# # pies <- nodepie(lik.anc, cols=1:2, color=c("red","blue"), alpha=0.8)
# # # basic tree
# # p12 <- ggtree(trpy) + geom_tiplab(offset=2.5, size=0.5) + xlim(0,1000) + ylim(0,1000)
# # p12 <- inset(p12, pies, width = 0.1, height = 0.1)
# # 
# # # This works, but the piecharts are still too big
# # pdf(file = "~/Desktop/phylo_tree_sleepy_fish_calibrated_ancestral.pdf", width = 20, height = 60)
# # p12
# # dev.off()
# # 
# # 
# # 
# # plotTree(trpy,fsize=0.4, type = "fan", ftype="i")
# # nodelabels(node=1:trpy$Nnode+Ntip(trpy),pie=fitER$lik.anc,piecol=cols,cex=0.2)
# # tiplabels(pie=to.matrix(diel.vector,sort(unique(diel.vector))),piecol=cols,cex=0.2)
# # add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=-max(nodeHeights(trpy)),fsize=0.8)
# # 
# # 
# # 
# # ## simulate single stochastic character map using empirical Bayes method
# # mtree <- make.simmap(trpy, diel.vector, model="ER")
# # 
# # plot(mtree,cols,type="fan",fsize=0.8,ftype="i")
# # add.simmap.legend(colors=cols, prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(mtree)), fsize=0.8)
# # 
# # 
# # mtrees<-make.simmap(trpy, diel.vector, model="ER", nsim=100)
# # par(mfrow=c(10,10))
# # null<-sapply(mtrees,plot,colors=cols,lwd=1,ftype="off")
# # 
# # pd<-summary(mtrees,plot=FALSE)
# # pd
# # 
# # plot(pd,fsize=0.6,ftype="i")
# # 
# # 
# # 
# # 
# # 
# # ## Calculate Pagel's lambda
# # 
# # pagel <- fitDiscrete(phy = trpy, dat = diel.vector)
# # 
# # 
# # ## Run fitDiscrete function using estimated lambda transformation 
# # ## ericson nocturnal 
# # EN.lambda<- fitDiscrete(trpy, diel.vector, model = c("ARD"), transform = c("lambda")) 
# # 
# # ## Run fitDiscrete function using the white-noise transformation 
# # EN.wn<- fitDiscrete(trpy, diel.vector, model = c("ARD"), transform = c("white")) 
# # 
# # 
# # 
# # 
# # # Calculate ancestral and make new phylograms
# # simmap <- make.simmap(trpy_cladogram, diel.vector, nsim = 10)
# # 
# # cols <- setNames(c("red", "blue"),c("diurnal", "nocturnal"))
# # plotSimmap(simmap, colors = cols, fsize = 0.4)
# # 
