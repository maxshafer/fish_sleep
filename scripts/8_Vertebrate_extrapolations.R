library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)
library(rfishbase)
library(ggrepel)
library(here)
library(xlsx)

setwd(here())

source(here("scripts", "fish_sleep_functions.R"))


## Load in the taxonomies for the different groups

# Mammals, from https://www.mammaldiversity.org/
# Downloaded on 17/08/2023

fish_tax <- as.data.frame(load_taxa(collect = T, version = "21.06"))
fish_tax <- fish_tax[,c("Species", "Order", "Family", "Genus")]
fish_tax$Species <- str_replace(fish_tax$Species, " ", "_")
colnames(fish_tax) <- c("species", "order", "family", "genus")
fish_tax$group <- "fish"

# Mammals, from https://www.mammaldiversity.org/
# Downloaded on 17/08/2023

mammal_tax <- read.csv(here("tetrapod_data", "MDD_v1.11_6649species.csv"))
mammal_tax <- mammal_tax[,c("sciName", "order", "family", "genus")]
colnames(mammal_tax) <- c("species", "order", "family", "genus")
mammal_tax$group <- "mammals"

# Amphibians, from https://amphibiaweb.org/index.html
# Downloaded on 17/08/2023

amphibian_tax <- read.csv(here("tetrapod_data", "amphib_names_17-08-2023.txt"), sep = "\t")
amphibian_tax <- amphibian_tax[,c("species", "order", "family", "genus")]
amphibian_tax$species <- paste(amphibian_tax$genus, amphibian_tax$species, sep = "_")
colnames(amphibian_tax) <- c("species", "order", "family", "genus")
amphibian_tax$group <- "amphibians"

# Sauropids (reptiles and birds), from http://www.reptile-database.org/ and https://www.birds.cornell.edu/clementschecklist/download/
# Downloaded on 17/08/2023

reptile_tax <- read.xlsx(here("tetrapod_data", "reptile_checklist_2023_07.xlsx"), sheetIndex = 1)
reptile_tax <- reptile_tax[,c("Species", "order", "Family")]
reptile_tax$genus <- str_split(reptile_tax$Species, " ")[[1]][[1]]
reptile_tax$Species <- str_replace(reptile_tax$Species, " ", "_")
colnames(reptile_tax) <- c("species", "order", "family", "genus")

aves_tax <- read.csv(here("tetrapod_data", "NEW_eBird-Clements-v2022-integrated-checklist-October-2022.csv"))
aves_tax <- aves_tax[,c("SCI_NAME", "ORDER1", "FAMILY")]
aves_tax$genus <- str_split(aves_tax$SCI_NAME, " ")[[1]][[1]]
colnames(aves_tax) <- c("species", "order", "family", "genus")
aves_tax$species <- str_replace(aves_tax$species, " ", "_")

sauropsida_tax <- rbind(reptile_tax, aves_tax)
sauropsida_tax$group <- "sauropsida"
# combine tax

taxonomies <- Reduce(rbind, list(fish_tax, mammal_tax, amphibian_tax, sauropsida_tax))
taxonomies <- taxonomies %>% mutate(., across(.cols = everything(), tolower))


## Load and combine data

fish_data <- loadTree(return = "trait_data", dataset = "fish", subset = "all")
colnames(fish_data) <- tolower(colnames(fish_data))
fish_data <- fish_data[,c("species", "diel1", "order")]
fish_data$family <- fish_tax$family[match(fish_data$species, fish_tax$species)]
fish_data$genus <- fish_tax$genus[match(fish_data$species, fish_tax$species)]
fish_data$group <- "fish"

mammal_data <- loadTree(return = "trait_data", dataset = "mammals", subset = "all")
colnames(mammal_data) <- tolower(colnames(mammal_data))
mammal_data <- mammal_data[,c("species", "diel1", "order","family","genus")]
mammal_data$group <- "mammals"

amphibian_data <- loadTree(return = "trait_data", dataset = "tetrapods", subset = "amphibians")
colnames(amphibian_data) <- tolower(colnames(amphibian_data))
amphibian_data <- amphibian_data[,c("species", "diel1", "order","family","genus")]
amphibian_data$group <- "amphibians"

sauropsida_data <- loadTree(return = "trait_data", dataset = "tetrapods", subset = "sauropsids")
colnames(sauropsida_data) <- tolower(colnames(sauropsida_data))
sauropsida_data <- sauropsida_data[,c("species", "diel1", "order","family","genus")]
sauropsida_data$group <- "sauropsida"

diel_data <- Reduce(rbind, list(fish_data, mammal_data, amphibian_data, sauropsida_data))
diel_data <- diel_data %>% mutate(., across(.cols = everything(), tolower))
colnames(diel_data) <- c("species", "diel", "order", "family", "genus", "group")

total_numb <- nrow(taxonomies)

tax_orders <- taxonomies %>% group_by(order) %>% summarise(n = n())
tax_family <- taxonomies %>% group_by(family) %>% summarise(n = n())
tax_genus <- taxonomies %>% group_by(genus) %>% summarise(n = n())
tax_group <- taxonomies %>% group_by(group) %>% summarise(n =n())


diel_data_2 <- diel_data[diel_data$species %in% taxonomies$species,]

diel_data_2$order <- taxonomies$order[match(diel_data_2$species, taxonomies$species)]
diel_data_2$family <- taxonomies$family[match(diel_data_2$species, taxonomies$species)]
diel_data_2$genus <- taxonomies$genus[match(diel_data_2$species, taxonomies$species)]
diel_data_2$group <- taxonomies$group[match(diel_data_2$species, taxonomies$species)]

diel_total_numb <- nrow(diel_data_2)
diel_order <- diel_data_2 %>% group_by(group, order, diel) %>% summarise(n = n()) %>% mutate(frac = n/sum(n))
diel_family <- diel_data_2 %>% group_by(group, family, diel) %>% summarise(n = n()) %>% mutate(frac = n/sum(n))
diel_genus <- diel_data_2 %>% group_by(group, genus, diel) %>% summarise(n = n()) %>% mutate(frac = n/sum(n))
diel_group <- diel_data_2 %>% group_by(group, diel) %>% summarise(n = n()) %>% mutate(frac = n/sum(n))

## Combine dfs

diel_order$tax_numb <- tax_orders$n[match(diel_order$order, tax_orders$order)]*diel_order$frac
diel_family$tax_numb <- tax_family$n[match(diel_family$family, tax_family$family)]*diel_family$frac
diel_genus$tax_numb <- tax_genus$n[match(diel_genus$genus, tax_genus$genus)]*diel_genus$frac
diel_group$tax_numb <- tax_group$n[match(diel_group$group, tax_group$group)]*diel_group$frac
diel_group$category <- "observed"

## Calculate totals
order_sum <- diel_order %>% group_by(group, diel) %>% summarise(n = sum(tax_numb)) %>% mutate(frac = n/sum(n), category = "order")
fam_sum <- diel_family %>% group_by(group, diel) %>% summarise(n = sum(tax_numb)) %>% mutate(frac = n/sum(n), category = "family")
genus_sum <- diel_genus %>% group_by(group, diel) %>% summarise(n = sum(tax_numb)) %>% mutate(frac = n/sum(n), category = "genus")
# group_sum <- diel_group %>% group_by(group, diel) %>% summarise(n = sum(tax_numb)) %>% mutate(frac = n/sum(n), category = "group")

summary <- Reduce(rbind, list(order_sum, fam_sum, genus_sum, diel_group[,c(1:4,6)]))

#summary %>% group_by(group, diel) %>% summarise(n_avg = mean(n))

## This makes for all vertebrates, and when shown next to mammals, it looks much better!
all_data <- summary %>% group_by(diel, category) %>% summarise(n = sum(n), group = "all")
all_data <- all_data %>% group_by(category) %>% mutate(frac = n/sum(n))
all_data <- all_data[,c("group", "diel", "n", "frac", "category")]
summary <- rbind(summary, all_data)

summary$group <- factor(summary$group, levels = c("fish", "amphibians", "mammals", "sauropsida", "all"))
summary$diel <- factor(summary$diel, levels = c("diurnal", "nocturnal"))

## This look convincing, but has a lot of wasted space
ggplot(summary, aes(x = diel, y = frac, colour = category, group = interaction(diel))) + geom_jitter() + theme_classic() + ylim(c(0,1)) + facet_wrap(~group, nrow = 1, scales = "fixed")
                                                

## This looks ok, bar plots
ggplot(summary, aes(y = interaction(category, diel, group), x = frac, fill = category, group = diel)) + geom_bar(stat="identity") + theme_classic() + xlim(c(0,1)) + facet_wrap(~group, ncol = 1, scales = "free")

man_cols <- scales::viridis_pal(option = "H")(4)[2:1]

## Can I do stacked?
summary$group <- factor(summary$group, levels = c("fish", "amphibians", "mammals", "sauropsida", "all"))
summary$category <- factor(summary$category, levels = c("order", "family", "genus", "observed"))

# Removing observed, because it doesn't look so convincing, or good?
predicted_plot <- ggplot(summary[summary$category %in% c("order", "family", "genus"),], aes(y = category, x = frac, fill = diel, group = diel)) + geom_bar(stat="identity") + theme_classic() + xlim(c(0,1)) + facet_wrap(~group, ncol = 1,  scales = "free") + scale_fill_manual(values = man_cols)

predicted_plot <- predicted_plot + theme(strip.background = element_blank(), strip.text = element_text(size = 0,margin = margin(0,0,0,0, "cm")))

pdf(file = "outs/Figures/Vertebrate_predictions.pdf", width = 4, height = 4)
predicted_plot
dev.off()



## Extract cetaceans and make Google searches

cetacean_fams <- c("Balaenidae", "Balaenopteridae", "Cetotheriidae", "Delphinidae", "Iniidae", "Kogiidae", "Lipotidae", "Monodontidae", "Phocoenidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Ziphiidae")

cetacean_species <- mammal_tax$species[tolower(mammal_tax$family) %in% tolower(cetacean_fams)]

new_df <- data.frame(species = cetacean_species)
new_df$species <- str_replace(new_df$species, "_", " ")

new_df$scholar <- paste("https://scholar.google.com/scholar?start=0&q=%22", new_df$species, "%22+AND+(%22circadian%22+OR+%22diel%22+OR+%22diurnal%22+OR+%22nocturnal%22+OR+%22crepuscular%22)&hl=en&as_sdt=0,5", sep = "")
new_df$google <- paste("https://www.google.com/search?q=%22", new_df$species, "%22+AND+%28circadian+OR+diel+OR+diurnal+OR+nocturnal+OR+crepuscular+%29&rlz=1C5CHFA_enCH822CH822&biw=1680&bih=831&sxsrf=AOaemvIvLVifugmAKccz0kD3lCJFQ9-laQ%3A1637151255625&ei=F_KUYYXQJYjjkgXZgb74BQ&oq=%22", new_df$species, "%22+AND+%28circadian+OR+diel+OR+diurnal+OR+nocturnal+OR+crepuscular+%29&gs_lcp=Cgdnd3Mtd2l6EANKBAhBGAFQ5QJY5QJgpARoAnAAeACAAWqIAWqSAQMwLjGYAQCgAQKgAQHAAQE&sclient=gws-wiz&ved=0ahUKEwjFt6PYr5_0AhWIsaQKHdmAD18Q4dUDCA4&uact=5", sep = "")

write_clip(new_df)

