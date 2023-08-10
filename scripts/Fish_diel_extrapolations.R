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

setwd(here())

source(here("scripts", "fish_sleep_functions.R"))
# Fetch the taxonomic levels from fishbase for all of the species and make it into a dataframe
fishbase_df <- load_taxa(collect = T, version = "21.06")
fishbase_df <- as.data.frame(fishbase_df)


resolved_names <- read.csv(file = "resolved_names_local.csv", row.names = "X")

dataset_variable <- "fish"
name_variable <- "all"

trpy_n <- loadTree(return = "tree", dataset = dataset_variable, subset = name_variable, custom_tips = c, only_model_data = FALSE)
trait.data_n <- loadTree(return = "trait_data", dataset = dataset_variable, subset = name_variable, only_model_data = FALSE)

trait.data_n$family <- fishbase_df$Family[match(str_replace(trait.data_n$species, "_", " "), fishbase_df$Species)]
trait.data_n$genus <- fishbase_df$Genus[match(str_replace(trait.data_n$species, "_", " "), fishbase_df$Species)]

######### Can we extrapolate the numbers based on proportions?

# Maybe Family, since order is a bit messed up in fishbase? Or can do for all


number_order <- table(trait.data_n$order, trait.data_n$diel2)
percent_order <- as.data.frame(number_order/rowSums(number_order))
phylum_numb <- table(fishbase_df$Order)
phylum_numb <- phylum_numb[names(phylum_numb) %in% percent_order$Var1]
predicted_numb_order <- percent_order
predicted_numb_order$Freq <- as.numeric(predicted_numb_order$Freq*as.vector(phylum_numb))

order_numb <- predicted_numb_order %>% group_by(Var2) %>% summarise(Freq = sum(Freq))

number_family <- table(trait.data_n$family, trait.data_n$diel2)
percent_family <- as.data.frame(number_family/rowSums(number_family))
phylum_numb <- table(fishbase_df$Family)[unique(fishbase_df$Family) %in% trait.data_n$family[!(is.na(trait.data_n$family))]]
# phylum_numb <- phylum_numb[names(phylum_numb) %in% percent_family$Var1]
predicted_numb_family <- percent_family
predicted_numb_family$Freq <- as.numeric(predicted_numb_family$Freq*as.vector(phylum_numb))
predicted_numb_family <- predicted_numb_family[!(is.na(predicted_numb_family$Freq)),]

family_numb <- predicted_numb_family %>% group_by(Var2) %>% summarise(Freq = sum(as.numeric(Freq)))

number_genus <- table(trait.data_n$genus, trait.data_n$diel2)
percent_genus <- as.data.frame(number_genus/rowSums(number_genus))
phylum_numb <- table(fishbase_df$Genus)[unique(fishbase_df$Genus) %in% trait.data_n$genus[!(is.na(trait.data_n$genus))]]
# phylum_numb <- phylum_numb[names(phylum_numb) %in% percent_genus$Var1]
predicted_numb_genus <- percent_genus
predicted_numb_genus$Freq <- as.numeric(predicted_numb_genus$Freq*as.vector(phylum_numb))

genus_numb <- predicted_numb_genus %>% group_by(Var2) %>% summarise(Freq = sum(Freq))


predicted_df <- data.frame(order = order_numb$Freq, family = family_numb$Freq, genus = genus_numb$Freq)
predicted_df$average <- rowMeans(predicted_df)

predicted_frac <- sweep(predicted_df,2,colSums(predicted_df),`/`)
predicted_frac$average <- rowMeans(predicted_frac)
predicted_frac$number <- predicted_frac$average*sum(predicted_df$order)

actual_frac <- data.frame(table(trait.data_n$diel2))
actual_frac$Perc <- actual_frac$Freq/sum(actual_frac$Freq)*100


### ID which clades are underrepresented (with the lowest % of total species richness)
database_numb <- table(trait.data_n$order)
phylum_numb <- table(fishbase_df$Order)
phylum_numb <- phylum_numb[names(phylum_numb) %in% names(database_numb)]

percent <- database_numb/phylum_numb

percent_order2 <- percent_order
percent_order2 <- percent_order2[percent_order2$Var2 == "nocturnal",]
percent_order2$percent_cov <- database_numb/phylum_numb
percent_order2$total_numb <- as.numeric(phylum_numb)

check_plot1 <- ggplot(percent_order2, aes(x = Freq, y = percent_cov, label = Var1, color = log(total_numb))) + geom_point(size = 3) + theme_classic() + geom_text_repel() + xlab("Frequency of temporal niche") + ylab("Fraction of sampled diversity")
check_plot2 <- ggplot(percent_order2, aes(x = log(total_numb), y = log(percent_cov), label = Var1, color = Freq)) + geom_point(size = 3) + theme_classic() + geom_text_repel() + xlab("log(Number of species per clade)") + ylab("log(Fraction of sampled diversity)") + theme(legend.position = c(0.8,0.8))

check_plot1 + check_plot2 + plot_layout(guides = "collect")

### Plot the proportions of temporal niches by order in barplot
plot.freq <- table(resolved_names$order, resolved_names$diel2)
plot.freq <- plot.freq/rowSums(plot.freq)
plot.df <- as.data.frame(table(resolved_names$order, resolved_names$diel2))
plot.df$Count <- plot.df$Freq
plot.df$Freq <- as.data.frame(plot.freq)[,3]
plot.df$Var2 <- factor(plot.df$Var2, levels = c("diurnal", "nocturnal", "crepuscular", "unknown"))

count_plot <- ggplot(plot.df, aes(y = Var1, x = log(Count))) + geom_bar(position = "stack", stat= "identity") + scale_fill_manual(values = c("firebrick3", "royalblue3", "goldenrod1", "lightgreen")) + theme_classic() + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")
bar_plot <- ggplot(plot.df, aes(y = Var1, x = Freq, fill = Var2)) + geom_bar(position = "stack", stat= "identity") + scale_fill_manual(values = c("firebrick3", "royalblue3", "goldenrod1", "lightgreen")) + theme_classic() + theme(axis.title.y = element_blank())

## Plot the estimated numbers
plot.df.predicted <- as.data.frame(predicted_numb_order)
plot.df.predicted$Var2 <- factor(plot.df.predicted$Var2, levels = c("diurnal", "nocturnal", "crepuscular", "unknown"))

count_plot_predicted <- ggplot(plot.df.predicted, aes(y = Var1, x = Freq, fill = Var2)) + geom_bar(position = "stack", stat= "identity") + scale_fill_manual(values = c("firebrick3", "royalblue3", "goldenrod1", "lightgreen")) + theme_classic()  + theme(axis.title.y = element_blank(), axis.text.y = element_blank())


layout <- "
ABC
DDD"

# Save pdf
pdf("outs/Figures/CountsAndEstimates_TemporalNiche_Orders.pdf", height = 15, width = 10)
bar_plot + count_plot + count_plot_predicted + check_plot2 + plot_layout(guides = "collect", widths = c(5,5,15), heights = c(7.5,5), design = layout)
dev.off()

pdf("outs/Figures/CountsAndEstimates_TemporalNiche_Orders_Inconsistancies.pdf", height = 8, width = 8)
check_plot2
dev.off()

## Maybe plot fraction sampled, vs fraction nocturnal?