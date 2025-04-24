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


## Maybe make some figures of the distribution of species by confidence and by order?
trait.data_n$diel2 <- factor(trait.data_n$diel2, levels = c("nocturnal", "diurnal", "crepuscular", "unclear"))

trait.data_n$diel_continuous <- paste(trait.data_n$confidence, trait.data_n$diel1, sep = "_")
trait.data_n$diel_continuous <- factor(trait.data_n$diel_continuous, levels = c("E_nocturnal", "D_nocturnal", "C_nocturnal", "B_nocturnal", "A_nocturnal", "A_diurnal", "B_diurnal", "C_diurnal", "D_diurnal", "E_diurnal"))

diel_conf_plot <- ggplot(trait.data_n[trait.data_n$diel1 %in% c("diurnal", "nocturnal"),], aes(x = diel_continuous, fill = factor(diel_continuous)), color = "black") + geom_bar(stat = "count") + theme_classic() + scale_fill_brewer(palette = "RdBu", direction = -1) + theme(legend.position = "none") + xlab("Nocturanl vs diurnal confidence") + ylab("# of species")
conf_by_diel_plot <- ggplot(trait.data_n, aes(x = confidence, fill = factor(diel2)), color = "black") + geom_bar(stat = "count") + theme_classic() + scale_fill_viridis_d(option = "H") + xlab("Diel pattern confidence") + ylab("# of species") + theme(legend.title = element_blank())

df <- data.frame(level = factor(c("Species", "Genus", "Family", "Order"), levels = c("Order", "Family", "Genus", "Species")), percent = c((length(unique(resolved_names$tips))/length(unique(fishbase_df$Species)))*100, (length(unique(resolved_names$genus))/length(unique(fishbase_df$Genus)))*100, (length(unique(resolved_names$family))/length(unique(fishbase_df$Family)))*100, (length(unique(resolved_names$order))/length(unique(fishbase_df$Order)))*100))
coverage_plot <- ggplot(df, aes(x = level, y = percent, fill = level), color = "black") + geom_bar(stat = "identity") + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5)) + xlab("Taxonomic rank") + ylab("% in database")


######### Can we extrapolate the numbers based on proportions?

## The frequency of each niche in my database
actual_frac <- data.frame(table(trait.data_n$diel2))
actual_frac$Perc <- actual_frac$Freq/sum(actual_frac$Freq)*100

# Maybe Family, since order is a bit messed up in fishbase? Or can do for all

# The number/percentage of species per order in my databasee
number_order <- table(trait.data_n$order, trait.data_n$diel2)
percent_order <- as.data.frame(number_order/rowSums(number_order))

# The number per taxonomic level in fishbase
phylum_numb <- table(fishbase_df$Order)
phylum_numb <- phylum_numb[names(phylum_numb) %in% percent_order$Var1]

# The predicted number per order for each diel niche
predicted_numb_order <- percent_order
predicted_numb_order$Freq <- as.numeric(predicted_numb_order$Freq*as.vector(phylum_numb))


## Same for family
number_family <- table(trait.data_n$family, trait.data_n$diel2)
percent_family <- as.data.frame(number_family/rowSums(number_family))
phylum_numb <- table(fishbase_df$Family)[unique(fishbase_df$Family) %in% trait.data_n$family[!(is.na(trait.data_n$family))]]
# phylum_numb <- phylum_numb[names(phylum_numb) %in% percent_family$Var1]
predicted_numb_family <- percent_family
predicted_numb_family$Freq <- as.numeric(predicted_numb_family$Freq*as.vector(phylum_numb))
predicted_numb_family <- predicted_numb_family[!(is.na(predicted_numb_family$Freq)),]

## Same for genus
number_genus <- table(trait.data_n$genus, trait.data_n$diel2)
percent_genus <- as.data.frame(number_genus/rowSums(number_genus))
phylum_numb <- table(fishbase_df$Genus)[unique(fishbase_df$Genus) %in% trait.data_n$genus[!(is.na(trait.data_n$genus))]]
# phylum_numb <- phylum_numb[names(phylum_numb) %in% percent_genus$Var1]
predicted_numb_genus <- percent_genus
predicted_numb_genus$Freq <- as.numeric(predicted_numb_genus$Freq*as.vector(phylum_numb))


# Summarised for each taxonomic level
# Combined and averaged across levels
order_numb <- predicted_numb_order %>% group_by(Var2) %>% summarise(Freq = sum(Freq))
family_numb <- predicted_numb_family %>% group_by(Var2) %>% summarise(Freq = sum(as.numeric(Freq)))
genus_numb <- predicted_numb_genus %>% group_by(Var2) %>% summarise(Freq = sum(Freq))


## Now predict the number
predicted_df <- data.frame(order = order_numb$Freq, family = family_numb$Freq, genus = genus_numb$Freq)
predicted_df <- sweep(predicted_df,2,colSums(predicted_df),`/`)
predicted_df$average <- rowMeans(predicted_df)
predicted_df$number <- predicted_df$average*sum(order_numb$Freq)






### ID which orders are underrepresented (with the lowest % of total species richness)
database_numb <- table(trait.data_n$order)
phylum_numb <- table(fishbase_df$Order)
phylum_numb <- phylum_numb[names(phylum_numb) %in% names(database_numb)]

## This is the frequency of nocturnal per order in my database, the percentage of each order covered, and the total # of species in per order in fishbase
percent_order2 <- percent_order[percent_order$Var2 == "nocturnal",]
percent_order2$percent_cov <- database_numb/phylum_numb
percent_order2$total_numb <- as.numeric(phylum_numb)

## I'm interested in ones that are low coverage, high nocturnal %
percent_order2$undersampled_numb <- (1-percent_order2$percent_cov)*percent_order2$Freq*percent_order2$total_numb

check_plot1 <- ggplot(percent_order2, aes(x = Freq, y = percent_cov, label = Var1, color = log(total_numb))) + geom_point(size = 3) + theme_classic() + geom_text_repel() + xlab("Frequency of temporal niche") + ylab("Fraction of sampled diversity")
check_plot2 <- ggplot(percent_order2, aes(x = log(total_numb), y = log(percent_cov), label = Var1, color = Freq)) + geom_point(size = 3) + theme_classic() + geom_text_repel() + xlab("log(Number of species per clade)") + ylab("log(Fraction of sampled diversity)") + theme(legend.position = c(0.8,0.8))

check_plot1 + check_plot2 + plot_layout(guides = "collect")



### Plot the proportions of temporal niches by order in barplot
plot.freq <- table(resolved_names$order, resolved_names$diel2)
plot.freq <- plot.freq/rowSums(plot.freq)
plot.df <- as.data.frame(table(resolved_names$order, resolved_names$diel2))
plot.df$Count <- plot.df$Freq
plot.df$Freq <- as.data.frame(plot.freq)[,3]

## Plot the estimated numbers
plot.df.predicted <- as.data.frame(predicted_numb_order)
plot.df.predicted$Var2 <- factor(plot.df.predicted$Var2, levels = c("nocturnal","diurnal","crepuscular","unclear"))

count_plot_predicted <- ggplot(plot.df.predicted, aes(y = Var1, x = Freq, fill = Var2)) + geom_bar(position = "stack", stat= "identity") + scale_fill_viridis_d(option = "H") + theme_classic()  + theme(axis.title.y = element_blank(), axis.text.y = element_blank())


plot.df$Var2 <- factor(plot.df$Var2, levels = c("nocturnal","diurnal","crepuscular","unknown"))
count_plot <- ggplot(plot.df[plot.df$Var1 %in% plot.df.predicted$Var1,], aes(y = Var1, x = log(Count))) + geom_bar(position = "stack", stat= "identity") + scale_fill_viridis_d(option = "H") + theme_classic() + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")
bar_plot <- ggplot(plot.df[plot.df$Var1 %in% plot.df.predicted$Var1,], aes(y = Var1, x = Freq, fill = Var2)) + geom_bar(position = "stack", stat= "identity") + scale_fill_viridis_d(option = "H") + theme_classic() + theme(axis.title.y = element_blank())


##############################
##############################

design <- "
ABC
DEF"

pdf("outs/Figures/Diel_database_statistics.pdf", height = 10, width = 10)
coverage_plot + conf_by_diel_plot + diel_conf_plot + bar_plot + count_plot + count_plot_predicted + plot_layout(design = design, heights = c(2,12), widths = c(5,5,10), guides = "collect")
dev.off()
 


pdf("outs/Figures/CountsAndEstimates_TemporalNiche_Orders_Inconsistancies.pdf", height = 8, width = 8)
check_plot2
dev.off()

## Maybe plot fraction sampled, vs fraction nocturnal?