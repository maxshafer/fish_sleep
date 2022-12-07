library(ape)
library(ggtree)
library(ggplot2)
library(viridis)

# Source this gist (From "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32"), which allows using multiple of the same aes for the same plot
# Here I use it to make multiple heatmaps for Annika's figure (but should work better than the super hack I used for Figure 5 of the Cavefish single-cell paper)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid_sleep_videos/_analysis2/R_plots/")
source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")

cichlids.diel <- read.csv("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/pheno_data/combined_cichlid_data_2022-08-04.csv")
lt_phylo <- read.nexus("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/05_BEAST_RAxML.tre")
cichlids.diel <- fv2

#colnames(cichlids.diel) <- c("y", "tips", "species", "peak_amplitude", "peak", "day_night_dif")
cichlids.diel$tips[grep("Aalcal", cichlids.diel$tips)] <- "Altcal"
cichlids.diel$tips[grep("Tel'sh", cichlids.diel$tips)] <- "Telshe"

cichlids.diel$tips <- cichlids.diel$six_letter_name_Ronco

cichlids.diel$x <- 11

subset <- keep.tip(lt_phylo, tip = cichlids.diel$tips[cichlids.diel$tips %in% lt_phylo$tip.label])
d <- fortify(subset)
cichlids.diel <- cichlids.diel[cichlids.diel$tips %in% subset$tip.label,]
cichlids.diel$tips <- factor(cichlids.diel$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
cichlids.diel$y <- as.numeric(cichlids.diel$tips)

cichlids.diel$diel1 <- ifelse(cichlids.diel$cluster == cichlids.diel$cluster[cichlids.diel$six_letter_name_Ronco == "Boumic"], "diurnal", ifelse(cichlids.diel$cluster == cichlids.diel$cluster[cichlids.diel$six_letter_name_Ronco == "Neotoa"], "nocturnal", NA))
cichlids.diel$diel2 <- ifelse(cichlids.diel$cluster %in% cichlids.diel$cluster[cichlids.diel$six_letter_name_Ronco %in% c("Neopul", "Petpol", "Neocau", "Julmrk", "Juldik")], "crepuscular", NA)

################################################################################################################################################################################################
#### Reconstruct ancestral states  #############################################################################################################################################################
################################################################################################################################################################################################

day_night_dif <- cichlids.diel$day_night_dif
names(day_night_dif) <- cichlids.diel$tips
day_night_dif_spd <- cichlids.diel$day_night_dif_spd
names(day_night_dif_spd) <- cichlids.diel$tips
total_rest <- cichlids.diel$total_rest
names(total_rest) <- cichlids.diel$tips
peak <- cichlids.diel$peak
names(peak) <- cichlids.diel$tips
peak_amp <- cichlids.diel$peak_amplitude
names(peak) <- cichlids.diel$tips

fitBM_day_night_dif <- ace(day_night_dif, subset, method = "ML", model = "BM", type = "continuous")
fitBM_day_night_dif_spd <- ace(day_night_dif_spd, subset, method = "ML", model = "BM", type = "continuous")
fitBM_total_rest <- ace(total_rest, subset, method = "ML", model = "BM", type = "continuous")
fitBM_peak <- ace(peak, subset, method = "ML", model = "BM", type = "continuous")
fitBM_peak_amp <- ace(peak_amp, subset, method = "ML", model = "BM", type = "continuous")

lik.anc <- data.frame(Day_vs_Night_activity = fitBM_day_night_dif$ace, Day_vs_Night_activity_spd = fitBM_day_night_dif_spd$ace, Total_Rest = fitBM_total_rest$ace, Peak_Percentage = fitBM_peak$ace, Peak_Amplitude = fitBM_peak_amp$ace)
lik.anc$node <- (1:nrow(lik.anc)) + length(subset$tip.label)

node.data <- data.frame(Day_vs_Night_activity = day_night_dif, Day_vs_Night_activity_spd = day_night_dif_spd, Total_Rest = total_rest, Peak_Percentage = peak, Peak_Amplitude = peak_amp, node = match(names(day_night_dif), subset$tip.label))
node.data <- rbind(node.data, lik.anc)


################################################################################################################################################################################################
#### Identify nodes for tribes  ################################################################################################################################################################
################################################################################################################################################################################################

ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab() + geom_text(aes(label = node))

lamps <- 63
ectos <- 104
trops <- 114
cyps <- 110
erets <- 52
limnos <- 42
cyphos <- 41
bouls <- 61
haplos <- 53

################################################################################################################################################################################################
#### Make Plots  ###############################################################################################################################################################################
################################################################################################################################################################################################

phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 5)

p3 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=x, fill = day_night_dif_spd, colour = diel1), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", limits = c(-50,50)) + scale_color_discrete(na.value = 'transparent')
p3 <- p3 + new_scale("size") + new_scale("colour") + geom_point(data = cichlids.diel, aes(y=y, x=x+1.75, size = total_rest, colour = total_rest)) + scale_colour_distiller(palette = "PRGn", direction = 1) 
p3 <- p3 + new_scale("fill") + new_scale("colour") + geom_tile(data = cichlids.diel, aes(y=y, x=x+3, fill = peak), width = 1, inherit.aes = FALSE) + scale_fill_distiller(palette = "BrBG", direction = 1, limits = c(0,1))
# p3 <- p3 + new_scale("fill") + geom_tile(data = cichlids.diel, aes(y=y, x=x+3, fill = peak), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "YlOrRd")
# p3 <- p3 + xlim(c(0,17))

pdf("Summary_figure_phylogram_circular.pdf", width = 10, height = 10)
p3
dev.off()

png("Summary_figure_phylogram_circular.png", width = 10, height = 10, units = "in", res = 500)
p3
dev.off()

## The below is just for a presentation, and doesn't colour a few species correctly by D/N difference (just for quick visuals)
phylo.plot <- ggtree(as.phylo(subset), size = 2) + theme_tree(bgcolor = NA) + geom_tiplab(offset = 2.5, size = 10)

p3 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=10.5, fill = day_night_dif), width = 1, height = 1, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", limits = c(-0.5,0.5), na.value = "red3")
p3 <- p3 + new_scale("colour") + geom_point(data = cichlids.diel, aes(y=y, x=x+0.5, colour = total_rest), size = 10) + scale_colour_distiller(palette = "PRGn", direction = 1) 
p3 <- p3 + xlim(c(0,15)) + theme(legend.position = c(0.1,0.9))

pdf("Summary_figure_phylogram_presentation.pdf", width = 10, height = 20)
p3
dev.off()

png("Summary_figure_phylogram_presentation.png", width = 10, height = 20, units = "in", res = 500)
p3
dev.off()

phylo.plot <- ggtree(as.phylo(subset), layout = "rectangular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 5)

p3 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=x, fill = day_night_dif_spd, colour = diel1), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", limits = c(-50,50)) + scale_color_discrete(na.value = 'transparent')
p3 <- p3 + new_scale("size") + new_scale("colour") + geom_point(data = cichlids.diel, aes(y=y, x=x+1.75, size = total_rest, colour = total_rest)) + scale_colour_distiller(palette = "PRGn", direction = 1) 
p3 <- p3 + new_scale("fill") + geom_tile(data = cichlids.diel, aes(y=y, x=x+3, fill = peak), width = 1, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "BrBG", direction = 1, limits = c(0,1))
# p3 <- p3 + new_scale("fill") + geom_tile(data = cichlids.diel, aes(y=y, x=x+3, fill = peak), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "YlOrRd")
#p3 <- p3 + xlim(c(0,17))

pdf("Summary_figure_phylogram_NEW.pdf", width = 5, height = 10)
p3
dev.off()

png("Summary_figure_phylogram_NEW.png", width = 5, height = 10, units = "in", res = 500)
p3
dev.off()

# Plot and save the tree with ancestral reconstruction colouring the edges (which are matched to a node)
#ancestral.plot <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "YlOrRd", direction = 1)


##### THIS NEEDS TO BE FIXED SO THE SCALES ARE CENTERED #######

ancestral.plot.day_night_dif <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", limits = c(-0.75,0.75))
ancestral.plot.total_rest <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Total_Rest) + scale_color_distiller(palette = "PRGn", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Total_Rest), shape = 16, size = 6) + scale_color_distiller(palette = "PRGn", direction = 1)
ancestral.plot.peak <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Peak_Percentage) + scale_color_distiller(palette = "BrBG", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Peak_Percentage), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = 1, limits = c(0,1))
ancestral.plot.peak_amp <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Peak_Amplitude) + scale_color_distiller(palette = "YlGn", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Peak_Amplitude), shape = 16, size = 6) + scale_color_distiller(palette = "YlGn", direction = 1)

ancestral.plot <- ancestral.plot.day_night_dif + ancestral.plot.total_rest + ancestral.plot.peak + ancestral.plot.peak_amp + plot_layout(nrow = 2)

################################################################################################################################################################################################
#### Save Plots  ###############################################################################################################################################################################
################################################################################################################################################################################################

pdf("Ancestral_phylograms_NEW.pdf", width = 20, height = 15)
ancestral.plot
dev.off()

png("Ancestral_phylograms_NEW.png", width = 20, height = 15, units = "in", res = 500)
ancestral.plot
dev.off()


################################################################################################################################################################################################
#### Plot Phenotypes on Big Tree  ##############################################################################################################################################################
################################################################################################################################################################################################

phylo.plot <- ggtree(as.phylo(lt_phylo), layout = "circular", right = T) + theme_tree(bgcolor = NA) + geom_tiplab(offset = 3, size = 2.5)

cichlids.diel.2 <- data.frame(X = 1:length(lt_phylo$tip.label), six_letter_name_Ronco = lt_phylo$tip.label, total_rest = NA, peak_amplitude = NA, peak = NA, day_night_dif = NA, cluster = NA, tips = as.character(lt_phylo$tip.label), x = 10.5, y = NA)
cichlids.diel.2$tips <- factor(cichlids.diel.2$tips, levels = rev(with(phylo.plot$data, label[order(y, decreasing=T)])))
cichlids.diel.2$y <- as.numeric(cichlids.diel.2$tips)
cichlids.diel.2$day_night_dif <- cichlids.diel$day_night_dif[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]
cichlids.diel.2$total_rest <- cichlids.diel$total_rest[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]
cichlids.diel.2$peak <- cichlids.diel$peak[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]
cichlids.diel.2$peak_amplitude <- cichlids.diel$peak_amplitude[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]

head(cichlids.diel.2)

p4 <- phylo.plot + geom_tile(data = cichlids.diel.2, aes(y=y, x=x, fill = day_night_dif), width = 1, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu", na.value = 'transparent', limits = c(-0.75,0.75))
p4 <- p4 + new_scale("size") + new_scale("colour") + geom_point(data = cichlids.diel.2, aes(y=y, x=x+1.25, size = total_rest, colour = total_rest)) + scale_colour_distiller(palette = "PRGn", direction = 1, na.value = 'transparent') 
p4 <- p4 + new_scale("colour") + geom_point(data = cichlids.diel.2, aes(y=y, x=x+2, colour = peak, size = 10), shape = 18) + scale_colour_distiller(palette = "BrBG", direction = 1, na.value = 'transparent', limits = c(0,1))
# p4 <- p4 + new_scale("fill") + geom_tile(data = cichlids.diel.2, aes(y=y, x=x+3, fill = peak), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "YlOrRd")
# p4 <- p4 + xlim(c(0,17))

pdf("Ancestral_phylograms_all_species.pdf", width = 12, height = 12)
p4
dev.off()

png("Ancestral_phylograms_all_species.png", width = 12, height = 12, units = "in", res = 750)
p4
dev.off()


################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

# crepuscular_ids <- c("Neosav", "Neopul", "Neooli", "Neokom", "Neohel", "Neobri")
# cathemeral_ids <- c("Neogra", "Neocyg",  "Neomar", "NeofaM", "Neocra", "Neofal")
# 
# fv3 <- cichlids.diel[cichlids.diel$six_letter_name_Ronco %in% c(crepuscular_ids, cathemeral_ids), c("six_letter_name_Ronco", "day_night_dif", "total_rest")]
# 
# isotope_data_3 <- isotope_data[isotope_data$SpeciesID %in% c(crepuscular_ids, cathemeral_ids),]
# isotope_data_3$type <- ifelse(isotope_data_3$SpeciesID %in% crepuscular_ids, "crepuscular", "non-crepuscular")
# 
# ggplot(isotope_data_3, aes(x = type, y = d15N)) + geom_jitter() + geom_boxplot() + theme_classic()

