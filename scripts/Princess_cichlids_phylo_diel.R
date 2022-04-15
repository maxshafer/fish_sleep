library(ape)
library(ggtree)
library(ggplot2)
library(viridis)

# Source this gist (From "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32"), which allows using multiple of the same aes for the same plot
# Here I use it to make multiple heatmaps for Annika's figure (but should work better than the super hack I used for Figure 5 of the Cavefish single-cell paper)
source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")

cichlids.diel <- read.csv("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid_sleep_videos/_analysis2/_combined_20210819/combined_diel_patterns_2021-08-20_dp.csv")
lt_phylo <- read.nexus("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/05_BEAST_RAxML.tre")
cichlids.diel <- fv2

#colnames(cichlids.diel) <- c("y", "tips", "species", "peak_amplitude", "peak", "day_night_dif")
cichlids.diel$tips[grep("Aalcal", cichlids.diel$tips)] <- "Altcal"
cichlids.diel$tips[grep("Tel'sh", cichlids.diel$tips)] <- "Telshe"

cichlids.diel$tips <- rownames(cichlids.diel)

cichlids.diel$x <- 11

subset <- keep.tip(lt_phylo, tip = cichlids.diel$tips[cichlids.diel$tips %in% lt_phylo$tip.label])
d <- fortify(subset)
cichlids.diel <- cichlids.diel[cichlids.diel$tips %in% subset$tip.label,]
cichlids.diel$tips <- factor(cichlids.diel$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
cichlids.diel$y <- as.numeric(cichlids.diel$tips)

phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 4)
phylo.plot <- ggtree(as.phylo(subset), layout = "rectangular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 4)




#############
# Reconstruct ancestral states
diel.vector <- cichlids.diel$move_mean_diff
names(diel.vector) <- cichlids.diel$tips


fitBM <- ace(diel.vector, subset, method = "ML", model = "BM", type = "continuous")

lik.anc <- data.frame(Day_vs_Night_activity = fitBM$ace)
lik.anc$node <- (1:nrow(lik.anc)) + length(subset$tip.label)

node.data <- data.frame(Day_vs_Night_activity = diel.vector, node = match(names(diel.vector), subset$tip.label))
node.data <- rbind(node.data, lik.anc)

# Plot and save the tree with ancestral reconstruction colouring the edges (which are matched to a node)
#ancestral.plot <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "YlOrRd", direction = 1)

ancestral.plot <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu")

dev.new()
ancestral.plot







p1 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=x, fill = move_mean_diff), width = 2, inherit.aes = F, color = "white") + scale_fill_distiller(palette = "RdBu")

# p1 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=x, fill = day_night_dif), width = 2, inherit.aes = FALSE, color = "black") + scale_fill_distiller(palette = "RdBu") + new_scale("fill") + geom_tile(data = cichlids.diel, aes(y=y, x=x+2, fill = peak_amplitude), width = 2, color = "black") + scale_fill_distiller(palette = "Purples", direction = 1)
# p1

p2 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=x, fill = move_mean_diff), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu") + new_scale("fill") + geom_point(data = cichlids.diel, aes(y=y, x=x+2, size = total_rest, colour = total_rest)) + scale_colour_distiller(palette = "Purples", direction = 1)
p2

pdf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid_sleep_videos/_analysis2/Summary_figure_phylogram_NEW.pdf", width = 10, height = 10)
p2
dev.off()

pdf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid_sleep_videos/_analysis2/Summary_figure_phylogram_ancestral_NEW.pdf", width = 10, height = 10)
ancestral.plot
dev.off()



crepuscular_ids <- c("Neosav", "Neopul.daffodil", "Neooli", "Neokom", "Neohel", "Neobri")
cathemeral_ids <- c("Neogra", "Neocyg",  "Neomar", "NeofaM", "Neocra", "Neofal")

fv3 <- fv2[c(crepuscular_ids, cathemeral_ids), c("move_mean_diff", "total_rest")]

isotope_data_3 <- isotope_data[isotope_data$SpeciesID %in% c(crepuscular_ids, cathemeral_ids),]
isotope_data_3$type <- ifelse(isotope_data_3$SpeciesID %in% crepuscular_ids, "crepuscular", "non-crepuscular")

ggplot(isotope_data_3, aes(x = type, y = d15N)) + geom_jitter() + geom_boxplot() + theme_classic()


# cichlids.diel <- data.frame(short = c("Trioto", "Neocyg", "Neogra", "Neocra", "Neomar", "Neobri", "Neohel", "Neosav", "Neooli", "Neopul", "Astbur", "Boumic", "Cypcol", "Cyplep", "Erecya", "Neocau", "Altcal", "Neobre", "Neomul", "Neonig", "Neotoa", "Ophboo", "Telvit", "Xenspi", "Lepatt", "Neodev", "Neolon", "Petpol", "Calple", "Enamel", "Neocyl"), 
#                             diel = c("Crepuscular", "Cathemeral", "Cathemeral", "Cathemeral", "Cathemeral", "Crepuscular", "Crepuscular", "Crepuscular", "Crepuscular", "Crepuscular", "Diurnal", "Diurnal", "Diurnal", "Diurnal", "Diurnal", "Diurnal", "Nocturnal", "Nocturnal", "Nocturnal", "Nocturnal", "Nocturnal", "Nocturnal", "Nocturnal", "Nocturnal", "Crepuscular", "Crepuscular", "Crepuscular", "Crepuscular", "Cathemeral", "Cathemeral", "Cathemeral"))


# phylo.plot <- ggtree(lt_phylo, layout = "circular") + theme_tree(bgcolor = NA)
# phylo.plot <- phylo.plot %<+% cichlids.diel + geom_tiplab(aes(colour = day_night_dif)) + scale_color_viridis()
# phylo.plot <- phylo.plot  + geom_hilight(node = 1, fill = "beige")+ geom_hilight(node = 269, fill = "brown") + geom_hilight(node = 275, fill = "grey") + geom_hilight(node = 284, fill = "pink") + geom_hilight(node = 392, fill = "yellow") + geom_hilight(node = 396, fill = "darkblue") + geom_hilight(node = 405, fill = "lightblue")  + geom_hilight(node = 447, fill = "darkred") + geom_hilight(node = 450, fill = "orange") + geom_hilight(node = 457, fill = "red") + geom_hilight(node = 469, fill = "purple") + geom_hilight(node = 491, fill = "green")
# 
# dev.new()
# phylo.plot
# 
# 
# phylo.plot <- phylo.plot %<+% cichlids.diel[,c("tips", "day_night_dif")] + geom_tiplab(aes(colour = day_night_dif)) + scale_colour_continuous(high = "red", low = "blue")
# phylo.plot
# 
# 
# phylo.plot <- ggtree(as.phylo(subset)) + theme_tree(bgcolor = NA)
# phylo.plot <- phylo.plot %<+% cichlids.diel[,c("tips", "peak_amplitude")] + geom_tiplab(aes(colour = peak_amplitude)) + scale_colour_continuous(high = "red", low = "blue")
# phylo.plot
# 
# 
# 
# # Make heatmap tiles
# 
# d <- fortify(subset)
# 
# cichlids.diel$tips <- factor(cichlids.diel$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
# 
# dnd <- ggplot(cichlids.diel, aes(y = tips, x = y, fill = day_night_dif)) + geom_tile() + theme_classic()
# pa <- ggplot(cichlids.diel, aes(y = tips, x = y, fill = peak_amplitude)) + geom_tile() + theme_classic()
# 
# library(patchwork)
# phylo.plot + dnd + pa
