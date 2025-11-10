
# Section 1: Diel pattern reconstructions ----------------------------------


#load in model file

filename <- "artiodactyla_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds"
#filename <- "whippomorpha_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds"
#filename <- "ruminants_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds"

all_model_results <- readRDS(here(paste0(filename)))

#separate the results by the model types we want to use (ER, SYM, ARD, bridge_only)
#uncomment the model you want to plot

# model_results <- all_model_results$ER_model
# model_name <- "ER"

# model_results <- all_model_results$SYM_model
# model_name <- "SYM"

# model_results <- all_model_results$CONSYM_model
# model_name <- "CONSYM"

# model_results <- all_model_results$ARD_model
# model_name <- "ARD"

model_results <- all_model_results$bridge_only
#model_results <- model_results$UNTITLED
model_name <- "bridge_only"

#from the model results file, tip states describes the trait states at the tips, states describes the trait states at the nodes
lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
phylo_tree <- model_results$phy

ancestral_plot <- ggtree(phylo_tree, layout = "circular", size = 2) + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_text(aes(label=node, colour = "red"), hjust=-.2, size = 3)
ancestral_plot

#associate each of these species and their trait states with its node
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

adjusted_white <- adjustcolor("white", alpha.f = 0.0)

ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + 
  aes(color = diurnal) + 
  scale_color_gradientn(colours = c(adjustcolor("white", alpha.f = 0.0), "red"))

#plot the ancestral reconstruction, displaying each of the three trait states (cathemeral, diurnal, nocturnal)
ancestral_plot_di <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_gradientn(colours = c(adjustcolor("white", alpha.f = 0.0), "white")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
ancestral_plot_di <- ancestral_plot_di + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_di

ancestral_plot_noc <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_gradientn(colours = c(adjustcolor("#00ffff", alpha.f = 0.0), "#00ffff")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
ancestral_plot_noc <- ancestral_plot_noc + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_noc

ancestral_plot_cath <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_gradientn(colours = c(adjustcolor("yellow", alpha.f = 0.0), "yellow")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
ancestral_plot_cath <- ancestral_plot_cath + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_cath

ancestral_plot_crep <- ggtree(phylo_tree, layout = "circular", size = 2) %<+% lik.anc + aes(color = crepuscular) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5) + scale_color_gradientn(colours = c(adjustcolor("#ff00ff", alpha.f = 0.0), "#ff00ff")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5)
ancestral_plot_crep <- ancestral_plot_crep + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_crep

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_diurnal_whippo_", model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_di
dev.off()

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_noc_whippo_",model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_noc
dev.off()

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_cath_whippo_", model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_cath
dev.off()

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_crep_whippo_", model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_crep
dev.off()

#plot the ancestral reconstruction, displaying each of the three trait states (cathemeral, diurnal, nocturnal)
ancestral_plot_di <- ggtree(phylo_tree, layout = "circular", size = 1) %<+% lik.anc + aes(color = diurnal) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5) + scale_color_gradientn(colours = c(adjustcolor("white", alpha.f = 0.0), "white")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = diurnal), shape = 16, size = 1.5)
ancestral_plot_di <- ancestral_plot_di + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_di

ancestral_plot_noc <- ggtree(phylo_tree, layout = "circular", size = 1) %<+% lik.anc + aes(color = nocturnal) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)+ scale_color_gradientn(colours = c(adjustcolor("#00ffff", alpha.f = 0.0), "#00ffff")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = nocturnal), shape = 16, size = 1.5)
ancestral_plot_noc <- ancestral_plot_noc + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_noc

ancestral_plot_cath <- ggtree(phylo_tree, layout = "circular", size = 1) %<+% lik.anc + aes(color = cathemeral) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5) + scale_color_gradientn(colours = c(adjustcolor("yellow", alpha.f = 0.0), "yellow")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = cathemeral), shape = 16, size = 1.5)
ancestral_plot_cath <- ancestral_plot_cath + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_cath

ancestral_plot_crep <- ggtree(phylo_tree, layout = "circular", size = 1) %<+% lik.anc + aes(color = crepuscular) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5) + scale_color_gradientn(colours = c(adjustcolor("#ff00ff", alpha.f = 0.0), "#ff00ff")) + geom_tiplab(size = 1.5, offset = 0.5) + geom_tippoint(aes(color = crepuscular), shape = 16, size = 1.5)
ancestral_plot_crep <- ancestral_plot_crep + theme(legend.position = "none", panel.background = element_rect(fill='transparent', colour = "transparent"), plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill='transparent'))
ancestral_plot_crep

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_diurnal_artio_", model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_di
dev.off()

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_noc_artio_",model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_noc
dev.off()

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_cath_artio_", model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_cath
dev.off()

pdf(paste("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_figures/", "ancestral_recon_crep_artio_", model_name, ".pdf", sep = ""), bg = "transparent")
ancestral_plot_crep
dev.off()

# Section 2: Pie chart reconstructions ----------------------------------
#load in model data
filename <- "artiodactyla_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds"
#filename <- "whippomorpha_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds"
#filename <- "ruminants_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds"

all_model_results <- readRDS(here(filename))
model_results <- all_model_results$bridge_only_model

phylo_tree <- model_results$phy

#rename column names for consistency in the next steps
colnames(model_results$data) <- c("tips", "Diel_Pattern")

#to make more clear we can colour the tips separately using geom_tipppoint 
#may have to adjust what trait data column is called in each
custom.colours <- c("#dd8ae7", "#EECBAD" ,"#FC8D62", "#66C2A5","grey")
base_tree <- ggtree(phylo_tree, layout = "rectangular") + geom_tiplab(size = 2, hjust = -0.1)
base_tree <- base_tree %<+% model_results$data[, c("tips", "Diel_Pattern")]
#base_tree <- base_tree + geom_tippoint(aes(color = Diel_Pattern), size = 3) 
base_tree <- base_tree + geom_tile(data = base_tree$data[1:length(phylo_tree$tip.label),], aes(x=x, y=y, fill = Diel_Pattern), inherit.aes = FALSE, colour = "transparent", width = 0.5) + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
base_tree

#make the dataframe of likelihoods at the internal nodes without the tips
lik.anc <- as.data.frame(model_results$states)

lik.anc$node <- c(1:nrow(lik.anc)) + nrow(model_results$data)

#get the pie charts from this database using nodepie
#the number of columns changes depending on how many trait states
pies <- nodepie(lik.anc, 1:(length(lik.anc)-1))
pie_tree <- base_tree + geom_inset(pies, width = .03, height = .03) + scale_colour_manual(values = custom.colours)
pie_tree

bars <- nodebar(lik.anc, 1:(length(lik.anc)-1))
pie_tree <- base_tree + geom_inset(bars, width = .01, height = .02) + scale_fill_manual(values = custom.colours)
pie_tree

base_tree <- ggtree(phylo_tree, layout = "dendrogram") + geom_tiplab(size = 2)
base_tree <- base_tree %<+% model_results$data[, c("tips", "Diel_Pattern")]
base_tree <- base_tree + geom_tippoint(aes(color = Diel_Pattern), size = 3) 
pie_tree <- base_tree + geom_inset(pies, width = .03, height = .03) 
pie_tree

#circular tree
base_tree <- ggtree(phylo_tree, layout = "rectangular")
base_tree <- base_tree %<+% model_results$data[, c("tips", "Diel_Pattern")]
base_tree <- base_tree + geom_plot(data=base_tree$data, mapping=aes(x=x,y=y, label=pie), vp.width=0.09, vp.height=0.09, hjust=0.5, vjust=0.5) + coord_geo()
base_tree

library(deeptime)
library(geiger)
set.seed(12345)
tree <- sim.bdtree(b = 0.1, d = 0, stop = "time", t = 20, seed = 12345)
trait <- rTraitDisc(tree)
trait_anc <- data.frame(rerootingMethod(tree, trait)$marginal.anc)
trait_anc$node <- row.names(trait_anc)
pies <- nodepie(trait_anc, cols = 1:2, color = c("darkorange1","blue"), alpha = 0.8)
df <- tibble::tibble(node=as.numeric(trait_anc$node), pies=pies)
p1 <- ggtree(tree)
p2 <- p1 %<+% df
library(ggpp)
p2 + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.09, vp.height=0.09, hjust=0.5, vjust=0.5) + coord_geo()
#this adds a the timescale for the entire tree
#pie_tree <- pie_tree + theme_tree2()
#reverses the timescale so it starts at 0mya at the tips and extends back to 50mya at ancestor
#not currently working
#pie_tree <- revts(pie_tree)

# Section 3: Take the average of all the reconstructions? Trees? ----------------------------------
#tbd


# # Section 4: Tutorial time ----------------------------------------------
#https://blog.phytools.org/2024/04/a-few-useful-demos-on-ancestral-state.html
## load tree and data
library(phytools)
data("primate.tree")
primate.tree
data("primate.data")
head(primate.data)

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- read.csv(here("Sleepy_artiodactyla_full.csv"))
trait.data <- trait.data[trait.data$tips %in% mam.tree$tip.label,]
artio.tree <- keep.tip(mam.tree, trait.data$tips)

activity2 <- setNames(as.factor(trait.data$max_crep), trait.data$tips)
activity<-setNames(primate.data$Activity_pattern,rownames(primate.data))
head(activity)

#run each of the mk models 
primate_er<-fitMk(primate.tree,activity, model="ER",pi="fitzjohn")
primate_sym<-fitMk(primate.tree,activity,model="SYM",pi="fitzjohn")
primate_ard<-fitMk(primate.tree,activity,model="ARD",pi="fitzjohn")

#run each of the mk models,not working 
artio_er<-fitMk(artio.tree,activity2, model="ER",pi="fitzjohn")
artio_sym<-fitMk(artio.tree,activity2,model="SYM",pi="fitzjohn")
artio_ard<-fitMk(artio.tree,activity2,model="ARD",pi="fitzjohn")

#After having fit each of these three models, let’s compare them using a generic anova call.
primate_aov<-anova(primate_er,primate_sym,primate_ard)

#Our analysis reveals that the strength of evidence supporting each of these three different 
#models is quite similar one to the other. As such, we might choose to use model averaging 
#in our ancestral reconstruction. What this entails is just computing a set of weighted 
#average ancestral states, in which the weights come from our Akaike weights. ancr will 
#to do this pass it our primate_aov object instead of any one of our original fitted models.

primate_ancr<-ancr(primate_aov)
primate_ancr

cols<-setNames(c("goldenrod","lightblue","black"),
               levels(activity))

node.cex<-apply(primate_ancr$ace,1,
                function(x) if(any(x>0.95)) 0.3 else 0.6)
plot(primate_ancr,
     args.plotTree=list(type="arc",arc_height=0.5,
                        fsize=0.5,offset=3),
     args.nodelabels=list(cex=node.cex,piecol=cols),
     args.tiplabels=list(cex=0.2,piecol=cols),
     legend=FALSE)
legend(0.2*max(nodeHeights(primate.tree)),
       0.3*max(nodeHeights(primate.tree)),
       levels(activity),pch=16,col=cols,
       horiz=FALSE,cex=0.8,bty="n",pt.cex=2,
       y.intersp=1.2)

# Section 5: Lineages through time plot -----------------------------------

all_model_results <- readRDS(here(paste0("artiodactyla_finalized_max_clade_cred_four_state_max_crep_traits_ER_SYM_CONSYM_ARD_bridge_only_models.rds")))
model_results <- all_model_results$bridge_only

#from the model results file, tip states describes the trait states at the tips, states describes the trait states at the nodes
lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("cathemeral", "crepuscular", "diurnal", "nocturnal")
phylo_tree <- model_results$phy
lik.anc$tips <- rownames(lik.anc)

#associate each of these species and their trait states with its node
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

lik.anc[, "max"] <- apply(lik.anc[, 1:4], 1, max)
#I don't know a faster way to do this
for(i in 1:nrow(lik.anc)){
  if(lik.anc[i, "max"] == lik.anc[i, "diurnal"]){
    lik.anc[i, "max"] <- "diurnal"}
  if(lik.anc[i, "max"] == lik.anc[i, "nocturnal"]){
    lik.anc[i, "max"] <- "nocturnal"}
  if(lik.anc[i, "max"] == lik.anc[i, "cathemeral"]){
    lik.anc[i, "max"] <- "cathemeral"}
  if(lik.anc[i, "max"] == lik.anc[i, "crepuscular"]){
    lik.anc[i, "max"] <- "crepuscular"}
}

#plot each activity patterns lineage over time
di_nodes <- lik.anc[lik.anc$max == "diurnal", ]
phylo_tree_di <- keep.tip(phylo_tree, tip = di_nodes$node)
par(bg = NA)
ltt(phylo_tree_di, plot = TRUE, col = "#FC8D62", lwd = 2, ylim = c(0,5), xlim = c(0, 65))
dev.copy(pdf, here("testplot1.pdf"))
dev.off()

noc_nodes <- lik.anc[lik.anc$max == "nocturnal", ]
phylo_tree_noc <- keep.tip(phylo_tree, tip = noc_nodes$node)
par(bg = NA)
ltt(phylo_tree_noc, plot = TRUE, col = "#66C2A5", lwd = 2, ylim = c(0,5), xlim = c(0, 65), xaxt = "n", yaxt = "n", ann = FALSE)
dev.copy(pdf, here("testplot2.pdf"))
dev.off()

cath_nodes <- lik.anc[lik.anc$max == "cathemeral", ]
phylo_tree_cath <- keep.tip(phylo_tree, tip = cath_nodes$node)
par(bg = NA)
ltt(phylo_tree_cath, plot = TRUE, col = "#dd8ae7", lwd = 2, ylim = c(0,5), xlim = c(0, 65), xaxt = "n", yaxt = "n", ann = FALSE)
dev.copy(pdf, here("testplot3.pdf"))
dev.off()

crep_nodes <- lik.anc[lik.anc$max == "crepuscular", ]
phylo_tree_crep <- keep.tip(phylo_tree, tip = crep_nodes$node)
par(bg = NA)
ltt(phylo_tree_crep, plot = TRUE, col = "#EECBAD", lwd = 2, ylim = c(0,5), xlim = c(0, 65), xaxt = "n", yaxt = "n", ann = FALSE)
dev.copy(pdf, here("testplot4.pdf"))
dev.off()


#create a multiphylo object
trees <- list(phylo_tree_di, phylo_tree_noc, phylo_tree_cath, phylo_tree_crep)
class(trees) <- "multiPhylo"

mltt.plot(trees) 

#plot with multiple lines
ltt(phylo_tree_di, lwd = 2)
for(i in 2:4) ltt.lines(trees[[i]])

ltt.lines(phylo_tree_noc, lty = 2)
