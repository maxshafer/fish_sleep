###############################################################################################################################################
### Run Hidden Rates Models for EXTINCTION ### 
###############################################################################################################################################

## We would want a 4state_2rate model, where 2 states are fixed for all nodes ("extinction_event", "not"), where they don't have 2 rates, and no transitions
## To test whether being in an "extinction_event" affects the rate states, and transitions (whether or not it puts a faster rate under extinction conditions)
phy <- trpy_n

# Include or not
node_heights <- nodeHeights(phy)
node.age <- node_heights[match(1:nrow(node_heights), phy$edge[,2]),1]
node.age[is.na(node.age)] <- 0
names(node.age) <- 1:nrow(node_heights)
node.age <- (node.age - max(node.age))*-1

extinction <- read.csv("Speciation_paleo.csv", header = F)
colnames(extinction) <- c("Time", "extinction")
table <- slice_max(extinction, extinction, n = 5)
table$end <- table$Time + 25

# Need to ID nodes that fall within an extinction event
# BUT this doesn't work, because the node should be labelled "diurnal_extinction", or "nocturnal_extinction", not just extinction, prob why it isn't able to reconstruct
extinction_nodes <- ifelse(node.age > table$Time[1] & node.age < table$end[1], TRUE, ifelse(node.age > table$Time[2] & node.age < table$end[2], TRUE, ifelse(node.age > table$Time[3] & node.age < table$end[3], TRUE, ifelse(node.age > table$Time[4] & node.age < table$end[4], TRUE, FALSE))))
label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
label.vector[extinction_nodes] <- 3
label.vector[!(extinction_nodes)] <- NA
#label.vector[(Ntip(phy)+1):(Ntip(phy)+Nnode(phy))] <- NA
phy$node.label <- label.vector[-c(1:Ntip(phy))]

ggtree(phy, layout = "circular") + geom_nodelab(colour = "blue", size = 5)

# Make the matrices
trait.data_n$extinction <- c(rep("radiation", length(phy$tip.label)-2), "extinction", "extinction") # UGH This makes two hagfish 'extinct', so it reconstructs the LCA as diurnal
#trait.data_n$extinction <- "radiation"
# DN_LegendAndRate <- getStateMat4Dat(data.frame(species = trait.data_n$species[1:4], diel1 = c("diurnal", "diurnal", "nocturnal", "nocturnal"), extinction = c("extinction", "radation", "extinction", "radation")))
DN_LegendAndRate <- getStateMat4Dat(trait.data_n[phy$tip.label, c("species", "diel1", "extinction")])

DN_R1 <- DN_LegendAndRate$rate.mat
# DN_R2 <- DN_LegendAndRate$rate.mat
# 
# DN_ObsStateClasses <- list(DN_R1, DN_R2)
# DN_RateClassMat <- getRateCatMat(2)
# 
# DN_FullMat <- getFullMat(DN_ObsStateClasses, DN_RateClassMat)
# plotMKmodel(DN_FullMat, rate.cat = 2, display = "square", text.scale = 0.9)

HMM_2state_1rate_FEN <- corHMM(phy = phy, data = trait.data_n[phy$tip.label, c("species", "diel1", "extinction")], rate.cat = 1, rate.mat = DN_R1, fixed.nodes = TRUE)

## The problem is, is that I'm dealing with 2 traits, where the internal nodes are fixed for one of them, I have to specify that they are nocturnal or diurnal AND extinct
## How can I code uncertainty in one trait?
## I don't think I can do it

plotMKmodel(HMM_2state_1rate_FEN)
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
### Plot the Hidden rate ancestral reconstructions (of rate and state) ### 
###############################################################################################################################################

model <- HMM_2state_1rate_FEN

lik.anc <- as.data.frame(rbind(model$tip.states, model$states))
colnames(lik.anc) <- c("diurnal_extinction", "diurnal_radiation", "nocturnal_extinction", "nocturnal_radiation")
lik.anc$node <- c(1:length(trpy_n$tip.label), (length(trpy_n$tip.label) + 1):(trpy_n$Nnode + length(trpy_n$tip.label)))

lik.anc$noc.sum <- rowSums(lik.anc[,c(3,4)])
lik.anc$di.sum <- rowSums(lik.anc[,c(1,2)])

lik.anc$Ext.sum <- rowSums(lik.anc[,c(1,3)])
lik.anc$Rad.sum <- rowSums(lik.anc[,c(2,4)])

# Plot just nocturnal vs diurnal, regardless of rate class (this is what I want to known anyway)
diel_2state_2rate <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = di.sum - noc.sum) + geom_tippoint(aes(color = di.sum - noc.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")

rate_plot <- ggtree(trpy_n, layout = "circular") %<+% lik.anc + aes(color = Ext.sum) + geom_tippoint(aes(color = Ext.sum), shape = 16, size = 1.5) + scale_color_distiller(palette = "PRGn") + scale_color_distiller(palette = "PRGn")
