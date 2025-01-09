library(ape)
library(geiger)
library(nlme)
library(phytools)

## Loading required package: maps
## Loading required package: rgl


# # Brownian motion tutorial -----------------------------------------------
#from tutorial http://lukejharmon.github.io/ilhabela/instruction/2015/07/02/simulating-Brownian-motion/
t<-0:100 # time = 100 generations
sig2<-0.01 #diffusion process/variance = 0.01 per generation (average rate the trait changes by)
## first, simulate a set of random deviates
x<-rnorm(n=length(t)-1,sd=sqrt(sig2)) #set x as a normal distribution, with a mean of 0, sd of sqrt(0.01) and length t (100 generations)
## now compute their cumulative sum
x<-c(0,cumsum(x)) #starting at zero add the value from a random sampling of the normal distribution we made
plot(t,x,type="l",ylim=c(-2,2)) #plot this sum over time

#each time you run the above code it will give a new stochastic random walk
#can run many simulations of this random walk at a time to give a better idea of how a trait may evolve
nsim<-100 #we'll run 100 simulations
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2)),nsim,length(t)-1) #create the same normal distribution as before but 100 times
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum))) #calculate the cumsum 100 times
plot(t,X[1,],xlab="time",ylab="phenotype",ylim=c(-2,2),type="l") #plot the result
apply(X[2:nsim,],1,function(x,t) lines(t,x),t=t) #plot the result for each of the 100 simulations

#can get a wider on narrower spread of possible phenotypes depending on if the σ2 is larger or smaller
#examples at link

#we can simulate this random walk along the branches of a phylogeny with three steps
#1. Set root state, 2. Draw random normal deviate for each branch 3. Add along path from root to each tip to get tip values
# on a tree closely related species tend to be similar (they covary)
#we can simulate  a discrete time phylogeny to do this
t <- 100 #total time
n <- 30 #total taxa
b <- (log(n) - log(2))/t
tree <- pbtree(b=b, n=n, t=t, type = "discrete") # use pbtree from phytools, works via reject sampling (ie 13 trees rejected before finding a tree)

plotTree(tree)

## simulate evolution along each edge, get the cumulative sum from a normal distribution and apply this over the length of each branch
X<-lapply(tree$edge.length,function(x) c(0,cumsum(rnorm(n=x,sd=sqrt(sig2)))))
## reorder the edges of the tree for pre-order traversal
cw<-reorder(tree)
## now simulate on the tree
ll<-tree$edge.length+1
for(i in 1:nrow(cw$edge)){
  pp<-which(cw$edge[,2]==cw$edge[i,1])
  if(length(pp)>0) X[[i]]<-X[[i]]+X[[pp]][ll[pp]]
  else X[[i]]<-X[[i]]+X[[1]][1]
}

## get the starting and ending points of each edge for plotting
H<-nodeHeights(tree)
## plot the simulation
plot(H[1,1],X[[1]][1],ylim=range(X),xlim=range(H),xlab="time",ylab="phenotype")
for(i in 1:length(X)) lines(H[i,1]:H[i,2],X[[i]])
## add tip labels if desired (labels each of the 30 tips/species)
yy<-sapply(1:length(tree$tip.label),function(x,y) which(x==y),y=tree$edge[,2])
yy<-sapply(yy,function(x,y) y[[x]][length(y[[x]])],y=X)
text(x=max(H),y=yy,tree$tip.label)

#this simulates continuous character evolution in the simplest set of conditions:
#discrete time (not continuous time), same rate of evolution across the entire tree, no covariance in traits due to shared ancestry
#link above gives examples of how to deal with more complex models


# PIC tutorial ------------------------------------------------------------
#PIC = phylogenetic independent contrasts
#analyzes the correlation between two traits while accounting for the phylogenetic relationships between their species
#PIC is like a statistical transformation that creates independent data points (trait data) by accounting for their relatedness
#there may be two reasons why traits are correlated, aka can we predict trait Y from trait X
#reason one: the species with traits X and Y are related, reason two: X and Y tend to evolve together (ie height and weight)
#normally, measurements from one species will not be independent from another if they share an ancestry
#tree structure may make it appear like two independent traits are co-evolving
#Felsenstein's worst-case scenario, we have a tree that separates into two large rapidly radiating clades
#the members of these clades will have similar X and Y values because they are all very closely related
#this will make it appear that X and Y are correlated when we 
#modelling brownian motion across this tree gives a significant regression slope between trait X and Y using Ordinary Least Squares

# #PGLS tutorial ----------------------------------------------------------


#tutorial from https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/

#import the df of trait data and the phylogenetic tree
anoleData <- read.csv("C:/Users/ameli/Downloads/anolisDataAppended.csv")
anoleTree <- read.tree("C:/Users/ameli/Downloads/anolis.phy")

#fix the trait data by setting the species names as rownames
rownames(anoleData) <- anoleData[,"X"]

#see what the tree looks like
plot(anoleTree)
#check that the names in the df are in the tree and vice versa (geiger function)
name.check(anoleTree, anoleData)

#this shows the correlation between two continuous (?) traits
plot(anoleData[, c("awesomeness", "hostility")])

#we can quantify how strong this correlation is with PIC
# Extract columns we want to test the correlation for
host <- anoleData[, "hostility"]
awe <- anoleData[, "awesomeness"]

#associate each of the entries with the correct row name (specific species)
names(host) <- names(awe) <- rownames(anoleData)

#calculate PICs
hPic <- pic(host, anoleTree)
aPic <- pic(awe, anoleTree)

#Make a model
picModel <- lm(hPic ~ aPic - 1)

#See if the correlation is significant
summary(picModel)
#yes is is significant, p vale of less than 2.2e-16, R-squared = 0.827

# plot results
plot(hPic ~ aPic)
abline(a = 0, b = coef(picModel))

# # PGLS tutorial ---------------------------------------------------------

#we can do this in less steps by using PGLS (phylogenetic generalized least squares)
#Brownian correlation because these are both continuous traits
#Brownian motion = stochastic "random walk" where the trait has a chance of slightly increasing or decreasing at each independent step
#method = ML, maximum likelihood
pglsModel <- gls(hostility ~ awesomeness, correlation = corBrownian(phy = anoleTree),
                 data = anoleData, method = "ML")

#summarize results of the pglsModel
summary(pglsModel)

#plot each species plus the correlation
coef(pglsModel)
plot(host ~ awe)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

#pGLS can be used with a continuous and discrete predictor (like ecomorph) as well

#for a discrete and continous variable, we still use brownian motion
pglsModel2 <- gls(hostility ~ ecomorph, correlation = corBrownian(phy = anoleTree),
                  data = anoleData, method = "ML")

#use an ANOVA to determine how hostility changes based on ecomorph
anova(pglsModel2)

coef(pglsModel2)

#we can also include multiple predictors, for example ecomorph, hostility, awesomeness
pglsModel3 <- gls(hostility ~ ecomorph * awesomeness, correlation = corBrownian(phy = anoleTree),
                  data = anoleData, method = "ML")
anova(pglsModel3)

#can also use OU rather than brownian motion (example at link above)


