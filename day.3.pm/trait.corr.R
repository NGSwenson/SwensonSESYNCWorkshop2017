setwd("/Users/nateswenson/sesync.workshop/data")

install.packages("phangorn")
install.packages("caper")

library(picante)
library(phytools)
library(phangorn)
library(gplots)
library(geiger)
library(ade4)
library(caper)

#### DATA

my.phylo = read.tree("toy.phylo.txt")
my.sample = readsample("toy.cdm.txt")
my.traits = read.csv("toy.traits.csv",row.names=1)



##########################################
##########################################
###########    CHAPTER 7  ################
##########################################
##########################################
##########################################
############################
############################
############################
############## PIC #########
############################
############################
############################
pic.x <- pic(my.traits[my.phylo$tip.label, 1], my.phylo)

pic.x

pic.y <- pic(my.traits[my.phylo$tip.label, 3], my.phylo)

pic.output <- lm(pic.y ~ pic.x - 1)

summary(pic.output)

plot(pic.y ~ pic.x)

abline(pic.output)

############################
############################
############################
############# pGLM #########
############################
############################
############################
source("pglm3.4.r")


vc = vcv(my.phylo)

pglm(trait1~trait2, data=my.traits, phylomat=vc)

pglmEstLambda(trait1~trait2, data=my.traits, phylomat=vc)

############################
############################
############################
############# PEigen #######
############################
############################
############################


p.dist.mat <- cophenetic(my.phylo)

phylo.pca <- princomp(p.dist.mat)

summary(phylo.pca)

pev.mod <- lm(my.traits[my.phylo$tip.label, 2] ~ my.traits[my.phylo$tip.label, 1] + phylo.pca$scores[ ,1] + phylo.pca$scores[ ,2])

summary(pev.mod)

#library(PVR)
#phy.eig = PVRdecomp(my.phylo)

#eigen.mod = PVR(phy.eig, , my.traits[,1], method = "moran", significance = TRUE, sig.treshold = 0.05, MI.treshold = 0.05)

#eigen.mod@Selection




############################
############################
############################
############# Prma   #######
############################
############################
############################



tr1 = my.traits[,1]
tr2 = my.traits[,2]

names(tr1)=rownames(my.traits)
names(tr2)=rownames(my.traits)

phyl.RMA(tr1,tr2,my.phylo, method="BM")


############################
############################
############################
############# pTtest #######
############################
############################
############################
phyl.pairedttest(my.phylo,my.traits[,1:2],lambda=1)


############################
############################
############################
############# Panova #######
############################
############################
############################



categorical.traits = as.data.frame(c("Green","Green","Green","Blue","Blue","Blue","Red","Red","Red","Red"))
rownames(categorical.traits)=my.phylo$tip.label
colnames(categorical.traits)="flower.color"
categorical.traits

cat.traits = categorical.traits$flower.color
names(cat.traits) = my.phylo$tip.label
phylANOVA(my.phylo, cat.traits,tr1)


############################
############################
############################
############################
############################
############################
############################
############################
############################
#install.packages(c("phyloclim","phylopars"))

install.packages("phyloclim")
library(phyloclim)

data(tree)
data(PNO)

# choose summer precipitation for analysis
clim <- PNO$PrecipitationWarmestQuarter

# estimate ancestral tolerances
ac <- anc.clim(target = tree, pno = clim, n = 100)




