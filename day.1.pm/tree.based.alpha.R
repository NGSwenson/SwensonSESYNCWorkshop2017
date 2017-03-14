setwd("/Users/nateswenson/sesync.workshop/data")

#install some packages
install.packages("GUniFrac")
install.packages("FD")
install.packages("geometry")
install.packages("fBasics", dependencies = T)
install.packages("SDMTools")

#Load packages
library(GUniFrac)
library(lefse)
library(picante)
library(phytools)
library(ade4)
library(phylobase)
library(fBasics)
library(geiger)
library(SDMTools)
library(fBasics)
library(geometry)
library(FD)

#### DATA

my.phylo = read.tree("toy.phylo.txt")
my.sample = readsample("toy.cdm.txt")
my.traits = read.csv("toy.traits.csv",row.names=1)

##########################################
##########################################
###########    CHAPTER 3 ALPHA DIVERSITY #
##########################################
##########################################
##########################################
##PHYLO DIVERSITY
##########################################

#plot your phylogeny to bring back good memories
plot(my.phylo)

##remove species present in your first community
my.sample[1, my.sample[1, ] > 0]

##calculate species richness
species.richness.com.1 = length(my.sample[1, my.sample[1, ] > 0])

species.richness.com.1


##transpose matrix to make it look like a 'trait matrix'
t(my.sample[1, my.sample[1, ] > 0])

##prune your phylo to only include species present from your community
treedata(my.phylo, t(my.sample[1, my.sample[1, ] > 0]))

##prune your phylo to only include species present from your community again keeping only the phylogeny and supressing warnings
pruned.tree <-treedata(my.phylo, t(my.sample[1, my.sample[1, ] > 0]), warnings = F)$phy

pruned.tree

quartz()

#plot phylo only containing species in community
plot(pruned.tree,cex=10)

#sum branch lengh
sum(pruned.tree$edge.length)

#WRITE OUR OWN PD METHOD
prune.sum.function <- function(x){

		tmp.tree <- treedata(my.phylo, x[x > 0],warnings=F)$phy
		sum(tmp.tree$edge.length)

	}


##CALCULATE PD
apply(my.sample, MARGIN = 1, prune.sum.function)

##PD FUNCTION FROM PICANTE
pd(my.sample, my.phylo)


com.1.phylo <- treedata(my.phylo, t(my.sample[1, my.sample[1, ] > 0]),warnings=F)$phy

branches <- matrix(NA, nrow = nrow(com.1.phylo$edge), ncol = 4)

branches[,1:2] <- com.1.phylo$edge

branches[,3] <- com.1.phylo$edge.length

head(branches)

plot.phylo(com.1.phylo, show.tip.label=F)

nodelabels(bg = "white")

tiplabels(bg = "white")

##FUNCTION FOR ABUNDANCE WEIGHTED FAITH'S INDEX
weighted.faith <-
function(my.phylo, my.sample){

weighted.faith.function = function(my.sub.sample){
		
			## extract the names of species in a community with an abundance greater than zero and use that information to make a pruned phylogeny for that community.
		tmp.tree = treedata(my.phylo, my.sub.sample[my.sub.sample > 0],warnings=F)$phy

	
			## Create empty branches matrix
		branches = matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)

			## Fill first two columns of the matrix with node numbers defining each edge.
		branches[,1:2] = tmp.tree$edge

			## Fill the third column with the length of each branch
		branches[,3] = tmp.tree$edge.length
	
				get.leaves<-function(x){
					leaves.node<-tips(tmp.tree,x[2])
				}
		
				## Apply the get.leaves() function to each row in the branches object. This will retrieve species names subtended by each branch (i.e. ## row) in the branches matrix
		leaves = apply(branches, MARGIN = 1, get.leaves)
	
					## Now quickly loop through each set of leaves to ## calculate the mean abundance (Ai) for those species.
				for(i in 1:length(leaves)){
					branches[i,4] = mean(my.sub.sample[leaves[[i]]], na.rm = T) 

				}

			## Lastly calculated the Weighted Faith’s Index
		nrow(tmp.tree$edge) * ((sum(branches[,3] * branches[,4])) / sum(branches[,4]))


	}
outt = apply(my.sample, MARGIN = 1, weighted.faith.function)
outt

}

## calculate weighted faith
weighted.faith(my.phylo,my.sample)
pd(my.sample, my.phylo)

my.sample2 = my.sample
my.sample2[1,1]=50

weighted.faith(my.phylo,my.sample2)


####################################################
####################################################
####################################################
####################################################
####################################################
###FUNCTIONAL ALPHA DENDROGRAM BASED################
####################################################
####################################################
####################################################
### MAKE A MULTIVARIATE TRAIT DISTANCE MATRIX
dist(my.traits,method="euclidean")


## MAKE A UPGMA HIERARCHICAL CLUSTERING DENDROGRAM AND PLOT IT
plot(hclust(dist(my.traits,method="euclidean"), method="average"))

## ASSIGN DENDROGRAM TO AN OBJECT CALLED MY.DENDRO
my.dendro = hclust(dist(my.traits,method="euclidean"), method="average")

## CALCULATE PETCHEY AND GASTON'S FD...WHAT WENT WRONG
pd(my.sample, my.dendro)

## CALCULATE PETCHEY AND GASTON'S FD
pd(my.sample, as.phylo(my.dendro))

## PLOT DENDRO AS A PHYLO
plot(as.phylo(my.dendro))

## CALCULATE ABUNDANCE WEIGHTED PETCHEY AND GASTON'S FD
weighted.faith(as.phylo(my.dendro),my.sample)

#### DO THE ABOVE WITH A SINGLE TRAIT

dist(my.traits[,1],method="euclidean")

#### PUT NAMES BACK ON DISTANCE MATRIX
uni.trait = as.matrix(dist(my.traits[,1],method="euclidean"))
rownames(uni.trait)=colnames(uni.trait)= rownames(my.traits)

uni.dendro = hclust(as.dist(uni.trait), method="average")

