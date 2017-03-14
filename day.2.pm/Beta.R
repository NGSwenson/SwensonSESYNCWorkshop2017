setwd("/Users/nateswenson/sesync.workshop/data")

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
###########    CHAPTER 5  ################
##########################################
##########################################
##########################################
############## PHYLOGENETIC BETA DIVERSITY
#### TREE BASED

## UNIFRAC PRES-ABS WEIGHTED

## PD OF COM 1
pd.1 <- pd(as.matrix(my.sample[1, ]), my.phylo)[, 1]
pd.1

## PD OF COM 2
pd.2 <- pd(as.matrix(my.sample[2, ]), my.phylo)[, 1]
pd.2

## TOTAL PD OF COM 1 + COM 2  ###GAMMA
pd.total <- pd(as.matrix(my.sample[1, ])+as.matrix(my.sample[2, ]), my.phylo)[, 1]

## CALCULATE PD OF SHARED BRANCHES
pd.shared <- (pd.1 + pd.2) - pd.total

## CALCULATE UNIFRAC
u.frac <- (pd.total - pd.shared) / pd.total
u.frac

## CALCULATE UNIFRAC VIA PICANTE
unifrac(my.sample, my.phylo)



## UNIFRAC ABUNDANCE WEIGHTED
## CALCULATE RELATIVE ABUNDANCE
my.ra.sample <- my.sample / rowSums(my.sample)

## EMPTY MATRIX WITH NROW = # OF INTERNAL NODES
node.bls <- matrix(NA, nrow = length(my.phylo$edge.length), ncol = 6)

## PLACE IN PARENT AND DAUGHTER NODES IN COLS 1 AND 2
node.bls[,1:2] <- my.phylo$edge

## PLOT PHYLO AND LOOK AT NODES
plot.phylo(my.phylo, show.tip.label = F)
nodelabels(bg = "white")
tiplabels(bg = "white")

## PLACE IN NODE LENGTHS
node.bls[,3] <- my.phylo$edge.length
head(node.bls)

## Begin the loop from row 1 to the number of rows in 
## the node.bls object
for(i in 1:dim(node.bls)[1]){
		
	## Using the terminal node number for each branch 
	## in column two of the node.bls object we pull 
	## out the taxa names subtended by the branch 
	## described in each row of the node.bls object.
		leaves.node <- tips(my.phylo, node.bls[i, 2])
		
	## Sum the total relative abundance of the 
	## species found in community 1, which is the 
	## first row in our my.ra.sample object. This is 
	## equal to the Ai/AT in the raw weighted UniFrac 
	## equation.
		node.bls[i, 4] <- sum(my.ra.sample[1, leaves.node])
		
	## Sum the total relative abundance of the 
	## species found in community 1, which is the 
	## first row in our my.ra.sample object. This is 
	## equal to the Bi/BT in the raw weighted UniFrac 
	## equation.
		node.bls[i, 5] <- sum(my.ra.sample[2, leaves.node])
	
	## Calculate the product of the branch length in 
	## column 3, the bi parameter) and the absolute 
	## value of the difference of the relative 
	## abundances in the two communities for the 
	## species subtended by the branch.	
		node.bls[i, 6] <- node.bls[i, 3] * (abs(node.bls[i, 4]-node.bls[i, 5]))

	## close loop
	}

## Calculate the raw weighted UniFrac by summing the 
## values in column 6.

raw.weighted.unifrac.output <- sum(node.bls[, 6])
raw.weighted.unifrac.output
scale.output <- matrix(NA, ncol = 3, nrow=ncol(my.ra.sample))

## Start the loop
for(j in 1:ncol(my.ra.sample)){

	## add the relative abundance of species j in 
	## both communities. This is the parenthetical 
	## calculation in the scaling factor equation.
		scale.output[j, 1] <- my.ra.sample[1, j] + my.ra.sample[2, j] 

	## Use the diagonal values produced from the 
	## phylogenetic variance-covariance matrix 
	## function vcv.phylo() to provide the distance 
	## from the root to the tip holding species j. 
	## This is parameter dj.
		scale.output[j, 2] <- vcv.phylo(my.phylo)[colnames(my.ra.sample)[j],colnames(my.ra.sample)[j]]
		
	## Multiply the branch length, dj, by the total 
	## relative abundance for species j in the two 
	## communities.
		scale.output[j, 3] <- scale.output[j, 2] * scale.output[j, 1] 

	## Close the loop.
}
	

## Calculate the scaling factor, D, by summing the 
## third column of the output.

scaling.factor.output <- sum(scale.output[, 3])
scaling.factor.output
raw.weighted.unifrac.output / scaling.factor.output

## CALCULATE USING GUNIFRAC
GUniFrac(my.sample, my.phylo, alpha = 1)

## CALCULATE PHYLOSOR BY YOURSELF
2 * (pd.shared / (pd.1 + pd.2))

## CALCULATE PHYLOSOR VIA PICANTE
phylosor(my.sample, my.phylo)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
#### DISTANCE BASED

## MAKE YOUR OWN PAIRWISE BETA MEASURE
new.comdist = function(my.sample, my.dist.mat){
	get.presents = function(x){
		names(x[x>0])
			}

	list.of.names = apply(my.sample, 1, get.presents)
	
	Dpw.apply.function = function(x){
		tmp.function = function(z){
			mean(my.dist.mat[x,z])
			}
		
		lapply(list.of.names,FUN=tmp.function)
		
		}

	dpw.output = lapply(list.of.names,Dpw.apply.function)

	outt = do.call(rbind,lapply(dpw.output,unlist))
	outt
}

## MEASURE PAIRWISE BETA MEASURE FOR PHYLO AND TRAITS
new.comdist(my.sample,cophenetic(my.phylo))
new.comdist(my.sample,as.matrix(dist(my.traits)))

## MEASURE PAIRWISE BETA MEASURE FOR PHYLO VIA PICANTE
comdist(my.sample, cophenetic(my.phylo))



## MAKE YOUR OWN PAIRWISE BETA MEASURE WITH ABUNDANCE INFORMATION
new.comdist.prime = function(my.sample, my.dist.mat){

	my.ra.sample = my.sample/rowSums(my.sample)

	get.presents = function(x){
		names(x[x>0])
		}

	get.weights = function(x){
		x[x>0]
		}

	list.of.weights = apply(my.ra.sample, 1, get.weights)

	Dpw.apply.function = function(x){
		tmp.function = function(z){
			dpw.output = sum(my.dist.mat[names(x),names(z)] * outer(x, z))			
			}
		
		lapply(list.of.weights, FUN=tmp.function)
		}

	dpw.output = lapply(list.of.weights, Dpw.apply.function)

	outt = do.call(rbind,lapply(dpw.output,unlist))

	outt
}

## MEASURE PAIRWISE BETA MEASURE FOR PHYLO AND TRAITS WITH ABUNDANCE
new.comdist.prime(my.sample,cophenetic(my.phylo))
new.comdist.prime(my.sample,as.matrix(dist(my.traits)))

## MEASURE PAIRWISE BETA MEASURE FOR PHYLO AND TRAITS WITH ABUNDANCE VIA PICANTE
comdist(my.sample, cophenetic(my.phylo),abundance.weighted=T)



## MAKE YOUR OWN NEAREST NEIGHBOR BETA MEASURE
new.comdistnn = function(my.sample, my.dist.mat){

	get.presents = function(x){
		names(x[x>0])
		}

	list.of.names = apply(my.sample, 1, get.presents)
	
	Dnn.apply.function = function(x){
		tmp.function = function(z){
			mean(c(apply(my.dist.mat[x,z], MARGIN=1,min, na.rm=T), apply(my.dist.mat[x,z], MARGIN=2, min, na.rm=T)))
			}
		
		lapply(list.of.names, FUN=tmp.function)
		}

		dnn.output = lapply(list.of.names,Dnn.apply.function)
		
		outt = do.call(rbind,lapply(dnn.output,unlist))

		outt

}

## MEASURE NEAREST NEIGHBOR BETA MEASURE FOR PHYLO AND TRAITS
new.comdistnn(my.sample,cophenetic(my.phylo))
new.comdistnn(my.sample,as.matrix(dist(my.traits)))

## MEASURE NEAREST NEIGHBOR BETA MEASURE FOR PHYLO AND TRAITS VIA PICANTE
comdistnt(my.sample, cophenetic(my.phylo),abundance.weighted=F)


## MAKE YOUR OWN NEAREST NEIGHBOR BETA MEASURE WITH ABUNDANCE INFORMATION
new.comdistnn.prime = function(my.sample, my.dist.mat){

my.ra.sample = my.sample/rowSums(my.sample)

	get.presents = function(x){
		names(x[x>0])
		}

	get.weights = function(x){
		x[x>0]
		}


	list.of.weights = apply(my.ra.sample, 1, get.weights)

	Dnn.apply.function = function(x){
		tmp.function = function(z){

			weighted.mean(c(apply(as.matrix(my.dist.mat[names(x),names(z)]), MARGIN=1,min, na.rm=T), apply(as.matrix(my.dist.mat[names(x),names(z)]), MARGIN=2, min, na.rm=T)),
			c(x,z)
			)		
			
			}
		
		lapply(list.of.weights, FUN=tmp.function)
		
		}

		dnn.output = lapply(list.of.weights, Dnn.apply.function)

		outt = do.call(rbind,lapply(dnn.output,unlist))
		outt
}


## MEASURE NEAREST NEIGHBOR BETA MEASURE FOR PHYLO AND TRAITS WITH ABUNDANCE
new.comdistnn.prime(my.sample,cophenetic(my.phylo))
new.comdistnn.prime(my.sample,as.matrix(dist(my.traits)))

## MEASURE NEAREST NEIGHBOR BETA MEASURE FOR PHYLO AND TRAITS WITH ABUNDANCE VIA PICANTE
comdistnt(my.sample, cophenetic(my.phylo),abundance.weighted=T)



plot(raoD(my.sample,my.phylo)$H, as.matrix(comdist(my.sample,cophenetic(my.phylo),abundance.weighted=T)))

plot(raoD(my.sample,my.phylo)$Dkl, as.matrix(comdist(my.sample,cophenetic(my.phylo),abundance.weighted=T)))



