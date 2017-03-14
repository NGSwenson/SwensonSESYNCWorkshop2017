setwd("/Users/nateswenson/sesync.workshop/data")

#install.packages("Rsundials")
library(Rsundials)

#Load packages
library(GUniFrac)
library(lefse)
library(picante)
library(phytools)
library(ade4)
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

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

comdist.is <- function(x){

as.matrix(comdist(randomizeMatrix(x, null.model = "independentswap"), cophenetic(my.phylo), abundance.weighted = F))
	}

nulls <- replicate(999, comdist.is(my.sample))

hist(nulls[1, 2, ])

nulls.means <- apply(nulls, c(1:2), mean, na.rm = T)

nulls.sds <- apply(nulls, c(1:2), sd, na.rm = T)

obs <- as.matrix(comdist(my.sample, cophenetic(my.phylo), abundance.weighted = F))

ses <- (obs - nulls.means) / nulls.sds

library(abind)
obs.nulls = abind(obs,nulls)

rank.out <- array(dim=dim(obs.nulls), t(apply(apply(obs.nulls,c(1,2),rank),3,t)))

#ranks/pvals
rank.out[,,1]


new.comdist.null = function(my.sample, my.dist.mat){

row.names(my.dist.mat) = sample(row.names(my.dist.mat), length(row.names(my.dist.mat)), replace=F)

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


dpw.output =do.call(rbind,lapply(dpw.output,unlist))

}



new.comdist.prime.null = function(my.sample, my.dist.mat){

row.names(my.dist.mat) = sample(row.names(my.dist.mat), length(row.names(my.dist.mat)), replace=F)

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


dpw.output =do.call(rbind,lapply(dpw.output,unlist))

}


new.comdistnn.null = function(my.sample, my.dist.mat){

row.names(my.dist.mat) = sample(row.names(my.dist.mat), length(row.names(my.dist.mat)), replace=F)

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


dnn.output =do.call(rbind,lapply(dnn.output,unlist))

}



new.comdistnn.prime.null = function(my.sample, my.dist.mat){

row.names(my.dist.mat) = sample(row.names(my.dist.mat), length(row.names(my.dist.mat)), replace=F)

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


dnn.output = do.call(rbind,lapply(dnn.output,unlist))

}



replicate(1,new.comdist.null(my.sample,cophenetic(my.phylo)))
replicate(1,new.comdist.prime.null(my.sample,cophenetic(my.phylo)))
replicate(1,new.comdistnn.null(my.sample,cophenetic(my.phylo)))
replicate(1,new.comdistnn.prime.null(my.sample,cophenetic(my.phylo)))