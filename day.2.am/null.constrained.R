






rand.lefse.traits <- lefse.traits

rand.lefse.traits

replicate(5, sample(rownames(rand.lefse.traits), length(rownames(rand.lefse.traits)), replace = F))

trait.shuffle.funk <- function(x){

		x <- x[colnames(lefse.sample), ]
	
		rownames(x) <- sample(rownames(x), length(rownames(x)), replace = F)

		mpd(lefse.sample, as.matrix(dist(x[colnames(lefse.sample), ])))
	
	}

replicate(10, trait.shuffle.funk(lefse.traits))


shuff.constrain.nn <- function(x){

		## Calculate the maximum and minimum value for 
		## trait 1 in community
		max.sam.1 <- max(lefse.traits[names(x[x > 0]), 1])
		
		min.sam.1 <- min(lefse.traits[names(x[x > 0]), 1])

		## Calculate the maximum and minimum value for 
		## trait 2 in community
		max.sam.2 <- max(lefse.traits[names(x[x > 0]), 2])
		
		min.sam.2 <- min(lefse.traits[names(x[x > 0]), 2])

		## Calculate the maximum and minimum value for 
		## trait 3 in community
		max.sam.3 <- max(lefse.traits[names(x[x > 0]), 3])
		
		min.sam.3 <- min(lefse.traits[names(x[x > 0]), 3])

		## Get the names of all the species in your trait 
		## data that are less or equal to the maximum 			
		## trait 1 value observed in community and names 
		## of species greater than or equal to the 
		## observed minimum trait 1 value in community.
		tr.1.a <- subset(rownames(lefse.traits),lefse.traits[ ,1] <= max.sam.1)
		tr.1.b <- subset(rownames(lefse.traits),lefse.traits[ ,1] >= min.sam.1)
		
		tr.2.a <- subset(rownames(lefse.traits),lefse.traits[ ,2] <= max.sam.2)
		tr.2.b <- subset(rownames(lefse.traits),lefse.traits[ ,2] >= min.sam.2)
		
		tr.3.a <- subset(rownames(lefse.traits),lefse.traits[ ,3] <= max.sam.3)
		tr.3.b <- subset(rownames(lefse.traits),lefse.traits[ ,3] >= min.sam.3)


		## Intersect all of the names from above such 
		## that the only remaining names are species that 
		## fall within the three dimensional range for 
		## the community. This is the new species pool.
		pruned.names <- Reduce(intersect, list(tr.1.a,tr.1.b, tr.2.a, tr.2.b, tr.3.a, tr.3.b))

		## Prune sample to only include species in new 
		## species pool.
		pruned.sam <- lefse.sample[ , pruned.names]

		## Prune trait matrix to only include species in 
		## new species pool.
		pruned.matrix <- lefse.traits[pruned.names, ]

		## Shuffle the species names on the pruned trait 
		## matrix
		rownames(pruned.matrix) <- 
		sample(rownames(pruned.matrix), 
		length(rownames(pruned.matrix)), replace = F)

		## Calculate mean nearest neighbor distance. 
		mntd(pruned.sam, as.matrix(dist(pruned.matrix)), 
		abundance.weighted = F)

	}
	
constrained.nulls = apply(replicate(999,apply(lefse.sample, 1, 	shuff.constrain.nn)), 3, function(x){ diag(x) })



pd.is <- function(x){

pd(randomizeMatrix(x, null.model = "independentswap"), my.phylo )[,1]

	}

obs.null.output <- cbind(pd(lefse.sample, my.phylo)[,1], replicate(999, pd.is(lefse.sample)))

obs.rank <- apply(obs.null.output, 1, rank)[,1]

obs.rank

ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / - apply(obs.null.output, 1, sd)

ses.value

ses.pd(lefse.sample, my.phylo, null.model = "independentswap", runs = 999, iterations = 1000)

pd.shuffle <- function(x){

		pd(lefse.sample, tipShuffle(x))[,1]

	}

obs.null.output <- cbind(pd(lefse.sample, my.phylo)[,1], replicate(999, pd.shuffle(my.phylo)))

ses.pd(lefse.sample, my.phylo, null.model = "taxa.labels", runs = 999, iterations = 1000)

mpd.is <- function(x){

	mpd(randomizeMatrix(x, null.model = "independentswap"), cophenetic(my.phylo))
	
}

obs.null.output <- cbind(mpd(lefse.sample, cophenetic(my.phylo)), replicate(999, mpd.is(lefse.sample)))

ses.mpd(lefse.sample, cophenetic(my.phylo), null.model = "independentswap", abundance.weighted = F,  runs = 999, iterations = 1000)

mpd.shuffle <- function(x){

		mpd(lefse.sample, cophenetic(tipShuffle(x)))

	}

mpd.shuff.trait <- function(x){

		rownames(x) <- sample(rownames(x))
		mpd(lefse.sample, as.matrix(dist(x)))

	}

obs.null.output <- cbind(mpd(lefse.sample, cophenetic(my.phylo)), replicate(999, mpd.shuffle(my.phylo)))

ses.mpd(lefse.sample, cophenetic(my.phylo), null.model = "taxa.labels", abundance.weighted = F,  runs = 999, iterations = 1000)

mntd.is <- function(x){

	mntd(randomizeMatrix(x, null.model = "independentswap"), cophenetic(my.phylo))
	
}

obs.null.output <- cbind(mntd(lefse.sample, cophenetic(my.phylo)), replicate(999, mntd.is(lefse.sample)))

ses.mntd(lefse.sample, cophenetic(my.phylo), null.model = "independentswap", abundance.weighted = F,  runs = 999, iterations = 1000)

mntd.shuffle <- function(x){
	
	mntd(lefse.sample, cophenetic(tipShuffle(x)))
	
}

mntd.shuff.trait <- function(x){
	
	rownames(x) = sample(rownames(x))
		mntd(lefse.sample, as.matrix(dist(x)))
	
}

obs.null.output <- cbind(mntd(lefse.sample, cophenetic(my.phylo)), replicate(999, mntd.shuffle(my.phylo)))

ses.mntd(lefse.sample, cophenetic(my.phylo), null.model = "taxa.labels", abundance.weighted = F,  runs = 999, iterations = 1000)

 dbfd.is <- function(x){

dbFD(lefse.traits[colnames(lefse.sample),], randomizeMatrix(lefse.sample, null.model = "independentswap")[] )$FDis
	
}

obs.null.output <- cbind(dbFD(lefse.traits[colnames(lefse.sample),], lefse.sample)$FDis, replicate(99, dbfd.is(lefse.traits)))

 dbfd.shuff <- function(x){

	rownames(x) <- sample(rownames(x), 	length(rownames(x)), replace = F)

		dbFD(x[colnames(lefse.sample),], lefse.sample)$FDis
	
}

obs.null.output <- cbind(dbFD(lefse.traits[colnames(lefse.sample),], lefse.sample)$FDis, replicate(99, dbfd.shuff(lefse.traits)))

comdist.is <- function(x){

as.matrix(comdist(randomizeMatrix(x, null.model = "independentswap"), cophenetic(my.phylo), abundance.weighted = F))
	}

nulls <- replicate(999, comdist.is(lefse.sample))

hist(nulls[1, 2, ])

nulls.means <- apply(nulls, c(1:2), mean, na.rm = T)

nulls.sds <- apply(nulls, c(1:2), sd, na.rm = T)

obs <- as.matrix(comdist(lefse.sample, cophenetic(my.phylo), abundance.weighted = F))

ses <- (obs - nulls.means) / nulls.sds





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



replicate(1,new.comdist.null(lefse.sample,cophenetic(my.phylo)))
replicate(1,new.comdist.prime.null(lefse.sample,cophenetic(my.phylo)))
replicate(1,new.comdistnn.null(lefse.sample,cophenetic(my.phylo)))
replicate(1,new.comdistnn.prime.null(lefse.sample,cophenetic(my.phylo)))