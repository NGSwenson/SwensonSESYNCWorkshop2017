setwd("/Users/nateswenson/sesync.workshop/data")

install.packages("Rsundials")
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


##########################################
##########################################
###########    CHAPTER 6  ################
##########################################
##########################################
##########################################
### NULL MODELS
### RANDOMIZE COMMUNITY DATA MATRICES

# MAKE A DUMMY COMMUNITY DATA MATRIX
com.data <- matrix(c(1,0,1,0,1,0,1,1,1, 1,1,0, 1,1,1, 0,0,1), nrow = 3)
rownames(com.data) = c("com1", "com2", "com3")
colnames(com.data) = c("A", "B", "C", "D", "E", "F")
com.data

# CALCULATE COMMUNITY SPECIES RICHNESS
rowSums(com.data)

# MAKE A RANDOM COMMUNITY OF 4 SPECIES
sample(colnames(com.data), size = 4, replace = F)

# MAKE A RANDOM COMMUNITY OF 4 SPECIES 10 TIMES - A 'RICHNESS' CONSTRAINED NULL MODEL
t(replicate(10, sample(colnames(com.data), size = 4, replace = F)))

# PERFORM 5 REPLICATES OF A RICHNESS NULL MODEL
reps <- replicate(5, randomizeMatrix(com.data, null.model = "richness"))
reps

# CALCULATE RICHNESS VALUES FOR RANDOM COMMUNITIES - CONSISTENT
apply(reps, MARGIN = 3, rowSums)

# CALCULATE OCCUPANCY RATES FOR RANDOM COMMUNITIES - NOT CONSISTENT
apply(reps, MARGIN = 3, colSums)

reps

# PERFORM 5 REPLICATES OF A INDEPENDENT SWAP NULL MODEL
reps.is <- replicate(5, randomizeMatrix(com.data, null.model = "independentswap"))

# CALCULATE RICHNESS VALUES FOR RANDOM COMMUNITIES - CONSISTENT
apply(reps.is, MARGIN = 3, rowSums)

# CALCULATE OCCUPANCY RATES FOR RANDOM COMMUNITIES - CONSISTENT
apply(reps.is, MARGIN = 3, colSums)

reps.is


# REDO INDEPENDENT SWAP WITH ABUNDANCE DATA
com.data.a <- matrix(c(3,0,2,0,10,0,4,2,3, 4,5,0, 6,1,2, 0,0,7), nrow = 3)
rownames(com.data.a) <- c("com1", "com2", "com3")
colnames(com.data.a) <- c("A", "B", "C", "D", "E", "F")
com.data.a 

reps.a.is <- replicate(5, randomizeMatrix(com.data.a, null.model = "independentswap"))
apply(reps.a.is, MARGIN = 3, rowSums)
apply(reps.a.is, MARGIN = 3, colSums)
reps.a.is



################################################################################
################################################################################
################################################################################
### RANDOMIZE PHYLO DATA

# PLOT YOUR ORIGINAL PHYLO
plot(my.phylo,cex=10)

# LOOK AT SPECIES NAMES ON PHYLO
my.phylo$tip.label

# RANDOMLY SAMPLE WITHOUT REPLACEMENT SPECIES NAMES ON PHYLO
sample(my.phylo$tip.label, length(my.phylo$tip.label), replace = F)

# MAKE TOY PHYLO TO SCRAMBLE
my.phylo.rand <- my.phylo

# RANDOMIZE NAMES ON PHYLO AND LOOK AT IT
my.phylo.rand$tip.label <- sample(my.phylo$tip.label, length(my.phylo$tip.label), replace = F)
quartz()
plot(my.phylo.rand,cex=10)

# USE PICANTE'S TIP SHUFFLE FUNCTIONN TO RANODMIZE NAMES ON PHYLO
my.phylo.rand.2 <- tipShuffle(my.phylo)

quartz()

plot(my.phylo.rand.2,cex=10)

# LET'S CALCULATE A RANDOM MEAN PAIRWISE DISTANCE MEASURE BY RANDOMIZING NAMES ON A PHYLOGENY
rand.mpd.fun <- function(x){
	
	tmp.phylo <- tipShuffle(x)
	
	## we use mpd() here, but this function itself 
	## is slow due to its use of for() loops through 
    ## rows in the community data matrix. See Chapter ## 3 for an alternative and faster approach.
	mpd(my.sample, cophenetic(tmp.phylo))

}

# RUN 10 ITERATIONS OF THE RANDOMIZATION
null.output <- replicate(10, rand.mpd.fun(my.phylo))
null.output

# PLOT YOUR NULL DISTRIBUTION FOR FIRST COMMUNITY
hist(null.output[1, ])

# CALCULATE THE 'REAL' OBSERVED VALUE FOR FIRST COMMUNITY
observed.mpd = mpd(my.sample, cophenetic(my.phylo))[1]

# PLOT WHERE OBSERVED LANDS IN NULL DISTRIBUTION
abline(v = mpd(my.sample, cophenetic(my.phylo))[1], col = "red", lwd = 2)

# CALCULATE A STANDARDIZED EFFECT SIZE FOR COMMUNITY 1
ses.1 <- (mpd(my.sample, cophenetic(my.phylo))[1] - mean(null.output[1,])) / sd(null.output[1,])
ses.1

# CALCULATE THE RANK
rank.1 <- rank(c(observed.mpd, null.output[1,]))[1]

# CALCULATE THE P VALUE
p.val.1 <- rank.1 / 11

# CALCULATE A STANDARDIZED EFFECT SIZE AND P'S FOR ALL COMMUNITIES
ses.all <- (mpd(my.sample, cophenetic(my.phylo)) - apply(null.output, MARGIN = 1, mean)) / apply(null.output, MARGIN = 1, sd)
p.val.all <- apply(cbind(mpd(my.sample, cophenetic(my.phylo)), null.output), MARGIN = 1, rank)[1,] / 11


# NOW DO THE SAME EXCEPT USING FAITH'S INDEX (PD)
observed.pd <- pd(my.sample, my.phylo, include.root = FALSE)
observed.pd

rand.pd.fun <- function(x){
	tmp.phylo <- tipShuffle(x)
	pd(my.sample, tmp.phylo)[,1]
}

null.output <- replicate(10, rand.pd.fun(my.phylo))
null.output

ses.all <- (observed.pd - apply(null.output, MARGIN = 1, mean))/apply(null.output, MARGIN = 1, sd)
p.val.all <- apply(cbind(observed.pd, null.output), MARGIN = 1, rank)[1,] / 11



################################################################################
################################################################################
################################################################################
### RANDOMIZE TRAIT DATA

# MAKE A TOY TRAIT MATRIX YOU CAN SCRAMBLE
rand.my.traits <- my.traits
rand.my.traits

# SCRAMBLE NAMES ON TRAIT MATRIX BY SAMPLING WITHOUT REPLACEMENT
replicate(5, sample(rownames(rand.my.traits), length(rownames(rand.my.traits)), replace = F))

# FUNCTION TO CALCULATE A RANDOM MEAN PAIRWISE DISTANCE
trait.shuffle.funk <- function(x){
		x <- x[colnames(my.sample), ]
		rownames(x) <- sample(rownames(x), length(rownames(x)), replace = F)
		mpd(my.sample, as.matrix(dist(x[colnames(my.sample), ])))
	
	}

# REPLICATE THE NULL MODEL FOR 10 ITERATIONS
replicate(10, trait.shuffle.funk(my.traits))


shuff.constrain.nn <- function(x){
		## Calculate the maximum and minimum value for 
		## trait 1 in community
		max.sam.1 <- max(my.traits[names(x[x > 0]), 1])
		
		min.sam.1 <- min(my.traits[names(x[x > 0]), 1])

		## Calculate the maximum and minimum value for 
		## trait 2 in community
		max.sam.2 <- max(my.traits[names(x[x > 0]), 2])
		
		min.sam.2 <- min(my.traits[names(x[x > 0]), 2])

		## Calculate the maximum and minimum value for 
		## trait 3 in community
		max.sam.3 <- max(my.traits[names(x[x > 0]), 3])
		
		min.sam.3 <- min(my.traits[names(x[x > 0]), 3])

		## Get the names of all the species in your trait 
		## data that are less or equal to the maximum 			
		## trait 1 value observed in community and names 
		## of species greater than or equal to the 
		## observed minimum trait 1 value in community.
		tr.1.a <- subset(rownames(my.traits),my.traits[ ,1] <= max.sam.1)
		tr.1.b <- subset(rownames(my.traits),my.traits[ ,1] >= min.sam.1)
		
		tr.2.a <- subset(rownames(my.traits),my.traits[ ,2] <= max.sam.2)
		tr.2.b <- subset(rownames(my.traits),my.traits[ ,2] >= min.sam.2)
		
		tr.3.a <- subset(rownames(my.traits),my.traits[ ,3] <= max.sam.3)
		tr.3.b <- subset(rownames(my.traits),my.traits[ ,3] >= min.sam.3)


		## Intersect all of the names from above such 
		## that the only remaining names are species that 
		## fall within the three dimensional range for 
		## the community. This is the new species pool.
		pruned.names <- Reduce(intersect, list(tr.1.a,tr.1.b, tr.2.a, tr.2.b, tr.3.a, tr.3.b))

		## Prune sample to only include species in new 
		## species pool.
		pruned.sam <- my.sample[ , pruned.names]

		## Prune trait matrix to only include species in 
		## new species pool.
		pruned.matrix <- my.traits[pruned.names, ]

		## Shuffle the species names on the pruned trait 
		## matrix
		rownames(pruned.matrix) <- 
		sample(rownames(pruned.matrix), 
		length(rownames(pruned.matrix)), replace = F)

		## Calculate mean nearest neighbor distance. 
		mntd(pruned.sam, as.matrix(dist(pruned.matrix)), 
		abundance.weighted = F)

	}
	
constrained.nulls = apply(replicate(999,apply(my.sample, 1, 	shuff.constrain.nn)), 3, function(x){ diag(x) })



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################## FUNCTIONS TO CALCULATE NULLS ################################################
################################################################################################
################################################################################################
################################################################################################

# PD INDEPENDENT SWAP
pd.is <- function(x){
		pd(randomizeMatrix(x, null.model = "independentswap"), my.phylo )[,1]
	}

obs.null.output <- cbind(pd(my.sample, my.phylo)[,1], replicate(999, pd.is(my.sample)))

obs.rank <- apply(obs.null.output, 1, rank)[,1]
obs.rank

ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / - apply(obs.null.output, 1, sd)
ses.value

# PD INDEPENDENT SWAP VIA PICANTE
ses.pd(my.sample, my.phylo, null.model = "independentswap", runs = 999, iterations = 1000)



# PD NAME SHUFFLING
pd.shuffle <- function(x){
		pd(my.sample, tipShuffle(x))[,1]
	}

obs.null.output <- cbind(pd(my.sample, my.phylo)[,1], replicate(999, pd.shuffle(my.phylo)))

# PD NAME SHUFFLING VIA PICANTE
ses.pd(my.sample, my.phylo, null.model = "taxa.labels", runs = 999, iterations = 1000)


# MPD INDEPENDENT SWAP
mpd.is <- function(x){
	mpd(randomizeMatrix(x, null.model = "independentswap"), cophenetic(my.phylo))
}

obs.null.output <- cbind(mpd(my.sample, cophenetic(my.phylo)), replicate(999, mpd.is(my.sample)))

# MPD INDEPENDENT SWAP VIA PICANTE
ses.mpd(my.sample, cophenetic(my.phylo), null.model = "independentswap", abundance.weighted = F,  runs = 999, iterations = 1000)


# MPD NAME SHUFFLING
mpd.shuffle <- function(x){
		mpd(my.sample, cophenetic(tipShuffle(x)))
	}

mpd.shuff.trait <- function(x){
		rownames(x) <- sample(rownames(x))
		mpd(my.sample, as.matrix(dist(x)))
	}

obs.null.output <- cbind(mpd(my.sample, cophenetic(my.phylo)), replicate(999, mpd.shuffle(my.phylo)))

# MPD NAME SHUFFLING VIA PICANTE
ses.mpd(my.sample, cophenetic(my.phylo), null.model = "taxa.labels", abundance.weighted = F,  runs = 999, iterations = 1000)

# MNTD INDEPENDENT SWAP
mntd.is <- function(x){
	mntd(randomizeMatrix(x, null.model = "independentswap"), cophenetic(my.phylo))
}

obs.null.output <- cbind(mntd(my.sample, cophenetic(my.phylo)), replicate(999, mntd.is(my.sample)))

# MNTD INDEPENDENT SWAP VIA PICANTE
ses.mntd(my.sample, cophenetic(my.phylo), null.model = "independentswap", abundance.weighted = F,  runs = 999, iterations = 1000)

# MNTD NAME SHUFFLING
mntd.shuffle <- function(x){	
	mntd(my.sample, cophenetic(tipShuffle(x)))
}

mntd.shuff.trait <- function(x){
	rownames(x) = sample(rownames(x))
		mntd(my.sample, as.matrix(dist(x)))
}

obs.null.output <- cbind(mntd(my.sample, cophenetic(my.phylo)), replicate(999, mntd.shuffle(my.phylo)))

# MNTD NAME SHUFFLING VIA PICANTE
ses.mntd(my.sample, cophenetic(my.phylo), null.model = "taxa.labels", abundance.weighted = F,  runs = 999, iterations = 1000)

# LALIBERTE & LEGENDRE METRICS INDEPENDENT SWAP
 dbfd.is <- function(x){
		dbFD(my.traits[colnames(my.sample),], randomizeMatrix(my.sample, null.model = "independentswap")[] )$FDis
}

obs.null.output <- cbind(dbFD(my.traits[colnames(my.sample),], my.sample)$FDis, replicate(99, dbfd.is(my.traits)))

# LALIBERTE & LEGENDRE METRICS NAME SHUFFLING
 dbfd.shuff <- function(x){
	rownames(x) <- sample(rownames(x), 	length(rownames(x)), replace = F)
		dbFD(x[colnames(my.sample),], my.sample)$FDis
}

obs.null.output <- cbind(dbFD(my.traits[colnames(my.sample),], my.sample)$FDis, replicate(99, dbfd.shuff(my.traits)))
