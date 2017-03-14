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
###########    CHAPTER 4  ################
##########################################
##########################################
##########################################

## GET SPECIES PRESENT IN FIRST COMMUNITY
spp <- names(my.sample[1, my.sample[1,] > 0])

## GET TRAITS FOR SPECIES PRESENT IN FIRST COMMUNITY
my.traits[spp, ]

my.traits[names(my.sample[1, my.sample[1,] > 0]),]

## CALCULATE MEAN TRAIT VALUE FOR FIRST TRAIT IN FIRST COMMUNITY
mean(my.traits[names(my.sample[1, my.sample[1,] >0]), 1], na.rm = T)

## CALCULATE MEAN TRAIT VALUE FOR SECOND TRAIT IN FIRST COMMUNITY
mean(my.traits[names(my.sample[1, my.sample[1,] >0]), 2], na.rm = T)

## CALCULATE STANDARD DEVIATION OF THE TRAIT VALUES FOR FIRST TRAIT IN FIRST COMMUNITY
sd(my.traits[names(my.sample[1, my.sample[1,] >0]), 1], na.rm = T)

## CALCULATE SKEWNESS OF THE TRAIT VALUES FOR FIRST TRAIT IN FIRST COMMUNITY
skewness(my.traits[names(my.sample[1, my.sample[1,] >0]), 1], method = "moment", na.rm = T)

## CALCULATE KURTOSIS OF THE TRAIT VALUES FOR FIRST TRAIT IN FIRST COMMUNITY
kurtosis(my.traits[names(my.sample[1, my.sample[1,] >0]), 1], method = "moment", na.rm = T)


## WRITE A FUNCTION THAT CALCULATES ALL TRAIT MOMENTS FOR ALL TRAITS IN ALL COMMUNITIES
com.trait.moments <-
function(my.sample, traits){

mean.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
sd.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
skew.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
kurt.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))

skewness = function(x){ 
			m3<-sum((x-mean(x))^3)/length(x)
			s3<-sqrt(var(x))^3
			m3/s3
			}

kurtosis = function(x){
			m4 <- sum((x-mean(x))^4)/length(x)
			s4 <- var(x)^2
			m4/s4 -3


}



	for(i in 1:ncol(traits)){

		kurt.funk = function(x){ 
			kurtosis(traits[names(x[x > 0]), i])
		}

		skew.funk = function(x){ 
			skewness(traits[names(x[x >0]), i])	
		}

		mean.funk = function(x){ 
			mean(traits[names(x[x >0]), i], na.rm = T)
		}

		sd.funk = function(x){ 
			sd(traits[names(x[x >0]), i], na.rm = T)
		}


mean.output[,i] = apply(my.sample, MARGIN = 1, mean.funk)	
sd.output[,i] = apply(my.sample, MARGIN = 1, sd.funk)
skew.output[,i] = apply(my.sample, MARGIN = 1, skew.funk)
kurt.output[,i] = apply(my.sample, MARGIN = 1, kurt.funk)

	}

output = cbind(mean.output, sd.output, skew.output, kurt.output)

colnames(output)=paste(c(rep("mean",ncol(traits)), rep("sd",ncol(traits)), rep("skew",ncol(traits)),rep( "kurtosis",ncol(traits)) ), rep(names(traits), 4), sep=".")
rownames(output) = rownames(my.sample)
output



}

## RUN THE FUNCTION
com.trait.moments(my.sample,my.traits)

##########################################################################################
############################## ABUNDANCE WEIGHTED MEAN AND SD ############################
##########################################################################################
##########################################################################################
##########################################################################################


#### LET'S SEE IF WE CAN WEIGHT THE MEAN AND SD BY ABUNDANCE. FIRST WRITE A FUNCTION
com.trait.weighted <-
function(my.sample, traits){

mean.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
sd.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))


	for(i in 1:ncol(traits)){
		weight.mean.funk = function(x){
			weighted.mean(traits[names(x[x>0]), i] , x[x > 0])
		}

		weight.sd.funk = function(x){
			wt.sd(traits[names(x[x>0]), i] , x[x>0])	
		}



		mean.output[,i] = apply(my.sample, MARGIN = 1, weight.mean.funk)	
		sd.output[,i] = apply(my.sample, MARGIN = 1, weight.sd.funk)
	}


output = cbind(mean.output, sd.output)
colnames(output)=paste(c(rep("mean",ncol(traits)), rep("sd",ncol(traits)) ), rep(names(traits), 2), sep=".")
rownames(output) = rownames(my.sample)
output

}

#### RUN THE FUNCTION
com.trait.weighted(my.sample,my.traits)


##########################################################################################
############################## SOME CONSIDERATION OF UNI- AND MULTI-VARIATE ##############
##########################################################################################
##########################################################################################
##########################################################################################

### GENERATE A DISTANCE MATRIX FOR A SINGLE TRAIT (HERE TRAIT 2)
trait.2 <- as.matrix(my.traits[,2])

rownames(trait.2) = rownames(my.traits)

trait.2

my.dist.mat.2 <- dist(trait.2, method = "euclidean")

dist(my.traits, method = "euclidean")


### SCALE YOUR TRAITS

my.traits.scaled <- apply(my.traits, MARGIN = 2, scale)

### CONSIDER REDUCING THE DATA
pc = princomp(my.traits.scaled)

summary(pc)

print(pc$loadings, cutoff = 0.001)

print(pc$scores, cutoff = 0.001)

pc.scores <- pc$scores[,1:3]

rownames(pc.scores) <- rownames(my.traits)

pc.scores

pc.dist.mat <- dist(pc.scores, method = "euclidean")

pc.dist.mat

### GOWER DISTANCE FOR MIXED VARIABLES AND (SOME) MISSING DATA

gowdis(my.traits)

####### MAKE YOUR PC DISTANCE OBJECT AS A MATRIX
square.dist.mat <- as.matrix(pc.dist.mat)

####### GET DIST MATRIX CONTAINING ONLY SPECIES PRESENT IN COMMUNITY 2 & TAKE A MEAN
mean(as.dist(square.dist.mat[names(my.sample[2, my.sample[2,] > 0]),names(my.sample[2, my.sample[2,] > 0])]))

###### MEAN PAIRWISE DISTANCE PRES-ABS WEIGHTED (SAME AS mpd(, abundance.weighted=F) IN PICANTE)
Fmpd(square.dist.mat, my.sample)

###### MEAN PAIRWISE DISTANCE ABUND WEIGHTED (SAME AS mpd(, abundance.weighted=T) IN PICANTE)
Fmpd.a(square.dist.mat, my.sample)

###### MEAN NEAREST NEIGHBOR DISTANCE PRES-ABS WEIGHTED (SAME AS mntd(, abundance.weighted=F) IN PICANTE)
Fmntd(square.dist.mat, my.sample)

###### MEAN NEAREST NEIGHBOR DISTANCE ABUND WEIGHTED (SAME AS mntd(, abundance.weighted=T) IN PICANTE)
Fmntd.a(square.dist.mat, my.sample)



################ RANGE OF TRAITS ONE AT A TIME

trait.range <-
function(my.sample, traits){

	range.sub = function(x){

	com.names = names(x[x > 0])
	
	apply(traits[com.names, ], MARGIN = 2, max,na.rm=T) - apply(traits[com.names, ], MARGIN = 2, min,na.rm=T)
	}

apply(my.sample, MARGIN = 1, range.sub)

}

trait.range(my.sample,my.traits)


################ MULTIVARIATE 'RANGE' A.K.A. CONVEX HULL VOLUME

convhulln(my.traits[names(my.sample[1, my.sample[1, ]> 0, 1]),])

convhulln(my.traits[names(my.sample[3, my.sample[3, ]> 0, 1]),])

hull.output = convhulln(my.traits[names(my.sample[3, my.sample[3, ]> 0, 1]),], options="FA")

names(hull.output)
hull.output$vol
hull.function <- function(x){


	## Get the names of the species present in our 
## community.
	com.names <- names(x[x > 0])

if(length(com.names)>3){

		## Calculate the hull using all my.traits.
	convhulln(my.traits[com.names, ], options = "FA")$vol

	}
	
	else{NA}
	}
	
	
apply(my.sample, MARGIN = 1, hull.function)


######## CALCULATE CONVEX HULL (FRic) USING FD PACKAGE
dbFD(my.traits[colnames(my.sample),], my.sample)$FRic

######## FUNCTIONAL EVENNESS (FEve) USING FD PACKAGE
dbFD(my.traits[colnames(my.sample),], my.sample, w.abun = T,  stand.x = T )$FEve

######## FUNCTIONAL DIVERGENCE (FDiv) USING FD PACKAGE
dbFD(my.traits[colnames(my.sample),], my.sample, w.abun = T,  stand.x = T )$FDiv

######## FUNCTIONAL DISPERSION (FDis) USING FD PACKAGE
dbFD(my.traits[colnames(my.sample),], my.sample, w.abun = T,  stand.x = T )$FDis

## CALCULATE SPECIES RICHNESS PER COMMUNITY
richness <- rowSums(decostand(my.sample, method = "pa", MARGIN = 1))

## PAIRWISE MEASURES VIA PICANTE
pw <- mpd(my.sample, as.matrix(dist(my.traits)), abundance.weighted = F)
pw.prime <- mpd(my.sample, as.matrix(dist(my.traits)), abundance.weighted = T)

## NEAREST NEIGHBOR MEASURES VIA PICANTE
nn <- mntd(my.sample, as.matrix(dist(my.traits)), abundance.weighted = F)
nn.prime <- mntd(my.sample, as.matrix(dist(my.traits)), abundance.weighted = T)

## LALIBERTE & LEGENDRE MEASURES VIA FD
fric <- dbFD(my.traits[colnames(my.sample), ], my.sample, w.abun = T,  stand.x = T )$FRic
feve <- dbFD(my.traits[colnames(my.sample), ], my.sample, w.abun = T,  stand.x = T )$FEve
fdiv <- dbFD(my.traits[colnames(my.sample), ], my.sample, w.abun = T,  stand.x = T )$FDiv
fdis <- dbFD(my.traits[colnames(my.sample), ], my.sample, w.abun = T,  stand.x = T )$FDis

## PUT ALL OUTPUTS TOGETHER
outputs <- as.data.frame(cbind(richness, pw, pw.prime, nn, nn.prime, fric, feve, fdiv, fdis))
names(outputs) <- c("richness", "pw", "pw.prime", "nn", "nn.prime", "FRic", "FEve", "FDiv", "Fdis" )

## PLOT IT
plot(outputs, pch = 16)

## CALCULATE CORRELATIONS
cor(outputs)


raoD(my.sample,my.phylo)

plot(raoD(my.sample,my.phylo)$Dkk, mpd(my.sample,cophenetic(my.phylo)))

plot(raoD(my.sample,my.phylo)$Dkk, mpd(my.sample,cophenetic(my.phylo),abundance.weighted=T))
abline(0,2)





#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
### PHYLOGENETIC DISTANCE BASED ALPHA DIVERSITY METRICS

## SIMPLE PAIRWISE MEASURES USING PICANTE FUNCTIONS
pw <- mpd(my.sample, cophenetic(my.phylo), abundance.weighted = F)
pw.prime <- mpd(my.sample, cophenetic(my.phylo), abundance.weighted = T)

## RAO'S PAIRWISE MEASURE VIA PICANTE
Rao.pw <- raoD(my.sample, my.phylo)$Dkk

## HELMUS MEASURES VIA PICANTE
helmus.psv <- psv(my.sample, my.phylo)[,1]
helmus.pse <- pse(my.sample, my.phylo)[,1]
helmus.psr <- psr(my.sample, my.phylo)[,1]

## SIMPLE NEAREST NEIGHBOR MEASURES USING PICANTE FUNCTIONS
nn <- mntd(my.sample, cophenetic(my.phylo), abundance.weighted = F)
nn.prime <- mntd(my.sample, cophenetic(my.phylo), abundance.weighted = T)

## PUT ALL OUTPUTS TOGETHER
outputs <- as.data.frame(cbind(pw, pw.prime, Rao.pw, helmus.psv, helmus.pse, helmus.psr, nn, nn.prime))
names(outputs) <- c("MPD", "MPD.abund", "Rao", "Helmus.PSV", "Helmus.PSE", "Helmus.PSR", "MNTD", "MNTD.abund" )

## PLOT IT
plot(outputs, pch = 16)

## CALCULATE CORRELATIONS
cor(outputs)

