
install.packages("Rsundials")
library(Rsundials)

##########################################################################################
###################OLIVIER HARDY'S CONSTRAINED NULL TO ACCOUNT FOR SIGNAL IN A ABUNDANCE##
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

## SET HARDY'S MYSTICAL K VALUE. NOTE THIS IS DIFFERENT FROM BLOMBERG'S K
hardy.K <- 3

## BIN YOUR DATA BY ABUNDANCE CATEGORIES
abund.bins <- gseq(from = min(colSums(my.sample)), to = max(colSums(my.sample)), by = hardy.K)
abund.bins

## ASSIGN YOUR SPECIES TO BINS
assigned.bins <- findInterval(colSums(my.sample) , abund.bins)
assigned.bins

names(assigned.bins) <- names(colSums(my.sample))
assigned.bins

## PUT SPECIES INTO THEIR BINS
split.bins <- split(names(assigned.bins), assigned.bins)
split.bins

## SHUFFLE NAMES WITHIN BINS NOT ACROSS
hardy.shuffle <- function(x){
	sample(x, length(x), replace = F)
	}

## APPLY ACROSS YOUR LIST/BINS
lapply(split.bins, hardy.shuffle)

## PUT NAMES BACK TOGETHER
unsplit(lapply(split.bins, hardy.shuffle), assigned.bins)

assigned.bins

## DO THE RANDOMIZATION 9 TIMES
null.names <- replicate(9, unsplit(lapply(split.bins, hardy.shuffle), assigned.bins))
null.names

## MAKE A TOY/RANDOM SMAPLE
tmp.sample <- my.sample

## FUNCTION FOR HARDY NULL USING PD
hardy.pd.null <- function(x){
	colnames(tmp.sample) <- x
	pd(tmp.sample, my.phylo)[,1]
}

## FUNCTION FOR HARDY NULL USING MPD
hardy.mpd.null <- function(x){
	colnames(tmp.sample) <- x
	mpd(tmp.sample, cophenetic(my.phylo))
}

## RUN THE NULLS
apply(null.names, MARGIN = 2,  hardy.pd.null)
apply(null.names, MARGIN = 2,  hardy.mpd.null)

