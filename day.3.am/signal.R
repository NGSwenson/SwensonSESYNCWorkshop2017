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

#################################################
#################################################
#################################################
#################################################

bmPlot(my.phylo,type="BM")


tmp.tree = rcoal(4)
bmPlot(tmp.tree,type="BM", anc=0, sig2=1/1000, ngen=1000)

##########################################
##########################################
###########    Signal No Phylo ###########
##########################################
##########################################
##########################################

library(nlme)
##read in fake trait matrix with columns indicating taxonomy
tra = read.csv("tr.partitioning.csv")

#look at data
head(tra)

#linear mixed effects model 
trait.lme<-lme(trait.value~1,random= ~1 | Family/Genus/sp, data=tra, na.action=na.omit)

#variance estimation per taxonomic level
trait.vcomp<-VarCorr(trait.lme)
trait.vcomp


##########################################
##########################################
###########    Signal Categorical Trait ##
##########################################
##########################################
##########################################



categorical.traits = as.data.frame(c("Green","Green","Green","Blue","Blue","Blue","Red","Red","Red","Red"))
rownames(categorical.traits)=my.phylo$tip.label
colnames(categorical.traits)="flower.color"
categorical.traits

OBS.null <- c(NA)
for(i in 1:999){
rand.tree = tipShuffle(my.phylo)
obs=t(data.frame(categorical.traits))
obs2<-phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
OBS.null[i]<-parsimony(rand.tree,obs2,method="sankoff")

}

obs.real = parsimony(my.phylo,obs2,method="sankoff")

p.value = (rank(c(obs.real,OBS.null))[1])/1000


cat.traits2 = cbind(as.matrix(c(rep(1,3),rep(2,3),rep(3,4)),ncol=1),1)
rownames(cat.traits2)=my.phylo$tip.label

heatmap.2(cat.traits2[,1:2],Rowv=as.dendrogram(as.hclust(my.phylo)))



##########################################
##########################################
###########    Signal Mantel Test ########
##########################################
##########################################
##########################################


p.dist.mat <- cophenetic(my.phylo)

my.traits <- my.traits[row.names(p.dist.mat), ]

trait.dist.mat <- dist(my.traits, method = "euclidean")

mantel(p.dist.mat, trait.dist.mat)


##########################################
##########################################
###########    PIC Variance Test #########
##########################################
##########################################
##########################################

pic(my.traits[my.phylo$tip.label, 1], my.phylo)

var(pic(my.traits[my.phylo$tip.label, 1], my.phylo))

pic.shuffle.var <- function(x){

## Shuffle the tip labels on the phylo.
x$tip.label <- sample(x$tip.label, length(x$tip.label), replace = FALSE)

## Calculate variance of the PICs with shuffled 
## names.
var(pic(my.traits[x$tip.label, 1], x))
	}

#replicate(999, pic.shuffle.var(my.phylo))

##calculate the rank of the observed data in a null distribution of 999. Divide output by 1000 to get 'pvalue'
rank(c(var(pic(my.traits[my.phylo$tip.label,1], my.phylo)), replicate(999, pic.shuffle.var (my.phylo))))[1]  


##########################################
##########################################
###########    PIC Mean Test #############
##########################################
##########################################
##########################################

mean(pic(my.traits[my.phylo$tip.label, 1], my.phylo))

pic.shuffle.mean <- function(x){

## Shuffle the tip labels on the phylo.
x$tip.label <- sample(x$tip.label, length(x$tip.label), replace = FALSE)

## Calculate mean PIC with shuffled names.
mean(pic(my.traits[x$tip.label, 1], x))
	}

#replicate(999, pic.shuffle.mean(my.phylo))

##calculate the rank of the observed data in a null distribution of 999. Divide output by 1000 to get 'pvalue'
rank(c(mean(pic(my.traits[my.phylo$tip.label,1], my.phylo)), replicate(999, pic.shuffle.mean(my.phylo))))[1] 
##########################################
##########################################
###########    Signal Blomberg's K Test ##
##########################################
##########################################
##########################################
trait.vector <- my.traits[my.phylo$tip.label , 1]

n <- length(trait.vector)

my.vcv <- vcv.phylo(my.phylo)

inv.vcv <- solve(my.vcv)

root.value <- sum(inv.vcv %*% trait.vector)/ sum(inv.vcv)

MSEo <- (t(trait.vector - root.value) %*% (trait.vector - root.value)) / (n - 1)

MSE <- (t(trait.vector - root.value) %*% inv.vcv %*% (trait.vector - root.value)) / (n - 1)

expected.MSEo.MSE <- (sum(diag(my.vcv)) - (n / sum(inv.vcv))) / (n - 1)

(MSEo / MSE) / expected.MSEo.MSE


phylosig(my.phylo, my.traits[my.phylo$tip.label, 1], method = "K", test = FALSE)

phylosig(my.phylo, my.traits[my.phylo$tip.label, 1], method = "K", test = TRUE, nsim = 1000)
phylosig(my.phylo, my.traits[my.phylo$tip.label, 2], method = "K", test = TRUE, nsim = 1000)
phylosig(my.phylo, my.traits[my.phylo$tip.label, 3], method = "K", test = TRUE, nsim = 1000)

phylosignal(my.traits[my.phylo$tip.label, 1],my.phylo)
phylosignal(my.traits[my.phylo$tip.label, 2],my.phylo)
phylosignal(my.traits[my.phylo$tip.label, 3],my.phylo)


phylosig(my.phylo, my.traits[my.phylo$tip.label, 1], method = "K", test = TRUE, nsim = 1000)
phylosignal(my.traits[my.phylo$tip.label, 1],my.phylo)

##########################################
##########################################
###########  Signal Pagel's Lambda Test ##
##########################################
##########################################
##########################################



lambda.phylo <- rescale(my.phylo, "lambda")

lambda.phylo(0.22)

par(mfrow = c(1,4))

plot(lambda.phylo(0))

title("Lambda = 0.00")

plot(lambda.phylo(0.25))

title("Lambda = 0.25")

plot(lambda.phylo(0.5))

title("Lambda = 0.50")

plot(lambda.phylo(1))

title("Lambda = 1.00")

phylosig(my.phylo, my.traits[my.phylo$tip.label, 1], method = "lambda", test = FALSE)

phylosig(my.phylo, my.traits[my.phylo$tip.label, 1], method = "lambda", test = TRUE)


##########################################
##########################################
###########  Phylogenetic Eigenvectors  ##
##########################################
##########################################
##########################################
##not supported on 3.3.2 yet apparently! ugh
#install.packages("PVR")
#library(PVR)

phylo.pca <- princomp(cophenetic(my.phylo))

pevs <- phylo.pca$scores[ , 1:3]
summary(lm(my.traits[row.names(pevs), 1 ] ~ pevs))


###PVR package code
#xx = PVRdecomp(my.phylo)
#yy = PSR(xx,my.traits[my.phylo$tip.label,])
#plot(yy)



table.phylog(as.data.frame(pevs),newick2phylog(write.tree(my.phylo)))
table.phylog(as.data.frame(cbind(pevs,my.traits[,1])),newick2phylog(write.tree(my.phylo)))


##########################################
##########################################
###########  Binary Traits    ############
##########################################
##########################################
##########################################
install.packages("caper")
library(caper)
binary.traits = as.data.frame(cbind( my.phylo$tip.label,c(1,1,1,1,1,1,0,0,0,0)))
colnames(binary.traits) = c("species", "trait")

phylo.d(binary.traits,my.phylo, names.col=species,binvar=trait)

####################################################################################
####################################################################################
####################################################################################
####################################################################################
##########################################
##########################################
###########  Comparing Models ############
##########################################
##########################################
##########################################
tr = my.traits[,1]
names(tr) = rownames(my.traits)

fitContinuous(my.phylo,tr,model="BM")
fitContinuous(my.phylo,tr,model="OU")

tmp.tra = sim.char(my.phylo,vcv(my.phylo),nsim=1,model="BM")[,1,]
fitContinuous(my.phylo,tmp.tra,model="BM")
fitContinuous(my.phylo,tmp.tra,model="OU")



####################################################################################
####################################################################################
####################################################################################
####################################################################################
##########################################
##########################################
##########################################
##########################################
###########  Node Signal with PIC ########
##########################################
##########################################
##########################################

pic.shuffle <- function(x){

## Shuffle the tip labels on the phylo.
x$tip.label <- sample(x$tip.label, length(x$tip.label), replace = FALSE)

## Calculate mean PIC with shuffled names.
pic(my.traits[x$tip.label, 1], x)
	}

#run null model of node pics
node.null <- replicate(999, pic.shuffle(my.phylo))

##combine obs and nulls
combined.obs.null <- cbind(pic(my.traits[my.phylo$tip.label, 1], my.phylo), node.null)

##calculate ranks 
apply(combined.obs.null, MARGIN = 1, rank)[1, ]

# put it in an object and plot it
my.ranks <- apply(combined.obs.null, MARGIN = 1, rank)[1, ]
plot(my.phylo)
nodelabels(my.ranks)

# get node ages/times
times <- branching.times(my.phylo)

## plot the ranks through time
plot(-1 * times, my.ranks, pch = 16)  


##########################################
##########################################
###########  Node Signal with DTT ########
##########################################
##########################################
##########################################

#calculate disparity through time with average squared distance and simulate 999 bm outcomes/nulls
my.dtt <- dtt(my.phylo, my.traits, index ="avg.sq", nsim = 999)

names(my.dtt)

my.dtt$dtt

my.dtt$times

my.dtt$sim[1:4,1:4]

apply(cbind(my.dtt$dtt, my.dtt$sim), MARGIN = 1, rank)[1, ]
