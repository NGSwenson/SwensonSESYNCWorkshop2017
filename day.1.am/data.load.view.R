setwd("/Users/nateswenson/sesync.workshop/data")
#Load packages
library(ape)
library(phytools)
library(picante)
library(ade4)
library(gplots)

####PHYLOGENETIC DATA
#load phylo data from lefse package
my.phylo <- read.tree("toy.phylo.txt")

#see what is in the object
names(my.phylo)

#what are the names of the tips (species) on the phylo
my.phylo$tip.label

#generically plot the phylogeny
plot(my.phylo)

#ask if the phylogeny is rooted
is.rooted(my.phylo)

#as if the phylogeny is ultrametric
is.ultrametric(my.phylo)

#look at the edge/branch descriptors
my.phylo$edge

#place node lables on plotted phylo
nodelabels()

#print out the lengths of each branch
my.phylo$edge.length

#sum up all branch lengths. This is called tree length.
tree.length = sum(my.phylo$edge.length)

#look at the length and nodes for each branch
cbind(my.phylo$edge,my.phylo$edge.length)

#plot the phylogeny again
plot(my.phylo)

#label just node 12 with a word.
nodelabels("howdy",12, frame="circle",bg="yellow")

#make one phylo object per clade in the big phylogeny
subs = subtrees(my.phylo)

#see how many clades. Should be equivalent to number of internal nodes
length(subs)

#plot a subtree
plot(subs[[1]])

#plot a subtree
plot(subs[[2]])

#get the time/date/relative time of each internal node
branching.times(my.phylo)

#plot the phylogeny as a fan/circular type tree
plot.phylo(my.phylo, type="fan",edge.color="red",show.tip.label=F)

#make a temporary phylo to manipulate
my.tmp.phylo = my.phylo

#shuffle names in the temporary phylogeny object
my.tmp.phylo$tip.label = sample(my.phylo$tip.label,length(my.phylo$tip.label), replace=F)

#plot the randomized phylogeny
plot(my.tmp.phylo)

#use the tipShuffle() function to randomize names and plot
plot(tipShuffle(my.tmp.phylo))

rand.names = function(ph){

ph$tip.label = sample(ph$tip.label,length(ph$tip.label), replace=F)

ph

}


## make 6 random phylos and plot next to each other
par(mfrow=c(2,3))

for(i in 1:6){
	plot(rand.names(my.tmp.phylo),cex=10)
	}

my.rand.phylos = list(1:6)
	
for(i in 1:6){
write.tree(rand.names(my.tmp.phylo), paste(i,".new.txt",sep=""))
	}


my.tree.lengths = matrix(NA, ncol=6,nrow=1)

for(i in 1:6){
my.tree.lengths[1,i] = sum(rand.names(my.tmp.phylo)$edge.length)
	}




#plot original phylo
plot(my.phylo)

#make a new phylo object that has pruned out species 1 and 3
pruned.phylo = drop.tip(my.phylo, c("t13","t1"))

#plot the pruned phylo
plot(pruned.phylo)

#randomly add a new extant species to your phylo
plot(add.random(my.phylo))

#make a temporary phylo
tmp.phylo = my.phylo

#add one species randomly 6 times until you have 6 new species
par(mfrow=c(2,3))

for(i in 1:6){
	tmp.phylo = add.random(tmp.phylo)
	plot(tmp.phylo)
		}

#add a new species in a specific location
plot(bind.tip(my.phylo,"sp11",where=12))

##make a clade in your phylogeny a star phylogeny 
plot(collapse.to.star(my.phylo,16))

#plot phylo and look at node labels
plot(my.phylo,cex=1-)
nodelabels()

#find the sister node to node 12
getSisters(my.phylo,12)

#get internal and terminal (tips) nodes descendant from node 12
getDescendants(my.phylo,12)

# get only terminal nodes
subset(getDescendants(my.phylo,12),getDescendants(my.phylo,12)<11)

sister = getSisters(my.phylo,12)
sister.div = length(subset(getDescendants(my.phylo,sister),getDescendants(my.phylo,sister)<11))
focal.div = length(subset(getDescendants(my.phylo,12),getDescendants(my.phylo,12)<11))


#make a fancy plot of your phylogeny
fancyTree(my.phylo)

#make a phylogenetic distance matrix
cophenetic(my.phylo)

#make a phylogenetic variance-covariance matrix
vcv(my.phylo)

#look at the branching times
branching.times(my.phylo)

##plot node labes
nodelabels()

#store branching times
my.times=branching.times(my.phylo)



#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
####COMMUNITY DATA

#get the community data matirx example
my.sample <- readsample("toy.cdm.txt")

# use the vegan function decostand() to divide all abundances by the toal abundance in community
my.sample.new = decostand(my.sample, method="total",MARGIN=1)

# use the vegan function decostand() to make your community data matrix into presence-absence data
decostand(my.sample.new, method="pa",MARGIN=1)

# calculate community species richness
rowSums(decostand(my.sample, method="pa"))

# calculate species occupancy rates in the system
colSums(decostand(my.sample, method="pa"))

# order community data matix by names on the phylogeny
my.sample[, my.phylo$tip.label]

# for the first community extract data for only species that are present
my.sample[1,my.sample[1,]>0]

# for the first community extract names of only species that are present
names(my.sample[1,my.sample[1,]>0])

#make a small function that could get the names of present species for a community (a vector representing one row in a community data matrix)
tmp.fun = function(my.sample){
	names(my.sample[my.sample>0])
	}

# apply the function to all rows in the community data matrix
present.names = apply(my.sample,1, tmp.fun)
present.names

#plot the phylogeny without the species found in community 1
plot(drop.tip(my.phylo,my.phylo$tip.label[match(present.names$com.1, my.phylo$tip.label)]))

#plot the phylogeny without the species found in community 3
plot(drop.tip(my.phylo,my.phylo$tip.label[match(present.names$com.13, my.phylo$tip.label)]))

sam = matrix2sample(my.sample)
sam
sample2matrix(sam)
writesample(my.sample,"phylocom.sample.txt")
new.samp = readsample("phylocom.sample.txt")
new.samp


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
####Traits data

#read in trait data from lefse package
trait.dat =read.csv("toy.traits.csv",row.names=1)

#look at trait data
trait.dat

#plot histograms of the data
par(mfrow=c(1,3))
for(i in 1:3){hist(trait.dat[,i])}

#plot all traits against each other
plot(trait.dat)

#summarize data
apply(trait.dat,2, mean,na.rm=T)
apply(trait.dat,2, sd,na.rm=T)
apply(trait.dat,2, max,na.rm=T)
apply(trait.dat,2, min,na.rm=T)

#z transform/scale data
scale(trait.dat)

#take out one trait and store as vectore
new.traits = trait.dat[,1]

#place species names on the values in the vector
names(new.traits)=rownames(trait.dat)

#take a look
new.traits

#plot a traitgram
traitgram(new.traits[my.phylo$tip.label],my.phylo)


#make new trait matrix with a column containing species names
new.traits = cbind(rownames(trait.dat),trait.dat)

#take a look
new.traits

#name the new column
names(new.traits)[1]="species"

#plot the phylo with species names colored by trait values for trait 1
color.plot.phylo(my.phylo, new.traits, "trait1","species")

#plot the phylo with species names colored by trait values for trait 1 but binning by 3 breaks
color.plot.phylo(my.phylo, new.traits, "trait1","species",num.breaks=3)

#simulate some traits
Y<- sim.corrs(my.phylo,vcv=matrix(c(1,0.75,0.75,1),2,2) )

#plot a 3d traitgram
fancyTree(my.phylo,type="traitgram3d",X=Y)

##make a new phylo object where all species from node 13 have the state "africa"
new.tree = paintSubTree(my.phylo, 13, "Africa", stem=FALSE)

#define a vector of colors that will match our states
cols=c("black","red")

#name the states in the colors matrix
names(cols)=c(1,"Africa")

#plot the new phylogeny with edges colored by states
plotSimmap(new.tree,cols)

#look at the node labels
nodelabels()

#id the stem of a node as having the state "Africa"
aa = paintBranches(new.tree,15,"Africa")

#plot the updated phylogeny
plotSimmap(aa, cols)

#plot your trait data in 2 dimensions linked by phylogenetic branches
phylomorphospace(my.phylo, trait.dat[,1:2],xlab="trait 1",ylab="trait 2")

#plot your trait data in 3 dimensions linked by phylogenetic branches
phylomorphospace3d(my.phylo, trait.dat[,1:3])


#new.phylo = phylo4d(my.phylo,trait.dat[,1:3])
#table.phylo4d(new.phylo)

#make a newick object by copying and pasting text from phylo txt file on hard drive surrounded by quotes
phylog.tre = "((((sp10:0.05647302435,sp5:0.05647302435):0.0549812367,sp8:0.1114542611):0.215108127,(sp2:0.1561272726,(sp7:0.07400684377,sp9:0.07400684377):0.08212042887):0.1704351154):0.8062464292,(sp1:0.1798632653,((sp6:0.05787762397,sp3:0.05787762397):0.001076379461,sp4:0.05895400343):0.1209092618):0.952945552);"

#transform newick to a phylog object
phylog.tre = newick2phylog(phylog.tre)

#plot the phylog
plot(phylog.tre)

#take a look at whats in the object
phylog.tre

 phylog.tre$nodes
 
 phylog.tre$Wdist
 
 phylog.tre$Wmat
 
 vcv(my.phylo)

#change a phylog object back to a phylo object 
 toPhylo = as.phylo(phylog.tre)
 
 #make a trait dendrogram and transform to a phylo object
 ff = as.phylo(hclust(dist(trait.dat)))

 ##make an association matrix
 matt = cbind(sort(my.phylo$tip.label), sort(ff$tip.label))
 
 #plot it
cophyloplot(my.phylo,ff,assoc=matt,  length.line = 4, space = 28, gap = 3)
  
 ##make an association matrix
 matt = cbind(my.phylo$tip.label, ff$tip.label)
 
 #plot it
 cophyloplot(my.phylo,ff,assoc=matt,  length.line = 4, space = 28, gap = 3)
 
 
 #heatmap
 heatmap.2(as.matrix(trait.dat[,1:3]),Rowv=as.dendrogram(as.hclust(my.phylo)), trace="none",tracecol="black", margins=c(10,5),col="heat.colors")
 
 
 
 
 