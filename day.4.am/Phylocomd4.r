################################
########### IN #################
#########   R      #############
################################
################################
################################

setwd("/Users/nateswenson/sesync.workshop/data")

library(picante)

#### DATA

my.phylo = read.tree("toy.phylo.txt")
my.sample = readsample("toy.cdm.txt")
my.traits = read.csv("toy.traits.csv",row.names=1)

setwd("/Users/nateswenson/phylocom-4.2/mac")


write.tree(my.phylo,"my.phylo.1.txt")
writesample(my.sample,"my.sample.1.txt")


################################
########### IN #################
######### PHYLOCOM #############
################################
################################
################################
./phylocom

./phylocom pd -f my.phylo.1.txt -s my.sample.1.txt

./phylocom comstruct -f my.phylo.1.txt -s my.sample.1.txt

./phylocom comdist -f my.phylo.1.txt -s my.sample.1.txt -n -m 0 -a 


################################
########### IN #################
#########   R      #############
################################
################################
################################

system("./phylocom pd -f my.phylo.1.txt -s my.sample.1.txt > output.txt")

################################
################################
################################
################################
################################
# Run something iteratively using both R and C program
# so...let's try and randomize a phylogeny in R
# put it in working directory (same wd as where your phylocom is) with a unique name
# run the pd command in phylocom and output a file uniquely named
# read output back into R & plot it

# telling R to get ready to plot 10 panels in 2 rows and 5 columns
par(mfrow=c(2,5))

# loop through 1 to 10
for(i in 1:10){

# shuffle names in phylo
tmp.tree = tipShuffle(my.phylo)

# write out tree to working directory using paste to put the number of the 
# iteration on the front of the file name to make it unique
write.tree(tmp.tree, paste(i,".tree.txt",sep=""))

# use system command to call phylocom and use past to call correct phylo and output unique file
system(paste("./phylocom pd -f ", i, ".tree.txt -s my.sample.1.txt > ", i, ".output.txt",sep=""))

# read in output
tmp.output = read.table(paste(i,".output.txt",sep=""), sep="\t", header=T)

# plot pd onto species richness
plot(tmp.output[,3] ~ tmp.output[,2],cex=20,pch=20, xlab="Species Richness", ylab="Faith's PD")

}


################################
################################
################################
################################
################################
## the goal here is to take a trait file and put it into phylocom format

## we are working with a 3 trait matirx with all continuous traits.
## so we put a new row with 3's to indicate continuous traits
## and a second row with names of traits
## remaining rows will be the trait data
my.traits2 = rbind(c(3,3,3), colnames(my.traits), my.traits)
my.traits2

## add names for the first two rows using phylocom format
rownames(my.traits2)[1:2]= c("type","name")
my.traits2

## write out file to wd
write.table(my.traits2, "my.traits4phylocom.txt",sep="\t",row.names=T,col.names=F,quote=F)

################################
########### IN #################
######### PHYLOCOM #############
################################
################################
################################
######### 1 = var ############## 
######### 2 = mean #############
######### 3 = mpd ##############
######### 4 = mntd #############
################################
################################

./phylocom comtrait -t my.traits4phylocom.txt -s my.sample.1.txt -x 1 > trait.output.txt


./phylocom aot -f my.phylo.1.txt -t my.traits4phylocom.txt

./phylocom nodesigl -f my.phylo.1.txt -s my.sample.1.txt > nodesig.output.txt



################################
########### IN #################
#########   R      #############
################################
################################
################################
setwd("/Users/nateswenson/sesync.workshop/data")
my.tax = read.csv("my.species.list.csv")


new.tax = cbind(my.tax[,1:2], paste(my.tax[,2], my.tax[,3],sep="_"))
new.tax

setwd("/Users/nateswenson/phylocom-4.2/mac")
write.table(new.tax, "my.taxa.txt",sep="/",row.names=F,col.names=F,quote=F)


system("./phylomatic -t my.taxa.txt -f R20100701.new")

new.tax2 = cbind(tolower(my.tax[,1]), tolower(my.tax[,2]), paste(tolower(my.tax[,2]), tolower(my.tax[,3]),sep="_"))
write.table(new.tax2, "my.taxa.2.txt",sep="/",row.names=F,col.names=F,quote=F)

system("./phylomatic -t my.taxa.2.txt -f R20100701.new")

new.tax2[6,1]="urticaceae"

write.table(new.tax2, "my.taxa.3.txt",sep="/",row.names=F,col.names=F,quote=F)

system("./phylomatic -t my.taxa.3.txt -f R20100701.new > phylomatic.output.txt")

phylomatic.tree = read.tree("phylomatic.output.txt")
plot(phylomatic.tree)

./phylocom cleanphy -f phylomatic.output.txt -e > phylomatic.output2.txt

./phylocom bladj -f phylomatic.output2.txt