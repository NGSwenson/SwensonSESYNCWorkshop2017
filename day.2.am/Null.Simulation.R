

tmp.tree = rcoal(100)
plot(tmp.tree)


outt=matrix(NA,nrow=98,ncol=100)

for(j in 1:100){
for(i in 1:98){
 
outt[i,j]  = sum(drop.tip(tmp.tree, sample(tmp.tree$tip.label, i))$edge.length)

}}

plot(rep(99:2,100), c(outt), xlab="Community Species Richness", ylab="Phylogenetic Diversity [Faith's Index]")




outt2 = matrix(NA, nrow=98, ncol=100)


for(j in 1:100){
for(i in 1:98){
 
outt2[i,j]  = mean(cophenetic(drop.tip(tmp.tree, sample(tmp.tree$tip.label, i))))

}}

plot(rep(99:2,100), c(outt2), xlab="Community Species Richness", ylab="Phylogenetic Diversity [Mean Pairwise Distance]")