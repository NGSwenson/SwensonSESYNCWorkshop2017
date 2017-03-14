tree.mess = read.tree("cluster10_1rr_MIortho2.tre")
tree.mess
plot(tree.mess)

nodelabels()

nodelabels(tree.mess$node.label)

tree.mess$tip.label


gsub("_","",tree.mess$tip.label)

gsub("._","",tree.mess$tip.label)

gsub(".._","",tree.mess$tip.label)

gsub(".*_","",tree.mess$tip.label)

gsub("@.","",tree.mess$tip.label)

gsub("@.*","",tree.mess$tip.label)

gsub("@.*","_savin.time",tree.mess$tip.label)

tree.clean = tree.mess 
tree.clean$tip.label = gsub("@.*","_savin.time",tree.mess$tip.label)

plot(tree.clean)

################################################################################################
################################################################################################
################################################################################################
################################################################################################

examp = c("yes","we","can")
examp

examp = gsub("es","ou",examp)
examp = gsub("we","are",examp)
gsub("can","great",examp)

################################################################################################
################################################################################################
################################################################################################
################################################################################################


horrible.names = c("spp_gr8", "spp_10", "spp_11r","spp_1")

gsub("\\d","",horrible.names)

gsub("\\D","",horrible.names)

horrible.names
gsub("_\\d","",horrible.names)
gsub("_\\d*","",horrible.names)
gsub("_\\d*.","",horrible.names)


horrible.names
gsub("_\\d*","",horrible.names)
gsub("_\\d*\\D*","",horrible.names)



horrible.names2 = c("Spp_gr8", "spp_10", "spp_11r","spp_1")

gsub("Spp","ack",horrible.names2)
gsub("Spp","ack",horrible.names2,ignore.case=T)


gsub("[[:upper:]]","i.hate.capital.letters.",horrible.names2)



################################################################################################
################################################################################################
################################################################################################
################################################################################################
#GREP
new.tax2

grep("aceae",new.tax2[,1],value=F)

grep("aceae",new.tax2[,1],value=T)

grep("fab",new.tax2[,1],value=T)

grep("fab",new.tax2[,1],value=F)

new.tax2[grep("fab",new.tax2[,1],value=F),]

new.tax2[grep("fab|mel",new.tax2[,1]),]

new.tax3 = new.tax2

new.tax3[grep("fab|mel",new.tax3[,1]),1] = "FAMILYaceae"





















