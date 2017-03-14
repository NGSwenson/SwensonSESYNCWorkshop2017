lseries = rlnorm(100000, meanlog = 0, sdlog = 1)

hist(log(lseries),ylab="log(Number of Species)",xlab= "Genus" , main="", col="light grey", lwd=2)

hist(lseries,ylab="Number of Species",xlab= "Genus" , main="", col="light grey", lwd=2)

tt = findInterval(lseries, c(0:62))

zz = cbind(tt, lseries)


outt = matrix(NA, nrow=100, ncol=1000)

for(i in 1:100){print(i)
for(j in 1:1000){

spe = i
gen = length(unique(sample(zz[,1], i, replace = F)))

outt[i,j] = gen/spe

}
}

write.table(outt, "GS.sim.csv",sep=",", row.names=F,col.names=F,quote=F)
write.table(zz, "GS.tab.csv",sep=",", row.names=F,col.names=F,quote=F)


bb = matrix(NA, nrow=100, ncol=2);aa = for(i in 1:100){ bb[i,2]=sd(outt[i,]); bb[i,1]=mean(outt[i,])}
cc= bb[,1]-bb[,2]
dd= bb[,1]+bb[,2]
bbb=cbind(1:100,bb,cc,dd)

bbb[,5]= ifelse(bbb[,5]>1, 1,bbb[,5])

plot(rep.int(1:100, 1000), c(outt), xlab="Community Species Richness", ylab = "Genus:Species", pch=21,cex=.1, bg="black")
 polygon(c(1:100,100:1), c(bbb[1:100,4], bbb[100:1,5]), col=adjustcolor("grey", alpha=.5) )
 lines(1:100, bbb[,2], lty=2, lwd=3)