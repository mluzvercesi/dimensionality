par(mfrow=c(2,2)) 
for (i in 1:4){
  orderedgenes <- genes_by_weight(pca,ncomp=i,retw = TRUE)
  plot(orderedgenes$weight, xlab="# gen", ylab="Peso", main=paste("PC",i))
  text(100, 0.004, round(sum(orderedgenes$weight[1:100]),digits=2),col="red")
  text(500, 0.003, round(sum(orderedgenes$weight[1:500]),digits=2),col="red")
  text(1000, 0.002, round(sum(orderedgenes$weight[1:1000]),digits=2),col="red")
}

for (i in 5:8){
  orderedgenes <- genes_by_weight(pca,ncomp=i,retw = TRUE)
  plot(orderedgenes$weight, xlab="# gen", ylab="Peso", main=paste("PC",i))
  text(100, 0.004, round(sum(orderedgenes$weight[1:100]),digits=2),col="red")
  text(500, 0.003, round(sum(orderedgenes$weight[1:500]),digits=2),col="red")
  text(1000, 0.002, round(sum(orderedgenes$weight[1:1000]),digits=2),col="red")
}

for (i in 97:100){
  orderedgenes <- genes_by_weight(pca,ncomp=i,retw = TRUE)
  plot(orderedgenes$weight, xlab="# gen", ylab="Peso", main=paste("PC",i))
  text(100, 0.004, round(sum(orderedgenes$weight[1:100]),digits=2),col="red")
  text(500, 0.003, round(sum(orderedgenes$weight[1:500]),digits=2),col="red")
  text(1000, 0.002, round(sum(orderedgenes$weight[1:1000]),digits=2),col="red")
}