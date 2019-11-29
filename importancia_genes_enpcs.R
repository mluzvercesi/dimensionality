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

genes_pesos <- function(results,ncomp=1){ 
  # Devuelve una lista ordenada de nombres de genes segun su peso sobre los componentes principales ncomp de X
  # Si retw = TRUE devuelve el peso ("importancia") de cada gen
  # Necesita matriz de scexp (X) y resultado de pca_svd o de prcomp (results)
  p <- results$rotation[,ncomp]
  lambdas <- (results$sdev)^2
  w <- lambdas[ncomp]/sum(lambdas)
  if (is.null(dim(p))){ # si solo hay un componente principal
    p <- p*w
    weights <- abs(p)
  }else{ # si hay mas de uno, es matriz
    p <- apply(p,1,function(x){x*w})
    weights <- apply(p,2,function(x){sum(abs(x))})
  }
  #weights <- weights/sum(weights)
  # los pesos pueden variar entre prcomp y pca_svd, pero el orden es el mismo
  
  indices <- order(weights,decreasing = TRUE) # orden de maximas proyecciones
  
  #nombre <- list("genes"=rownames(results$rotation)[indices],"weights"=weights[indices])
  return(weights[indices])
}

par(mfrow=c(1,1)) 
orderedgenes1 <- genes_pesos(pca,ncomp=1)
plot(orderedgenes, xlab="# gen", ylab="Peso", main=paste("PC 1 vs 100"))
orderedgenes100 <- genes_pesos(pca,ncomp=100)
points(orderedgenes, col="blue")
