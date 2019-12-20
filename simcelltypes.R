edgecoms2 <- apply(edges,2,function(x){celltypes[x]})

bordes2 <- edges[edgecoms2[,1]!=edgecoms2[,2],]
bordesim2 <- apply(bordes2, 1, function(x){matriz[x[1],x[2]]})

interior2 <- edges[edgecoms2[,1]==edgecoms2[,2],]
interiorsim2 <- apply(interior2, 1, function(x){matriz[x[1],x[2]]})


hist(matriz[upper.tri(matriz, diag=FALSE)], freq=FALSE, main="similitud en toda la red",sub="celltype", xlim=c(0.5,1), ylim=c(0,12))
hist(bordesim2, freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))

hist(edgesim, freq=FALSE, main= "similitud entre nodos enlazados", sub="celltype", xlim=c(0.5,1), ylim=c(0,12))
hist(bordesim2, freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))

#por comunidad
bordescoms2 <- edgecoms2[edgecoms2[,1]!=edgecoms2[,2],]
borde_por_com2 <- rep( list(list()), length(unique(bordescoms2[,1])))
for (c in unique(bordescoms2[,1])){
  idx <- which(apply(bordescoms2,1,function(x){if (x[1]==c|x[2]==c){1} else {0}})==1)
  borde_por_com2[[c]] <- apply(bordes2[idx,], 1, function(x){matriz[x[1],x[2]]})
  hist(edgesim, freq=FALSE, main= paste("todos los bordes vs bordes de comunidad \n Comunidad",c, "tamaño", table(celltypes)[c]),
       xlim=c(0.5,1), ylim=c(0,12))
  hist(borde_por_com2[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
}

#comparar enlaces de borde vs internos
interiorcoms2 <- edgecoms2[edgecoms2[,1]==edgecoms2[,2],]
int_por_com2 <- rep( list(list()), length(unique(interiorcoms2[,1])))
for (c in unique(interiorcoms2[,1])){
  idx <- which(interiorcoms2[,1]==c)
  int_por_com2[[c]] <- apply(interior2[idx,], 1, function(x){matriz[x[1],x[2]]})
  hist(int_por_com2[[c]], freq=FALSE, main=paste("enlaces internos vs bordes por comunidad \n Comunidad",c, "tamaño", table(celltypes)[c]),
       xlim=c(0.5,1), ylim=c(0,12))
  hist(borde_por_com2[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
}
