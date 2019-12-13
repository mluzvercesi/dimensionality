#Analisis de comunidades: participacion de los nodos, coreness de enlaces, bordes de clusters

# Asortatividad vectorial por comunidad
matriz <- cos_sim(nes0)

for (i in 1:length(comms)){
  nombres <- comms[[i]]
  nombres <- nombres[nombres %in% colnames(matriz)]
  red <- induced_subgraph(knn, nombres)
  r <- assortativity_vect(red, matriz[nombres,nombres])
  cat("r = ", r, "Comunidad", i, "tamaño", table(membMCL)[i], "\n")
}

for (i in 1:length(unique(celltypes))){
  nombres <- names(which(celltypes[lccnm]==unique(celltypes)[i]))
  nombres <- nombres[nombres %in% colnames(matriz)]
  red <- induced_subgraph(knn, nombres)
  r <- assortativity_vect(red, matriz[nombres,nombres])
  cat("r = ", r, "Comunidad", unique(celltypes)[i], "tamaño", table(celltypes[lccnm])[i],"\n")
}


enlaces <- as_edgelist(knn.lcc, names = TRUE)

nodos0 <- colnames(matriz)[colnames(matriz) %in% vertex_attr(knn.lcc, name = "name")] #nodos en nes0 y en el lcc
knn.lcc0 <- induced_subgraph(knn.lcc, nodos0)

edges <- as_edgelist(knn.lcc0, names = TRUE)
edgesim <- apply(edges, 1, function(x){matriz[x[1],x[2]]})
edgecoms <- apply(edges,2,function(x){membMCL[x]})

bordes <- edges[edgecoms[,1]!=edgecoms[,2],]
bordesim <- apply(bordes, 1, function(x){matriz[x[1],x[2]]})

interior <- edges[edgecoms[,1]==edgecoms[,2],]
interiorsim <- apply(interior, 1, function(x){matriz[x[1],x[2]]})


hist(matriz[upper.tri(matriz, diag=FALSE)], freq=FALSE, main="similitud en toda la red", xlim=c(0.5,1), ylim=c(0,12))
hist(bordesim, freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))

hist(edgesim, freq=FALSE, main= "similitud entre nodos enlazados", xlim=c(0.5,1), ylim=c(0,12))
hist(bordesim, freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))

#por comunidad
bordescoms <- edgecoms[edgecoms[,1]!=edgecoms[,2],]
borde_por_com <- rep( list(list()), length(unique(bordescoms[,1])))
for (c in sort(unique(bordescoms[,1]))){
  idx <- which(apply(bordescoms,1,function(x){if (x[1]==c|x[2]==c){1} else {0}})==1)
  borde_por_com[[c]] <- apply(bordes[idx,], 1, function(x){matriz[x[1],x[2]]})
  hist(edgesim, freq=FALSE, main= paste("todos los bordes vs bordes de comunidad \n Comunidad",c, "tamaño", table(membMCL)[c]),
       xlim=c(0.5,1), ylim=c(0,12))
  hist(borde_por_com[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
}

#comparar enlaces de borde vs internos
interiorcoms <- edgecoms[edgecoms[,1]==edgecoms[,2],]
int_por_com <- rep( list(list()), length(unique(interiorcoms[,1])))
for (c in sort(unique(interiorcoms[,1]))){
  idx <- which(interiorcoms[,1]==c)
  int_por_com[[c]] <- apply(interior[idx,], 1, function(x){matriz[x[1],x[2]]})
  hist(int_por_com[[c]], freq=FALSE, main=paste("enlaces internos vs bordes por comunidad \n Comunidad",c, "tamaño", table(membMCL)[c]),
       xlim=c(0.5,1), ylim=c(0,12))
  hist(borde_por_com[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
}



# atributo para enlaces: vecinos en comun
frac_comun <- apply(enlaces, 1, function(x){vecinoscomun(knn.lcc,x) })

#Centralidad de los enlaces
kcoreness <- coreness(knn.lcc)
enlaces_k <- enlaces
enlaces_k[,1] <- kcoreness[enlaces[,1]]
enlaces_k[,2] <- kcoreness[enlaces[,2]]
enlaces_k <- matrix(as.integer(enlaces_k), nrow=dim(enlaces_k)[1], ncol=dim(enlaces_k)[2])

edge_core <- apply(enlaces_k, 1, min)

sim.nes <- cos_sim(nes0)
knn.nes <- delete_vertices(knn.lcc, v=colnames(cells.nes)[is.nan(apply(cells.nes, 2, min))])
e <- as_edgelist(knn.nes, names = TRUE)
edge_sim <- apply(e, 1, function(x){ sim.nes[x[1], x[2]] })

mean_edge_sim <- rep(0, length(unique(edge_core)))
for (i in 1:length(unique(edge_core))){
  mean_edge_sim[i] <- edge_sim[edge_core=unique(edge_core)[i]]
}

par(mfrow=c(1,1))
plot(edge_core, edge_sim)
points(unique(edge_core), mean_edge_sim, col="green", pch=16)

corelist <- sort(unique(edge_core), decreasing = TRUE)

g <- graph_from_edgelist(enlaces[edge_core>corelist[7],], directed = FALSE)
is.connected(g)
alphavec <- kcoreness[V(g)]/max(kcoreness[V(g)])
plot.igraph(g, vertex.label=NA, vertex.size=5, edge.width=0.5,
            vertex.color = rgb(0,1,0,alphavec))

#------------------------------------------------------------------------

#para participacion y z-score (calculate_toproles.R)
membership <- membMCL[lccnm]
g_users <- knn.lcc

plot(kcoreness, p_i)