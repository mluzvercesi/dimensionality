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
enlacescoms <- apply(enlaces,2,function(x){membMCL[x]})

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

hist(edgesim, freq=FALSE, main= "similitud entre nodos enlazados", xlim=c(0.5,1), ylim=c(0,12), col=rgb(0,0,0.5,0.5))
hist(bordesim, freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))

#por comunidad
bordescoms <- edgecoms[edgecoms[,1]!=edgecoms[,2],]
borde_por_com <- rep( list(list()), length(unique(bordescoms[,1])))
for (c in sort(unique(bordescoms[,1]))){
  idx <- which(apply(bordescoms,1,function(x){if (x[1]==c|x[2]==c){1} else {0}})==1)
  borde_por_com[[c]] <- apply(bordes[idx,], 1, function(x){matriz[x[1],x[2]]})
  hist(edgesim, freq=FALSE, main= paste("todos los bordes vs bordes de comunidad \n Comunidad",c, "tamaño", table(membMCL)[c]),
       xlim=c(0.5,1), ylim=c(0,12), col=rgb(0,0,0.5,0.5))
  hist(borde_por_com[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
}

#comparar enlaces de borde vs internos
interiorcoms <- edgecoms[edgecoms[,1]==edgecoms[,2],]
int_por_com <- rep( list(list()), length(unique(interiorcoms[,1])))
for (c in sort(unique(interiorcoms[,1]))){
  idx <- which(interiorcoms[,1]==c)
  int_por_com[[c]] <- apply(interior[idx,], 1, function(x){matriz[x[1],x[2]]})
  hist(int_por_com[[c]], freq=FALSE, main=paste("enlaces internos vs bordes por comunidad \n Comunidad",c, "tamaño", table(membMCL)[c]),
       xlim=c(0.5,1), ylim=c(0,12), col=rgb(0,0,0.5,0.5))
  hist(borde_por_com[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
}

#comparacion
options(bitmapType="cairo")
for (c in sort(unique(bordescoms[,1]))){
  png(paste0("fig/coms_mcl_",c,".png"), width = 988, height = 494)
  par(mfrow=c(1,2))
  hist(edgesim, freq=FALSE, main= paste("todos los bordes vs bordes de comunidad \n Comunidad",c, "tamaño", table(membMCL)[c]),
       xlim=c(0.5,1), ylim=c(0,12), col=rgb(0.8,0.1,0,0.5))
  hist(borde_por_com[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
  hist(int_por_com[[c]], freq=FALSE, main=paste("enlaces internos vs bordes por comunidad \n Comunidad",c, "tamaño", table(membMCL)[c]),
       xlim=c(0.5,1), ylim=c(0,12), col=rgb(0,0,0.5,0.5))
  hist(borde_por_com[[c]], freq=FALSE, add=TRUE, col=rgb(0,0.6,0.5,0.5))
  dev.off()
}


boxplot.default(c(borde_por_com[1:13], int_por_com[1:13]), at=c(seq(1,26,2),seq(2,26,2)), col=c(rep("orange",13),rep("cyan",13)))
legend('bottomright', legend = c('borde','interior'), col = c('orange', 'cyan'), pch=15)


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

#-----
#grafo de comunidades como nodos
b <- apply(enlacescoms, 1, sort)
g <- make_graph(as.vector(b), directed = FALSE)

E(g)$weight <- 1
g2 <- simplify(g, edge.attr.comb=list(weight="sum"), remove.loops = FALSE)

nodesizes <- as.integer(table(membMCL)[1:vcount(g2)])
g2 <- set_vertex_attr(g2, "size", value=2*log(nodesizes))

#veo las comunidades con loops
E(g2)[which_loop(g2)]
loopweights <- get.edge.attribute(g2, name = "weight", index = which(which_loop(g2)))

gMCL <- simplify(g2)
plot(gMCL)

#COMPARAR SIMILITUD ENTRE COMUNIDADES "LEJANAS"



#------------------------------------------------------------------------

#para participacion y z-score (calculate_toproles.R)
membership <- membMCL[lccnm]
g_users <- knn.lcc

nodos_perif <- unique(as.list(enlaces[enlacescoms[,1]!=enlacescoms[,2],]))
n <- unique(as.list(enlaces[enlacescoms[,1]==enlacescoms[,2],]))
nodos_centr <- n[!(n %in% nodos_perif)]
nodos_centr <- do.call(c, nodos_centr)
nodos_perif <- do.call(c, nodos_perif)

pz <- list(names(connector_hubs), names(connector_nonhubs), names(provincial_hubs), names(provincial_nonhubs))
names(pz) <- c("con_h", "con_nh", "prov_h", "prov_nh")
a <- reverseSplit(pz)
a <- unlist(a)
centralidad <- list(nodos_centr,nodos_perif)
names(centralidad) <- c("centr", "perif")
b <- reverseSplit(centralidad)
b <- unlist(b)
table(a,b)


membership <- celltypes_nro[lccnm]

enlacesct <- apply(enlaces,2,function(x){celltypes[x]})
nodos_perif2 <- unique(as.list(enlaces[enlacesct[,1]!=enlacesct[,2],]))
n2 <- unique(as.list(enlaces[enlacesct[,1]==enlacesct[,2],]))
nodos_centr2 <- n2[!(n2 %in% nodos_perif2)]
nodos_centr2 <- do.call(c, nodos_centr2)
nodos_perif2 <- do.call(c, nodos_perif2)

pz2 <- list(names(connector_hubs), names(connector_nonhubs), names(provincial_hubs), names(provincial_nonhubs))
names(pz2) <- c("con_h", "con_nh", "prov_h", "prov_nh")
c <- reverseSplit(pz2)
c <- unlist(c)
centralidad2 <- list(nodos_centr2,nodos_perif2)
names(centralidad2) <- c("centr", "perif")
d <- reverseSplit(centralidad2)
d <- unlist(d)
table(c,d)


plot(p_i1[nodos_perif], z_i1[nodos_perif], col="blue",main="Comunidades MCL", 
     xlim=c(min(p_i),max(p_i)), ylim = c(min(z_i), max(z_i)), xlab="participacion", ylab="grado en comunidad")
points(p_i1[nodos_centr], z_i1[nodos_centr], col="cyan")
abline(v=mean(p_i1),lwd=2, lty=2,col='green')
abline(h=0,lwd=2, lty=2,col='green')
legend('topright', legend=c('perifericos','centrales'),col=c('blue','cyan'),pch=c(1,1))

plot(p_i2[nodos_perif2], z_i2[nodos_perif2], col="blue",main="Comunidades celltype", 
     xlim=c(min(p_i2),max(p_i2)), ylim = c(min(z_i2), max(z_i2)), xlab="participacion", ylab="grado en comunidad")
points(p_i2[nodos_centr2], z_i2[nodos_centr2], col="cyan")
abline(v=mean(p_i2),lwd=2, lty=2,col='green')
abline(h=0,lwd=2, lty=2,col='green')
legend('topright', legend=c('perifericos','centrales'),col=c('blue','cyan'),pch=c(1,1))


plot(kcoreness, p_i)