#Analisis de comunidades: similitud, participacion de los nodos, coreness de enlaces, bordes de clusters

Ngen <- dim(dataAsub)[1]
enlaces <- as_edgelist(knn.lcc)
#cuantos genes tienen 0 en ambos nodos del enlace
a <- apply(enlaces, 1, function(x){sum((dataAsub[,x[1]]+dataAsub[,x[2]])==0)})
dropoutperc <- a/Ngen

similitud <- cos_sim(dataAsub)
similitudpares <- apply(enlaces, 1, function(x){similitud[x[1],x[2]]})
par(mfrow=c(1,1))
plot(similitudpares, dropoutperc)

intracom <- membMCL[enlaces[,1]]
idx <- apply(enlaces, 1, function(x){membMCL[x[1]]!=membMCL[x[2]]})
intracom[idx] <- 0
intracom[idx] <- max(intracom)+1
colores <- rainbow(length(unique(intracom))-1)
colores <- c(colores, 'black')

plot(similitudpares, dropoutperc, col=colores[intracom], 
     main="Similitud entre pares de nodos enlazados vs \n porcentaje de dropouts en ambos nodos del enlace", 
     xlab="similitud", ylab="porcentaje de dropouts")
legend('bottomleft', legend=c(seq(1:18),'dif'), col=colores, pch=1, cex=0.7)

#por comunidad
xs <- par("usr")[1:2]
ys <- par("usr")[3:4]

par(mfrow=c(2,3))
for (i in 1:12){
  idx1 <- apply(enlaces, 1, function(x){if (membMCL[x[1]]==i | membMCL[x[2]]==i){TRUE} else {FALSE}})
  plot(similitudpares[idx1], dropoutperc[idx1], col=colores[intracom[idx1]], main=paste('Comunidad',i), xlim=xs, ylim=ys,
       xlab="similitud", ylab="porcentaje dropout")
}


#---- REVISAR A PARTIR DE ACA
#distribucion de similitud entre pares por comunidad (comparo con azar pero sin incluir la propia comunidad)
ady <- as_adjacency_matrix(knn, sparse = FALSE)
matriz <- cos_sim(nes0)
N <- 100
par(mfrow=c(3,2))

for (i in 1:6){
  icomm <- which(membMCL[lccnm]==i)
  idx <- colnames(matriz) %in% names(membMCL[lccnm])[icomm]
  A <- ady[idx,idx]
  enlacesvect <- A[upper.tri(A)]
  enlacesvect <- ifelse(enlacesvect, TRUE, FALSE)
  X <- matriz[idx,idx]
  similvect <- X[upper.tri(X)]
  
  hist(similvect[enlacesvect], freq=FALSE,
       main=paste("C",i,"tama�o",table(membMCL[lccnm])[i]), xlab = "", ylim=c(0,10))
  
  npar <- (length(X)-dim(X)[1])/2
  
  X <- matriz[!idx,!idx]
  pares <- X[upper.tri(X, diag = FALSE)]
  A <- matrix(0, nrow = N, ncol = npar)
  for (j in 1:dim(A)[1]){
    A[j,] <- sample(pares, npar)
  }
  hist(A, freq=FALSE, col=rgb(1,0,0,0.5), add=TRUE)
}


#----
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