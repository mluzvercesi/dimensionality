#Analisis de comunidades: similitud, participacion de los nodos, coreness de enlaces, bordes de clusters

#grafo simplificado de comunidades como nodos
gMCL <- grafosimple(knn, membMCL, plt='lcc')


#COMPARAR SIMILITUD ENTRE COMUNIDADES "LEJANAS" por ej: 1-3 1-6 2-6 3-4


# similitud vs dropouts----
Ngen <- dim(dataAsub)[1]
enlaces <- as_edgelist(lcc)
#cuantos genes tienen 0 en ambos nodos del enlace
X <- dataAsub==0
a <- apply(enlaces, 1, function(x){sum(X[,x[1]]*X[,x[2]])})
dropoutperc <- a/Ngen

similitud <- cos_sim(dataAsub)
similitudpares <- apply(enlaces, 1, function(x){similitud[x[1],x[2]]})
par(mfrow=c(1,1))
plot(similitudpares, dropoutperc)

intracom <- membMCL[enlaces[,1]]
idx <- apply(enlaces, 1, function(x){membMCL[x[1]]!=membMCL[x[2]]})
#intracom[idx] <- 0
intracom[idx] <- max(intracom)+1 
colores <- rainbow(length(unique(intracom))-1)
colores <- c(colores, 'black')

plot(similitudpares, dropoutperc, col=colores[intracom], 
     main="Similitud entre pares de nodos enlazados vs \n porcentaje de dropouts en ambos nodos del enlace", 
     xlab="similitud", ylab="porcentaje de dropouts")
legend('bottomleft', legend=c(seq(1:18),'dif'), col=colores, pch=1, cex=0.7, ncol = 2)

#por comunidad
xs <- par("usr")[1:2]
ys <- par("usr")[3:4]

par(mfrow=c(2,3))
for (i in 1:12){
  idx1 <- apply(enlaces, 1, function(x){if (membMCL[x[1]]==i | membMCL[x[2]]==i){TRUE} else {FALSE}})
  plot(similitudpares[idx1], dropoutperc[idx1], col=colores[intracom[idx1]], main=paste('Comunidad',i), xlim=xs, ylim=ys,
       xlab="similitud", ylab="porcentaje dropout")
}


# Asortatividad vectorial por comunidad----
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
  rhoch[i] <- r
  cat("r = ", r, "Comunidad", unique(celltypes)[i], "tamaño", table(celltypes[lccnm])[i],"\n")
}

#----
lccnm <- vertex_attr(lcc, name = "name") #todos estan en lcc y en nes0

edges <- as_edgelist(lcc, names = TRUE)
edgesim <- apply(edges, 1, function(x){matriz[x[1],x[2]]})
edgecoms <- apply(edges,2,function(x){membMCL[x]})

interior_idx <- edgecoms[,1]==edgecoms[,2]

#comparar enlaces de borde vs internos
hist(edgesim[interior_idx], freq=FALSE, main= "Distribucion de similitud de NES>0", col=rgb(0,0.6,0.4,0.5), xlab="similitud")
hist(edgesim[!interior_idx], freq=FALSE, add=TRUE, col=rgb(0,0,0.5,0.5))
legend("topleft", legend=c("interior","bordes"), fill=c(rgb(0,0.6,0.4,0.5), rgb(0,0,0.5,0.5)))

ccint <- sort(unique(as.vector(edgecoms[interior_idx,]))) #comunidades con enlaces internos


#Centralidad de los enlaces
kcoreness <- coreness(lcc)
enlaces_k <- enlaces
enlaces_k[,1] <- kcoreness[enlaces[,1]]
enlaces_k[,2] <- kcoreness[enlaces[,2]]
enlaces_k <- matrix(as.integer(enlaces_k), nrow=dim(enlaces_k)[1], ncol=dim(enlaces_k)[2])

edgecore <- apply(enlaces_k, 1, min)

mean_edgesim <- rep(0, length(unique(edgecore)))
for (i in 1:length(mean_edgesim)){
  mean_edgesim[i] <- mean(edgesim[edgecore==unique(edgecore)[i]])
}
par(mfrow=c(1,1))
plot(edgecore, edgesim, xlab="coreness de enlace", ylab="similitud de enlace")
points(unique(edgecore), mean_edgesim, col="green", pch=16)

corelist <- sort(unique(edgecore), decreasing = TRUE)
#grafo donde solo quedan los enlaces con coreness mayor a cierto umbral
g <- graph_from_edgelist(enlaces[edgecore>corelist[7],], directed = FALSE)
is.connected(g)
alphavec <- kcoreness[V(g)]/max(kcoreness[V(g)])
plot.igraph(g, vertex.label=NA, vertex.size=5, edge.width=0.5,
            vertex.color = rgb(0,1,0,alphavec))


#por comunidad
lista <- comms[as.character(ccint)]
colores <- rainbow(length(lista))
for (i in 1:length(lista)){
  nombres <- lista[[i]]
  nombres <- nombres[nombres %in% colnames(matriz)]
  g <- induced_subgraph(lcc, nombres)
  e <- as_edgelist(g)
  kcore <- coreness(g)
  enlaces_k <- e
  enlaces_k[,1] <- kcore[e[,1]]
  enlaces_k[,2] <- kcore[e[,2]]
  enlaces_k <- matrix(as.integer(enlaces_k), nrow=dim(enlaces_k)[1], ncol=dim(enlaces_k)[2])
  
  ecore <- apply(enlaces_k, 1, min)
  
  edgesim <- apply(e, 1, function(x){matriz[x[1],x[2]]})
  mean_edgesim <- rep(0, length(unique(ecore)))
  for (j in 1:length(mean_edgesim)){
    mean_edgesim[j] <- mean(edgesim[ecore==unique(ecore)[j]])
  }
  png(filename=paste0("fig/corenesscomm",names(lista)[i],".png"), width = 1600, height = 600, res = 150)
  par(mfrow=c(1,2))
  plot(ecore, edgesim, xlab="coreness de enlace", ylab="similitud de enlace", xlim=c(0,24),
       main=paste('Comunidad', names(lista)[i]))
  
  clr <- col2rgb(colores[i])
  points(unique(ecore), mean_edgesim, col=rgb(t(clr), maxColorValue = 255), pch=16)
  
  corelist <- sort(unique(ecore), decreasing = TRUE)
  #grafo donde solo quedan los enlaces con coreness mayor a cierto umbral
  #g1 <- graph_from_edgelist(e[ecore>corelist[7],], directed = FALSE)
  #is.connected(g1)
  alphav <- kcore[V(g)]/max(kcore[V(g)])
  plot.igraph(g, vertex.label=NA, vertex.size=5, edge.width=0.5,
              vertex.color = rgb(t(clr)/255, alpha = alphav), cex.main=0.5, 
              main = paste(vcount(g), 'nodos y', ecount(g), 'enlaces'))
  dev.off()
}

listasims <- lista
listacores <- lista
for (i in 1:length(lista)){
  nombres <- lista[[i]]
  nombres <- nombres[nombres %in% colnames(matriz)]
  g <- induced_subgraph(lcc, nombres)
  e <- as_edgelist(g)
  kcore <- coreness(g)
  enlaces_k <- e
  enlaces_k[,1] <- kcore[e[,1]]
  enlaces_k[,2] <- kcore[e[,2]]
  enlaces_k <- matrix(as.integer(enlaces_k), nrow=dim(enlaces_k)[1], ncol=dim(enlaces_k)[2])
  
  listacores[[i]] <- apply(enlaces_k, 1, min)
  listasims[[i]] <- apply(e, 1, function(x){matriz[x[1],x[2]]})
}

plot(listacores[[1]], listasims[[1]], col=colores[1], xlim=c(1,24), ylim=c(0.1,1), 
     xlab='coreness de enlace', ylab='similitud de enlace')
for (i in 2:length(listacores)){
  points(listacores[[i]], listasims[[i]], col=colores[i])
}


# atributo para enlaces: vecinos en comun
frac_comun <- apply(enlaces, 1, function(x){vecinoscomun(lcc,x) })


#---- REVISAR A PARTIR DE ACA
#distribucion de similitud entre pares por comunidad (comparo con azar pero sin incluir la propia comunidad)
ady <- as_adjacency_matrix(knn, sparse = FALSE)
matriz <- cos_sim(nes0)
N <- 100
par(mfrow=c(3,2))

for (i in 1:12){
  icomm <- which(membMCL[lccnm]==i)
  idx <- colnames(matriz) %in% names(membMCL[lccnm])[icomm]
  A <- ady[idx,idx]
  enlacesvect <- A[upper.tri(A)]
  enlacesvect <- ifelse(enlacesvect, TRUE, FALSE)
  X <- matriz[idx,idx]
  similvect <- X[upper.tri(X)]
  
  hist(similvect[enlacesvect], freq=FALSE, col=rgb(0,1,0,0.4), xlim=c(0,1), 
       main=paste("C",i,"tamaño",table(membMCL[lccnm])[i]), xlab = "")
  
  npar <- (length(X)-dim(X)[1])/2
  
  X <- matriz[!idx,!idx]
  pares <- X[upper.tri(X, diag = FALSE)]
  A <- matrix(0, nrow = N, ncol = npar)
  for (j in 1:dim(A)[1]){
    A[j,] <- sample(pares, npar)
  }
  hist(A, freq=FALSE, col=rgb(0,0,1,0.5), add=TRUE)
}


#------------------------------------------------------------------------

#para participacion y z-score (calculate_toproles.R)
membership <- membMCL[lccnm]
g_users <- lcc

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
