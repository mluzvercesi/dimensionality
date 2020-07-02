wd <- getwd()
setwd(paste0(wd,"/Tesis/dimensionality"))
load("results/Asubconjunto.RData")
make.knn <- function(expr, k){
  # grafo knn mutuos con correlacion a partir de expr
  # una matriz de expresion donde las columnas son celulas
  require(igraph)
  
  X <- cor(expr) # se podria hacer con similarity en vez de correlacion
  diag(X) <- 0
  
  # knn: ordeno cada fila (modulo), me quedo con los k mayores (k es umbral)
  ordenk <- apply(abs(X),1, function(x){order(x, decreasing = TRUE)})[1:k,]
  # las columnas de ordenk corresponden al orden de las filas de X
  
  A <- matrix(0L, nrow = dim(X)[1], ncol = dim(X)[2])
  rownames(A) <- rownames(X)
  colnames(A) <- colnames(X)
  for (i in 1:dim(X)[1]){
    A[i,ordenk[,i]] <- 1 #X[ordenk[,i],i]
  }
  ady <- (A+t(A))/2
  rm(i, A, X)
  ady <- floor(ady) #matriz de adyacencia de k vecinos cercanos mutuos
  
  #los nodos estan en el mismo orden que en la matriz
  knn <- graph_from_adjacency_matrix(ady, mode="undirected", diag = FALSE)
  return(knn)
}
cos_sim <- function(x){
  # calcula matriz de similarity entre columnas de matriz x
  y <- colSums(x^2)
  coseno <- (t(x) %*% x) / (sqrt( y %*% t(y) ))
  coseno[is.nan(coseno)] <- 0
  nombres <- colnames(x)
  colnames(coseno) <- nombres
  rownames(coseno) <- nombres
  return(coseno)
}
load("nodos_clusters.RData")
load("filtrosGO.RData")
load("pca_corregido_2filtros.RData")
load("gsea_filtrolineal_zmeans.RData")

cells.nes <- cells.nes.zmeans
cells.padj <- cells.padj.zmeans
rm(wd,pca_5000, pcaAsub, pcaAsub_scaled, cells.es.zmeans, cells.nes.zmeans, cells.pval.zmeans, cells.padj.zmeans)

assortativity_vect <- function(network, X){
  # network puede ser la red o la matriz de adyacencia
  # X matriz de similitud o correlacion entre propiedades vectoriales de los nodos
  # devuelve el coeficiente de asortatividad
  
  if (is.igraph(network)){ # si tengo grafo, armo matriz de adyacencia
    A <- as.matrix(as_adjacency_matrix(network))
  }else if (is.matrix(network)){
    A <- network
  }else{stop("No se reconoce como red")}
  
  diag(A) <- 0 # sin auto enlaces
  m <-  sum(A)/2 # enlaces totales
  
  if (m==0){
    r<-0
  }else{
    X <- X[colnames(A),colnames(A)]
    r <- sum(diag(A %*% X)) / (2*m)
  }
  return(r)
}


pca <- pca_lin
cells_pc <- t(pca$x[,1:100]) # celulas (filas de x) en 100 primeros pcs, transpuesto
knn <- make.knn(cells_pc, k = 40)

#Hipotesis nula: comparo con azar cambiando el orden de las columnas. armo distribucion (bootstrap)
N <- 100
N <- 0
nfilt <- nrow(GOidx)
ra_nesp <- matrix(0,ncol=N+1, nrow=nfilt)
for (i in 1:nfilt){
  nes <- cells.nes[colnames(GOidx)[GOidx[i,]],]
  nes <- ifelse(nes<0,0,nes)
  sim <- cos_sim(nes)
  ra_nesp[i,N+1] <- assortativity_vect(knn, sim)
  # for (j in 1:N){
  #   cat("\r", sprintf("%.1f", 100*i*j/(N*nfilt)),"%")
  #   azar <- sample(colnames(nes))
  #   sim_azar <- sim[azar,azar]
  #   colnames(sim_azar) <- colnames(nes)
  #   rownames(sim_azar) <- colnames(nes)
  #   ra[i,j] <- assortativity_vect(knn, sim_azar)
  # }
}

par(mfrow=c(3,3))
for (i in 1:nrow(ra)){
  hist(ra[i,1:100], xlab="r", xlim=c(min(ra[i,]), max(ra[i,])), main=i)
  abline(v=ra[i,101], col="green")
}


clusters <- nodos[,-c(3,4)]
r_particion <- apply(clusters,2,function(x){0})
for (i in 1:length(r_particion)){
  r_particion[i] <- assortativity.nominal(knn, types=clusters[,i], directed=F)
}

library(fossil)
rands <- apply(clusters,2,function(x){0})
rands <- rands[-1]
rands <- rbind(rands, rands)
rownames(rands) <- c("ri", "ri.ajustado")
celltypes_nro <- as.numeric(factor(clusters[,"cell_type"]))
for (i in 1:ncol(rands)){
  rands[1,i] <- rand.index(clusters[,i], celltypes_nro)
  rands[2,i] <- adj.rand.index(clusters[,i], celltypes_nro)
}




# silhouette
library(cluster)
# [1] "cell_type"     "postnatal_day" "id"            "label"         "membMCL"       "NES1"          "NESp1"        
# [8] "NES2"          "NESp2"         "NES3"          "NESp3"         "NES4"          "NESp4"         "NES5"         
# [15] "NESp5"         "NES6"          "NESp6"         "NES7"          "NESp7"         "NES8"          "NESp8"        
# [22] "NES9"          "NESp9" 
i <- 23
X <- nodos[,i]
names(X) <- rownames(nodos)
table(X)
sil <- silhouette(X[X<23], dmatrix = 1-sim[names(X)[X<23],names(X)[X<23]], full=TRUE) #sim = cos_sim(nes)
x11()
plot(sil)


#comparacion entre particiones----
contingencia <- table(membMCL,celltypes)
cntg_rows <- rowSums(contingencia)
cntg_cols <- colSums(contingencia)

cntg_uni <- contingencia
cntg_max <- contingencia
cntg_min <- contingencia

for (i in 1:dim(contingencia)[1]){
  for (j in 1:dim(contingencia)[2]){
    cntg_uni[i,j] <- contingencia[i,j]/(cntg_rows[i]+cntg_cols[j]-contingencia[i,j])
    cntg_max[i,j] <- contingencia[i,j]/max(cntg_rows[i],cntg_cols[j])
    cntg_min[i,j] <- contingencia[i,j]/min(cntg_rows[i],cntg_cols[j])
  }
}
rm(i,j)

#-----
#me quedo con los filtros lindos: 2 y 9
GOidx_buenos <- GOidx[c(2,9),]
clusters_buenos <- nodos[,c("cell_type", "postnatal_day","membMCL","NES2", "NESp2","NES9","NESp9")]
#                           cell_type postnatal_day membMCL NES2 NESp2 NES9 NESp9
# X10X43_1_AACGTTCTTAACCG.1      nIPC            12       3    8    10    7    13
# X10X43_1_ACTTCCCTACCCTC.1      nIPC            12       3    8    13    7    19

lccnms <- vertex_attr(decompose(knn, max.comps = 1, min.vertices = 2)[[1]], name = "name")

#para el NES2
membnes2 <- (clusters_buenos[lccnms,"NES9"])[clusters_buenos[lccnms,"NES9"] %in% which(table(clusters_buenos[lccnms,"NES9"])>20)]
names(membnes2) <- lccnms[clusters_buenos[lccnms,"NES9"] %in% which(table(clusters_buenos[lccnms,"NES9"])>20)]

nes <- cells.nes[colnames(GOidx)[GOidx[9,]],]

nes2_commean <- matrix(0, nrow=nrow(nes), ncol=length(unique(membnes2)))
colnames(nes2_commean) <- sort(unique(membnes2))
rownames(nes2_commean) <- rownames(nes)
for (m in unique(membnes2)){
  nes2_commean[,m] <- apply(nes[, names(which(membnes2==m))], 1, mean)
}

nes2_comsdev <- nes2_commean
for (m in unique(membnes2)){
  nes2_comsdev[,m] <- apply(nes[, names(which(membnes2==m))], 1, sd)
}

bps_nes2_top20 <- apply(nes2_commean, 2, function(x){rownames(nes2_commean)[order(x, decreasing = T)[1:20]]})

bps_nes2_f <- matrix(0, nrow=nrow(nes), ncol=length(unique(membnes2)))
colnames(bps_nes2_f) <- sort(unique(membnes2))
rownames(bps_nes2_f) <- rownames(nes)
fracciones <- rep(list(bps_nes2_f), 3)
names(fracciones) <- c("0", ".5", "1")
for (m in unique(membnes2)){
  X <- nes[, names(which(membnes2==m))]
  ncels <- ncol(X)
  for (fs in names(fracciones)){
    fracciones[[fs]][,m] <- apply(X,1,function(x){sum(x>as.numeric(fs))/ncels})
  }
}

bp_repr_xcom <- vector(mode = "list", length = length(unique(membnes2)))
names(bp_repr_xcom) <- sort(unique(membnes2))
for (m in unique(membnes2)){
  a <- (rownames(fracciones[[1]])[order(fracciones[[1]][,m], decreasing = T)])[1:20]
  b <- (rownames(fracciones[[2]])[order(fracciones[[2]][,m], decreasing = T)])[1:20]
  c <- (rownames(fracciones[[3]])[order(fracciones[[3]][,m], decreasing = T)])[1:20]
  d <- bps_nes2_top20[,m]
  #bp_repr_xcom[[m]] <- Term(Reduce(intersect, list(a,b,c,d)))
  bp_repr_xcom[[m]] <- Term(intersect(a,d))
}
unlist(lapply(bp_repr_xcom, length))

plot(nes2_commean[,1], nes2_comsdev[,1], ylim=c(.08, 1.4), xlim=c(-2,2))
for (i in 1:ncol(nes2_commean)){
  points(nes2_commean[,i], nes2_comsdev[,i], col=rainbow(ncol(nes2_commean))[i])
}
