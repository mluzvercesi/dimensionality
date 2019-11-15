library(igraph)
library(hbm) #para clustering MCL

#------------------------------------------------------------------------
# Topologia
# Primero armo el grafo knn
# Cell connected to its K (5-100) most similar cells, obtained using Euclidean distances on the PC-reduced expression space

load("~/Documents/dimensionality/results/Asubconjunto.RData")
rm(pcaAsub_scaled, dataAsub)
cells_pc <- pcaAsub$x[,1:10] # celulas (filas) representadas en los 10 primeros pcs
cor_cells <- cor(t(cells_pc))

X <- cor_cells # se podria hacer con similarity en vez de correlacion pero da las mismas comunidades
diag(X) <- 0
# knn: ordeno cada fila (modulo), me quedo con los k mayores (k es umbral)
ordenk <- apply(abs(X),1, function(x){order(x, decreasing = TRUE)})
# las columnas de ordenk corresponden al orden de las filas de X
k <- 40
ordenk <- ordenk[1:k,]

A <- matrix(0L, nrow = dim(X)[1], ncol = dim(X)[2])
rownames(A) <- rownames(X)
colnames(A) <- colnames(X)
for (i in 1:dim(X)[1]){
  A[i,ordenk[,i]] <- 1 #X[ordenk[,i],i]
}
rm(i)
ady <- (A+t(A))/2
ady <- floor(ady) #matriz de adyacencia de k vecinos cercanos mutuos

#los nodos estan en el mismo orden que en la matriz
knn_sub <- graph_from_adjacency_matrix(ady, mode="undirected", diag = FALSE)

jacc <- similarity(knn_sub, method = "jaccard")
knn.jac_sub <-  graph_from_adjacency_matrix(jacc, mode="undirected", weighted = TRUE, diag = FALSE)
com_jac_sub  <- mcl(jacc, infl = 1.25, iter = 300, verbose = TRUE)
# A sub: 11 coms, maxsize=1089, 39 iteraciones, 4 coms de 1
# A: 43 coms, maxsize=2743, 60 iteraciones, 15 coms de 1. NOTA: hay 19 coms con <5 nodos

# I = 6,    2774 comunidades, max size = 7, 8 iteraciones
# I = 2,    739 comunidades, max size= 37, 16 iteraciones
# I = 1.6,  128 comunidades, max size= 95, 25 iteraciones, 13 coms de 1
# I = 1.4,  50 comunidades, max size=283, 32 iteraciones, 8 coms de 1
# I = 1.25, 19 comunidades, max size=691, 56 iteraciones, 7 coms de 1

names(com_jac_sub) <- rownames(ady)
nro_com <- unique(com_jac_sub)

#corrijo los numeros de las comunidades para que queden en orden
for (i in 1:length(nro_com)){
  com_jac_sub[com_jac_sub==nro_com[i]] <- i
}
rm(i)

# me tengo que quedar con el largest connected component ahora o entre knn y jaccard. miro que coincidan:
components(knn_sub)$csize
components(knn.jac_sub)$csize
table(com_jac_sub)
which(components(knn.jac_sub)$membership>1)
which(components(knn_sub)$membership>1)
com_jac_sub[which(components(knn_sub)$membership>1)]

components(knn)$csize
components(knn.jac)$csize
com_jac[components(knn)$membership %in% which(components(knn)$csize==1)]

lcc <- components(knn.jac)$membership %in% which(components(knn.jac)$csize>1)
lcc_sub <- components(knn.jac_sub)$membership %in% which(components(knn.jac_sub)$csize>1)

knn_lcc <- induced_subgraph(knn, lcc, impl = "copy_and_delete")
knn_sub_lcc <- induced_subgraph(knn_sub, lcc_sub, impl = "copy_and_delete")

#------------------------------------------------------------------------
load("~/Documents/dimensionality/results/A_knn40.RData") #tiene com_jac, com_jac_sub, knn, knn_sub, knn.jac, knn.jac_sub

metadata <- read.csv("~/Documents/dimensionality/Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt", 
                     sep="\t", header=TRUE, row.names=1)
rownames(metadata) <- sub("-",".",rownames(metadata))
rownames(metadata) <- paste0("X",rownames(metadata))

celltypes <- as.character(metadata[names(com_jac),"cell_type"])
names(celltypes) <- names(com_jac)

celltypes_sub <- as.character(metadata[names(com_jac_sub),"cell_type"])
names(celltypes_sub) <- names(com_jac_sub)

celltypes_nro <- as.numeric(factor(celltypes))
names(celltypes_nro) <- names(celltypes)


# Indice Rand
indice.rand(com_jac, celltypes_nro)

library(fossil)
rand.index(com_jac, celltypes_nro)
adj.rand.index(com_jac, celltypes_nro)

rand.index(com_jac_sub, celltypes_nro[names(com_jac_sub)])
adj.rand.index(com_jac_sub, as.numeric(factor(celltypes_nro[names(com_jac_sub)])))


# Asortatividad de etiquetas MCL y celltype
assortativity_nominal(knn.jac, types=com_jac, directed=FALSE)
assortativity_nominal(knn.jac_sub, types=com_jac_sub, directed=FALSE)

assortativity_nominal(knn.jac, types=celltypes_nro, directed=FALSE)
assortativity_nominal(knn.jac_sub, types=celltypes_nro[names(com_jac_sub)], directed=FALSE)

#------------------------------------------------------------------------
# Similitud GSEA
load("~/Documents/dimensionality/results/fgseaA_res.RData") # cells.es, cells.nes, cells.padj, cells.pval, asort
# asort es la tabla de valores de asortatividad calculado con todas las distintas medidas

logp <- -log(cells.padj)

# NOTA: Hay NaN en los resultados de NES, 42 valores en 9 celulas (obs: en cada uno de ellos, pval=1). ??
# NES: enrichment score normalized to mean enrichment of random samples of the same size
nes <- cells.nes[,!is.nan(apply(cells.nes, 2, min))] #tiro las celulas problematicas
es <- cells.es[,!is.nan(apply(cells.nes, 2, min))]

sim.es <- cos_sim(cells.es)

# OJO! Puede que muchas cosas den similares porque se parecen en que ninguna de las dos esta enriquecida en cierta GO
# Ver: hist(as.dist(sim.es, diag = FALSE, upper = FALSE))
# O: hist(sim.es[upper.tri(sim.es)])

# Convertir en 0 los negativos de ES y NES (son los mismos en ambas)
es0 <- ifelse(cells.es < 0, 0, cells.es)
nes0 <- ifelse(nes < 0, 0, nes)


# Filtrar ES y NES por padj: tomo -log y me fijo los significativos
hist(-log(cells.padj))
abline(v=-log(0.05), col="blue")
abline(v=-log(0.1), col="green")
abline(v=-log(0.15), col="red")


# Asortatividad vectorial por comunidad
ncoms <- unique(com_jac_sub[lcc_sub])

#prueba------
matriz <- cos_sim(cells.es)
a <- ifelse(cells.padj[,!is.nan(apply(cells.nes, 2, min))]<0.15, nes, 0)
a <- ifelse(a>0, a, 0)
matriz <- cos_sim(a)


for (i in 1:length(ncoms)){
  nombres <- names(which(com_jac_sub[lcc_sub]==ncoms[i]))
  nombres <- nombres[nombres %in% colnames(matriz)]
  red <- induced_subgraph(knn_sub, nombres)
  r <- assortativity_vect(red, matriz[nombres,nombres])
  cat("r = ", r, "Comunidad", ncoms[i], "tamaño", table(com_jac_sub[lcc_sub])[i], "\n")
}

for (i in 1:length(unique(celltypes_sub))){
  nombres <- names(which(celltypes_sub[lcc_sub]==unique(celltypes_sub)[i]))
  nombres <- nombres[nombres %in% colnames(matriz)]
  red <- induced_subgraph(knn_sub, nombres)
  r <- assortativity_vect(red, matriz[nombres,nombres])
  cat("r = ", r, "Comunidad", unique(celltypes_sub)[i], "tamaño", table(celltypes_sub[lcc_sub])[i],"\n")
}
#----

#Hipotesis nula: comparo con el azar cambiando el orden de las columnas
for (i in 1:length(ncoms)){
  nombres <- which(com_jac_sub[lcc_sub]==ncoms[i])
  red <- induced_subgraph(knn_sub, nombres)
  r <- assortativity_vect(red, sim.es[nombres,nombres])
  cat("r = ", r, "Comunidad", ncoms[i], "tamaño", table(com_jac_sub[lcc_sub])[i], "\n")
}

for (i in 1:length(unique(celltypes_sub))){
  nombres <- which(celltypes_sub[lcc_sub]==unique(celltypes_sub)[i])
  red <- induced_subgraph(knn_sub, nombres)
  r <- assortativity_vect(red, sim.es[nombres,nombres])
  cat("r = ", r, "Comunidad", unique(celltypes_sub)[i], "tamaño", table(celltypes_sub[lcc_sub])[i],"\n")
}
ra <- 0
N <- 100
for (j in 1:N){
  azar <- cells.es[,sample(ncol(cells.es))]
  colnames(azar) <- colnames(cells.es)
  sim.es_azar <- cos_sim(azar)
  ra <- ra + assortativity_vect(knn_sub, sim.es_azar)
}
ra <- ra/N #es casi lo mismo usar la red completa que las comunidades
cat("r azar = ",r)


# silhouette
library(cluster)
sil <- silhouette(com_jac_sub, dmatrix = 1-sim.es, full=TRUE)
x11()
plot(sil)


#distribucion (histograma) de similitud entre pares, por comunidad
pares <- as.dist(sim.es, diag = FALSE, upper = FALSE)
com_ind <- c(1, 4, 5, 6)
par(mfrow=c(2,2))
for (i in 1:length(com_ind)){
  cuales <- which(com_jac_sub[lcc_sub]==ncoms[com_ind[i]])
  X <- sim.es[cuales, cuales]
  diag(X) <- NaN
  hist(as.dist(X, diag = FALSE, upper = FALSE), freq=FALSE,
       main=paste("C",ncoms[com_ind[i]],"tamaño",table(com_jac_sub[lcc_sub])[com_ind[i]]), xlab = "")
  
  A <- matrix(0, nrow = N, ncol = (length(X)-dim(X)[1])/2)
  for (j in 1:dim(A)[1]){
    A[j,] <- sample(pares, (length(X)-dim(X)[1])/2)
  }
  hist(A, freq=FALSE, col=rgb(1,0,0,0.5), add=TRUE)
}


#Centralidad de los enlaces
kcoreness <- coreness(knn_sub_lcc)
enlaces <- as_edgelist(knn_sub_lcc, names = TRUE)
enlaces_k <- enlaces
enlaces_k[,1] <- kcoreness[enlaces[,1]]
enlaces_k[,2] <- kcoreness[enlaces[,2]]
enlaces_k <- matrix(as.integer(enlaces_k), nrow=dim(enlaces_k)[1], ncol=dim(enlaces_k)[2])

edge_core <- apply(enlaces_k, 1, min)
edge_sim <- apply(enlaces, 1, function(x){ sim.es[x[1], x[2]] })

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

# participacion vs coreness de cada nodo
g_users <- knn_sub_lcc
membership <- com_jac_sub[lcc_sub]
#usar calculate_toproles.R
plot(kcoreness, p_i)


contingencia <- table(com_jac,celltypes)
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


#Entropia
S_rows <- apply(cntg_uni,1,function(x){-log(x)*x})
S_rows <- t(S_rows)
S_rows[cntg_uni==0] <- 0
M <- dim(cntg_uni)[2]
S_row_max <- -log(1/M)

S_cols <- apply(cntg_uni,2,function(x){-log(x)*x})
S_cols[cntg_uni==0] <- 0
N <- dim(cntg_uni)[1]
S_col_max <- -log(1/N)

#------------------------------------------------------------------------
nro_com <- unique(com_jac)
single_com <- nro_com[table(com_jac)==1]
vecinos <- adjacent_vertices(knn, com_jac %in% single_com)
vecinos[as.numeric(lapply(vecinos,length))>0]

#para participacion y z-score (calculate_toproles.R)
membership <- com_jac_sub[lcc_sub]
g_users <- knn_sub_lcc

#------------------------------------------------------------------------
#Comparacion rapida: grafo knn vs similarity pval
comm.sorted <- sort(sizes(comunidades))
comm.sorted <- comm.sorted[comm.sorted>2]
comm.sorted <- as.list(comm.sorted)

for (i in 1:length(comm.sorted)){
  n <- as.integer(names(comm.sorted)[i])
  ind <- comunidades$names[comunidades$membership==n]
  M <- sim.log[ind,ind]
  diag(M) <- NaN
  hist(M,main=n)
}
rm(i)
