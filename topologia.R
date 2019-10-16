library(igraph)
library(hbm) #para clustering MCL

#------------------------------------------------------------------------
# Similitud GSEA
logp <- -log(cells.padj)
sim.log <- cos_sim(logp)
cor.log <- cor(logp)

sim.nes <- cos_sim(cells.nes) 
# NOTA: Hay NaN en los resultados de NES, 42 valores en 9 celulas (obs: en cada uno de ellos, pval=1). ??

load("~/Documents/dimensionality/results/fgseaA_res.RData") # cells.es, cells.nes, cells.padj, cells.pval, sim.log
#------------------------------------------------------------------------
# Topologia
# Primero armo el grafo knn
# Cell connected to its K (5-100) most similar cells, obtained using Euclidean distances on the PC-reduced expression space

load("~/Documents/dimensionality/results/Asubconjunto.RData")
rm(pcaAsub_scaled)
cells_pc <- pcaAsub$x[,1:10] # celulas (filas) representadas en los 10 primeros pcs
cor_cells <- cor(t(cells_pc))

X <- cor_cells # se podria hacer con similarity en vez de correlacion pero da las mismas comunidades
diag(X) <- 0
# knn: ordeno cada fila (modulo), me quedo con los k mayores (k es umbral)
ordenk <- apply(abs(X),1, function(x){order(x, decreasing = TRUE)})
# las columnas de ordenk corresponden al orden de las filas de X

k <- 40
A <- matrix(0L, nrow = dim(X)[1], ncol = dim(X)[2])
rownames(A) <- rownames(X)
colnames(A) <- colnames(X)
for (i in 1:dim(X)[1]){
  A[i,ordenk[1:k,i]] <- 1
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
which(components(knn.jac_sub)$membership>1)
which(components(knn_sub)$membership>1)
table(com_jac_sub)
com_jac_sub[which(components(knn_sub)$membership>1)]

lcc_sub <- which(components(knn_sub)$membership==1)

components(knn)$csize
components(knn.jac)$csize
com_jac[components(knn)$membership %in% which(components(knn)$csize==1)]

lcc <- components(knn)$membership %in% which(components(knn)$csize==1)
lcc <- !lcc

knn_lcc <- induced_subgraph(knn, lcc, impl = "copy_and_delete")
knn_sub_lcc <- induced_subgraph(knn_sub, lcc_sub, impl = "copy_and_delete")

#------------------------------------------------------------------------
load("~/Documents/dimensionality/results/A_knnjacc_k40.RData") #tiene com_jac y com_jac_sub, knn.jac y knn.jac_sub

metadata <- read.csv("~/Documents/dimensionality/Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt", 
                     sep="\t", header=TRUE, row.names=1)
rownames(metadata) <- sub("-",".",rownames(metadata))
rownames(metadata) <- paste0("X",rownames(metadata))

celltypes <- as.character(metadata[names(com_jac),"cell_type"])
names(celltypes) <- names(com_jac)

celltypes_nro <- rep(0, length(celltypes))
names(celltypes_nro) <- names(celltypes)
for (i in 1:length(unique(celltypes))){
  celltypes_nro[celltypes==unique(celltypes)[i]] <- i
}
rm(i)


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


indice.rand(com_jac, celltypes_nro)
library(fossil)
rand.index(com_jac, celltypes_nro)
adj.rand.index(com_jac, celltypes_nro)


# Asortatividad
# tengo que agregar los atributos mcl y cell al grafo
set_vertex_attr(knn.jac, "mcl", value=com_jac) #?
assortativity_nominal(knn.jac, types=V(knn.jac)$mcl)

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
membership <- com_jac
g_users <- knn

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
