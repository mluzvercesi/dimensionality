#------------------------------------------------------------------------
# Similitud
logp <- -log(cells.padj)
sim.log <- cos_sim(logp)
cor.log <- cor(logp)

sim.nes <- cos_sim(cells.nes) 
# NOTA: Hay NaN en los resultados de NES, 42 valores en 9 celulas (obs: en cada uno de ellos, pval=1). ??

load("~/Documents/dimensionality/results/fgseaA_res.RData")
# cells.es, cells.nes, cells.padj, cells.pval, sim.log

#------------------------------------------------------------------------
# Topologia
# Cell connected to its K (5-100) most similar cells, obtained using Euclidean distances on the PC-reduced expression space
library(igraph)

load("~/Documents/dimensionality/results/Asubconjunto.RData")
rm(pcaAsub_scaled)
cells_pc <- pcaAsub$x[,1:10] # celulas (filas) representadas en los 10 primeros pcs
sim_cells <- cos_sim(t(cells_pc))
cor_cells <- cor(t(cells_pc))

X <- cor_cells #con similarity da las mismas comunidades
diag(X) <- 0
# knn: ordeno cada fila (abs), me quedo con los k mayores (k es umbral)
ordenk <- apply(abs(X),1, function(x){order(x, decreasing = TRUE)})
# las columnas de ordenk corresponden al orden de las filas de X

k <- 20
A <- matrix(0L, nrow = dim(X)[1], ncol = dim(X)[2])
rownames(A) <- rownames(X)
colnames(A) <- colnames(X)
for (i in 1:dim(X)[1]){
  A[i,ordenk[1:k,i]] <- 1
}
rm(i)
ady <- (A+t(A))/2
ady <- floor(ady) #matriz de adyacencia de k vecinos mutuos cercanos

#los nodos estan en el mismo orden que en la matriz
G <- graph_from_adjacency_matrix(ady, mode="undirected", diag = FALSE)

library(hbm)
# I = 6,    2774 comunidades, max size = 7, 8 iteraciones
# I = 2,    739 comunidades, max size= 37, 16 iteraciones
# I = 1.6,  128 comunidades, max size= 95, 25 iteraciones, 13 coms de 1
# I = 1.4,  50 comunidades, max size=283, 32 iteraciones, 8 coms de 1
# I = 1.25, 19 comunidades, max size=691, 56 iteraciones, 7 coms de 1

jacc <- similarity(G, method = "jaccard")
knn <-  graph_from_adjacency_matrix(jacc, mode="undirected", weighted = TRUE, diag = FALSE)
com_jac  <- mcl(jacc, infl = 1.25, iter = 300, verbose = TRUE) # I=1.25
# A sub: 28 coms, maxsize=1051, 49 iteraciones, 7 coms de 1
# A: 96 coms, maxsize=989, 53 iteraciones, 45 coms de 1. NOTA: hay 54 coms con <5 nodos

names(com_jac) <- rownames(ady)
nro_com <- unique(com_jac)

#corrijo los numeros de las comunidades para que queden en orden
for (i in 1:length(nro_com)){
  com_jac[com_jac==nro_com[i]] <- i
}
rm(i)

#------------------------------------------------------------------------
load("~/Documents/dimensionality/results/A_knnjacc.RData") #tiene com_jac y com_jac_sub, knn y knn_sub

nro_com <- unique(com_jac)
single_com <- nro_com[table(com_jac)==1]
vecinos <- adjacent_vertices(knn, com_jac %in% single_com)
vecinos[as.numeric(lapply(vecinos,length))>0]


#para participacion y z-score (calculate_toproles.R)
membership <- com_jac
g_users <- knn

#------------------------------------------------------------------------
# ejemplos de grafos para ver como es el grafico zP de cada uno
g_bara <- sample_pa(100, directed = FALSE)
g_erd <- erdos.renyi.game(100, 1/100)

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
