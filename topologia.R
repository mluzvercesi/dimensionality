library(igraph)
library(Biobase)

#library(hbm) # NO! para clustering MCL usar el software posta
#com  <- mcl(jacc, infl = 1.25, iter = 300, verbose = TRUE)
# A sub: 11 coms, maxsize=1089, 39 iteraciones, 4 coms de 1
# A: 43 coms, maxsize=2743, 60 iteraciones, 15 coms de 1. NOTA: hay 19 coms con <5 nodos
# I = 6,    2774 comunidades, max size = 7, 8 iteraciones
# I = 2,    739 comunidades, max size= 37, 16 iteraciones
# I = 1.6,  128 comunidades, max size= 95, 25 iteraciones, 13 coms de 1
# I = 1.4,  50 comunidades, max size=283, 32 iteraciones, 8 coms de 1
# I = 1.25, 19 comunidades, max size=691, 56 iteraciones, 7 coms de 1

#------------------------------------------------------------------------
# Topologia
# Primero armo el grafo knn
# Cell connected to its K (5-100) most similar cells, obtained using Euclidean distances on the PC-reduced expression space

load("~/Documents/dimensionality/results/Asubconjunto.RData")
rm(pcaAsub_scaled, dataAsub)
cells_pc <- pcaAsub$x[,1:100] # celulas (filas) representadas en los 100 primeros pcs
cor_cells <- cor(t(cells_pc))

X <- cor_cells # se podria hacer con similarity en vez de correlacion pero da las mismas comunidades
diag(X) <- 0
k <- 40
# knn: ordeno cada fila (modulo), me quedo con los k mayores (k es umbral)
ordenk <- apply(abs(X),1, function(x){order(x, decreasing = TRUE)})
# las columnas de ordenk corresponden al orden de las filas de X
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
knn <- graph_from_adjacency_matrix(ady, mode="undirected", diag = FALSE)

jaccmat <- similarity(knn, method = "jaccard")
rownames(jaccmat) <- vertex_attr(knn, name = "name")
colnames(jaccmat) <- vertex_attr(knn, name = "name")
jaccgr <-  graph_from_adjacency_matrix(jaccmat, mode="undirected", weighted = TRUE, diag = FALSE)


# MCL (Ojo: mcl indexa desde 0)
df <- jaccgr
outfile<- "tmp.txt"

cmd    <- paste("mcl - --abc -I 1.25 -o - >",outfile,sep="")
pw     <- pipe(cmd,open="wb")
write.graph(df,file=pw,format="ncol")
close(pw)

pr<-pipe(paste("cat ",outfile,sep=""),open="r")
lines<-readLines(pr)
comms<-strsplit(split="\t", lines)
close(pr)

names(comms)<-seq_along(comms) #comms es una lista (de comunidades) de listas (de genes)
a<-reverseSplit((comms))
membMCL<-as.numeric(unlist(a))
names(membMCL)<-names(a) #membMCL es un vector de genes con membership

#----
# me tengo que quedar con el largest connected component ahora o entre knn y jaccard. miro que coincidan:
components(knn)$csize
components(jaccgr)$csize
table(membMCL)

lcc <- components(knn)$membership %in% which.max(components(knn)$csize)
knn.lcc <- induced_subgraph(knn, lcc, impl = "copy_and_delete")

lccnm <- vertex_attr(knn.lcc, name = "name")

sum(lccnm %in% names(membMCL)) #todos los nodos del lcc estan en el resultado de MCL
table(membMCL[lccnm])
table(membMCL) #los unicos nodos que faltan en el lcc son de dos pequeñas comunidades

#------------------------------------------------------------------------
load("~/Documents/dimensionality/results/A_knn40.RData") #tiene com_jac, com_jac_sub, knn, knn_sub, knn.jac, knn.jac_sub

metadata <- read.csv("~/Documents/dimensionality/Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt", 
                     sep="\t", header=TRUE, row.names=1)
rownames(metadata) <- sub("-",".",rownames(metadata))
rownames(metadata) <- paste0("X",rownames(metadata))

celltypes <- as.character(metadata[names(membMCL),"cell_type"])
names(celltypes) <- names(membMCL)

celltypes_nro <- as.numeric(factor(celltypes))
names(celltypes_nro) <- names(celltypes)


# Indice Rand
indice.rand(membMCL, celltypes_nro)

library(fossil)
rand.index(membMCL, celltypes_nro)
adj.rand.index(membMCL, celltypes_nro)


# Asortatividad de etiquetas MCL y celltype
assortativity_nominal(induced_subgraph(jaccgr, names(membMCL)), types=membMCL, directed=FALSE)
assortativity_nominal(induced_subgraph(jaccgr, names(celltypes_nro)), types=celltypes_nro, directed=FALSE)

assortativity_nominal(induced_subgraph(knn, names(membMCL)), types=membMCL, directed=FALSE)
assortativity_nominal(induced_subgraph(knn, names(celltypes_nro)), types=celltypes_nro, directed=FALSE)

#------------------------------------------------------------------------
# Similitud GSEA
load("~/Documents/dimensionality/results/fgseaA_res.RData") # cells.es, cells.nes, cells.padj, cells.pval, asort
load("~/Documents/dimensionality/results/fgseaAsub100pcs500genes.RData")

# NOTA: Hay NaN en los resultados de NES, en 14 celulas (obs: en cada uno de ellos, pval=1). ??
# NES: enrichment score normalized to mean enrichment of random samples of the same size
nes <- cells.nes[,!is.nan(apply(cells.nes, 2, min))] #tiro las celulas problematicas

# Puede que muchas cosas den similares porque se parecen en que ninguna esta enriquecida en cierta GO
# Convertir los negativos de ES y NES (son los mismos en ambas) en 0s
es0 <- ifelse(cells.es < 0, 0, cells.es)
nes0 <- ifelse(nes < 0, 0, nes)

# Filtrar ES y NES por padj: tomo -log y me fijo los significativos
logpadj <- -log(cells.padj)
logpadj <- ifelse(cells.es < 0, 0, logpadj)
#logpadj <- ifelse(logpadj > -log(0.15), logpadj, 0)
# CUALQUIER filtro para pval y padj empeora la asortatividad de toda la red y tambien por comunidad


# Asortatividad vectorial
asort <- read.table("~/Documents/dimensionality/results/resumen_asort.txt", sep = "\t")

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
#----

#Hipotesis nula: comparo con el azar cambiando el orden de las columnas
# armo distribucion (bootstrap)
N <- 500
ra <- rep(0,N)
for (j in 1:N){
  azar <- nes0[,sample(ncol(nes0))]
  colnames(azar) <- colnames(nes0)
  sim.nes_azar <- cos_sim(azar)
  red <- induced_subgraph(knn, colnames(azar))
  ra[j] <- assortativity_vect(red, sim.nes_azar)
}

hist(ra,xlim=c(0.7,0.8))
abline(v=asort[3,2],col="green")


#hacer lo mismo para azar con correlacion
N <- 100
ra.corr <- rep(0,N)
for (j in 1:N){
  azar <- nes0[,sample(ncol(nes0))]
  colnames(azar) <- colnames(nes0)
  nes_azar <- cor(azar)
  red <- induced_subgraph(knn, colnames(azar))
  ra.corr[j] <- assortativity_vect(red, nes_azar)
}

hist(ra.corr)#,xlim=c(0.7,0.8))
abline(v=asort[4,2],col="green")


#distribucion (histograma) de similitud entre pares, por comunidad
# comparo con azar pero sin incluir la propia comunidad
matriz <- cos_sim(nes0)
N <- 100
par(mfrow=c(3,2))
for (i in 1:6){
  cuales <- which(membMCL[lccnm]==i)
  idx <- colnames(matriz) %in% names(membMCL[lccnm])[cuales]
  X <- matriz[idx,idx]
  hist(X[upper.tri(X, diag = FALSE)], freq=FALSE,
       main=paste("C",i,"tamaño",table(membMCL[lccnm])[i]), xlab = "")
  
  npar <- (length(X)-dim(X)[1])/2
  
  X <- matriz[!idx,!idx]
  pares <- X[upper.tri(X, diag = FALSE)]
  A <- matrix(0, nrow = N, ncol = npar)
  for (j in 1:dim(A)[1]){
    A[j,] <- sample(pares, npar)
  }
  hist(A, freq=FALSE, col=rgb(1,0,0,0.5), add=TRUE)
}


# silhouette
library(cluster)
a <- names(membMCL) %in% colnames(nes0)
b <- colnames(nes0) %in% names(membMCL)
X <- membMCL[a]
Z <- nes0[,b]
sil <- silhouette(X[X<23], dmatrix = 1-cos_sim(Z[,names(X[X<23])]), full=TRUE)
x11()
plot(sil)


# atributo para enlaces: vecinos en comun
enlaces <- as_edgelist(knn.lcc, names = TRUE)
frac_comun <- apply(enlaces, 1, function(x){vecinoscomun(knn.lcc,x) })


#Centralidad de los enlaces
kcoreness <- coreness(knn.lcc)
enlaces <- as_edgelist(knn.lcc, names = TRUE)
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


#comparacion entre distintas particiones
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

#para participacion y z-score (calculate_toproles.R)
membership <- membMCL[lccnm]
g_users <- knn.lcc

plot(kcoreness, p_i)
