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

# Topologia----
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
ady <- (A+t(A))/2
rm(i, A, X)
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

rm(df, outfile, cmd, a, lines, pr, pw)
#----
# me tengo que quedar con el largest connected component ahora o entre knn y jaccard. miro que coincidan:
components(knn)$csize
components(jaccgr)$csize

lcc <- components(knn)$membership %in% which.max(components(knn)$csize)
knn.lcc <- induced_subgraph(knn, lcc, impl = "copy_and_delete")

lccnm <- vertex_attr(knn.lcc, name = "name")

sum(lccnm %in% names(membMCL)) #todos los nodos del lcc estan en el resultado de MCL
table(membMCL[lccnm])
table(membMCL) #los unicos nodos que faltan en el lcc son de dos pequeñas comunidades

#celltype----
#load("~/Documents/dimensionality/results/A_knn40.RData") #tiene com_jac, com_jac_sub, knn, knn_sub, knn.jac, knn.jac_sub

metadata <- read.csv("~/Documents/dimensionality/Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt", 
                     sep="\t", header=TRUE, row.names=1)
rownames(metadata) <- sub("-",".",rownames(metadata))
rownames(metadata) <- paste0("X",rownames(metadata))

#
a <- read.csv(file = "membMCL.txt",header = TRUE, sep=" ", row.names = 1)
membMCL <- as.array(t(a))
names(membMCL) <- rownames(a)
comms <- split(names(membMCL), membMCL)
rm(a)
#

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
cat("Asortatividad de comunidades MCL:", 
    assortativity_nominal(induced_subgraph(knn, names(membMCL)), types=membMCL, directed=FALSE))
cat("Asortatividad de celltypes:", 
    assortativity_nominal(induced_subgraph(knn, names(celltypes_nro)), types=celltypes_nro, directed=FALSE))


# Similitud GSEA----
#load("~/Documents/dimensionality/results/fgseaA_res.RData") # cells.es, cells.nes, cells.padj, cells.pval, asort
load("results/fgseaAsub100pcs500genes.RData")

# NOTA: Hay NaN en los resultados de NES (enrichment score normalized to mean enrichment of random samples of the same size)
# en 14 celulas (obs: en cada uno de ellos, pval=1). ??
nes <- cells.nes[,!is.nan(apply(cells.nes, 2, min))] #tiro las celulas problematicas

# Convertir los negativos de ES (NES) en 0s para evitar similitud en falta de enriquecimiento
es0 <- ifelse(cells.es < 0, 0, cells.es) #tal vez conviene eliminar esas 9 celulas para que tengan la misma dimension
nes0 <- ifelse(nes < 0, 0, nes)

# Filtrar ES y NES por padj: tomo -log y me fijo los significativos
logpadj <- -log(cells.padj)
logpadj <- ifelse(cells.es < 0, 0, logpadj)
#logpadj <- ifelse(logpadj > -log(0.15), logpadj, 0)
# CUALQUIER filtro para pval y padj empeora la asortatividad de toda la red y tambien por comunidad


# Asortatividad vectorial
asort <- read.table("results/resumen_asort.txt", sep = "\t")

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
hist(ra,xlim=c(0.7,0.8), main = "Distribucion de asortatividad \n en una red al azar", xlab="r")
abline(v=asort["sim (ES>0)","NES"],col="green")


#hacer lo mismo para azar con correlacion (tarda mas)
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
abline(v=asort["corr (ES>0)","NES"],col="green")


#distribucion de similitud entre pares por comunidad (comparo con azar pero sin incluir la propia comunidad)
matriz <- cos_sim(nes0)
N <- 100
par(mfrow=c(3,2))
for (i in 1:6){
  cuales <- which(membMCL[lccnm]==i)
  idx <- colnames(matriz) %in% names(membMCL[lccnm])[cuales]
  X <- matriz[idx,idx]
  hist(X[upper.tri(X, diag = FALSE)], freq=FALSE,
       main=paste("C",i,"tamaño",table(membMCL[lccnm])[i]), xlab = "", ylim=c(0,10))
  
  npar <- (length(X)-dim(X)[1])/2
  
  X <- matriz[!idx,!idx]
  pares <- X[upper.tri(X, diag = FALSE)]
  A <- matrix(0, nrow = N, ncol = npar)
  for (j in 1:dim(A)[1]){
    A[j,] <- sample(pares, npar)
  }
  hist(A, freq=FALSE, col=rgb(1,0,0,0.5), add=TRUE)
}

#repito lo mismo, pero para pares que pertenezcan a distintas comunidades, sea cual sea
m <- membMCL[rownames(matriz)[rownames(matriz) %in% names(membMCL)]]
ccomp <- matrix(rep(0, length(m)^2), nrow=length(m), ncol=length(m),
               dimnames = list(names(m),names(m)))
for (i in unique(m)){
    nombres <- names(m)[m==i]
    ccomp[nombres,nombres] <- 1
}
rm(i,nombres)
idx <- ccomp[upper.tri(ccomp, diag = FALSE)]
X <- matriz[rownames(ccomp),colnames(ccomp)]
pares <- X[upper.tri(X, diag = FALSE)]
pares <- pares[!idx]

par(mfrow=c(3,2))
for (i in 1:6){
  cuales <- which(membMCL[lccnm]==i)
  idx <- colnames(matriz) %in% names(membMCL[lccnm])[cuales]
  X <- matriz[idx,idx]
  hist(X[upper.tri(X, diag = FALSE)], freq=FALSE,
       main=paste("C",i,"tamaño",table(membMCL[lccnm])[i]), xlab = "", ylim=c(0,10))
  
  npar <- (length(X)-dim(X)[1])/2
  
  A <- matrix(0, nrow = N, ncol = npar)
  for (j in 1:dim(A)[1]){
    A[j,] <- sample(pares, npar)
  }
  hist(A, freq=FALSE, col=rgb(1,0,0,0.5), add=TRUE)
}


# silhouette
library(cluster)
a <- names(membMCL) %in% colnames(nes0)
X <- membMCL[a]

b <- colnames(nes0) %in% names(membMCL)
Z <- nes0[,b]
sil <- silhouette(X[X<23], dmatrix = 1-cos_sim(Z[,names(X[X<23])]), full=TRUE)
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

#----
Ngen <- dim(dataAsub)[1]
enlaces <- as_edgelist(knn.lcc)
a <- apply(enlaces, 1, function(x){sum((dataAsub[,x[1]]+dataAsub[,x[2]])==0)})
dropoutperc <- a/Ngen

similitud <- cos_sim(dataAsub)
similitudpares <- apply(enlaces, 1, function(x){similitud[x[1],x[2]]})
plot(similitudpares, dropoutperc)

intracom <- membMCL[enlaces[,1]]
idx <- apply(enlaces, 1, function(x){membMCL[x[1]]!=membMCL[x[2]]})
intracom[idx] <- 0
intracom[idx] <- max(intracom)+1
colores <- rainbow(length(unique(intracom))-1)
colores <- c(colores, 'black')

plot(similitudpares, dropoutperc, col=colores[intracom])
legend('bottomleft', legend=c(seq(1:18),'dif'), col=colores, pch=1, cex=0.7)

#por comunidad
xs <- par("usr")[1:2]
ys <- par("usr")[3:4]

par(mfrow=c(2,3))
for (i in 1:12){
  idx1 <- apply(enlaces, 1, function(x){if (membMCL[x[1]]==i | membMCL[x[2]]==i){TRUE} else {FALSE}})
  plot(similitudpares[idx1], dropoutperc[idx1], col=colores[intracom[idx1]], main=paste('Comunidad',i), xlim=xs, ylim=ys)
}

# calcular similitud pero solo con las partes no nulas del vector?
