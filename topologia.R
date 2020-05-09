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

make.knn <- function(expr, k){
  # grafo knn mutuos con correlacion a partir de expr
  # una matriz de expresion donde las columnas son celulas
  
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


load("results/Asubconjunto.RData")
rm(pcaAsub_scaled)
cells_pc <- t(pcaAsub$x[,1:100]) # celulas (filas de x) en 100 primeros pcs, transpuesto
knn <- make.knn(cells_pc, k = 40)

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

#
a <- read.csv(file = "membMCL.txt",header = TRUE, sep=" ", row.names = 1)
membMCL <- a[,1]
names(membMCL) <- rownames(a)
comms <- split(names(membMCL), membMCL)
rm(a)
#

#----
# me tengo que quedar con el largest connected component ahora o entre knn y jaccard. miro que coincidan:
components(knn)$csize
components(jaccgr)$csize

lcc <- decompose(knn, max.comps = 1, min.vertices = 2)[[1]]
lccnm <- vertex_attr(lcc, name = "name")

sum(lccnm %in% names(membMCL)) #todos los nodos del lcc estan en el resultado de MCL
table(membMCL[lccnm])
table(membMCL) #los unicos nodos que faltan en el lcc son de dos pequenias comunidades

#celltype----
metadata <- read.csv("Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt", sep="\t", header=TRUE, row.names=1)
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
cat("Asortatividad de comunidades MCL:", 
    assortativity_nominal(induced_subgraph(knn, names(membMCL)), types=membMCL, directed=FALSE))
cat("Asortatividad de celltypes:", 
    assortativity_nominal(induced_subgraph(knn, names(celltypes_nro)), types=celltypes_nro, directed=FALSE))


# Similitud GSEA----
#load("results/fgseaA_res.RData") # cells.es, cells.nes, cells.padj, cells.pval, asort
load("results/fgseaAsub100pcs500genes.RData")

#Filtro por p ajustado: solo me quedo con la categoria GO si tengo mas de dos entradas con p menor que .05
GOidx <- apply(cells.padj, 1, function(x){if (sum(x<0.05)>2) {TRUE} else {FALSE}})
# si hay NA en los resultados ver gsea.R

nes <- cells.nes[GOidx,]
es <- cells.es[GOidx,]
logpadj <- -log(cells.padj[GOidx,])
logp <- -log(cells.pval[GOidx,])

# Convertir los negativos de ES (NES) en 0s para evitar similitud en falta de enriquecimiento
es0 <- ifelse(es < 0, 0, es)
nes0 <- ifelse(nes < 0, 0, nes)
logpadj0 <- ifelse(es < 0, 0, logpadj)
logp0 <- ifelse(es < 0, 0, logp)

# Filtrar ES y NES por padj significativos? CUALQUIER filtro empeora la asortatividad (tambien por comunidad)

# Asortatividad vectorial
asort <- read.csv("results/resumen_asort.txt", row.names = 1)

A <- es0
assortativity_vect(induced_subgraph(knn, colnames(A)), cos_sim(A))

#Hipotesis nula: comparo con azar cambiando el orden de las columnas. armo distribucion (bootstrap)
N <- 500
ra <- rep(0,N)
A <- logpadj0
for (j in 1:N){
  azar <- A[,sample(ncol(A))]
  colnames(azar) <- colnames(A)
  sim_azar <- cos_sim(azar)
  red <- induced_subgraph(knn, colnames(azar))
  ra[j] <- assortativity_vect(red, sim_azar)
  n <- 100*j/N
  if (n%%2==0){
    cat(paste0(n,"% "))
  }
}
r1 <- assortativity_vect(induced_subgraph(knn, colnames(A)), cos_sim(A))
hist(ra, main = "Distribucion de asortatividad \n en una red al azar", xlab="r", xlim=c(min(ra), r1))
abline(v=r1, col="green")


#hacer lo mismo para azar con correlacion (tarda mas)
N <- 100
ra.corr <- rep(0,N)
A <- nes0
for (j in 1:N){
  azar <- A[,sample(ncol(A))]
  colnames(azar) <- colnames(A)
  cor_azar <- cor(azar)
  red <- induced_subgraph(knn, colnames(azar))
  ra.corr[j] <- assortativity_vect(red, cor_azar)
  n <- 100*j/N
  if (n%%2==0){
    cat(paste0(n,"% "))
  }
}
hist(ra.corr)
abline(v=assortativity_vect(induced_subgraph(knn, colnames(A)), cor(A)),col="green")

#distribución de similitud en enlaces vs no-enlaces
ady <- as_adjacency_matrix(knn, sparse = FALSE)
enlacesvect <- ady[upper.tri(ady)]
enlacesvect <- ifelse(enlacesvect, TRUE, FALSE)
matriz <- cos_sim(nes0)
similvect <- matriz[upper.tri(matriz)]
hist(similvect[enlacesvect], freq=F, col=rgb(0,0.6,0.4,0.5), xlim=c(0,1), 
     main="Distribución de similitud de NES>0", xlab="similitud")
hist(similvect[!enlacesvect], freq=F, add=T, col=rgb(0,0,0.5,0.5)) # el total es casi igual a este
legend("topleft", legend=c("enlazados","no enlazados"), fill=c(rgb(0,0.6,0.4,0.5), rgb(0,0,0.5,0.5)))


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


#Entropia - REVISAR
S_rows <- apply(cntg_uni,1,function(x){-log(x)*x})
S_rows <- t(S_rows)
S_rows[cntg_uni==0] <- 0
M <- dim(cntg_uni)[2]
S_row_max <- -log(1/M)

S_cols <- apply(cntg_uni,2,function(x){-log(x)*x})
S_cols[cntg_uni==0] <- 0
N <- dim(cntg_uni)[1]
S_col_max <- -log(1/N)

