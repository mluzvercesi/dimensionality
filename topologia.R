library(igraph)
library(Biobase)
#library(hbm) # NO! para clustering MCL usar el software posta

# Topologia----
# Primero armo el grafo knn
# Cell connected to its K (5-100) most similar cells, obtained using Euclidean distances on the PC-reduced expression space
load("resultados/datasetA.RData")
load("resultados/pca_2filtros.RData")
rm(dataA, pca_5000)
pca <- pca_lin
cells_pc <- t(pca$x[,1:100]) # celulas (filas de x) en 100 primeros pcs, transpuesto
knn <- make.knn(cells_pc, k = 40) #tambien puedo usar matriz de NES

jaccmat <- similarity(knn, method = "jaccard")
rownames(jaccmat) <- vertex_attr(knn, name = "name")
colnames(jaccmat) <- vertex_attr(knn, name = "name")
jaccgr <-  graph_from_adjacency_matrix(jaccmat, mode="undirected", weighted = TRUE, diag = FALSE)

# MCL
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
# coinciden el largest connected component entre knn y jaccard
components(knn)$csize[1]==components(jaccgr)$csize[1]

lcc <- decompose(knn, max.comps = 1, min.vertices = 2)[[1]]
lccnm <- vertex_attr(lcc, name = "name")

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


# Asortatividad de clusters MCL y celltype
cat("Asortatividad de comunidades MCL:", 
    assortativity_nominal(induced_subgraph(knn, names(membMCL)), types=membMCL, directed=FALSE))
cat("Asortatividad de celltypes:", 
    assortativity_nominal(induced_subgraph(knn, names(celltypes_nro)), types=celltypes_nro, directed=FALSE))


# Similitud GSEA----
load("resultados/gsea_zmeans.RData")
cells.nes <- cells.nes.zmeans
cells.padj <- cells.padj.zmeans
narow <- is.na(apply(cells.padj, 1, min))
cells.nes <- cells.nes[!narow,]
cells.padj <- cells.padj[!narow,]

load("resultados/filtrosGO.RData")
GOidx[,colnames(GOidx) %in% nervdev_gos] <- TRUE
idx <- GOidx[1,]
nes <- cells.nes[idx,]
nes0 <- ifelse(nes < 0, 0, nes) #para evitar similitud en falta de enriquecimiento

# Asortatividad vectorial
rvectNES <- assortativity_vect(knn, cos_sim(nes))
rvectNESp <- assortativity_vect(knn, cos_sim(nes0))
cat("Asortatividad vectorial NES:", rvectNES, "\nAsortatividad vectorial NES+:", rvectNESp)


#Hipotesis nula: comparo con azar cambiando el orden de las columnas. armo distribucion (bootstrap)
N <- 500
ra <- rep(0,N+1)
A <- nes
similitud <- cos_sim(A)
ra[N+1] <- assortativity_vect(knn, similitud)
for (j in 1:N){
  azar <- sample(colnames(A))
  sim_azar <- similitud[azar,azar]
  colnames(sim_azar) <- colnames(A)
  rownames(sim_azar) <- colnames(A)
  ra[j] <- assortativity_vect(knn, sim_azar)
  n <- 100*j/N
  if (n%%2==0){
    cat(paste0(n,"% "))
  }
}
hist(ra[1:N], main = "Distribucion de asortatividad \n en una red al azar", xlab="r", xlim=c(min(ra), max(ra)))
abline(v=ra[N+1], col="green")

#hacer lo mismo para azar con correlacion (tarda mas)
N <- 100
ra.corr <- rep(0,N+1)
correlacion <- cor(A)
ra.corr[N+1] <- assortativity_vect(knn, correlacion)
for (j in 1:N){
  azar <- sample(colnames(A))
  cor_azar <- correlacion[azar,azar]
  colnames(cor_azar) <- colnames(A)
  rownames(cor_azar) <- colnames(A)
  ra.corr[j] <- assortativity_vect(knn, cor_azar)
  n <- 100*j/N
  if (n%%2==0){
    cat(paste0(n,"% "))
  }
}
hist(ra.corr[1:N], main = "Distribucion de asortatividad \n en una red al azar", xlab="r", xlim=c(min(ra.corr), max(ra.corr)))
abline(v=ra.corr[N+1],col="green")


#distribución de similitud en enlaces vs no-enlaces
ady <- as_adjacency_matrix(knn, sparse = FALSE)
enlacesvect <- ady[upper.tri(ady)]
enlacesvect <- ifelse(enlacesvect, TRUE, FALSE)
similvect <- similitud[upper.tri(similitud)]
hist(similvect[enlacesvect], freq=F, col=rgb(0,0.6,0.4,0.5), xlim=c(-1,1), 
     main="Distribución de similitud de NES>0", xlab="similitud")
hist(similvect[!enlacesvect], freq=F, add=T, col=rgb(0,0,0.5,0.5)) # el total es casi igual a este
legend("topleft", legend=c("enlazados","no enlazados"), fill=c(rgb(0,0.6,0.4,0.5), rgb(0,0,0.5,0.5)))


# silhouette
library(cluster)
a <- names(membMCL) %in% colnames(nes)
X <- membMCL[a]

b <- colnames(nes0) %in% names(membMCL)
Z <- nes0[,b]
sil <- silhouette(X[X<23], dmatrix = 1-cos_sim(Z[,names(X[X<23])]), full=TRUE)
x11()
plot(sil)

celltype <- as.numeric(factor(droplevels(metadata[rownames(pca$x),"cell_type"])))
dissimilarity <- 1-cos_sim(t(pca$x))
sil <- silhouette(celltype, dmatrix = dissimilarity, full=TRUE)
x11()
plot(sil, col=rainbow(length(unique(celltype))))

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
