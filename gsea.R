library(fgsea) 
# fgsea: Update old packages: 'AnnotationDbi', 'foreign', 'IRanges', 'nlme', 'Rcpp', 'RSQLite', 'sys'
library(org.Mm.eg.db)
#------------------------------------------------------------------------
# Pre procesamiento
load("results/Asubconjunto.RData")
rm(pcaAsub_scaled)
X <- dataAsub

N_gen  <- dim(X)[1]
M_cell <- dim(X)[2]

n0 <- as.matrix(apply(X,2,function(x){sum(x==0)}))/N_gen 
#max0s <- c(0.92,0.88,0.95) # para los subconjuntos
X0 <- X[,(n0<0.92)] # Quedan solo las que no superen el umbral

# Normalizar a 10000 cuentas por celula en A (5000 en C?)
X0 <- apply(X0,2,function(x){x*10000/sum(x)})

# Eliminar genes con poca variacion
mean_gen <- apply(X0,1,function(x){mean(x)})
sdev_gen <- apply(X0,1,function(x){norm_vect(x-mean(x))})
cv <- sdev_gen/mean_gen
cv[sdev_gen==0] <- 0 # corrijo los ceros que dieron nan
cv2 <- cv^2

# Valor medio de cv2 con bineado logaritmico
cv2mean <- mean_logbin(cv2,mean_gen)

plot(mean_gen,cv2,log="xy",
     xlab=expression(mu),ylab=expression(CV^2),
     main=paste('Dataset',dataset),sub=paste0('(',N_gen,' genes, ',dim(X0)[2],' celulas)'))
points(cv2mean,col='green')

# ajuste lineal con bins "centrales" dataset A sub
x <- cv2mean$x[10:75]
y <- cv2mean$y[10:75]

fit <- nls(log(y) ~ a*log(x)+b, start = c(a = -1, b = 10))
cv2pred <- predict(fit, data.frame(x=cv2mean$x))
lines(cv2mean$x, exp(cv2pred), col="green")

# Me quedo con los genes que esten por encima de este ajuste lineal
pred <- predict(fit, data.frame(x=mean_gen))
filtro_gen <- (cv2>exp(pred))

X1 <- X0[filtro_gen,]
rm(cv,cv2,cv2mean,filtro_gen,fit,mean_gen,mean_logbin,pred,sdev_gen,x,y,cv2pred, n0)

N_gen  <- dim(X1)[1]
M_cell <- dim(X1)[2]

Xlog <- X1
Xlog[X1 != 0] <- log2(X1[X1 != 0])


# cuenta original de fold change, aplicado a todas las celulas por igual
Z <- apply(Xlog, 1, function(x){(x-mean(x))/sd(x)})
Z <- t(Z) #filas = genes; columnas = celulas
rownames(Z) <- sym2eg(rownames(Z))


# fold change sin comparar con los primeros k vecinos de cada celula
cells_pc <- t(pcaAsub$x[,1:100])
knn <- make.knn(cells_pc, k = 40)

for (cel in colnames(Xlog)){
  vecinos <- neighbors(knn, cel)$name
  media <- apply(Xlog[,c(cel,vecinos)], 1, mean)
  z
}

#------------------------------------------------------------------------
# GSEA
# Cada columna de la matriz Z es una celula, es decir un vector de ranks
load("results/gseaA.RData") #AHORA NO SIRVE Z,gopcsA

gonames <- rownames(gopcs)
Ngo <- length(gonames)
Ncells <- dim(Z)[2]
go_list <- list() #lista (de goids) de listas (de genes en cada goid)
go_terms <- lapply(gonames,function(x){ mappedLkeys(org.Mm.egGO2ALLEGS[x]) }) # genes de cada goid
go_list[gonames] <- go_terms

# ejemplo de grafico de ES
#plotEnrichment(pathway = go_list[[1]], stats = Z[,1], gseaParam = 1, ticksSize = 0.2)

# Matrices: p value, p BH-adjusted, ES, NES
cells.pval <- matrix(rep(0.0, Ngo*Ncells), nrow=Ngo, ncol=Ncells)
rownames(cells.pval) <- gonames
colnames(cells.pval) <- colnames(Z)
rm(gonames,go_terms,Ngo)

cells.padj <- cells.pval
cells.es <- cells.pval
cells.nes <- cells.pval

# fgsea devuelve tabla de 8 variables: pathway, pval, padj, ES, NES, nMoreExtreme, size, leadingEdge
dospor <- round(0.02*Ncells)
start_time <- Sys.time()
for (i in 1:Ncells){
  if(i%%dospor==0){
    p <- 100*i/Ncells
    cat(sprintf("%.1f ", p))
  }
  zranks <- Z[,i]
  gsea_res <- fgsea(pathways = go_list, stats = zranks, nperm = 1000)
  cells.pval[,i] <- gsea_res$pval
  cells.padj[,i] <- gsea_res$padj
  cells.es[,i] <- gsea_res$ES
  cells.nes[,i] <- gsea_res$NES
}
rm(i,p,dospor,zranks,gsea_res)
end_time <- Sys.time()
print(end_time-start_time)
rm(start_time,end_time)


load("results/fgseaAsub100pcs500genes.RData") #cells.es, cells.nes, cells.pval, cells.padj

sum(is.na(cells.padj))
nacol <- is.na(apply(cells.padj, 2, min))
narow <- is.na(apply(cells.padj, 1, min))
sum(nacol)
sum(narow)
# me quedo con el menor de estos dos, para eliminar la menor cantidad de info posible
# me fijo si eliminando esas entradas desaparecen los NA:
sum(is.na(cells.padj[,!nacol]))

cells.es <- cells.es[,!nacol]
cells.nes <- cells.nes[,!nacol]
cells.pval <- cells.pval[,!nacol]
cells.padj <- cells.padj[,!nacol]

GOidx <- apply(cells.padj, 1, function(x){if (sum(x<0.05)>2) {TRUE} else {FALSE}})
dim(cells.padj[GOidx,])
