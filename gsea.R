load("~/Documents/dimensionality/results/Asubconjunto.RData")
library(fgsea) 
# fgsea: Update old packages: 'AnnotationDbi', 'foreign', 'IRanges', 'nlme', 'Rcpp', 'RSQLite', 'sys'
library(org.Mm.eg.db)
#------------------------------------------------------------------------
# Pre procesamiento
load("~/Documents/dimensionality/results/Asubconjunto.RData")
rm(pcaAsub,pcaAsub_scaled)
X <- dataAsub

N_gen  <- dim(X)[1]
M_cell <- dim(X)[2]

n0 <- as.matrix(apply(X,2,function(x){sum(x==0)}))/N_gen 
#max0s <- c(0.92,0.88,0.95) # para los subconjuntos
X_0 <- X[,(n0<0.92)] # Quedan solo las que no superen el umbral

# Normalizar a 10000 cuentas por celula en A (5000 en C?)
X_0 <- apply(X_0,2,function(x){x*10000/sum(x)})

# Eliminar genes con poca variacion
mean_gen <- apply(X_0,1,function(x){mean(x)})
sdev_gen <- apply(X_0,1,function(x){norm_vect(x-mean(x))})
cv <- sdev_gen/mean_gen
cv[sdev_gen==0] <- 0 # corrijo los ceros que dieron nan
cv2 <- cv^2

# Valor medio de cv2 con bineado logaritmico
cv2mean <- mean_logbin(cv2,mean_gen)

plot(mean_gen,cv2,log="xy",
     xlab=expression(mu),ylab=expression(CV^2),
     main=paste('Dataset',dataset),sub=paste0('(',N_gen,' genes, ',dim(X_0)[2],' células)'))
points(cv2mean,col='green')

# ajuste lineal con bins "centrales"
#dataset A sub
x <- cv2mean$x[10:75]
y <- cv2mean$y[10:75]

fit_cv2mean <- nls(log(y) ~ a*log(x)+b, start = c(a = -1, b = 10))
y_cv2mean <- predict(fit_cv2mean, data.frame(x=cv2mean$x))
lines(cv2mean$x, exp(y_cv2mean), col="green")

# Me quedo con los genes que esten por encima de este ajuste lineal
pred_lin <- predict(fit_cv2mean, data.frame(x=mean_gen))
filtro_gen <- (cv2>exp(pred_lin))

X_1 <- X_0[filtro_gen,]
rm(cv,cv2,cv2mean,filtro_gen,fit_cv2mean,mean_gen,mean_logbin,pred_lin,sdev_gen,x,y,y_cv2mean)

N_gen  <- dim(X_1)[1]
M_cell <- dim(X_1)[2]

Xlog <- X_1
Xlog[X_1 != 0] <- log2(X_1[X_1 != 0])

Z <- apply(Xlog, 1, function(x){(x-mean(x))/norm_vect(x-mean(x))})
Z <- t(Z) #filas = genes; columnas = celulas
rownames(Z) <- sym2eg(rownames(Z))

#------------------------------------------------------------------------
# GSEA
# Cada columna de la matriz Z es una celula, es decir un vector de ranks
# save(Z,gopcsA,file="~/Documents/dimensionality/results/gseaA.RData")
load("~/Documents/dimensionality/results/gseaA.RData")

gonames <- rownames(gopcs)
Ngo <- length(gonames)
Ncells <- dim(Z)[2]
go_list <- list() #lista (de goids) de listas (de genes en cada goid)
go_terms <- lapply(gonames,function(x){ mappedLkeys(org.Mm.egGO2ALLEGS[x]) }) # genes de cada goid
go_list[gonames] <- go_terms

# ejemplo de grafico de ES
# plotEnrichment(pathway = go_list[[1]], stats = Z[,1], gseaParam = 1, ticksSize = 0.2)

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
save(cells.es, cells.nes, cells.pval, cells.padj, file="~/Documents/dimensionality/results/fgseaAsub100pcs500genes.RData")
# 1.968373 hs, 1553 de 2932 con 10k permutaciones, el resto con 1k (res_mix.RData)

#------------------------------------------------------------------------
# Similitud
logp <- -log(cells.padj)
sim.log <- cos_sim(logp)
cor.log <- cor(logp)

load("~/Documents/dimensionality/results/fgseaA_res.RData")
# cells.es, cells.nes, cells.padj, cells.pval, sim.log

sim.nes <- cos_sim(cells.nes) 
# NOTA: Hay NaN en los resultados de NES, 42 valores en 9 celulas (obs: en cada uno de ellos, pval=1). ??
