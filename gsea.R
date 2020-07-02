library(fgsea)
library(org.Mm.eg.db)

# Pre procesamiento----
#lo mismo que para pca
load("resultados/datasetA.RData")
load("resultados/pca_2filtros.RData")
rm(dataA, pca_5000)
X <- dataAsub
pca <- pca_lin

igen <- (apply(X,1,sum)>20) # genes que aparezcan al menos 20 veces en todas las celulas
X <- X[igen,]
igen <- (apply(X,1, function(x){sum(x=!0)})/dim(X)[2])<0.6 # genes que aparezcan en menos del 60% de las celulas
X <- X[igen,]

nGen  <- dim(X)[1]
minGen <- 600 # es la cantidad minima de genes que deben expresarse
max0cell <- 1-minGen/nGen

n0 <- as.matrix(apply(X,2,function(x){sum(x==0)}))/nGen 
X0 <- X[,(n0<max0cell)] # Quedan solo las que no superen el umbral

# Eliminar genes con poca variacion
mean_gen <- apply(X0,1, mean)
sdev_gen <- apply(X0,1, sd)
cv <- sdev_gen/mean_gen
cv[sdev_gen==0] <- 0 # corrijo los ceros que dieron nan
cv2 <- cv^2

plot(mean_gen,cv2,log="xy",xlab=expression(mu),ylab=expression(CV^2),
     sub=paste0('(',nGen,' genes, ',dim(X0)[2],' celulas)'))

x <- mean_gen
y <- cv2
fit <- nls(log(y) ~ a*log(x)+b, start = c(a = -1, b = 10))
pred <- predict(fit, data.frame(x=mean_gen))
filtro_gen <- (cv2>exp(pred)) # solo los genes por encima de este ajuste lineal

X1 <- X0[filtro_gen,]
rm(cv,cv2,filtro_gen,fit,mean_gen,pred,sdev_gen,x,y,n0)

# Normalizar a 10000 cuentas por celula en A (5000 en C?)
X1 <- apply(X1,2,function(x){x*10000/sum(x)})
X1 <- round(X1)

Xlog <- X1
Xlog[X1 != 0] <- log2(X1[X1 != 0])

# Z: ranks para gsea (como fold change)----
# cuenta original aplicado a todas las celulas por igual
Z <- apply(Xlog, 1, function(x){(x-mean(x))/sd(x)})
Z <- t(Z) #para volver a filas = genes; columnas = celulas
rownames(Z) <- sym2eg(rownames(Z))
zraw <- Z


# sin comparar con los primeros k vecinos de cada celula
cells_pc <- t(pca$x[,1:100])
knn <- make.knn(cells_pc, k = 40)

#agregar verificacion de que nodos del knn = columnas de Xlog
nms_intersect <- intersect(get.vertex.attribute(knn, name='name'), colnames(Xlog))
length(nms_intersect)/length(union(get.vertex.attribute(knn, name='name'), colnames(Xlog)))
knn_intersect <- induced.subgraph(knn, vids=nms_intersect)
Xlog_intersect <- Xlog[, nms_intersect]
zmeans <- Xlog_intersect
m <- dim(zmeans)[2]

for (i in 1:m){
  cat("\r", sprintf("%.1f", 100*i/m),"%")
  cel <- colnames(Xlog_intersect)[i]
  vecinos <- neighbors(knn_intersect, cel)$name
  if (isEmpty(vecinos)){
    media <- Xlog_intersect[,cel]
    resto <- Xlog_intersect
  }else{
    media <- apply(Xlog_intersect[,c(cel,vecinos)], 1, mean)
    idxrest <- !(colnames(Xlog_intersect) %in% c(cel,vecinos))
    resto <- Xlog_intersect[,idxrest]
    resto <- cbind(resto,media)
  }
  mediaresto <- apply(resto, 1, mean)
  sdresto <- apply(resto, 1, sd)
  zmeans[,i] <- (media-mediaresto)/sdresto
}
rownames(zmeans) <- sym2eg(rownames(zmeans))

#puedo cambiar los infinitos por el maximo
if (sum(is.infinite(zmeans))){
  zmeans[is.infinite(zmeans)] <- ceiling(max(zmeans[!(is.infinite(zmeans))]))
}

# #CORREGIR esto
infomapcoms <- infomap.community(knn)
infomapmemb <- membership(infomapcoms)
zinfomap <- matrix(NA, ncol=length(unique(infomapmemb)), nrow=nrow(Xlog))
rownames(zinfomap) <- sym2eg(rownames(Xlog))
for (i in 1:ncol(zinfomap)){
  cidx <- names(infomapmemb)[infomapmemb==i]
  if (length(cidx)>1){
    zinfomap[,i] <- apply(Z[,cidx], 1, mean)
  }else{
    zinfomap[,i] <- Z[,cidx]
  }
}


# GSEA----
# Cada columna de la matriz Z es una celula, es decir un vector de ranks
load("resultados/hgs.RData")

gonames <- rownames(gopcs) #ver hyperg.R
Ngo <- length(gonames)
go_list <- list() #lista (de goids) de listas (de genes en cada goid)
go_terms <- lapply(gonames,function(x){ mappedLkeys(org.Mm.egGO2ALLEGS[x]) }) # genes de cada goid
go_list[gonames] <- go_terms

Z <- zmeans
Ncells <- dim(Z)[2]

# ejemplo de grafico de ES
#plotEnrichment(pathway = go_list[[1]], stats = Z[,1], gseaParam = 1, ticksSize = 0.2)

# Matrices: p value, p BH-adjusted, ES, NES
cells.pval <- matrix(0, nrow=Ngo, ncol=Ncells)
rownames(cells.pval) <- gonames
colnames(cells.pval) <- colnames(Z)
rm(gonames,go_terms,Ngo)

cells.padj <- cells.pval
cells.es <- cells.pval
cells.nes <- cells.pval

# fgsea devuelve tabla de 8 variables: pathway, pval, padj, ES, NES, nMoreExtreme, size, leadingEdge
start_time <- Sys.time()
for (i in 1:Ncells){
  cat("\r", sprintf("%.1f", 100*i/Ncells), "%")
  zranks <- Z[,i]
  gsea_res <- fgsea(pathways = go_list, stats = zranks, nperm = 1000)
  cells.pval[,i] <- gsea_res$pval
  cells.padj[,i] <- gsea_res$padj
  cells.es[,i] <- gsea_res$ES
  cells.nes[,i] <- gsea_res$NES
}
end_time <- Sys.time()
print(end_time-start_time)
rm(i,p,dospor,zranks,gsea_res,start_time,end_time)


load("resultados/gsea_zmeans.RData")
cells.nes <- cells.nes.zmeans
cells.padj <- cells.padj.zmeans
sum(is.na(cells.padj))
nacol <- is.na(apply(cells.padj, 2, min))
narow <- is.na(apply(cells.padj, 1, min))
sum(nacol)
sum(narow)
# me quedo con el menor de estos dos, para eliminar la menor cantidad de info posible
# o siempre elimino filas (BPs) para no perder celulas
cells.nes <- cells.nes[!narow,]
cells.padj <- cells.padj[!narow,]


# filtro de BPs (reduccion de dimension)----
cellsperbp <- apply(cells.padj, 1, function(x){sum(x<.05)})
hist(cellsperbp)
# filtros: entre mediana y 3er cuartil
quantile(cellsperbp)
GOidx <- cellsperbp>as.numeric(quantile(cellsperbp)[3]) & cellsperbp<as.numeric(quantile(cellsperbp)[4])

# en ese intervalo hay mas BPs de desarrollo de sist nervioso
library(GO.db)
nervdev_gos <- get("GO:0007399",GOBPOFFSPRING)
hist(cellsperbp[nervdev_gos[nervdev_gos %in% rownames(cells.padj)]], breaks= as.numeric(quantile(cellsperbp)), plot=F)$counts

#a este filtro tambien le puedo agregar una condicion sobre information content
library(GOSemSim)
mmGO <- godata("org.Mm.eg.db", ont="BP", computeIC=T)
# select <- AnnotationDbi::select
mmGO@IC <- computarIC("org.Mm.eg.db", keytype = "ENTREZID", ont = "BP")
plot(cellsperbp, mmGO@IC[names(cellsperbp)])
GOidx_IC <- (GOidx & mmGO@IC[names(cellsperbp)]>6)

GOidx <- rbind(GOidx, GOidx_IC)
GOidx[,colnames(GOidx) %in% nervdev_gos] <- TRUE


#otros filtros que probe----
load("resultados/filtrosGO_viejos.RData")
ncelmin <- as.numeric(names(which(cumsum(rev(table(cellsperbp)))>100)[1]))
GOidx0 <- cellsperbp>=ncelmin # elijo N tq resulte mas cerca de 100 BPs
GOidx0 <- cellsperbp>50 & cellsperbp<ncol(cells.padj)/2 & mmGO@IC[names(cellsperbp)]>6
load("resultados/BPs_prop.v2.RData")
bp.prop <- bp.prop[,c("IC", "ncelulas","N0lcc","N","Nlcc","Ncom")]
fraccion <- bp.prop[,"N0lcc"]/bp.prop[,"ncelulas"]
GOidx0 <- fraccion>.8
GOidx0 <- fraccion>.9
GOidx0 <- fraccion>.8 & bp.prop[,"ncelulas"]>20
GOidx0 <- fraccion>.9 & bp.prop[,"ncelulas"]>20
GOidx0 <- fraccion>.8 & bp.prop[,"IC"]>6


#grafico de IC vs celulas por BP, con cuartiles y BPs de desarrollo del sist nerv
library(plotly)
vline <- function(x = 0, color = "red") {
  list(type = "line", 
       y0 = 0, y1 = 1, yref = "paper",x0 = x, x1 = x, 
       line = list(color = color, dash = 'dash'))
}
f <- list(family = "Courier New, monospace",size = 18,color = "#7f7f7f")
x <- list(title = "celulas por BP",titlefont = f)
y <- list(title = "IC",titlefont = f)
fig <- plot_ly(type = 'scatter', mode = 'markers')
fig <- fig %>%
  add_trace(
    x = cellsperbp[!names(cellsperbp) %in% nervdev_gos], 
    y = mmGO@IC[names(cellsperbp)[!names(cellsperbp) %in% nervdev_gos]],
    text = Term(names(cellsperbp[!names(cellsperbp) %in% nervdev_gos])),
    hoverinfo = 'text',
    marker = list(color='black'),
    showlegend = F
  )
fig <- fig %>%
  add_trace(
    x = cellsperbp[names(cellsperbp) %in% nervdev_gos], 
    y = mmGO@IC[names(cellsperbp)[names(cellsperbp) %in% nervdev_gos]],
    text = Term(names(cellsperbp[names(cellsperbp) %in% nervdev_gos])),
    hoverinfo = 'text',
    marker = list(color='green'),
    showlegend = F
  )
fig <- fig %>%
  layout(shapes = list(vline(as.numeric(quantile(cellsperbp))[2]),
                       vline(as.numeric(quantile(cellsperbp))[3]),
                       vline(as.numeric(quantile(cellsperbp))[4])))
fig <- fig %>% layout(xaxis = x, yaxis = y)
fig
