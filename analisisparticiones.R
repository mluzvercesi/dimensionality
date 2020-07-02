load("resultados/datasetA.RData")
load("resultados/filtrosGO.RData")
load("resultados/pca_2filtros.RData")
load("resultados/gsea_zmeans.RData")
#funciones: cos_sim, make.knn, assortativity_vect

cells.nes <- cells.nes.zmeans
cells.padj <- cells.padj.zmeans
pca <- pca_lin

rm(pca_5000, pca_lin, cells.es.zmeans, cells.nes.zmeans, cells.pval.zmeans, cells.padj.zmeans)

cells_pc <- t(pca$x[,1:100]) # celulas (filas de x) en 100 primeros pcs, transpuesto
knn <- make.knn(cells_pc, k = 40)


if(FALSE){
  #Hipotesis nula: comparo con azar cambiando el orden de las columnas
  #armo distribucion para cada filtro de BPs
  N <- 100
  nfilt <- nrow(GOidx)
  ra_nesp <- matrix(0,ncol=N+1, nrow=nfilt)
  for (i in 1:nfilt){
    nes <- cells.nes[colnames(GOidx)[GOidx[i,]],]
    #nes <- ifelse(nes<0,0,nes)
    sim <- cos_sim(nes)
    ra_nesp[i,N+1] <- assortativity_vect(knn, sim)
    for (j in 1:N){
      cat("\r", sprintf("%.1f", 100*i*j/(N*nfilt)),"%")
      azar <- sample(colnames(nes))
      sim_azar <- sim[azar,azar]
      colnames(sim_azar) <- colnames(nes)
      rownames(sim_azar) <- colnames(nes)
      ra[i,j] <- assortativity_vect(knn, sim_azar)
    }#end loop repeticiones
  }#end loop filtros BPs
}#todos los filtros de BPs dan muy por encima del azar en asortatividad vect NES

par(mfrow=c(1,2))
for (i in 1:nrow(ra)){
  hist(ra[i,1:100], xlab="r", xlim=c(min(ra[i,]), max(ra[i,])), main=i)
  abline(v=ra[i,101], col="green")
}


load("resultados/nodos_clusters.RData")
clusters <- nodos[,-c(3,4)]
r_particion <- apply(clusters,2,function(x){0})
for (i in 1:length(r_particion)){
  r_particion[i] <- assortativity.nominal(knn, types=clusters[,i], directed=F)
}

library(fossil)
rands <- apply(clusters,2,function(x){0})
rands <- rands[-1]
rands <- rbind(rands, rands)
rownames(rands) <- c("ri", "ri.ajustado")
celltypes_nro <- as.numeric(factor(clusters[,"cell_type"]))
for (i in 1:ncol(rands)){
  rands[1,i] <- rand.index(clusters[,i], celltypes_nro)
  rands[2,i] <- adj.rand.index(clusters[,i], celltypes_nro)
}

# silhouette
library(cluster)
# [1] "cell_type"     "postnatal_day" "id"            "label"         "membMCL"       "NES1"          "NESp1"        
# [8] "NES2"          "NESp2"         "NES3"          "NESp3"         "NES4"          "NESp4"         "NES5"         
# [15] "NESp5"         "NES6"          "NESp6"         "NES7"          "NESp7"         "NES8"          "NESp8"        
# [22] "NES9"          "NESp9" 
i <- 23
X <- nodos[,i]
names(X) <- rownames(nodos)
table(X)
sil <- silhouette(X[X<23], dmatrix = 1-sim[names(X)[X<23],names(X)[X<23]], full=TRUE) #sim = cos_sim(nes)
x11()
plot(sil)


#por comunidad-----
load("resultados/particiones.RData")
lccnms <- vertex_attr(decompose(knn, max.comps = 1, min.vertices = 2)[[1]], name = "name")

#me quedo con los clusters mas grandes (mas de 20 celulas)
membnes <- (particiones[lccnms,"nes"])[particiones[lccnms,"nes"] %in% which(table(particiones[lccnms,"nes"])>20)]
names(membnes) <- lccnms[particiones[lccnms,"nes"] %in% which(table(particiones[lccnms,"nes"])>20)]

nes <- cells.nes[colnames(GOidx)[GOidx[1,]],]

nes_commean <- matrix(0, nrow=nrow(nes), ncol=length(unique(membnes)))
colnames(nes_commean) <- sort(unique(membnes))
rownames(nes_commean) <- rownames(nes)
nes_comsdev <- nes_commean
for (m in unique(membnes)){
  nes_commean[,m] <- apply(nes[, names(which(membnes==m))], 1, mean)
  nes_comsdev[,m] <- apply(nes[, names(which(membnes==m))], 1, sd)
}
plot(nes_commean, nes_comsdev, xlab="enriquecimiento medio por comunidad", ylab="desvio de enriquecimiento por comunidad")
#los de mayor enriquecimiento son los de menor desvio

bps_top20 <- apply(nes_commean, 2, function(x){rownames(nes_commean)[order(x, decreasing = T)[1:20]]})

#grafo simplificado de comunidades como nodos
gcoms <- grafo.coms(knn, membnes)
plot(simplify(gcoms, remove.loops=TRUE),
     vertex.size = sqrt(V(gcoms)$size),
     vertex.color = rgb(0,.5,.5, alpha = V(gcoms)$alpha), edge.arrow.size = .25)
table(droplevels(particiones[names(membnes),"cell_type"]), membnes)

Term(bps_top20[,1]) #para ver la descripcion de los terminos por comunidad
