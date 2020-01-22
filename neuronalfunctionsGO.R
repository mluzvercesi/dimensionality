library(GO.db)
GOfuns <- read.table(file="GOfuns.txt", header=FALSE, sep="\n", row.names = NULL)
GOfuns <- as.character(GOfuns[,1])

tbl=toTable(GOTERM)
tbl[grep(faltan[4], tbl$Term),] -> x

allnervdev <- get("GO:0007399",GOBPOFFSPRING) # nervous system development (incluye neurogen)
allneurogn <- get("GO:0022008",GOBPOFFSPRING) # neurogenesis

get(GOBPid, org.Mm.egGO2ALLEGS)


#porcentajes de los genes de cada regulon que pertenecen a cada funcion neuronal GO
idx <- tbl$Term %in% GOfuns
GOfunsID <- unique(tbl$go_id[idx])
porcentajes <- matrix(rep(0, length(TFlist)*length(GOfunsID)), nrow=length(TFlist), ncol=length(GOfunsID))
rownames(porcentajes) <- names(TFlist)
colnames(porcentajes) <- GOfunsID

for (i in 1:dim(porcentajes)[1]){
  lista <- TFlist[[i]]
  L <- length(lista)
  for (j in 1:dim(porcentajes)[2]){
    goid <- colnames(porcentajes)[j]
    porcentajes[i,j] <- length(genesinbp(mapIds(org.Mm.eg.db, keys=lista, keytype = "ENSEMBL", column="ENTREZID"), goid))/L
  }
}


# sobrerrepresentacion GO de cada regulon
library(GOstats)

hgCutoff <- 1 #usar cutoff 1 para que devuelva todos los pvalues

hgs <- rep(NULL,8)

gene_universe <- as.character(unique(TFnet[,"Target"]))
gene_universe_id <- mapIds(org.Mm.eg.db, keys=gene_universe, keytype = "ENSEMBL", column="ENTREZID")
for (i in 11:112){
  selected_genes <- TFlist[[i]]
  L <- length(selected_genes)
  if (L<500 & L>10){
    selected_genes_id <- mapIds(org.Mm.eg.db, keys=selected_genes, keytype = "ENSEMBL", column="ENTREZID")
    params <- new("GOHyperGParams", geneIds=selected_genes_id, universeGeneIds=gene_universe_id, annotation="org.Mm.eg.db",
                  ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
    hyperg_results <- hyperGTest(params)
    hyperg_df <- summary(hyperg_results)
    # devuelve GOBPID (ids GO de c/categoria), pvalue, oddsratio, expcount, 
    # count (cantidad de selected genes en c/cat), size, term
    
    # corrijo los pvalues
    hyperg_df[,'Pvalue'] <- p.adjust(hyperg_df[,'Pvalue'],method='fdr')
    
    index <- hyperg_df[,"Pvalue"]<0.15
    # cambiar este pvalue si quiero abarcar mas o menos, para no tener que armar la matriz cada vez
    
    # agrego el PC solo si hay pvalue menor
    if (sum(index)>0){
      hg <- hyperg_df[index,]
      hg <- hg[hg[,'Size']<200,]
      hg <- hg[hg[,'Size']>10,]
      
      if (dim(hg)[1]>0){ # agrego el PC solo si tiene elementos
        hg <- cbind(hg,rep(i,dim(hg)[1]))
        names(hg)[8] <- "PC"
        hgs <- rbind(hgs,hg)
      }#end if size
    }#end if pvalue
  }#end loop L
}#end loop i (pcs)
#colnames(hgs)[8] <- "PC

hgs_redux <- rep(NULL,8)
for (i in 1:length(unique(hgs[,"Regulon"]))){
  r <- unique(hgs[,"Regulon"])[i]
  filas <- hgs[,"Regulon"]==r
  X <- hgs[filas,]
  filas <- hgs[filas,"Pvalue"]==min(hgs[filas,"Pvalue"])
  hgs_redux <- rbind(hgs_redux, X[filas,])
}

X <- hgs_redux[,c(T,T,F,F,F,F,F,T)]
idx <- which(X[,"GOBPID"] %in% colnames(porcentajes))
X <- X[idx,]
X <- cbind(X, rep(0, dim(X)[1]))
X <- cbind(X, rep(0, dim(X)[1]))
colnames(X)[4:5] <- c("Porcentaje","IsMax")


for (i in 1:dim(X)[1]){
  p <- porcentajes[as.numeric(X[i,"Regulon"]), X[i,"GOBPID"]]
  X[i, "Porcentaje"] <- p
  X[i, "IsMax"] <- maximos[as.numeric(X[i,"Regulon"])]==p
}
