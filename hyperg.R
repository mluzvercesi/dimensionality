source("~/Documents/funciones.R")
#------------------------------------------------------------------------
load("~/Documents/results/Cpca.RData")
load("~/Documents/results/Csubconjunto.RData")
#------------------------------------------------------------------------
# require("package") or "package" %in% rownames(installed.packages())
library(GO.db)
library(GOstats)
library(org.Mm.eg.db)

ls("package:org.Mm.eg.db")
columns(org.Mm.eg.db)
# get("17708", org.Mm.egGENENAME) # para un solo gen, por ejemplo

#---------
#podria ser util para corregir algunos errores en los nombres
a<-read.table("Documents/Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt",sep="\t",header=TRUE)
table(a[,2])

a<-as.list(org.Mm.egALIAS2EG)
as.numeric(names(table(grep("Rik",names(a)))))
#---------

pca <- pcaCsub

gene_universe <- rownames(pca$rotation)
gene_universe_id <- sym2eg(gene_universe)
if (sum(is.na(gene_universe_id))/length(gene_universe_id) > 0.1){
  cat(sprintf("Fraccion del universo de genes sin mapear: %.2f", sum(is.na(selected_genes_id))/length(selected_genes_id)))
}

hgCutoff <- 1 #usar cutoff 1 para que devuelva todos los pvalues

# hgs es una matriz con los resultados del test hipergeometrico con ciertos pvalue para cada PC
hgs <- rep(NULL,8)

start_time <- Sys.time()
for (i in 1:10){
  orderedgenes <- genes_by_weight(pca,ncomp=i)
  selected_genes <- orderedgenes[1:100]
  selected_genes_id <- sym2eg(selected_genes)
  if (sum(is.na(selected_genes_id))/length(selected_genes_id)>0.1){
    cat(sprintf("Fraccion de genes sin mapear: %.2f", sum(is.na(selected_genes_id))/length(selected_genes_id)))
  }
  
  params <- new("GOHyperGParams", geneIds=selected_genes_id, universeGeneIds=gene_universe_id, annotation="org.Mm.eg.db",ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
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
      hgs <- rbind(hgs,hg)
    }#end if size
  }#end if pvalue
}#end loop i (pcs)
colnames(hgs)[8] <- "PC"
end_time <- Sys.time()
cat(sprintf("Tardó %.2f minutos\n", end_time-start_time))
#------------------------------------------------------------------------
load("~/Documents/results/Chg.RData")
neuro_offs <- get("GO:0022008",GOBPOFFSPRING) # neurogenesis

hgs <- hgsA
cat <- unique(hgs[,"GOBPID"]) #categorias GO
gopcs <- matrix(0L, nrow = length(cat), ncol = 11)
rownames(gopcs) <- cat
cols <- c(seq(1,10,1),"Neuro")
colnames(gopcs) <- cols

for (i in 1:10){
  # indices (booleanos) de las categorias que tienen PC=i
  ind_cat <- rownames(gopcs) %in% hgs[hgs[,"PC"]==i,"GOBPID"]
  # indices (con orden)
  ind_pval <- match(rownames(gopcs)[ind_cat],hgs[hgs[,"PC"]==i,"GOBPID"])
  # reemplazo de pvalues
  gopcs[ind_cat,i] <- hgs[hgs[,"PC"]==i,"Pvalue"][ind_pval]
}

ind_neuro <- rownames(gopcs) %in% neuro_offs
gopcs[ind_neuro,"Neuro"] <- .001

if(FALSE){
  a<-gopcs
  a[a==0]<-1 #hacer una matriz binaria
  a<-ifelse(a<.05,0,1)
  i11<-which(apply(a,1,sum)<ncol(a)) # GOs que aparecen
  
  heatmap.2(a[i11,])
  
  for(ipca in 1:11){
    igo<-which(a[,ipca]==0)
    print(Term(names(igo)))
    readline("Press ENTER")     
  }
  
}

heatmap(gopcs) # todo
heatmap(gopcs, Rowv=NA, Colv=NA) #sin reordenar con dendrograma
heatmap(sqrt(gopcs[gopcs[,"Neuro"]==1,])) # solo las neuro reescaladas


require(gplots) # para tener referencias
heatmap.2(sqrt(gopcs), trace="none")

#------------------------------------------------------------------------
# que genes aparecen en cada BP
genes_enbps <- lapply(hyperg_df[,"GOBPID"],function(x){genesinbp(selected_genes_id,x)})
names(genes_enbps) <- hyperg_df[,"GOBPID"]

# representacion visual de GO (los resultados, no todo GO)
# comparar distribucion de genes importantes en BP vs dist de todos los genes que no estan en BP
