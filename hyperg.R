source("~/Documents/funciones.R")
#------------------------------------------------------------------------
load("~/Documents/dimensionality/results/Apca.RData")
load("~/Documents/dimensionality/results/Asubconjunto.RData")

library(GO.db)
library(GOstats)
library(org.Mm.eg.db)

require(gplots) # para tener heatmap con referencias
#------------------------------------------------------------------------
# para ver que tiene el paquete de Mm
ls("package:org.Mm.eg.db")
columns(org.Mm.eg.db)
get("17708", org.Mm.egGENENAME) # para un solo gen, por ejemplo

#---------
#podria ser util para corregir algunos errores en los nombres:
a<-read.table("Documents/Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt",sep="\t",header=TRUE)
table(a[,2])

a<-as.list(org.Mm.egALIAS2EG)
as.numeric(names(table(grep("Rik",names(a)))))

#TEST HIPERGEOMETRICO PARA SOBRERREPRESENTACION--------------------------
# hgs es una matriz con resultados del test: pvalue para cada PC

pca <- pcaCsub

gene_universe <- rownames(pca$rotation)
gene_universe_id <- sym2eg(gene_universe)
if (sum(is.na(gene_universe_id))/length(gene_universe_id) > 0.1){
  cat(sprintf("Fraccion del universo de genes sin mapear: %.2f", sum(is.na(selected_genes_id))/length(selected_genes_id)))
}

hgCutoff <- 1 #usar cutoff 1 para que devuelva todos los pvalues

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
#REPRESENTACION VISUAL---------------------------------------------------
load("~/Documents/dimensionality/results/Ahg.RData")

neuro_offs <- get("GO:0022008",GOBPOFFSPRING) # neurogenesis
nervdev_offs <- get("GO:0007399",GOBPOFFSPRING) # nervous system development (incluye neurogen)

hgs <- hgsA
cat <- unique(hgs[,"GOBPID"]) #categorias GO
cols <- unique(hgs[,"PC"])
Npcs <- length(cols)
gopcs <- matrix(0L, nrow = length(cat), ncol = Npcs+1)
rownames(gopcs) <- cat
colnames(gopcs) <- c(cols,"Neuro")
rm(cat)

for (i in cols){
  # indices (booleanos) de las categorias que tienen PC=i
  ind_cat <- rownames(gopcs) %in% hgs[hgs[,"PC"]==i,"GOBPID"]
  # indices (con orden)
  ind_pval <- match(rownames(gopcs)[ind_cat],hgs[hgs[,"PC"]==i,"GOBPID"])
  # reemplazo de pvalues
  gopcs[ind_cat,i] <- hgs[hgs[,"PC"]==i,"Pvalue"][ind_pval]
}
rm(i,cols,ind_cat,ind_pval,Npcs)


ind_neuro <- rownames(gopcs) %in% neuro_offs
ind_nerv <- rownames(gopcs) %in% nervdev_offs
gopcs[ind_neuro,"Neuro"] <- .001
gopcs[ind_nerv,"Neuro"] <- .001
rm(ind_neuro, ind_nerv)

heatmap.2(gopcs) # todo
heatmap.2(sqrt(gopcs), trace="none") #reescalada
heatmap(sqrt(gopcs[gopcs[,"Neuro"]==.001,])) # solo las neuro, reescaladas


if(FALSE){
  gopcs01 <- gopcs
  gopcs01[gopcs01==0]<-1 #hacer una matriz binaria
  gopcs01 <- ifelse(gopcs01<.05,0,1)
  i11 <- which(apply(gopcs01,1,sum)<ncol(gopcs01)) # GOs que aparecen
  
  heatmap.2(gopcs01[i11,],trace="none")
  
  for(ipca in 1:11){
    igo<-which(gopcs01[,ipca]==0)
    print(Term(names(igo)))
    readline("Press ENTER")     
  }
}

#------------------------------------------------------------------------
# que genes aparecen en cada BP
genes_enbps <- lapply(hyperg_df[,"GOBPID"],function(x){genesinbp(selected_genes_id,x)})
names(genes_enbps) <- hyperg_df[,"GOBPID"]
# comparar distribucion de genes importantes en BP vs dist de todos los genes que no estan en BP

#------------------------------------------------------------------------
#SEMANTIC SIMILARITY-----------------------------------------------------
library(GOSemSim)
library("igraph")

mmGO <- godata("org.Mm.eg.db", ont="BP")

# ALGUNOS GO NO APARECEN EN LA BASE DE DATOS. Hay que agregar manualmente:
# mmGO@keys, mmGO@IC, mmGO@geneAnno (ENTREZID GO EVIDENCE ONTOLOGY)
indices <- rownames(gopcs01) %in% names(mmGO@IC)
faltan <- rownames(gopcs01)[!indices]


get(faltan[3], org.Mm.egGO2ALLEGS)
mapIds(org.Mm.eg.db, keys=faltan, column=c("ENTREZID"), keytype="GO")
universo_BP <- get("GO:0008150", org.Mm.egGO2ALLEGS)
universo_BP_GO <- eg2go(universo_BP)
a <- get("GO:0048731", org.Mm.egGO2ALLEGS)
a_GO <- eg2go(a)
-log2(length(unique(a_GO))/length(unique(universo_BP_GO))) # NO FUNCIONA

#ademas, sacar cosas que sobran de la lista. Usar unicamente evidencia experimental:
# Inferred from Experiment (EXP)
# Inferred from Direct Assay (IDA)
# Inferred from Physical Interaction (IPI)
# Inferred from Mutant Phenotype (IMP)
# Inferred from Genetic Interaction (IGI)
# Inferred from Expression Pattern (IEP)
# Inferred from High Throughput Experiment (HTP)
# Inferred from High Throughput Direct Assay (HDA)
# Inferred from High Throughput Mutant Phenotype (HMP)
# Inferred from High Throughput Genetic Interaction (HGI)
# Inferred from High Throughput Expression Pattern (HEP)

go_list <- rownames(gopcs01)[indices]
sim <- mgoSim(go_list,go_list,semData = mmGO, measure ="Resnik",combine=NULL)
sim2 <- mgoSim(go_list,go_list,semData = mmGO, measure ="Wang",combine=NULL)

#grafo completo
G <- graph_from_adjacency_matrix(adjmatrix = sim, weighted = TRUE, diag = FALSE, mode = "undirected")
coord <- layout_with_fr(G)
l <- norm_coords(coord, ymin=-1, ymax=1, xmin=-1, xmax=1)
max_weight <- max(E(G)$weight)

lista_nerv <- V(G)$name[V(G)$name %in% nervdev_offs]
colores2 <- c("blue","lightblue")

#si pertenece a neuro tipo=1, sino tipo=2
V(G)$type <- 2
V(G)$type[V(G)$name %in% lista_nerv] <- 1

png("~/Documents/dimensionality/fig/grafo_resnik_completo.png")
plot.igraph(G, layout=l, rescale=FALSE, vertex.size=10, vertex.color=adjustcolor(colores2[V(G)$type], alpha.f=.5),
            vertex.label=NA,edge.width=E(G)$weight/max_weight,
            main="Categorias GO (BP)")
legend('topright', legend=c('Neuro','No neuro'), pch=21, col="black", pt.bg=colores2)
dev.off()


colores3 <- c("blue","gray","lightblue")

for (i in 1:(length(colnames(gopcs))-1)){ #para cada pc
  lista_pc <- V(G)$name[V(G)$name %in% rownames(gopcs01)[gopcs01[,i]==0]]
  
  subG <- subgraph(G,unique(c(lista_pc,lista_nerv))) #subgrafo con ambos (pc y nerv)
  maximo <- max(E(subG)$weight)
  
  indices_nerv <- V(subG)$name %in% lista_nerv
  indices_pc <- V(subG)$name %in% lista_pc
  
  V(subG)$type <- 0
  V(subG)$type[indices_nerv] <- V(subG)$type[indices_nerv]+1
  V(subG)$type[indices_pc] <- V(subG)$type[indices_pc]+2
  
  ind_coord <- match(V(subG)$name,V(G)$name)
  if (length(V(subG)$name) != sum(V(subG)$name == V(G)$name[ind_coord]) ){
    print("Paso algo con los indices")
  }
  
  fn <- paste0("~/Documents/dimensionality/fig/grafo_resnik_pc",i)
  png(fn)
  plot.igraph(subG, layout=l[ind_coord,], rescale=FALSE, vertex.size=10,
              vertex.color=adjustcolor(colores3[V(subG)$type], alpha.f=.5), vertex.label=NA,
              edge.width=E(subG)$weight/maximo,
              main=paste("PC #",colnames(gopcs)[i]))
  legend('topright', legend=c('Neuro','PC','Ambos'), pch=21, col="black", pt.bg=colores3)
  dev.off()
}
rm(i, lista_pc, subG, maximo, indices_nerv,indices_pc,ind_coord,fn)
