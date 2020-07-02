library(GO.db)
library(GOstats)
library(org.Mm.eg.db)

# para ver que tiene el paquete de Mm
ls("package:org.Mm.eg.db")
columns(org.Mm.eg.db)
get("17708", org.Mm.egGENENAME) # por ejemplo, nombre de un gen

#----
#podria ser util para corregir algunos errores en los nombres:
a<-read.table("Neurogenesis/Linnarson_NatNeuro2018/GSE95315/GSE95315_metadata.txt",sep="\t",header=TRUE)
table(a[,2])
a<-as.list(org.Mm.egALIAS2EG)
as.numeric(names(table(grep("Rik",names(a)))))

# TEST HIPERGEOMETRICO PARA SOBRERREPRESENTACION----
# hgs es una matriz con resultados del test: pvalue para cada PC
load("resultados/datasetA.RData")
load("resultados/pca_2filtros.RData")
rm(dataA,pca_5000)

pca <- pca_lin

gene_universe <- rownames(pca$rotation)
gene_universe_id <- sym2eg(gene_universe)
if (sum(is.na(gene_universe_id))/length(gene_universe_id) > 0.1){
  cat(sprintf("Fraccion del universo de genes sin mapear: %.2f", sum(is.na(selected_genes_id))/length(selected_genes_id)))
}

hgCutoff <- 1 #usar cutoff 1 para que devuelva todos los pvalues, filtro despues

nmaxpcs <- 100 #cantidad de pcs que considera
ngenes <- 500 #genes por pc
hgs <- rep(NULL,8)

if(FALSE){
  #como elegir cantidad de genes? veo donde dejan de aportar peso
  orderedgenes <- genes_by_weight(pca,ncomp=1, retw = T)
  selected_genes <- orderedgenes[1:1000]
  plot(selected_genes, col=rainbow(nmaxpcs)[1])
  for (i in 2:nmaxpcs){
    orderedgenes <- genes_by_weight(pca,ncomp=i, retw=T)
    selected_genes <- orderedgenes[1:1000]
    points(selected_genes, col=rainbow(nmaxpcs)[i])
  }
}

start_time <- Sys.time()
for (i in 1:nmaxpcs){
  cat("\r", sprintf("%.1f", 100*i/nmaxpcs),"%")
  orderedgenes <- genes_by_weight(pca,ncomp=i)
  selected_genes <- orderedgenes[1:ngenes]
  selected_genes_id <- sym2eg(selected_genes)
  if (sum(is.na(selected_genes_id))/length(selected_genes_id)>0.1){
    cat(sprintf("Fraccion de genes sin mapear: %.2f", sum(is.na(selected_genes_id))/length(selected_genes_id)))
  }
  
  params <- new("GOHyperGParams", geneIds=selected_genes_id, universeGeneIds=gene_universe_id, annotation="org.Mm.eg.db",
                ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
  hyperg_results <- hyperGTest(params)
  hyperg_df <- summary(hyperg_results)
  # devuelve GOBPID, pvalue, oddsratio, expcount, 
  # count (cantidad de selected genes en c/cat), size, term
  
  hyperg_df[,'Pvalue'] <- p.adjust(hyperg_df[,'Pvalue'],method='fdr') # corrijo los pvalues
  
  index <- hyperg_df[,"Pvalue"]<0.15
  # cambiar este pvalue si quiero abarcar mas o menos, para no tener que armar la matriz cada vez
  
  # agrego el PC solo si hay algun pvalue menor
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
}#end loop i (pcs)
end_time <- Sys.time()
print(end_time-start_time)
rm(i, hg, hyperg_df, index, hyperg_results, params, selected_genes_id, selected_genes, orderedgenes,end_time, start_time)

# REPRESENTACION VISUAL----
require(gplots)
load("resultados/hgs.RData")

nervdev_gos <- get("GO:0007399",GOBPOFFSPRING) # nervous system development (incluye neurogenesis)
neurogen_gos <- get("GO:0022008",GOBPOFFSPRING) # neurogenesis

gonms <- unique(hgs[,"GOBPID"]) #categorias GO
pcnms <- unique(hgs[,"PC"])
gopcs <- matrix(1L, nrow = length(gonms), ncol = length(pcnms)+1)
rownames(gopcs) <- gonms
colnames(gopcs) <- c(pcnms,"Neuro")

for (i in pcnms){
  # submatriz con PC=i
  X <- hgs[hgs[,"PC"]==i,]
  # reemplazo de pvalues
  gopcs[X[,"GOBPID"],as.character(i)] <- X[,"Pvalue"]
}
igos <- rownames(gopcs) %in% nervdev_gos # neuro_gos
gopcs[igos,"Neuro"] <- .Machine$double.eps # el menor pvalue posible distinto de 0
rm(i,gonms,pcnms,X,igos)

heatmap.2(-log(gopcs), trace="none")
heatmap(-log(gopcs[gopcs[,"Neuro"]!=1,1:100])) # solo las neuro

if(FALSE){
  gopcs01 <- ifelse(gopcs<.05,0,1) #hacer una matriz binaria
  i_gos <- which(apply(gopcs01,1,sum)<ncol(gopcs01)) # GOs que aparecen con este filtro
  
  heatmap.2(gopcs01[i_gos,],trace="none")
  
  for(ipca in 1:ncol(gopcs01)){
    igo<-which(gopcs01[,ipca]==0)
    print(Term(names(igo)))
    readline("Press ENTER")     
  }
}

# DISTANCIA ENTRE PCS----
coseno <- cos_sim(-log(gopcs)) #le resto la media?
distcos <- 1-coseno #deberia calcular otra distancia?
loc <- cmdscale(distcos)
x <- loc[, 1]
y <- loc[, 2]
nombres <- colnames(gopcs)
for (i in 1:(length(nombres)-1)){
  nombres[i] <- paste0("PC",nombres[i])
}
nombres[length(nombres)] <- "nerv.syst. development"

summaries <- NULL
for (i in 1:length(nombres)){
  igo <- which(gopcs[,i]<1)
  
  indices_nerv <- names(igo) %in% nervdev_gos
  terms_si <- as.character(Term(names(igo[indices_nerv])))
  terms_si <- paste(terms_si, collapse="</br>")
  terms_si <- paste0("</br>",terms_si)
  
  terms_no <- as.character(Term(names(igo[!indices_nerv])))
  if (length(terms_no)>10){
    terms_no <- paste(terms_no[1:10], collapse="</br>")
  }else{
    terms_no <- paste(terms_no, collapse="</br>")
  }
  terms_no <- paste0("</br><i>",terms_no, "</i>")
  
  go_terms_order <- paste(terms_si, terms_no)
  
  summaries <- c(summaries, go_terms_order)
}
rm(i, igo, go_terms_order, terms_no, terms_si)

# solo el grafico de distancias
plot(x, y, type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE)
text(x, y, nombres, cex = 0.6)

plot(c(x[1:10],x[length(x)]), c(y[1:10], y[length(y)]), type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE)
text(c(x[1:10],x[length(x)]), c(y[1:10], y[length(y)]), c(nombres[1:10],nombres[length(nombres)]), cex = 0.6)

# el grafico con las descripciones de cada PC
library(plotly)
ax <- list(title = "", zeroline = F, showline = F, showticklabels = F, showgrid = F)
plot_ly (x=x,y=y, type='scatter', mode='text',
         text=nombres, hovertext=summaries, hoverinfo='text', showlegend=F, 
         hoverlabel=list(namelength=-1, font=list(size=9)))%>%
  layout(xaxis = ax, yaxis = ax)


h <- hclust(as.dist(distcos),method = "ward.D")
plot(h)


# SEMANTIC SIMILARITY----
library(GOSemSim)
library(igraph)

# OJO CON LA FUNCION godata
mmGO <- godata("org.Mm.eg.db", ont="BP", computeIC=FALSE)
mmGO@IC <- computarIC("org.Mm.eg.db", keytype = "ENTREZID", ont = "BP")
# Arme una funcion de info content porque al crear al environment global,
# no usa GO.db entonces faltan algunos GOIDs

# Por que hay elementos que no existen en Mm, si hgs fue creado usando Mm?

#ademas, sacar cosas que sobran de la lista. Usar unicamente evidencia experimental
# fuente: http://geneontology.org/docs/guide-go-evidence-codes/
evidencia <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP")
ind_ev <- mmGO@geneAnno$EVIDENCE %in% evidencia
go_exp <- mmGO@geneAnno[ind_ev,"GO"]

go_list <- intersect(rownames(gopcs), go_exp)

sim <- mgoSim(go_list,go_list,semData = mmGO, measure ="Resnik",combine=NULL)
sim2 <- mgoSim(go_list,go_list,semData = mmGO, measure ="Wang",combine=NULL)

#solo faltaria ver como calcula Resnik. Ver que dan distinto:
goSim("GO:0007399", "GO:0007399", semData=mmGO, measure="Resnik") #ejemplo: nerv sys dev
mmGO@IC["GO:0007399"]

#normaliza pero no se con que:
diagonal <- diag(sim)
ic_gopcs <- mmGO@IC[go_list]
plot(diagonal,ic_gopcs)

#grafo completo
G <- graph_from_adjacency_matrix(adjmatrix = sim, weighted = TRUE, diag = FALSE, mode = "undirected")
coord <- layout_with_fr(G)
l <- norm_coords(coord, ymin=-1, ymax=1, xmin=-1, xmax=1)
max_weight <- max(E(G)$weight)

lista_nerv <- V(G)$name[V(G)$name %in% nervdev_gos]
colores2 <- c("blue","lightblue")

#si pertenece a neuro tipo=1, sino tipo=2
V(G)$type <- 2
V(G)$type[V(G)$name %in% lista_nerv] <- 1

png("fig/grafo_resnik_completo.png")
plot.igraph(G, layout=l, rescale=FALSE, vertex.size=5, vertex.color=adjustcolor(colores2[V(G)$type], alpha.f=.5),
            vertex.label=NA,edge.width=E(G)$weight/max_weight,
            main="Categorias GO (BP)")
legend('topright', legend=c('Neuro','No neuro'), pch=21, col="black", pt.bg=colores2)
dev.off()


colores3 <- c("blue","gray","lightblue")

for (i in 1:10){ #para cada uno de los 10 primeros pcs
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
  
  fn <- paste0("fig/grafo_resnik_pc",i,".png")
  png(fn)
  plot.igraph(subG, layout=l[ind_coord,], rescale=FALSE, vertex.size=10,
              vertex.color=adjustcolor(colores3[V(subG)$type], alpha.f=.5), vertex.label=NA,
              edge.width=E(subG)$weight/maximo,
              main=paste("PC #",colnames(gopcs)[i]))
  legend('topright', legend=c('Neuro','PC','Ambos'), pch=21, col="black", pt.bg=colores3)
  dev.off()
}
rm(i, lista_pc, subG, maximo, indices_nerv,indices_pc,ind_coord,fn)
