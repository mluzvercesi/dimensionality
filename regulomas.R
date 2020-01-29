expresion <- read.table(file = "aracne/MatrizRemoveBatchARACNe.txt", header = TRUE, row.names = 1)

TFnet0 <- read.table(file = "aracne/network.txt", header=TRUE)

tgrm <- TFnet0[,"Target"] %in% unique(TFnet0[,"Regulator"])
TFnet1 <- TFnet0[!tgrm,]

#filtro por MI
TFnet <- TFnet1[TFnet1[,"MI"]>0.5,]

#para usar la base de datos CHEA3, tengo que cambiar los nombres
library(org.Mm.eg.db)

TFnm.ens <-  as.character(unique(TFnet[,"Regulator"]))
TFnm.hugo <- as.character(mapIds(org.Mm.eg.db, keys = TFnm.ens, column="SYMBOL", keytype = "ENSEMBL"))
TFlist <- rep( list(list()), length(TFnm.ens))
names(TFlist) <- TFnm.ens
TFranks <- TFlist
for (nm in names(TFlist)){
  idx <- which(TFnet[,"Regulator"]==nm)
  TFlist[[nm]] <- as.character(TFnet[idx,"Target"])
}
#o sino
TFnm.ens1 <-  as.character(unique(TFnet1[,"Regulator"]))
TFlist1 <- rep( list(list()), length(TFnm.ens1))
names(TFlist1) <- TFnm.ens1
for (nm in names(TFlist1)){
  idx <- which(TFnet1[,"Regulator"]==nm)
  submatrix <- TFnet1[idx,]
  orden <- order(submatrix[,"MI"], decreasing = TRUE)
  TFlist1[[nm]] <- submatrix[orden[1:50],2:3]
}


#creo que no es necesario filtrar por tamaño de regulon, porque el filtro de MI se deshace de los mas grandes
#lo mismo para cantidad de reguladores por target

lista <- as.character(unique(TFnet[,"Target"]))
hugolist <- mapIds(org.Mm.eg.db, keys = lista, column="SYMBOL", keytype = "ENSEMBL")
writeLines(as.character(hugolist), con="lista_hugo.txt")

library(httr)
library(jsonlite)
url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
encode = "json"

#Tomo cada regulon como un gene set, y me fijo que rank tiene el regulador en cada base de datos de CHEA
for (i in length(TFnm.ens):1){
  i <- 8
  geneset <- TFlist[[i]]
  hugo <- as.character(mapIds(org.Mm.eg.db, keys = geneset, column="SYMBOL", keytype = "ENSEMBL"))
  hugo <- hugo[!is.na(hugo)]
  
  if (length(hugo)>1){
    payload = list(query_name = TFnm.hugo[i], gene_set = hugo)
    response = POST(url = url, body = payload, encode = encode) #POST to ChEA3 server
    json = httr::content(response, "text") # CUIDADO! esta funcion es de httr, no de bioconductor
    results_chea = fromJSON(json) #results as list of R dataframes
    
    for (j in 1:length(results_chea)){
      idx <- which(results_chea[[j]][,"TF"]==toupper(TFnm.hugo[i]))
      if (isEmpty(idx)){
        results_chea[[j]] <- paste("El FT", TFnm.hugo[i], "no aparece rankeado")
      }
      else{
        results_chea[[j]] <- results_chea[[j]][idx,]
      }
    }#loop j libraries ChEA3
    TFranks[[i]] <- results_chea
    rm(results_chea,j, payload, response)
  }#if tamaño >1
  else{
    TFranks[[i]] <- paste(TFnm.hugo[i], "tiene < 2 elementos")
  }
  rm(geneset, hugo)
  Sys.sleep(0.4)
}#loop i en regulones

load("aracne/TFranks.RData")

summary_reg <- matrix(rep(0, length(TFranks)*8), nrow = length(TFranks), ncol = 8)
rownames(summary_reg) <- TFnm.hugo
colnames(summary_reg) <- names(TFranks[[26]])

for (i in 1:length(TFnm.ens)){
  if (length(TFranks[[i]])<2){ #si es muy grande o muy chico (<2 elementos)
    summary_reg[i,] <- rep(NA, 8)
  }else{
    for (j in 1:8){
      if (length(TFranks[[i]][[j]])==1){ #si no esta rankeado
        summary_reg[i,j] <- 0
      }else{
        summary_reg[i,j] <- as.numeric(TFranks[[i]][[j]][["Rank"]])
      }
    }
  }
}


TFcorr <- function(target, tflist, expr){
  idx <- unlist(lapply(tflist, function(x){target %in% x}))
  cortargets <- unique(unlist(tflist[idx]))
  
  correlation <- apply(expr[cortargets,], 1, cor, y = as.numeric(expr[target,]))
  # no quiero TODAS las correlaciones, solo entre mi target y los que esten agrupados con el
  
  return(correlation)
}

#funcion de Tomas----
generar_reguloma <- function(grafo){
  reguladores    <- unique(get.edgelist(grafo)[,1]) #Vector de reguladores que están mi grafo
  data           <-  get.edgelist(grafo)
  data           <- as.data.frame(data)
  reguloma <- list()
  for(i in 1:length(reguladores)){
    regulon      <- data[which(reguladores[i] == data$V1), ]
    regulon      <- regulon[which(!(regulon$V2 %in% reguladores)),] #Me saco de encima los FT como targets
    regulon      <- as.character(regulon$V2)
    regulon      <- c(regulon,reguladores[i]) #Le agrego el regulador de ese reguloma.
    reguloma[[i]] <- regulon
  }
  names(reguloma) <- reguladores
  return(reguloma)
}