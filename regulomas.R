expresion <- read.table(file = "aracne/MatrizRemoveBatchARACNe.txt", header = TRUE, row.names = 1)

TFnet0 <- read.table(file = "aracne/network.txt", header=TRUE)

tgrm <- TFnet0[,"Target"] %in% unique(TFnet0[,"Regulator"])
TFnet <- TFnet0[!tgrm,]

#filtro por MI
TFnet <- TFnet[TFnet[,"MI"]>0.5,]

TFnm.ens <-  as.character(unique(TFnet[,"Regulator"]))
TFnm.hugo <- as.character(mapIds(org.Mm.eg.db, keys = TFnm.ens, column="SYMBOL", keytype = "ENSEMBL"))
TFlist <- rep( list(list()), length(TFnm.ens))
names(TFlist) <- TFnm.ens
TFranks <- TFlist
for (nm in names(TFlist)){
  idx <- which(TFnet[,"Regulator"]==nm)
  TFlist[[nm]] <- as.character(TFnet[idx,"Target"])
}

#creo que no es necesario filtrar por tamaño de regulon, porque el filtro de MI se deshace de los mas grandes
#lo mismo para cantidad de reguladores por target

#para usar la base de datos CHEA3, tengo que cambiar los nombres
library(org.Mm.eg.db)
url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
encode = "json"

#Tomo cada regulon como un gene set, y me fijo que rank tiene el regulador en cada base de datos de CHEA
for (i in length(TFnm.ens):1){
  #i <- 75
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
    }
    TFranks[[i]] <- results_chea
  }
}

lista <- as.character(unique(TFnet[,"Target"]))
hugolist <- mapIds(org.Mm.eg.db, keys = lista, column="SYMBOL", keytype = "ENSEMBL")
writeLines(as.character(hugo), con="lista_hugo.txt")

library(httr)
library(jsonlite)

genes = c("SMAD9","FOXO1","MYC","STAT1",'STAT3',"SMAD3") #ejemplo
genes <- as.character(hugolist)
url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes)

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = httr::content(response, "text") # CUIDADO! esta funcion es de httr, no de bioconductor

#results as list of R dataframes
results_chea = fromJSON(json)


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