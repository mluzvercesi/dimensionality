expresion <- read.table(file = "aracne/MatrizRemoveBatchARACNe.txt", header = TRUE, row.names = 1)

TFnet0 <- read.table(file = "aracne/network.txt", header=TRUE)

tgrm <- TFnet0[,"Target"] %in% unique(TFnet0[,"Regulator"])
TFnet <- TFnet0[!tgrm,]

#filtro por MI
TFnet <- TFnet[TFnet[,"MI"]>0.5,]

TFlist <- rep( list(list()), length(unique(TFnet[,"Regulator"])))
names(TFlist) <- as.character(unique(TFnet[,"Regulator"]))
for (nm in names(TFlist)){
  idx <- which(TFnet[,"Regulator"]==nm)
  TFlist[[nm]] <- as.character(TFnet[idx,"Target"])
}

#creo que no es necesario filtrar por tamaño de regulon, porque el filtro de MI se deshace de los mas grandes
#lo mismo para cantidad de reguladores por target

TFcorr <- function(target, tflist, expr){
  idx <- unlist(lapply(tflist, function(x){target %in% x}))
  cortargets <- unique(unlist(tflist[idx]))

  correlation <- apply(expr[cortargets,], 1, cor, y = as.numeric(expr[target,]))
  # no quiero TODAS las correlaciones, solo entre mi target y los que esten agrupados con el
  
  return(correlation)
}


#para usar la base de datos CHEA3, tengo que cambiar los nombres
library(org.Mm.eg.db)
lista <- as.character(unique(TFnet[,"Target"]))
lista <- TFlist[[28]]
hugolist <- mapIds(org.Mm.eg.db, keys = lista, column="SYMBOL", keytype = "ENSEMBL")

#hay que sacarle las comillas a los nombres
#sub("-",".",rownames(metadata))
write.table(hugolist, file="listahugo.txt",  sep="\n", row.names = FALSE)

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