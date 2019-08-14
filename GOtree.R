# Funcion para visualizar grafo de GOs
library(igraph)
library(GO.db)

addedge2parents <- function(edges,term){
  parents <- get(term,GOBPPARENTS)
  for (i in 1:length(parents)){
    edge <- c(term,parents[i])
    edges <- rbind(edges,edge)
  }
  return(edges)
}

ipca <- 4
igo <- which(gopcs01[,ipca]==0)
lista <- names(igo)

# d es una matriz de dos columnas: from y to
d <- rep(NULL,2)

#primero agrego solo los parents de la lista
for (termino in lista){
  d <- addedge2parents(d, termino)
}


g = graph_from_data_frame(d, directed = TRUE, vertices = NULL)