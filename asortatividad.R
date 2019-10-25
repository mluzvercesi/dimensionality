library(igraph)
A <- as.matrix(as_adjacency_matrix(knn_sub))
sim.es <- cos_sim(cells.es)
X <- sim.es


assortativity_vect <- function(network, X){
  # network puede ser la red o la matriz de adyacencia
  # X matriz de similitud o correlacion entre propiedades vectoriales de los nodos
  # devuelve el coeficiente de asortatividad
  
  if (is.igraph(network)){ # si tengo grafo, armo matriz de adyacencia
    A <- as.matrix(as_adjacency_matrix(network))
  }else if (is.matrix(network)){
    A <- network
  }else{stop("No es una red")}
  
  diag(A) <- 0 # sin auto enlaces
  m <-  sum(A)/2 # enlaces totales
  
  if (m==0){
    r<-0
  }else{
    r <- sum(diag(A %*% X)) / (2*m)
  }
  
  # si X es un vector, es decir la propiedad es escalar:
  #  media <- sum(A %*% X)/sum(A)
  #  r <- ((X-media) %*% A %*% (X-media))/sum(A%*% (X-media)^2)
  return(r)
}
