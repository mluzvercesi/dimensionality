library(igraph)
A <- as.matrix(as_adjacency_matrix(knn_sub))
sim.es <- cos_sim(cells.es)
X <- sim.es


#function(network, X)
# red o matriz de adyacencia
# vector o matriz de etiquetas para asortatividad

# si tengo grafo, armo matriz de adyacencia
if (is.igraph(network)){
  A <- as.matrix(as_adjacency_matrix(network))
}else if (is.matrix(network)){
  A <- network
}else{stop("No es una red")}

diag(A) <- 0 # sin auto enlaces
k <- rowSums(A)
m <-  sum(k)/2 #solo es asi si el grafo es no dirigido (la matriz, simetrica)

if (is.matrix(X)){ #???
  n1 <- sum(A %*% X)
  d1 <- k^2 %*% diag(X)
  c <- (t(k) %*% X %*% k)/(2*m)
  r <- (n1 - c)/(d1 - c)
}else{
  # si x es un vector, es decir la propiedad es escalar:
  media <- sum(A %*% x)/sum(A)
  r <- ((x-media) %*% A %*% (x-media))/sum(A%*% (x-media)^2) #funciona
}

# sum( (A - (k %*% t(k))/(2*m)) %*% X )/sum( (diag(k) - (k %*% t(k))/(2*m)) %*% X) #no funciona asi
