A <- as.matrix(as_adjacency_matrix(knn_sub))
X <- sim.es

#function(network, X)
# red o matriz de adyacencia
# vector o matriz de etiquetas para asortatividad

# si tengo grafo, armo matriz de adyacencia
if (is.igraph(network)){
  A <- as.matrix(as_adjacency_matrix(network))
}else if (is.matrix(network)){
  A <- network
}

diag(A) <- 0 # sin loops
k <- rowSums(A)
m <-  sum(k)/2

if (is.matrix(X)){ #???
  n1 <- sum(A %*% X)
  d1 <- k %*% diag(X)
  c <- (t(k) %*% X %*% k)/(2*m)
}else{
  #d1 <- 
}

r <- (n1 - c)/(d1 - c)
