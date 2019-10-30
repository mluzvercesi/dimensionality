# Funcion horrible, faltar definir un monton de cosas, solo funciona en topologia.R

comparar_similitud <- function(x, title, hist_por_comunidad=FALSE){
  #Hipotesis nula: comparo con el azar cambiando el orden de las columnas
  simil <- cos_sim(x)
  
  #Asortatividad vectorial
  print(title)
  for (i in 1:length(ncoms)){
    nombres <- which(com_jac_sub[lcc_sub]==ncoms[i])
    red <- induced_subgraph(knn_sub, nombres)
    r <- assortativity_vect(red, simil[nombres,nombres])
    cat("r = ", r, "Comunidad", ncoms[i], "tamaño", table(com_jac_sub[lcc_sub])[i], "\n")
  }
  
  for (i in 1:length(unique(celltypes_sub))){
    nombres <- which(celltypes_sub[lcc_sub]==unique(celltypes_sub)[i])
    red <- induced_subgraph(knn_sub, nombres)
    r <- assortativity_vect(red, simil[nombres,nombres])
    cat("r = ", r, "Comunidad", unique(celltypes_sub)[i], "tamaño", table(celltypes_sub[lcc_sub])[i],"\n")
  }
  ra <- 0
  N <- 100
  for (j in 1:N){
    azar <- x[,sample(ncol(x))]
    colnames(azar) <- colnames(x)
    ra <- ra + assortativity_vect(knn_sub, cos_sim(azar))
  }
  ra <- ra/N #es casi lo mismo usar la red completa que las comunidades
  cat("r (azar)",ra)
  
  par(mfrow=c(1,1))
  hist(as.dist(simil, diag=FALSE, upper=FALSE),
       main=paste("Red entera"),sub="(todos los pares)", xlab = paste("Similitud de",title))
  
  if (hist_por_comunidad){
    #distribucion (histograma) de similitud entre pares, por comunidad
    pares <- as.dist(simil, diag = FALSE, upper = FALSE)
    par(mfrow=c(2,2))
    for (i in 1:length(com_ind)){
      cuales <- which(com_jac_sub[lcc_sub]==ncoms[com_ind[i]])
      X <- simil[cuales, cuales]
      hist(as.dist(X, diag = FALSE, upper = FALSE), freq=FALSE,
           main=paste("C",ncoms[com_ind[i]],"tamaño",table(com_jac_sub[lcc_sub])[com_ind[i]]), xlab = paste("Similitud de",title))
      
      A <- matrix(0, nrow = N, ncol = (length(X)-dim(X)[1])/2)
      for (j in 1:dim(A)[1]){
        A[j,] <- sample(pares, (length(X)-dim(X)[1])/2)
      }
      hist(A, freq=FALSE, col=rgb(1,0,0,0.5), add=TRUE)
    }
  }
}
