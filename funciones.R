wd <- getwd()
setwd(paste0(wd,"/Documents/dimensionality"))

#generales----
norm_vect <- function(x) {sqrt(sum(x^2))}

cos_sim <- function(x){
  # calcula matriz de similarity entre columnas de matriz x
  y <- colSums(x^2)
  coseno <- (t(x) %*% x) / (sqrt( y %*% t(y) ))
  coseno[is.nan(coseno)] <- 0
  nombres <- colnames(x)
  colnames(coseno) <- nombres
  rownames(coseno) <- nombres
  return(coseno)
}

mean_logbin <- function(y,x){
  # Valores medios de y con bineado logaritmico, con limites de x
  bins <- seq(log(min(x[x>0])),log(max(x)), length=100)
  bins <- exp(bins)
  bins <- c(0,bins)
  histo <- hist(x, breaks = bins, plot=FALSE)
  
  centros <- rep(0, length(bins)-1)
  y_mean <- rep(0, length(centros))
  # centros de bins
  for (i in 1:length(centros)){
    centros[i] <- (bins[i+1]+bins[i])/2
  }
  # calculo los valores medios en estos bins
  y_order <- y[order(x)]
  k1 <- 1
  k2 <- 0
  for (i in 1:length(y_mean)){
    c <- histo$counts[i]
    if (c==0){
      y_mean[i] <- NaN
    }
    else{
      k2 <- k2 + c
      y_mean[i] <- sum(y_order[k1:k2])/c
      k1 <- k1 + c
    }
  }
  result <- list('x'=centros,'y'=y_mean)
  return(result)
}

pca_svd <- function(X,L=NULL){
  # Funcion para encontrar los L primeros componentes principales de X usando SVD
  # Devuelve desvios standard (raiz de autovalores de t(X)*X) y componentes principales (autovect)
  
  # PARA SC-EXP: que X sea la matriz transpuesta (filas:celulas, columnas:genes)
  
  X <- apply(X,2,function(x){x-mean(x)}) # columnas centradas en 0
  M <- dim(X)[1] # cantidad de filas (celulas)
  N <- dim(X)[2] # cantidad de columnas (genes)
  
  # cantidad de componentes principales: minimo entre lados de la matriz y L
  if (is.null(L)){
    k <- min(N,M)
  }else{
    k <- min(N,M,L)
  }
  
  # Singular value decomposition de la matriz
  d_u_v <- svd(X)
  d <- d_u_v$d # vector de singular values (sigmas)
  u <- d_u_v$u # matriz de MxM cuyas columnas son lsv (autovect de X*Xt)
  v <- d_u_v$v # matriz de NxN cuyas columnas son rsv (autovect de Xt*X)
  
  # Devuelve vector de desvios standard, la matriz v cuyas columnas son PCs y
  # matriz cuyas columnas son las nuevas mediciones rotadas
  result <- list('sdev'=d_u_v$d[1:k],'rotation'=d_u_v$v[,1:k],'x'=d_u_v$u[,1:k])
  return(result)
}

#grafos----
make.knn <- function(expr, k){
  # grafo knn mutuos con correlacion a partir de expr
  # una matriz de expresion donde las columnas son celulas
  require(igraph)
  
  X <- cor(expr) # se podria hacer con similarity en vez de correlacion
  diag(X) <- 0
  
  # knn: ordeno cada fila (modulo), me quedo con los k mayores (k es umbral)
  ordenk <- apply(abs(X),1, function(x){order(x, decreasing = TRUE)})[1:k,]
  # las columnas de ordenk corresponden al orden de las filas de X
  
  A <- matrix(0L, nrow = dim(X)[1], ncol = dim(X)[2])
  rownames(A) <- rownames(X)
  colnames(A) <- colnames(X)
  for (i in 1:dim(X)[1]){
    A[i,ordenk[,i]] <- 1 #X[ordenk[,i],i]
  }
  ady <- (A+t(A))/2
  rm(i, A, X)
  ady <- floor(ady) #matriz de adyacencia de k vecinos cercanos mutuos
  
  #los nodos estan en el mismo orden que en la matriz
  knn <- graph_from_adjacency_matrix(ady, mode="undirected", diag = FALSE)
  return(knn)
}

grafo.simplificado <- function(g, memb, plt = c('none', 'lcc', 'all')){
  # Funcion que toma un grafo y una particion, y genera un nuevo grafo
  # mas chico donde cada nodo es una comunidad, con enlaces pesados no dirigidos
  plt <- match.arg(plt)
  
  if (vcount(g)!=length(memb)){
    interseccion <- intersect(names(memb), vertex_attr(g, name='name'))
    stopifnot(length(interseccion)>0)
    g <- induced_subgraph(g, interseccion)
    memb <- memb[interseccion]
    print('Se utilizo el subconjunto interseccion')
  }
  
  enlaces <- as_edgelist(g)
  ecoms <- apply(enlaces, 2, function(x){memb[x]})
  G <- graph_from_edgelist(ecoms, directed = FALSE)
  E(G)$weight <- 1
  G <- simplify(G, edge.attr.comb=list(weight="sum"), remove.loops=FALSE)
  
  #la transparencia es la cantidad de loops
  alphas <- rep(0, vcount(G))
  idx <- as_edgelist(G)[which_loop(G),1]
  alphas[idx] <- log(E(G)$weight[which_loop(G)]+1)
  alphas <- alphas/max(alphas)
  V(G)$alpha <- alphas
  V(G)$size <- as.numeric(table(memb))
  
  if (plt=='all'){
    plot(
      simplify(G, remove.loops=TRUE), # para que no grafique loops
      vertex.size     = log(V(G)$size, base=min(V(G)$size)), # el tamanio de la comunidad
      vertex.label    = NA,
      edge.arrow.size = .25,
      vertex.color    = rgb(0,.5,.5, alpha = V(G)$alpha)
    )
  } else if (plt=='lcc'){
    lcc <- decompose(G, max.comps = 1, min.vertices = 2)
    if (length(lcc)==0){
      print('No hay lcc con mas de un nodo')
    }else{
      lcc <- lcc[[1]]
      plot(
        simplify(lcc, remove.loops=TRUE),
        vertex.size     = log(V(lcc)$size, base=min(V(lcc)$size)), # el tamanio de la comunidad
        vertex.label    = NA,
        edge.arrow.size = .25,
        vertex.color    = rgb(0,.5,.5, alpha = V(lcc)$alpha)
      )
    } #end if hay lcc
  } #end if plt='lcc'
  
  return(G)
}

vecinoscomun <- function(g, par){
  vec <- adjacent_vertices(g, v=par)
  inter <- length(intersect(vec[[1]],vec[[2]]))
  union <- length(union(vec[[1]],vec[[2]]))
  if(union==0){
    return(0)
  }else{
    return(inter/union)
  }
}

assortativity_vect <- function(network, X){
  # network puede ser la red o la matriz de adyacencia
  # X matriz de similitud o correlacion entre propiedades vectoriales de los nodos
  # devuelve el coeficiente de asortatividad
  
  if (is.igraph(network)){ # si tengo grafo, armo matriz de adyacencia
    A <- as.matrix(as_adjacency_matrix(network))
  }else if (is.matrix(network)){
    A <- network
  }else{stop("No se reconoce como red")}
  
  diag(A) <- 0 # sin auto enlaces
  m <-  sum(A)/2 # enlaces totales
  
  if (m==0){
    r<-0
  }else{
    X <- X[colnames(A),colnames(A)]
    r <- sum(diag(A %*% X)) / (2*m)
  }
  
  # si X es un vector, es decir la propiedad es escalar:
  #  media <- sum(A %*% X)/sum(A)
  #  r <- ((X-media) %*% A %*% (X-media))/sum(A%*% (X-media)^2)
  return(r)
}


indice.rand <- function (group1, group2){
  # Devuelve:
  # r = indice Rand
  # a = cantidad de pares juntos en ambas particiones y su fraccion
  # b = cantidad de pares separados en ambas particiones y su fraccion
  # cd= cantidad de pares juntos en un grupo y separados en el otro y su fraccion
  if (length(group1)==length(group2)){
    x <- abs(sapply(group1, function(x) x - group1))
    x[x > 1] <- 1
    y <- abs(sapply(group2, function(x) x - group2))
    y[y > 1] <- 1
    cuentas <- table(x+y)
    a <- (cuentas[[1]]-length(group1))/2 #cantidad de 0s sin diagonal
    b <- cuentas[[3]]/2 #cantidad de 2s
    cd <- cuentas[[2]]/2 #cantidad de 1s
    n <- a+b+cd
    res <- list('r'=(a+b)/n,'a'=c(a, a/n),'b'=c(b, b/n),'cd'=c(cd, cd/n))
    return(res)
  }else{
    print("Los grupos son diferentes")
  }
}

concatparticiones <- function(x, particion){
  # x es un vector o una matriz en la cual
  # las filas son vertices y cada columna es una particion.
  # particion es un vector o una lista de vectores para agregar a x
  # (simplifica guardar nodos para gephi)
  
  x <- as.data.frame(x)
  
  oldcolnms <- colnames(x)
  
  if (typeof(particion)!="list"){
    nm <- deparse(substitute(particion))
    particion <- list(particion)
    names(particion) <- nm
  }
  if (is.null(names(particion))){
    names(particion) <- paste0('part',seq(length(particion)))
  }
  for (z in particion){
    if(sum(names(z) %in% rownames(x))!=length(z)){
      print("Algunos elementos de la particion no se encuentran en la matriz y seran descartados")
      z <- z[names(z) %in% rownames(x)]
    }
    x <- cbind(x, rep(NA, n=nrow(x)))
    x[names(z),ncol(x)] <- z
    faltan <- seq(sum(is.na(x[,ncol(x)]))) + max(z)
    x[is.na(x[,ncol(x)]),ncol(x)] <- faltan
  }
  colnames(x) <- c(oldcolnms,names(particion))
  return(x)
}

#genes----
genes_by_weight <- function(results,ncomp=1,retw=FALSE){ 
  # Devuelve una lista ordenada de nombres de genes segun su peso sobre los componentes principales ncomp de X
  # Si retw = TRUE devuelve el peso ("importancia") de cada gen
  # Necesita matriz de scexp (X) y resultado de pca_svd o de prcomp (results)
  p <- results$rotation[,ncomp]
  lambdas <- (results$sdev)^2
  w <- lambdas[ncomp]/sum(lambdas)
  if (is.null(dim(p))){ # si solo hay un componente principal
    p <- p*w
    weights <- abs(p)
  }else{ # si hay mas de uno, es matriz
    p <- apply(p,1,function(x){x*w})
    weights <- apply(p,2,function(x){sum(abs(x))})
  }
  weights <- weights/sum(weights)
  # los pesos pueden variar entre prcomp y pca_svd, pero el orden es el mismo
  
  indices <- order(weights,decreasing = TRUE) # orden de maximas proyecciones
  
  if (retw){
    nombre <- weights[indices]
  }else{
    nombre <- rownames(results$rotation)[indices]
  }
  return(nombre)
}

##' @importFrom AnnotationDbi as.list
##' @importFrom GO.db GOBPOFFSPRING
# fuente: GOSemSim https://rdrr.io/bioc/GOSemSim/src/R/computeIC.R
computarIC <- function(OrgDb, keytype = "ENTREZID", ont) {
  ont <- toupper(ont)
  ont <- match.arg(ont, c("BP", "CC", "MF"))
  
  OrgDb <- load_OrgDb(OrgDb)
  kk <- keys(OrgDb, keytype=keytype)
  goAnno <- suppressMessages(
    select(OrgDb, keys=kk, keytype=keytype,
           columns=c("GO", "ONTOLOGY")))
  
  goAnno <- goAnno[!is.na(goAnno$GO), ]
  goAnno <- goAnno[goAnno$ONTOLOGY == ont,]
  
  goids <- select(GO.db, keys="BP", columns=c("GOID"), keytype="ONTOLOGY")
  goids <- goids$GOID
  ## all GO terms appearing in an given ontology ###
  goterms=goAnno$GO
  gocount <- table(goterms)
  ## goid of specific organism and selected category.
  goname  <- names(gocount) 
  
  ## ensure goterms not appearing in the specific annotation have 0 frequency
  go.diff        <- setdiff(goids, goname) #esta en goids pero no en goname
  m              <- double(length(go.diff))
  names(m)       <- go.diff
  gocount        <- as.vector(gocount)
  names(gocount) <- goname
  gocount        <- c(gocount, m)
  
  Offsprings <- switch(ont,
                       MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                       BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                       CC = AnnotationDbi::as.list(GOCCOFFSPRING))
  
  cnt <- gocount[goids] + sapply(goids, function(i) sum(gocount[Offsprings[[i]]], na.rm=TRUE))
  names(cnt) <- goids
  # CORRECCION DE ESTAS DOS ULTIMAS LINEAS:
  #cnt <- gocount + sapply(names(gocount), function(i) sum(gocount[Offsprings[[i]]], na.rm=TRUE))
  #names(cnt) <- names(gocount)
  
  ## the probabilities of occurrence of GO terms in a specific corpus.
  p <- cnt/sum(gocount)
  ## IC of GO terms was quantified as the negative log likelihood.
  IC <- -log(p)
  return(IC)
}


#GO----
sym2eg <- function(genes,inverse=FALSE){
  # Toma una lista de genes como simbolos (nombres) y devuelve sus egIDs, o viceversa
  # mapIds va de key a columna
  # NOTA: algunos elementos de la lista (mt-Co1, mt-Cytb, mt-Nd1) tienen el simbolo oficial pero 
  # no aparecen asi en la base de datos, ni como alias. Otros elementos (3110021N24Rik) se desconocen.
  require(org.Mm.eg.db)
  
  if (inverse){
    genes_t <- mapIds(org.Mm.eg.db, keys=genes, column=c("ALIAS"), keytype="ENTREZID")
  }else{
    genes_t <- mapIds(org.Mm.eg.db, keys=genes, column=c("ENTREZID"), keytype="ALIAS")
    
    # algunos nombres son secuencias:
    genes_na <- names(genes_t[is.na(genes_t)])
    genes_na_id <- try(mapIds(org.Mm.eg.db, keys=genes_na, column=c("ENTREZID"), keytype="ACCNUM"),silent=TRUE)
    genes_t[names(genes_na_id)] <- genes_na_id
    
    #paste0(substring("FRMPD2",1,1),tolower(substring("FRMPD2",2)))
  }
  return(genes_t)
}

eg2go <- function(genes_id){
  # Toma una lista de genes como egIDs y devuelve sus GO IDs
  genes_go <- mapIds(org.Mm.eg.db, keys=genes_id, column=c("GO"), keytype="ENTREZID")
  return(genes_go)
}

genesinbp <- function(genes_id,GOBPid){
  # Devuelve una lista de todos los genes de genes_id que estan en GO BP
  genesBP <- genes_id[genes_id %in% get(GOBPid, org.Mm.egGO2ALLEGS)]
  #genesBPsym <- sym2eg(genesBP,inverse=TRUE)
  #return(list('eg'=genesBP,'symbol'=genesBPsym))
  return(genesBP)
}


#Fuera de uso------------------------------------------------------------
# pca_1c <- function(X,tol=0.001){
#   # Funcion para encontrar el primer componente principal de X
#   # Devuelve el maximo autovalor y su autovector (comp principal)
#   N <- dim(X)[1]
#   M <- dim(X)[2]
#   error <- 1 # error inicial
#   
#   # Resta la media a cada variable X[i,]
#   X <- apply(X,1,function(x){x-mean(x)})
#   X <- t(X)
#   
#   r <- runif(M, -1.0, 1.0) # vector random de longitud M
#   r <- r/norm_vect(r) # vector random normalizado
#   
#   # Maximiza la suma de los cuadrados de las proyecciones
#   while (error>tol){
#     s <- rep(0,M)
#     for (i in 1:N){
#       s <- s + (X[i,]%*%r) * X[i,]
#     }
#     ev <- r%*%s
#     error <- norm_vect(ev*r - s)
#     r <- s/norm_vect(s)
#   }
#   result <- list("autovalor"=ev,"autovector"=r)
#   return(result)
# }
# 
# pca_iter <- function(X,tol=0.001,L=10){
#   # Funcion iterativa para encontrar los L primeros componentes principales de X
#   # Devuelve sigmas (raices de autovalores) y autovectores como columnas de matriz
#   # Necesita funcion pca_1c
#   N <- dim(X)[1] # cantidad de variables (filas)
#   M <- dim(X)[2] # cantidad de mediciones (columnas)
#   
#   # cantidad de componentes principales, puede ser como maximo el minimo lado de la matriz
#   k <- min(N,M,L)
#   
#   # Autovalores y autovectores
#   aval <- rep(0,k)
#   avec <- matrix(0L, nrow = M, ncol = k) #las columnas son los autovectores
#   
#   # Iteracion
#   for (i in 1:k){
#     C <- pca_1c(X,tol)
#     aval[i] <- C[[1]]
#     avec[,i] <- C[[2]]
#     
#     # Saca componentes de X de manera sucesiva
#     X <- X - X%*%C[[2]]%*%C[[2]]
#   }
#   # Devuelve un vector de desvios standard y la matriz cuyas L columnas son PCs
#   result <- list("sdev"=sqrt(aval),"x"=avec)
#   return(result)
# }

# PARA PROBAR------------------------------------------------------------
# pca_forzado <- function(X,nombre){
#   # Funcion que permite forzar un gen "nombre" a aparecer en PCA como PC
#   # Devuelve las componentes principales como columnas de la matriz pcs y la lista de genes importantes
#   # Usa funcion pca_svd y genes_proy
#   
#   # Si es uno solo lo transformo en componente principal
#   ind_gen = match(nombre, rownames(X)) # indice del gen forzado
#   pc1 <- rep(0, dim(X)[1]) # gen forzado
#   pc1[ind_gen] <- 1
#   
#   sdev1 <- sqrt(sum(X[ind_gen,]^2))
#   
#   #simplifico esta cuenta para el espacio ortogonal: X - X*pc1*pc1
#   X_reducida <- X
#   X_reducida[ind_gen,] <- 0
#   # Saco el gen de la matriz, aplico PCA al resto
#   results_reducida <- pca_svd(t(X_reducida))
#   
#   sdev <- c(sdev1,results_reducida$sdev)
#   pcs <- cbind(pc1,results_reducida$rotation)
#   lista <- c(nombre, genes_proy(X_reducida,results_reducida,10)$genes)
#   
#   nombre <- list('sdev'=sdev,'rotation'=pcs, 'genes'=lista)
#   
#   # Si es una lista de varios genes. Opciones: 
#   # - cada uno es una PC
#   # - se pueden combinar: hacer PCA en el subespacio de la lista
#   return(nombre)
# }

# pca_svd_zeros <- function(X,L=NULL){
#   # Funcion para encontrar los L primeros componentes principales de X usando SVD
#   # Devuelve desvios standard (raiz de autovalores de t(X)*X) y componentes principales (autovect)
#   
#   # Forzando que se mantengan los ceros
# 
#   # PARA SC-EXP: que X sea la matriz transpuesta (filas:celulas, columnas:genes)
#   zeros <- (X==0)
#   X <- apply(X,2,function(x){x-mean(x)}) # columnas centradas en 0
#   X[zeros] <- 0
#   M <- dim(X)[1] # cantidad de filas (celulas)
#   N <- dim(X)[2] # cantidad de columnas (genes)
#   
#   # cantidad de componentes principales: minimo entre lados de la matriz y L
#   if (is.null(L)){
#     k <- min(N,M)
#   }else{
#     k <- min(N,M,L)
#   }
#   
#   # Singular value decomposition de la matriz
#   d_u_v <- svd(X)
#   d <- d_u_v$d # vector de singular values (sigmas)
#   u <- d_u_v$u # matriz de MxM cuyas columnas son lsv (autovect de X*Xt)
#   v <- d_u_v$v # matriz de NxN cuyas columnas son rsv (autovect de Xt*X)
#   
#   # Devuelve vector de desvios standard, la matriz v cuyas columnas son PCs y
#   # matriz cuyas columnas son las nuevas mediciones rotadas
#   result <- list('sdev'=d_u_v$d[1:k],'rotation'=d_u_v$v[,1:k],'x'=d_u_v$u[,1:k])
#   return(result)
# }
# 
# pca_svd_min <- function(X,L=NULL){
#   # Funcion para encontrar los L primeros componentes principales de X usando SVD
#   # Devuelve desvios standard (raiz de autovalores de t(X)*X) y componentes principales (autovect)
#   
#   # PARA SC-EXP: que X sea la matriz transpuesta (filas:celulas, columnas:genes)
# 
#   # Forzando que los minimos sean ceros. Hay que reescalar despues?
#   
#   X <- apply(X,2,function(x){x-min(x)}) # columnas llevadas al 0
#   M <- dim(X)[1] # cantidad de filas (celulas)
#   N <- dim(X)[2] # cantidad de columnas (genes)
#   
#   # cantidad de componentes principales: minimo entre lados de la matriz y L
#   if (is.null(L)){
#     k <- min(N,M)
#   }else{
#     k <- min(N,M,L)
#   }
#   
#   # Singular value decomposition de la matriz
#   d_u_v <- svd(X)
#   d <- d_u_v$d # vector de singular values (sigmas)
#   u <- d_u_v$u # matriz de MxM cuyas columnas son lsv (autovect de X*Xt)
#   v <- d_u_v$v # matriz de NxN cuyas columnas son rsv (autovect de Xt*X)
#   
#   # Devuelve vector de desvios standard, la matriz v cuyas columnas son PCs y
#   # matriz cuyas columnas son las nuevas mediciones rotadas
#   result <- list('sdev'=d_u_v$d[1:k],'rotation'=d_u_v$v[,1:k],'x'=d_u_v$u[,1:k])
#   return(result)
# }
