# load("~/Documents/results/myworkspace.RData")
# save(object1,object2,file="~/Documents/results/myworkspace.RData")
# require("package") or "package" %in% rownames(installed.packages())

norm_vect <- function(x) {sqrt(sum(x^2))}

#------------------------------------------------------------------------
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
  ## all GO terms appearing in an given ontology ###########
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

#------------------------------------------------------------------------
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

#------------------------------------------------------------------------
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

#------------------------------------------------------------------------
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
		nombre <- list("genes"=rownames(results$rotation)[indices],"weights"=weights[indices])
	}else{
	  nombre <- rownames(results$rotation)[indices]
	}
  return(nombre)
}

#------------------------------------------------------------------------
# FUNCIONES PARA GO
sym2eg <- function(genes,inverse=FALSE){
  # Toma una lista de genes como simbolos (nombres) y devuelve sus egIDs, o viceversa
  # mapIds va de key a columna
  # NOTA: algunos elementos de la lista (mt-Co1, mt-Cytb, mt-Nd1) tienen el simbolo oficial pero 
  # no aparecen asi en la base de datos, ni como alias. Otros elementos (3110021N24Rik) se desconocen.
  
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
  # Devuelve una lista de nombres de todos los genes en genes_id de una categoria GO de BP
  genesBP <- genes_id[genes_id %in% get(GOBPid, org.Mm.egGO2ALLEGS)]
  #genesBPsym <- sym2eg(genesBP,inverse=TRUE)
  #return(list('eg'=genesBP,'symbol'=genesBPsym))
  return(genesBP)
}


#------------------------------------------------------------------------
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

#------------------------------------------------------------------------
#FUNCIONES PARA PROBAR---------------------------------------------------
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

# #------------------------------------------------------------------------
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
