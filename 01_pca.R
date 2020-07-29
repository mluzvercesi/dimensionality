# Para leer los archivos por primera vez
dataset <- 1 # (1,2,3) son (A,B,C)
datasets <- list("GSE95315/GSE95315_10X_expression_data.tab.gz",
                 "GSE95752/GSE95752_C1_expression_data.tab.gz",
                 "GSE104323/GSE104323_10X_expression_data_V2.tab.gz") # en orden, A B y C
file <- paste0("Neurogenesis/Linnarson_NatNeuro2018/",datasets[dataset])
X <- as.matrix(read.table(file, header=TRUE, row.names = 1,as.is=TRUE))

# O para cargarlos
load("resultados/datasetA.RData")

X <- dataAsub

# genes que aparezcan al menos 20 veces en todas las celulas
igen <- (apply(X,1,sum)>20)
X <- X[igen,]
# genes que aparezcan en menos del 60% de las celulas
igen <- (apply(X,1, function(x){sum(x=!0)})/dim(X)[2])<0.6
X <- X[igen,]

#----
counts <- apply(X, 2, sum) #counts por sample o celula
features <- apply(X, 2, function(x){sum(x!=0)}) # genes por celula
#----
nGen  <- dim(X)[1]
minGen <- 600 # es la cantidad minima de genes que deben expresarse
max0cell <- 1-minGen/nGen

# Eliminar celulas con poca expresion
n0 <- as.matrix(apply(X,2,function(x){sum(x==0)}))/nGen # porcentaje de ceros por celula

# Para ver porcentaje maximo de 0s por celula aceptado:
hist(n0,xlab='Fraccion de ceros por celula',sub=paste0('(',nGen,' genes, ',dim(X)[2],' celulas)'))
abline(v=max0cell,lwd=2, lty=2,col='red')

X0 <- X[,(n0<max0cell)] # Quedan solo las que no superen el umbral


# Eliminar genes con poca variacion
mean_gen <- apply(X0,1,mean)
sdev_gen <- apply(X0,1,sd)
cv <- sdev_gen/mean_gen
cv[sdev_gen==0] <- 0 # corrijo los ceros que dieron nan
cv2 <- cv^2

plot(mean_gen,cv2,log="xy",xlab=expression(mu),ylab=expression(CV^2),
     sub=paste0('(',nGen,' genes, ',dim(X0)[2],' celulas)'))

# ajuste lineal 
fit <- nls(log(cv2) ~ a*log(mean_gen)+b, start = c(a = -1, b = 10))

# Valor medio de cv2 con bineado logaritmico
cv2mean <- mean_logbin(cv2,mean_gen)
points(cv2mean,col='green')
cv2pred <- predict(fit, data.frame(x=cv2mean$x))
lines(cv2mean$x, exp(cv2pred), col="green")
legend('topright', legend=c('data',expression(paste(mu,'(bin)')),'ajuste'),
       col=c('black','green','green'),lty=c(NA,NA,1),pch=c(1,1,NA))

# Me quedo con los genes que esten por encima del ajuste lineal
pred <- predict(fit, data.frame(x=mean_gen))
filtro_gen <- (cv2>exp(pred))
X1 <- X0[filtro_gen,]
# # en Hochgerner: 5000 genes con mas distancia al ajuste
# o <- order(cv2-exp(pred), decreasing = TRUE)
# X5000 <- X0[o[1:5000],]

# Normalizar a 10000 cuentas por celula en dataset A
X1 <- apply(X1,2,function(x){x*10000/sum(x)})
X1 <- round(X1)

Xlog <- X1
Xlog[X1 != 0] <- log2(X1[X1 != 0])

# PCA -----
#results_iter <- pca_iter(t(X1))
#results_svd <- pca_svd(t(Xlog))

#puede ser con variables logaritmicas o escalado con std dev
pca <- prcomp(t(Xlog),retx = TRUE)
#pca <- prcomp(t(X1),retx = TRUE, scale=TRUE) #para probar escalado

lambdas <- (pca$sdev)^2
w <- lambdas/sum(lambdas)
plot(w[1:100],xlab="# PC",ylab="% desvio estandar")
plot(cumsum(w[1:100]),xlab="# PC",ylab="% desvio estandar acumulado")


#analisis muy basico de genes en PCs----
# los genes ordenados segun sus proyecciones en los 10 primeros componentes principales
orderedgenes <- genes_by_weight(pca,ncomp=1:100,retw = TRUE)

# proyeccion maxima del primer gen importante en algun PC
g <- names(orderedgenes)[1]
m <- max(abs(pca$rotation[g,]))
pci <- colnames(pca$rotation)[which.max(abs(pca$rotation[g,]))]
cat('El gen mas importante', g, 'de peso', as.numeric(orderedgenes[1]), 'tiene una proyeccion maxima de', m, 'en el', pci)

# proyeccion maxima del primer PC o de cualquier PC en algun gen
m <- max(abs(pca$rotation[,1]))
pci <- colnames(pca$rotation)[1]
g <- rownames(pca$rotation)[which.max(abs(pca$rotation[,1]))]
cat(pci, 'tiene una proyeccion maxima de', m, 'del gen', g)
m <- max(abs(pca$rotation))
pci <- colnames(pca$rotation)[which.max(apply(abs(pca$rotation),2,max))]
g <- rownames(pca$rotation)[which.max(apply(abs(pca$rotation),1,max))]
cat('La maxima proyeccion de', m, 'es del gen', g, 'sobre el', pci)

par(mfrow=c(2,2)) 
for (i in c(1,10,50,100)){ #los pesos relativos tienen comportamiento similar en todas las PCs
    orderedgenes <- genes_by_weight(pca,ncomp=i,retw = T)
    plot(orderedgenes, xlab="# gen", ylab="Peso", main=paste("PC",i), ylim=c(0,.005))
    text(100, 0.004, round(sum(orderedgenes[1:100]),digits=2),col="red")
    text(500, 0.003, round(sum(orderedgenes[1:500]),digits=2),col="red")
    text(1000, 0.002, round(sum(orderedgenes[1:1000]),digits=2),col="red")
}

par(mfrow=c(1,1))
orderedgenes1 <- genes_by_weight(pca,ncomp=1, retw = T)
plot(orderedgenes1, xlab="# gen", ylab="Peso normalizado", main=paste("PC 1 vs 100"))
orderedgenes100 <- genes_by_weight(pca,ncomp=100, retw = T)
points(orderedgenes100, col="blue")


dropouts <- apply(X1,1,function(x){sum(x==0)/length(x)})
dropouts <- dropouts[order(dropouts,decreasing=FALSE)]

# para ver si los dropouts tienen algo que ver (son inversos?) con la importancia
plot(dropouts[names(orderedgenes)], xlab='indice de gen (en orden de peso total)', ylab='dropout')
plot(dropouts[names(orderedgenes)], orderedgenes, ylab='peso del gen en PCs', xlab='dropout')

#---
# comparacion filtros de varianza de genes----

# Ajuste de cv^2: overdispersion
fit_overdisp <- nls(cv2 ~ a/mean_gen + b, start = c(a = 1, b = 1))
x_mu <- seq(min(mean_gen), max(mean_gen), length=length(mean_gen)) # para el grafico, tiene que estar ordenado
y_overdisp <- predict(fit_overdisp, data.frame(mean_gen=x_mu))

# proyecto PCs con filtro lineal contra PCs con filtro overdispersion
proyecciones <- matrix(0L, nrow = 10, ncol = 10)
for (i in 1:dim(proyecciones)[1]){
  for (j in 1:dim(proyecciones)[2]){
    proyecciones[i,j] <- (results_svd[[2]][,i])%*%(results_overdisp[[2]][,j])
  }
}# las filas son de la lineal, las columnas son de overdisp

# puedo ver si hay correspondencia entre algunos de los PCs
apply(proyecciones,2,function(x){max(abs(x))}) # que tan parecidos son ambos filtros
apply(proyecciones,2,function(x){which(abs(x)==max(abs(x)))}) #a que indices del lineal corresponden
