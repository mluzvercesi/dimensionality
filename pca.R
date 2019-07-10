# FUNCIONES--------------------------------------------------------------
source("~/Documents/funciones.R")
#------------------------------------------------------------------------
# MAIN-------------------------------------------------------------------

dataset <- 3 # (1,2,3) son (A,B,C)

# Para leer los archivos
datasets <- list("GSE95315/GSE95315_10X_expression_data.tab.gz",
                 "GSE95752/GSE95752_C1_expression_data.tab.gz",
                 "GSE104323/GSE104323_10X_expression_data_V2.tab.gz") # en orden, A B y C
file <- paste0("Neurogenesis/Linnarson_NatNeuro2018/",datasets[dataset])
X <- as.matrix(read.table(file, header=TRUE, row.names = 1,as.is=TRUE))

# Para cargarlos
datasets <- list("A","B","C")
load(paste0("~/Documents/dataset",datasets[dataset],".RData"))
load(paste0("~/Documents/results/",datasets[dataset],"subconjunto.RData"))
X <- dataCsub


N_gen  <- dim(X)[1]
M_cell <- dim(X)[2]

# Eliminar celulas con poca expresion
n0 <- as.matrix(apply(X,2,function(x){sum(x==0)}))/N_gen # porcentaje de ceros por columna

# Para elegir porcentaje maximo de 0s por celula aceptado:
hist(n0,main=paste('Dataset',dataset),xlab='Fracción de ceros por célula',sub=paste0('(',N_gen,' genes, ',M_cell,' células)'))
max0s <- c(0.95,0.9,0.95)
#max0s <- c(0.92,0.88,0.95) # para los subconjuntos
max0cell <- max0s[dataset]
abline(v=max0cell,lwd=2, lty=2,col='red')
N_gen*(1-max0cell) # es la cantidad minima de genes que deben expresarse

X_0 <- X[,(n0<max0cell)] # Quedan solo las que no superen el umbral


# Normalizar a 10000 cuentas por celula en A (5000 en C?)
X_0 <- apply(X_0,2,function(x){x*10000/sum(x)})


# Eliminar genes con poca variacion
mean_gen <- apply(X_0,1,function(x){mean(x)})
sdev_gen <- apply(X_0,1,function(x){norm_vect(x-mean(x))})
cv <- sdev_gen/mean_gen
cv[sdev_gen==0] <- 0 # corrijo los ceros que dieron nan
cv2 <- cv^2

# Valor medio de cv2 con bineado logaritmico
cv2mean <- mean_logbin(cv2,mean_gen)

plot(mean_gen,cv2,log="xy",
     xlab=expression(mu),ylab=expression(CV^2),
     main=paste('Dataset',dataset),sub=paste0('(',N_gen,' genes, ',dim(X_0)[2],' células)'))
points(cv2mean,col='green')

# ajuste lineal con bins "centrales"
#dataset A
x <- cv2mean$x[1:65]
y <- cv2mean$y[1:65]
#dataset A sub
x <- cv2mean$x[10:75]
y <- cv2mean$y[10:75]
#dataset B
x <- cv2mean$x[5:65]
y <- cv2mean$y[5:65]
#dataset C
x <- cv2mean$x[5:70]
y <- cv2mean$y[5:70]
#dataset C sub
x <- cv2mean$x[10:70]
y <- cv2mean$y[10:70]

fit_cv2mean <- nls(log(y) ~ a*log(x)+b, start = c(a = -1, b = 10))
y_cv2mean <- predict(fit_cv2mean, data.frame(x=cv2mean$x))

lines(cv2mean$x, exp(y_cv2mean), col="green")
legend('topright', legend=c('data','medios (bin log)','ajuste lineal'),
       col=c('black','green','green'),lty=c(NA,NA,1),pch=c(1,1,NA))


# Me quedo con los genes que esten por encima de este ajuste lineal
pred_lin <- predict(fit_cv2mean, data.frame(x=mean_gen))
filtro_gen <- (cv2>exp(pred_lin))
X_1 <- X_0[filtro_gen,]

Xlog <- X_1
Xlog[X_1 != 0] <- log2(X_1[X_1 != 0])

# PCA -------------------------------------------------------------------
#results_iter <- pca_iter(t(X_1))
#results_svd <- pca_svd(t(Xlog))

#puede ser con variables logaritmicas o escalado con std dev
pca <- prcomp(t(Xlog),retx = TRUE)
#pca <- prcomp(t(X_1),retx = TRUE, scale=TRUE)

lambdas <- (pca$sdev)^2
w <- lambdas/sum(lambdas)
plot(w[1:20],main=paste('Dataset',dataset),xlab="# PC",ylab="% desvio estandar")


# los genes ordenados segun sus proyecciones en los 10 primeros componentes principales
orderedgenes <- genes_by_weight(pca,ncomp=1:10,retw = TRUE)

dropouts <- apply(X_1,1,function(x){sum(x==0)/length(x)})
dropouts <- dropouts[order(dropouts,decreasing=FALSE)]


# proyeccion maxima del primer gen importante en algun PC
max(abs(pca$rotation[match(orderedgenes$genes[1],rownames(Xlog)),]))

# proyeccion maxima del primer PC o de cualquier PC en algun gen
max(pca$rotation[1])
max(apply(pca$rotation,2,max))

# para ver si los dropouts tienen algo que ver (son inversos?) con la importancia
points(dropouts[orderedgenes$genes[1:20]],col="blue")

#------------------------------------------------------------------------
# comparacion filtros de varianza de genes-------------------------------

# Ajuste de cv^2: overdispersion
cv2 <- cv^2
fit_overdisp <- nls(cv2 ~ a/mean_gen + b, start = c(a = 1, b = 1))
# para el grafico, tiene que estar ordenado
x_mu <- seq(min(mean_gen), max(mean_gen), length=length(mean_gen))
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
#------------------------------------------------------------------------