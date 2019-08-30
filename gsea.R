load("~/Documents/dimensionality/results/Asubconjunto.RData")
library(fgsea)
#------------------------------------------------------------------------
# Pre procesamiento
load(paste0("~/Documents/dimensionality/results/Asubconjunto.RData"))
rm(pcaAsub,pcaAsub_scaled)
X <- dataAsub

N_gen  <- dim(X)[1]
M_cell <- dim(X)[2]

#max0s <- c(0.92,0.88,0.95) # para los subconjuntos
X_0 <- X[,(n0<0.92)] # Quedan solo las que no superen el umbral

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
#dataset A sub
x <- cv2mean$x[10:75]
y <- cv2mean$y[10:75]

fit_cv2mean <- nls(log(y) ~ a*log(x)+b, start = c(a = -1, b = 10))
y_cv2mean <- predict(fit_cv2mean, data.frame(x=cv2mean$x))
lines(cv2mean$x, exp(y_cv2mean), col="green")

# Me quedo con los genes que esten por encima de este ajuste lineal
pred_lin <- predict(fit_cv2mean, data.frame(x=mean_gen))
filtro_gen <- (cv2>exp(pred_lin))

rm(cv,cv2,cv2mean,filtro_gen,fit_cv2mean,mean_gen,mean_logbin,pred_lin,sdev_gen,x,y,y_cv2mean)

X_1 <- X_0[filtro_gen,]
N_gen  <- dim(X_1)[1]
M_cell <- dim(X_1)[2]

Xlog <- X_1
Xlog[X_1 != 0] <- log2(X_1[X_1 != 0])

Z <- apply(Xlog, 1, function(x){(x-mean(x))/norm_vect(x-mean(x))})
Z <- t(Z) #filas = genes; columnas = celulas

#------------------------------------------------------------------------
# GSEA
# Cada columna de la matriz Z es una celula, es decir un vector de ranks
# save(Z,gopcsA,file="~/Documents/dimensionality/results/gseaA.RData")
load("~/Documents/dimensionality/results/gseaA.RData")

zranks <- Z[,1]

gonames <- rownames(gopcsA)
Ngo <- length(gonames)

# falta hacer una lista (de goids) de listas (de genes de cada categoria)
goids <- gopcsA[,1]

gsea_res <- fgsea(pathways = goids, stats = zranks, nperm = 10000)
