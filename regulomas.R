expresion <- read.table(file = "aracne/MatrizRemoveBatchARACNe.txt", header = TRUE, row.names = 1)

TFnet0 <- read.table(file = "aracne/network.txt", header=TRUE)

tgrm <- TFnet0[,"Target"] %in% unique(TFnet0[,"Regulator"])
TFnet <- TFnet0[!tgrm,]

TFlist <- rep( list(list()), length(unique(TFnet[,"Regulator"])))
names(TFlist) <- as.character(unique(TFnet[,"Regulator"]))

for (nm in names(TFlist)){
  idx <- which(TFnet[,"Regulator"]==nm)
  TFlist[[nm]] <- as.character(TFnet[idx,"Target"])
}

TFcorr <- function(target, expr){
  idx <- which(TFnet[,"Target"]==target)
  regs <- as.character(TFnet[idx,"Regulator"])
  targets <- unique(as.character(unlist(TFlist[regs])))
  
  correlation <- apply(expr[targets,], 1, cor, y = as.numeric(expresion[target,]))
  # no quiero TODAS las correlaciones, solo entre mi target y los que esten agrupados con el
  
  return(correlation)
}
