# scImpute
# 
# GEO GSE60361
# ERCC spike-ins
# 27.1% zeros
# Fluidigm C1
# Ref 19 (Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq, Zeisel)
# 
# ArrayExpress GSE60361
# Cell cycle
# 22.6% zeros
# Fluidigm C1
# Ref 20
# 
# GEO GSE45719
# mouse embryo
# 61.0% zeros
# custom ?
# Ref 22


file <- ("~/Documents/dimensionality/Neurogenesis/otros/GSE60361_C1-3005-Expression.txt.gz")

tabla <- read.table(file, header=TRUE, as.is=TRUE)
# tiene 2 filas repetidas: "Mar-02" "Mar-01"

genes <- tabla[,"cell_id"]
a <- duplicated(genes)
index <- which(genes %in% genes[a])

genes[index]
# revisar a mano cuales son los duplicados, por ej:
X <- tabla[-18089,]
X <- X[-10773,]
filas <- X[,1]
X <- X[,-1]
rownames(X) <- filas

X <- as.matrix(X)

mar2 <- tabla[10773,-1]
mar1 <- tabla[18089,-1]
mar1 <- array(as.numeric(unlist(mar1)))
mar2 <- array(as.numeric(unlist(mar2)))
X[1261,] <- X[1261,] + mar2
X[4665,] <- X[4665,] + mar1
