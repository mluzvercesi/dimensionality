datasets <- list("GSE95315/GSE95315_metadata.txt",
                 "GSE95752/GSE95752_metadata.txt",
                 "GSE104323/GSE104323_metadata_barcodes_24185cells.txt") # en orden, A B y C
#---------------
#para dataset A:
dataset <- 1 # (1,2,3) son (A,B,C)
metadata <- read.csv(paste0("~/Documents/Neurogenesis/Linnarson_NatNeuro2018/",datasets[dataset]), 
                     sep="\t", header=TRUE, row.names=1)

rownames(metadata) <- sub("-",".",rownames(metadata))
rownames(metadata) <- paste0("X",rownames(metadata))

table(metadata[,"cell_type"])

# Uno cada categoria glia-like con la mas diferenciada
metadata$cell_type <- sub(" glia-like","",metadata$cell_type)
# excepto los neuroblast 1 y 2, uno todos los que tengan un mismo nombre con != nro
metadata$cell_type[!grepl("Neuroblast", metadata$cell_type)] <- sub("\\d+$", "", metadata$cell_type[!grepl("Neuroblast", metadata$cell_type)])

subconjunto <- rbind(metadata[metadata$cell_type=="nIPC",],
                     metadata[metadata$cell_type=="Neuroblast",],
                     metadata[metadata$cell_type=="Neuroblast2",],
                     metadata[metadata$cell_type=="Granule-immature",],
                     metadata[metadata$cell_type=="Granule-mature",]
                     )

table(metadata[,"cell_type"])
table(subconjunto[,"cell_type"])

dataAsub <- dataA[,rownames(subconjunto)]

#---------------
#para dataset B:
dataset <- 2
metadata <- read.csv(paste0("~/Documents/Neurogenesis/Linnarson_NatNeuro2018/",datasets[dataset]), 
                     sep="\t", header=TRUE, row.names=1)

rownames(metadata) <- paste0("X",rownames(metadata))

table(metadata[,"cell_type"])

subconjunto <- rbind(metadata[metadata$cell_type=="nIPC",],
                     metadata[metadata$cell_type=="Neuroblast",],
                     metadata[metadata$cell_type=="Granule-immature",],
                     metadata[metadata$cell_type=="Granule-mature",]
)
table(subconjunto[,"cell_type"])

dataBsub <- dataB[,rownames(subconjunto)]

#---------------
#para dataset C:
dataset <- 3
metadata <- read.csv(paste0("~/Documents/Neurogenesis/Linnarson_NatNeuro2018/",datasets[dataset]), 
                     sep="\t", header=TRUE, row.names=NULL)
metadata <- metadata[!apply(metadata=="", 1, all),] # para eliminar las ultimas 30 filas vacias

rownames(metadata) <- metadata$Sample.name..24185.single.cells.
rownames(metadata) <- sub("-",".",rownames(metadata))
rownames(metadata) <- paste0("X",rownames(metadata))

table(metadata[,"characteristics..cell.cluster"])

subconjunto <- rbind(metadata[metadata$characteristics..cell.cluster=="nIPC",],
                     metadata[metadata$characteristics..cell.cluster=="Neuroblast",],
                     metadata[metadata$characteristics..cell.cluster=="Immature-GC",],
                     metadata[metadata$characteristics..cell.cluster=="GC-juv",],
                     metadata[metadata$characteristics..cell.cluster=="GC-adult",],
                     metadata[metadata$characteristics..cell.cluster=="nIPC-perin",])
table(subconjunto[,"characteristics..cell.cluster"])

dataCsub <- dataC[,rownames(subconjunto)]
