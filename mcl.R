d100 <- as.matrix(dist(pcaAsub$x[,1:100]))

# KNN
kk <- 40
knn_B <- apply(d100,1,function(x){
  saux <- names(sort(rank(x),decreasing=FALSE))[2:(kk+1)]
  return(saux) #devuelve el nombre de k_primeros_vecinos
})
#Calculo el grado de cada nodo (sin armar el grafo todavia)
k<-rep(kk,ncol(knn_B))
names(k)<-colnames(knn_B)
for(i in 1:ncol(knn_B)){
  resto <- colnames(knn_B)[!colnames(knn_B)%in%knn_B[,i]]
  k[i]<-k[i]+sum(as.vector(knn_B[,resto])%in%colnames(knn_B)[i]) #cuantas veces aparece el i en el resto de las columnas
}


#Calculo TOM
res<-c()
for(i in 1:ncol(knn_B)){
  for(j in (i+1):ncol(knn_B)){
    if(j>ncol(knn_B)) next
    if(!(colnames(knn_B)[i]%in%knn_B[,j] | colnames(knn_B)[j]%in%knn_B[,i])) next
    ii<-intersect(knn_B[,i],knn_B[,j])
    if(length(ii)>0){
      a   <- ifelse(colnames(knn_B)[i]%in%knn_B[,j],1,0)
      tom <- (length(ii)+a)/(min(k[i],k[j])-a+1)
      
      res<-rbind(res,c(a=colnames(knn_B)[i],
                       b=colnames(knn_B)[j],
                       w=tom))
    }  
  } 
}
toml<-data.frame(A=res[,1],B=res[,2],weight=as.numeric(res[,3]),stringsAsFactors=FALSE)
#save(toml,file="~/Documents/dimensionality/results/toml.Rdata")


require(igraph)
g<-graph_from_data_frame(toml,directed=FALSE)

#MCL
df <- knn.jac_sub
outfile<- "tmp.txt"

cmd    <- paste("mcl - --abc -I 1.25 -o - >",outfile,sep="")
pw     <- pipe(cmd,open="wb")
write.graph(df,file=pw,format="ncol")
close(pw)

pr<-pipe(paste("cat ",outfile,sep=""),open="r")
lines<-readLines(pr)
comms<-strsplit(split="\t", lines)
close(pr)

names(comms)<-seq_along(comms)
a<-reverseSplit((comms))
membMCL<-as.numeric(unlist(a))
names(membMCL)<-names(a)
table(samples[,"cell_type"],membMCL[samples[,"cell_id"]]) #samples seria mi metadata, pero no se para que sirve esta linea

#Nota: la red tenía 2932 nodos, pero las comunidades tienen 2846. No se si elimina las comunidades de 1
#Los nombres estan bien (llegan hasta 2931)

#comms es una lista (de comunidades) de listas (de genes)
#si quiero un vector de genes con nro de comunidades:
comms_array <- rep(0L, length(unique(purrr::flatten(comms))))
names(comms_array) <- sort(as.array(as.numeric(unique(purrr::flatten(comms)))))

for(nm in names(comms)){
  comms_array[comms[[nm]]] <- as.numeric(nm)
}
rm(nm)
names(comms_array) <- as.numeric(names(comms_array))+1

#----
comunidades <- rep(0, length(com_jac_sub_100))
names(comunidades) <- names(com_jac_sub_100)
for (i in 1:length(comms_mias)){
  comunidades[as.integer(comms_mias[[i]])] <- i
}

# Asortatividad vectorial por comunidad
ncoms <- unique(comunidades[lcc_sub])
comms_array[lcc_sub[as.numeric(names(comms_array))]]


#prueba------
matriz <- cos_sim(cells.es)

a <- ifelse(cells.padj[,!is.nan(apply(cells.nes, 2, min))]<0.15, nes, 0)
a <- ifelse(a>0, a, 0)
matriz <- cos_sim(a)
commslcc <- comms_array[lcc_sub[as.numeric(names(comms_array))]] #vector de comunidades del lcc

for (i in 1:length(comms)){
  v <- as.numeric(names(commslcc[commslcc==i]))
  red <- induced_subgraph(knn_sub, v)
  r <- assortativity_vect(red, matriz[v,v])
  cat("r = ", r, "Comunidad", i, "tamaño", table(commslcc)[i], "\n")
}

#------
#si quiero ponerle un filtro a knn
thres <- (sort((jacc[upper.tri(jacc)])[jacc[upper.tri(jacc)]>0]))[round(0.82*length((jacc[upper.tri(jacc)])[jacc[upper.tri(jacc)]>0])):length((jacc[upper.tri(jacc)])[jacc[upper.tri(jacc)]>0])][1]
jacc_thres <- ifelse(jacc>thres,jacc,0)
knnjacsubthres <-  graph_from_adjacency_matrix(jacc_thres, mode="undirected", weighted = TRUE, diag = FALSE)

#MCL
df <- knnjacsubthres
outfile<- "tmp1.txt"

cmd    <- paste("mcl - --abc -I 1.25 -o - >",outfile,sep="")
pw     <- pipe(cmd,open="wb")
write.graph(df,file=pw,format="ncol")
close(pw)

pr<-pipe(paste("cat ",outfile,sep=""),open="r")
lines<-readLines(pr)
comms_thres<-strsplit(split="\t", lines)
close(pr)

names(comms_thres)<-seq_along(comms_thres)

#comms es una lista (de comunidades) de listas (de genes)
#si quiero un vector de genes con nro de comunidades:
comms_thres_array <- rep(0L, length(unique(purrr::flatten(comms_thres))))
names(comms_thres_array) <- sort(as.array(as.numeric(unique(purrr::flatten(comms_thres)))))

for(nm in names(comms_thres)){
  comms_thres_array[comms_thres[[nm]]] <- as.numeric(nm)
}
rm(nm)
names(comms_thres_array) <- as.numeric(names(comms_thres_array))+1

#----