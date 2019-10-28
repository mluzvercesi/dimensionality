g_users = gph

#Ahora calculo los módulos de alguna forma.
#No puedo usar cluster_optimal() Guimera (2005) porque CRAN ya no supportea GLPK.
membership = membership(cluster_louvain(g_users))#Pruebo con Louvain (rapida).
horizontal_member_mat = as.array(rep(1, length(V(g_users)))) %o% membership #La matriz tiene las memberships en direccion horizontal y las repite en cada fila.
comembership_mat = as.matrix(horizontal_member_mat == t(horizontal_member_mat)) #Si lo comparo con su transpuesta hago una matriz donde la componente (i,j) me dice si i y j pertenecen al mismo modulo.
comembership_mat[which(comembership_mat)] = 1 #Paso los booleanos a valores numericos.
comembership_mat[which(!comembership_mat)] = 0

#Calculo las cantidades que usa Guimera (2005) de forma matricial.
adyacencia = as.matrix(as_adjacency_matrix(g_users))
k_i = rowSums(adyacencia * comembership_mat) #El k_i es el grado del nodo i dentro de su mismo modulo.
k_medio_si = (comembership_mat %*% k_i) / rowSums(comembership_mat) #Calculo la media del k_i para el modulo al cual pertenece cada nodo.
k_var_si = (comembership_mat %*% ((k_i - k_medio_si) ** 2)) / rowSums(comembership_mat) #Lo mismo con la varianza.
z_i = as.numeric((k_i - k_medio_si) / sqrt(k_var_si)) #Finalmente obtengo el z-score (que tan bien conectado esta a su modulo).
z_i[which(is.nan(z_i))] = 0 #Los que tenian varianza cero por ser todos identicos los nodos en el cluster, les pongo cero z-score.

#Ahora planteo una matriz que mapee del nodo i al modulo m.
modu_map = matrix(nrow=length(unique(membership)), ncol=length(V(g_users)))
modu_map[,] = 0 #Seteo default todo a cero.
for(i in 1:length(V(g_users))){
    modu_map[unique(membership)==membership[i],i] = 1 #Senalo a que modulo pertenece el nodo i.
}
module_links = modu_map %*% adyacencia #Si mapeo la adyacencia por los modulos tengo la cantidad de links del nodo i a cada modulo.

p_i = 1 - (colSums(module_links ** 2) / (degree(g_users) ** 2)) #Calculo el p-score de Guimera (2005).
p_i[degree(g_users)==0] <- 0

#Calculo cuales cumplen cada rol.
small_z = which(z_i < 0)
big_z = which(z_i >= 0)
small_p = which(p_i < mean(p_i))
big_p = which(p_i >= mean(p_i))

connector_hubs = big_p[which(big_p %in% big_z)]
connector_nonhubs = big_p[which(big_p %in% small_z)]
provincial_hubs = small_p[which(small_p %in% big_z)]
provincial_nonhubs = small_p[which(small_p %in% small_z)]


plot(p_i,z_i)
abline(v=mean(p_i),lwd=2, lty=2,col='green')
abline(h=0,lwd=2, lty=2,col='green')

# grafico por comunidades
plot(p_i, z_i, col = rainbow(length(unique(membership)))[membership])
abline(v=mean(p_i),lwd=2, lty=2)
abline(h=0,lwd=2, lty=2)

coms1 <- which(table(membership)>10)

for(i in 1:length(coms1)){
  plot(p_i[membership==coms1[i]], z_i[membership==coms1[i]], col = rainbow(length(coms1))[i], xlab = "p", ylab = "z",
       xlim=c(min(p_i),max(p_i)),ylim=c(min(z_i),max(z_i)), main=paste0("C",coms1[i], " (",table(membership)[coms1[i]]," nodos)"))
  abline(v=mean(p_i),lwd=2, lty=2)
  abline(h=0,lwd=2, lty=2)
  readline("Press ENTER")     
}
