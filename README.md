Los scripts numerados son los del pipeline principal:

01_pca.R tiene filtros de QC, PCA y un analisis muy basico de genes importantes en componentes principales.

02_hyperg.R se basa en el resultado de PCA para hacer un test hipergeometrico para sobrerrepresentacion de procesos biologicos GO.
Tiene algunas representaciones visuales (heatmaps; distancias entre PCs).

03_gsea.R se basa en el resultado del test hipergeometrico para hacer un test GSEA.
Tiene distintos filtros para procesos biologicos segun el p-valor y el IC.

04_topologia.R tiene grafo knn, clusterings MCL y tipo celular, asortatividad.

05_analisisparticiones.R compara las particiones (MCL, tipo celular, NES)

###
Otros scripts:

calculate_toproles.R es el analisis de participacion de Guimera.

subconjuntos.R tiene los tipos celulares que nos interesan de cada dataset.

GOfuns.txt es la lista de procesos biologicos que podrian ser relevantes en nuestros experimentos.
