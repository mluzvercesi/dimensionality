##' @importFrom AnnotationDbi as.list
##' @importFrom GO.db GOBPOFFSPRING
# fuente: https://rdrr.io/bioc/GOSemSim/src/R/computeIC.R

computeIC <- function(OrgDb, keytype = "ENTREZID", ont) {
  ont <- toupper(ont)
  ont <- match.arg(ont, c("BP", "CC", "MF"))

  OrgDb <- load_OrgDb(OrgDb)
  kk <- keys(OrgDb, keytype=keytype)
  goAnno <- suppressMessages(
    select(OrgDb, keys=kk, keytype=keytype,
           columns=c("GO", "ONTOLOGY")))
  
  goAnno <- goAnno[!is.na(goAnno$GO), ]
  goAnno <- goAnno[goAnno$ONTOLOGY == ont,]

  if (!exists(".GOSemSimEnv")) .initial()
  .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
  godata <- get("gotbl", envir=.GOSemSimEnv)
  
  goids <- unique(godata[godata$Ontology == ont, "go_id"])
  ## all GO terms appearing in an given ontology ###########
  goterms=goAnno$GO
  gocount <- table(goterms)
  ## goid of specific organism and selected category.
  goname  <- names(gocount) 
  
  ## ensure goterms not appearing in the specific annotation have 0 frequency..
  go.diff        <- setdiff(goids, goname)
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
  
  ## the probabilities of occurrence of GO terms in a specific corpus.
  p <- cnt/sum(gocount)
  ## IC of GO terms was quantified as the negative log likelihood.
  IC <- -log(p)
  return(IC)
}