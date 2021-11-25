# intomics lookup through https://inbio-discover.com/

genelist <- function(protein)
{
  prot <- subset(gap.datasets::inf1,prot==protein)
  s <- subset(ann,Protein==protein,SYMBOL!=prot[["gene"]])
  c(s[["SYMBOL"]], s[["NEAREST"]])
}

INF <- Sys.getenv("INF")
ann <- read.delim(file.path(INF,"annotate","INF1.merge-annotate.tsv"))

for (protein in c("IL.12B","SCF","TRAIL"))
{
  genes <- genelist(protein)
  write.table(setdiff(genes,"-"),file=paste0(file.path("INF","annotate",protein),".gene"),col.names=FALSE,row.names=FALSE,quote=FALSE)
}
