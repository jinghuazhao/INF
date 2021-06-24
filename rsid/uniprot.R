uniprotid <- "ACC+ID"
uniprotids <- noquote(gap::inf1[["uniprot"]])
query <- paste(uniprotids,collapse=" ")
for (to in c("ENSEMBL","HGNC","KEGG","UCSC","CHEMBL","DRUGBANK"))
{
  r <- pQTLtools::uniprot2ids(uniprotid,paste0(to,"_ID"),query)
  cat(r,file=paste("INF1.merge",tolower(to),sep="."))
}
