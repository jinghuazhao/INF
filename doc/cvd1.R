# 28-9-2018 JHZ

d <- read.delim("cvd1.tsv", as.is=TRUE)
proteins <- d[c("Olink_name", "gene", "Uniprot")]
dim(proteins)
head(proteins, 40)
