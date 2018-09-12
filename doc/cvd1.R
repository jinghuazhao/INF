# 17-8-2018 JHZ

d <- read.delim("cvd1", as.is=TRUE)
proteins <- d[c("Olink_name", "gene", "Uniprot")]
dim(proteins)
head(proteins, 40)
