options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_merge <- read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE)

INF1_aggr <- INF1_merge %>%
  select(CHR,POS,prot,Chrom,SNP) %>%
  group_by(CHR,POS,Chrom,SNP) %>%
  summarise(nassocs = n(), prots = paste(prot, collapse = "; ")) %>%
  data.frame()

rsid <- INF1_aggr[["SNP"]]
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37

r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)

save(INF1_aggr,r,file=file.path(INF,"work","INF1.merge.GWAS"))
