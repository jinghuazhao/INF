pqtltools <- function()
  snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)

ps <- function()
# manually with two lists of SNPs and phenoscanner
{
  r1 <- phenoscanner(rsid[1:100], catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  r2 <- phenoscanner(rsid[101:162], catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  list(snps=rbind(r1$snps,r2$snps),results=rbind(r1$results,r2$results))
}

catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37
INF <- Sys.getenv("INF")
INF1_merge <- read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE)
ord <- with(INF1_merge,order(CHR,POS))

library(dplyr)
library(phenoscanner)
library(pQTLtools)
INF1_aggr <- INF1_merge %>%
  select(CHR,POS,prot,Chrom,SNP) %>%
  group_by(CHR,POS,Chrom,SNP) %>%
  summarise(nassocs = n(), prots = paste(prot, collapse = "; ")) %>%
  data.frame()
rsid <- INF1_aggr[["SNP"]]
# chrpos <- with(INF1_aggr,paste(Chrom,POS,sep=":"))

r <- pqtltools()
save(INF1_aggr,r,file=file.path(INF,"work","INF1.merge.GWAS"))
