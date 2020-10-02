pqtltools <- function()
# chr:pos in order and one go
{
  chrpos <- unique(with(INF1_merge[ord,c("Chrom","POS")],paste(Chrom,POS,sep=":")))
  pQTLtools::snpqueries(chrpos, catalogue="GWAS", proxies="EUR", p=5e-8, r2=0.8, build=37)
}

ps <- function()
# manually with two lists of SNPs and phenoscanner
{
  rsid <- unique(INF1_merge[ord,"SNP"])
  r1 <- phenoscanner::phenoscanner(rsid[1:100], catalogue="GWAS", proxies="EUR", p=5e-10, r2=0.8, build=37)
  r2 <- phenoscanner::phenoscanner(rsid[101:162], catalogue="GWAS", proxies="EUR", p=5e-8, r2=0.8, build=37)
  list(snps=rbind(r1$snps,r2$snps),results=rbind(r1$results,r2$results))
}

INF <- Sys.getenv("INF")
INF1_merge <- read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE)
ord <- with(INF1_merge,order(CHR,POS))

r1 <- pqtltools()
r2 <- ps()
save(r1,r2,file=file.path(INF,"work","INF1.merge.GWAS"))
