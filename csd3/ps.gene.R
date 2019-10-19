# 19-10-2019 JHZ

options(width=500)
glist <- scan("work/INF1.merge.gene",what="")
genes <- data.frame()
results <- data.frame()
require(phenoscanner)
for(s in 1:6)
{
  gset <- ifelse(s==6, 61:68, (s-1)*10+1:10)
  g <- phenoscanner(genequery=glist[gset], catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  genes <- rbind(genes,with(g,genes))
  results <- rbind(results,with(g,results))
}
r <- list(genes=genes,results=results)
lapply(r,dim)
save(r,file="work/INF1.merge.genes")
attach(r)
results <- within(results,{
   ref_a1 <- a1
   ref_a2 <- a2
   a1 <- ref_a1
   a2 <- ref_a2
   swap <- ref_a1 > ref_a2
   a1[swap] <- ref_a2[swap]
   a2[swap] <- ref_a1[swap]
   ref_snpid <- paste0(hg19_coordinates,":",a1,"_",a2)
})
sink("gene.genes")
genes
sink()
for(d in unique(with(results,dataset)))
{
  cat(d,"\n")
  sink(paste("gene",d,sep="."))
  s <- subset(results[c("rsid","ref_snpid","rsid","p","n","trait","dataset","pmid")],dataset==d)
  print(s)
 sink()
}
detach(r)
