# 19-10-2019 JHZ

# m <- read.delim("work/INF1.merge",as.is=TRUE)[c("prot","MarkerName")]
# names(m) <- c("prot","ref_snpid")

require(phenoscanner)
catalogue <- Sys.getenv("catalogue")
rsid <- with(read.table("work/INF1.merge.rsid",as.is=TRUE,col.names=c("snpid","rsid")), rsid)
r1 <- phenoscanner(snpquery=rsid[1:100], catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
lapply(r1,dim)
r2 <- phenoscanner(snpquery=rsid[101:162], catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
lapply(r2,dim)
r <- list(snps=rbind(with(r1,snps),with(r2,snps)),results=rbind(with(r1,results),with(r2,results)))
lapply(r,dim)
save(r,file=paste0("work/INF1.merge.",catalogue))
options(width=500)
load(paste0("work/INF1.merge.",catalogue))
attach(r)
results <- within(results,{
   a1 <- ref_a1
   a2 <- ref_a2
   swap <- ref_a1 > ref_a2
   a1[swap] <- ref_a2[swap]
   a2[swap] <- ref_a1[swap]
   ref_snpid <- paste0(ref_hg19_coordinates,":",a1,"_",a2)
})
for(d in unique(with(results,dataset)))
{
  cat(d,"\n")
  sink(paste(catalogue,d,sep="."))
  s <- subset(results[c("ref_rsid","ref_snpid","rsid","r2","p","trait","dataset","pmid")],dataset==d)
  print(s)
  sink()
}
detach(r)
