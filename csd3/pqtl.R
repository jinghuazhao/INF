# 19-3-2020 JHZ

catalogue <- Sys.getenv("catalogue")
ps <- function(rsid)
      phenoscanner::phenoscanner(snpquery=rsid, catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
allpQTL <- scan("work/INF1.merge.snp",what="")
knownpQTL <- c("rs2228145","rs704", "rs8064426", "rs113010081", "rs6827617")
rsid <- setdiff(allpQTL,knownpQTL)
m <- length(rsid)
r1 <- ps(rsid[1:100]); lapply(r1,dim)
r2 <- ps(rsid[101:m]); lapply(r2,dim)
snps <- rbind(with(r1,snps),with(r2,snps))
results <- rbind(with(r1,results),with(r2,results))
r <- list(snps=snps,results=results); lapply(r,dim)
save(r,file=paste0("work/INF1.merge.",catalogue))
options(width=500)
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
  s <- subset(results,dataset==d)
  print(s[c("ref_rsid","ref_snpid","rsid","r2","p","trait")])
  sink()
}
