# 21-3-2020 JHZ

catalogue <- Sys.getenv("catalogue")
ps <- function(rsid)
      phenoscanner::phenoscanner(snpquery=rsid, catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
allpQTL <- scan("work/INF1.merge.snp",what="")
sun18 <- c("rs2228145","rs704", "rs8064426", "rs113010081", "rs6827617", "rs635634", "rs579459", "rs2247769", "rs4734879", "rs9469127",
           "rs6993770", "rs705379", "rs757973", "rs13229619", "rs6921438", "rs9272226", "rs11759846", "rs10733789", "rs10822155",
           "rs7090111", "rs16840522")
rsid <- setdiff(allpQTL,sun18)
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
   ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
})
for(d in unique(with(results,dataset)))
{
  cat(d,"\n")
  sink(paste(catalogue,d,sep="."))
  s <- subset(results,dataset==d)
  print(s[c("ref_rsid","ref_snpid","rsid","r2","p","trait")])
  sink()
}

# Novel pQTL
# rs579459 -- 
# rs4734879 -- 
# rs6993770 -- CXCL5
# rs705379 -- SCF
# rs13229619 -- FGF.21
# rs9272226 -- CDCP1
# rs11759846 -- IL.1.alpha
# rs10733789 -- CCL11
# rs10822155 -- VEGF.A
# rs16840522 -- TRAIL
