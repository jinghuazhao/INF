# 13-4-2020 JHZ

require(phenoscanner)
catalogue <- "pQTL"
rsid <- scan("work/INF1.snp",what="")
batches <- split(rsid, ceiling(seq_along(rsid)/100))
s <- t <- list()
for(i in 1:length(batches))
{
  cat("Block ",i,"\n")
  q <- phenoscanner(snpquery=batches[[i]], catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  s[[i]] <- with(q,snps)
  t[[i]] <- with(q,results)
}
snps <- do.call(rbind,s)
results <- do.call(rbind,t)
results <- within(results,{
   a1 <- ref_a1
   a2 <- ref_a2
   swap <- ref_a1 > ref_a2
   a1[swap] <- ref_a2[swap]
   a2[swap] <- ref_a1[swap]
   ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
})
r <- list(snps=snps,results=results)
save(r,file=paste0("work/INF1.merge.",catalogue))

options(width=500)
attach(r)
for(d in unique(with(results,dataset)))
{
  cat(d,"\n")
  vars <- c("ref_rsid","ref_snpid","rsid","r2","p","trait","dataset","pmid")
  s <- subset(results[vars],dataset==d)
  if (d=="Sun-B_pQTL_EUR_2017")
  {
     gs <-read.delim("INTERVAL_box.tsv",as.is=TRUE)
     m <- merge(s,gs,by.x="trait",by.y="TargetFullName")
     s <- m[c("ref_rsid","ref_snpid","rsid","r2","p","trait","UniProt","UniProts","symbol")]
  }
  sink(paste(catalogue,d,sep="."))
  print(s)
  sink()
}
detach(r)
