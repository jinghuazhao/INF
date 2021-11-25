#!/usr/bin/bash

require(phenoscanner)
catalogue <- Sys.getenv("catalogue")
pvalue <- Sys.getenv("pvalue")
r2 <- Sys.getenv("r2")
rsid <- Sys.getenv("rsid")
rsid <- scan(rsid,what="")
batches <- split(rsid, ceiling(seq_along(rsid)/100))
s <- t <- list()
for(i in 1:length(batches))
{
  cat("Block ",i,"\n")
  q <- phenoscanner(snpquery=batches[[i]], catalogue=catalogue, proxies="EUR", pvalue=as.numeric(pvalue), r2=as.numeric(r2), build=37)
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
save(r,file=paste0("ps/INF1.merge.",catalogue))

snps_results <- with(r,dplyr::right_join(snps,results))
for(d in unique(with(snps_results,dataset)))
{
  cat(d,"\n")
  vars <- c("ref_rsid","ref_snpid","rsid","ref_chr","ref_pos_hg19","chr","pos_hg19","proxy","r2","p","trait","dataset","pmid","study")
  s <- subset(snps_results[vars],dataset==d)
  HOME <- Sys.getenv("HOME")
  if (d=="Sun-B_pQTL_EUR_2017")
  {
     gs <-read.delim(file.path(HOME,"SomaLogic","doc","SOMALOGIC_Master_Table_160410_1129info.tsv"),as.is=TRUE)
     m <- merge(s,gs,by.x="trait",by.y="TargetFullName")
     s <- m[c(vars,"UniProt")]
  }
  write.table(s,file.path(HOME,"INF/ps",paste(catalogue,d,"tsv",sep=".")),quote=FALSE,row.names=FALSE,sep="\t")
}
