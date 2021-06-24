#!/bin/bash
# 31-5-2019 JHZ

export rt=/scratch/jhz22/INF/snps/clump0
cut -d' ' -f3 $rt/INF1.clumped | \
awk 'NR>1' | \
sort -k1,1 | \
uniq | \
join work/INTERVAL.rsid - > INF1.rsid

R --no-save <<END
  prot <- MarkerName <- NA
  tbl <- read.table("/scratch/jhz22/INF/snps/clump0/INF1.clumped",as.is=TRUE,header=TRUE)
  tbl <- tbl[setdiff(names(tbl), c("NSIG", "S05", "S01", "S001", "S0001", "SP2"))]
  rsid <- read.table("INF1.rsid",as.is=TRUE,col.names=c("SNP","rsid"))
  requireNamespace("dplyr")
  d <- dplyr::nest_join(tbl,rsid)
  dy <- d["y"]
  rsid <- unlist(dy)
  d <- d[setdiff(names(d),"y")]
  t <- cbind(d,rsid)
  ord <- with(t,order(prot,CHR))
  tbl <- t[ord,]
  write.table(tbl,file="INF1.clumped",col.names=FALSE,row.names=FALSE,quote=FALSE)
END

export rt=/scratch/jhz22/INF/snps/cojo
(
  sort -k1,1 -k2,2n ${rt}/INF1.lz | \
  parallel --env rt -j1 -C' ' '
    awk "NR>1 && \$2 <=5e-10" $rt/lz/{1}-{2}.lz | \
    sort -k2,2g | \
    awk -vOFS="\t" -vprot={1} -vregion={2} "NR==1{print "prot", "region", \$0}" 
  '
) > sentinels.txt
