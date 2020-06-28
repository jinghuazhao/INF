#!/usr/bin/bash

# count genomic regions for all sentinels
(
  awk -v OFS='\t' '(NR==1){print $1,$2,$3,$5,$6,$8,$9}' work/INF1.merge 
  awk -v OFS='\t' '(NR>1){print $1,$2,$3,$5,$6,$8,$9}' work/INF1.merge | \
  sort -k6,6n -k2,2n
) | \
bedtools merge | \
wc -l

# ST2
R --no-save -q <<END
  cvt<- read.table("work/INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  m <- read.delim("work/INF1.merge",as.is=TRUE)
  m <- merge(m[c("prot","MarkerName","Start","End")],cvt[c("uniprot","prot","SNP","cis.trans")],by.x=c("prot","MarkerName"),by.y=c("prot","SNP"))
  load("ds/latest/fp/INF1.rda")
  e <- with(tbl,prot=="CCL25" & MarkerName=="chr19:49206145_C_G")
  tbl <- subset(tbl,!e)
  m <- merge(tbl,m,by=c("prot","MarkerName"))
  m <- merge(rsid,m,by="MarkerName")
  s <- setdiff(names(m),c("FreqSE","MinFreq","MaxFreq"))
  write.table(m[s],file="work/INF1.METAL",quote=FALSE,row.names=FALSE,sep="\t")
END

R --no-save -q <<END
  tbl <- read.delim("work/INF1.tbl",as.is=TRUE)
  attach(tbl)
  MAF <- Freq1
  repl <- MAF > 1-MAF
  MAF[repl] <- 1-MAF[repl]
  png("work/b-maf.png",width=10,height=8,units="in",pointsize=4,res=300)
  plot(MAF,abs(Effect),cex=0.85,xlab="MAF",ylab="Effect size")
  dev.off()
END

R --no-save -q <<END
  cvt <- read.table("work/INF1.merge.out",as.is=TRUE,header=TRUE,nrows=70)
  H <- with(cvt,table(total))
  M <- names(H)
  png(file = "work/signals_by_protein.png",width=6,height=5,units="in",pointsize=4,res=300)
  barplot(H,names.arg=M,xlab="Protein",ylab="Signals",col="blue",
  main="Number of signals by protein",border="red")
  dev.off()
END
