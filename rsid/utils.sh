#!/usr/bin/bash

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
  load("ds/latest/fp/INF1.rda")
  e <- with(tbl,prot=="CCL25" & MarkerName=="chr19:49206145_C_G")
  tbl <- subset(tbl,!e)
  m <- merge(rsid,tbl,by="MarkerName")
  s <- setdiff(names(m),c("FreqSE","MinFreq","MaxFreq"))
  write.table(m[s],file="work/INF1.METAL",quote=FALSE,row.names=FALSE,sep="\t")
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
