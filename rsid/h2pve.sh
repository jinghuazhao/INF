#!/usr/bin/bash

Rscript -e '
  INF <- Sys.getenv("INF")
  require(gap)
  t <- read.delim(file.path(INF,"work","INF1.tbl"),as.is=TRUE)
  tbl <- within(t, {
      prot <- sapply(strsplit(Chromosome,":"),"[",1)
      Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
  })
## to obtain variance explained
  tbl <- within(tbl,
  {
    x2 <- (Effect/StdErr)^2
    r2 <- x2 / (N - 2 + x2)
    v <- 2 / (N - 1)
  })
  s <- with(tbl, aggregate(r2,list(prot),sum))
  names(s) <- c("prot", "pve")
  se2 <- with(tbl, aggregate(v,list(prot),sum))
  names(se2) <- c("p1","se")
  m <- with(tbl, aggregate(r2,list(prot),length))
  names(m) <- c("p2","m")
  pve <- within(cbind(s,se2,m),{se=sqrt(se)})
  ord <- with(pve, order(pve))
  sink(file.path(INF,"h2","pve.dat"))
  print(pve[ord, c("prot","pve","se","m")], row.names=FALSE)
  sink()
  write.csv(tbl,file=file.path(INF,"h2","INF1.csv"),quote=FALSE,row.names=FALSE)
'

join <(sort -k1,1 ${INF}/h2/h2.tsv) <(sed '1d' ${INF}/h2/pve.dat | sort -k1,1) > ${INF}/h2/h2pve.dat

Rscript -e '
  INF <- Sys.getenv("INF")
  png(file.path(INF,"h2","h2pve.png"), res=300, units="cm", width=40, height=20)
  par(mfrow=c(1,2))
  INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE)
  t.test(Effect~as.factor(cis.trans),data=INF1_METAL)
  summary(lm(Effect~as.factor(cis.trans),data=INF1_METAL))
# boxplot(Effect~as.factor(cis.trans),data=INF1_METAL)
  with(INF1_METAL,
  {
    MAF <- Freq1
    repl <- MAF > 1-MAF
    MAF[repl] <- 1-MAF[repl]
    v <- 2*MAF*(1-MAF)*Effect^2
    col <- c("blue","red")[1+(cis.trans=="cis")]
#   plot(MAF,abs(Effect),cex.axis=1.5,cex.lab=1.5,pch=19,main="a",xlab="MAF",ylab="Effect size",col=col)
    library(scatterplot3d)
    scatterplot3d(MAF,abs(Effect),v,color=c("blue","red")[(cis.trans=="cis")+1], 
                  main="a",cex.axis=1.5,cex.lab=1.5,
                  box=TRUE, pch=16, type="h", xlab="MAF", ylab="Effect", zlab=expression(italic(2*MAF(1-MAF)*Effect^2)))
    legend("right", legend=levels(as.factor(cis.trans)), box.lwd=0, col=c("red", "blue"), pch=16)
  })
  h2pve <- read.table(file.path(INF,"h2","h2pve.dat"),col.names=c("prot","h2","h2se","pve","sepve","m"))
  with(h2pve,summary(h2))
  with(h2pve,cor(pve,h2))
  subset(h2pve, h2>0.3 & pve > 0.3)
  attach(h2pve)
  cor(h2,pve)
  plot(h2,pve,pch=19,cex.axis=1.5,cex.lab=1.5,main="b",xlab=expression(italic(h^2)),ylab="")
  mtext("PVE", side=2, at=0.13, line=3, cex=1.5)
  reg <- lm(pve~h2)
  summary(reg)
  abline(reg, col="red")
# lines(lowess(h2,pve), col="blue")
  detach(h2pve)
  dev.off()
'

function merge()
{
Rscript -e '
# This part needs to be run inside R
  options(width=200)
  INF <- Sys.getenv("INF")
  png(file.path(INF,"h2","h2-pve.png"), res=300, units="cm", width=40, height=40)
  par(mfrow=c(3,1))
  h2 <- read.table(file.path(INF,"h2","h2.tsv"),as.is=TRUE,col.names=c("prot","h2","se"))
  summary(h2)
  ord <- with(h2, order(h2))
  sink(file.path(INF,"h2","h2.dat"))
  print(h2[ord, c("prot","h2","se")], row.names=FALSE)
  sink()
  np <- nrow(h2)
  with(h2[ord,], {
    plot(h2, cex=0.8, pch=16, axes=FALSE, main="a", xlab="", cex.lab=2)
    xy <- xy.coords(h2)
    l <- h2-1.96*se
    l[l<0] <- 0
    u <- h2+1.96*se
    segments(xy$x,l, xy$x,u)
    xtick <- seq(1, np, by=1)
    axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
    axis(side=2, cex.axis=2)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=1.2)
  })
  ldak <- read.table(file.path(INF,"ldak","INF1.ldak.h2"),as.is=TRUE,skip=1)
  names(ldak) <- c("prot","h2","se","inf","inf_se")
  summary(ldak)
  ord <- with(ldak,order(h2))
  sink(file.path(INF,"ldak","ldak.dat"))
  print(ldak[ord, c("prot","h2","se")], row.names=FALSE)
  sink()
  with(ldak[ord,],{
    plot(h2, cex=0.8, pch=16, axes=FALSE, main="b", xlab="", cex.lab=2)
    xy <- xy.coords(h2)
    l <- h2-1.96*se
    l[l<0] <- 0
    u <- h2+1.96*se
    segments(xy$x,l, xy$x,u)
    xtick <- seq(1, np, by=1)
    axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
    axis(side=2, cex.axis=2)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=1.2)
  })
  pve <- read.table(file.path(INF,"h2","pve.dat"),as.is=TRUE,header=TRUE)
  summary(pve)
  np <- nrow(pve)
  with(pve, {
      plot(pve, cex=0.8, pch=16, axes=FALSE, main="c", xlab="Protein", cex.lab=2)
      xy <- xy.coords(pve)
      segments(xy$x, pve-1.96*se, xy$x, pve+1.96*se)
      xtick <- seq(1, np, by=1)
      axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
      axis(side=2, cex.axis=2)
      text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=1.2)
  })
  dev.off()
  names(h2) <- c("prot","h2_interval","SE_h2_interval")
  names(ldak) <- c("prot","h2_scallop","SE_h2_scallop", "inf", "se_inf")
  names(pve) <- c("prot","pve","SE_pve","m")
  h2_ldak_pve <- merge(merge(h2,ldak,by="prot"),pve,by="prot",all=TRUE)
  write.csv(h2_ldak_pve,file=file.path(INF,"ldak","h2-ldak-pve.csv"),quote=FALSE,row.names=FALSE)
'
}
