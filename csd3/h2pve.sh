#!/usr/bin/bash

join <(sort -k1,1 h2.tsv) <(sed '1d' pve.dat | sort -k1,1) > h2pve.dat
R --no-save -q <<END
  png("h2pve.png", res=300, units="cm", width=40, height=20)
  par(mfrow=c(1,2))
  with(read.delim("INF1.METAL",as.is=TRUE),
  {
    MAF <- Freq1
    repl <- MAF > 1-MAF
    MAF[repl] <- 1-MAF[repl]
    plot(MAF,abs(Effect),cex.axis=2,cex.lab=2,pch=19,main="a",xlab="MAF",
         ylab="Effect size",col=c("red","blue")[2-(cis.trans=="cis")])
  })
  h2pve <- read.table("h2pve.dat",col.names=c("prot","h2","h2se","pve","sepve","m"))
  with(h2pve,summary(h2))
  subset(h2pve, h2>0.3 & pve > 0.3)
  attach(h2pve)
  cor(h2,pve)
  plot(h2,pve,pch=19,cex.axis=2,cex.lab=2,main="b",xlab="h2",ylab="pve")
  abline(lm(pve~h2), col="red")
  lines(lowess(h2,pve), col="blue")
  detach(h2pve)
  dev.off()
END

R --no-save -q <<END
  png("h2-pve.png", res=300, units="cm", width=40, height=40)
  par(mfrow=c(3,1))
  h2 <- read.table("h2.tsv",as.is=TRUE,col.names=c("prot","h2","se"))
  summary(h2)
  ord <- with(h2, order(h2))
  sink("h2.dat")
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
  ldak <- read.table("ldak/INF1.ldak.h2",as.is=TRUE,skip=1)
  names(ldak) <- c("prot","h2","se","inf","inf_se")
  summary(ldak)
  ord <- with(ldak,order(h2))
  sink("ldak.dat")
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
  pve <- read.table("pve.dat",as.is=TRUE,header=TRUE)
  summary(pve)
  np <- nrow(pve)
  with(pve, {
      plot(pve, cex=0.8, pch=16, axes=FALSE, main="c", xlab="Protein", cex.lab=2)
      xy <- xy.coords(pve)
      se <- sqrt(v)
      segments(xy$x, pve-1.96*se, xy$x, pve+1.96*se)
      xtick <- seq(1, np, by=1)
      axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
      axis(side=2, cex.axis=2)
      text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=1.2)
  })
  dev.off()
END
