#!/usr/bin/bash

join <(sort -k1,1 h2.tsv) <(sed '1d' pve.dat | sort -k1,1) > h2pve.dat
R --no-save -q <<END
  h2pve <- read.table("h2pve.dat",col.names=c("prot","h2","h2se","pve","sepve","m"))
  with(h2pve,summary(h2))
  subset(h2pve, h2>0.3 & pve > 0.3)
  png("h2pve.png",width=6,height=5,units="cm",pointsize=4,res=300)
  attach(h2pve)
  cor(h2,pve)
  plot(h2,pve,pch=19,cex.axis=2,cex.lab=2,xlab="h2",ylab="pve")
  abline(lm(pve~h2), col="red")
  lines(lowess(h2,pve), col="blue")
  detach(h2pve)
  dev.off()
END

R --no-save -q <<END
  h2 <- read.table("h2.tsv",as.is=TRUE,col.names=c("prot","h2","se"))
  ord <- with(h2, order(h2))
  sink("h2.dat")
  print(h2[ord, c("prot","h2","se")], row.names=FALSE)
  sink()
  png("h2-pve.png", res=300, units="cm", width=40, height=40)
  par(mfrow=c(2,1))
  np <- nrow(h2)
  with(h2[ord,], {
    plot(h2, cex=0.8, pch=16, axes=FALSE, main="a")
    xy <- xy.coords(h2)
    l <- h2-1.96*se
    l[l<0] <- 0
    u <- h2+1.96*se
    segments(xy$x,l, xy$x,u)
    xtick <- seq(1, np, by=1)
    axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
    axis(side=2)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.85)
  })
  pve <- read.table("pve.dat",as.is=TRUE,header=TRUE)
  np <- nrow(pve)
  with(pve, {
      plot(pve, cex=0.8, pch=16, axes=FALSE, main="b")
      xy <- xy.coords(pve)
      se <- sqrt(v)
      segments(xy$x, pve-1.96*se, xy$x, pve+1.96*se)
      xtick <- seq(1, np, by=1)
      axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
      axis(side=2)
      text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.85)
  })
  dev.off()
END
