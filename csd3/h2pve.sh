#!/usr/bin/bash

join <(sort -k1,1 h2.tsv) <(sed '1d' pve.dat | sort -k1,1) > h2pve.dat
R --no-save -q <<END
  h2pve <- read.table("h2pve.dat",col.names=c("prot","h2","h2se","pve","sepve","m"))
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
