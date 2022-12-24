#!/bin/bash

qctool -g ${INF}/h2/#.bed -s ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_1.sample -filetype binary_ped -ofiletype binary_ped -og ${INF}/h2/INTERVAL
plink --bfile ${INF}/h2/INTERVAL --make-grm-bin --threads 10 --out ${INF}/h2/INTERVAL
awk -vOFS="\t" '{print $2,$2}' ${INF}/h2/INTERVAL.grm.id > ${INF}/h2/INTERVAL.grm.id2
cp ${INF}/h2/INTERVAL.grm.id2 ${INF}/h2/INTERVAL.grm.id

export s=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
cut -d' ' -f1-2,5-7 ${s} | awk 'NR>3' > ${INF}/h2/INTERVAL.covar
cut -d' ' -f1-2,4,8-28 ${s} | awk 'NR>3' > ${INF}/h2/INTERVAL.qcovar
cut -d' ' -f1-2,29- ${s} | awk 'NR>3' > ${INF}/h2/INTERVAL.pheno

cd ${INF}/h2
grep V\(G\) *hsq | grep Vp | sed 's|.hsq:V(G)/Vp||g' > h2.tsv
R --no-save -q <<END
  h2 <- read.table("h2.tsv",as.is=TRUE,col.names=c("prot","h2","se"))
  ord <- with(h2, order(h2))
  sink("h2.dat")
  print(h2[ord, c("prot","h2","se")], row.names=FALSE)
  sink()
  png("h2.png", res=300, units="in", width=6, height=4)
  np <- nrow(h2)
  with(h2[ord,], {
    plot(h2, cex=0.4, pch=16, xaxt="n", xlab="protein", ylab=expression(h^2))
    xtick <- seq(1, np, by=1)
    axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.01)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.3)
  })
  title("(a)")
  dev.off()
END
cd -
