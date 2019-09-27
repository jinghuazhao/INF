#!/bin/bash

export TMPDIR=/rds/user/jhz22/hpc-work/work
export rt=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/INTERVAL
export s=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/o5000-inf1-outlier_in-r2.sample

cut -d' ' -f1-2,5-7 ${s} | awk 'NR>2' > ${rt}.covar
cut -d' ' -f1-2,4,8-28 ${s} | awk 'NR>2' > ${rt}.qcovar
cut -d' ' -f1-2,29- ${s} | awk 'NR>2' > ${rt}.pheno

plink --bfile ${rt} \
      --exclude range tryggve/high-LD-regions-hg19.txt \
      --indep-pairwise 500kb 1 0.80 --maf 0.01 --out $rt
plink --bfile ${rt} --extract ${rt}.prune.in --make-bed --out ${rt}.prune
plink --bfile ${rt}.prune --make-grm-bin --threads 2 --out ${rt}

cut -d' ' -f29- ${s} | head -1 | sed 's/ /\n/g' | awk '{split($1,a,"__"); print a[1]}' > h2.list

sbatch --wait csd3/h2.sb

cd work
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
    xy <- xy.coords(h2)
    segments(xy$x,h2-1.96*se, xy$x,h2+1.96*se)
    xtick <- seq(1, np, by=1)
    axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.01)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.3)
  })
  dev.off()
END
cd -
