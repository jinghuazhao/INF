#!/usr/bin/bash

export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107/
export A2=COVID19_HGI_A2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export B2=COVID19_HGI_B2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export C2=COVID19_HGI_C2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export A2=COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.txt.gz
export B2=COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.txt.gz
export C2=COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.txt.gz

for trait in A2 B2 C2
do
  echo ${trait}
  (
    echo "SNP A1 A2 freq b se p N"
    gunzip -c ${HGI}/COVID19_HGI_${trait}_ALL_eur_leave_23andme_20210107.b37.txt.gz | \
    awk '
    {
      CHR=$1
      POS=$2
      a1=$4
      a2=$3
      if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
      else snpid="chr" CHR ":" POS "_" a1 "_" a2
      if (NR>1) print snpid, a1, a2, $12, $7, $8, $9, $11
    }'
  ) | \
  awk 'a[$1]++==0' | \
  gzip -f > ${INF}/HGI/gsmr_${trait}.txt.gz
done

R --no-save -q <<END
  INF <- Sys.getenv("INF")
  gsmr <- within(read.table(file.path(INF,"HGI","A2-B2-C2.txt"),header=TRUE),{col="black"})
  recolor <- with(gsmr,prot=="LIF.R")
  gsmr[recolor,"col"] <- "red"
  n.prot <- nrow(gsmr)
  xtick <- seq(1,n.prot)
  png(file.path(INF,"HGI","A2-B2-C2.gsmr.png"), res=300, units="cm", width=40, height=20)
  with(gsmr, {
      par(mfrow=c(3,1))
      plot(-log10(p_A2), col=col, cex=2, pch=16, axes=FALSE, main="Critical illness (A2)", xlab="", ylab="", cex.lab=1)
      axis(side=2, cex.axis=2)
      text(0.5,3.5,"-log10(p)", cex=1.5)
      plot(-log10(p_B2), col=col, cex=2, pch=17, axes=FALSE, main="Hospitalisation (B2)", xlab="", ylab="",cex.lab=1)
      axis(side=2, cex.axis=2)
      text(0.5,1.5,"-log10(p)", cex=1.5)
      plot(-log10(p_C2), col=col, cex=2, pch=18, axes=FALSE, main="Reported infection (C2)", xlab="", ylab="",cex.lab=1)
      axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
      axis(side=2, cex.axis=2)
      text(0.5,6,"-log10(p)", cex=1.5)
      text(x=xtick, y=-1, col=col, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=1.2)
  })
  dev.off()
# source("http://cnsgenomics.com/software/gcta/res/gsmr_plot.r")
  source(file.path(INF,"csd3","gsmr_plot.r"))
  read_gsmr_by_trait_p <- function(trait,p)
  {
      gsmr_data <- read_gsmr_data(file.path(INF,"HGI",paste0("gsmr_",trait,"_",p,".eff_plot.gz")))
      plot_gsmr_effect(gsmr_data, p, trait, colors()[75])
      snp_effect <- with(gsmr_data,snp_effect[,-c(2:6)])
      colnames(snp_effect) <- c("snpid",paste0(trait,"_",c("b1","se1","b2","se2")))
      as.data.frame(snp_effect)
  }
  png(file.path(INF,"HGI","A2-B2-C2.mr.png"),res=300, units="cm", width=40, height=20)
  par(mfrow=c(1,3))
  A2 <- read_gsmr_by_trait_p("A2","LIF.R")
  B2 <- read_gsmr_by_trait_p("B2","LIF.R")
  C2 <- read_gsmr_by_trait_p("C2","LIF.R")
  dev.off()
  A2_B2_C2 <- merge(merge(A2,B2,by="snpid",all.y=TRUE),C2,by="snpid")
  snp_effect_id <- read.table(file.path(INF,"HGI","A2-B2-C2.snp_effects"))[,1:2]
  snp_effects <- data.frame(snpid=snp_effect_id[["V2"]],apply(A2_B2_C2[,-1],2,as.numeric))
  png(file.path(INF,"HGI","A2-B2-C2.ESplot.png"),res=300, units="cm", width=40, height=20)
  par(mfrow=c(1,3))
  gap::ESplot(snp_effects[c("snpid","A2_b2","A2_se2")], lty=2, v=0, xlim=c(-0.4,0.4), logscale=FALSE)
  snp_effects["snpid"] <- ""
  gap::ESplot(snp_effects[c("snpid","B2_b2","B2_se2")], lty=2, v=0, xlim=c(-0.3,0.4), logscale=FALSE)
  gap::ESplot(snp_effects[c("snpid","C2_b2","C2_se2")], lty=2, v=0, xlim=c(-0.3,0.4), logscale=FALSE)
  dev.off()
END

# gunzip -c $HGI/$C2 | head -1 | tr '\t' '\n' | awk '{print "#" NR,$1}'
# export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
# gunzip -c $HGI/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37.txt.gz | \
