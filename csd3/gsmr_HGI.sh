#!/usr/bin/bash

export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107/
export A2=COVID19_HGI_A2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export B2=COVID19_HGI_B2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export C2=COVID19_HGI_C2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz

function trait()
{
  echo ${trait}
  (
    echo "SNP A1 A2 freq b se p N"
    gunzip -c ${HGI}/COVID19_HGI_${trait}_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz | \
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
}

function collect()
{
  if [ -f ${INF}/HGI/INF1_${trait}.gsmr ]; then rm ${INF}/HGI/INF1_${trait}.gsmr; fi
  (
    cat ${INF}/HGI/gsmr_${trait}*.gsmr | \
    head -1
    ls ${INF}/HGI/gsmr_${trait}*gsmr | \
    parallel -j1 -C' ' '
      if [ -f {} ]; then
         awk "NR>1" {}
      fi
    '
  ) | \
  grep -v nan > ${INF}/HGI/INF1_${trait}.gsmr
}

for trait in A2 B2 C2
do
  trait
  collect
done

(
  echo prot uniprot A2 b_A2 se_A2 p_A2 n_A2 B2 b_B2 se_B2 p_B2 n_B2 C2 b_C2 se_C2 p_C2 n_C2
  join -a1 -e "NA"   -o1.1,1.2,2.2,2.3,2.4,2.5,2.6 <(sort -k1,1 ${INF}/work/inf1.tmp) <(sed '1d' ${INF}/HGI/INF1_A2.gsmr | cut -f1-6 | sort -k1,1) | \
  join -a1 -e "NA" - -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6 <(sed '1d' ${INF}/HGI/INF1_B2.gsmr | cut -f1-6 | sort -k1,1 ) | \
  join -a1 -e "NA" - -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3,2.4,2.5,2.6 \
                     <(sed '1d' ${INF}/HGI/INF1_C2.gsmr | cut -f1-6 | sort -k1,1) | \
  awk '$3!="NA"'
) > ${INF}/HGI/A2-B2-C2.txt

join -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8 \
     <(gunzip -c ${INF}/HGI/gsmr_B2_LIF.R.eff_plot.gz | \
       awk '/effect_begin/,/effect_end/' | \
       grep -v effect) \
     <(gunzip -c ${INF}/HGI/gsmr_C2_LIF.R.eff_plot.gz | \
       awk '/effect_begin/,/effect_end/' | \
       grep -v effect) | \
join -a1 -e "NA" -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16 \
       - \
       <(gunzip -c ${INF}/HGI/gsmr_A2_LIF.R.eff_plot.gz | \
       awk '/effect_begin/,/effect_end/' | \
       grep -v effect) | \
awk '{$1=$9};1' | \
join work/INTERVAL.rsid - > ${INF}/HGI/A2-B2-C2.snp_effects

R --no-save -q <<END
  INF <- Sys.getenv("INF")
  gsmr <- within(read.table(file.path(INF,"HGI","A2-B2-C2.txt"),header=TRUE),{col="black"})
  recolor <- with(gsmr,prot=="LIF.R")
  gsmr[recolor,"col"] <- "red"
  n.prot <- nrow(gsmr)
  xtick <- seq(1,n.prot)
  png(file.path(INF,"HGI","A2-B2-C2.png"), res=300, units="cm", width=40, height=20)
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
  source("http://cnsgenomics.com/software/gcta/res/gsmr_plot.r")
  read_gsmr_by_trait_p <- function(trait,p)
  {
      gsmr_data <- read_gsmr_data(file.path(INF,"HGI",paste0("gsmr_",trait,"_",p,".eff_plot.gz")))
      plot_gsmr_effect(gsmr_data, p, trait, colors()[75])
      snp_effect <- with(gsmr_data,snp_effect[,-c(2:6)])
      colnames(snp_effect) <- c("snpid",paste0(trait,"_",c("b1","se1","b2","se2")))
      as.data.frame(snp_effect)
  }
  A2 <- read_gsmr_by_trait_p("A2","LIF.R")
  B2 <- read_gsmr_by_trait_p("B2","LIF.R")
  C2 <- read_gsmr_by_trait_p("C2","LIF.R")
  A2_B2_C2 <- merge(merge(A2,B2,by="snpid",all.y=TRUE),C2,by="snpid")
  snp_effect_id <- read.table(file.path(INF,"HGI","A2-B2-C2.snp_effects"))[,1:2]
  snp_effects <- data.frame(snpid=snp_effect_id[["V2"]],apply(A2_B2_C2[,-1],2,as.numeric))
  par(mfrow=c(1,3))
  gap::ESplot(snp_effects[c("snpid","A2_b2","A2_se2")], logscale=FALSE)
  gap::ESplot(snp_effects[c("snpid","B2_b2","B2_se2")], logscale=FALSE)
  gap::ESplot(snp_effects[c("snpid","C2_b2","C2_se2")], logscale=FALSE)
END

# gunzip -c $HGI/$C2 | head -1 | tr '\t' '\n' | awk '{print "#" NR,$1}'
# export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
# gunzip -c $HGI/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37.txt.gz | \
