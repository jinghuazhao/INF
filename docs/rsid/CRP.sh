#!/usr/bin/bash

function CRP()
{
export f=work/INF1.crp${tag}
(
  echo snpid chr pos allele1 b1 se1 log10p1 a1 b2 se2 p2
  cat work/INF1.merge.cis | \
  sed 's/chr//;s/_/:/g' | \
  zgrep -f - ukb/30710${tag}.gwas.imputed_v3.both_sexes.tsv.bgz | \
  cut -f1,2,8,9,11 | \
  awk '
  {
    split($1,a,":")
    if (a[3] < a[4]) {a1=a[3];a2=a[4]} else {a1=a[4];a2=a[3]}
    snpid=a[1] ":" a[2] ":" a1 ":" a2
    print snpid,toupper(a[4]),$3,$4,$5
  }' | \
  sort -k1,1 | \
  join <(awk '$NF=="cis"' work/INF1.METAL | cut -f1,4-6,9-11 | sed 's/chr//;s/_/:/g' | sort -k1,1) - | \
  sort -k2,2n -k3,3n | \
  awk '{$4=toupper($4)};1'
) > ${f}

R --no-save -q <<END
  options(width=200)
  tag <- Sys.getenv("tag")
  f <- Sys.getenv("f")
  png(paste0(f,".png"), res=300, units="cm", width=40, height=40)
  crp <- read.table(f,as.is=TRUE,header=TRUE)
  swap <- with(crp,allele1!=a1)
  crp
  with(crp,cor.test(b1,b2))
  crp[swap,"b2"] <- -crp[swap,"b2"]
  with(crp,cor.test(b1,b2))
  with(crp,plot(b1,b2,pch=19,xlab="INF effect size",ylab=paste("CRP effect size (",tag,")")))
  dev.off()
END
}

export tag=_raw
CRP
export tag=_irnt
CRP
