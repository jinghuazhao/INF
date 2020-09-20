#!/usr/bin/bash

(
  echo snpid chr pos allele1 b1 se1 log10p1 a1 b2 se2 p2
  cut -f6 work/INF1.merge | \
  sed 's/chr//;s/_/:/g' | \
  zgrep -f - ukb/30710_irnt.gwas.imputed_v3.both_sexes.tsv.bgz | \
  cut -f1,2,8,9,11 | \
  awk '
  {
    split($1,a,":")
    if (a[3] < a[4]) {a1=a[3];a2=a[4]} else {a1=a[4];a2=a[3]}
    snpid=a[1] ":" a[2] ":" a1 ":" a2
    print snpid,toupper($2),$3,$4,$5
  }' | \
  sort -k1,1 | \
  join <(cut -f1,4-6,9-11 work/INF1.METAL | sed 's/chr//;s/_/:/g' | sort -k1,1) - | \
  sort -k2,2n -k3,3n | \
  awk '{
    $4=toupper($4)
    if($4!=$8) $9=-$9
    print
  }'
) > work/INF1.merge.crp

R --no-save <<END
  options(width=200)
  png("work/INF1.crp.png", res=300, units="cm", width=40, height=20)
  crp <- read.table("work/INF1.merge.crp",as.is=TRUE,header=TRUE)
  crp
  with(crp,plot(b1,b2,pch=19,xlab="INF effect size",ylab="CRP effect size"))
  dev.off()
END
