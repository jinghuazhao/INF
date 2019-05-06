#!/bin/bash

. /etc/profile.d/modules.sh
module load default-cardio
module load use.own

# PhenoScanner uses old R
module load R/3.4.2 phenoscanner/phenoscanner_v1.1
export R_LIBS=

cd ps
for i in $(cut -f1 INF1.jma | awk 'NR>1' | uniq )
do
  echo $i
  awk -vprot=$i '$1==prot{print "chr" $2 ":" $4}' INF1.jma > $i.ps
# default -r 0.6
  phenoscanner -c All -l 1000G -p 0.0000001 -i $i.ps -o $i
done
cd -

function MR()
{
# Fine with R 3.6.0 but only possible with 100 SNPs as query from the Web
R --no-save -q <<END
  library(MendelianRandomization)
# SNPID not possible
  INF1 <- read.table("INF1.ps",as.is=TRUE)
# batches in 100 only
  INF1 <- with(read.csv("INF1_PhenoScanner_SNP_Info.csv"),rsID)
  r1 <- phenoscanner(snpquery=  INF1[1:100], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r2 <- phenoscanner(snpquery=INF1[101:200], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r3 <- phenoscanner(snpquery=INF1[201:300], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r4 <- phenoscanner(snpquery=INF1[301:376], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  save(r1,r2,r3,r4,file="INF1r.rda",version=2)
END
}
