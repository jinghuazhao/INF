#!/bin/bash

export TMPDIR=/scratch/jhz22/tmp
export INF=$HOME/INF

#### PhenoScanner

. /etc/profile.d/modules.sh
module load default-cardio
module load use.own

module load R/3.4.2 phenoscanner/phenoscanner_v2
export R_LIBS=

#### snpstats

function rsid()
{
  gunzip -c $INF/work/INTERVAL.rsid.gz | \
  awk '$2!="."' > INTERVAL.rsid
  cut -f4 $INF/aild/cojo/INF1.jma | \
  awk 'NR>1' | \
  uniq | \
  sort -k1,1 | \
  join - INTERVAL.rsid > INF1.rsid
  cut -d' ' -f2 INF1.rsid > INF1.ps
}

cd $INF/ps

phenoscanner -s T -c All -x EUR -p 0.0000001 -r 0.6 -i INF1.ps -o INF1
rsid
for i in $(cut -f1 INF1.jma | awk 'NR>1' | uniq )
do
  echo $i;
  awk -vprot=$i '$1==prot{print $4}' INF1.jma | \
  sort -k1,1 | \
  join INF1.rsid - | \
  awk '{print $2}'> $i.ps;
  phenoscanner -s T -c All -x EUR -p 0.0000001 -r 0.6 -i $i.ps -o $i
done

module load gcc/5.2.0

function MR()
{
# Fine with R 3.6.0 but only possible with 100 SNPs as query from the Web
R --no-save -q <<END
  library(MendelianRandomization)
# rsid prepared from INTERVAL data
  rsid <- with(read.table("INF1.ps",as.is=TRUE,col.names=c("rsid")), rsid)
  r1 <- phenoscanner(snpquery=  rsid[1:100], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r2 <- phenoscanner(snpquery=rsid[101:200], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r3 <- phenoscanner(snpquery=rsid[201:300], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r4 <- phenoscanner(snpquery=rsid[301:400], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r5 <- phenoscanner(snpquery=rsid[401:437], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  save(r1,r2,r3,r4,r5,file="INF1r.rda",version=2)
END
}

function ps_v1.1 ()
# This version works with chr:pos
{
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
}
