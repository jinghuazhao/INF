#!/usr/bin/bash

export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107/
export A2=COVID19_HGI_A2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export B2=COVID19_HGI_B2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz
export C2=COVID19_HGI_C2_ALL_eur_leave_ukbb_23andme_20210107.b37.txt.gz

(
  echo "SNP A1 A2 freq b se p N"
  gunzip -c $HGI/$C2 | \
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
gzip -f > ${INF}/HGI/gsmr_C2.txt.gz

if [ -f ${INF}/HGI/INF1_C2.gsmr ]; then rm ${INF}/HGI/INF1_C2.gsmr; fi
(
  cat ${INF}/HGI/gsmr_C2*.gsmr | \
  head -1
  ls ${INF}/HGI/gsmr_C2*gsmr | \
  parallel -j1 -C' ' '
    if [ -f {} ]; then
       awk "NR>1" {}
    fi
  '
) | \
grep -v nan > ${INF}/HGI/INF1_C2.gsmr

# gunzip -c $HGI/$C2 | head -1 | tr '\t' '\n' | awk '{print "#" NR,$1}'
# export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
# gunzip -c $HGI/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37.txt.gz | \
