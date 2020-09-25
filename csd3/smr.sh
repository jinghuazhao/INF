#!/usr/bin/bash

# setup
smr --beqtl-summary ~/COVID-19/smr/cage_eqtl_data_lite_hg19/CAGE.sparse.lite --extract-snp-p 5e-8 --add-n 2765 --make-besd --out CAGE
smr --beqtl-summary CAGE --show-n 
ln -sf CAGE.besd CAGE_snpid.besd
ln -sf CAGE.epi CAGE_snpid.epi
awk '
{
  chr=$1
  pos=$4
  a1=$5
  a2=$6
  if (a1<a2) snpid="chr" chr ":" pos "_" a1 "_" a2;
  else snpid="chr" chr ":" pos "_" a2 "_" a1;
  $2=snpid
  print
}' OFS='\t' CAGE.esi > CAGE_snpid.esi

# SNP query
sed '1d' ${INF}/work/INF1.merge | \
cut -f6 | \
sort | \
uniq > INF1.merge.MarkerName
smr --beqtl-summary CAGE_snpid --extract-snp INF1.merge.MarkerName --query 5.0e-8 --out INF1.merge
cut -f5,6 ${INF}/work/INF1.merge | \
awk 'NR>1{print $2,$1}' | \
sort -k1,1 | \
join - <(sed '1d' INF1.merge.txt | sort -k1,1) | \
awk '$2==$11' > INF1.merge.coloc
awk '{print $2,$1}' INF1.merge.coloc | \
sort -k1,2 | \
uniq | \
parallel -j1 -C' ' 'awk -vm={2} -vp={1} "\$5==p && \$6==m" ${INF}/work/INF1.merge'

# formal SMR
echo 6 25392021 33392022 1 > HLA.hg19
sed '1d' ${INF}/work/INF1.merge | \
tr '\t' ' ' | \
parallel -j1 -C' ' '
# PLINK file
  gunzip -c ${INF}/work/INF1.merge.*-{5}-{6}.gz | \
  cut -f3 | \
  sed "1d" > {5}-{6}.snps
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --chr {8} --exclude range HLA.hg19 --extract {5}-{6}.snps --make-bed --out {5}-{6}
# BLD format
  smr --bfile {5}-{6} --make-bld --r --ld-wind 4000 --out {5}-{6}
# SMR and HEIDI test
  smr --bld {5}-{6} --gwas-summary ${INF}/work/{5}.ma --beqtl-summary CAGE_snpid --out {5} --thread-num 10 
'
