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
  smr --bld {5}-{6} --gwas-summary ${INF}/work/{5}.ma --beqtl-summary CAGE_snpid --out {5}-{6} --thread-num 10 
'

(
  cat *.smr | \
  head -1 | \
  awk -vOFS="\t" "{print \"prot\",\"MarkerName\",\"gene\",\"cistrans\",\$0}"
  sed '1d' ${INF}/work/INF1.merge.cis.vs.trans | \
  cut -d, -f2,5,10,14 | \
  tr ',' ' ' | \
  parallel -j1 -C' ' '
    if [ -f {1}-{2}.smr ]; then
       sed "1d" {1}-{2}.smr | \
       awk -vprot={1} -vsnpid={2} -vgene={3} -vcistrans={4} -vOFS="\t" "{print prot,snpid,gene,cistrans,\$0}"
    fi
  '
) > INF1.merge.smr

(
  awk 'NR==1 {$1="topSNPid toprsid MarkerName rsid";$2="prot";$9="";print}' INF1.merge.smr | \
  awk '{$1=$1};1'
  join -22 ${INF}/work/INTERVAL.rsid <(sed '1d' INF1.merge.smr | sort -k2) | \
  sort -k10 | \
  join -210 ${INF}/work/INTERVAL.rsid - | \
  awk '$6==$10'
) > INF1.merge.coloc

function plotSMR()
{
  export bfile=$1
  export ma=$2
  export probe=$3
  export out=$4
  smr --bfile ${bfile} --gwas-summary ${ma} --beqtl-summary CAGE_snpid --gene-list plot/glist_hg19_refseq.txt \
      --probe ${probe} --probe-wind 500 \
      --plot --out ${out}
}
