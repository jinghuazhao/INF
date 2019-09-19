#!/bin/bash

export job=132
export TMPDIR=/data/jinhua/work
export INF=/home/jinhua/INF
export list=${INF}/work/INF1_nold.sentinels
export p=$(awk -v job=${job} 'NR==job+1{print $1}' ${list})
export chr=$(awk -v job=${job} 'NR==job+1{print $2}' ${list})
export pos=$(awk -v job=${job} 'NR==job+1{print $3}' ${list})
export r=$(awk -v job=${job} 'NR==job+1{print $4}' ${list})
export pr=${p}-${r}
export sumstats=${INF}/METAL/${p}-1.tbl.gz
export flanking=1e6
export start=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
export end=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{print pos+flanking}')
export bfile=${INF}/EUR
export sample=${INF}/work/o5000-inf1-outlier_in-r2.sample
export study=INTERVAL
export N=4994

cd work
# z0
(
  zcat ${sumstats} | head -1 | awk '{print $3,$1,$2,$5,$4,"MAF",$10,$11}'
  zcat ${sumstats} | \
  awk 'NR > 1' | \
  sort -k3,3 | \
  awk -vchr=${chr} -vstart=${start} -vend=${end} '
  {
    if ($1==chr && $2 >= start && $2 < end) {
    if ($6 < 0.5) maf = $6; else maf = 1-$6
    if (maf > 0 && maf <= 0.5 && $10 != "NA" && $11 != "NA") print $3, $1, $2, toupper($5), toupper($4), maf, $10, $11
  }
  } ' | \
  join ${pr} -
) > ${pr}.z0
awk 'NR > 1{print $1} ' ${pr}.z0 > ${pr}.incl

# bgen
module load gcc/5.4.0 lapack/3.8.0 qctool/2.0.1
qctool -g ${bfile}.bed -og ${pr}.bgen -ofiletype bgen -incl-rsids ${pr}.incl

# bgi
module load bgen/20180807
bgenix -g ${pr}.bgen -index -clobber
ln -sf ${pr}.bgen.bgi ${pr}.bgi

# z
(
  join <(awk 'NR > 1' ${pr}.z0 | cut -d' ' -f1,4,5  | sort -k1,1) \
       <(bgenix -g ${pr}.bgen -list 2>&1 | awk 'NR>9 && NF==7'| cut -f2,6,7 | sort -k1,1) | \
  awk '{print $1, ($2!=$4)}'
) > ${pr}.flip
(
  awk 'BEGIN {print "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "flip"}'
  join <(awk 'NR > 1' ${pr}.z0) ${pr}.flip | awk '{if($9==1) {t=$4;$4=$5;$5=t};$7=-$7; print}'
) > ${pr}.z

function ld_finemap()
# ldstore 1.1 and finemap 1.3.1.
{
  (
    echo "z;ld;snp;config;cred;log;n_samples"
    echo "${pr}.z;${pr}.ld;${pr}.snp;${pr}.config;${pr}.cred;${pr}.log;$N"
  ) > ${pr}.master
  ldstore_v1.1 --bcor ${pr}-1 --bgen ${pr}.bgen --n-threads 1
  ldstore_v1.1 --bcor ${pr}-1 --merge 1
  ldstore_v1.1 --bcor ${pr}-1 --matrix ${pr}.ld
  rm ${pr}-1_*
  mv ${pr}-1 ${pr}.bcor
  rm -rf ${pr}.cred* ${pr}.dose ${pr}.snp ${pr}.config
  finemap_v1.3.1 --sss --in-files ${pr}.master --log --n-causal-snps 10 \
               --corr-config 0.95 --corr-group 0.9999 --group-snps
}

function ld_finemap2()
# ldstore 2.0b
{
  (
    echo "z;ld;bgen;bgi;bcor;bdose;snp;config;cred;log;n_samples"
    echo "${pr}.z;${pr}.ld;${pr}.bgen;${pr}.bgi;${pr}.bcor;${pr}.bdose;${pr}.snp;${pr}.config;${pr}.cred;${pr}.log;$N"
  ) > ${pr}.master2
  ldstore_v2.0b --in-files ${pr}.master2 --write-bcor
  ldstore_v2.0b --in-files ${pr}.master2 --bcor-to-text
  finemap_v1.4 --sss --in-files ${pr}.master2 --log --n-causal-snps 10 \
               --corr-config 0.95 --corr-group 0.9999 --group-snps
}}

# xlsx
gzip -f ${pr}.ld
R -q --no-save < ${INF}/csd3/finemap.R > ${pr}-finemap.log

cd -
