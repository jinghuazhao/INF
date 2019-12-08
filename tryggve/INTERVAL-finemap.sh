#!/bin/bash

export TMPDIR=/data/jinhua/work
export INF=/home/jinhua/INF
export list=${INF}/merge/INF1.merge
export n=10000
export N=4994
export LDREF=INTERVAL

module load gcc/5.4.0 lapack/3.8.0 qctool/2.0.1
module load bgen/20180807
cd work
for job in 97 # `seq 180`
do
  # variables
  export p=$(awk -v job=${job} 'NR==job+1{print $5}' ${list})
  export r=$(awk -v job=${job} 'NR==job+1{print $6}' ${list})
  export chr=$(awk -v job=${job} 'NR==job+1{print $8}' ${list})
  export pos=$(awk -v job=${job} 'NR==job+1{print $9}' ${list})
  export start=$(awk -v job=${job} 'NR==job+1{print $2}' ${list})
  export end=$(awk -v job=${job} 'NR==job+1{print $3}' ${list})
  export flanking=1e6
  export start=$(awk -vpos=${start} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
  export end=$(awk -vpos=${end} -vflanking=${flanking} 'BEGIN{print pos+flanking}')
  export pr=${p}-${r}
  export sumstats=${INF}/sumstats/INTERVAL/INTERVAL.${p}.gz

  # z0
  # SNPID   CHR     POS     STRAND  N       EFFECT_ALLELE   REFERENCE_ALLELE        CODE_ALL_FQ     BETA    SE      PVAL    RSQ     RSQ_IMP IMP
  (
    echo rsid chromosome position allele1 allele2 maf beta se
    zcat ${sumstats} | cut -f1-3,6-10 | sed '1d;s/\t/ /g' | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} '
    {
      if ($6 < 0.5) maf = $6; else maf = 1-$6
      if ($2==chr && $3 >= start && $3 < end && maf > 0 && maf <= 0.5 && $7 != "NA" && $8 != "NA") print $1, $2, $3, toupper($5), toupper($4), maf, $7, $8
    } ' | \
    sort -k1,1 | \
    join ${INF}/bgen/${pr} -
  ) > ${pr}.z0

  awk 'NR > 1{print $1} ' ${pr}.z0 > ${pr}.incl

  # bgen
  qctool -g ${INF}/bgen/${pr}.bgen -og ${pr}.bgen -ofiletype bgen -incl-rsids ${pr}.incl

  # bgi
  bgenix -g ${pr}.bgen -index -clobber
  ln -sf ${pr}.bgen.bgi ${pr}.bgi

  # z
  (
    join <(awk 'NR > 1' ${pr}.z0 | cut -d' ' -f1,4,5  | sort -k1,1) \
         <(bgenix -g ${pr}.bgen -list 2>&1 | awk 'NR>9 && NF==7' | cut -f2,6,7 | sort -k1,1) | \
    awk '{print $1, ($2!=$4)}'
  ) > ${pr}.flip
  (
    awk 'BEGIN {print "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "flip"}'
    join <(awk 'NR > 1' ${pr}.z0) ${pr}.flip | awk '{if($9==1) {t=$4;$4=$5;$5=t};$7=-$7; print}'
  ) > ${pr}.z

  # ldstore 1.1 and finemap 1.3.1.
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
                 --corr-config 0.95 --corr-group 0.99 --group-snps

# xlsx
  gzip -f ${pr}.ld
  R -q --no-save < ${INF}/csd3/finemap.R > ${pr}-finemap.log

done
cd -
