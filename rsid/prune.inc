#!/bin/bash

export job=$1
export TMPDIR=${HPC_WORK}/work
export list=${INF}/work/INF1.merge-rsid
export Start=$(awk -v job=${job} 'NR==job+1{print $2}' ${list})
export End=$(awk -v job=${job} 'NR==job+1{print $3}' ${list})
export p=$(awk -v job=${job} 'NR==job+1{print $5}' ${list})
export r=$(awk -v job=${job} 'NR==job+1{print $6}' ${list})
export pr=${p}-${r}
export chr=$(awk -v job=${job} 'NR==job+1{print $8}' ${list})
export pos=$(awk -v job=${job} 'NR==job+1{print $9}' ${list})
export flanking=1e6
export start=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
export end=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{print pos+flanking}')

function prune()
{
  module load plink/2.00-alpha
  plink2 --bfile ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr} \
         --chr ${chr} --from-bp ${start} --to-bp ${end} \
         --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.8 --out ${pr}
  if [ $(grep -w ${r} ${pr}.prune.in | wc -l) -eq 0 ]; then
     export i=$(grep -w -f ${pr}.prune.in ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr}.bim | \
     awk -vpos=${pos} 'function abs(x) {if (x<0) return -x; else return x;}
         {d=$4-pos;print $1, $2, $4, d}' | sort -r -k4,4n | awk 'NR==1 {print $2}')
     sed -i 's/'"$i"'/'"$r"'/g' ${pr}.prune.in
  fi
  if [ ${chr} -eq 19 ]; then
     sort ${pr}.prune.in | \
     join -v1 - ${INF}/work/NLRP2 > ${pr}.prune
  else
     sort ${pr}.prune.in > ${pr}.prune
  fi
  rm ${pr}.prune.in ${pr}.prune.out
}

cd ${INF}/prune
prune
cd -
