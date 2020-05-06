#!/usr/bin/bash

module load ceuadmin/stata/15

export inf1=$(ls METAL/*-1.tbl.gz | sed 's|METAL/||g;s/-1.tbl.gz//g;s/\*//g')

for p in ${inf1[@]}
do
   export p=${p}
   echo ${p}
   (
     zcat METAL/${p}-1.tbl.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
     zcat sentinels/${p}.p.gz | \
     awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}'
   ) > ${p}.a
   export p0=$(awk 'NR>1&&10^$13==0' ${p}.a | wc -l | cut -d' ' -f1)
   (
      echo SNP A1 A2 freq b se p N
      if [ ${p0} -eq 0 ]; then
         awk 'NR>1{print $4,toupper($5),toupper($6),$7,$11,$12,10^$13,$19}' ${p}.a
         rm ${p}.a
      else
         awk -vOFS="\t" '{if(NR==1) print $0, "N"; else print $0,NR-1}' work/${p}.merged > ${p}.b
         bedtools intersect -wb -a ${p}.a -b ${p}.b > ${p}.ab
         stata -b -q csd3/ma.do 
         awk -f csd3/ma.awk ${p}.txt
         rm ${p}.a ${p}.b ${p}.ab ${p}.txt ${p}.dta
      fi
   ) > sentinels/${p}.ma
    export lines=$(wc -l sentinels/${p}.ma | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}.ma with $lines lines
      rm sentinels/${p}.ma
    fi
done
rm ma.log

# The Stata code implements an old trick (Sun et al. 2018) such that only one p==0 is kept and others reset to 1.
# However, the latter version of GCTA uses effect size and this is unnecessary.
