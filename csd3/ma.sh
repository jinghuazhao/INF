#!/usr/bin/bash
# 6-11-2019 JHZ

module load ceuadmin/stata/15

export tag=_nold
export inf1=$(ls METAL/*-1.tbl.gz | sed 's|METAL/||g;s/-1.tbl.gz//g;s/\*//g')

for p in ${inf1[@]}
do echo ${p}
   (
     zcat METAL/${p}-1.tbl.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
     zcat sentinels/${p}.p.gz | \
     awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
     awk '!($1 == "chr6" && $3 >= 25392021 && $3 < 33392022)'
   ) > ${p}.a
   awk -vOFS="\t" '{if(NR==1) print $0, "N"; else print $0,NR-1}' work/${p}.merged > ${p}.b
   bedtools intersect -wb -a ${p}.a -b ${p}.b > ${p}.ab
   stata -b -q csd3/ma.do 
   (
     echo SNP A1 A2 freq b se p N
     awk -f csd3/ma.awk ${p}.txt
   ) > sentinels/${p}.ma
   rm ${p}.a ${p}.b ${p}.ab ${p}.txt ${p}.dta ${p}0.dta ${p}.log
   export lines=$(wc -l sentinels/${p}.ma | cut -d' ' -f1)
   if [ $lines -eq 1 ]; then
      echo removing ${p}.ma with $lines lines
      rm sentinels/${p}.ma
   fi
done
