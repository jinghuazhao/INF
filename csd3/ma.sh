#!/usr/bin/bash
# 5-11-2019 JHZ

export tag=_nold

export inf1=$(ls METAL/*-1.tbl.gz | sed 's|METAL/||g;s/-1.tbl.gz//g;s/\*//g')
for p in ${inf1[@]}
do echo $p
  (
    echo SNP A1 A2 freq b se p N
    (
      zcat METAL/${p}-1.tbl.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      zcat sentinels/${p}.p.gz | \
      awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
      awk '!($1 == "chr6" && $3 >= 25392021 && $3 < 33392022)'
    ) | \
    awk 'NR>1{print $4, toupper($5), toupper($6), $7, $11, $12, 10^$13, $19}'
  ) > sentinels/${p}.ma
  export lines=$(wc -l sentinels/${p}.ma | cut -d' ' -f1)
  if [ $lines -eq 1 ]; then
     echo removing ${p}.ma with $lines lines
     rm sentinels/${p}.ma
  fi
done
