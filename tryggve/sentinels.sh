# 3-6-2019 JHZ

module load bedtools/2.27.1

# extract all significant SNPs
ls METAL/*-1.tbl.gz | \
sed 's|METAL/||g;s/-1.tbl.gz//g' | \
parallel -j3 -C' ' '
(
# zcat METAL/{}-1.tbl.gz | head -1
  zcat METAL/{}-1.tbl.gz | awk "NR>1 && length(\$4)==1 && length(\$5)==1 && \$12<5e-10" | sort -k1,1n -k2,2n
) | gzip -f > work/{}.p.gz 

# removing those in high LD regions'
for p in $(ls METAL/*-1.tbl.gz | sed 's|METAL/||g;s/-1.tbl.gz//g')
do
  (
    zcat METAL/${p}-1.tbl.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
    zcat work/${p}.p.gz | \
    awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}'
  ) | bedtools subtract -header -a - -b tryggve/high-LD-regions-hg19.bed > work/${p}.p
# echo $(zcat ${p}.p.gz | wc -l) $(wc -l ${p}.p)
  export lines=$(wc -l work/${p}.p|cut -d' ' -f1)
  if [ $lines -eq 1 ]; then
    echo removing ${p} with $lines lines
    rm work/${p}.p
  fi
done

# find sentinels
for prot in $(ls work/*.p | sed 's|work/||g;s|\.p||g')
do 
  export prot=${prot}
  echo ${prot}
  R --no-save -q < tryggve/sentinels.R > work/${prot}.o
done
cd work
(
  awk -vOFS="," 'BEGIN{print "prot","CHR","BP","SNP","l","u","d","log10p","Groupid", "Type"}'
  awk -vFS="," -vOFS="," '!/option/{
      SNPID=$2
      split(SNPID,a,":")
      split(a[2],b,"_")
      gsub(/chr/,"",a[1])
      $1=$1 "," a[1] "," b[1]
      print
  }' *.o
) > INF1.sentinels
cd -
