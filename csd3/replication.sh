# 9-5-2020 JHZ

(
  gunzip -c sumstats/ARISTOTLE/*txt.gz | \
  head -1 | \
  awk -vOFS="\t" '{print "Protein",$0}'
  cut -f5,6 work/INF1.merge | \
  sed '1d' | \
  tr '\t' ' ' | \
  parallel -j8 -C' ' '
    echo {1}-{2}
    zgrep -w {2} sumstats/ARISTOTLE/ARISTOTLE.{1}.txt.gz;
    awk -vOFS="\t" -vprot={1} "{print prot,\$0}"
  '
) > work/INF1.merge.replication.txt

# snpid --> rsid
cd work
for f in INF1.merge.replication.txt
do
  awk -vOFS="\t" '{if(NF==1) printf $1 OFS; else print}' ${f} > ${f}-rsid
  (
  cat INF1.merge.rsid | \
  parallel --dry-run -C' ' "
    export s={1};
    export r={2};
    sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}-rsid
  "
  ) | bash
done
awk -vOFS="\t" '{if(NR>1) {split($1,a,"-");$1=a[1]};print}' INF1.merge.replication.txt-rsid | xsel -i
for p in 5e-10 1e-5 5e-2
do
  echo ${p}
  cut -f12 work/INF1.merge.replication.txt-rsid | awk -v p=${p} '$1<p{print $1}' | wc -l
done
cd -
