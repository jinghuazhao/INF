# 4-4-2020 JHZ

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
