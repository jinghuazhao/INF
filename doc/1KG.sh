# 16-11-2018 JHZ

(
  echo -e "snpid\trsid"
  for i in $(seq 22)
  do
    zgrep -v -w CHROM ${i}.csv.gz | \
    awk -vFS="," -vOFS="\t" '
    {
      OFS="\t"
      CHR=$1
      POS=$3
      a1=$4
      a2=$5
      rsid=$2
      if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
      else snpid="chr" CHR ":" POS "_" a1 "_" a2
      print snpid, rsid
   }'
  done
) | \
gzip -f > 1KG.tsv.gz

# x.csv.gz
# 1 CHR
# 2 SNP
# 3 POS
# 4 A1
# 5 A2
