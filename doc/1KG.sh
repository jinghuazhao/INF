# 2-12-2018 JHZ

# to take on results from 1KG.sb
(
  awk -vOFS="\t" 'BEGIN{print "SNP","CHR","POS","MINOR","MAJOR","MAF"}'
  for chr in $(seq 22)
  do
    zgrep -v -w CHR 1KGp3v5-${chr}.txt.gz 
  done
) | \
gzip -f > 1KGp3v5.txt.gz

R --no-save -q <<END
z <- gzfile("1KGp3v5.txt.gz")
allele_ref_std <- read.table(z,header=TRUE,as.is=TRUE)
save(allele_ref_std,file="1KGp3v3.RData")
END

# to work on data from LocusZoom 1.4
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
