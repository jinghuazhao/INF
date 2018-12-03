# 3-12-2018 JHZ

# allele frequencies from 1KG phase 3
sbatch -wait doc/1KG.slurm
(
  awk -vOFS="\t" 'BEGIN{print "SNP","CHR","POS","MINOR","MAJOR","MAF"}'
  for chr in $(seq 22); do zgrep -v -w CHR 1KGp3v5-${chr}.txt.gz ; done
) | \
gzip -f > 1KGp3v5.tsv.gz

R --no-save -q <<END
  z <- gzfile("1KGp3v5.tsv.gz")
  allele_ref_std <- read.table(z,header=TRUE,as.is=TRUE)
  save(allele_ref_std,file="1KGp3v5.RData")
END

# snpid-rsid from genotype data
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
gunzip -c  /data/jinhua//1KGp3/1KG.tsv.gz | cut -f1 | awk 'NR>1' | sort > 1KG.snpid

# x.csv.gz
# 1 CHR
# 2 SNP
# 3 POS
# 4 A1
# 5 A2

# to pick SNPIDs not in the reference panel

gunzip -c /data/jinhua/data/1KG/1KG.tsv.gz | \
cut -f1 | \
awk 'NR>1' | \
sort -k1,1 > 1KG.snpid

export protein=CD6
grep -w $protein METAL/METAL.tmp | \
parallel -j1 -C' ' '
   gunzip -c {1} | \
   cut -f1 | \
   awk "NR>1"| \
   sort -k1,1 | \
   join -j1 -v1 - 1KG.snpid > {2}.v1
' 

