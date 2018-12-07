# 7-12-2018 JHZ

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

echo "--> 100Genomes phase 3 version 5 at Cardio"

export p3v5=/scratch/public_databases/1000_Genomes_phase3v5a
export TMPDIR=/scratch/jhz22/tmp
if [ ! -f EUR.list ]; then
  grep EUR $p3v5/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > EUR.list
fi
bcftools view -S EUR.list -O v $p3v5/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz | \
vcftools --vcf - --freq --stdout | \
awk -vOFS="\t" '
{
  if(NR==1) print "SNP","CHR","POS","MINOR","MAJOR","MAF";
  else
  {
     split($5,x,":"); a=x[1]; af=x[2];
     split($6,y,":"); b=y[1]; bf=y[2];
     if (a>b) snpid="chr" $1 ":" $2 "_" b "_" a;
         else snpid="chr" $1 ":" $2 "_" a "_" b
     if (af<bf) {minor=a;major=b;maf=af}
           else {minor=b;major=a;maf=bf}
     if (maf>0) print snpid, $1, $2, minor, major, maf
  }
}' | \
gzip -f > $HOME/INF/work/1KGp3v5-X.txt.gz
# bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT'
