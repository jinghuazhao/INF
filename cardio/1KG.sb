#!/bin/bash
#SBATCH -J 1KGp3v5
#SBATCH -e 1KGp3v5.log
#SBATCH -o 1KGp3v5.%A
#SBATCH -a 1-22
#SBATCH -p long
#SBATCH -t 4-0:0
#SBATCH --export ALL

export p3v5=/scratch/public_databases/1000_Genomes_phase3v5a
export TMPDIR=/scratch/jhz22/tmp

if [ ! -f EUR.list ]; then
  grep EUR $p3v5/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > EUR.list
fi
bcftools view -S EUR.list -O v $p3v5/ALL.chr${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
vcftools --vcf - --freq --stdout | \
awk -vOFS="\t" '
# input
#10	60515	2	5008	C:0.9998	T:0.000199681
# output
#         SNP CHR    POS MINOR MAJOR MAF SOURCE               DATE_ADDED
#1 rs10399749   1  45162     T     C   0 hapmap Thu Nov 22 16:08:25 2018
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
gzip -f > $HOME/INF/work/1KGp3v5-${SLURM_ARRAY_TASK_ID}.txt.gz

function info()
{
  cd $p3v5
# 503 EUR sample
  grep EUR integrated_call_samples_v3.20130502.ALL.panel | wc -l
# 84802133 sites
  gunzip -c ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | wc -l
}
