#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
module load gcc/6
function prototype()
{
R --no-save <<END
  options(width=200)
  library(gwasvcf)
  library(VariantAnnotation)
  library(gap)
  INF <- Sys.getenv("INF")
  prots <- unique(with(inf1,prot))
  for (prot in prots)
  {
    cat(prot,"\n")
    metal <- read.table(file.path(INF,"METAL",paste0(prot,"-1.tbl.gz")),as.is=TRUE,header=TRUE)
    out <- with(subset(metal,!grepl("<",Allele1)&grepl("<",Allele2)),
               create_vcf(chrom=Chromosome, pos=Position, nea=toupper(Allele2), ea=toupper(Allele1), snp=MarkerName, ea_af=Freq1,
                                effect=Effect, se=StdErr, pval=10^log.P., n=as.integer(N), name=prot)
               )
    writeVcf(out, file=file.path(INF,"METAL",paste0(prot,".vcf")),index=TRUE)
  }
END
)

for prot in $(cut -f1 work/inf1.tmp | grep -v BDNF | sort)
do
  echo ${prot}
  export prot=${prot}
  (
    head -13 ${INF}/rsid/gwasvcf.hdr | \
    sed 's|NAME|'"$prot"'|g'
    gunzip -c ${INF}/METAL/${prot}-1.tbl | \
    awk -vOFS="\t" 'NR>1{print $1,$2,$3,toupper($4),toupper($5),".","PASS",".","AF:ES:SE:LP:SS",$6 ":" $10 ":" $11 ":" (-$12) ":" $18}' | \
    grep -v "<" | \
    sort -k1,1n -k2,2n
  ) | \
  bgzip -f > ${INF}/METAL/vcf/${prot}.vcf.gz
  tabix -f ${INF}/METAL/vcf/${prot}.vcf.gz
done

# sbatch ${INF}/rsid/gwas2vcf.sb
