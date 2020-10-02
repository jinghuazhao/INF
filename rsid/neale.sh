#!/usr/bin/bash

# This script extracts signals from a compressed ${src} file for a given ${pval} -- tuned for UKB sumstats from Neale lab
# The result file is ${pval}/p.signals
# Adapted from counterparts for multiple traits (proteins) with a trait variable in the output as a historical hallmark

export pval=5e-8
export src=ukb/30750_raw.gwas.imputed_v3.both_sexes.tsv.bgz
export tag=_nold

if [ ! -d ${pval} ]; then mkdir ${pval}; fi

function pgz()
# 1. extract all significant SNPs
{
  (
    zcat ${src} | head -1 | awk -vOFS="\t" '{print "chr","pos",$0}'
    zcat ${src} | awk -v p=${pval} -vOFS="\t" 'NR>1 && $11 <= p {split($1,a,":");print a[1],a[2],$0}' | sort -k1,1n -k2,2n
  ) | gzip -f > ${pval}/p.gz
}

function _HLA()
# 2. handling HLA
{
  (
    zcat ${src} | head -1 | awk -vOFS="\t" '{$1="Chrom" OFS "Start" OFS "End" OFS $1;print}'
    (
      zcat ${pval}/p.gz | \
      awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
      awk 'NR>1 && !($1 == "6" && $3 >= 25392021 && $3 < 33392022)'
      zcat ${pval}/p.gz | \
      awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
      awk '$1 == "6" && $3 >= 25392021 && $3 < 33392022' | \
      sort -k13,13g | \
      awk 'NR==1'
    ) | \
    sort -k1,1n -k2,2n -k3,3n | \
    awk -v OFS="\t" '{$1="chr" $1};1'
  ) > ${pval}/${tag}.p
  export lines=$(wc -l ${pval}/${tag}.p | cut -d' ' -f1)
  if [ $lines -eq 1 ]; then
    echo removing ${tag} with $lines lines
    rm ${pval}/${tag}.p
  fi
}

for cmd in pgz _HLA; do $cmd; done

(
  mergeBed -i ${pval}/${tag}.p -d 1000000 -c 14 -o min | \
  awk -v OFS="\t" -v trait=HbA1c '
  {
    if(NR==1) print "Chrom", "Start", "End", "P", "trait"
    print $0, trait
  }'
) > ${pval}/p.merged
(
  cut -f1-4,14 ${pval}/_nold.p | \
  bedtools intersect -a ${pval}/p.merged -b - -wa -wb | \
  awk '$4==$10' | \
  cut -f1-6,8-10 | \
  awk -v OFS="\t" '
  {
    if(NR==1) print "Chrom", "Start", "End", "P", "trait", "MarkerName", "CHR", "POS", "SNP", "P_check"
    $5=$5 OFS $6 ":" $7
    gsub(/chr/,"",$6)
    print
  }'
) | uniq > ${pval}/p.sentinels

module load gcc/6
R --no-save -q <<END
  pval <- Sys.getenv("pval")
  f <- file.path(pval,"p.sentinels")
  m <- read.table(f,header=TRUE,as.is=TRUE)
  dim(m)
  head(m)
  library(dplyr)
  t <- m %>% group_by(trait,Chrom,Start,End) %>% slice(which.min(P))
  t
  write.table(t,file=paste(pval,"p.merge",sep="/"),quote=FALSE,row.names=FALSE,sep='\t')
  P <- with(m,P)
  p <- table(P)[table(P)>1]
  print(p)
  merge <- read.delim(file.path(pval,"p.merge"),as.is=TRUE)
  m <- subset(merge,MarkerName!=".")
  cols <- c(1:5,9)
  write.table(m[,cols],file=file.path(pval,"p.signals"),row.names=FALSE,quote=FALSE,sep="\t")
END
