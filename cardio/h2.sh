#!/bin/bash

export TMPDIR=/scratch/jhz22/tmp
export rt=/scratch/jhz22/data/INTERVAL/INTERVAL
export s=/home/jp549/post-doc/genetics/r2-test/sample_files/o5000-inf1-outlier_out-r2.sample

cut -d' ' -f1-2,5-7 ${s} | awk 'NR>3' > ${rt}.covar
cut -d' ' -f1-2,4,8-28 ${s} | awk 'NR>3' > ${rt}.qcovar
cut -d' ' -f1-2,29- ${s} | awk 'NR>3' > ${rt}.pheno

plink --bfile ${rt} --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out $rt
plink --bfile ${rt} --extract ${rt}.prune.in --make-bed --out ${rt}.prune

plink --bfile ${rt}.prine \
      --make-grm-bin \
      --threads 2 \
      --out ${rt}
