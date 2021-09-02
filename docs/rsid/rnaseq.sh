#!/usr/bin/bash

Rscript -e 'write.table(pQTLtools::inf1,file=file.path(Sys.getenv("INF"),"work","inf1.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")'
export rnaseq=${INF}/work/tensorqtl_allSNPs_MAF0.005_merged_annotated.csv
export TMPDIR=${HPC_WORK}/work
join -a1 <(grep -f ${INF}/sentinels/INF1.jma-rsid.cis ${INF}/sentinels/INF1.jma-rsid | cut -f1,3-6 --output-delimiter=' ' | sort -k1,1) \
         <(cut -f2,5 --output-delimiter=' ' ${INF}/work/inf1.txt | sed '1d' | sort -k1,1) | \
awk '{print $2":"$4,$3,$6}' | \
sort -k1,1 | \
join - <(cut -d, -f1,3,5,6,13 ${rnaseq} | awk -vFS="," 'NR>1{print $3":"$4,$1,$5,$2}' | sort -k1,1) > ${INF}/work/rnaseq.dat
awk '$2==$4 && $3==$5 && $6<5e-8'
