#!/usr/bin/bash

export pvalue=5e-8
export r2=0.8
export rsid=${INF}/ps/INF1.merge.snp

awk 'NR>1{print $8,$9}' ${INF}/work/INF1.merge | \
sort -k1,1n -k2,2n | \
uniq | \
awk '{print "chr" $1 ":" $2}' > ${rsid}
for catalogue in eQTL pQTL mQTL methQTL GWAS
do
  export catalogue=${catalogue}
  R --no-save -q < ${INF}/csd3/ps.R > ${INF}/ps/INF1.merge.${catalogue}.log
done

# gene
R --no-save -q < ${INF}/csd3/ps.gene.R > ${INF}/ps/INF1.merge.gene.log

# CLI
phenoscanner -s T -c All -x EUR -p ${pvalue} --r2 ${r2} -i ${rsid} -o ${INF}/ps/INF1
