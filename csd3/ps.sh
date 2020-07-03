#!/usr/bin/bash

export rt="work/INF1.merge"
awk 'NR>1{print $8,$9}' ${rt}.merge | \
sort -k1,1n -k2,2n | \
uniq | \
awk '{print "chr" $1 ":" $2}' > ${rt}.snp

export rsid=${rt}.snp
export pvalue=1.5e-11
export r2=0.7
for catalogue in eQTL pQTL mQTL methQTL GWAS
do
  export catalogue=${catalogue}
  R --no-save -q < csd3/ps.R > work/INF1.merge.${catalogue}.log
done

# gene
R --no-save -q < csd3/ps.gene.R > work/INF1.merge.gene.log

phenoscanner -s T -c All -x EUR -p 0.0000001 --r2 0.6 -i INF1.merge.snp -o INF1
