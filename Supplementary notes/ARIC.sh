#!/usr/bin/bash

export sumstats=~/rds/results/public/proteomics/ARIC

if [ -d ${INF}/ARIC ]; then mkdir ${INF}/ARIC; fi

cat <(sed '1d' ${sumstats}/seqid.txt | cut -f2) <(echo P12034) <(echo P30203) | \
grep -f - doc/olink.inf.panel.annot.tsv | \
cut -f3 | \
sed 's\"\\g' | \
grep -f - work/inf1.tmp | \
grep -v BDNF | \
cut -f1 | \
grep -f - work/INF1.b38

cat <(echo pos38) <(cut -f2 ${INF}/work/INF1.b38) | \
paste ${INF}/work/INF1.METAL - | \
sed '1d;s/chr[0-9]*:[0-9]*//' | \
awk '{$1="chr"$4":"$22$1;print $1,$2,$3,$20}' | \
sort -k4,4 | \
join -12 -24 <(sed '1d' ${sumstats}/seqid.txt | cut -f1,2 | sort -k2,2) - | \
cut -d' ' -f1 --complement | \
parallel -j12 -C' ' --env sumstats '
  grep -w {3} ${sumstats}/EA/{1}.PHENO1.glm.linear | \
  awk -vseqid={1} -vsnpid={2} -vrsid={3} -vprot={4} -vOFS="\t" "\$13<=5e-8{print seqid,snpid,rsid,prot,\$0}"
'
