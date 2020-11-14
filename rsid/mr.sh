#!/usr/bin/bash

if [ ! -f work/INF1.merge.genes ]; then
   cut -f3,8,9,10 doc/olink.inf.panel.annot.tsv | grep -v BDNF | sed 's/"//g' | sort -k1,1 | join -12 work/inf1.tmp - > work/INF1.merge.genes
fi
if [ ! -d work/mr ]; then mkdir work/mr; fi

function MR_dat()
{
cut -f3 work/INF1.METAL | sed '1d' | sort | uniq | grep -w -f - work/INF1.merge.genes | \
parallel -j5 -C' ' '
  echo --- {2} ---
  gunzip -c METAL/{2}-1.tbl.gz | \
  cut -f1-6,10-12,18 | \
  awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vlogp=-5.45131 -vsuffix=${suffix} "
        (suffix==\"cis\" && \$1==chr && \$2>=start-M && \$2 <= end+M && \$9<=logp) || (suffix==\"pan\" && \$9<=logp)
      " > work/mr/{2}-${suffix}.mri
  (
    echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
    awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" work/mr/{2}-${suffix}.mri | \
    sort -k3,3 | \
    join -23 work/INTERVAL.rsid - | \
    awk -v prot={2} "{\$1=prot;print}" | \
    tr " " "\t"
  ) | gzip -f > work/mr/{2}-${suffix}.mrx
'
}
for type in cis pan; do export suffix=${type}; MR_dat; done
R --no-save < ~/INF/rsid/mr.R > ~/INF/work/mr.log

# entirely independent of TwoSampleMR:
# cut -f3 work/mr/{2}-${suffix}.mri > work/mr/{2}-${suffix}.mrs
# plink --bfile INTERVAL/cardio/INTERVAL --extract work/mr/{2}-${suffix}.mrs \
#       --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.01 --out work/mr/{2}-${suffix}
#   grep -w -f work/mr/{2}-${suffix}.prune.in work/mr/{2}-${suffix}.mri | \
#   join -23 <(zgrep "chr{3}" ${SUMSTATS}/snp150.snpid_rsid.gz) - | \
#   awk "{\$3=\"chr\"\$1\":\"\$2;print}" | \
#   join -23 -12 work/snp_pos - | \
