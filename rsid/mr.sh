#!/usr/bin/bash

function MR_dat()
{
cut -f3,8,9,10 doc/olink.inf.panel.annot.tsv | grep -v BDNF | sed 's/"//g' | sort -k1,1 | join -12 work/inf1.tmp - | \
parallel -j1 -C' ' '
  echo --- {2} ---
  gunzip -c METAL/{2}-1.tbl.gz | \
  cut -f1-6,10-12,18 | \
  awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vlogp=-5.45131 "
      {
        if(ENVIRON["suffix"]=="cis") if(\$1==chr && \$2>=start-M && \$2 <= end+M && \$9<=logp) print; else if(\$9<=logp) print
      }" > work/{2}.mri-${suffix}
  cut -f3 work/{2}.mri-${suffix} > work/{2}.mrs-${suffix}
  plink --bfile INTERVAL/cardio/INTERVAL --extract work/{2}.mrs-${suffix} --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.1 \
        --out work/{2}-${suffix}
  if [ -f work/{2}-${suffix} ]; then
  (
    echo -e "rsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\nN"
    grep -w -f work/{2}-${suffix}.prune.in work/{2}.mri-${suffix} | \
    awk "{\$3=\"chr\"\$1\":\"\$2;\$8=-\$8;print}" | \
    sort -k3,3 | \
    join -23 -12 work/snp_pos - | \
    cut -d" " -f1 --complement | \
    tr " " "\t"
  ) | gzip -f > work/{2}.mrx-${suffix}
  fi
'
}

for type in cis pan; do export suffix=${type}; MR_dat; done

export outcomes="ukb-b-20208"
cut -f3,8,9,10 doc/olink.inf.panel.annot.tsv | grep -v BDNF | sed 's/"//g' | sort -k1,1 | join -12 work/inf1.tmp - | \
parallel -j1 -C' ' 'export prot={2}; R --no-save <rsid/mr.R'
