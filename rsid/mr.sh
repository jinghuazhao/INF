#!/usr/bin/bash

cd work
if [ ! -f INF1.merge.genes ]; then
   cut -f3,8,9,10 ${INF}/doc/olink.inf.panel.annot.tsv | grep -v BDNF | sed 's/"//g' | sort -k1,1 | join -12 inf1.tmp - > INF1.merge.genes
fi
if [ ! -d mr ]; then mkdir mr; fi

function MR_dat()
{
cut -f3 INF1.METAL | sed '1d' | sort | uniq | grep -w -f - INF1.merge.genes | \
parallel -j5 -C' ' '
  echo --- {2} ---
  gunzip -c ${INF}/METAL/{2}-1.tbl.gz | \
  cut -f1-6,10-12,18 | \
  awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vlogp=-5.45131 -vsuffix=${suffix} "
        (suffix==\"cis\" && \$1==chr && \$2>=start-M && \$2 <= end+M && \$9<=logp) || (suffix==\"pan\" && \$9<=logp)
      " > mr/{2}-${suffix}.mri
  (
    echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
    awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" mr/{2}-${suffix}.mri | \
    sort -k3,3 | \
    join -23 INTERVAL.rsid - | \
    awk -v prot={2} "{\$1=prot;print}" | \
    tr " " "\t"
  ) | gzip -f > mr/{2}-${suffix}.mrx
'
}
for type in cis pan; do export suffix=${type}; MR_dat; done
R --no-save <<END
  url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
  efo <- subset(openxlsx::read.xlsx(url, sheet="EFO", colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:4), rows=c(1:78)),!is.na(MRBASEID))
  write.table(efo, file="efo.txt",quote=FALSE,row.names=FALSE,sep="\t")
END
parallel --env INF -C' ' '
  export MRBASEID={1}; 
  export prot={2}; 
  export type={3}; 
  export prefix={1}-{2}-{3};
  echo ${prefix}
  R --no-save <${INF}/rsid/mr.R>/dev/null
  for f in result loo single; do cut -f1,2,5,6 --complement mr/${prefix}-${f}.txt | awk -vFS="\t" "NR==1||\$5<0.05" > mr/${prefix}-${f}.sig; done
  for f in result loo single; do export l=$(wc -l mr/${prefix}-${f}.sig | cut -d" " -f1); if [ ${l} -le 1 ]; then rm mr/${prefix}-${f}.sig; fi; done
' ::: $(awk -vFS="\t" 'NR>1 {print $4}' efo.txt) ::: $(sed '1d' INF1.merge | cut -f5 | sort -k1,1 | uniq) ::: cis pan
cd -

# uncomment if clumping outside TwoSampleMR:
# cut -f3 mr/{2}-${suffix}.mri > mr/{2}-${suffix}.mrs
# plink --bfile INTERVAL/cardio/INTERVAL --extract mr/{2}-${suffix}.mrs \
#       --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.01 --out mr/{2}-${suffix}
#   grep -w -f mr/{2}-${suffix}.prune.in mr/{2}-${suffix}.mri | \
#   join -23 <(zgrep "chr{3}" ${SUMSTATS}/snp150.snpid_rsid.gz) - | \
#   awk "{\$3=\"chr\"\$1\":\"\$2;print}" | \
#   join -23 -12 snp_pos - | \
