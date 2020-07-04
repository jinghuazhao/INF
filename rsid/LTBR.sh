#!/usr/bin/bash

export rnaseq=tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated.csv
export rsid=rs2364485

function init()
{
  grep -w ${rsid} ${rnaseq}
  zgrep ENSG00000256433 ${INF}/work/ensGtp.txt.gz | \
  cut -f2 | \
  zgrep -f - ${INF}/work/ensemblToGeneName.txt.gz
}

cd work
init
export flank_kb=2000
grep -w ${rsid} ${rnaseq} | \
grep LTBR | \
awk -vFS="," -vd=$((${flank_kb}*1000)) '{print $5,$6-d,$6+d}' > st.tmp
echo LTBR | \
parallel --env rsid --env flank_kb -j1 -C' ' '
   read chrom start end < st.tmp
   (
     awk -vOFS="\t" "BEGIN{print \"MarkerName\",\"P-value\", \"Weight\"}"
     awk -vFS="," -vOFS="\t" -vchr=$chrom -vstart=$start -vend=$end \
         "(\$5 == chr && \$6 >= start && \$6 <= end)" ${rnaseq} | \
     sort -k6,6n | \
     awk "{print \$1,\$3,1}"
   )  > {}.lz
   rm -f ld_cache.db
   locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
             --markercol MarkerName --pvalcol P-value --refsnp ${rsid} --flank ${flank_kb}kb \
             --no-date --plotonly --prefix={} --rundir .
   pdftopng -r 300 {}_${refsnp}..pdf {}_${refsnp}
   mv {}_${refsnp}-000001.png {}_${refsnp}-1.png
   mv {}_${refsnp}-000002.png {}_${refsnp}-2.png
'
cd -
