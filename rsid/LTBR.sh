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
         "(\$5 == chr && \$6 >= start && \$6 <= end) {print \$1,\$3,5000}" ${rnaseq}
   )  > {}.lz
   rm -f ld_cache.db
   locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
             --markercol MarkerName --pvalcol P-value --refsnp ${rsid} --flank ${flank_kb}kb \
             --no-date --plotonly --prefix={} --rundir .
   mv {}-chr${chrom}_${start}-${end}.pdf {}.lz.pdf
   pdftopng -r 300 {}.lz.pdf {}
   mv {}-000001.png {}.lz-1.png
   mv {}-000002.png {}.lz-2.png
'
cd -
