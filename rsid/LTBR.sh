#!/usr/bin/bash

export rnaseq=tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated.csv
export rsid=rs2364485
export flank_kb=15000

function interval()
{
# init
  grep -w ${rsid} ${rnaseq}
  zgrep ENSG00000256433 ${INF}/work/ensGtp.txt.gz | \
  cut -f2 | \
  zgrep -f - ${INF}/work/ensemblToGeneName.txt.gz
# region
  grep -w ${rsid} ${rnaseq} | \
  grep LTBR | \
  awk -vFS="," -vd=$((${flank_kb}*1000)) '{print $5,$6-d,$6+d}' > st.tmp
# LocusZoom plot
  echo LTBR | \
  parallel --env rsid --env flank_kb -j1 -C' ' '
     read chrom start end < st.tmp
     (
       awk -vOFS="\t" "BEGIN{print \"MarkerName\",\"P-value\", \"Weight\"}"
       awk -vFS="," -vchr=$chrom -vstart=$start -vend=$end \
           "(\$5 == chr && \$6 >= start && \$6 <= end)" ${rnaseq} | \
       tr "," " " | \
       sort -k6,6n | \
       awk -vOFS="\t" "{print \$1,\$3,1}"
     )  > {}.lz
     rm -f ld_cache.db
     locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
               --markercol MarkerName --pvalcol P-value --refsnp ${rsid} --flank ${flank_kb}kb \
               --no-date --plotonly --prefix={} --rundir .
     export f={}_${rsid}
     pdftopng -r 300 ${f}.pdf ${f}
     mv ${f}-000001.png ${f}-1.png
     mv ${f}-000002.png ${f}-2.png
  '
}

# https://www.eqtlgen.org/trans-eqtls.html
# https://www.eqtlgen.org/cis-eqtls.html

function eQTLGen()
{
  export cis=2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz  
  export trans=2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
  zgrep -w ${rsid} ${cis} | \
  grep LTBR | \
  awk -vd=$((${flank_kb}*1000)) '{print $3,$4-d,$4+d}' > st.tmp
  read chrom start end < st.tmp
  awk -vOFS="\t" "BEGIN{print \"MarkerName\",\"P-value\", \"Weight\"}" > eQTLGen.lz
  (
    gunzip -c ${cis} | \
    awk -vchr=$chrom -vstart=$start -vend=$end '$3==chrom && $4>=start && $4<=end'
    gunzip -c ${trans} | \
    awk -vchr=$chrom -vstart=$start -vend=$end '$3==chrom && $4>=start && $4<=end'
  ) | \
  sort -k4,4n | \
  awk -v OFS="\t" '{print $2,$1,$13}' >> eQTLGen.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal eQTLGen.lz \
            --markercol MarkerName --pvalcol P-value --refsnp ${rsid} --flank ${flank_kb}kb \
            --no-date --plotonly --prefix=eQTLGen --rundir .
  export f=eQTLGen_${rsid}
  pdftopng -r 300 ${f}.pdf ${f}
  mv ${f}-000001.png ${f}-1.png
  mv ${f}-000002.png ${f}-2.png
}

cd work
eQTLGen
cd -
