#!/usr/bin/bash

function INTERVAL()
{
# init
  export rnaseq=tensorqtl_allSNPs_MAF0.005_merged_annotated.csv
  export rnaseq=tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated.csv
  grep -w ${rsid} ${rnaseq}
  zgrep ENSG00000256433 ${INF}/work/ensGtp.txt.gz | \
  cut -f2 | \
  zgrep -f - ${INF}/work/ensemblToGeneName.txt.gz
# region
  awk -vchr=${chr} -vpos=${pos} -vd=$((${flank_kb}*1000)) 'BEGIN{print chr,pos-d,pos+d}' > st.tmp
# LocusZoom plot
  read chr start end < st.tmp
  awk -vFS="," -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} 'NR==1 || ($5==chr && $6>=start && $6<=end && index($0,gene)>0)' ${rnaseq} | \
  tr "," "\t" > LTBR.lz
# export b1=$(cut -f6 LTBR.lz| sed '1d' | sort -k1,1n | awk 'NR==1')
# export b2=$(cut -f6 LTBR.lz| sed '1d' | sort -k1,1n | awk 'END{print}')
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal LTBR.lz --delim tab title="INTERVAL-LTBR" \
            --markercol variant_id --pvalcol pval --chr ${chr} --start ${b1} --end ${b2} \
            --no-date --plotonly --prefix=INTERVAL --rundir .
  mv INTERVAL_chr${chr}_${bracket}.pdf INTERVAL-LTBR-cis.pdf
}

function eQTLGen()
# https://www.eqtlgen.org/trans-eqtls.html
# https://www.eqtlgen.org/cis-eqtls.html
{
  export cis=2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz  
  export trans=2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
  zgrep -w ${rsid} ${cis} | \
  grep LTBR | \
  awk -vchr=${chr} -vpos=${pos} -vd=$((${flank_kb}*1000)) 'BEGIN{print chr,pos-d,pos+d}' > st.tmp
  read chr start end < st.tmp
  (
    gunzip -c ${cis} | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} 'NR==1 || ($3==chr && $4>=start && $4<=end && index($0,gene)>0)'
    gunzip -c ${trans} | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} 'NR==1 || ($3==chr && $4>=start && $4<=end && index($0,gene)>0)'
  ) > eQTLGen.lz
# export b1=$(cut -f4 eQTLGen.lz| sed '1d' | sort -k1,1n | awk 'NR==1')
# export b2=$(cut -f4 eQTLGen.lz| sed '1d' | sort -k1,1n | awk 'END{print}')
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal eQTLGen.lz --delim tab title="eQTLGen-LTBR" \
            --markercol SNP --pvalcol Pvalue --chr ${chr} --start ${b1} --end ${b2} \
            --no-date --plotonly --prefix=eQTLGen --rundir .
  mv eQTLGen_chr${chr}_${bracket}.pdf eQTLGen-LTBR-cis.pdf
}

function SCALLOP()
{
  awk -vchr=${chr} -vpos=${pos} -vd=$((${flank_kb}*1000)) 'BEGIN{print chr, pos-d, pos+d}' > st.tmp
  read chr start end < st.tmp
  (
    awk -vOFS="\t" 'BEGIN{print "MarkerName","P-value","Weight"}'
    gunzip -c ${INF}/METAL/TNFB-1.tbl.gz | \
    awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} '$1 == chr && $2 >= start && $2 <= end {print $3,10^$12,$18}' | \
    sort -k1,1 | \
    join <(awk -vchrom=chr${chr} 'index($0,chrom)>0' INTERVAL.rsid) - | \
    awk -vOFS="\t" '{print $2, $3, $4}'
  ) > TNFB.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal TNFB.lz --delim tab title="SCALLOP-TNFB" \
            --chr ${chr} --start ${b1} --end ${b2} --gwas-cat whole-cat_significant-only \
            --no-date --plotonly --prefix=TNFB --rundir .
  mv TNFB_chr${chr}_${bracket}.pdf SCALLOP-TNFB-cis.pdf
}

cd work
export chr=12
export pos=6514963
export gene=LTBR
export rsid=rs2364485
export flank_kb=1000
export b1=6300000
export b2=6700000
export bracket=${b1}-${b2}

tabix ${INF}/METAL/gwas2vcf/TNFB.tsv.gz ${chr}:${bracket} > TNFB.tbx

module load python/2.7

INTERVAL
eQTLGen
SCALLOP
cd -
