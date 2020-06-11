#!/usr/bin/bash

export trait=$1
export tissue=Brain_Hypothalamus.P01
export eQTL=${HOME}/genetics/VCF-liftover
export FUSION=${HPC_WORK}/fusion_twas

cd work

(
  awk -v OFS='\t' 'BEGIN{print "#chrom","chromStart","chromEnd","MarkerName","z"}'
  gunzip -c ${INF}/METAL/${trait}-1.tbl.gz | \
  awk -v OFS='\t' 'NR>1 {print "chr" $1,$2-1,$2,$3,$10/$11}'
) | \
bedtools intersect -a ${HOME}/FM-pipeline/1KG/EUR.bed -b - -wa -wb | \
awk -v OFS='\t' 'NR>1{print $8,$4,$9}' | \
sort -k1,1 | \
join <(cat INTERVAL.rsid | tr ' ' '\t') ${trait}.torus.zval -t$'\t'| \
cut -f2,3,4 | \
awk '{gsub(/region/,"",$2)};1' | \
sort -k2,2n | \
awk -vOFS='\t' '{print $1,"loc" $2,$3}' | \
gzip -f > ${trait}.torus.zval.gz

torus -d ${trait}.torus.zval.gz --load_zval -dump_pip ${trait}.gwas.pip
gzip ${trait}.gwas.pip

fastenloc -eqtl ${eQTL}/gtex_v8.eqtl_annot_rsid.hg19.vcf.gz -gwas ${trait}.gwas.pip.gz -thread 4 -prefix ${trait}
sort -grk6 ${trait}.enloc.sig.out

# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensemblToGeneName.txt.gz
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensGtp.txt.gz

for o in sig snp
do
  export o=${o}
  cp ${trait}.enloc.${o}.out ${trait}-${o}.out
  join -12 <(awk 'NR>1' ${trait}.enloc.${o}.out | awk '{split($1,a,":");print a[1]}' | zgrep -w -f - ensGtp.txt.gz | cut -f1,2 | sort -k2,2) \
           <(gunzip -c ensemblToGeneName.txt.gz | sort -k1,1) | \
  cut -d' ' -f2,3 | \
  parallel --env o -C' ' 'sed -i "s/{1}/{1}-{2}/g" ${trait}-${o}.out'
done

gunzip -c ${INF}/METAL/${trait}-1.tbl.gz | \
(
  awk -v OFS='\t' 'BEGIN{print "SNP","A1","A2","N","CHISQ","Z"}'
  awk -v OFS='\t' 'NR>1{print $3,toupper($4),toupper($5),$18,($10/$11)^2,$10/$11}'
) > ${trait}.sumstats

(
  awk -v OFS='\t' 'BEGIN{print "PANEL","FILE","ID","CHR","P0","P1","HSQ","BEST.GWAS.ID","BEST.GWAS.Z",
                  "EQTL.ID","EQTL.R2","EQTL.Z","EQTL.GWAS.Z","NSNP","NWGT",
                  "MODEL","MODELCV.R2","MODELCV.PV","TWAS.Z","TWAS.P",
                  "COLOC.PP0","COLOC.PP1","COLOC.PP2","COLOC.PP3","COLOC.PP4"}'
  seq 22 | \
  parallel --env FUSION --env trait --env tissue -C' ' '
  Rscript ${FUSION}/FUSION.assoc_test.R \
          --sumstats ${trait}.sumstats \
          --weights ${FUSION}/WEIGHTS/${tissue}.pos \
          --weights_dir ${FUSION}/WEIGHTS \
          --ref_ld_chr ${FUSION}/LDREF/1000G.EUR. \
          --chr {} \
          --coloc_P 5e-8 \
          --GWASN 15335 \
          --out ${trait}-${tissue}-{}.dat
  cat ${trait}-${tissue}-{}.dat | \
  if [ {} -eq 1 ]; then cat; else awk "NR>1"; fi
  rm ${trait}-${tissue}-{}.dat
  '
) | \
(
  grep -v -e WARNING -e skipped -e complete -e consider
) > ${trait}-${tissue}.dat
cd -
