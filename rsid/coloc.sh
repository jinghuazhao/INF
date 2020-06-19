#!/usr/bin/bash

export trait=$1
export tissue=Whole_Blood.P01
export P01=$(echo ${tissue} | sed 's/.P01//')
export FUSION=${HPC_WORK}/fusion_twas

cd work
(
  awk -v OFS='\t' 'BEGIN{print "SNP","A1","A2","N","CHISQ","Z"}'
  gunzip -c ${INF}/METAL/${trait}-1.tbl.gz | \
  awk -v OFS='\t' 'NR>1{print $3,toupper($4),toupper($5),$18,($10/$11)^2,$10/$11,$1,$2}' | \
  sort -k1,1 | \
  join <(cat INTERVAL.rsid | tr ' ' '\t') - -t$'\t' | \
  sort -k8,8n -k9,9n | \
  cut -d' ' -f2-7 | \
) > ${trait}.sumstats
ln -sf ${FUSION}/WEIGHTS WEIGHTS
(
  awk -v OFS='\t' 'BEGIN{print "PANEL","FILE","ID","CHR","P0","P1","HSQ","BEST.GWAS.ID","BEST.GWAS.Z",
                  "EQTL.ID","EQTL.R2","EQTL.Z","EQTL.GWAS.Z","NSNP","NWGT",
                  "MODEL","MODELCV.R2","MODELCV.PV","TWAS.Z","TWAS.P",
                  "COLOC.PP0","COLOC.PP1","COLOC.PP2","COLOC.PP3","COLOC.PP4"}'
  (
    seq 22 | \
    parallel --env FUSION --env trait --env tissue -C' ' '
    Rscript ${FUSION}/FUSION.assoc_test.R \
            --sumstats ${trait}.sumstats \
            --weights ${FUSION}/WEIGHTS/${tissue}.pos \
            --weights_dir WEIGHTS \
            --ref_ld_chr ${FUSION}/LDREF/1000G.EUR. \
            --chr {} \
            --coloc_P 1e-4 \
            --GWASN 15335 \
            --out ${trait}-${tissue}-{}.dat
    cat ${trait}-${tissue}-{}.dat | \
    if [ {} -eq 1 ]; then cat; else awk "NR>1"; fi
    rm ${trait}-${tissue}-{}.dat
    '
  ) | \
  grep -v -e WARNING -e skipped -e complete -e consider -e MHC -e TWAS | \
  sed "s|WEIGHTS/${P01}/${P01}.||g;s/.wgt.RDat//g" | \
  sort -k4,4n -k5,5n
) | \
awk '$NF!="NA"' > ${trait}-${P01}-coloc.dat
cd -
