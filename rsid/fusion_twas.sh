#!/usr/bin/bash

function init()
{
  if [ ! -d ${INF}/FUSION ]; then mkdir ${INF}/FUSION; fi
  parallel --env FUSION --env INF -C' ' 'awk -vOFS="\t" "
  {
    if (\$5 < \$6) snpid=\"chr\"\$1\":\"\$4\"_\"\$5\"_\"\$6;
    else snpid=\"chr\"\$1\":\"\$4\"_\"\$6\"_\"\$5
    print snpid, \$2
  }" ${FUSION}/LDREF/1000G.EUR.{}.bim > ${INF}/FUSION/1000G.EUR.{}.snpid' ::: $(seq 22)
  Rscript -e "write.table(pQTLtools::inf1,file=file.path(Sys.getenv('INF'),'FUSION','INF1.tsv'),row.names=FALSE,quote=FALSE,sep='\t')"
}
export GCTA=${HPC_WORK}/bin/gcta64
export PLINK="${HPC_WORK}/bin/plink-1.9 --allow-no-sex"
export GEMMA=${HPC_WORK}/bin/gemma
export FUSION=${HPC_WORK}/fusion_twas
sed '1d' ${INF}/FUSION/INF1.tsv | cut -f2,5-8 --output-delimiter=' ' | sort -k2,2n | \
parallel --env FUSION --env INF -C' ' '
  export prot={1}
  export gene={2}
  export chr={3}
  export P0=$(awk -v start={4} "BEGIN{if(start<1e6) print 1; else print start-1e6}")
  export P1=$(expr {5} + 1000000)
  export col=$(awk "\$1 == ENVIRON[\"prot\"] {print NR}" ${INF}/csd3/prot.list)
  $PLINK --bfile ${INF}/INTERVAL/cardio/INTERVAL --make-bed --out ${INF}/FUSION/${prot}-${gene} \
         --chr ${chr} --from-bp ${P0} --to-bp ${P1} --maf 0.01 > /dev/null
  cut -d" " -f1,2,${col} ${INF}/INTERVAL/cardio/INTERVAL.pheno > ${INF}/FUSION/${prot}-${gene}.pheno
  Rscript ${FUSION}/FUSION.compute_weights.R --verbose 2 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA \
          --bfile ${INF}/FUSION/${prot}-${gene} \
          --pheno ${INF}/FUSION/${prot}-${gene}.pheno --covar ${INF}/INTERVAL/cardio/INTERVAL.qcovar \
          --tmp ${INF}/FUSION/${prot}-${gene} --out ${INF}/FUSION/${prot}-${gene} \
          --models top1,blup,bslmm,lasso,enet
'
sed '1d' ${INF}/FUSION/INF1.tsv | cut -f2,5-8 --output-delimiter=' ' | sort -k2,2 | \
join -o1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4 -e NA -a2 -12 -24 - <(sort -k4,4 ${INF}/csd3/glist-hg19)
