# 1-11-2018 JHZ

export INTERVAL=/scratch/jp549/olink-merged-output

function CD6()
{
  gunzip -c $INTERVAL/INTERVAL_inf1_CD6___Q8WWJ7_chr_merged.gz | \
  awk -vOFS="\t" '(NR>1){
       SNP=$2;
       chr=$3; sub(/0/,"",chr);
       pos=$4;
       a1=$5;
       a2=$6;
       N=$18;
       EAF=(0.5*$15+$16)/N;
       beta=$24;
       se=$25;
       p=$22;
       print SNP,a2,a1,EAF,beta,se,p,N,chr,pos;
  }' | \
  sort -k9,9n -k10,10n > CD6
}

export BGEN_DIR=/scratch/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_typed_only
export BGEN=$BGEN_DIR/interval_olink_subset_unrelated_4994_pihat_0.1875_autosomal_typed_only

qctool -g $BGEN.bgen -s $BGEN.sample -snp-stats -osnp INTERVAL.snpstats

export SCRIPT=/scratch/jp549/analyses/interval_subset_olink/inf1/r2/outlier_in/pcs1_3
export BS=/scratch/jp549/apps/bram-scripts
export INF=/scratch/curated_genetic_data/phenotypes/interval/high_dimensional_data/Olink_proteomics_inf/gwasqc/olink_qcgwas_inf.csv
export REF=/scratch/curated_genetic_data/reference_files/interval/

function NOTE()
{
ls $INTERVAL/*gz | \
grep inf1 | \
xargs -l basename | \
sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g' | \
cut -d'_' --output-delimiter=' ' -f1,4 | \
grep -w CD6 | \
parallel -j2 --env INTERVAL -C' ' '
  gunzip -c $INTERVAL/INTERVAL_inf1_{1}___{2}_chr_merged.gz
'
}
