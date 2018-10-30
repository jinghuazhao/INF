export SRC=/scratch/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_typed_only
export BGEN_SAMPLE=$SRC/interval_olink_subset_unrelated_4994_pihat_0.1875_autosomal_typed_only
export SCRIPT=/scratch/jp549/analyses/interval_subset_olink/inf1/r2/outlier_in/pcs1_3
export INTERVAL=/scratch/jp549/olink-merged-output
export BS=/scratch/jp549/apps/bram-scripts

ls $INTERVAL/*gz | \
grep inf1 | \
xargs -l basename | \
sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g' | \
cut -d'_' --output-delimiter=$'\t' -f1,4
