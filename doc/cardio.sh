export INTERVAL=/scratch/jp549/olink-merged-output

ls $INTERVAL/*gz | \
grep inf1 | \
xargs -l basename | \
sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g' | \
cut -d'_' --output-delimiter=$'\t' -f1,4

