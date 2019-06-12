# coloc

# prototype coloc pipeline
# https://bitbucket.org/mgloud/production_coloc_pipeline/

# GWAS
(
  awk -vOFS="\t" 'BEGIN{print "chr", "snp_pos", "alt", "ref", "beta", "se", "pvalue"}'
  gunzip -c /scratch/jhz22/BMI/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz | \
  awk -vOFS="\t" 'NR>1' | \
  cut -f1,2,4,5,7,8,9 | \
  sort -k1,1n -k2,2n
) > bmi.sumstats
bgzip -f bmi.sumstats
tabix -f -S 1 -s 1 -b 2 -e 2 bmi.sumstats.gz

# eQTL
(
  awk -vOFS="\t" 'BEGIN{print "chr", "snp_pos", "alt", "ref", "chisq", "pvalue"}'
  gunzip -c /scratch/jhz22/INF/METAL/VEGF.A-1.tbl.gz | \
  awk 'NR>1' | \
  awk -vOFS="\t" '(NR>1 && $11!=0) {print $1,$2,$4,$5,($10/$11)^2,$12}' | \
  sort -k1,1n -k2,2n
) > VEGF.A.sumstats
bgzip -f VEGF.A.sumstats
tabix -f -S 1 -s 1 -b 2 -e 2 VEGF.A.sumstats.gz

# coloc
export p3v5=/scratch/public_databases/1000_Genomes_phase3v5a/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
python scripts/dispatch.py -g bmi.sumstats.gz -e VEGF.A.sumstats.gz -v ${p3v5} -N 2504 -c 6 -s 43737921
