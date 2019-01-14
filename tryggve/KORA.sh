# 14-1-2019 JHZ

# concatenate .info and filter .info>=0.4
export info=/data/jinhua/data/KORA/impute_out_info/info

ls $info | \
sed 's/KORAS4F4_AffyAxiom_N3788_chr//g;s/_imputation.impute2_info\*//g;s/_/ /g' | \
sort -t' ' -k1,1g -k2,2g | \
parallel -j1 --env info -C' '  '
(
  head -1 $info/KORAS4F4_AffyAxiom_N3788_chr{1}_1_imputation.impute2_info
  awk "NR>1" $info/KORAS4F4_AffyAxiom_N3788_chr{1}_{2}_imputation.impute2_info
) > chr{1}.info
'
seq 22 | \
parallel -j1 -C' ' 'awk "NR>1 && \$7>=0.4 {print \$2}" chr{}.info > chr{}.rs_id'

# SNP pruning
export vcf="/data/jinhua/data/KORA/Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220"

seq 22 | \
parallel -j3 -C' ' 'ln -sf Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220/chr{}_N3775_1000GPhase3.dose.vcf.gz chr{}.vcf.gz'

module load plink2/1.90beta5.4
seq 22 | \
parallel -j3 -C' ' '
plink --vcf chr{}.vcf.gz --extract chr{}.rs_id --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out chr{}
plink --vcf chr{}.vcf.gz --extract chr{}.prune.in --make-bed --out chr{}'
