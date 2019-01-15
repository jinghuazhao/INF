# 15-1-2019 JHZ

## concatenate .info and filter .info>=0.4
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
parallel -j1 -C' ' 'awk -vchr={} "NR>1 && \$7>=0.4 {print chr \":\" \$3}" chr{}.info > chr{}.id'

## SNP filtering for .vcf.gz
export vcf="/data/jinhua/data/KORA/Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220"

seq 22 | \
parallel -j3 -C' ' 'ln -sf Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220/chr{}_N3775_1000GPhase3.dose.vcf.gz chr{}.vcf.gz'

module load bcftools/1.9
# bcftools query -f"%CHROM\t%POS\t%REF\t%ALT\t%ID\n" chr22.vcf.gz

## covariates and proteins
bcftools query -l chr22.vcf.gz | \
sort > genotyped.id
export llod=/data/jampet/KORA/kora.below.llod.normalised.prot.txt
awk -vOFS="\t" '{
  $1=$1 OFS $1
  if($3=="M") $3=0; else if($3=="F") $3=1
  print
}' $llod | \
awk -vOFS="\t" '{if(NR==1) {$1="FID"; $2="IID"}};1' > phenocovar.txt
cut -f1-2 phenocovar.txt | \
awk 'NR>1' | \
sort -k1,1 | \
join genotyped.id - | \
cut -d' ' -f1 > protein.id
join -v2 genotype.id protein.id | \
awk -vOFS="\t" '{print $1,$1}' > remove.id

module load plink2/1.90beta5.4

## good quality SNPs
seq 22 | \
parallel -j3 -C' ' '
plink --vcf chr{}.vcf.gz --list-duplicate-vars require-same-ref --out chr{}
awk "NR>1{split(\$NF,dupids,\" \");print dupids[1]}" chr{}.dupvar > chr{}.dupid
plink --vcf chr{}.vcf.gz --exclude chr{}.dupid --remove remove.id --make-bed --out chr{} --threads 1
'
(
  seq 22 | \
  parallel -j3 -C' ' 'bcftools query -i "MAF>0.01 && R2>=0.4" -f"%ID\n" chr{}.vcf.gz > chr{}.id'
) > MAFR2.id
seq 22 | \
awk -vp=chr '{print p $1}' > merge-list
plink --merge-list merge-list --extract MAFR2.id --make-bed --out KORA

function prune()
# for .info files only ~1M variants left with info>=0.4, so preferably skipped
{
  seq 22 | \
  parallel -j3 -C' ' '
  plink --vcf chr{}.vcf.gz --extract chr{}.id --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out chr{}
  plink --vcf chr{}.vcf.gz --extract chr{}.prune.in --make-bed --out chr{}'
}

## association analysis
module load bolt-lmm/2.3.2

# https://data.broadinstitute.org/alkesgroup/BOLT-LMM/#x1-220005.1.2

parallel -j2 -C' ' '
  --bfile KORA \
  --phenoFile=phenocovar.txt --phenoCol UH_O_{} \
  --covarFile=phenocovar.txt --covarCol sex age \
  --remove remove.id \
  --lmm --statsFileImpute2Snps={}-snp --statsFile={}-stats > 2>&1 | tee {}.log' ::: OPG TNFSF14
'
