# 22-1-2019 JHZ

module load bcftools/1.9
module load plink2/1.90beta5.4
module load bolt-lmm/2.3.2

function info()
{
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
  parallel -j1 -C' ' 'awk -vchr={} "NR>1 && \$7>=0.4 {print chr \":\" \$3}" chr{}.info > chr{}.id
# only ~1M variants left with info>=0.4, so preferably skipped
  plink --vcf chr{}.vcf.gz --extract chr{}.id --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out chr{}
  plink --vcf chr{}.vcf.gz --extract chr{}.prune.in --make-bed --out chr{}'
}

function ln()
{
  ## SNP filtering for .vcf.gz
  export vcf="/data/jinhua/data/KORA/Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220"

  seq 22 | \
  parallel -j3 -C' ' 'ln -sf Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220/chr{}_N3775_1000GPhase3.dose.vcf.gz chr{}.vcf.gz'
}

function id()
{
  bcftools query -l chr22.vcf.gz | \
  sort > genotype.id
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
  join genotype.id - | \
  cut -d' ' -f1 > protein.id
  join -v1 genotype.id protein.id | \
  awk -vOFS="\t" '{print $1,$1}' > remove.id
}

function snp()
{
  seq 22 | \
  parallel -j1 -C' ' '
# bcftools annotate --set-id "chr%CHROM\:%POS\_%REF\_%ALT" chr{}.vcf.gz -O z -o KORA{}.vcf.gz
  plink --vcf KORA{}.vcf.gz --list-duplicate-vars --out chr{}
  awk "NR>1{split(\$NF,dupids,\" \");print dupids[1]}" chr{}.dupvar > chr{}.dupid
  bcftools query -i "MAF>0.01 && R2>=0.4" -f"%ID\n" KORA{}.vcf.gz | \
  join -v2 chr{}.dupid - > chr{}.mafr2
  plink --vcf KORA{}.vcf.gz --extract chr{}.mafr2 --remove remove.id --make-bed --out nodup{}
  awk -vOFS="\t" "
  {
    CHR=\$1
    POS=\$4
    a1=\$5
    a2=\$6
    if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
    else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
    print snpid, \$2
  }" nodup{}.bim > nodup{}.snpid
  plink --bfile nodup{} --update-name nodup{}.snpid 1 2 --make-bed --out KORA{}
  '
  seq 22 | \
  awk -vp=KORA '{print p NR}' > merge-list
  plink --merge-list merge-list --make-bed --out KORA
  plink --bfile KORA --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out KORA
  plink --bfile KORA --extract KORA.prune.in --make-bed --out KORA.prune'
  seq 22 | \
  parallel -j3 -C' ' 'bcftools convert --samples-file protein.id KORA{}.vcf.gz -g protein{}'
  parallel -j1 'echo {} KORA{}.gen.gz' > KORA.list
  awk -vOFS="\t" '{print $1, $1}' protein.id > KORA.id
}

function assoc()
{
## association analysis
# https://data.broadinstitute.org/alkesgroup/BOLT-LMM/#x1-220005.1.2
  seq 22 | \
  parallel -j2 -C' ' '
  bolt \
  --bfile KORA.prune \
  --impute2FileList=KORA.list \
  --impute2FidIidFile=KORA.id \
  --LDscoresUseChip \
  --phenoFile=phenocovar.txt --phenoCol UH_O_{} \
  --covarFile=phenocovar.txt --covarCol sex --covarCol age \
  --remove remove.id \
  --lmm --statsFileImpute2Snps={}-snp --statsFile={}-stats 2>&1 | tee {}.log' ::: OPG TNFSF14
}

cd KORA
assoc
