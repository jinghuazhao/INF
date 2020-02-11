# 11-2-2020 JHZ

export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
cut -f2,16 work/crp.ukb > work/crp.SNP.pvalue
(
  echo "0.001 0 0.001"
  echo "0.05 0 0.05"
  echo "0.1 0 0.1"
  echo "0.2 0 0.2"
  echo "0.3 0 0.3"
  echo "0.4 0 0.4"
  echo "0.5 0 0.5"
) > work/crp.range_list

awk '
{
  # SNP=$2, A1=$5, A2=$6, $minor_allele=$7, minor_AF=$8, P=$16
  if (NR==1) print "SNP", "EAF", "P";
  else {
    EA=$5
    if($5==$7) MAF <- $8;
    else {EA=$6; MAF=1-$7}
    print $2,EA,$16
  }
}' FS="\t" OFS="\t" work/crp.ukb > work/crp.SNP.MAF.pvalue
cut -f2 work/crp.ukb | sed '1d' > work/crp.SNP

# extract SNPs?
qctool \
    -g ${UKB}/ukb_imp_chr#_v3.bgen -s $UKB/ukb_BP_imp_v3.sample \
    -og work/crp.bgen -ofiletype bgen_v1.1 \
    -incl-rsids work/INF1.merge.snp

plink \
    --bgen work/crp.bgen --sample $UKB/ukb_BP_imp_v3.sample \
    --score work/crp.SNP.MAF.pvalue 1 2 3 header \
    --q-score-range work/crp.range_list work/crp.SNP.pvalue \
    --covar work/crp.cov --pheno work/crp.cvd --pheno-name cv \
    --out work/crp.score

R --no-save <<END
  for(p in c(0.001,0.05,0.1,0.2,0.3,0.4,0.5))
  {
    profile <- paste0("work/crp.score.",p,".profile")
    print(profile)
    d <- read.table(profile,as.is=TRUE,header=TRUE)
    s <- subset(d,PHENO!=-999 & PHENO!=-9)
    g <-glm(PHENO~SCORE,family="binomial",data=s)
    r <- summary(g)
    print(r)
  }
END
