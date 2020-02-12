# 12-2-2020 JHZ

join work/INTERVAL.rsid <( \
sed '1d' work/crp.ukb | awk '
{
 # SNP=$2, A1=$5, A2=$6, minor_allele=$7, minor_AF=$8, P=$16
   if($5==$7) {EA=$5; score=$13} else {EA=$6; score=-$13}
   print $1,EA,score,$16
}') | \
cut -d ' ' -f2-5 > work/crp.SNP.EA.beta.pvalue
export p="0.001 0.05 0.1  0.2  0.3  0.4  0.5"
awk 'BEGIN{split(ENVIRON["p"],a);for(i=1;i<=length(a);i++) print a[i],0,a[i]}' > work/crp.range_list

# extract SNPs?
export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K

function extract()
{
qctool \
    -g ${UKB}/ukb_imp_chr#_v3.bgen -s $UKB/ukb_BP_imp_v3.sample \
    -ofiletype bgen_v1.1 -og work/crp.bgen \
    -incl-rsids work/INF1.merge.snp
}

plink \
    --bgen work/crp.bgen --sample $UKB/ukb_BP_imp_v3.sample \
    --score work/crp.SNP.EA.beta.pvalue 1 2 3 double-dosage \
    --q-score-range work/crp.range_list work/crp.SNP.EA.beta.pvalue 1 4 \
    --pheno work/crp.cvd --pheno-name cv \
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
