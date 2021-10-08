# 8-10-2021 JHZ

function setup()
# irnt, raw for base data, CRP in mg/L
{
  export suffix=$1
  export SUMSTATS=${INF}/ukb/30710_${suffix}.gwas.imputed_v3.both_sexes.tsv.bgz
}

setup raw
(
  gunzip -c ${SUMSTATS} | \
  head -1 | \
  awk -vOFS="\t" '{$1="snpid\tchr\tpos\tA1\tA2"};1'
  gunzip -c ${SUMSTATS} | \
  awk '{
     OFS="\t"
     if (NR>1)
     {
       split($1,a,":")
       CHR=a[1]
       POS=a[2]
       a1=a[4]
       a2=a[3]
       if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
       else snpid="chr" CHR ":" POS "_" a1 "_" a2
       $1=snpid "\t" a[1] "\t" a[2] "\t" a[4] "\t" a[3]
     }
     print
  }' | \
  sort -k1,1 | \
  join -t$'\t' - <(sed '1d' work/INF1.merge | cut -f6 | sort -k1,1 | uniq)
) > ${INF}/crp/crp.${suffix}

# Target data, UKB
export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
join <(awk 'NR>1 {print $1}' ${INF}/crp/crp.${suffix} | sort -k1,1) \
     <(awk 'NR > 1 {print $1,$2}' ${INF}/work/INF1.merge.ukbsnp | sort -k1,1) | \
cut -d' ' -f2 > ${INF}/crp/INF1.merge.ukbsnpid
(
  awk -v OFS="\t" 'NR == 1 {$1=$1 "\t" "id"; print}' ${INF}/crp/crp.${suffix}
  join -t$'\t' <(awk -v OFS="\t" 'NR > 1 {print $1,$2}' ${INF}/work/INF1.merge.ukbsnp | sort -k1,1) \
               <(awk -v OFS="\t" 'NR > 1' ${INF}/crp/crp.${suffix})
) > ${INF}/crp/crp-${suffix}.ukb
sed 's/ID_1/FID/g;s/ID_2/IID/g;2d' ${UKB}/ukb_BP_imp_v3.sample | cut -d' ' -f1,2,4 > ${INF}/crp/crp.sample

module load ceuadmin/stata
stata <<END
  local INF : env INF
  insheet using `INF'/crp/crp.sample, case clear delim(" ")
  sort FID
  save `INF'/crp/crp, replace
  gzuse `INF'/ukb/analysis, clear
  gen chd=ep1_chd
  gen cv=ep1_cv
  destring idno, gen(FID)
  merge 1:1 FID using `INF'/crp/crp, gen(m1)
  drop if FID==. | FID <0 | IID==. | IID<0
  replace chd=. if ep1_chd==6
  replace cv=. if ep1_cv==6
  merge 1:1 FID using `INF'/crp/ukb.pca.dta, gen(m2)
  keep FID IID chd cv sex ages PC1-PC50
  save `INF'/crp/ukb, replace
  mvencode _all, mv(-999)
  outsheet FID IID chd cv using `INF'/crp/crp.cvd, nolabel noquote replace
  outsheet FID IID sex ages PC1-PC50 using `INF'/crp/crp.cov, nolabel noquote replace
  outsheet FID IID using `INF'/crp/crp-chd.excl if chd==-999 | sex==-999 | ages==-999, nolabel noquote replace
  outsheet FID IID using `INF'/crp/crp-cv.excl if cv==-999 | sex==-999 | ages==-999, nolabel noquote replace
END

# PRSice analysis
## .all.score
PRSice --base ${INF}/crp/crp-${suffix}.ukb --snp id --chr chr --bp pos --A1 A1 --A2 A2 --beta beta --pvalue pval \
       --target ${UKB}/ukb_imp_chr#_v3 --type bgen --pheno ${INF}/crp/crp.sample \
       --extract ${INF}/crp/INF1.merge.ukbsnpid --model add --no-regress --score avg \
       --out ${INF}/crp/crp
## .best
for pheno in chd cv
do
  PRSice --base ${INF}/crp/crp-${suffix}.ukb --snp id --chr chr --bp pos --A1 A1 --A2 A2 --beta beta --pvalue pval \
         --target ${UKB}/ukb_imp_chr#_v3 --type bgen --pheno ${INF}/crp/crp.cvd --pheno-col ${pheno} \
                                                     --cov ${INF}/crp/crp.cov --remove `INF'/crp/crp-${pheno}.excl --binary-target T \
         --extract ${INF}/crp/INF1.merge.snp --model add --score avg \
         --out ${INF}/crp/crp-${pheno}
done
stata -b -q <<END
foreach v in "chd" "cv" {
  local INF : env INF
  insheet using `INF'/crp/crp-`v'.best, case clear delim(" ")
  sort FID
  merge 1:1 FID using `INF'/crp/ukb
  logit `v' sex ages PRS
  stset ages, failure(`v')
  stcox sex PRS if `v'!=.
}
END

setup irnt
(
  gunzip -c ${SUMSTATS} | \
  head -1 | \
  awk -vOFS="\t" '{$1="snpid\tchr\tpos\tA1\tA2\tvariant"};1'
  gunzip -c ${SUMSTATS} | \
  awk '{
     OFS="\t"
     if (NR>1)
     {
       split($1,a,":")
       CHR=a[1]
       POS=a[2]
       a1=a[4]
       a2=a[3]
       if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
       else snpid="chr" CHR ":" POS "_" a1 "_" a2
       $1=snpid "\t" a[1] "\t" a[2] "\t" a[4] "\t" a[3] "\t" $1
     }
     print
  }' | \
  sort -k1,1 | \
  join -t$'\t' - <(sed '1d' work/INF1.merge | cut -f6 | sort -k1,1 | uniq)
) > ${INF}/crp/crp-ext.${suffix}

(
  join -14 -21 <(head -1 ${INF}/crp/crp-ext.${suffix} | sort -k6,6 | cut -f1,4-6,13-16) \
               <(gunzip -c ${INF}/ukb/variants.tsv.bgz | \
                 cut -f1-10 | \
                 head -1 | \
                 sort -k1,1)
  join -14 -21 <(sed '1d' ${INF}/crp/crp-ext.${suffix} | sort -k6,6 | cut -f1,4-6,13-16) \
               <(gunzip -c ${INF}/ukb/variants.tsv.bgz | \
                 cut -f1-10 | \
                 sed '1d' | \
                 sort -k1,1)
) > ${INF}/crp/crp-ext-ann.${suffix}

module load gcc/6
R --no-save -q <<END
  INF <- Sys.getenv("INF")
  suffix <- Sys.getenv("suffix")
  f <- file.path(INF,"crp",paste0("crp-ext-ann.",suffix))
  m <- read.table(f,header=TRUE,as.is=TRUE)
  options(width=200)
  library(dplyr)
  t <- m %>% group_by(snpid) %>% slice(which.min(pval))
  m <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
       select(MarkerName,prot,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans) %>%
       mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
       group_by(MarkerName) %>%
       slice(which.min(log.P.)) %>%
       left_join(m, by=c('MarkerName'='snpid')) %>%
       mutate(beta=if_else(Allele1==alt,beta,-beta)) %>%
       select(prot,MarkerName,rsid,Effect,StdErr,log.P.,beta,se,pval,cis.trans,consequence,consequence_category,info)
  write.table(m,file=file.path(INF,"crp","crp-ext-ann.txt"),row.names=FALSE,quote=FALSE,sep="\t")
  with(filter(m,!is.na(beta)),cor(Effect,beta,use="everything"))
  with(filter(m,!is.na(beta)&cis.trans=="cis"),cor(Effect,beta,use="everything"))
  filter(m,consequence_category=="ptv")
  with(m,plot(Effect,beta))
END

