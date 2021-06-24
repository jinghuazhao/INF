# 10-5-2020 JHZ

# irnt for CRP
export suffix=irnt
# Base data, CRP in mg/L
export suffix=raw
export SUMSTATS=ukb/30710_${suffix}.gwas.imputed_v3.both_sexes.tsv.bgz
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
) > work/crp.${suffix}

# Target data, UKB
export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
join <(awk 'NR>1 {print $1}' work/crp.${suffix} | sort -k1,1) \
     <(awk 'NR > 1 {print $1,$2}' work/INF1.merge.ukbsnp | sort -k1,1) | \
cut -d' ' -f2 > work/INF1.merge.ukbsnpid
(
  awk -v OFS="\t" 'NR == 1 {$1=$1 "\t" "id"; print}' work/crp.${suffix}
  join -t$'\t' <(awk -v OFS="\t" 'NR > 1 {print $1,$2}' work/INF1.merge.ukbsnp | sort -k1,1) \
               <(awk -v OFS="\t" 'NR > 1' work/crp.${suffix})
) > work/crp-${suffix}.ukb
sed 's/ID_1/FID/g;s/ID_2/IID/g;2d' ${UKB}/ukb_BP_imp_v3.sample | cut -d' ' -f1,2,4 > work/crp.sample

module load ceuadmin/stata
stata <<END
  insheet using work/crp.sample, case clear delim(" ")
  sort FID
  save work/crp, replace
  gzuse ukb/analysis, clear
  gen chd=ep1_chd
  gen cv=ep1_cv
  destring idno, gen(FID)
  merge 1:1 FID using work/crp, gen(m1)
  drop if FID==. | FID <0 | IID==. | IID<0
  replace chd=. if ep1_chd==6
  replace cv=. if ep1_cv==6
  merge 1:1 FID using work/ukb.pca.dta, gen(m2)
  keep FID IID chd cv sex ages PC1-PC50
  save work/ukb, replace
  mvencode _all, mv(-999)
  outsheet FID IID chd cv using work/crp.cvd, nolabel noquote replace
  outsheet FID IID sex ages PC1-PC50 using work/crp.cov, nolabel noquote replace
  outsheet FID IID using work/crp-chd.excl if chd==-999 | sex==-999 | ages==-999, nolabel noquote replace
  outsheet FID IID using work/crp-cv.excl if cv==-999 | sex==-999 | ages==-999, nolabel noquote replace
END

# PRSice analysis
## .all.score
PRSice --base work/crp-${suffix}.ukb --snp id --chr chr --bp pos --A1 A1 --A2 A2 --beta beta --pvalue pval \
       --target ${UKB}/ukb_imp_chr#_v3 --type bgen --pheno work/crp.sample \
       --extract work/INF1.merge.ukbsnpid --model add --no-regress --score avg \
       --out work/crp
## .best
for pheno in chd cv
do
  PRSice --base work/crp-${suffix}.ukb --snp id --chr chr --bp pos --A1 A1 --A2 A2 --beta beta --pvalue pval \
         --target ${UKB}/ukb_imp_chr#_v3 --type bgen --pheno work/crp.cvd --pheno-col ${pheno} \
                                                     --cov work/crp.cov --remove work/crp-${pheno}.excl --binary-target T \
         --extract work/INF1.merge.snp --model add --score avg \
         --out work/crp-${pheno}
done
stata -b -q <<END
foreach v in "chd" "cv" {
  insheet using work/crp-`v'.best, case clear delim(" ")
  sort FID
  merge 1:1 FID using work/ukb
  logit `v' sex ages PRS
  stset ages, failure(`v')
  stcox sex PRS if `v'!=.
}
END
