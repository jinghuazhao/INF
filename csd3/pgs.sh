# 12-3-2020 JHZ

export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
qctool -g ${UKB}/ukb_imp_chr#_v3.bgen -s $UKB/ukb_BP_imp_v3.sample \
       -ofiletype bgen_v1.1 -og work/crp.score.bgen -incl-snpids work/INF1.merge.ukbsnpid

grep -f work/INF1.merge.ukbsnpid -w work/crp.ukb | \
awk '
{
 # SNP=$2, A1=$5, A2=$6, minor_allele=$7, minor_AF=$8, P=$16
   if($5==$7) {EA=$5; score=$13} else {EA=$6; score=-$13}
   print $1,$2,EA,score,$16
}' | join work/INTERVAL.rsid - > work/crp.SNP.EA.beta.pvalue

export p="0.001 0.05 0.1 0.2 0.3 0.4 0.5"
awk 'BEGIN{split(ENVIRON["p"],a);for(i=1;i<=length(a);i++) print a[i],0,a[i]}' > work/crp.range_list

plink --bgen work/crp.score.bgen --sample $UKB/ukb_BP_imp_v3.sample \
      --out work/crp.score \
      --pheno work/crp.cvd --pheno-name cv \
      --q-score-range work/crp.range_list work/crp.SNP.EA.beta.pvalue 2 6 \
      --score work/crp.SNP.EA.beta.pvalue 2 4 5

R --no-save <<END
   p <- Sys.getenv("p")
   sink('work/crp-cv.out')
   for(pval in as.numeric(unlist(strsplit(p," "))))
   {
     profile <- paste0("work/crp.score.",pval,".profile")
     d <- read.table(profile,as.is=TRUE,header=TRUE)
     s <- subset(d,PHENO!=-999 & PHENO!=-9)
     cov <- read.delim("work/crp.cov",as.is=TRUE)
     scov <- merge(s,cov,by="FID")
     f <- paste0("PHENO~SCORE+sex+ages+",paste0("PC",1:50,collapse="+"))
     g <-glm(as.formula(f),family="binomial",data=scov)
     print(f)
     r <- summary(g)
     print(profile)
     print(r)
   }
   sink()
END

(
  gunzip -c ${INF}/METAL/*gz | \
  head -1 | \
  awk -v OFS="\t" '{$1="prot" OFS $1};1'
  cut -f5,6 --output-delimiter=' ' work/INF1.merge | \
  sed '1d' | \
  parallel -C' ' 'zgrep -H -w {2} ${INF}/METAL/{1}-1.tbl.gz | sed "s|${INF}/METAL/||g;s/-1.tbl.gz:/\t/g"'
) > work/INF1.merge.tbl

awk 'a[$4]++==0' work/INF1.merge.tbl > work/INF1.merge.1st

R --no-save -q <<END
  inf <- read.delim("work/INF1.merge.1st",as.is=TRUE)
  inf <- within(inf, {Allele1=toupper(Allele1);Allele2=toupper(Allele2)})
  ukb <- read.table("work/crp.SNP.EA.beta.pvalue",as.is=TRUE,col.names=c("MarkerName","rsid","ukbid","EA","beta","p"))
  vars <- c("prot","MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "log.P.", "Direction","N")
  inf_ukb <- merge(inf[vars],ukb,by="MarkerName")
  swap <- with(inf_ukb,Allele1!=EA)
  inf_ukb[swap,"beta"] <- -inf_ukb[swap,"beta"]
  with(inf_ukb,{
    plot(beta,Effect,cex=0.4)
    r <- lm(beta~Effect)
    summary(r)
  })
END
