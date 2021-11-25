# 7-10-2021 JHZ

export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
export p="0.001 0.05 0.1 0.2 0.3 0.4 0.5"
export suffix=raw

function ukb_bgen()
{
qctool -g ${UKB}/ukb_imp_chr#_v3.bgen -s ${UKB}/ukb_BP_imp_v3.sample \
       -ofiletype bgen_v1.1 -og ${INF}/crp/crp.score.bgen -incl-snpids ${INF}/crp/INF1.merge.ukbsnpid
}
awk 'BEGIN{split(ENVIRON["p"],a);for(i=1;i<=length(a);i++) print a[i],0,a[i]}' > ${INF}/crp/crp.range_list

## UKB

function qqman()
{
  R --no-save -q <<\ \ END
  INF <- Sys.getenv("INF")
  suffix <- Sys.getenv("suffix")
  require(qqman)
  gz <- file.path(INF,"ukb",paste0("30710_",suffix,".gwas.imputed_v3.both_sexes.tsv.bgz"))
  tbl <- read.delim(gz,as.is=TRUE)
  tbl <- within(tbl,{
     SNP <- variant
     SNP_split <- sapply(tbl[["variant"]],strsplit,":",simplify=TRUE)
     CHR <- as.integer(unlist(lapply(SNP_split,"[",1)))
     BP <- as.integer(unlist(lapply(SNP_split,"[",2)))
     P <- pval
  })
  tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
  qq <- paste0("crp-",suffix,"_qq.png");
  png(qq,width=12,height=10,units="in",pointsize=4,res=300)
  qq(with(tbl,P))
  dev.off()
  manhattan <- paste0("crp-",suffix,"_manhattan.png");
  png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
  manhattan(tbl,main="CRP",genomewideline=-log10(5e-8),suggestiveline=FALSE,ylim=c(0,50));
  dev.off();
  END
}

function UKB()
{
grep -f ${INF}/crp/INF1.merge.ukbsnpid -w ${INF}/crp/crp-${suffix}.ukb | \
awk '
{
 # SNPID=$1, id=$2, A1=$5, A2=$6, minor_allele=$7, minor_AF=$8, P=$16
   if($5==$7) {EA=$5; beta=$13} else {EA=$6; beta=-$13}
   print $1,$2,EA,beta,$16
}' | join ${INF}/work/INTERVAL.rsid - > ${INF}/crp/crp.${suffix}.EA.beta.pvalue

awk '{$1=ENVIRON["suffix"] "-" $1};1' ${INF}/crp/crp.range_list > ${INF}/crp/crp-${suffix}.range_list

plink --bgen ${INF}/crp/crp.score.bgen --sample ${UKB}/ukb_BP_imp_v3.sample \
      --out ${INF}/crp/crp.score \
      --pheno ${INF}/crp/crp.cvd --pheno-name cv \
      --q-score-range ${INF}/crp/crp-${suffix}.range_list ${INF}/crp/crp.${suffix}.EA.beta.pvalue 2 6 \
      --score ${INF}/crp/crp.SNP.EA.beta.pvalue 2 4 5

R --no-save <<END
   INF <- Sys.getenv("INF")
   p <- Sys.getenv("p")
   suffix <- Sys.getenv("suffix")
   sink(file.path(INF,"crp",paste0("crp-cv-",suffix,".out")))
   for(pval in as.numeric(unlist(strsplit(p," "))))
   {
     profile <- file.path(INF,"crp",paste0("crp.score.",suffix,"-",pval,".profile"))
     d <- read.table(profile,as.is=TRUE,header=TRUE)
     s <- subset(d,PHENO!=-999 & PHENO!=-9)
     cov <- read.delim(file.path(INF,"crp","crp.cov"),as.is=TRUE)
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
}

## INF

function INF()
{
cut -d' ' -f1-3 ${INF}/crp/crp.${suffix}.EA.beta.pvalue > ${INF}/crp/crp.${suffix}.id
cut -f4,5,11,13 ${INF}/crp/INF1.merge.1st | sed '1d' | sort -k1,1 | awk '{$2=toupper($2);$4=10^$4};1' | \
join ${INF}/crp/crp.${suffix}.id - > ${INF}/crp/crp.INF.EA.beta.pvalue

awk '{$1="INF-" $1};1' ${INF}/crp/crp.range_list > ${INF}/crp/crp-INF.range_list

plink --bgen ${INF}/crp/crp.score.bgen --sample ${UKB}/ukb_BP_imp_v3.sample \
      --out ${INF}/crp/crp.score \
      --pheno ${INF}/crp/crp.cvd --pheno-name cv \
      --q-score-range work/crp-INF.range_list ${INF}/crp/crp.INF.EA.beta.pvalue 2 6 \
      --score ${INF}/crp/crp.INF.EA.beta.pvalue 2 4 5

R --no-save <<END
   INF <- Sys.getenv("INF")
   p <- Sys.getenv("p")
   sink(paste0("work/crp-cv-INF.out"))
   for(pval in as.numeric(unlist(strsplit(p," "))))
   {
     profile <- paste0("work/crp.score.INF-",pval,".profile")
     d <- read.table(profile,as.is=TRUE,header=TRUE)
     s <- subset(d,PHENO!=-999 & PHENO!=-9)
     cov <- read.delim(file.path(INF,"crp","crp.cov"),as.is=TRUE)
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
}

# UKB_INF

function UKB_INF()
{
(
  gunzip -c ${INF}/METAL/*gz | \
  head -1 | \
  awk -v OFS="\t" '{$1="prot" OFS $1};1'
  cut -f5,6 --output-delimiter=' ' ${INF}/work/INF1.merge | \
  sed '1d' | \
  parallel -C' ' 'zgrep -H -w {2} ${INF}/METAL/{1}-1.tbl.gz | sed "s|${INF}/METAL/||g;s/-1.tbl.gz:/\t/g"'
) > ${INF}/crp/INF1.merge.tbl

awk 'a[$4]++==0' ${INF}/crp/INF1.merge.tbl > ${INF}/crp/INF1.merge.1st

R --no-save -q <<END
  INF <- Sys.getenv("INF")
  suffix <- Sys.getenv("suffix")
  inf <- read.delim(file.path(INF,"crp","INF1.merge.1st"),as.is=TRUE)
  inf <- within(inf, {Allele1=toupper(Allele1);Allele2=toupper(Allele2)})
  vars <- c("prot","MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "log.P.", "Direction","N")
  ukb <- read.table(file.path(INF,"crp",paste0("crp.",suffix,".EA.beta.pvalue")),as.is=TRUE,
                    col.names=c("MarkerName","rsid","id","EA","beta","p"))
  inf_ukb <- merge(inf[vars],ukb,by="MarkerName")
  swap <- with(inf_ukb,Allele1!=EA)
  inf_ukb[swap,"beta"] <- -inf_ukb[swap,"beta"]
  inf_ukb[c("Effect","beta")]
  summary(inf_ukb)
  with(inf_ukb,{
    pdf(file.path(INF,"crp",paste0("INF1.ukb-",suffix,".pdf")))
    plot(Effect,beta,xlab="INF",ylab="UKB",cex=0.4)
    title("Effect sizes in INF vs ukb")
    print(cor(Effect,beta))
    r <- lm(beta~Effect)
    summary(r)
    dev.off()
  })
END
}

# ukb_bgen
# UKB
  INF
  UKB_INF
