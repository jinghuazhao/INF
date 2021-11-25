#!/bin/bash
# 11-11-2019 JHZ

function jma_cojo()
# IL.10R1 data: snpid, rsid, geno, phenocovar
{
  cut -f2 ${1}/${2} | awk 'NR>1' | sort -k1,1 > snpid
  awk /chr2:/ work/INTERVAL.rsid | sort -k1,1 | join - snpid > snpid-rsid
  cut -d' ' -f2 snpid-rsid > rsid
  export s=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/o5000-inf1-outlier_in-r2.sample
  qctool  -filetype bgen -g INTERVAL/per_chr/interval.imputed.olink.chr_2.bgen \
          -ofiletype dosage -og IL18R1.gen.gz -s ${s} \
          -incl-rsids rsid
  zcat IL18R1.gen.gz | \
  awk '
  {
    for (i=1;i<=NF;i++) a[NR,i]=$i
    if(big <= NF) big=NF
  }
  END {
     for(i=1;i<=big;i++)
     {
       for(j=1;j<=NR;j++) printf OFS a[j,i]; printf "\n"
     }
  }' | \
  awk 'NR != 1 && NR !=3 && NR != 4 && NR != 5 && NR != 6' > geno
  awk 'NR != 2' $s > phenocovar
  (
  R --no-save -q <<\ \ END
    g <- read.table("geno",as.is=TRUE,header=TRUE)
    pc <-read.table("phenocovar",as.is=TRUE,header=TRUE)
    gpc <- merge(g,pc,by.x="SNPID",by.y="ID_1")
    s <- read.table("rsid",as.is=TRUE)
    require(stringr)
    nog <- c("age", "sexPulse", "season", "plate", "bleed_to_process_time", paste0("PC",1:10), t(s))
    rhs <- str_c(nog, collapse="+")
    m1 <- with(gpc, lm(as.formula(paste0("IL.18R1___Q13478~",rhs))))
    summary(m1)
  END
  ) > IL.18R1.log
}

function INF1()
# INF1: previously 102810080
{
  jma_cojo work IL.18R1-chr2:102992675_C_T.jma.cojo
}

function INTERVAL()
# INTERVAL: solely the variant
{
  jma_cojo INTERVAL IL.18R1-chr2:102992675_C_T.jma.cojo
}

INF1

# rm snpid snpid-rsid rsid phenocovar geno IL18R1.gen.gz

function cond_assoc_test()
## SNPTEST v2.5.2 conditioning on independent variants from GCTA
{
  qctool  -filetype bgen -g INTERVAL/per_chr/interval.imputed.olink.chr_2.bgen \
          -ofiletype gen -og IL18R1.gen.gz -incl-range 101992675-103992675 -assume-chromosome 2 -s ${s}
  snptest_v2.5.2 -data IL18R1.gen.gz ${s} -pheno IL.18R1___Q13478 -full_parameter_estimates \
                 -condition_on $(cat rsid) \
                 -frequentist 1 -method em -use_raw_phenotypes -use_raw_covariates -cov_all -o IL.18R1.out
}

# assoc_test
