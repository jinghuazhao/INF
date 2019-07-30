# 30-7-2019 JHZ

# IL.10R1 data

function jma_cojo()
{
  cut -f2 ${1}/IL.18R1-chr2:102810080_A_G.jma.cojo | awk 'NR>1' | sort -k1,1 > l
  awk /chr2:/ work/INTERVAL.rsid | sort -k1,1 | join - l > ll
  cut -d' ' -f2 ll > lll
  export s=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/o5000-inf1-outlier_in-r2.sample
  qctool  -filetype bgen -g INTERVAL/per_chr/interval.imputed.olink.chr_2.bgen \
          -ofiletype dosage -og IL18R1.gen.gz -s ${s} \
          -incl-rsids lll
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

  R --no-save -q <<\ \ END
    g <- read.table("geno",as.is=TRUE,header=TRUE)
    pc <-read.table("phenocovar",as.is=TRUE,header=TRUE)
    gpc <- merge(g,pc,by.x="SNPID",by.y="ID_1")
    s <- read.table("lll",as.is=TRUE)
    require(stringr)
    nog <- c("age", "sexPulse", "season", "plate", "bleed_to_process_time", paste0("PC",1:20), t(s))
    rhs <- str_c(nog, collapse="+")
    m1 <- with(gpc, lm(as.formula(paste0("IL.18R1___Q13478~",rhs))))
    summary(m1)
  END
}

function assoc_test()
## SNPTEST v2.5.2 for CKDGen-type finemap analysis
{
  qctool  -filetype bgen -g INTERVAL/per_chr/interval.imputed.olink.chr_2.bgen \
          -ofiletype gen -og IL18R1.gen.gz -incl-range 101810080-103810080 -assume-chromosome 2 -s ${s}
  snptest_v2.5.2 -data IL18R1.gen.gz ${s} -pheno IL.18R1___Q13478 -full_parameter_estimates \
                 -condition_on $(cat lll) \
                 -frequentist 1 -method em -use_raw_phenotypes -use_raw_covariates -cov_all -o IL.18R1.out
}

jma_cojo work/test
jma_cojo work
assoc_test

# rm l ll lll phenocovar geno IL18R1.gen.gz
