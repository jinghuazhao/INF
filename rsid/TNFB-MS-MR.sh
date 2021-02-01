#!/usr/bin/bash

export prot=TNFB

function pgz()
{
# all significant SNPs
  (
    zcat ${INF}/METAL/${prot}-1.tbl.gz | head -1
    zcat ${INF}/METAL/${prot}-1.tbl.gz | \
    awk 'NR>1 && length($4)==1 && length($5)==1 && $12<=-5' | \
    sort -k1,1n -k2,2n | \
    gzip -f > ${INF}/MS/${prot}.p.gz
# HLA
    zcat ${INF}/MS/${prot}.p.gz | \
    awk 'NR>1 && !($1 == 6 && $2 >= 25392021 && $2 < 33392022)'
    zcat {INF}/MS/${prot}.p.gz | \
    awk '$1 == 6 && $2 >= 25392021 && $2 < 33392022' | \
    sort -k12,12g | \
    awk 'NR==1'
  ) > ${INF}/MS/${prot}.p
  export lines=$(expr $(wc -l ${INF}/MS/${prot}.p | cut -d' ' -f1) - 1)
  if [ $lines -eq 0 ]; then rm ${INF}/MS/${prot}.p; fi
  (
    awk -vprot=${prot} -vOFS="\t" 'BEGIN{print prot, "rsid", "chr", "pos", "beta", "se", "snpid", "A1", "A2", "EAF", "P", "N"}'
    awk -vFS="\t" '{split($3,a,"_"); print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,-$12,$18}' ${prot}.p | \
    sort -k1,1 | \
    join -12 -21 ${INF}/work/snp_pos - | \
    awk 'a[$7]++==0' | \
    awk -vprot=${prot} -vOFS="\t" '{print prot, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | \
    sort -k3,3n -k4,4n
  ) > ${INF}/MS/${prot}-pQTL.dat
  R --no-save -q <<\ \ END
    options(echo=FALSE,width=200)
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rsid <- Sys.getenv("rsid")
    pQTL <- file.path(INF,"MS",paste0(prot,"-pQTL.dat"))
    library(TwoSampleMR)
    x <- read_exposure_data(pQTL,
         clump = TRUE,
         sep = "\t",
         phenotype_col = "TNFB",
         snp_col = "rsid",
         beta_col = "beta",
         se_col = "se",
         eaf_col = "EAF",
         effect_allele_col = "A1",
         other_allele_col = "A2",
         pval_col = "P",
         samplesize_col = "N",
         id_col = "rsid",
         log_pval = TRUE)
    for(id in c("ieu-b-18","ieu-a-1025","ukb-b-17670","finn-a-G6_MS"))
    {
      cat("--",pQTL,"-",id,"--\n")
      y <- extract_outcome_data(with(x,SNP), id, proxies = TRUE, rsq = 0.8)
      xy <- mr(harmonise_data(x, y))
      print(xy)
    }
  END
}

pgz > ${INF}/MS/${prot}-MS-MR.log

function pqtl_qtl_rsid()
{
  export start=$(expr ${pos} - ${M})
  export end=$(expr ${pos} + ${M})
  if [ ${get_data} == "yes" ]; then
    (
      awk -vprot=${prot} -vOFS="\t" 'BEGIN{print prot, "rsid", "chr", "pos", "beta", "se", "snpid", "A1", "A2", "EAF", "P", "N"}'
      gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
      awk -vFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} '(NR>1 && $1 == chr && $2 >= start && $2 <= end) {
          split($3,a,"_"); print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,-$12,$18
      }' | \
      sort -k1,1 | \
      join -12 -21 ${INF}/work/snp_pos - | \
      awk 'a[$7]++==0' | \
      awk -vprot=${prot} -vOFS="\t" '{print prot, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | \
      sort -k3,3n -k4,4n
    ) > ${INF}/work/${prot}-pQTL-${rsid}.dat
  fi
  R --no-save -q <<\ \ END
    options(echo=FALSE,width=200)
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rsid <- Sys.getenv("rsid")
    pQTL <- file.path(INF,"work",paste0(prot,"-pQTL-",rsid,".dat"))
    library(TwoSampleMR)
    x <- subset(read_exposure_data(pQTL,
                clump = TRUE,
                sep = "\t",
                phenotype_col = "TNFB",
                snp_col = "rsid",
                beta_col = "beta",
                se_col = "se",
                eaf_col = "EAF",
                effect_allele_col = "A1",
                other_allele_col = "A2",
                pval_col = "P",
                samplesize_col = "N",
                id_col = "rsid",
                log_pval = TRUE), pval<1e-5)
    for(id in c("ieu-b-18","ieu-a-1025","ukb-b-17670","finn-a-G6_MS"))
    {
      cat("##",pQTL,"-",id,"##\n")
      y <- extract_outcome_data(with(x,SNP), id, proxies = TRUE, rsq = 0.8)
      xy <- mr(harmonise_data(x, y))
      print(xy)
    }
  END
  rm work/${prot}-*-${rsid}.dat
}

function pqtl()
# pQTLs
{
  (
     awk -vOFS="\t" 'BEGIN{print "Phenotype", "SNP", "chr", "pos", "beta", "se", "snpid", "effect_allele", "other_allele", "eaf", "pval", "N"}'
     zgrep -e chr12:6514963_A_C -e chr12:6440009_C_T -e chr6:31540757_A_C -e chr6:31073047_A_G ${INF}/METAL/${prot}-1.tbl.gz | \
     awk -vFS="\t" '{split($3,a,"_"); print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,10^$12,$18}' | \
     sort -k1,1 | \
     join -12 -21 ${INF}/work/snp_pos - | \
     awk 'a[$7]++==0' | \
     awk -vprot=${prot} -vOFS="\t" '{print prot, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | \
     sort -k3,3n -k4,4n
  ) > ${INF}/work/${prot}-pQTL.dat
  R --no-save -q <<\ \ END
    library(pQTLtools)
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    pqtl <- file.path(INF,"work",paste0(prot,"-pQTL.dat"))
    ivs <- read.table(pqtl,as.is=TRUE,header=TRUE)
    setwd(file.path(INF,"MS"))
    for(id in c("ieu-b-18","ieu-a-1025","ukb-b-17670","finn-a-G6_MS"))
    {
      pqtlMR(subset(ivs,SNP%in%c("rs2364485","rs1800693")),id,prefix=paste0(prot,"-",id,"-trans"))
      pqtlMR(subset(ivs,SNP%in%c("rs2229092","rs9263621")),id,prefix=paste0(prot,"-",id,"-cis"))
    }
  END
}

function pqtl_flanking()
# +/- 0.5Mb
{
  (
  # cis pQTLs
  # chr6:31540757_A_C rs2229092
  # chr6:31073047_A_G rs9263621
    export chr=6
    export rsid=rs2229092
    export pos=31540757
    pqtl_qtl_rsid
    export rsid=rs9263621
    export pos=31073047
    pqtl_qtl_rsid
  # trans pQTL
  # chr12:6514963_A_C rs2364485
  # chr12:6440009_C_T rs1800693
  # r2=0.0029
    export chr=12
    export rsid=rs2364485
    export pos=6514963
    pqtl_qtl_rsid
    export rsid=rs1800693
    export pos=6440009
    pqtl_qtl_rsid
  ) >> ${INF}/MS/${prot}-MS-MR.log
}

export M=1000000
export get_data=no

pqtl_flanking
pqtl

R --no-save -q <<END
  info <- function()
  {
    # https://gwas.mrcieu.ac.uk/datasets/?trait__icontains=multiple%20sclerosis
      ms_ids <- c("ukb-b-17670","ieu-b-18","ieu-a-1025","ebi-a-GCST005531","ieu-a-1024","ebi-a-GCST001198","ieu-a-820","ieu-a-821","finn-a-G6_MS")
      ms_gwasinfo <- ieugwasr::gwasinfo(ms_ids)
      cols <- c("id","year","ncase","ncontrol","nsnp","build","author","pmid","population")
      ms_gwasinfo[cols]
#                 id year ncase ncontrol     nsnp       build         author     pmid population
# 1      ukb-b-17670 2018  1679   461254  9851867 HG19/GRCh37   Ben Elsworth       NA   European
# 2         ieu-b-18 2019 47429    68374  6304359 HG19/GRCh37 Patsopoulos NA 31604244   European
# 3       ieu-a-1025 2013 14498    24091   156632 HG19/GRCh37        Beecham 24076602   European
# 4 ebi-a-GCST005531 2013 14498    24091   132089 HG19/GRCh37     Beecham AH 24076602   European
# 5       ieu-a-1024 2011  9722    17376   465435 HG19/GRCh37       Sawcer S 21833088   European
# 6 ebi-a-GCST001198 2011  9772    16849   463040 HG19/GRCh37       Sawcer S 21833088   European
# 7        ieu-a-820 2007   931     1862   327095 HG19/GRCh37      Hafler DA 17660530   European
# 8        ieu-a-821 2009   978      883   514572 HG19/GRCh37   Baranzini SE 19010793   European
# 9     finn-a-G6_MS 2020   378    44677 16152119 HG19/GRCh37             NA       NA   European
# 2,9 with no VCF
# https://gwas.mrcieu.ac.uk/files/ukb-b-17670/ukb-b-17670.vcf.gz
      mr_info <- as.data.frame(epigraphdb::mr(outcome_trait = "Multiple sclerosis", pval_threshold = 1e-8))
      select <- c("ukb-a-100","ukb-a-104","ieu-a-294","ieu-a-295","ieu-a-971")
      subset(mr_info,exposure.id%in%select,select=-c(outcome.trait,mr.selection,mr.method,mr.moescore))
# exposure.id                                             exposure.trait outcome.id       mr.b      mr.se      mr.pval
#   ieu-a-971                                         Ulcerative colitis ieu-a-1025 -0.4839367 0.06111146 2.395840e-15
#   ieu-a-295                                 Inflammatory bowel disease ieu-a-1025 -0.2593652 0.03799687 8.733802e-12
#   ieu-a-294                                 Inflammatory bowel disease ieu-a-1025  0.1232999 0.01697710 1.375433e-10
#   ukb-a-104 Non-cancer illness code  self-reported: ulcerative colitis ieu-a-1025 20.4956226 3.23034883 2.228468e-10
#   ukb-a-100          Non-cancer illness code  self-reported: psoriasis ieu-a-1025 11.0119963 1.86965907 3.865659e-09
  }
END
