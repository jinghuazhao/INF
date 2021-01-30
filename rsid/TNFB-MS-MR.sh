#!/usr/bin/bash

R --no-save <<END
options(width=200)
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
END

function pqtl_qtl_mr()
{
  (
    awk -vOFS="\t" 'BEGIN{print "pheno", "rsid", "chr", "pos", "beta", "se", "snpid", "A1", "A2", "EAF", "P", "N"}'
    gunzip -c METAL/${prot}-1.tbl.gz | \
    sed '1d' | \
    awk -vFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
    {
      if ($1 == chr && $2 >= start-M && $2 <= end+M)
      {
        split($3,a,"_")
        print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,10^$12,$18
      }
    }' | \
    sort -k1,1 | \
    join -12 -21 work/snp_pos - | \
    awk 'a[$7]++==0' | \
    awk -vOFS="\t" '{print "TNFB", $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'
  ) > work/${prot}-pQTL-${rsid}.2s
  gunzip -c ~/rds/results/public/gwas/multiple_sclerosis/discovery_metav3.0.meta.gz | \
  awk 'NR>1 && $1==ENVIRON["chr"] && $2>=ENVIRON["start"] && $2<=ENVIRON["end"] {$3="chr" $1 ":" $2;print}' | \
  sort -k3,3 | \
  join -12 -23 work/snp_pos - | \
  awk 'a[$1]++==0' > work/${prot}-QTL-${rsid}.2s
  R --no-save -q <<\ \ END
    prot <- Sys.getenv("prot")
    rsid <- Sys.getenv("rsid")
    ms <- within(read.table(paste0("work/",prot,"-",rsid,"-QTL.2s"),as.is=TRUE,
                 col.names=c("chrpos","rsid","chr","pos","A1","A2","N","P","OR")),
    {
      beta <- log(P)
      se <- abs(beta/qnorm(P/2))
    })
    library(TwoSampleMR)
    x <- read_exposure_data(paste0("work/",prot,"-",rsid,"-pQTL.2s"),
         clump = FALSE,
         sep = "\t",
         phenotype_col = "pheno",
         snp_col = "rsid",
         beta_col = "beta",
         se_col = "se",
         eaf_col = "EAF",
         effect_allele_col = "A1",
         other_allele_col = "A2",
         pval_col = "P",
         samplesize_col = "N",
         gene_col = "gene",
         id_col = "rsid",
         log_pval = FALSE
       )
    y <- extract_outcome_data(
         with(ms,rsid),
         "ieu-b-18",
         proxies = TRUE,
         rsq = 1,
         align_alleles = 1,
         palindromes = 1,
         maf_threshold = 0.3,
         access_token = ieugwasr::check_access_token(),
         splitsize = 10000,
         proxy_splitsize = 500
       )
    subset(y,SNP%in%c("rs2364485"))
    h <- harmonise_data(x, y, action = 2)
    xy <- mr(h)
    print(xy)
  END
# rm work/${prot}-pQTL-${rsid}.2s work/${prot}-QTL-${rsid}.2s
}

export prot=TNFB
export gene=LTBR
# trans pQTL
# chr12:6514963_A_C rs2364485
# chr12:6440009_C_T rs1800693
# r2=0.0029
export rsid=rs2364485
export chr=12
export start=6400000
export end=6520000
export M=0
pqtl_qtl_mr
export rsid=rs1800693
export start=6434009
export end=6446009
pqtl_qtl_mr

# cis pQTLs
# chr6:31540757_A_C rs2229092
export rsid=rs2229092
export chr=6
export start=31534757
export end=31546757
pqtl_qtl_mr

# chr6:31073047_A_G rs9263621
export rsid=rs9263621
export start=31067047
export end=31079047
pqtl_qtl_mr
