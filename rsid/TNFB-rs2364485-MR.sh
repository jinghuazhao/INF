#!/usr/bin/bash

export chr=12
export start=6400000
export end=6520000
export M=0
export gene=LTBR
export prot=TNFB

echo SCALLOP/INF
# chr12:6514963_A_C
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
) > work/${prot}-pQTL.2s

gunzip -c data/discovery_metav3.0.meta.gz | \
awk 'NR>1 && $1==ENVIRON["chr"] && $2>=ENVIRON["start"] && $2<=ENVIRON["end"] {$3="chr" $1 ":" $2;print}' | \
sort -k3,3 | \
join -12 -23 work/snp_pos - | \
awk 'a[$1]++==0' > work/${prot}-QTL.2s

R --no-save -q <<END
prot <- Sys.getenv("prot")
ms <- within(read.table(paste0("work/",prot,"-QTL.2s"),as.is=TRUE,
             col.names=c("chrpos","rsid","chr","pos","A1","A2","N","P","OR")),
{
  beta <- log(P)
  se <- abs(beta/qnorm(P/2))
})
library(TwoSampleMR)
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
x <- read_exposure_data(paste0("work/",prot,"-pQTL.2s"),
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
yr <- extract_outcome_data(
       with(x,SNP),
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
# both y and yr miss rs1800693
h <- harmonise_data(x, y, action = 2)
subset(y,SNP%in%c("rs1800693","rs2364485"))
subset(yr,SNP%in%c("rs1800693","rs2364485"))
subset(h,SNP%in%c("rs1800693","rs2364485"))
xy <- mr(h)
END
