options(width=500)
library(dplyr)
library(pQTLtools)
library(TwoSampleMR)
INF <- Sys.getenv("INF")

pqtlMR <- function()
{
  f <- file.path(INF,"work","HGI","INF.ins")
  ivs <- read.delim(f,sep=" ")
  for(row in 1:nrow(ivs))
  {
    prot <- ivs[row,"Phenotype"]
    rsid <- ivs[row,"SNP"]
    e <- format_data(ivs[row,], type="exposure", header = TRUE, snp_col = "SNP",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "EAF", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                     samplesize_col = "N")
    d <- read.delim(file.path(INF,paste0(prot,"-",rsid)))
    o <- format_data(d, type="outcome", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "ALT", other_allele_col = "REF",
                     eaf_col = "all_meta_AF", 
                     beta_col = "all_inv_var_meta_beta",
                     se_col = "all_inv_var_meta_sebeta",
                     pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                     samplesize_col = "all_meta_sample_N")
    outcome_dat <- clump_data(o)
    if (nrow(e) != nrow(subset(o,SNP==ivs[row,"SNP"]))) next
    dat <- harmonise_data(e, subset(o,SNP==ivs[row,"SNP"]), action = 1)
    directionality <- directionality_test(dat)
    result <- mr(dat)
    heterogeneity <- mr_heterogeneity(dat)
    pleiotropy <- mr_pleiotropy_test(dat)
    single <- mr_singlesnp(dat)
    loo <- mr_leaveoneout(dat)
    prefix <- paste0("pqtlMR-",prot,"-",rsid)
    invisible(lapply(c("directionality","result","heterogeneity","pleiotropy","single","loo"), function(x) {
                    v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                    if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"work","HGI",paste0(prefix,"-",x,".txt")),
                                                 quote=FALSE,row.names=FALSE,sep="\t")
               }))
  }
}

MR <- function()
{
  f <- file.path(INF,"work","HGI","INF.ins")
  ivs <- read.delim(f,sep=" ")
  for(row in 1:nrow(ivs))
  {
    prot <- ivs[row,"Phenotype"]
    rsid <- ivs[row,"SNP"]
    d <- read.delim(file.path(INF,"work","HGI",paste0(prot,"-",rsid,".mrx")))
    d <- within(d,{P <- 10^logP})
    e <- format_data(d, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                     samplesize_col = "N")
    d <- read.delim(file.path(INF,"work","HGI",paste0(prot,"-",rsid,"_1e-5")))
    if (nrow(d)==0) next
    o <- format_data(d, type="outcome", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "ALT", other_allele_col = "REF",
                     eaf_col = "all_meta_AF", 
                     beta_col = "all_inv_var_meta_beta",
                     se_col = "all_inv_var_meta_sebeta",
                     pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                     samplesize_col = "all_meta_sample_N")
    outcome_dat <- clump_data(o)
    d <- merge(e,o,by="SNP")
    e <- subset(e,SNP%in%with(d,SNP))
    o <- subset(o,SNP%in%with(d,SNP))
    if (nrow(e)==0 | nrow(o)==0) next
    dat <- harmonise_data(e, o, action = 1)
    directionality <- directionality_test(dat)
    result <- mr(dat)
    heterogeneity <- mr_heterogeneity(dat)
    pleiotropy <- mr_pleiotropy_test(dat)
    single <- mr_singlesnp(dat)
    loo <- mr_leaveoneout(dat)
    prefix <- paste0("MR-",prot,"-",rsid)
    invisible(lapply(c("directionality","result","heterogeneity","pleiotropy","single","loo"), function(x) {
                    v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                    if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"work","HGI",paste0(prefix,"-",x,".txt")),
                                                 quote=FALSE,row.names=FALSE,sep="\t")
               }))
  }
}

pqtlMR()
MR()
