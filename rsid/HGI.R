options(width=500)
library(dplyr)
library(pQTLtools)
library(TwoSampleMR)

INF <- Sys.getenv("INF")
trait <- Sys.getenv("trait")

MR <- function(clumping=FALSE)
{
  f <- file.path(INF,"HGI","mr","INF.ins")
  ivs <- read.delim(f,sep=" ")
  for(row in 1:nrow(ivs))
  {
    prot <- ivs[row,"Phenotype"]
    rsid <- ivs[row,"SNP"]
    cat(row, trait, "-", prot,"-",rsid,"\n")
    d <- read.delim(file.path(INF,"HGI","mr",paste0(trait,"-",prot,"-",rsid,".mrx")))
    d <- within(d,{P <- 10^logP})
    e <- format_data(d, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                     samplesize_col = "N")
    d <- read.delim(file.path(INF,"HGI","mr",paste0(trait,"-",prot,"-",rsid)))
    if (nrow(d)<=1) next
    o <- format_data(d, type="outcome", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "ALT", other_allele_col = "REF",
                     eaf_col = "all_meta_AF", 
                     beta_col = "all_inv_var_meta_beta",
                     se_col = "all_inv_var_meta_sebeta",
                     pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                     samplesize_col = "all_meta_sample_N")
    if (clumping)
    {
      outcome_dat <- clump_data(o,clump_r2 = 0.01)
      if (nrow(outcome_dat)==0) next
      d <- merge(e,outcome_dat,by="SNP")
      if (nrow(d)<=1) next
      e <- subset(e,SNP%in%with(d,SNP))
      o <- subset(outcome_dat,SNP%in%with(d,SNP))
      if (nrow(e)==0 | nrow(outcome_dat)==0) next
      dat <- harmonise_data(e, outcome_dat, action = 1)
    } else {
      d <- merge(e,o,by="SNP")
      if (nrow(d)<=1) next
      e <- subset(e,SNP%in%with(d,SNP))
      o <- subset(o,SNP%in%with(d,SNP))
      if (nrow(e)==0 | nrow(o)==0) next
      dat <- harmonise_data(e, o, action = 1)
    }
    cat(nrow(dat), "\n")
    directionality <- directionality_test(dat)
    result <- mr(dat)
    heterogeneity <- mr_heterogeneity(dat)
    pleiotropy <- mr_pleiotropy_test(dat)
    single <- mr_singlesnp(dat)
    loo <- mr_leaveoneout(dat)
    prefix <- paste0("MR-",trait,"-",prot,"-",rsid)
    invisible(lapply(c("directionality","result","heterogeneity","pleiotropy","single","loo"), function(x) {
                    v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                    if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"HGI","mr",paste0(prefix,"-",x,".txt")),
                                                 quote=FALSE,row.names=FALSE,sep="\t")
               }))
  }
}
MR()

pqtlMR <- function()
{
  f <- file.path(INF,"HGI","mr","INF.ins")
  ivs <- read.delim(f,sep=" ")
  for(row in 1:nrow(ivs))
  {
    prot <- ivs[row,"Phenotype"]
    rsid <- ivs[row,"SNP"]
    cat(row, trait, "-", prot,"-",rsid,"\n")
    e <- format_data(ivs[row,], type="exposure", header = TRUE, snp_col = "SNP",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "EAF", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                     samplesize_col = "N")
    d <- read.delim(file.path(INF,"HGI","mr",paste0(trait,"-",prot,"-",rsid)))
    o <- format_data(d, type="outcome", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "ALT", other_allele_col = "REF",
                     eaf_col = "all_meta_AF", 
                     beta_col = "all_inv_var_meta_beta",
                     se_col = "all_inv_var_meta_sebeta",
                     pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                     samplesize_col = "all_meta_sample_N")
    outcome_dat <- clump_data(o,clump_r2 = 0.01)
    if (nrow(e) != nrow(subset(o,SNP==ivs[row,"SNP"]))) next
    dat <- harmonise_data(e, subset(o,SNP==ivs[row,"SNP"]), action = 1)
    directionality <- directionality_test(dat)
    result <- mr(dat)
    heterogeneity <- mr_heterogeneity(dat)
    pleiotropy <- mr_pleiotropy_test(dat)
    single <- mr_singlesnp(dat)
    loo <- mr_leaveoneout(dat)
    prefix <- paste0("pqtlMR-",trait,"-",prot,"-",rsid)
    invisible(lapply(c("directionality","result","heterogeneity","pleiotropy","single","loo"), function(x) {
                    v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                    if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"HGI","mr",paste0(prefix,"-",x,".txt")),
                                                 quote=FALSE,row.names=FALSE,sep="\t")
               }))
  }
}
pqtlMR()
