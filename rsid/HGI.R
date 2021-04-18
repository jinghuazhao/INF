options(width=500, echo=FALSE)
library(dplyr)
library(pQTLtools)
library(TwoSampleMR)

INF <- Sys.getenv("INF")
trait <- Sys.getenv("trait"); prot <- Sys.getenv("prot"); rsid <- Sys.getenv("rsid");
id3 <- paste(trait,prot,rsid,sep="-")
cat(id3,"\n",sep="")

Allele1 <- Sys.getenv("Allele1")
Allele2 <- Sys.getenv("Allele2")
EAF <- Sys.getenv("EAF")
Effect <- Sys.getenv("Effect")
StdErr <- Sys.getenv("StdErr")
P <- Sys.getenv("P")
N <- Sys.getenv("N")

MR <- function(clumping=FALSE, r2=0.01)
{
  y <- within(read.delim(file.path(INF,"HGI","mr",id3)), {outcome=prot})
  if (nrow(y)<=1) return
  o <- format_data(y, type="outcome", header = TRUE, phenotype_col = "outcome", snp_col = "rsid",
                   effect_allele_col = "ALT", other_allele_col = "REF",
                   eaf_col = "all_meta_AF", 
                   beta_col = "all_inv_var_meta_beta",
                   se_col = "all_inv_var_meta_sebeta",
                   pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                   samplesize_col = "all_meta_sample_N")
  x <- read.delim(file.path(INF,"HGI","mr",paste0(id3,".gz")))
  e <- format_data(x, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                   effect_allele_col = "Allele1", other_allele_col = "Allele2",
                   eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "logP", log_pval = TRUE,
                   samplesize_col = "N") %>% mutate(pval.exposure=pval)
  if (clumping)
  {
    e <- clump_data(e,clump_r2 = r2)
    if (nrow(e)==0) return
  }
  d <- merge(e,o,by="SNP")
  if (nrow(d)<=1) return
  e <- subset(e,SNP%in%d[["SNP"]])
  o <- subset(o,SNP%in%d[["SNP"]])
  if (nrow(e)==0 | nrow(o)==0) return
  dat <- harmonise_data(e, o, action = 1)
  cat(nrow(dat), "\n")
  result <- mr(dat)
  heterogeneity <- mr_heterogeneity(dat)
  pleiotropy <- mr_pleiotropy_test(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)
  prefix <- paste0("MR-",id3)
  invisible(lapply(c("result","heterogeneity","pleiotropy","single","loo"), function(x) {
                 v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                 if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"HGI","mr",paste0(prefix,"-",x,".txt")),
                                              quote=FALSE,row.names=FALSE,sep="\t")
         }))
}

pqtlMR <- function(r2=0.001)
{
  cat(id3,"\n")
  e <- format_data(data.frame(SNP=rsid,Allele1,Allele2,EAF,Effect,StdErr,P,N), type="exposure", header = TRUE, snp_col = "SNP",
                   effect_allele_col = "Allele1", other_allele_col = "Allele2",
                   eaf_col = "EAF", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                   samplesize_col = "N")
  d <- within(read.delim(file.path(INF,"HGI","mr",id3)),{outcome=prot})
  o <- format_data(d, type="outcome", header = TRUE,  phenotype_col = "outcome", snp_col = "rsid",
                   effect_allele_col = "ALT", other_allele_col = "REF",
                   eaf_col = "all_meta_AF", 
                   beta_col = "all_inv_var_meta_beta",
                   se_col = "all_inv_var_meta_sebeta",
                   pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                   samplesize_col = "all_meta_sample_N")
  outcome_dat <- clump_data(o,clump_r2 = r2)
  if (nrow(e) != nrow(subset(o,SNP==rsid))) return
  dat <- harmonise_data(e, subset(o,SNP==rsid), action = 1)
  directionality <- directionality_test(dat)
  result <- mr(dat)
  heterogeneity <- mr_heterogeneity(dat)
  pleiotropy <- mr_pleiotropy_test(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)
  prefix <- paste0("pqtlMR-",id3)
  invisible(lapply(c("directionality","result","heterogeneity","pleiotropy","single","loo"), function(x) {
                  v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                  if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"HGI","mr",paste0(prefix,"-",x,".txt")),
                                               quote=FALSE,row.names=FALSE,sep="\t")
             }))
}

# MR(clumping=TRUE,r2=0.001)
pqtlMR()
