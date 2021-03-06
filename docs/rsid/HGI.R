pqtlMR <- function(r2=0.001)
{
  cat(id3,"\n")
  y <- within(read.delim(file.path(INF,"HGI","mr",paste0(id3,".gz"))),{outcome=prot})
  o <- format_data(y, type="outcome", header = TRUE, phenotype_col = "outcome", snp_col = "rsid",
                   effect_allele_col = "ALT", other_allele_col = "REF",
                   eaf_col = "all_meta_AF", 
                   beta_col = "all_inv_var_meta_beta",
                   se_col = "all_inv_var_meta_sebeta",
                   pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                   samplesize_col = "all_meta_sample_N")
  Allele1 <- Sys.getenv("Allele1")
  Allele2 <- Sys.getenv("Allele2")
  EAF <- as.numeric(Sys.getenv("EAF"))
  Effect <- as.numeric(Sys.getenv("Effect"))
  StdErr <- as.numeric(Sys.getenv("StdErr"))
  P <- as.numeric(Sys.getenv("P"))
  N <- as.integer(Sys.getenv("N"))
  x <- data.frame(prot,SNP=rsid,Allele1,Allele2,EAF,Effect,StdErr,P,N)
  e <- format_data(x, type="exposure", header = TRUE, phenotype_col="prot", snp_col = "SNP",
                   effect_allele_col = "Allele1", other_allele_col = "Allele2",
                   eaf_col = "EAF", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                   samplesize_col = "N")
  print(head(e))
  print(head(o))
  dat <- harmonise_data(e, o, action = 1)
  result <- mr(dat)
  heterogeneity <- mr_heterogeneity(dat)
  pleiotropy <- mr_pleiotropy_test(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)
  prefix <- paste0("pqtlMR-",id3)
  invisible(lapply(c("result","heterogeneity","pleiotropy","single","loo"), function(x) {
                  v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                  if (!is.null(v)) write.table(format(v,digits=3),file=file.path(INF,"HGI","mr",paste0(prefix,"-",x,".txt")),
                                               quote=FALSE,row.names=FALSE,sep="\t")
             }))
}

MR <- function(clumping=TRUE, r2=0.001)
{
  y <- within(read.delim(file.path(INF,"HGI","mr",paste0(id3,".gz"))), {outcome=prot})
  if (nrow(y)<=1) return (-1)
  o <- format_data(y, type="outcome", header = TRUE, phenotype_col = "outcome", snp_col = "rsid",
                   effect_allele_col = "ALT", other_allele_col = "REF",
                   eaf_col = "all_meta_AF", 
                   beta_col = "all_inv_var_meta_beta",
                   se_col = "all_inv_var_meta_sebeta",
                   pval_col = "all_inv_var_meta_p", log_pval = FALSE,
                   samplesize_col = "all_meta_sample_N")
  x <- read.delim(file.path(INF,"HGI","mr","tsv",paste0(id3,".tsv.gz")))
  e <- format_data(x, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                   effect_allele_col = "Allele1", other_allele_col = "Allele2",
                   eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                   samplesize_col = "N")
  if (clumping)
  {
    e <- clump_data(e,clump_r2 = r2)
    if (nrow(e)==0) return (-2)
  }
  d <- merge(e,o,by="SNP")
  if (nrow(d)<=1) return (-3)
  e <- subset(e,SNP%in%d[["SNP"]])
  o <- subset(o,SNP%in%d[["SNP"]])
  if (nrow(e)==0 | nrow(o)==0) return (-4)
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

options(width=500, echo=FALSE)
library(dplyr)
library(pQTLtools)
library(TwoSampleMR)

INF <- Sys.getenv("INF")
trait <- Sys.getenv("trait"); prot <- Sys.getenv("prot"); rsid <- Sys.getenv("rsid");
id3 <- paste("r6",trait,prot,rsid,sep="-")
cat(id3,"\n",sep="")

MR()
pqtlMR()
