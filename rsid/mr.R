library(TwoSampleMR)
library(pQTLtools)
type <- Sys.getenv("type")
prot <- Sys.getenv("prot")
outcomes <- Sys.getenv("MRBASEID")
INF <- Sys.getenv("INF")
outdir <- file.path(INF,"mr",type)
gz <- gzfile(file.path(outdir,paste0(prot,"-",type,".mrx")))
d <- lapply(gz, function(x) tryCatch(read.delim(gz,as.is=TRUE), error=function(e) NULL))[[1]]
d <- within(d,{P <- 10^logP})
e <- format_data(d, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                 effect_allele_col = "Allele1", other_allele_col = "Allele2",
                 eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                 samplesize_col = "N")
exposure_dat <- clump_data(e)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1,
                                    maf_threshold = 0.5)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
prefix <- file.path(outdir,paste0(outcomes,"-",prot,"-",type))
pdf(paste0(prefix,".pdf"))
run_TwoSampleMR(list(exposure=e, outcome=outcome_dat, clump=exposure_dat, harmonise=dat), plot="New", prefix=prefix)
dev.off()
