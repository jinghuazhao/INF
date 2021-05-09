library(TwoSampleMR)
library(pQTLtools)
library(dplyr)
type <- Sys.getenv("type")
prot <- Sys.getenv("prot")
outcomes <- Sys.getenv("MRBASEID")
INF <- Sys.getenv("INF")
outdir <- file.path(INF,"mr",type)
gz <- gzfile(file.path(outdir,paste0(prot,"-",type,".mrx")))
d <- lapply(gz, function(x) tryCatch(read.delim(gz,as.is=TRUE), error=function(e) NULL))[[1]] %>%
     group_by(rsid) %>%
     slice(which.max(abs(Effect/StdErr)))
e <- format_data(within(d, {P <- 10^logP}), type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                 effect_allele_col = "Allele1", other_allele_col = "Allele2",
                 eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                 samplesize_col = "N")
exposure <- clump_data(e,clump_r2 = 0.01)
outcome <- extract_outcome_data(exposure$SNP, outcomes, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.5)
harmonise <- harmonise_data(exposure, outcome, action = 2)
prefix <- file.path(outdir,paste0(outcomes,"-",prot,"-",type))
pdf(paste0(prefix,".pdf"))
run_TwoSampleMR(harmonise, mr_plot="New", prefix=prefix)
dev.off()
