type <- Sys.getenv("type")
prot <- Sys.getenv("prot")
outcomes <- Sys.getenv("MRBASEID")
outdir <- "mr/"
gz <- gzfile(paste0(outdir,prot,"-",type,".mrx"))
d <- lapply(gz, function(x) tryCatch(read.delim(gz,as.is=TRUE), error=function(e) NULL))[[1]]
if (nrow(d)!=0)
{
  library(TwoSampleMR)
  d <- within(d,{P <- 10^logP})
  e <- format_data(d, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                   effect_allele_col = "Allele1", other_allele_col = "Allele2",
                   eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                   samplesize_col = "N")
  exposure_dat <- clump_data(e)
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1,
                                      maf_threshold = 0.5)
  if(!is.null(outcome_dat))
  {
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    if (nrow(dat)!=0)
    {
      mr_result <- mr(dat)
      mr_heterogeneity(dat)
      mr_pleiotropy_test(dat)
      mr_single <- mr_singlesnp(dat)
      mr_loo <- mr_leaveoneout(dat)
      prefix <- paste0(outdir,prot,"-",outcomes,"-",type)
      invisible(lapply(paste0("mr_",c("result","heterogeneity","pleiotropy_test","single","loo")), function(x) {
                      v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                      if (!is.null(v)) write.table(v,file=paste0(prefix,"-",x,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
                    }))
      save(dat,mr_result,mr_single,mr_loo,file=paste0(prefix,".mro"))
      pdf(paste0(prefix,".pdf"))
      mr_scatter_plot(mr_result, dat)
      mr_forest_plot(mr_single)
      mr_leaveoneout_plot(mr_loo)
      mr_funnel_plot(mr_single)
      dev.off()
    }
  }
}
