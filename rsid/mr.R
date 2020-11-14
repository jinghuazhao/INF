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
      result <- mr(dat)
      heterogeneity <- mr_heterogeneity(dat)
      pleiotropy <- mr_pleiotropy_test(dat)
      single <- mr_singlesnp(dat)
      loo <- mr_leaveoneout(dat)
      prefix <- paste0(outdir,outcomes,"-",prot,"-",type)
      invisible(lapply(c("result","heterogeneity","pleiotropy","single","loo"), function(x) {
                      v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                      if (!is.null(v)) write.table(format(v,digits=3),file=paste0(prefix,"-",x,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
                    }))
      pdf(paste0(prefix,".pdf"))
      mr_scatter_plot(result, dat)
      mr_forest_plot(single)
      mr_leaveoneout_plot(loo)
      mr_funnel_plot(single)
      dev.off()
    }
  }
}
