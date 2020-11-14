options(echo=TRUE,width=200)
INF <- Sys.getenv("INF")
source(paste0(INF,"/rsid/efo_inf.R"))
INF1_merge <- read.delim(paste0(INF,"/work/INF1.merge"),as.is=TRUE)
outdir <- paste0(INF,"/work/mr/")
library(TwoSampleMR)
for (type in c("cis","pan"))
{
  for(prot in with(INF1_merge,unique(prot)))
  {
    cat(type,prot,"\n")
    gz <- gzfile(paste0(outdir,prot,"-",type,".mrx"))
    d <- lapply(gz, function(x) tryCatch(read.delim(gz,as.is=TRUE), error=function(e) NULL))[[1]]
    if(nrow(d)==0) next
    d <- within(d,{P <- 10^logP})
    e <- format_data(d, type="exposure", phenotype_col="prot", header = TRUE, snp_col = "rsid",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                     samplesize_col = "N")
    exposure_dat <- clump_data(e)
    print(exposure_dat)
    for(outcomes in with(efo,MRBASEID))
    {
      cat(prot,"-",outcomes,"-",type,"\n")
      outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1,
                                          maf_threshold = 0.5)
      if(is.null(outcome_dat)) next
      dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
      if(nrow(dat)==0) next
      mr_result <- mr(dat)
      print(mr_result)
      mr_heterogeneity(dat)
      mr_pleiotropy_test(dat)
      mr_single <- mr_singlesnp(dat)
      mr_loo <- mr_leaveoneout(dat)
      save(dat,mr_result,mr_single,mr_loo,file=paste0(outdir,prot,"-",outcomes,"-",type,".mro"))
      pdf(paste0(outdir,prot,"-",outcomes,"-",type,".pdf"))
      mr_scatter_plot(mr_result, dat)
      mr_forest_plot(mr_single)
      mr_leaveoneout_plot(mr_loo)
      mr_funnel_plot(mr_single)
      dev.off()
    }
  }
}
