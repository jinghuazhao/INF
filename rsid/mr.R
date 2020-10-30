prot <- Sys.getenv("prot")
suffix <- Sys.getenv("suffix")
outcomes <- Sys.getenv("outcomes")

library(TwoSampleMR)
gz <- gzfile(paste0("work/",prot,".mrx",suffix))
d <- within(read.delim(gz,as.is=TRUE),{P <- 10^-logP})
exposure_dat <- format_data(d, type="exposure", snp_col = "rsid", effect_allele_col = "Allele1", other_allele_col = "Allele2",
                            eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", pval_col = "P", log_pval = FALSE,
                            samplesize_col = "N")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1,
                                    maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res_mr <- mr(dat)
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)
save(res_mr,res_single,res_loo,file=paste0("work/",prot,".mro",suffix),quote=FALSE,row.names=FALSE)
pdf(paste0("work/",prot,suffix,".pdf"))
mr_scatter_plot(res_mr, dat)
mr_forest_plot(res_single)
mr_leaveoneout_plot(res_loo)
mr_funnel_plot(res_single)
dev.off()

heightMR <- function()
{
  gz <- gzfile("~/hpc-work/results/BMI/Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt.gz")
  d <- read.delim(gz,as.is=TRUE)
  exposure_dat <- format_data(d, type="exposure", snp_col = "SNP", effect_allele_col = "Tested_Allele", other_allele_col = "Other_Allele",
                              eaf_col = "Freq_Tested_Allele_in_HRS", beta_col = "BETA_COJO", se_col = "SE_COJO", pval_col = "P_COJO", 
                              samplesize_col = "N")
  ao <- available_outcomes()
  subset(ao,consortium=="DIAGRAM")
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, 23, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
  dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
  res_mr <- mr(dat)
  mr_heterogeneity(dat)
  mr_pleiotropy_test(dat)
  res_single <- mr_singlesnp(dat)
  res_loo <- mr_leaveoneout(dat)
  pdf("heightMR.pdf")
  mr_scatter_plot(res_mr, dat)
  mr_forest_plot(res_single)
  mr_leaveoneout_plot(res_loo)
  mr_funnel_plot(res_single)

  library(MendelianRandomization)
  MRInputObject <- with(dat, mr_input(bx = beta.exposure, bxse = se.exposure, by = beta.outcome, byse = se.outcome,
                                      exposure = "Body mass index", outcome = "T2D", snps = SNP))
  mr_ivw(MRInputObject, model = "default", robust = FALSE, penalized = FALSE, weights = "simple", distribution = "normal", alpha = 0.05)
  mr_egger(MRInputObject, robust = FALSE, penalized = FALSE, distribution = "normal", alpha = 0.05)
  mr_maxlik(MRInputObject, model = "default", distribution = "normal", alpha = 0.05)
  mr_median(MRInputObject, weighting = "weighted", distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
  mr_allmethods(MRInputObject, method = "all")
  mr_plot(MRInputObject, error = TRUE, orientate = FALSE, interactive = TRUE, labels = TRUE, line = "ivw")
  dev.off()
}
