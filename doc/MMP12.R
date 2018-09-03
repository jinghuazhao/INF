# 3/9/2018 JHZ

setwd("u:/work")
mmp12 <- read.table("MMP12.dat",as.is=TRUE, col.names=c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "logP"))
mmp12 <- within(mmp12, {phen <- "MMP12";P <- 10^logP;N <- 3400})

library(TwoSampleMR)
exposure_dat <- format_data(cad_mmp12, type="exposure", snp_col = "markername", effect_allele_col = "Allele1", other_allele_col = "Allele2",
                            eaf_col = "effect_allele_freq", beta_col = "Effect", se_col = "StdErr", pval_col = "P", samplesize_col = "N")
ao <- available_outcomes()
subset(ao,consortium=="CARDIoGRAMplusC4D")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 7, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 
0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res_mr <- mr(dat)
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)
pdf("MMP12-MR.pdf")
mr_scatter_plot(res_mr, dat)
mr_forest_plot(res_single)
mr_leaveoneout_plot(res_loo)
mr_funnel_plot(res_single)

library(MendelianRandomization)
MRInputObject <- with(dat, mr_input(bx = beta.exposure, bxse = se.exposure, by = beta.outcome, byse = se.outcome,
                                    exposure = "MMP-12", outcome = "Coronary heart disease", snps = SNP))
mr_ivw(MRInputObject, model = "default", robust = FALSE, penalized = FALSE, weights = "simple", distribution = "normal", alpha = 0.05)
mr_egger(MRInputObject, robust = FALSE, penalized = FALSE, distribution = "normal", alpha = 0.05)
mr_maxlik(MRInputObject, model = "default", distribution = "normal", alpha = 0.05)
mr_median(MRInputObject, weighting = "weighted", distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
mr_allmethods(MRInputObject, method = "all")
mr_plot(MRInputObject, error = TRUE, orientate = FALSE, interactive = TRUE, labels = TRUE, line = "ivw")
dev.off()
