# 31/8/2018 JHZ

setwd("u:/work")
mmp12 <- read.delim("MMP12.txt",as.is=TRUE)
mmp12 <- within(mmp12, {phen <- "MMP12";P <- 10^log.P.;N <- 3400})
dim(mmp12)
cad <- read.delim("CAD.txt",as.is=TRUE)
dim(cad)
cad_mmp12 <- merge(cad,mmp12,by.x="bp_hg19",by.y="position")
dim(cad_mmp12)
names(cad_mmp12)

library(TwoSampleMR)
# less ideal to obtain RSid's from cad and PhenoScanner actually does better
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
