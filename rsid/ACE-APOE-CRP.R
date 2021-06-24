# 31/8/2018 JHZ

# location for mrbase.oauth and pdf
setwd("u:/work")

library(MRInstruments)
d <- subset(proteomic_qtls,analyte%in%c("ACE","ApoE","CRP"))
d <- within(d, {N=5000})

library(TwoSampleMR)
exposure_dat <- format_data(d, type="exposure", snp_col = "SNP", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                            eaf_col = "eaf", beta_col = "beta", se_col = "se", pval_col = "pval", samplesize_col = "N")
ao <- available_outcomes()
subset(ao,consortium=="CARDIoGRAMplusC4D")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 7, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res_mr <- mr(dat)
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)
pdf("ACE-ApoE-CRP.pdf")
mr_scatter_plot(res_mr, dat)
mr_forest_plot(res_single)
mr_leaveoneout_plot(res_loo)
mr_funnel_plot(res_single)

library(MendelianRandomization)
MRInputObject <- with(dat, mr_input(bx = beta.exposure, bxse = se.exposure, by = beta.outcome, byse = se.outcome,
                                    exposure = "ACE-ApoE-CRP", outcome = "Coronary heart disease", snps = SNP))
mr_ivw(MRInputObject, model = "default", robust = FALSE, penalized = FALSE, weights = "simple", distribution = "normal", alpha = 0.05)
mr_egger(MRInputObject, robust = FALSE, penalized = FALSE, distribution = "normal", alpha = 0.05)
mr_maxlik(MRInputObject, model = "default", distribution = "normal", alpha = 0.05)
mr_median(MRInputObject, weighting = "weighted", distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
mr_allmethods(MRInputObject, method = "all")
mr_plot(MRInputObject, error = TRUE, orientate = FALSE, interactive = TRUE, labels = TRUE, line = "ivw")
dev.off()

# documentation example
MRMVInputObject <- mr_mvinput(bx = cbind(ldlc, hdlc, trig),
                              bxse = cbind(ldlcse, hdlcse, trigse),
                              by = chdlodds, 
                              byse = chdloddsse)

MRMVInputObject
MRMVObject <- mr_mvivw(MRMVInputObject, 
                          model = "default",
			  correl = FALSE,
                          distribution = "normal",
                          alpha = 0.05)

MRMVObject <- mr_mvivw(MRMVInputObject)
MRMVObject

# PhenoScanner
path.noproxy <- system.file("extdata", "vitD_snps_PhenoScanner.csv", package = "MendelianRandomization")
path.proxies <- system.file("extdata", "vitD_snps_PhenoScanner_proxies.csv", package = "MendelianRandomization")
extract.pheno.csv(exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
                  outcome = "Tanner stage", pmidO = 24770850, ancestryO = "European", file = path.noproxy)
extract.pheno.csv(exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
                  outcome = "Tanner stage", pmidO = 24770850, ancestryO = "European", rsq.proxy = 0.6, file = path.proxies)
