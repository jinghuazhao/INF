# 2-12-2018 JHZ

z <- gzfile("/data/jinhua/1KGp3/1KGp3v5.txt.gz")
allfreq_ref_std <- read.table(z,header=TRUE,as.is=TRUE)
save(allfreq_ref_std,file="1KGp3v3.RData")

library(QCGWAS)
QCresults <- QC_GWAS("sumstats/INTERVAL/INTERVAL.IFN.gamma.gz",
		header_translations = "tryggve/header_translations.txt",
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.01,
		HQfilter_imp = 0.3,
		QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
		QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
		NAfilter = TRUE,
		allele_ref_std = "1KGp3v5.RData",
		allele_name_std = "OneKG",
		remove_mismatches = TRUE,
		allele_ref_alt = "ref_alternative.RData",
		allele_name_alt = "alternative",
		update_alt = TRUE,
		update_as_rdata = TRUE,
		backup_alt = TRUE)

QC_series(	data_files= c("data1.txt","data2.txt","data3.txt"),
		output_filenames = c("output1.txt","output2.txt","output3.txt"),
		dir_data = "preQC",
		dir_output = "work",
		dir_references = ".",
		header_translations = "tryggve/header_translations.txt",
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.01,
		HQfilter_imp = 0.3,
		QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
		QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
		NAfilter = TRUE,
		allele_ref_std = "ref_hapmap.RData",
		allele_name_std = "HapMap",
		remove_mismatches = TRUE,
		allele_ref_alt = "ref_alternative.RData",
		allele_name_alt = "alternative",
		update_alt = TRUE, update_as_rdata = TRUE, backup_alt = TRUE)
