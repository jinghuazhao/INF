# 3-12-2018 JHZ

library(QCGWAS)
QC_series(data_files=c("INTERVAL.ARTN.gz","ORCADES.ARTN.gz","STABILITY.ARTN.gz"),
	dir_data="sumstats/work",
	dir_output="work",
	dir_references="/data/jinhua/data/1KG",
	output_filenames = c("INTERVAL.ARTN.gz","ORCADES.ARTN.gz","STABILITY.ARTN.gz"),
	header_translations = "tryggve/header_translations.txt",
	save_final_dataset = TRUE,
	HQfilter_FRQ = 0.01,
	HQfilter_imp = 0.3,
	QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
	QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
	NAfilter = TRUE,
	allele_ref_std = "1KGp3v5.RData",
	allele_name_std = "OneKG",
	remove_mismatches = TRUE)
