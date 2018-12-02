# 2-12-2018 JHZ

header_translations <- read.delim("tryggve/header_translations.txt",as.is=TRUE)

library(QCGWAS)
QC_series(data_files=c("INTERVAL.ARTN.gz","INTERVAL.IFN.gamma.gz","INTERVAL.IL.13.gz"),
	dir_data="sumstats/INTERVAL",
	dir_output="work",
	dir_references=".",
	header_translations = header_translations,
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
