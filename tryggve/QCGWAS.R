# 3-12-2018 JHZ

library(QCGWAS)

prot <- Sys.getenv("protein")
src <- "INTERVAL"
src_file <- paste(src, prot, "gz", sep=".")
EGCUT <- paste0("EGCUT_INF",c("_autosomal","_male","_female"))
studies <- c(EGCUT, "INTERVAL", "NSPHS_INF", "ORCADES", "STABILITY", "STANLEY", "VIS")
QC_files <- paste(studies, prot, "gz", sep=".")
header_translations <- read.delim("tryggve/header_translations.tsv",as.is=TRUE)

QC_GWAS(src_file,
        dir_data = "sumstats/work",
        dir_output = "work",
        dir_references = "/data/jinhua/data/EasyQC",
        header_translations = header_translations,
        save_final_dataset = FALSE,
        HQfilter_FRQ = 0.01,
        HQfilter_imp = 0.3,
        QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
        QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
        NAfilter = TRUE,
        allele_ref_std = "1KGp3v5.RData",
        allele_name_std = "KG",
        remove_mismatches = TRUE,
        allele_ref_alt = NULL,
        allele_name_alt = "alternative",
        update_alt = TRUE,    
        update_savename = "ref_alternative",
        update_as_rdata = TRUE)

QC_series(data_files = QC_files,
	dir_data = "sumstats/work",
	dir_output = "work",
	dir_references = "/data/jinhua/data/EasyQC",
	output_filenames = QC_files,
	header_translations = header_translations,
	save_final_dataset = TRUE,
	HQfilter_FRQ = 0.01,
	HQfilter_imp = 0.3,
	QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
	QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
	NAfilter = TRUE,
	allele_ref_std = "1KGp3v5.RData",
	allele_name_std = "KG",
	remove_mismatches = TRUE,
	allele_ref_alt = "ref_alternative.RData",
	allele_name_alt = "alternative",
	update_alt = TRUE,
	update_savename = "ref_alternative",
	update_as_rdata = TRUE)
