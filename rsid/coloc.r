prot <- Sys.getenv("prot")
pQTL <- Sys.getenv("pQTL")
M <- Sys.getenv("M")

library(gap)
isnpid <- within(inv_chr_pos_a1_a2(pQTL),
{
  chr <- gsub("chr","",chr)
  pos <- as.integer(pos)
  start <- pos-M
  if (start<0) start <- 0
  end <- pos+M
})
region <- with(isnpid,paste0(chr,":",start,"-",end))

library(pQTLtools)
lapply(c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer"), require, character.only = TRUE)
tfp <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(tfp, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
tfpi <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_imported.tsv")
imported_tabix_paths <- read.delim(tfpi, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

pdf(file.path(INF,"coloc",paste0(prot,"-",pQTL,".pdf")))
# Association at a specific locus
platelet_df <- dplyr::filter(tabix_paths, study == "CEDAR", tissue_label == "platelet")
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
column_names <- names(read.delim(hdr))
summary_stats <- import_eQTLCatalogue(platelet_df$ftp_path, region, selected_gene_id = "ENSG00000163947", column_names)
summary_stats
ggplot(summary_stats, aes(x = position, y = -log(pvalue, 10))) + geom_point()
# GWAS sumstat from the same region
gwasvcf::set_bcftools("/rds/user/jhz22/hpc-work/bin/bcftools")
gwas_stats <- gwasvcf::query_gwas(file.path(INF,"METAL/vcf",paste0(prot,".vcf.gz"), chrompos = region)
gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats) %>% keepSeqlevels("3") %>% renameSeqlevels("chr3")
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)
gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
  unlist() %>%
  renameSeqlevels("3") %>%
  dplyr::as_tibble() %>%
  dplyr::transmute(chromosome = seqnames, position = start, AF, ES, SE, LP, SS) %>%
  dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
  dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(row_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(row_count == 1)
ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
# Colocalisation
res <- run_coloc(summary_stats, gwas_stats_hg38)

# a. all other eQTL datasets
microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
column_names <- names(read.delim(hdr))
summary_list <- purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region, selected_gene_id = "ENSG00000163947", column_names))
coloc_df_microarray <- purrr::map_df(summary_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

# b. Uniformly processed RNA-seq datasets
rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
column_names <- names(read.delim(hdr))
safe_import <- purrr::safely(import_eQTLCatalogue)
summary_list <- purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = "ENSG00000163947", column_names))
result_list <- purrr::map(summary_list, ~.$result)
result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
coloc_df_rnaseq <- purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

# c. GTEx_v8 imported eQTL datasets
rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
column_names <- names(read.delim(hdr))
safe_import <- purrr::safely(import_eQTLCatalogue)
summary_list <- purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = "ENSG00000163947", column_names))
result_list <- purrr::map(summary_list, ~.$result)
result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
result_filtered <- purrr::map(result_list, ~dplyr::filter(., !is.na(se)))
coloc_df_imported <- purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")

coloc_df = dplyr::bind_rows(coloc_df_microarray, coloc_df_rnaseq, coloc_df_imported)
dplyr::arrange(coloc_df, -PP.H4.abf)
ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
dev.off()
