#!/rds/user/jhz22/hpc-work/bin/Rscript --vanilla

liftRegion <- function(x,chain)
{
  require(GenomicRanges)
  gr <- with(x,GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end)))
  seqlevelsStyle(gr) <- "UCSC"
  gr38 <- rtracklayer::liftOver(gr, chain)
  chr <- colnames(table(seqnames(gr38)))
  chr <- gsub("chr","",chr)
  start <- min(unlist(start(gr38)))
  end <- max(unlist(end(gr38)))
# chr <- as.character(seqnames(gr38))
# start <- as.character(start(gr38))
# end <- as.character(end(gr38))
  invisible(list(chr=chr,start=start,end=end,region=paste0(chr,":",start,"-",end)))
}

coloc_sumstats <- function(prot)
{
  cat("GWAS sumstats\n")
  gwas_stats <- gwasvcf::query_gwas(file.path(INF,"METAL/vcf",paste0(prot,".vcf.gz")), chrompos = region37)
  gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats) %>% keepSeqlevels(chr) %>% renameSeqlevels(paste0("chr",chr))
  gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
    unlist() %>%
    renameSeqlevels(chr) %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(chromosome = seqnames, position = start, AF, ES, SE, LP, SS) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
    dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(row_count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(row_count == 1)
  ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  gwas_stats_hg38
}

coloc_a <- function(gwas_stats_hg38)
{
  cat("a. eQTL datasets\n")
  microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
  column_names <- names(read.delim(hdr))
  summary_list <- purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region38, selected_gene_id = ensGene, column_names))
  coloc_df_microarray <- purrr::map_df(summary_list[lapply(summary_list,nrow)!=0], ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

coloc_b <- function(gwas_stats_hg38)
{
  cat("b. Uniformly processed RNA-seq datasets\n")
  rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list, ~safe_import(., region38, selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list[lapply(result_list,nrow)!=0], ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  coloc_df_rnaseq <- purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

coloc_c <- function(gwas_stats_hg38)
{
  cat("c. GTEx_v8 imported eQTL datasets\n")
  rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>% dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list, ~safe_import(., region38, selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list, ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0], ~dplyr::filter(., !is.na(se)))
  coloc_df_imported <- purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

mixed_coloc <- function(prot,chr,ensGene,chain,region37,region38,out,run_all=FALSE)
{
  pdf(paste0(out,".pdf"))
  gwas_stats_hg38 <- coloc_sumstats(prot)
  if (run_all)
  {
    coloc_df_microarray <- coloc_a(gwas_stats_hg38)
    coloc_df_rnaseq <- coloc_b(gwas_stats_hg38)
    coloc_df_imported <- coloc_c(gwas_stats_hg38)
    if (exists("coloc_df_microarray") & exits("coloc_df_rnaseq") & exists("coloc_df_imported"))
    {
      coloc_df = dplyr::bind_rows(coloc_df_microarray, coloc_df_rnaseq, coloc_df_imported)
      saveRDS(coloc_df, file=paste(out,".RDS"))
      dplyr::arrange(coloc_df, -PP.H4.abf)
      ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
    }
  } else {
    coloc_df_imported <- coloc_c(gwas_stats_hg38)
    if (exists("coloc_df_imported")
    {
      saveRDS(coloc_df_imported,file=paste(out,".RDS"))
      dplyr::arrange(coloc_df_imported, -PP.H4.abf)
      ggplot(coloc_df_imported, aes(x = PP.H4.abf)) + geom_histogram()
    }
  }
  dev.off()
}

single_run <- function(r)
{
  sentinel <- sentinels[r,]
  isnpid <- within(gap::inv_chr_pos_a1_a2(sentinel[["SNP"]]),
  {
    chr <- gsub("chr","",chr)
    pos <- as.integer(pos)
    start <- pos-M
    if (start<0) start <- 0
    end <- pos+M
  })
  chr <- with(isnpid,chr)
  region37 <- with(isnpid, paste0(chr,":",start,"-",end))
  ensRegion37 <- with(subset(inf1,prot==sentinel[["prot"]]),paste0(chr,":",start,"-",end))
  region38 <- with(liftRegion(isnpid,chain),region)
  ensGene <- subset(inf1,prot==sentinel[["prot"]])[["ensembl_gene_id"]]
  ensRegion38 <- with(liftRegion(subset(inf1,prot==sentinel[["prot"]]),chain),region)
  f <- file.path(INF,"coloc",with(sentinel,paste0(prot,"-",SNP)))
  cat(chr,region37,region38,ensGene,ensRegion37,ensRegion38,"\n")
  mixed_coloc(sentinel[["prot"]],chr,ensGene,chain,region37,region38,f)
# mixed_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
}

options(width=200)
library(pQTLtools)
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)
invisible(lapply(c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer"), require, character.only = TRUE))
gwasvcf::set_bcftools("/rds/user/jhz22/hpc-work/bin/bcftools")
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_imported.tsv")
imported_tabix_paths <- read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
INF <- Sys.getenv("INF")
M <- 1e6
sentinels <- subset(read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans")),cis)

# uncomment for all runs inside R
# for (r in 1:nrow(sentinels))
{
  r <- as.integer(Sys.getenv("r"))
  single_run(r)
}
