liftRegion <- function(x,chain,flanking=1e6)
{
  require(GenomicRanges)
  gr <- with(x,GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end))+flanking)
  seqlevelsStyle(gr) <- "UCSC"
  gr38 <- rtracklayer::liftOver(gr, chain)
  chr <- gsub("chr","",colnames(table(seqnames(gr38))))
  start <- min(unlist(start(gr38)))
  end <- max(unlist(end(gr38)))
  invisible(list(chr=chr,start=start,end=end,region=paste0(chr,":",start,"-",end)))
}

sumstats <- function(prot,chr,region37)
{
  cat("GWAS sumstats\n")
  vcf <- file.path(INF,"METAL/gwas2vcf",paste0(prot,".vcf.gz"))
  gwas_stats <- gwasvcf::query_gwas(vcf, chrompos = region37)
  gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats) %>%
                keepSeqlevels(chr) %>%
                renameSeqlevels(paste0("chr",chr))
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
}

microarray <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("a. eQTL datasets\n")
  microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>%
                   dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
  column_names <- names(read.delim(hdr))
  summary_list <- purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region38,
                             selected_gene_id = ensGene, column_names))
  purrr::map_df(summary_list[lapply(summary_list,nrow)!=0],
                ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

rnaseq <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("b. Uniformly processed RNA-seq datasets\n")
  rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>%
               dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list, ~safe_import(., region38,
                             selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list[lapply(result_list,nrow)!=0], ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

gtex <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("c. GTEx_v8 imported eQTL datasets\n")
  rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
               dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list,
                             ~safe_import(., region38, selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list, ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0],
                                ~dplyr::filter(., !is.na(se)))
  purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

coloc <- function(prot,chr,ensGene,chain,region37,region38,out,run_all=FALSE)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37)
  if (run_all)
  {
    df_microarray <- microarray(gwas_stats_hg38,ensGene,region38)
    df_rnaseq <- rnaseq(gwas_stats_hg38,ensGene,region38)
    df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
    if (exists("df_microarray") & exits("df_rnaseq") & exists("df_gtex"))
    {
      coloc_df = dplyr::bind_rows(df_microarray, df_rnaseq, df_gtex)
      saveRDS(coloc_df, file=paste0(out,".RDS"))
      dplyr::arrange(coloc_df, -PP.H4.abf)
      p <- ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
    }
  } else {
    df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
    if (exists("df_gtex"))
    {
      saveRDS(df_gtex,file=paste0(out,".RDS"))
#     dplyr::arrange(df_gtex, -PP.H4.abf)
#     p <- ggplot(df_gtex, aes(x = PP.H4.abf)) + geom_histogram()
    }
  }
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
# ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
#        height = 15, width = 15, units = "cm", dpi = 300)
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
  ensRegion37 <- with(subset(inf1,prot==sentinel[["prot"]]),
                      {
                        start <- start-M
                        if (start<0) start <- 0
                        end <- end+M
                        paste0(chr,":",start,"-",end)
                      })
  region38 <- with(liftRegion(isnpid,chain),region)
  ensGene <- subset(inf1,prot==sentinel[["prot"]])[["ensembl_gene_id"]]
  ensRegion38 <- with(liftRegion(subset(inf1,prot==sentinel[["prot"]]),chain),region)
  f <- file.path(INF,"coloc",with(sentinel,paste0(prot,"-",SNP)))
  cat(chr,region37,region38,ensGene,ensRegion37,ensRegion38,"\n")
# coloc(sentinel[["prot"]],chr,ensGene,chain,region37,region38,f)
  coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
}

# slow with the following loop:
loop <- function() for (r in 1:nrow(sentinels)) single_run(r)

library(pQTLtools)
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)
pkgs <- c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer")
invisible(lapply(pkgs, require, character.only = TRUE))
HPC_WORK <- Sys.getenv("HPC_WORK")
gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
HOME <- Sys.getenv("HOME")
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_imported.tsv")
imported_tabix_paths <- within(read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
      {ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx_V8/ge",
                        paste0(HOME,"/rds/public_databases/GTEx/csv"),ftp_path)})
options(width=200)
library(dplyr)
INF <- Sys.getenv("INF")
M <- 1e6
sentinels <- subset(read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans")),cis)
cvt_rsid <- file.path(INF,"work","INF1.merge.cis.vs.trans-rsid")
prot_rsid <- subset(read.delim(cvt_rsid,sep=" "),cis,select=c(prot,SNP))
# loop()
# Faster with parallel Bash runs.
r <- as.integer(Sys.getenv("r"))
single_run(r)

# with the same setup, to collect results when all single runs are done
collect <- function()
{
  df_coloc <- data.frame()
  for(r in 1:nrow(sentinels))
  {
    prot <- sentinels[["prot"]][r]
    snpid <- sentinels[["SNP"]][r]
    rsid <- prot_rsid[["SNP"]][r]
    f <- file.path(INF,"coloc",paste0(prot,"-",snpid,".RDS"))
    if (!file.exists(f)) next
    cat(prot,"-",rsid,"\n")
    rds <- readRDS(f)
    if (nrow(rds)==0) next
    df_coloc <- rbind(df_coloc,data.frame(prot=prot,rsid=rsid,snpid=snpid,rds))
  }
  df_coloc <- within(df_coloc,{qtl_id <- gsub("GTEx_V8_","",qtl_id)}) %>%
              rename(H0=PP.H0.abf,H1=PP.H1.abf,H2=PP.H2.abf,H3=PP.H3.abf,H4=PP.H4.abf)
  write.table(subset(df_coloc,H3+H4>=0.9&H4/H3>=3),file=file.path(INF,"coloc","GTEx.tsv"),
              quote=FALSE,row.names=FALSE,sep="\t")
  write.table(df_coloc,file=file.path(INF,"coloc","GTEx-all.tsv"),
              quote=FALSE,row.names=FALSE,sep="\t")
}
