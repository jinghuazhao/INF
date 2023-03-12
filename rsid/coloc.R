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
    dplyr::transmute(chromosome = seqnames,
                     position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
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
  f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_gtex.tsv")
  imported_tabix_paths <- within(read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
        {ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx_V8/ge",
                          paste0(HOME,"/rds/public_databases/GTEx/csv"),ftp_path)})
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

ge <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("d. eQTL datasets\n")
  f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_ge.tsv")
  imported_tabix_paths <- within(read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
        {ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Alasoo_2018/ge",
                          paste0(HOME,"/rds/public_databases/eQTLCatalogue"),ftp_path)})
  ftp_path_list <- setNames(as.list(imported_tabix_paths$ftp_path), imported_tabix_paths$unique_id)
  hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list,
                             ~safe_import(., region38, selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list, ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0],
                                ~dplyr::filter(., !is.na(se)))
  purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "unique_id")
}

gtex_coloc <- function(prot,chr,ensGene,chain,region37,region38,out,run_all=FALSE)
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

ge_coloc <- function(prot,chr,ensGene,chain,region37,region38,out,run_all=FALSE)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37)
  if (run_all)
  {
    df_microarray <- microarray(gwas_stats_hg38,ensGene,region38)
    df_rnaseq <- rnaseq(gwas_stats_hg38,ensGene,region38)
    df_ge <- ge(gwas_stats_hg38,ensGene,region38)
    if (exists("df_microarray") & exits("df_rnaseq") & exists("df_gtex"))
    {
      coloc_df = dplyr::bind_rows(df_microarray, df_rnaseq, df_gtex)
      saveRDS(coloc_df, file=paste0(out,".RDS"))
      dplyr::arrange(coloc_df, -PP.H4.abf)
      p <- ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
    }
  } else {
    df_ge <- ge(gwas_stats_hg38,ensGene,region38)
    if (exists("df_ge"))
    {
      saveRDS(df_ge,file=paste0(out,".RDS"))
#     dplyr::arrange(df_ge, -PP.H4.abf)
#     p <- ggplot(df_ge, aes(x = PP.H4.abf)) + geom_histogram()
    }
  }
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
# ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
#        height = 15, width = 15, units = "cm", dpi = 300)
}

single_run <- function(r, batch="GTEx")
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
  cat(chr,region37,region38,ensGene,ensRegion37,ensRegion38,"\n")
  if (batch=="GTEx")
  {
    f <- file.path(INF,"coloc",with(sentinel,paste0(prot,"-",SNP)))
    gtex_coloc(sentinel[["prot"]],chr,ensGene,chain,region37,region38,f)
  # gtex_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
  } else {
    f <- file.path(INF,"eQTLCatalogue",with(sentinel,paste0(prot,"-",SNP)))
    ge_coloc(sentinel[["prot"]],chr,ensGene,chain,region37,region38,f)
  # ge_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
  }
}

collect <- function(coloc_dir="eQTLCatalogue")
# to collect results when all single runs are done
{
  df_coloc <- data.frame()
  for(r in 1:nrow(sentinels))
  {
    prot <- sentinels[["prot"]][r]
    snpid <- sentinels[["SNP"]][r]
    rsid <- prot_rsid[["SNP"]][r]
    f <- file.path(INF,coloc_dir,paste0(prot,"-",snpid,".RDS"))
    if (!file.exists(f)) next
    cat(prot,"-",rsid,"\n")
    rds <- readRDS(f)
    if (nrow(rds)==0) next
    df_coloc <- rbind(df_coloc,data.frame(prot=prot,rsid=rsid,snpid=snpid,rds))
  }
  df <- dplyr::rename(df_coloc,H0=PP.H0.abf,H1=PP.H1.abf,H2=PP.H2.abf,H3=PP.H3.abf,H4=PP.H4.abf)
  if (coloc_dir=="coloc") {
    df_coloc <- within(df,{qtl_id <- gsub("GTEx_V8_","",qtl_id)})
    write.table(subset(df,H3+H4>=0.9&H4/H3>=3),file=file.path(INF,coloc_dir,"GTEx.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    write.table(df,file=file.path(INF,coloc_dir,"GTEx-all.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
  } else {
    write.table(subset(df,H3+H4>=0.9&H4/H3>=3),file=file.path(INF,coloc_dir,"eQTLCatalogue.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    write.table(df,file=file.path(INF,coloc_dir,"eQTLCatalogue-all.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
  }
}

loop_slowly <- function() for (r in 1:nrow(sentinels)) single_run(r)

# Environmental variables

options(width=200)
HOME <- Sys.getenv("HOME")
HPC_WORK <- Sys.getenv("HPC_WORK")
INF <- Sys.getenv("INF")
M <- 1e6

pkgs <- c("dplyr", "gap", "ggplot2", "readr", "coloc", "GenomicRanges","pQTLtools","seqminer")
invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))

f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)
gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
sentinels <- subset(read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans")),cis)
cvt_rsid <- file.path(INF,"work","INF1.merge.cis.vs.trans-rsid")
prot_rsid <- subset(read.delim(cvt_rsid,sep=" "),cis,select=c(prot,SNP))

r <- as.integer(Sys.getenv("r"))
single_run(r)
single_run(r,batch="eQTLCatalogue")
