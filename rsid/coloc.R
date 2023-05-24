liftRegion <- function(x,chain,flanking=1e6)
{
  require(GenomicRanges)
  gr <- with(x,GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end))+flanking)
  seqlevelsStyle(gr) <- "UCSC"
  gr38 <- rtracklayer::liftOver(gr, chain)
  chr <- gsub("chr","",colnames(table(seqnames(gr38))))
  start <- min(unlist(start(gr38)))
  end <- max(unlist(end(gr38)))
  invisible(list(chr=chr[1],start=start,end=end,region=paste0(chr[1],":",start,"-",end)))
}

sumstats <- function(prot,chr,region37)
{
  cat("GWAS sumstats\n")
  vcf <- file.path(INF,"METAL/gwas2vcf",paste0(prot,".vcf.gz"))
  gwas_stats <- gwasvcf::query_gwas(vcf, chrompos = region37) %>%
                gwasvcf::vcf_to_granges() %>%
                keepSeqlevels(chr) %>%
                renameSeqlevels(paste0("chr",chr))
  gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
    unlist() %>%
#   renameSeqlevels(chr) %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(chromosome = seqnames,
                     position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
    dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(row_count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(row_count == 1) %>%
    mutate(chromosome=gsub("chr","",chromosome))
}

microarray <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("a. eQTL datasets\n")
  microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>%
                   dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
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
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
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
  fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_gtex.tsv")
  imported_tabix_paths <- within(read.delim(fp, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
        {
          f <- lapply(strsplit(ftp_path,"/csv/|/ge/"),"[",3);
          ftp_path <- paste0("~/rds/public_databases/GTEx/csv/",f)
        })
  imported_tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  gtex_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
             dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(gtex_df$ftp_path), gtex_df$qtl_id)
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
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
  fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_ge.tsv")
  imported_tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
  imported_tabix_paths <- within(read.delim(fp, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
        {
          f <- lapply(strsplit(ftp_path,"/csv/|/ge/"),"[",3)
          ftp_path <- paste0("~/rds/public_databases/eQTLCatalogue/",f)
        })
  ftp_path_list <- setNames(as.list(imported_tabix_paths$ftp_path), imported_tabix_paths$unique_id)
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
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
  df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
  if (!exists("df_gtex")) return
  saveRDS(df_gtex,file=paste0(out,".RDS"))
  dplyr::arrange(df_gtex, -PP.H4.abf)
  p <- ggplot(df_gtex, aes(x = PP.H4.abf)) + geom_histogram()
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
  ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
}

ge_coloc <- function(prot,chr,ensGene,chain,region37,region38,out,run_all=FALSE)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37)
  df_ge <- ge(gwas_stats_hg38,ensGene,region38)
  if (!exists("df_ge")) return
  saveRDS(df_ge,file=paste0(out,".RDS"))
  dplyr::arrange(df_ge, -PP.H4.abf)
  p <- ggplot(df_ge, aes(x = PP.H4.abf)) + geom_histogram()
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
  ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
}

all_coloc <- function(prot,chr,ensGene,chain,region37,region38,out,run_all=FALSE)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37)
  df_microarray <- microarray(gwas_stats_hg38,ensGene,region38)
  df_rnaseq <- rnaseq(gwas_stats_hg38,ensGene,region38)
  df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
  df_ge <- ge(gwas_stats_hg38,ensGene,region38)
  if (exists("df_microarray") & exits("df_rnaseq") & exists("df_gtex") & exists("df_ge"))
  {
    coloc_df = dplyr::bind_rows(df_microarray, df_rnaseq, df_gtex, df_ge)
    saveRDS(coloc_df, file=paste0(out,".RDS"))
    dplyr::arrange(coloc_df, -PP.H4.abf)
    p <- ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
  }
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
  ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
}

single_run <- function(r, batch="GTEx")
{
  sentinel <- sentinels[r,]
  chr <- with(sentinel,Chr)
  ss <- subset(inf1,prot==sentinel[["prot"]])
  ensRegion37 <- with(ss,
                      {
                        start <- start-M
                        if (start<0) start <- 0
                        end <- end+M
                        paste0(chr,":",start,"-",end)
                      })
  ensGene <- ss[["ensembl_gene_id"]]
  ensRegion38 <- with(liftRegion(ss,chain),region)
  if (!is.na(ensGene)) ensRegion38 <- with(ss,paste0(chr,":",start38-M,"-",end38+M))
  cat(chr,ensGene,ensRegion37,ensRegion38,"\n")
  if (batch=="GTEx")
  {
    f <- file.path(INF,"coloc",with(sentinel,paste0(prot,"-",SNP)))
    gtex_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
  } else {
    f <- file.path(INF,"eQTLCatalogue",with(sentinel,paste0(prot,"-",SNP)))
    ge_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
  }
}

collect <- function(coloc_dir="coloc")
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
    write.table(subset(df,+H4>=0.8),file=file.path(INF,coloc_dir,"GTEx.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    write.table(df,file=file.path(INF,coloc_dir,"GTEx-all.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
  } else {
    write.table(subset(df,H4>=0.8),file=file.path(INF,coloc_dir,"eQTLCatalogue.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    write.table(df,file=file.path(INF,coloc_dir,"eQTLCatalogue-all.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
  }
}

loop_slowly <- function() for (r in 1:nrow(sentinels)) single_run(r)

# Environmental variables

pkgs <- c("dplyr", "gap", "ggplot2", "readr", "coloc", "GenomicRanges","pQTLtools","seqminer")
invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))

options(width=200)
HOME <- Sys.getenv("HOME")
HPC_WORK <- Sys.getenv("HPC_WORK")
INF <- Sys.getenv("INF")
M <- 1e6

sevens <- "
ENSG00000131142 - CCL25 19 8052318 8062660
ENSG00000125735 - TNFSF14 19 6661253 6670588
ENSG00000275302 - CCL4 17 36103827 36105621
ENSG00000274736 - CCL23 17 36013056 36017972
ENSG00000013725 - CD6 11 60971680 61020377
ENSG00000138675 - FGF5 4 80266639 80336680
ENSG00000277632 - CCL3 17 36088256 36090169
"
updates <- as.data.frame(scan(file=textConnection(sevens),what=list("","","",0,0,0))) %>%
           setNames(c("ensGene","dash","gene","chromosome","start38","end38"))
inf1 <- left_join(pQTLdata::inf1,updates) %>%
        mutate(ensembl_gene_id=if_else(!is.na(ensGene),ensGene,ensembl_gene_id))
sentinels <- subset(read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans")),cis)
cvt_rsid <- file.path(INF,"work","INF1.merge.cis.vs.trans-rsid")
prot_rsid <- subset(read.delim(cvt_rsid,sep=" "),cis,select=c(prot,SNP))

f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
chain <- rtracklayer::import.chain(f)
gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# r <- as.integer(Sys.getenv("r"))
# single_run(r)
# single_run(r,batch="eQTLCatalogue")

collect()
collect(coloc_dir="eQTLCatalogue/ensGene")
