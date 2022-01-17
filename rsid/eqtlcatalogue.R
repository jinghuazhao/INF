fetch_eqtl <- function(id="Alasoo_2018.macrophage_naive",ensembl_gene_id,chr=1,start=158175353,end=160525679,flank=1e6)
{
 catalogueR::fetch_restAPI(
             unique_id=id,
             quant_method = "ge",
             infer_region = TRUE,
             gwas_data = NULL,
             chrom = chr,
             bp_lower = start-flank,
             bp_upper = end+flank,
             is_gwas = FALSE,
             size = NULL,
             verbose = TRUE
           ) %>%
           filter(ensembl_gene_id==gene_id.QTL)
}

test <- function()
{
  for (i in 1:nrow(hgnc))
  {
    with(hgnc[i,], {
      cat(gene,":")
      for (j in 1:nrow(mmeta))
      {
        id <- with(mmeta[j,], unique_id)
        cat(" ", id)
        f <- paste0(prot,"-",id)
        t <- fetch_eqtl(id,ensembl_gene_id,chr,start,end)
        save(t,file=file.path(INF,"eQTLCatalogue",paste0(f,".rda")))
      }
      cat("\n")
    })
  }
}

suppressMessages(library(dplyr))
suppressMessages(library(catalogueR))

INF <- Sys.getenv("INF")
HOME <- Sys.getenv("HOME")
hgnc <- read.table(file.path(INF,"work","INF1.merge-cis.genes"),col.names=c("prot","uniprot","chr","start","end","cis.trans")) %>%
        left_join(pQTLtools::inf1) %>%
        select(-cis.trans) %>%
        arrange(chr,start,end)
if(!dir.exists(file.path(INF,"eQTLCatalogue"))) dir.create(file.path(INF,"eQTLCatalogue"))
mmeta <- subset(meta, study!="GTEx" & quant_method=="ge") %>%
         select(unique_id,tissue_ontology_id,sample_size,ftp_path)
save(hgnc,mmeta,INF,file=file.path(INF,"eQTLCatalogue","setup.rda"))
write.table(hgnc,file=file.path(INF,"eQTLCatalogue","hgnc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(mmeta,file=file.path(HOME,"pQTLtools/inst/eQTL-Catalogue","tabix_ftp_paths_ge.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
