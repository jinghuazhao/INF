setup <- function()
{
  hgnc <- read.table(file.path(INF,"work","INF1.merge-cis.genes"),col.names=c("prot","uniprot","chr","start","end","cis.trans")) %>%
          left_join(pQTLtools::inf1) %>%
          select(-cis.trans) %>%
          arrange(chr,start,end)
  if(!dir.exists(file.path(INF,"eQTLCatalogue"))) dir.create(file.path(INF,"eQTLCatalogue"))
  data(meta)
  mmeta <- filter(meta, quant_method=="ge")
  save(hgnc,mmeta,INF,file=file.path(INF,"eQTLCatalogue","setup.rda"))
  write.table(hgnc,file=file.path(INF,"eQTLCatalogue","hgnc.tsv"),sep="\t",quote=FALSE)
  write.table(mmeta,file=file.path(INF,"eQTLCatalogue","paths.tsv"),sep="\t",quote=FALSE)
}

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
setup()

(
  sort -k6,6n -k2,2n
) | \
bedtools merge | \
wc -l

function region()
(
  grep -f work/INF1.cis work/INF1.merge | \
  awk -vd=1e6 -v OFS='\t' '
  {
    if($3-$2<=2) {$2=$2-d;$3=$3+d}
    if ($2<0) $2=0
    print $1,$2,$3,$5,$6,$8,$9
  }' | \
  sort -k6,6n -k2,2n
)

region | \
bedtools merge | \
wc -l

region | \

curl -X GET http://www.ebi.ac.uk/eqtl/api/chromosomes/4/associations?paginate=False&study=Alasoo_2018&qtl_group=macrophage_naive&quant_method=ge&bp_lower=14737349&bp_upper=16737284
