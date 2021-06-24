post <- function(regions)
{
  library(rGREAT)
  job = submitGreatJob(get(regions), species="hg19", version="3.0.0")
  et = getEnrichmentTables(job,download_by = 'tsv')
  tb <- do.call('rbind',et)
  write.table(tb,file=file.path(INF,"GREAT",paste0(regions,".tsv")),
              quote=FALSE,row.names=FALSE,sep="\t")
  invisible(list(job=job,tb=tb))
}

all_regions <- function()
{
  cistrans.post <- post("cistrans")
  job <- with(cistrans.post,job)
  png(file.path(INF,"GREAT","cistrans.png"), res=300, units="cm", width=30, height=20)
  plotRegionGeneAssociationGraphs(job,type=c(1,3))
  dev.off()
  availableOntologies(job)
  pdf(file.path(INF,"GREAT","cistrans.pdf"),width=12, height=8)
    par(mfcol=c(3,1))
    plotRegionGeneAssociationGraphs(job, ontology="GO Molecular Function", termID="GO:0005126", type=c(1,3))
    plotRegionGeneAssociationGraphs(job, ontology="GO Biological Process", termID="GO:0009611", type=c(1,3))
    plotRegionGeneAssociationGraphs(job, ontology="GO Cellular Component", termID="GO:0005615", type=c(1,3))
  dev.off()
}

three_regions <- function()
{
  tb_all <- data.frame()
  for (r in c("IL12B","KITLG","TNFSF10"))
  {
    r.post <- post(r)
    tb_all <- rbind(tb_all,data.frame(gene=r,with(r.post,tb)))
  }
  write.table(tb_all,file=file.path(INF,"GREAT",paste0("IL12B-KITLG-TNFSF10.tsv")),
              quote=FALSE,row.names=FALSE,sep="\t")
}

library(dplyr)
M <- 1e+6
INF <- Sys.getenv("INF")
INF1_merge <- read.delim(file.path(INF,"work","INF1.merge")) %>%
              mutate(chr=Chrom, start=POS-M, end=POS+M) %>%
              mutate(start=if_else(start<0,0,start)) %>%
              select(prot,MarkerName,chr,start,end)
cistrans <- INF1_merge %>% select(chr,start,end) %>% arrange(chr,start,end) %>% distinct()

# Specific regions
IL12B <- filter(INF1_merge,prot=="IL.12B") %>% select(chr,start,end)
KITLG <- filter(INF1_merge,prot=="SCF") %>% select(chr,start,end)
TNFSF10 <- filter(INF1_merge,prot=="TRAIL") %>% select(chr,start,end)

all_regions()
