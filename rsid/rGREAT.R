go <- function(i,tb,condition)
{
  id <- subset(tb[[i]],eval(parse(text=condition)))
  col <- unlist(strsplit(condition,"<"))[1]
  ord <- order(with(id,eval(parse(text=col))))
# Hyper_Adjp_BH without download_by option
# Augmented Exploration of the GO with Interactive Simulations
  id[ord,]
}

post <- function(regions,name,toplot,condition)
{
  flank <- 5e+5
  regions["Start"] <- regions["POS"] - flank
  regions["End"] <- regions["POS"] + flank
  reset <- with(regions,Start < 0)
  regions[reset,"Start"] <- 0
  library(rGREAT)
  job = submitGreatJob(regions, species="hg19", version="3.0.0")
  tb = getEnrichmentTables(job,download_by = 'tsv')
  write.table(do.call('rbind',tb),file=file.path(INF,"GREAT",paste0(name,"-all.tsv")),quote=FALSE,row.names=FALSE,sep="\t")
  go1 <- go(1,tb,condition)
  go2 <- go(2,tb,condition)
  go3 <- go(3,tb,condition)
  write.table(rbind(go1,go2,go3),file=file.path(INF,"GREAT",paste0(name,".tsv")),quote=FALSE,row.names=FALSE,sep="\t")
  if (toplot)
  {
    png(file.path(INF,"GREAT",paste0(name,".png")), res=300, units="cm", width=30, height=20)
    plotRegionGeneAssociationGraphs(job,type=c(1,3))
    dev.off()
    availableOntologies(job)
    pdf(file.path(INF,"GREAT",paste0(name,".pdf")),width=12, height=8)
    par(mfcol=c(3,1))
    plotRegionGeneAssociationGraphs(job, ontology="GO Molecular Function", termID="GO:0005126", type=c(1,3))
    plotRegionGeneAssociationGraphs(job, ontology="GO Biological Process", termID="GO:0009611", type=c(1,3))
    plotRegionGeneAssociationGraphs(job, ontology="GO Cellular Component", termID="GO:0005615", type=c(1,3))
    dev.off()
  }
  invisible(list(job=job,tb=tb,go1=go1,go2=go2,go3=go3))
}

INF <- Sys.getenv("INF")
INF1_merge <- read.table(file.path(INF,"work","INF1.merge"),as.is=TRUE,header=TRUE)

cis <- scan(file.path(INF,"work","INF1.merge.cis"),what="")
INF1_merge_cis <- subset(INF1_merge,MarkerName%in%cis)
cis.regions <- INF1_merge_cis[c("Chrom","Start","End","POS")]
IL12B <- subset(INF1_merge,prot=="IL.12B",select=c(Chrom,Start,End,POS))
KITLG <- subset(INF1_merge,prot=="SCF",select=c(Chrom,Start,End,POS))
TNFSF10 <- subset(INF1_merge,prot=="TRAIL",select=c(Chrom,Start,End,POS))

post(cis.regions,"GREAT",TRUE,"HyperFdrQ<=1e-5")
post(IL12B,"IL12B",FALSE,"BinomP<=0.0001")
post(KITLG,"KITLG",FALSE,"BinomP<=0.0001")
post(TNFSF10,"TNFSF10",FALSE,"BinomP<=0.0001")
