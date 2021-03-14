go <- function(i)
{
  id <- subset(tb[[i]],HyperFdrQ<1e-5)
  ord <- order(with(id,HyperFdrQ))
# Hyper_Adjp_BH without download_by option
# Augmented Exploration of the GO with Interactive Simulations
  id[ord,]
}

INF <- Sys.getenv("INF")
INF1_merge <- read.table(file.path(INF,"work","INF1.merge"),as.is=TRUE,header=TRUE)
cis <- scan(file.path(INF,"work","INF1.merge.cis"),what="")
INF1_merge_cis <- subset(INF1_merge,MarkerName%in%cis)
regions <- INF1_merge_cis[c("Chrom","Start","End","POS")]
flank <- 5e+5
regions["Start"] <- regions["POS"] - flank
regions["End"] <- regions["POS"] + flank
reset <- with(regions,Start < 0)
regions[reset,"Start"] <- 0
library(rGREAT)
job = submitGreatJob(regions, species="hg19", version="3.0.0")
tb = getEnrichmentTables(job,download_by = 'tsv')
class(tb)
names(tb)
go1 <- go(1)
go2 <- go(2)
go3 <- go(3)
write.table(rbind(go1,go2,go3),file=file.path(INF,"GREAT","great.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
png(file.path(INF,"GREAT","rGREAT.png"), res=300, units="cm", width=30, height=20)
plotRegionGeneAssociationGraphs(job,type=c(1,3))
dev.off()
availableOntologies(job)
pdf(file.path(INF,"GREAT","rGREAT-GO-top.pdf"),width=12, height=8)
par(mfcol=c(3,1))
plotRegionGeneAssociationGraphs(job, ontology="GO Molecular Function", termID="GO:0005126", type=c(1,3))
plotRegionGeneAssociationGraphs(job, ontology="GO Biological Process", termID="GO:0009611", type=c(1,3))
plotRegionGeneAssociationGraphs(job, ontology="GO Cellular Component", termID="GO:0005615", type=c(1,3))
dev.off()
