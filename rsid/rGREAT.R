go <- function(i,lines=15)
{
  ord <- order(with(tb[[i]],Hyper_Adjp_BH))
  print(head(tb[[i]][ord1,],lines))
  id <- with(subset(tb[[i]],Hyper_Adjp_BH<1e-8),ID)
# Augmented Exploration of the GO with Interactive Simulations
  write.table(id,file=paste0("id",i,".tsv"),col.names=FALSE,quote=FALSE,row.names=FALSE)
}

INF1_merge <- read.table("work/INF1.merge",as.is=TRUE,header=TRUE)
regions <- INF1_merge[c("Chrom","Start","End")]
singletons <- with(regions, Start-End<=2)
flank <- 5e+2
regions[singletons,"Start"] <- regions[singletons,"Start"] - flank
regions[singletons,"End"] <- regions[singletons,"End"] + flank
reset <- with(regions,Start < 0)
regions[reset,"Start"] <- 0
library(rGREAT)
job = submitGreatJob(regions, species="hg19", version="3.0.0")
tb = getEnrichmentTables(job)
class(tb)
names(tb)
go(1)
go(2)
go(3)
png("work/rGREAT.png", res=300, units="cm", width=30, height=20)
plotRegionGeneAssociationGraphs(job,type=c(1,3))
dev.off()
availableOntologies(job)
pdf("work/rGREAT-GO-top.pdf",width=12, height=8)
par(mfcol=c(3,1))
plotRegionGeneAssociationGraphs(job, ontology="GO Molecular Function", termID="GO:0005126", type=c(1,3))
plotRegionGeneAssociationGraphs(job, ontology="GO Biological Process", termID="GO:0009611", type=c(1,3))
plotRegionGeneAssociationGraphs(job, ontology="GO Cellular Component", termID="GO:0005615", type=c(1,3))
dev.off()
