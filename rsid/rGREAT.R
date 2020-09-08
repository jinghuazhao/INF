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
lapply(tb,head,15)
png("work/rGREAT-1.png", res=300, units="cm", width=30, height=30)
plotRegionGeneAssociationGraphs(job,type=1)
dev.off()
png("work/rGREAT-3.png", res=300, units="cm", width=30, height=30)
plotRegionGeneAssociationGraphs(job,type=3)
dev.off()
availableOntologies(job)
pdf("work/rGREAT-GO-0005126.pdf")
plotRegionGeneAssociationGraphs(job, ontology="GO Molecular Function", termID="GO:0005126", type=1)
plotRegionGeneAssociationGraphs(job, ontology="GO Molecular Function", termID="GO:0005126", type=3)
dev.off()
