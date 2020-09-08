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
png("work/rGREAT-1.png", res=300, units="cm", width=30, height=30)
res = plotRegionGeneAssociationGraphs(job,type=1)
dev.off()
png("work/rGREAT-3.png", res=300, units="cm", width=30, height=30)
res = plotRegionGeneAssociationGraphs(job,type=3)
dev.off()
# term only in hg38
res = plotRegionGeneAssociationGraphs(job, ontology="GO_Molecular_Function", termID="GO:0004984")
