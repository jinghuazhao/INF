# 6-9-2018 JHZ

# title("Missense variant rs28929474:T in SERPINA1 is a trans pQTL hotspot", line = -33, cex = 0.4)
xlsx <- "https://github.com/jinghuazhao/INF/blob/master/doc/SOMAscan.xlsx?raw=true"
t <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, cols=c(3,5,7:16,23,24), rows=c(5,1019:1037))
hgTables <- read.delim("hgTables.txt",as.is=TRUE)
hgTables <- within(hgTables, UniProt <- unlist(lapply(strsplit(hgTables$name,"-"),"[",1)))
b <- merge(t,hgTables,by="UniProt")
b1 <- with(b, data.frame(chr=paste0("chr",Chr),start=94844947-1,end=94844947,value1=1,gene="SERPINA1"))
b2 <- with(b, data.frame(chr=X.chrom,start=chromStart,end=chromEnd,value1=t[["Meta-analysis"]]/X14,gene=geneName,Target=Target,UniProt=UniProt))
b1
b2
SERPINA1 <- data.frame(chr="chr14", start=94843083, end=94857029, value1=1, gene="SERPINA1", Target="P01009-1", UniProt="P01009")
b12 <- rbind(b2,SERPINA1)
pdf("SERPINA1.pdf")
library(circlize)
circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
circos.genomicLabels(b12,labels.column = 5, side="inside")
circos.genomicLink(b1, b2, col = 10, border = 10, lwd = 2)
circos.clear()
dev.off()
