# 6-9-2018 JHZ

library(circlize)

xlsx <- "https://github.com/jinghuazhao/INF/blob/master/doc/SOMAscan.xlsx?raw=true"
t <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, cols=c(3,5,7:16,23,24), rows=c(5,1019:1037))
hgTables <- read.delim("hgTables.txt",as.is=TRUE)
hgTables <- within(hgTables, UniProt <- unlist(lapply(strsplit(hgTables$name,"-"),"[",1)))
b <- merge(t,hgTables,by="UniProt")
b1 <- with(b, data.frame(chr=paste0("chr",Chr),start=94844947-1,end=94844947,value1=1,gene="SERPINA1"))
b2 <- with(b, data.frame(chr=X.chrom,start=chromStart,end=chromEnd,value1=t[["Meta-analysis"]]/X14,Target=Target,UniProt=UniProt))
b1
b2
pdf("circlize.pdf")
circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()
circos.genomicLink(b1, b2, col = 10, border = 1, lwd = 15)
circos.clear()
dev.off()
