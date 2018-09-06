# 6-9-2018 JHZ

library(circlize)
# toy example
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 1), ]
bed1 <- rbind(bed1,bed1,bed1,bed1,bed1,bed1)
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 6), ]
circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = sample(1:6, 6, replace = TRUE), border = NA)
circos.clear()

# SOMAscan case, the b2 below should not be wxclusively chr14
xlsx <- "https://github.com/jinghuazhao/INF/blob/master/doc/SOMAscan.xlsx?raw=true"
tabs <- "ST4 - pQTL summary"
t <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, cols=c(3,5,7:16,23,24), rows=c(5,1019:1037))
hgTables <- read.delim("hgTables.txt",as.is=TRUE)
hgTables <- within(hgTables, UniProt <- unlist(lapply(strsplit(hgTables$name,"-"),"[",1)))
b <- merge(t,hgTables,by="UniProt")
b1 <- with(b, data.frame(chr=paste0("chr",Chr),start=Pos-1,end=Pos,value1=1))
b2 <- with(b, data.frame(chr=X.chrom,start=chromStart,end=chromEnd,value2=t["Meta-analysis"]/X14,Target=Target,UniProt=UniProt))
circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()
circos.genomicLink(b1, b2, col = sample(1:6, 19, replace = TRUE), border = NA)
circos.clear()

