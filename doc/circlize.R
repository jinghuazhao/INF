#28-11-2018 JHZ

inf1 <- read.delim("olink.inf.panel.annot.tsv", as.is=TRUE)
inf1[with(inf1, uniprot=="Q8NF90"),"hgnc_symbol"] <- "FGF5"
inf1[with(inf1, uniprot=="Q8WWJ7"),"hgnc_symbol"] <- "CD6"
prot <- read.table("inf1.list",col.names=c("prot","uniprot"),as.is=TRUE,sep="\t")
p <- merge(inf1,prot,by="uniprot")[c("chromosome_name","start_position","end_position","hgnc_symbol","prot","uniprot")]
names(p) <- c("chr","start","end","gene","prot","uniprot")
# INTERVAL
clumped <- read.table("INTERVAL.clumped.dat",as.is=TRUE,header=TRUE)
hits <- merge(clumped[c("CHR","BP","SNP","prot")],p[c("prot","uniprot")],by="prot")
names(hits) <- c("prot","Chr","bp","SNP","uniprot")
require(gap)
cvt <- cis.vs.trans.classification(hits,p)
with(cvt,summary(data))
b1 <- with(cvt,data[c("Chr","bp")])
b1 <- within(b1,{Chr=paste0("chr",Chr);start=bp-1})
names(b1) <- c("chr","end","start")
b2 <- with(cvt,data[c("p.chr","cis.start","cis.end","p.gene","p.prot")])
b2 <- within(b2,{p.chr=paste0("chr",p.chr)})
names(b2) <- c("chr","start","end","gene","prot")
ann <- read.table("st.bed",as.is=TRUE,header=TRUE)
names(ann)[1] <- "chr"
ann <- within(ann, {chr=paste0("chr",chr);start=start-1e6;end <- end+1e6})
ann[with(ann,start<0),"start"] <- 0
library(circlize)
pdf("INTERVAL.circlize.pdf")
circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
circos.genomicLabels(ann,labels.column = 4, side="inside")
circos.genomicLink(b1, b2, col = 10, border = 10, lwd = 2)
circos.clear()
dev.off()

#INF1
clumped <- read.table("INF1.clumped",as.is=TRUE,header=TRUE)
hits <- merge(clumped[c("CHR","BP","SNP","prot")],p[c("prot","uniprot")],by="prot")
names(hits) <- c("prot","Chr","bp","SNP","uniprot")
hits_excl <- c( "ADA", "CCL25", "CD6", "CST5", "FGF.5", "IFN.gamma",
                "IL.13", "IL.18R1", "IL.1.alpha", "IL.20", "IL.20RA", "IL.22.RA1",
                "IL.24", "IL.2RB", "IL.33", "LIF", "MCP.2", "NRTN", "TSLP",
                "IL.10RA", "IL.5", "TNF")
hits_inc <- subset(hits, !(prot%in%hits_excl))
cvt <- cis.vs.trans.classification(hits_inc,p)
with(cvt,summary(data))
b1 <- with(cvt,data[c("Chr","bp")])
b1 <- within(b1,{Chr=paste0("chr",Chr);start=bp-1})
names(b1) <- c("chr","end","start")
b2 <- with(cvt,data[c("p.chr","cis.start","cis.end","p.gene","p.prot")])
b2 <- within(b2,{p.chr=paste0("chr",p.chr)})
names(b2) <- c("chr","start","end","gene","prot")
ann <- read.table("st.bed",as.is=TRUE,header=TRUE)
names(ann)[1] <- "chr"
ann <- within(ann, {chr=paste0("chr",chr);start=start-1e6;end <- end+1e6})
ann[with(ann,start<0),"start"] <- 0
library(circlize)
pdf("INF1.circlize.pdf")
circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
circos.genomicLabels(ann,labels.column = 4, side="inside")
circos.genomicLink(b1, b2, col = 10, border = 10, lwd = 2)
circos.clear()
dev.off()


#> head(cistrans$data)
#  uniprot   prot Chr        bp                SNP p.chr   p.start     p.end
#1  O00300    OPG   8 120081031 chr8:120081031_C_T     8 119935796 119964439
#2  O00300    OPG  17  26612996 chr17:26612996_C_T     8 119935796 119964439
#3  O14625 CXCL11   4  76964956  chr4:76964956_A_G     4  76954835  76962568
#4  O14788 TRANCE   3 172274232 chr3:172274232_A_C    13  43136872  43182149
#5  O14788 TRANCE   8 120201029 chr8:120201029_C_T    13  43136872  43182149
#6  O14788 TRANCE  13  42999096 chr13:42999096_G_T    13  43136872  43182149
#     p.gene p.prot cis.trans   cis   cis.end cis.start
#1 TNFRSF11B    OPG       cis  TRUE 120964439 118935796
#2 TNFRSF11B    OPG     trans FALSE 120964439 118935796
#3    CXCL11 CXCL11       cis  TRUE  77962568  75954835
#4   TNFSF11 TRANCE     trans FALSE  44182149  42136872
#5   TNFSF11 TRANCE     trans FALSE  44182149  42136872
#6   TNFSF11 TRANCE       cis  TRUE  44182149  42136872
