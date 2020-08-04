# 4-8-2020 JHZ

library(regioneR)
INF1_merge <- read.delim("work/INF1.merge")[c("Chrom","Start","End","prot","MarkerName")]
singletons <- with(INF1_merge,End-Start==1)
INF1_merge[singletons,"Start"] <- INF1_merge[singletons,"Start"] - 1e6
INF1_merge[singletons,"End"] <- INF1_merge[singletons,"End"] + 1e6
small <- with(INF1_merge, Start<0)
INF1_merge[small,"Start"] <- 0
cvt <- read.table("work/INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
INF1_merge_cvt <- merge(INF1_merge,cvt,by.x=c("prot","MarkerName"),by.y=c("prot","SNP"))
ord <- with(INF1_merge_cvt,order(Chr,bp))
INF1_merge_cvt <- INF1_merge_cvt[ord,]

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
bm <- with(INF1_merge_cvt, getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
           filters="chromosomal_region", values=paste0(Chr,":",bp,":",bp), mart=mart))
genes <- with(bm,toGRanges(chromosome_name,start_position,end_position,labels=hgnc_symbol))

library(karyoploteR)

png("INF1.merge.png",width=16,height=14,units="cm",res=300)
attach(INF1_merge_cvt)
sentinels <- toGRanges(Chr,bp-1,bp,labels=p.gene)
cis.regions <- toGRanges(Chr,cis.start,cis.end)
loci <- toGRanges(Chr,Start,End)
panel <- toGRanges(p.chr,p.start,p.end,labels=p.gene)
colors <- c("red","blue")
seqlevelsStyle(sentinels) <- "UCSC"
kp <- plotKaryotype(genome="hg19",chromosomes=levels(seqnames(sentinels)))
kpAddBaseNumbers(kp)
kpPlotRegions(kp,data=loci,r0=0.05,r1=0.15,border="black")
kpPlotMarkers(kp, data=sentinels, labels=p.gene, text.orientation="vertical",
              cex=0.45, y=0.3*seq_along(p.gene)/length(p.gene), srt=30, ignore.chromosome.ends=TRUE,
              adjust.label.position=TRUE, label.color=colors[2-cis], label.dist=0.002,
              cex.axis=3, cex.lab=3)
legend("bottomright", bty="n", pch=c(19,19), col=colors, pt.cex=0.4, legend=c("cis", "trans"), text.col=colors, cex=0.8, horiz=FALSE)
kpPlotLinks(kp, data=loci, data2=panel, col=colors[2-cis])
detach(INF1_merge_cvt)
dev.off()
q('no')

kp <- plotKaryotype(genome="hg19",chromosomes=unique(paste0("chr",seqnames(genes))))
kpAddBaseNumbers(kp)
kpPlotRegions(kp,data=genes,r0=0.05,r1=0.15,border="black")
attach(bm)
kpPlotMarkers(kp, data=genes, labels=bm$hgnc_symbol, text.orientation="vertical",
              cex=0.6, y=0.3*seq_along(bm$hgnc_symbol)/length(bm$hgnc_symbol), srt=30, ignore.chromosome.ends=TRUE,
              adjust.label.position=TRUE, label.dist=0.002,
              cex.axis=3, cex.lab=3)
legend("bottomright", pch=c(19,19), col=colors, pt.cex=0.4, legend=c("cis", "trans"), text.col=colors, cex=1, horiz=TRUE)
detach(bm)

init <- function()
{
# https://bernatgel.github.io/karyoploter_tutorial//Examples/PlotGenes/PlotGenes.html
  library(biomaRt)
  library(regioneR)
  gene.symbols <- c("AKT", "APC", "BCR", "BIRC3", "BRAF", "BRCA1", "BRCA2", "CDKN2C", "FEV", "TP53", "PTEN", "RB1")
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                 filters = 'hgnc_symbol', values = gene.symbols, mart = ensembl))
  seqlevelsStyle(genes) <- "UCSC"
  library(karyoploteR) 
  kp <- plotKaryotype(genome="hg38")
  kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "horizontal",
                r1=0.5, cex=0.8, adjust.label.position = FALSE)

# https://www.bioconductor.org/packages/devel/bioc/vignettes/gtrellis/inst/doc/gtrellis.html

  library(ComplexHeatmap)
  library(circlize)
  bed = generateRandomBed(nr = 10000)
  bed = bed[sample(10000, 100), ]
  col_fun = colorRamp2(c(-1, 0, 1), c("green", "yellow", "red"))
  cm = ColorMapping(col_fun = col_fun)
  lgd = color_mapping_legend(cm, plot = FALSE, title = "Value")

  library(gtrellis)
  gtrellis_layout(n_track = 1, ncol = 1, track_axis = FALSE, xpadding = c(0.1, 0),
      gap = unit(4, "mm"), border = FALSE, asist_ticks = FALSE, add_ideogram_track = TRUE, 
      ideogram_track_height = unit(2, "mm"), legend = lgd)
  add_track(bed, panel_fun = function(gr) {
      grid.rect((gr[[2]] + gr[[3]])/2, unit(0.2, "npc"), unit(1, "mm"), unit(0.8, "npc"), 
          hjust = 0, vjust = 0, default.units = "native", 
          gp = gpar(fill = col_fun(gr[[4]]), col = NA))    
  })
  add_track(track = 2, clip = FALSE, panel_fun = function(gr) {
      chr = get_cell_meta_data("name")
      if(chr == "chrY") {
         grid.lines(get_cell_meta_data("xlim"), unit(c(0, 0), "npc"), default.units = "native")
      }
      grid.text(chr, x = 0, y = 0, just = c("left", "bottom"))
  })
}
