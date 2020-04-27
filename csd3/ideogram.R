# 27-4-2020 JHZ

library(biomaRt)
library(regioneR)
library(karyoploteR)
cvt <- read.table("work/INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
attach(cvt)
genes <- toGRanges(Chr,bp-1,bp,p.gene)
seqlevelsStyle(genes) <- "UCSC"
colors <- c("red","blue")
png("INF1.merge.png",width=12,height=10,units="in",pointsize=4,res=300)
kp <- plotKaryotype(genome="hg19",chromosomes="autosomal")
kpPlotMarkers(kp, data=genes, labels=p.gene, text.orientation = "vertical",
              r1=0.5, cex=0.8, srt=45, adjust.label.position = TRUE, label.color=colors[cis+1])
legend("bottomright", legend=c("cis", "trans"), box.lty=0, text.col=c("red", "blue"), cex=0.8)
dev.off()
detach(cvt)

init <- function()
{
# https://bernatgel.github.io/karyoploter_tutorial//Examples/PlotGenes/PlotGenes.html
  library(biomaRt)
  library(regioneR)
  gene.symbols <- c("AKT", "APC", "BCR", "BIRC3", "BRAF", "BRCA1", "BRCA2", "CDKN2C", "FEV", "TP53", "PTEN", "RB1")
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                 filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
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
