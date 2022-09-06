# Revised version for circlize plot

setup <- function()
{
  annotate <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE) %>%
              select(-cis,-cis.end,-cis.start,-p.prot,-p.target) %>%
              rename(MarkerName=SNP) %>%
              mutate(chr=Chr,chrom=paste0("hs",Chr),start=bp,end=bp,p.chrom=paste0("hs",p.chr),value=if_else(-log10p>150,150,-log10p),
                     fcolor=ifelse(cis.trans=="cis","color=vdred","color=vdblue"),
                     lcolor=ifelse(cis.trans=="cis","color=lred","color=lblue"),
                     chrbp=paste(Chr,bp,sep=":"),gene=p.gene)
  require(openxlsx)
  f <- file.path(INF,"NG","trans-pQTL_annotation.xlsx")
  annotation_trans <- read.xlsx(f,sheet=1,startRow = 6,colNames=FALSE,cols=c(1:3,14),skipEmptyRows=TRUE) %>%
                      setNames(c("No","rsid","gene.encoding","gene.causal"))
  metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
           select(MarkerName,rsid,prot,Chromosome,Position,log.P.,uniprot,cis.trans) %>%
           rename(log10p=log.P.)
  inf1 <- rename(gap.datasets::inf1,p.target.short=target.short,p.gene=gene,p.chr=chr,p.start=start,p.end=end)
  annotate <- left_join(metal,annotation_trans) %>%
              left_join(inf1) %>%
              mutate(gene=if_else(cis.trans=="cis",p.gene,gene.encoding)) %>%
              rename(Chr=Chromosome,bp=Position) %>%
              arrange(Chr,bp) %>%
              mutate(chr=Chr,chrom=paste0("hs",Chr),start=bp,end=bp,chrbp=paste0(Chr,bp,sep=":"),
                     p.chrom=paste0("hs",p.chr),value=if_else(-log10p>150,150,-log10p),
                     fcolor=ifelse(cis.trans=="cis","color=vdred","color=vdblue"),
                     lcolor=ifelse(cis.trans=="cis","color=lred","color=lblue"))
  subset(annotate,prot=="IL.12B")
  subset(annotate,prot=="TRAIL")
  write.table(annotate,file=file.path(INF,"circos","annotate.txt"),row.names=FALSE)
  f <- file.path(INF,"work","INF1.merge")
  INF1_merge <- read.delim(f)[c("Chrom","Start","End","prot","MarkerName")]
  INF1_merge_cvt <- merge(INF1_merge,annotate,by=c("prot","MarkerName")) %>%
                    left_join(gap.datasets::inf1[c("prot","target.short")]) %>%
                    arrange(Chr,bp)
  pQTLs <- select(INF1_merge_cvt,chrom,start,end,value,fcolor)
  write.table(pQTLs,file=file.path(INF,"circos","pQTLs.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
  label1 <- filter(annotate,cis.trans=="cis") %>%
            mutate(chr=paste0("hs",Chr),start=bp,end=bp,
                   value1=gene,
                   value2=paste0("color=vd","red")) %>%
            select(chr,start,end,value1,value2)
  label2 <- filter(annotate,cis.trans=="trans") %>%
            select(rsid,Chr,bp,gene.encoding,gene.causal) %>%
            distinct() %>%
            mutate(chr=paste0("hs",Chr),start=bp,end=bp,gene.encoding=gsub("; |, ",";",gene.encoding),gene.causal=gsub("; |, ",";",gene.causal),
                   value1=if_else(gene.causal=="-",gene.encoding,paste0(gene.encoding," [",gene.causal,"]")),
                   value2=paste0("color=vd","blue")) %>%
            select(chr,start,end,value1,value2)
  label <- bind_rows(label1,label2) %>%
           arrange(chr,start)
  write.table(label,file=file.path(INF,"circos","pQTL_labels.txt"),
              col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  pQTL_links <- filter(INF1_merge_cvt,cis.trans=="trans") %>%
                select(p.chrom,p.start,p.end,chrom,start,end,lcolor)
  write.table(pQTL_links,file=file.path(INF,"circos","pQTL_links.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
  pQTL_source_genes <- filter(INF1_merge_cvt,cis.trans=="cis") %>%
                       select(p.chr,p.chrom,,p.start,p.end,p.gene) %>%
                       rename(chr=p.chr,chrom=p.chrom,start=p.start,end=p.end,gene=p.gene) %>%
                       distinct()
  pQTL_target_genes <- filter(INF1_merge_cvt,cis.trans=="trans") %>%
                       select(chr,chrom,start,end,gene) %>%
                       mutate(gene=gsub("; ",";",gene))
  pQTL_genes <- rbind(pQTL_source_genes,pQTL_target_genes) %>%
                filter(gene!=".") %>%
                arrange(chr,start) %>%
                select(-chr) %>%
                distinct()
  write.table(pQTL_genes,file=file.path(INF,"circos","pQTL_genes.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
}

circlize <- function()
# circlize
{
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(gridBase))
  suppressMessages(library("circlize"))
  setEPS()
  postscript(file=file.path(INF,"circos","circlize.eps"), width=7.08, height=8.07, horizontal=FALSE, paper="special", colormodel="rgb")
  col_fun <- colorRamp2(c(-1, 1), c("red", "blue"))
  circle_size <- unit(1, "snpc")
  pushViewport(viewport(x=0.5, y=1, width=circle_size, height=circle_size, just=c("center", "bottom")))
  pQTLs <- read.table(file.path(INF,"circos","pQTLs.txt"),col.names=c("chr","start","end","value1","value2")) %>%
           mutate(chr=gsub("hs","chr",chr),value2=gsub("color=vd","",value2))
  pQTL_labels <- read.table(file.path(INF,"circos","pQTL_labels.txt"),col.names=c("chr","start","end","value1","value2"),sep="\t") %>%
           mutate(chr=gsub("hs","chr",chr),value2=gsub("color=vd","",value2))
  pQTL_links <- read.table(file.path(INF,"circos","pQTL_links.txt"),col.names=c("chr1","start1","end1","chr2","start2","end2","color")) %>%
                mutate(chr1=gsub("hs","chr",chr1),chr2=gsub("hs","chr",chr2),color=gsub("color=l","",color))
  pQTL_genes <- read.table(file.path(INF,"circos","pQTL_genes.txt"),col.names=c("chr","start","end","gene")) %>%
                mutate(chr=gsub("hs","chr",chr))
  circos.clear()
  circos.par(start.degree=90, gap.degree=c(rep(c(0.7), 21), 8), track.margin=c(0.005, 0.005), cell.padding=c(0.001, 0.01, 0.01, 0.001))
  circos.initializeWithIdeogram(cytoband=file.path(INF,"circos","cytoband.txt"),plotType=NULL)
  circos.genomicLabels(pQTL_labels, labels.column=4, side="outside", cex=0.45, line_lwd=0.8,
                       labels_height=min(c(cm_h(3.15), max(strwidth(pQTL_labels, cex = 0.4, font = par("font"))))),
                       connection_height=convert_height(8, "mm"), col=pQTL_labels[[5]], line_col=pQTL_labels[[5]])
  circos.track(ylim=c(0,1), track.height=0.05, bg.border=NA,
               panel.fun=function(x, y) {
                 chr=gsub("chr", CELL_META$sector.index, replace="")
                 xlim=CELL_META$xlim
                 ylim=CELL_META$ylim
                 circos.rect(xlim[1], 0, xlim[2], 1, col="white", cex=0.2, lwd=0.5)
                 circos.text(mean(xlim), mean(ylim), chr, cex=0.4, col="black", facing="inside", niceFacing=TRUE)
               })
  circos.track(ylim=c(0,1), track.height=0.03, bg.border=NA, panel.fun=function(x, y) {
               chr=gsub("chr", CELL_META$sector.index, replace="")
               xlim=CELL_META$xlim
               ylim=CELL_META$ylim
               circos.genomicAxis(h="top", direction="inside", labels.cex=0.2, major.at=seq(0,1e10,5e7))})
  circos.genomicTrackPlotRegion(pQTLs, track.height=0.15, bg.border=NA, bg.col="#FFFFFF", ylim=c(0, 150),
                                panel.fun=function(region, value, ...) circos.genomicPoints(region, value, pch=16, col=value[,2], cex=0.3))
  circos.yaxis(side="left", at=seq(0, 150, 50), labels=seq(0, 150, 50), sector.index=get.all.sector.index()[1], labels.cex=0.3, lwd=0.3,
               tick.length=0.5*(convert_x(1, "mm", get.cell.meta.data("sector.index"), get.cell.meta.data("track.index"))))
  circos.genomicText(data.frame(start=1,end=1),sector.index=get.all.sector.index()[1],
                     labels = "-log10(P)",
                     h = "bottom",
                     cex = 0.4,
                     y = 90,
                     adj = c(0.2, 1.5),
                     facing = "clockwise")
  circos.genomicLink(pQTL_links[,1:3], pQTL_links[,4:6], col=pQTL_links[[7]], border=NA, directional=1, arr.length=0.05,
                     arr.width=0.03, arr.lwd=0.05)
  dev.off()
  system("convert -density 300 ${INF}/circos/circlize.eps ${INF}/circos/circlize.png")
}

options(width=200)
INF <- Sys.getenv("INF")
suppressMessages(library(dplyr))

setup()
circlize()
