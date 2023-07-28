rsid2gene <- function()
# ieugwasr annotation
{
  geneinfo <- vector("character")
  for(i in 1:180) geneinfo[i] <- with(ieugwasr::variants_chrpos(cvt$chrbp[i],300),geneinfo)
  gene <- gsub(":([0-9])*","",geneinfo)
  gene <- stringr::str_split(gene,"[|]",simplify=TRUE)[,1]
}

setup <- function(simplify=TRUE)
{
# nearest <- rsid2gene()
  vep <- read.delim(file.path(INF,"annotate","INF1.merge-annotate.tsv")) %>%
         select(Protein,Location,NEAREST) %>%
         rename(prot=Protein,nearest=NEAREST) %>%
         arrange(Location)
# only VEP annotated genes
  annotate <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE) %>%
              select(-cis,-cis.end,-cis.start,-p.prot,-p.target) %>%
              rename(MarkerName=SNP) %>%
              mutate(Location=paste(Chr,bp,sep=":")) %>%
              left_join(vep) %>%
              mutate(chr=Chr,chrom=paste0("hs",Chr),start=bp,end=bp,p.chrom=paste0("hs",p.chr),value=if_else(-log10p>150,150,-log10p),
                     fcolor=ifelse(cis.trans=="cis","color=vdred","color=vdblue"),
                     lcolor=ifelse(cis.trans=="cis","color=lred","color=lblue"),
                     chrbp=paste(Chr,bp,sep=":"),gene=if_else(cis.trans=="cis",p.gene,nearest))
# causal genes
  require(openxlsx)
  f <- file.path(INF,"NG","trans-pQTL_annotation.xlsx")
  annotation_trans <- read.xlsx(f,sheet=1,startRow = 6,colNames=FALSE,cols=c(1:3,14),skipEmptyRows=TRUE)
  names(annotation_trans) <- c("No","rsid","gene.encoding","gene.causal")
  metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
           select(MarkerName,rsid,prot,Chromosome,Position,log.P.,uniprot,cis.trans) %>%
           rename(log10p=log.P.)
  inf1 <- rename(gap.datasets::inf1,p.target.short=target.short,p.gene=gene,p.chr=chr,p.start=start,p.end=end)
  annotate <- left_join(metal,annotation_trans) %>%
              left_join(inf1) %>%
              mutate(Location=paste0(Chromosome,":",Position)) %>%
              left_join(vep) %>%
              mutate(gene=if_else(cis.trans=="cis",p.gene,gsub("; |, ",";",gene.causal))) %>%
              mutate(gene=if_else(gene=="-",nearest,gene)) %>%
              mutate(gene=if_else(rsid=="rs7612912","ACKR2",gene)) %>%
              rename(Chr=Chromosome,bp=Position,chrbp=Location) %>%
              arrange(Chr,bp) %>%
              mutate(chr=Chr,chrom=paste0("hs",Chr),start=bp,end=bp,
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
  if (!simplify)
  {
    pQTL_labels <- filter(INF1_merge_cvt,gene!=".") %>%
                   select(chrom,start,end,gene,fcolor)
    write.table(pQTL_labels,file=file.path(INF,"circos","pQTL_labels.txt"),col.names=FALSE,row.names=FALSE,quote=TRUE,sep="\t")
  } else {
    collapse <- annotate %>%
                group_by(gene,cis.trans) %>%
                summarise(nprots=n(),chr=paste(chr,collapse=","),pos=paste(bp,collapse=","),
                          chrpos=paste(chrbp,collapse=";"),
                          prots=paste(p.target.short,collapse=";"),
                          log10p=paste(value,collapse=";"),
                          cistrans=paste(cis.trans,collapse=";")) %>%
                mutate(CHR=unlist(lapply(lapply(strsplit(chr,","),as.integer),median)),
                       POS=round(unlist(lapply(lapply(strsplit(pos,","),as.integer),median)))) %>%
                arrange(CHR,POS) %>%
                select(-chr,-pos) %>%
                mutate(chr=paste0("hs",CHR),start=POS,end=POS,
                       value1=paste0(gene," [",prots,"]"),
                       value2=paste0("color=vd",if_else(cistrans=="cis","red","blue")))
    write.table(collapse[c("chr","start","end","value1","value2")],file=file.path(INF,"circos","pQTL_labels.txt"),
                col.names=FALSE,row.names=FALSE,quote=TRUE,sep="\t")
  }
  pQTL_links <- filter(INF1_merge_cvt,cis.trans=="trans") %>%
                select(p.chrom,p.start,p.end,chrom,start,end,lcolor)
  write.table(pQTL_links,file=file.path(INF,"circos","pQTL_links.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
  pQTL_source_genes <- filter(INF1_merge_cvt,cis.trans=="cis") %>%
                       select(p.chr,p.chrom,,p.start,p.end,p.gene) %>%
                       rename(chr=p.chr,chrom=p.chrom,start=p.start,end=p.end,gene=p.gene) %>%
                       distinct()
  pQTL_target_genes <- filter(INF1_merge_cvt,cis.trans=="trans") %>%
                       select(chr,chrom,start,end,gene)
  pQTL_genes <- rbind(pQTL_source_genes,pQTL_target_genes) %>%
                filter(gene!=".") %>%
                arrange(chr,start) %>%
                select(-chr) %>%
                distinct()
  write.table(pQTL_genes,file=file.path(INF,"circos","pQTL_genes.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
}

## circlize

circlize <- function()
{
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(gridBase))
  suppressMessages(library("circlize"))
  setEPS()
  postscript(file=file.path(INF,"circlize-old.ps"), width=7.08, height=7.08, horizontal=FALSE, paper="special", colormodel="rgb")
  col_fun <- colorRamp2(c(-1, 1), c("red", "blue"))
  circle_size <- unit(1, "snpc")
  llabels <- Legend(at=c("cis","trans"), type="points", legend_gp=gpar(col=c("red","blue")), title_position="topleft",
                    title="1. Gene [target proteins]", nrow=1)
  lpoints <- Legend(at=c("cis","trans"), type="points", legend_gp=gpar(col=c("red","blue")), title_position="topleft", 
                    title="2. -log10(P) [ceiling=150]", nrow=1)
  llinks <- Legend(at="trans", type="lines", legend_gp=gpar(col="blue", lwd=2), title_position="topleft",
                    title="3. trans-pQTL connections", nrow=1)
  cis.trans <- Legend(at=c("cis","trans"), type="points", legend_gp=gpar(col=c("red","blue")), title_position="leftcenter",
                      title="cis/trans", nrow=1)
  llist_vertical = packLegend(llabels, lpoints, llinks, direction = "vertical")
  llist_horizontal = packLegend(cis.trans, direction = "horizontal")
  plot.new()
  pushViewport(viewport(x=0.5, y=1, width=circle_size, height=circle_size, just=c("center", "top")))
  pQTLs <- read.table(file.path(INF,"circos","pQTLs.txt"),col.names=c("chr","start","end","value1","value2")) %>%
           mutate(chr=gsub("hs","chr",chr),value2=gsub("color=vd","",value2))
  pQTL_labels <- read.table(file.path(INF,"circos","pQTL_labels.txt"),col.names=c("chr","start","end","value1","value2"),sep="\t") %>%
           mutate(chr=gsub("hs","chr",chr),value2=gsub("color=vd","",value2))
  pQTL_links <- read.table(file.path(INF,"circos","pQTL_links.txt"),col.names=c("chr1","start1","end1","chr2","start2","end2","color")) %>%
                mutate(chr1=gsub("hs","chr",chr1),chr2=gsub("hs","chr",chr2),color=gsub("color=l","",color))
  pQTL_genes <- read.table(file.path(INF,"circos","pQTL_genes.txt"),col.names=c("chr","start","end","gene")) %>%
                mutate(chr=gsub("hs","chr",chr))
  par(omi=gridOMI(), new=TRUE)
  circos.clear()
  circos.par(start.degree=90, gap.degree=c(rep(c(0.7), 21), 8), track.margin=c(0.005, 0.005), cell.padding=c(0.001, 0.01, 0.01, 0.001))
  circos.initializeWithIdeogram(cytoband=file.path(INF,"circos","cytoband.txt"),plotType=NULL)
  circos.genomicLabels(pQTL_labels, labels.column=4, side="outside", cex=0.45, line_lwd=0.8,
                       connection_height=convert_height(8, "mm"), col=pQTL_labels[[5]], line_col=pQTL_labels[[5]])
  circos.track(ylim=c(0, 1), track.height=0.05, bg.border=NA,
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
  circos.genomicTrackPlotRegion(pQTLs, track.height=0.25, bg.border=NA, bg.col="#FFFFFF", ylim=c(0, 150),
                                panel.fun=function(region, value, ...) circos.genomicPoints(region, value, pch=16, col=value[,2], cex=0.3))
  circos.yaxis(side="left", at=seq(0, 150, 50), labels=seq(0, 150, 50), sector.index=get.all.sector.index()[1], labels.cex=0.3, lwd=0.3,
               tick.length=0.5*(convert_x(1, "mm", get.cell.meta.data("sector.index"), get.cell.meta.data("track.index"))))
  circos.genomicLink(pQTL_links[,1:3], pQTL_links[,4:6], col=pQTL_links[[7]], border=NA, directional=1, arr.length=0.05,
                     arr.width=0.03, arr.lwd=0.05)
  upViewport()
  draw(llist_horizontal, x=circle_size*0.96, y=circle_size/12, just="right")
  dev.off()
  system("ps2pdf ${INF}/circlize-old.ps")
  system("convert -density 300 ${INF}/circlize-old.ps ${INF}/circlize-old.png")
}

# https://www.rapidtables.com/web/color/RGB_Color.html

options(width=200)
INF <- Sys.getenv("INF")
suppressMessages(library(dplyr))

# target.short --> prot and quote=FALSE for circos

setup()
circlize()
