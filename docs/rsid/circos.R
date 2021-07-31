circos <- function()
{
f <- file.path(INF,"work","INF1.merge")
INF1_merge <- read.delim(f)[c("Chrom","Start","End","prot","MarkerName")]

cvt <- read.csv(file.path(INF,"work","/INF1.merge.cis.vs.trans"),as.is=TRUE) %>%
                rename(MarkerName=SNP) %>% 
                mutate(chr=Chr,chrom=paste0("hs",Chr),start=bp,end=bp,p.chrom=paste0("hs",p.chr),value=-log10p,
                       fcolor=ifelse(cis,"color=vdred","color=vdblue"),
                       lcolor=ifelse(cis,"color=lred","color=lblue"),
                       chrbp=paste(Chr,bp,sep=":"))

geneinfo <- gene <- vector("character")
for(i in 1:180)
{
  geneinfo[i] <- with(ieugwasr::variants_chrpos(cvt$chrbp[i],5000),geneinfo)
  gene[i]=gsub(":([0-9])*","",geneinfo[i])
}
annotate <- within(cvt,{gene=gene})
is.cis <- with(annotate,cis)
annotate[is.cis,"gene"] <- annotate[is.cis,"p.gene"]
INF1_merge_cvt <- merge(INF1_merge,annotate,by=c("prot","MarkerName")) %>%
                  left_join(gap.datasets::inf1[c("prot","target.short")]) %>%
                  arrange(Chr,bp)

pQTLs <- select(INF1_merge_cvt,chrom,start,end,value,fcolor)
write.table(pQTLs,file="pQTLs.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
pQTL_labels <- filter(INF1_merge_cvt,gene!=".") %>%
               mutate(gene=paste0(gene,"/",target.short)) %>%
               select(chrom,start,end,gene,fcolor)
write.table(pQTL_labels,file="pQTL_labels.txt",col.names=FALSE,row.names=FALSE,quote=TRUE)
pQTL_links <- filter(INF1_merge_cvt,!cis) %>%
              select(p.chrom,p.start,p.end,chrom,start,end,lcolor)
write.table(pQTL_links,file="pQTL_links.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

## not used:

pQTL_source_genes <- filter(INF1_merge_cvt,cis) %>%
                     select(p.chr,p.chrom,,p.start,p.end,p.gene) %>%
                     rename(chr=p.chr,chrom=p.chrom,start=p.start,end=p.end,gene=p.gene) %>%
                     distinct()
pQTL_target_genes <- filter(INF1_merge_cvt,!cis) %>%
                     select(chr,chrom,start,end,gene)
pQTL_genes <- rbind(pQTL_source_genes,pQTL_target_genes) %>%
              filter(gene!=".") %>%
              arrange(chr,start) %>%
              select(-chr) %>%
              distinct()
write.table(pQTL_genes,file="pQTL_genes.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

## circlize

circlize <- function()
{
  pQTLs <- read.table("pQTLs.txt",col.names=c("chr","start","end","value1","value2")) %>%
           mutate(chr=gsub("hs","chr",chr),value2=gsub("color=vd","",value2))
  pQTL_labels <- read.table("pQTL_labels.txt",col.names=c("chr","start","end","value1","value2")) %>%
           mutate(chr=gsub("hs","chr",chr),value2=gsub("color=vd","",value2))
  pQTL_links <- read.table("pQTL_links.txt",col.names=c("chr1","start1","end1","chr2","start2","end2","color")) %>%
                mutate(chr1=gsub("hs","chr",chr1),chr2=gsub("hs","chr",chr2),color=gsub("color=l","",color))
  library("circlize")
  setEPS()
  postscript(file = "circlize.eps", width = 7.08, height = 7.08, horizontal = FALSE, paper = "special", colormodel="rgb")
  circos.clear()
  circos.par("start.degree" = 90, gap.degree = c(rep(c(0.7), 21), 8), track.margin = c(0.005, 0.005), cell.padding = c(0.001, 0.01, 0.01, 0.001))
# circos.initializeWithIdeogram(plotType = NULL, species = "hg19", chromosome.index = paste0("chr", 1:22))
  circos.initializeWithIdeogram(cytoband=file.path(INF,"circos","cytoband.txt"),plotType=NULL)
  circos.genomicLabels(pQTL_labels, labels.column = 4, side = "outside", cex = 0.45, line_lwd = 0.8,
                       connection_height = convert_height(8, "mm"), col=pQTL_labels[[5]], line_col=pQTL_labels[[5]])
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
               panel.fun = function(x, y) {
                 chr  = gsub("chr", CELL_META$sector.index, replace = "")
                 xlim = CELL_META$xlim
                 ylim = CELL_META$ylim
                 circos.rect(xlim[1], 0, xlim[2], 1, col = "white", cex = 0.2, lwd = 0.5)
                 circos.text(mean(xlim), mean(ylim), chr, cex = 0.4, col = "black", facing = "inside", niceFacing = TRUE)
               })
  circos.track(ylim=c(0,1), track.height = 0.03, bg.border = NA, panel.fun=function(x, y) {
               chr  = gsub("chr", CELL_META$sector.index, replace = "")
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               circos.genomicAxis(h = "top", direction = "inside", labels.cex=0.2, major.at=seq(0,1e10,5e7))})
  circos.genomicTrackPlotRegion(pQTLs, track.height = 0.25, bg.border = NA, bg.col = "#FFFFFF", ylim = c(0, 150),
                                panel.fun = function(region, value, ...) circos.genomicPoints(region, value, pch = 16, col = value[,2], cex = 0.3))
  circos.yaxis(side = "left", at = seq(0, 150, 50), labels = seq(0, 150, 50), sector.index = get.all.sector.index()[1], labels.cex = 0.3, lwd = 0.3,
               tick.length = 0.5*(convert_x(1, "mm", get.cell.meta.data("sector.index"), get.cell.meta.data("track.index"))))
  circos.genomicLink(pQTL_links[,1:3], pQTL_links[,4:6], col = pQTL_links[[7]], border = NA, directional = 1, arr.length = 0.05,
                     arr.width = 0.03, arr.lwd=0.05)
  dev.off()
}
# https://www.rapidtables.com/web/color/RGB_Color.html

INF <- Sys.getenv("INF")
library(dplyr)

# target.short --> prot and quote=FALSE for circos
circos()
circlize()
