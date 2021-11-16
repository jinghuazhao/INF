# 16-11-2021 JHZ

export HOTSPOT=$1
export cvt=${INF}/work/INF1.merge.cis.vs.trans

bedtools intersect -a <(grep -w ${HOTSPOT} ${cvt} | awk -v FS="," -v OFS="\t" '{print "chr" $3,$4-1, $4}') \
                   -b <(sort -k1,1n -k2,2n ${INF}/csd3/glist-hg19 | grep -v -w -e X -e Y -e XY | awk '{$1="chr" $1;print}' | sed 's/ /\t/g') \
                   -wa -wb -loj | \
cut -f1-3,7 > hotspot.txt

Rscript -e '
  library(dplyr)
  INF <- Sys.getenv("INF")
  HOTSPOT <- Sys.getenv("HOTSPOT")
  b <- read.csv(Sys.getenv("cvt"),as.is=TRUE) %>% filter(SNP==HOTSPOT) %>%
       mutate(p.chr=paste0("chr",p.chr)) %>% rename(chr=p.chr,start=p.start,end=p.end,gene=p.gene,cistrans=cis.trans)
  cols <- rep(12,nrow(b))
  cols[b[["cis"]]] <- 10
  d <- read.table("hotspot.txt",as.is=TRUE,col.names=c("chr","start","end", "gene"))
  u <- unique(d)
  a <- aggregate(u,by=list(with(u,chr),with(u,start),with(u,end)),FUN="paste",collapse=" ") %>% select(-c(chr,start,end))
  names(a) <- c("chr","start","end", "gene")
  labels <- bind_rows(select(b,chr,start,end,gene),a)
  suppressMessages(library(circlize))
  pdf(file.path(INF,"hotspots",paste0("hotspot-",HOTSPOT,".pdf")))
  circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
  circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  circos.genomicLabels(labels,labels.column = 4, side="inside")
  circos.genomicLink(d, b[c("chr","start","end")], col = cols, directional=1, border = 10, lwd = 2)
  circos.clear()
  dev.off()
  unlink("hotspot.txt")
'

pdftopng -r 300 ${INF}/hotspots/hotspot-${HOTSPOT}.pdf ${HOTSPOT}
mv ${HOTSPOT}-000001.png ${INF}/hotspots/hotspot-${HOTSPOT}.png

# chr12:11058117_C_T
# chr17:34326215_A_C
# chr19:49206145_C_G
