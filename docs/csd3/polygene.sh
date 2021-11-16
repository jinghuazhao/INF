# 16-11-2021 JHZ

export POLYGENE=$1
export cvt=${INF}/work/INF1.merge.cis.vs.trans

bedtools intersect -a <(cut -d, -f3,4,10,14 ${cvt} | \
                        tr ',' '\t' | \
                        grep -w ${POLYGENE} | \
                        awk -vFS="\t" -vOFS="\t" '{print "chr" $1,$2-1,$2}') \
                   -b <(sort -k1,1n -k2,2n ${INF}/csd3/glist-hg19 | \
                        grep -v -w -e X -e Y -e XY | \
                        awk '{$1="chr" $1;print}' | \
                        tr ' ' '\t') -wa -wb -loj | cut -f1-3,7 > polygene.txt
Rscript -e '
  library(dplyr)
  INF <- Sys.getenv("INF")
  POLYGENE <- Sys.getenv("POLYGENE")
  cvt <- read.csv(Sys.getenv("cvt"),as.is=TRUE) %>% filter(p.gene==POLYGENE) %>% mutate(p.chr=paste0("chr",p.chr))
  b <- cvt[c("p.chr","p.start","p.end","p.gene","cis.trans")]
  names(b) <- c("chr","start","end", "gene", "cistrans")
  d <- read.table("polygene.txt",as.is=TRUE,col.names=c("chr","start","end", "gene"))
  cols <- rep(12,nrow(b))
  cols[b[["cis"]]] <- 12
  a <- aggregate(d,by=list(with(d,chr),with(d,start),with(d,end)),FUN="c") %>% select(-c(chr,start,end))
  names(a) <- c("chr","start","end", "gene")
  labels <- bind_rows(a,distinct(select(b,-cistrans))) %>% filter(gene!=".")
  suppressMessages(library(circlize))
  pdf(file.path(INF,"hotspots",paste0("polygene-",POLYGENE,".pdf")))
  circos.par(start.degree = 90, track.height=0.1, cell.padding=c(0, 0, 0, 0))
  circos.initializeWithIdeogram(species="hg19", track.height=0.05, ideogram.height=0.06)
  circos.genomicLabels(labels,labels.column=4,side="inside")
  circos.genomicLink(a, b, col = cols, directional=2, border = 10, lwd = 2)
  circos.clear()
  dev.off()
  unlink("polygene.txt")
'

pdftopng -r 300 ${INF}/hotspots/polygene-${POLYGENE}.pdf ${POLYGENE}
mv ${POLYGENE}-000001.png ${INF}/hotspots/polygene-${POLYGENE}.png
