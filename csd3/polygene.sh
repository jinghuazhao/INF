# 15-10-2019 JHZ

export HOTSPOT=$1
export cvt=work/INF1.merge.cis.vs.trans
export DEBUG=0
(
  grep ${HOTSPOT} ${cvt} | \
  cut -d' ' -f3,4 | \
  awk -v OFS="\t" '{print "chr" $1,$2-1, $2}'
) > a1

(
  sort -k1,1n -k2,2n csd3/glist-hg19 | \
  grep -v X | \
  grep -v Y | \
  awk '{$1="chr" $1;print}' | \
  sed 's/ /\t/g'
) > a2

bedtools intersect -a a1 -b a2 -wa -wb -loj | \
cut  -f1-3,7 > a

(
  grep ${HOTSPOT} ${cvt} | \
  awk -v OFS="\t" '{print "chr" $14,$15,$16,$13,$17}' | \
  sed 's/\"//g'
) > b

R --no-save -q <<END
  HOTSPOT <- Sys.getenv("HOTSPOT")
  cvt <- Sys.getenv("cvt")
  cvt <- subset(read.table(cvt,as.is=TRUE,header=TRUE), p.gene==HOTSPOT)
  cvt <- within(cvt,p.chr <- paste0("chr",p.chr))
  write.table(cvt[c("p.chr","p.start","p.end","p.gene","cis.trans")],file="b",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  library(circlize)
  a <- read.table("b",as.is=TRUE,col.names=c("chr","start","end", "gene", "cistrans"))
  b <- read.table("a",as.is=TRUE,col.names=c("chr","start","end", "gene"))
  ab <- rbind(data.frame(b,cistrans="."),data.frame(unique(a[,-5]),cistrans="."))
  HOTSPOT <- Sys.getenv("HOTSPOT")
  pdf(paste0(HOTSPOT,".pdf"))
  circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
  circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  circos.genomicLabels(ab,labels.column = 4, side="inside")
  cols <- rep(12,nrow(a))
  cols[a["cistrans"]=="trans"] <- 10
  circos.genomicLink(b, a, col = cols, directional=1, border = 10, lwd = 2)
  circos.clear()
  dev.off()
END

pdftopng -r 300 ${HOTSPOT}.pdf ${HOTSPOT}
mv ${HOTSPOT}-000001.png ${HOTSPOT}.png

rm a1 a2 a b
