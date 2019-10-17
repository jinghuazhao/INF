# 17-10-2019 JHZ

export POLYGENE=$1
export cvt=work/INF1.merge.cis.vs.trans
export DEBUG=0
(
  grep -w ${POLYGENE} ${cvt} | \
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
  grep ${POLYGENE} ${cvt} | \
  awk -v OFS="\t" '{print "chr" $14,$15,$16,$13,$17}' | \
  sed 's/\"//g'
) > b

R --no-save -q <<END
  POLYGENE <- Sys.getenv("POLYGENE")
  cvt <- Sys.getenv("cvt")
  cvt <- subset(read.table(cvt,as.is=TRUE,header=TRUE), p.gene==POLYGENE)
  cvt <- within(cvt,p.chr <- paste0("chr",p.chr))
  b <- cvt[c("p.chr","p.start","p.end","p.gene","cis.trans")]
  names(b) <- c("chr","start","end", "gene", "cistrans")
  library(circlize)
  d <- read.table("a",as.is=TRUE,col.names=c("chr","start","end", "gene"))
  cols <- rep(10,nrow(b))
  a <- aggregate(d,by=list(with(d,chr),with(d,start),with(d,end)),FUN="c")[-c(4:6)]
  if (class(a[,4]) != "list") {
     a <- d
     cols[with(b,cistrans)=="cis"] <- 12
  } else {
     a[,4] <- as.character(a[,4])
     names(a) <- c("chr","start","end", "gene")
  }
  labels <- rbind(a,unique(b[,-5]))
  pdf(paste0(POLYGENE,".pdf"))
  circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
  circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  circos.genomicLabels(labels,labels.column = 4, side="inside")
  circos.genomicLink(a, b, col = cols, directional=1, border = 10, lwd = 2)
  circos.clear()
  dev.off()
END

pdftopng -r 300 ${POLYGENE}.pdf ${POLYGENE}
mv ${POLYGENE}-000001.png ${POLYGENE}.png

rm a1 a2 a
