# 17-10-2019 JHZ

export HOTSPOT=$1
export cvt=work/INF1.merge.cis.vs.trans

(
  grep -w ${HOTSPOT} ${cvt} | \
  awk -v OFS="\t" '{print "chr" $3,$4-1, $4}'
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

R --no-save -q <<END
  HOTSPOT <- Sys.getenv("HOTSPOT")
  cvt <- Sys.getenv("cvt")
  cvt <- subset(read.table(cvt,as.is=TRUE,header=TRUE), SNP==HOTSPOT)
  cvt <- within(cvt,p.chr <- paste0("chr",p.chr))
  write.table(cvt[c("p.chr","p.start","p.end","p.gene","cis.trans")],file="b",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  library(circlize)
  b <- read.table("b",as.is=TRUE,col.names=c("chr","start","end", "gene", "cistrans"))
  cols <- rep(10,nrow(b))
  d <- read.table("a",as.is=TRUE,col.names=c("chr","start","end", "gene"))
  a <- aggregate(d,by=list(with(d,chr),with(d,start),with(d,end)),FUN="paste")[-c(4:6)]
  if (class(a[,4]) != "matrix") {
    a <- d
    cols[b["cistrans"]=="cis"] <- 12
  } else {
    a[,4] <- apply(a[,4],1,"paste",collapse=",")
    names(a) <- c("chr","start","end","gene")
  }
  labels <- rbind(b,data.frame(unique(a),cistrans="."))
  pdf(paste0(HOTSPOT,".pdf"))
  circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
  circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  circos.genomicLabels(labels,labels.column = 4, side="inside")
  circos.genomicLink(a, b, col = cols, directional=1, border = 10, lwd = 2)
  circos.clear()
  dev.off()
END

pdftopng -r 300 ${HOTSPOT}.pdf ${HOTSPOT}
mv ${HOTSPOT}-000001.png ${HOTSPOT}.png

rm a1 a2 a b
