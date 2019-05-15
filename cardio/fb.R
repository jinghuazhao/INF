# 15-5-2019 JHZ

jma <- read.delim("snps/cojo/INF1.jma",as.is=TRUE)
pdf("work/fb.pdf")
with(jma, {
  plot(freq,b,cex=0.6)
  plot(b,bJ,cex=0.6)
})
dev.off()
