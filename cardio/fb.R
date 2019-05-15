# 15-5-2019 JHZ

jma <- read.delim("snps/cojo/INF1.jma",as.is=TRUE)
jma <- within(jma, {
       F <- ifelse(b<0, 1-freq, freq)
       B <- abs(b)
       MAF <- ifelse(freq < 0.5, freq, 1-freq)
})
pdf("work/fb.pdf")
with(jma, {
  print(cor(b,bJ))
  plot(MAF,B,cex=0.6)
  plot(b,bJ,cex=0.6)
})
dev.off()
