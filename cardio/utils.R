# 16-5-2019 JHZ

# Effect size ~ MAF plot

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

# MAF ~ N correspondence

MAF <- seq(0.01,0.1,0.01)
N <- seq(10,100,10)
PROD <- (MAF^2+MAF*(1-MAF))%o%N
colnames(PROD) <- N
rownames(PROD) <- MAF
PROD

# Comparison of GC lambda's between INTERVAL and INF1

INTERVAL <- read.table("work/INTERVAL.lambda.dat",as.is=TRUE)
INF1 <- read.table("work/INF1.lambda.dat",as.is=TRUE)
names(INTERVAL) <- names(INF1) <- c("prot","lambda")
lambda <- merge(INTERVAL,INF1,by="prot")
pdf("work/gc.lambda.pdf")
with(lambda,plot(lambda.x,lambda.y,cex=0.6))
dev.off()
